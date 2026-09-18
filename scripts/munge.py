import argparse,os,gzip,subprocess,shlex
import numpy as np
import pandas as pd
from utils import file_exists,make_sure_path_exists,get_path_info,fix_header,identify_separator,mapcount_gzip,tmp_bash,mapcount,pretty_print,load_rsid_mapping,return_header,progressBar
from pathlib import Path

"""
Pandas-vectorized rewrite of munge_new.py. Same CLI, same intermediate file names/columns
(rsid_file / chrompos_file / rejected_1 / rejected_2), same lift.py contract, same final
5-column SNP A1 A2 BETA STAT output -- the row-by-row python loop over the sumstat is replaced
by a single bulk read + vectorized pandas/numpy operations, which is where the real cost lives
for GWAS-sized files (millions of rows).

Two deliberate, flagged behavior differences from munge_new.py:
 - an OR <= 0 (or otherwise producing a non-finite log) is now actually rejected as
   'or_problems'. In munge_new.py this case silently writes 'nan'/'inf' into BETA instead,
   because np.log doesn't raise -- the try/except there never actually fires for that case.
 - a sumstat row with the wrong number of columns is dropped by the C parser (on_bad_lines) and
   is not attributed to any rejection category, whereas munge_new.py's per-row try/except would
   catch it as 'format'. Only the alternate_parse-style digit-extraction failure is caught here.
"""


def _bulk_reject(handle,file_root,reason,df,cols):
    """
    Writes file_root/reason/<cols...> rows in bulk instead of one write() call per row.
    Column-wise string concatenation (not a row-wise .agg/apply, which is a python-level
    loop in disguise and was the actual bottleneck on reject-heavy chunks).
    """
    if df.empty: return
    lines = f"{file_root}\t{reason}\t" + df[cols[0]].astype(str)
    for col in cols[1:]:
        lines = lines + '\t' + df[col].astype(str)
    handle.write('\n'.join(lines) + '\n')


def _process_chunk(chunk,rename,use_pos_cols,effect_type):
    """
    Vectorized equivalent of one pass through munge_new.py's per-row loop, applied to a chunk
    of rows read out of the sumstat. Returns (rsid_out,chrompos_out,rejected_effect,rejected_or,
    rejected_format), each already in the column shape needed for writing/rejecting.
    """
    df = chunk.rename(columns = rename)

    df['a1_up'] = df['a1'].str.upper()
    df['a2_up'] = df['a2'].str.upper()

    # pd.to_numeric(errors='coerce') is only used to find which rows parse at all -- in this
    # pandas version it silently truncates precision on the actual values (verified against
    # plain float()/astype(float), which are exact), so the real values below are computed via
    # astype(float) on the already-validated subset instead of trusting to_numeric's output.
    effect_check = pd.to_numeric(df['effect'],errors = 'coerce')
    effect_valid = effect_check.notna() & np.isfinite(effect_check)
    rejected_effect = df.loc[~effect_valid]

    effect_num = pd.Series(np.nan,index = df.index,dtype = 'float64')
    effect_num.loc[effect_valid] = df.loc[effect_valid,'effect'].astype(float)

    is_rsid = df['variant'].str.contains('rs',regex = False,na = False)
    rsid_mask = is_rsid & effect_valid
    chrompos_mask = (~is_rsid) & effect_valid

    beta = np.log(effect_num) if effect_type == 'OR' else effect_num
    beta_valid = np.isfinite(beta)
    # OR <= 0 (non-finite log) is rejected here explicitly -- see module docstring
    or_problem_mask = (rsid_mask | chrompos_mask) & ~beta_valid
    rejected_or = df.loc[or_problem_mask]
    rsid_mask = rsid_mask & beta_valid
    chrompos_mask = chrompos_mask & beta_valid
    beta_str = beta.astype(str)

    rsid_out = pd.DataFrame({
        'snp': df.loc[rsid_mask,'variant'],
        'a1': df.loc[rsid_mask,'a1_up'],
        'a2': df.loc[rsid_mask,'a2_up'],
        'beta': beta_str.loc[rsid_mask],
        'stat': df.loc[rsid_mask,'stat'],
    })

    sub = df.loc[chrompos_mask]
    if use_pos_cols:
        chrom_vals,pos_vals = sub['in_chrom'],sub['in_pos']
        format_ok = pd.Series(True,index = sub.index)
    else:
        # same idea as alternate_parse: take the first two runs of digits out of the variant id
        extracted = sub['variant'].str.extract(r'(\d+)\D+(\d+)')
        chrom_vals,pos_vals = extracted[0],extracted[1]
        format_ok = chrom_vals.notna() & pos_vals.notna()
    rejected_format = sub.loc[~format_ok]

    chrom_ok,pos_ok = chrom_vals[format_ok],pos_vals[format_ok]
    chrompos_out = pd.DataFrame({
        'chr': chrom_ok,
        'snp': chrom_ok + '_' + pos_ok,
        'a1': sub.loc[format_ok,'a1_up'],
        'a2': sub.loc[format_ok,'a2_up'],
        'pos': pos_ok,
        'beta': beta_str.loc[chrompos_mask][format_ok],
        'stat': sub.loc[format_ok,'stat'],
    })

    return rsid_out,chrompos_out,rejected_effect,rejected_or,rejected_format


def parse_file(args):
    """
    Splits the sumstat into rsid_file (rsid-labeled rows, already PRScs-ready) and chrompos_file
    (chrom_pos-labeled rows, still needing liftover + a FinnGen rsid assigned in merge_files).
    """

    tmp_path = os.path.join(args.out,'tmp_parse')
    make_sure_path_exists(tmp_path)
    rej_path = os.path.join(tmp_path,'rejected_variants')
    make_sure_path_exists(rej_path)

    file_path,file_root,file_extension = get_path_info(args.ss)
    pretty_print(f"{file_root}",l=50)

    lines  = os.path.join(tmp_path,f'{file_root}.variantcount')
    if not os.path.isfile(lines) or not mapcount(lines): tmp_bash(f'zcat {args.ss} | wc -l > {lines}')
    total_lines  =  int(open(lines).read()) -1

    rsid_file = os.path.join(tmp_path,f'rsid_{file_root}.gz')
    chrompos_file = os.path.join(tmp_path,f'chrompos_{file_root}.gz')
    rej_log = os.path.join(tmp_path,'rejected_variants',f'rejected_1_{file_root}.gz')

    if os.path.isfile(rsid_file) and os.path.isfile(chrompos_file) and not args.force:
        print(str(total_lines) + ' variants already parsed')
        args.force = False
    else:
        args.force = True

    if args.force:
        header_fix = fix_header(args.ss)

        needed = [args.variant,args.ref,args.alt,args.effect,args.statistic]
        if not all(elem in header_fix for elem in needed):
            raise Exception(f"Missing columns in header: {[elem for elem in needed if elem not in header_fix]}")

        print(f"Using statistic type: {args.statistic_type} for column '{args.statistic}'")

        use_pos_cols = all(elem in header_fix for elem in [args.chrom,args.pos]) and all(elem != "NA" for elem in [args.chrom,args.pos])
        usecols = list(needed) + ([args.chrom,args.pos] if use_pos_cols else [])
        args.print(f'columns to parse: {usecols}, regular parse: {use_pos_cols}')

        sep = identify_separator(args.ss)
        read_kwargs = dict(header = None,skiprows = 1,names = header_fix,usecols = usecols,dtype = object,
                            keep_default_na = False,na_values = [''],on_bad_lines = 'skip')
        read_kwargs['sep'] = r'\s+' if sep == ' ' else sep

        rename = {args.variant:'variant',args.ref:'a1',args.alt:'a2',args.effect:'effect',args.statistic:'stat'}
        if use_pos_cols: rename.update({args.chrom:'in_chrom',args.pos:'in_pos'})
        reject_cols = ['variant','a1','a2','effect','stat'] + (['in_chrom','in_pos'] if use_pos_cols else [])

        if args.test:
            chunks = [pd.read_csv(args.ss,nrows = 10,**read_kwargs)]
        else:
            chunks = pd.read_csv(args.ss,chunksize = args.chunksize,**read_kwargs)
            n_chunks = -(-total_lines // args.chunksize)
            print(f'reading in {n_chunks} chunk(s) of up to {args.chunksize} rows')

        rsid_count = chrompos_count = rejected_count = processed = 0
        with gzip.open(rsid_file,'wt') as r,gzip.open(chrompos_file,'wt') as c,gzip.open(rej_log,'wt') as rej:
            r.write('\t'.join(['snp','a1','a2','beta','stat']) + '\n')
            c.write('\t'.join(['chr','snp','a1','a2','pos','beta','stat']) + '\n')

            for chunk in chunks:
                rsid_out,chrompos_out,rejected_effect,rejected_or,rejected_format = _process_chunk(chunk,rename,use_pos_cols,args.effect_type)

                rsid_out.to_csv(r,sep = '\t',index = False,header = False)
                chrompos_out.to_csv(c,sep = '\t',index = False,header = False)
                _bulk_reject(rej,file_root,'effect_missing',rejected_effect,reject_cols)
                _bulk_reject(rej,file_root,'or_problems',rejected_or,reject_cols)
                _bulk_reject(rej,file_root,'format',rejected_format,reject_cols)

                rsid_count += len(rsid_out)
                chrompos_count += len(chrompos_out)
                rejected_count += len(rejected_effect) + len(rejected_or) + len(rejected_format)
                processed += len(chunk)
                if not args.test: progressBar(processed,total_lines)

        if not args.test: print()
        print('done.')
        print(f'{rsid_count} variants labeled with rsid, {chrompos_count} labeled with chrompos, {rejected_count} rejected while parsing (of {total_lines} total)')
        if not args.test and rsid_count + chrompos_count + rejected_count == total_lines: print('SUCCESS: number of variant matches')

    if args.force:
        lift(args,chrompos_file)
        lifted_snps = mapcount_gzip(f"{chrompos_file}.lifted.gz") - 1
        print(f'{lifted_snps} chrompos variants survived liftover')


def merge_files(args,chrompos_to_rsid,fg_rsids):
    """
    Filters both branches down to FinnGen's set and writes the final 5-column
    SNP A1 A2 BETA STAT munged file.
    """
    pretty_print("MERGING",l = 20)
    tmp_path = os.path.join(args.out,'tmp_parse')
    file_path,file_root,file_extension = get_path_info(args.ss)

    rsid_file = os.path.join(tmp_path,f'rsid_{file_root}.gz')
    chrompos_file = os.path.join(tmp_path,f'chrompos_{file_root}.gz.lifted.gz')
    rej_log = os.path.join(tmp_path,'rejected_variants',f'rejected_2_{file_root}.gz')

    if args.prefix: args.prefix += "_"
    out_file = os.path.join(args.out,f"{args.prefix}{file_root}.munged.gz")

    if os.path.isfile(out_file) and not args.force:
        print(f'{out_file} already munged')
        return
    print(f"generating {out_file}")

    stat_col = 'P' if args.statistic_type == 'PVAL' else 'SE'

    if args.lift: column_names = ['beta','stat','lift_chr','lift_pos','REF','ALT']
    else: column_names = ['beta','stat','chr','pos','a1','a2']
    header = return_header(chrompos_file)

    if args.test:
        rsid_chunks = [pd.read_csv(rsid_file,sep = '\t',dtype = object,nrows = 30)]
        chrompos_chunks = [pd.read_csv(chrompos_file,sep = '\t',dtype = object,header = None,skiprows = 1,
                                        names = header,usecols = column_names,nrows = 30)]
    else:
        rsid_chunks = pd.read_csv(rsid_file,sep = '\t',dtype = object,chunksize = args.chunksize)
        chrompos_chunks = pd.read_csv(chrompos_file,sep = '\t',dtype = object,header = None,skiprows = 1,
                                       names = header,usecols = column_names,chunksize = args.chunksize)

    # no upfront total here on purpose: getting one would mean fully decompressing rsid_file/
    # chrompos_file an extra time before the merge even starts, which is expensive on real-sized
    # files -- just report a running count as chunks come in instead.
    final_variants = rsid_missing = no_rsid_mapping = processed = chunk_idx = 0
    with gzip.open(out_file,'wt') as o,gzip.open(rej_log,'wt') as rej:
        o.write('\t'.join(['SNP','A1','A2','BETA',stat_col]) + '\n')

        for chunk in rsid_chunks:
            # .map(set.__contains__) instead of .isin(fg_rsids): pandas' isin rebuilds a hashtable
            # from fg_rsids on every call, which is fine for a small set but disastrous when
            # fg_rsids has millions of entries (FinnGen's real rsid map) and this runs once per
            # chunk -- map() reuses the set's own O(1) lookup with no rebuild.
            keep = chunk['snp'].map(fg_rsids.__contains__)
            accepted,rejected = chunk.loc[keep],chunk.loc[~keep]
            accepted.rename(columns = {'snp':'SNP','a1':'A1','a2':'A2','beta':'BETA','stat':stat_col})[['SNP','A1','A2','BETA',stat_col]] \
                     .to_csv(o,sep = '\t',index = False,header = False)
            _bulk_reject(rej,file_root,'rsid_missing',rejected,['snp','a1','a2','beta','stat'])
            final_variants += len(accepted); rsid_missing += len(rejected); processed += len(chunk); chunk_idx += 1
            if not args.test: print(f'\rmerged {processed} rows so far (chunk {chunk_idx})',end = '',flush = True)

        for chunk in chrompos_chunks:
            chunk = chunk.rename(columns = dict(zip(column_names,['beta','stat','chrom','pos','a1','a2'])))
            chunk['key'] = chunk['chrom'].str.replace(r'\D','',regex = True) + '_' + chunk['pos']
            chunk['rsid'] = chunk['key'].map(chrompos_to_rsid)
            got_rsid = chunk['rsid'].notna() & (chunk['rsid'] != '')
            accepted,rejected = chunk.loc[got_rsid],chunk.loc[~got_rsid]
            accepted.rename(columns = {'rsid':'SNP','a1':'A1','a2':'A2','beta':'BETA','stat':stat_col})[['SNP','A1','A2','BETA',stat_col]] \
                     .to_csv(o,sep = '\t',index = False,header = False)
            _bulk_reject(rej,file_root,'no_rsid_mapping',rejected,['key','a1','a2','beta','stat'])
            final_variants += len(accepted); no_rsid_mapping += len(rejected); processed += len(chunk); chunk_idx += 1
            if not args.test: print(f'\rmerged {processed} rows so far (chunk {chunk_idx})',end = '',flush = True)

    if not args.test: print()
    original_variants = int(open(os.path.join(tmp_path,f'{file_root}.variantcount')).read()) - 1
    print(f'merge rejections: {rsid_missing} rsid_missing, {no_rsid_mapping} no_rsid_mapping')
    print(f'{final_variants}/{original_variants} variants retained in the final munged file ({final_variants/original_variants:.2%})')


def lift(args,chrompos_file):

    pretty_print("LIFTOVER",l=20)
    file_path,*_ = get_path_info(chrompos_file)
    lifted_file = f"{chrompos_file}.lifted.gz"
    if not os.path.isfile(lifted_file) or args.force:
        if args.lift:
            cmd = f"python3 {os.path.join(args.root_path,'lift','lift.py')} {chrompos_file} --chainfile {args.chainfile} --info chr pos a1 a2 --out {file_path}"
        else:
            cmd = f"cp {chrompos_file} {lifted_file}"
            print('chainfile missing or empty, no lifting will take place')
        print(cmd)
        subprocess.call(shlex.split(cmd))
    else:
        print('already lifted file')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description ="Munge a GWAS summary stat (pandas-vectorized).")

    parser.add_argument("-o",'--out',type = str, help = "folder in which to save the results", required = True)
    parser.add_argument("--ss", help = "Path to gwas summary stat.", required = True,type = file_exists)
    parser.add_argument("--rsid-map", help = "Path to rsid to chrompos tsv mapping (or its pre-built pickle).", required = True,type = file_exists)
    parser.add_argument("--chainfile", help = "Path to liftover chainfile.")

    parser.add_argument('--test',action = 'store_true',help = 'Flag for testing purposes.')
    parser.add_argument('--force',action = 'store_true',help = 'Flag for forcing re-run.')
    parser.add_argument('--chunksize',type = int,default = 500_000,help = 'Rows read/processed per chunk (progress bar granularity).')

    parser.add_argument('--effect_type',  type = lambda s : s.upper(), choices=['BETA','OR'])
    parser.add_argument('--variant', type=str,required=True,help = 'Column entry of variant id')
    parser.add_argument('--ref', type=str,required=True,help='Column entry of ref (effect)')
    parser.add_argument('--alt', type=str,required=True,help = 'Column entry of other allel')
    parser.add_argument('--effect', type=str,required=True,help='Column entry of effect column (beta/OR)')
    parser.add_argument('--statistic', type=str,required=True,help='Column entry of statistic (p-value or SE)')
    parser.add_argument('--statistic-type', type = lambda s : s.upper(), choices=['SE','PVAL'],required=True,help='Whether --statistic is a standard error or a p-value')
    parser.add_argument('--chrom', type=str,help='Column entry of chrom')
    parser.add_argument('--pos', type=str,help='Column entry of position')
    parser.add_argument('--prefix',type = str,help = "string to prepend to output",default = "")

    args = parser.parse_args()
    args.ss = os.path.abspath(args.ss)

    make_sure_path_exists(args.out)

    args.lift = True
    if not args.chainfile or os.path.getsize(args.chainfile) == 0:
        args.lift = False

    if args.test:
        def vprint(x):
            print(x)
    else:
        vprint = lambda *a: None

    args.print = vprint
    args.print(args)
    args.root_path  = Path(os.path.realpath(__file__)).parent.absolute()

    chrompos_to_rsid = load_rsid_mapping(args.rsid_map,inverse = True)
    fg_rsids = set(chrompos_to_rsid.values())

    parse_file(args)
    merge_files(args,chrompos_to_rsid,fg_rsids)
