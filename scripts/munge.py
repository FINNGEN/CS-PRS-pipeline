import argparse,os,pickle,gzip,itertools,subprocess,shlex
import numpy as np
from utils import file_exists,make_sure_path_exists,basic_iterator,get_path_info,fix_header,isfloat,mapcount_gzip,progressBar,tmp_bash,mapcount,pretty_print,load_rsid_mapping,load_pos_mapping,map_alleles
from functools import partial
from pathlib import Path
from collections import defaultdict as dd



def merge_files(args):
    """
    Function that takes the rsid and chrompos lifted file and prints out the final variants, checking if the variant exists in our dataset (also checking for strand flip).
    """
    tmp_path = os.path.join(args.out,'tmp_parse')
    rej_path = os.path.join(tmp_path,'rejected_variants')
    file_path,file_root,file_extension = get_path_info(args.ss)

    rsid_file = os.path.join(tmp_path,f'rsid_{file_root}.gz')
    chrompos_file =  os.path.join(tmp_path,f'chrompos_{file_root}.gz.lifted.gz')
    rej_log = os.path.join(tmp_path,'rejected_variants',f'rejected_2_{file_root}.gz')
    
    out_file = os.path.join(args.out,f"{file_root}.munged.gz")
    if os.path.isfile(out_file) and not args.force:
        print(f'{out_file} already munged')
        return

    pos_dict = load_pos_mapping(args.chrompos_map,args.out)
    rsid_positions = len(pos_dict.keys())
    # starting final merge pass
    final_variants = 0.
    with gzip.open(rej_log,'wt') as rej,gzip.open(out_file,'wt') as o:
        out_header = '\t'.join(['CHR','SNP','A1','A2','BP','OR','P'])
        o.write(out_header + '\n')
        #looping of rsid file
        iterator = basic_iterator(rsid_file,skiprows = 1)
        loop = itertools.islice(iterator,30) if args.test else iterator
        for entry in loop:
            args.print(entry)
            chrom,rsid,a1,a2,pos,OR,pval = entry
            pass_bool,out_line = process_variant(pos_dict,chrom,pos,a1,a2,OR,pval,file_root,rsid)
            final_variants += pass_bool
            out_file = o if pass_bool else rej
            out_file.write(out_line)  
        #looping of lifted file
        iterator = basic_iterator(chrompos_file,skiprows = 1)
        loop = itertools.islice(iterator,30) if args.test else iterator
        for entry in loop:
            args.print(entry)
            *_,OR,pval,chrom,pos,a1,a2 = entry
            chrom = ''.join([s for s in chrom if s.isdigit()]) # extract integer from chrom field
            variant_id =  f"{chrom}_{pos}_{a1}_{a2}"
            pass_bool,out_line = process_variant(pos_dict,chrom,pos,a1,a2,OR,pval,file_root,variant_id)
            final_variants += pass_bool
            out_file = o if pass_bool else rej
            out_file.write(out_line)            

        # count lines
        if not args.test:
            original_variants  =  float(open(os.path.join(tmp_path,f'{file_root}.variantcount')).read()) -1  
            print('compared to input sumstat:',original_variants,final_variants,final_variants/original_variants)
            print('compared to rsid positions in FG:',rsid_positions,final_variants,final_variants/rsid_positions)
            
def process_variant(pos_dict,chrom,pos,a1,a2,OR,pval,file_name,variant_id):
    """
    Function that parses the lines of the chrompos and rsid files. It makes sure that the variant in each line is (potentially) the same as Finngen's by trying all possible combinations of strand flip and direction.
    """
    
    finngen_variant = pos_dict[f"{chrom}_{pos}"]
    # the position exists in finngen data
    if finngen_variant:
        # standard out line in case of missing position
        rej_line = '\t'.join([file_name,'wrong_alleles',variant_id,a1,a2,OR,pval,chrom,pos])+'\n'
        variant_check = map_alleles(a1,a2)
        for finngen_ref,finngen_alt in finngen_variant:
            if map_alleles(finngen_ref,finngen_alt) == variant_check:
                finngen_snp = f"chr{chrom}_{pos}"
                out_line = '\t'.join([chrom,finngen_snp,a1,a2,pos,OR,pval]) + '\n'
                return True,out_line
       
    # in case of position missing
    else:
        rej_line = '\t'.join([file_name,'missing_position',variant_id,a1,a2,OR,pval,chrom,pos]) + '\n'
           
    return False,rej_line



#####################################
#-----------FIRST PASS -------------#
#####################################
def parse_file(args):
    """
    Function that loops through a gwas assoc file. Based on the fields that are passed the file is parsed and split into two subfiles:
    1) containing the variants labeled with rsids and if the rsid is in the finngen data (via rsid_dict)
    2) containing the variants for which we have chrom_pos notation and need lifting
    In both cases the effect is mapped to OR if needed. Variants with no beta/or are skipped, as well as ones whose rsid is not in our data or whose format is non standard.
    """

    tmp_path = os.path.join(args.out,'tmp_parse')
    make_sure_path_exists(tmp_path)
    rej_path = os.path.join(tmp_path,'rejected_variants')
    make_sure_path_exists(rej_path) 

    file_path,file_root,file_extension = get_path_info(args.ss)
    pretty_print(f"{file_root}")

    lines  = os.path.join(tmp_path,f'{file_root}.variantcount')
    if not os.path.isfile(lines) or not mapcount(lines): tmp_bash(f'zcat {args.ss} | wc -l > {lines}')
    total_lines  =  int(open(lines).read()) -1  
    
    # define output files (chrompos,rsid,rejected)
    rsid_file = os.path.join(tmp_path,f'rsid_{file_root}.gz')
    chrompos_file = os.path.join(tmp_path,f'chrompos_{file_root}.gz')
    rej_log = os.path.join(tmp_path,'rejected_variants',f'rejected_1_{file_root}.gz')
    # if files already exist, check if number of variants matches.
    if os.path.isfile(rsid_file) and os.path.isfile(chrompos_file) and not args.force:
        print(str(total_lines) + ' variants already parsed')
        counts = [mapcount_gzip(elem) for elem in [rsid_file,chrompos_file,rej_log]]
        if np.sum(counts) - 2 == total_lines:
            print('SUCCESS: number of variant matches')
            args.force = False
        else:
            print(counts,np.sum(counts))
            args.force = True
            print('mismatch')
    else:
        args.force = True
        
    if args.force:
        # check if stats have OR or BETA
        if args.effect_type  == 'OR':
            def or_func(x): return x
        elif args.effect_type =='BETA':
            def or_func(x):
                return str(float(np.exp(np.float128(x))))

        # fix headers that have extra spaces
        header_fix = fix_header(args.ss)
        args.print(f'header: {header_fix}')
        
        # relevant columns to parse
        columns = [args.variant,args.ref,args.alt,args.effect,args.pval]
        args.print(f'columns to parse: {columns}')
        if not all([elem in header_fix for elem in columns]):
            raise Exception(f"Missing columns in header: {[elem for elem in columns if elem not in header_fix]}")
       
        indexes  = [header_fix.index(elem) for elem in columns]
        # define parsing function based on inputs
        pos_indexes = [header_fix.index(elem) for elem in [args.chrom,args.pos] if elem]
        if len(pos_indexes) ==2:
            indexes += pos_indexes
            parse_func = partial(regular_parse,or_func = or_func)
        else:
            parse_func = partial(alternate_parse,or_func = or_func)
        
            args.print(f'indexes: {indexes}')

            # WE ARE NOW READY TO PARSE THE FILE
            rsid_dict = load_rsid_mapping(args.rsid_map,args.out)
            
        with gzip.open(rsid_file,'wt') as r,gzip.open(chrompos_file,'wt') as c,gzip.open(rej_log,'wt') as rej:
            out_header = '\t'.join(['chr','snp','a1','a2','pos','or','p'])
            c.write(out_header + '\n')
            r.write(out_header + '\n')

            # loop through file, but only do first 30 lines if in test mode
            iterator = basic_iterator(args.ss,skiprows = 1,columns = indexes)
            loop = itertools.islice(iterator,10) if args.test else iterator
            for i,info in enumerate(loop):
                if not i % 1000:progressBar(i,total_lines)

                variant,a1,a2,effect,pval,*chrompos = info            
                args.print(f"Reading columns: {info}")

                # check if effect can be mapped to float (might be missing)
                if not isfloat(effect):
                    args.print(effect)
                    out_line,out_file = '\t'.join([file_root,'effect_missing'] + info),rej

                else:
                    if 'rs' in variant: # RSID FORMAT --> check if rsid exsists in our data. Update position!
                        chrompos = rsid_dict[variant]
                        if chrompos:
                            chrom,pos = chrompos.split('_')
                            OR = str(or_func(float(effect)))
                            out_line,out_file = '\t'.join([chrom,variant,a1.upper(),a2.upper(),pos,OR,pval]),r
                        else:
                            out_line,out_file = '\t'.join([file_root,'rsid_missing'] + info),rej

                    else: # NOT RSID --> get chrom/pos info if available, else retrieve it from snpid
                        try:
                            out_line,out_file  = '\t'.join(parse_func(info)),c
                        except:
                            out_line,out_file = '\t'.join([file_root,'format'] + info),rej

                args.print(f"output:{out_line}")
                out_file.write(out_line +'\n')
        print('done.')
        if not args.test:
            counts = [mapcount_gzip(elem) for elem in [rsid_file,chrompos_file,rej_log]]
            if np.sum(counts) - 2 == total_lines: print('SUCCESS: number of variant matches')

    lift(chrompos_file,args.chainfile,args.force)


def lift(chrompos_file,chainfile,force):

    file_path,*_ = get_path_info(chrompos_file)
    if not os.path.isfile(f"{chrompos_file}.lifted.gz") or force:
        cmd = f"python3 {os.path.join(args.root_path,'lift','lift.py')} {chrompos_file} --chainfile {chainfile} --info chr pos a1 a2 --out {file_path}"
        subprocess.call(shlex.split(cmd))
    else:
        print('already lifted file')
        
def regular_parse(info,or_func):
    """
    Standard parsing function when chrom & pos are provided. It just converts OR to BETA if needed.
    """
    variant,a1,a2,effect,pval,chrom,pos = info
    OR = or_func(effect)
    return chrom,f"{chrom}_{pos}",a1.upper(),a2.upper(),pos,OR,pval
    
def alternate_parse(info,or_func):
    """
    Parsing if chrom & pos are missing, extract the first two integers from the variant string
    """

    variant,a1,a2,effect,pval,*_  = info
    chrom,pos,*_ = ''.join((ch if ch.isdigit() else ' ') for ch in variant).split()
    OR = or_func(effect)
    return chrom,f"{chrom}_{pos}",a1.upper(),a2.upper(),pos,OR,pval

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description ="Munge a GWAS summary stat.")


    #Basic inputs
    parser.add_argument("-o",'--out',type = str, help = "folder in which to save the results", required = True)
    parser.add_argument("--ss", help = "Path to gwas summary stat.", required = True,type = file_exists)
    parser.add_argument("--rsid-map", help = "Path to rsid to chrompos tsv mapping.", required = True,type = file_exists)
    parser.add_argument("--chrompos-map", help = "Path to gwas summary stat.", required = True,type = file_exists)
    parser.add_argument("--chainfile", help = "Path to liftover chainfile.", required = True,type = file_exists)

    parser.add_argument('--test',action = 'store_true',help = 'Flag for testing purposes.')
    parser.add_argument('--force',action = 'store_true',help = 'Flag for forcing re-run.')

    # COLUMNS OF RELEVANT FIELDS
    parser.add_argument('--effect_type', type=str, choices=['BETA','OR'])
    parser.add_argument('--variant', type=str,required=True,help = 'Column entry of variant id')
    parser.add_argument('--ref', type=str,required=True,help='Column entry of ref (effect)')
    parser.add_argument('--alt', type=str,required=True,help = 'Column entry of other allel')
    parser.add_argument('--effect', type=str,required=True,help='Column entry of effect column (beta/OR)')
    parser.add_argument('--pval', type=str,required=True,help='Column entry of pvalue')
    parser.add_argument('--chrom', type=str,help='Column entry of chrom')
    parser.add_argument('--pos', type=str,help='Column entry of position')

    
    
    args = parser.parse_args()
    args.ss = os.path.abspath(args.ss)
    
    make_sure_path_exists(args.out)
    
    # test cases
    if args.test:
        def vprint(x):
            print(x)
    else:
        vprint = lambda *a: None  # do-nothing function

    args.print = vprint
    args.print(args)
    args.root_path  = Path(os.path.realpath(__file__)).parent.absolute()
    parse_file(args)
    merge_files(args)
