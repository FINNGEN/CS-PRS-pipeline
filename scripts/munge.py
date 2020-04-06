import argparse,os,pickle,gzip,itertools
import numpy as np
from utils import file_exists,make_sure_path_exists,basic_iterator,return_header,get_path_info,return_open_func,fix_header,isfloat,mapcount_gzip,progressBar,tmp_bash,mapcount,pretty_print,load_rsid_mapping
from functools import partial
from collections import defaultdict as dd


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
            return
        else:
            print('mismatch')

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
                if 'rs' in variant: # RSID FORMAT --> check if rsid exsists in our data
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

    counts = [mapcount_gzip(elem) for elem in [rsid_file,chrompos_file,rej_log]]
    if np.sum(counts) - 2 == total_lines: print('SUCCESS: number of variant matches')



    #### ADD LIFTING HERE ###
    # after lifting, i just need to "merge" the two files

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
    parse_file(args)
