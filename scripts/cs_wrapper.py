#!/usr/bin/env python3
import argparse,os.path,shlex,subprocess,sys,random
from collections import defaultdict as dd
from itertools import product
from utils import file_exists,make_sure_path_exists,tmp_bash,get_path_info,timing_function,natural_sort,get_filepaths,basic_iterator,pretty_print,mapcount,merge_files

root_path = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'
data_path = os.path.join(root_path,'data')

convert = os.path.join(root_path,'scripts','convert_rsids.py')
PRScs =  os.path.join(root_path,'PRScs','PRScs.py')

allele_map = {'A':'T','C':'G','G':'C','T':'A'}
allele_couples = [elem for elem in list(product(allele_map.keys(),allele_map.keys())) if elem[0] not in (elem[1],allele_map[elem[1]])] 

allele_couple_dict = {}
for ac in allele_couples:
    allele_couple_dict[ac] = [ac,ac[::-1],[allele_map[a] for a in ac],[allele_map[a] for a in ac[::-1]]]       

@timing_function
def to_rsid(args):
    """
    Fixes input sum stat to match cs_prs format.
    Converts chrom_pos_ref_alt summary stats into rsid.
    """
    munge_path = os.path.join(args.out,'munge')
    make_sure_path_exists(munge_path)
    out_file = os.path.join(munge_path,os.path.basename(args.sum_stats).replace('.gz','.rsid'))
    if os.path.isfile(out_file) and not args.force:
        print(f"{out_file} already esists")        
    else:
        cmd = f"python3 {convert} -f {args.sum_stats} -o {munge_path} --map {args.map}  --to-rsid --metadata SNP  --columns SNP A1 A2 BETA P "
        subprocess.call(shlex.split(cmd))

    args.sum_stats = out_file 


def weights(args):
    """
    Calculates weights with PRS-CS
    """
    pretty_print("WEIGHTS")
    log_path = os.path.join(args.out,'logs')
    make_sure_path_exists(log_path)
    # manipulate ref dir
    ref_dir,*_ = get_path_info(args.ref_file)
    print(f"ref dir: {ref_dir}")
    #manipulate bim file
    bim_dir,bim_root,_ = get_path_info(args.bim_file)
    bim_prefix = os.path.join(bim_dir,bim_root)
    print(f"bim file: {bim_prefix}")
    #mainpulate out dir
    args.ss_root= os.path.basename(args.sum_stats).split('.munged')[0]
    print(f"ss_root: {args.ss_root}")
    args.weights_path = os.path.join(args.out,'weights')
    make_sure_path_exists(args.weights_path)
    
    args.chrom_requested = sorted(args.chrom) if args.chrom else list(map(str,range(1,23)))
    print(f"requested chrom list :{' '.join(args.chrom_requested)}")
    if not args.force:
        #chromosomes to check to run
        chrom_list = []
        weights = natural_sort(get_filepaths(args.weights_path))
        for i in args.chrom_requested:
            # for each chrom to run check if weight already exists
            if not any(weight.endswith(f"chr{i}.txt") and os.path.basename(weight).startswith(args.ss_root) for weight in weights):
                chrom_list.append(str(i))
        if not chrom_list:
            print('All chromosomes ran')
            return
    else:
        chrom_list = args.chrom_requested

    print(f"final chrom list :{' '.join(chrom_list)}")
    if args.test:
        args.kwargs += " --n_iter=100"

    if args.parallel > 1:random.shuffle(chrom_list)
    
    logfile = os.path.join(log_path,f'{args.ss_root}.weights.log')
    log_template = os.path.join(log_path,f'{args.ss_root}.{{}}.weights.log')
    cmd = f"""parallel -j {args.parallel} "python3 -u {PRScs} --ref_dir {ref_dir} --bim_prefix {bim_prefix} --sst_file {args.sum_stats} --n_gwas {args.N} --out_dir {os.path.join(args.weights_path,args.ss_root)} {args.kwargs} --chrom {{}} >  {log_template} && echo {{}} " :::   {' '.join(chrom_list)} """
    print(cmd)
    tmp_bash(cmd)
    logs = [elem for elem in natural_sort(get_filepaths(log_path))]
    merge_files(logfile,logs)


#        cmd = f'python2.7 -u {PRScs} --ref_dir {ref_dir} --bim_prefix {bim_prefix} --sst_file {args.sum_stats} --n_gwas {args.N} --out_dir {os.path.join(args.weights_path,args.ss_root)} --chrom {",".join(chrom_list)} {args.kwargs} |& tee {logfile}'
 #       print(cmd)

    args.force = True
    
def to_chrompos(args):
    """
    Writes out weights in chrom_pos_ref_alt so we can use finngen data. 
    In order to do so, each rsid + ref/alt is "split" into 4 possible options based on whatever combination of min/maj allele in both strands. 
    """
    pretty_print("RSID FIX")
    print(args.chrom_requested)
    #chromosomes to check to run
    file_list = []
    weights = [elem for elem in natural_sort(get_filepaths(args.weights_path))]
    out_file = os.path.join(args.out,args.ss_root + '.weights.txt')
    if os.path.isfile(out_file) and mapcount(out_file) > 0 and not args.force:
        print(f"{out_file} already generated")
        return
    else:
        print(f"Saving to {out_file}")
        
    for weight in weights:
        # for each chrom to run check if weight already exists and if fixed weights exist
        if any(weight.endswith(f"chr{i}.txt") for i in args.chrom_requested):file_list.append(weight)

    
    with open(out_file,'wt') as o:
        for f in file_list:
            # convert from rsid to chrompos based on our data!
            chrompos_path,chrompos_root,_ = get_path_info(f)
            chrompos_file = os.path.join(chrompos_path,chrompos_root + '.chrompos')
            if not os.path.isfile(chrompos_file):
                cmd = f"python3 {convert} -f {f} --map {args.map} -o {args.weights_path} -m 1 3 4 --to-chrompos --no-header"
                subprocess.call(shlex.split(cmd))

            # create 4x entries with all possible REF_ALT combinations, also updating pos to build 38
            line_iterator = basic_iterator(chrompos_file)
            for line in line_iterator:
                snp = line[1]
                chrom,pos,ref,alt = snp.split('_')
                for a1,a2 in allele_couple_dict[(ref,alt)] :
                    line[2] = pos
                    line[1] = '_'.join([chrom,pos,a1,a2])
                    o.write('\t'.join(line) + '\n')
                                           
    return

  
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description ="Calculation of PRS for summary stats.")

    parser.add_argument('--ref-file',type = file_exists,help = 'Path to properly formatted hdf3 hm3 files',required = True)
    parser.add_argument('--bim-file',type = file_exists,help = 'Path to bim file',required = True)
    parser.add_argument('--sum-stats',type = file_exists,help = 'Path to summary stats',required = True)
    parser.add_argument('--N',type = int, help = 'Number of samples in the GWAS',required = True)
    parser.add_argument('--out','-o',help = "Out path",required = True)
    parser.add_argument('--kwargs',type = str,help = "Other args to pass to PRScs",default = "")
    parser.add_argument('--chrom',action = 'store',type = str,nargs = '*')
    parser.add_argument('--force',action = 'store_true')
    parser.add_argument('--test',action = 'store_true')
    parser.add_argument('--parallel',type = int,default =1)
    parser.add_argument('--map',type = file_exists,help = 'File that maps to/from rsids',required = True)
                        
    args = parser.parse_args()
    to_rsid(args)
    weights(args)
    to_chrompos(args)
