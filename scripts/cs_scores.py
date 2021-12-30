#!/usr/bin/env python3

import argparse,os,subprocess,shlex
from utils import file_exists,make_sure_path_exists,tmp_bash,pretty_print,get_filepaths,basic_iterator,mem_mib,cpus,get_path_info,mapcount

mem_mib = int(mem_mib)

def scores(args):

    scores_path = os.path.join(args.out,'scores')
    make_sure_path_exists(scores_path)

    # get final list of weights to run
    if args.weight:
        weight_iterator = [args.weight]
    else:
        # read list of files
        weight_iterator = basic_iterator(args.weight_list,columns= 0,separator = '\t')
    
    weight_files = []
    for weight in weight_iterator:
        if os.path.isfile(weight):
            root_name = os.path.basename(weight).split('.weights')[0]
            score_file = os.path.join(scores_path, root_name + args.suffix)
            weight_files.append((weight,score_file))  
  
    # fix plink command for bed file
    if args.bed:
        plink_root = args.bed.split('.bed')[0]
        plink_cmd = f"plink2 --bfile {plink_root} "
    elif args.pgen:
        plink_root = args.pgen.split('.pgen')[0]
        plink_cmd = f"plink2 --pfile {plink_root} "            
    else:
        raise TypeError("wrong inputs")

    # check if freq file is provided
    freq_file = plink_root + '.afreq'
    if os.path.isfile(freq_file):
        print('freq file present')
        plink_cmd += f' --read-freq {freq_file} '


    # now run for all weight_files
    for i,entry in enumerate(weight_files):
        weight_file,score_file = entry
        print(f"{i+1}/{len(weight_files)} {weight_file}")
        
        cmd = plink_cmd +  f" --out {score_file} --score {weight_file} 2 4 6  header center list-variants ignore-dup-ids  --memory {mem_mib} --silent "
        if not os.path.isfile(score_file + '.sscore'):
            subprocess.call(shlex.split(cmd))
        else:print(f'{score_file} already generated')

        if args.region and mapcount(args.region):
            _,region_root,file_extension = get_path_info(args.region)
            score_file += '.no_' + region_root
            if not os.path.isfile(score_file + '.sscore'):
                
                cmd =  plink_cmd +  f" --memory {mem_mib} --score {weight_file} 2 4 6  center list-variants --exclude range {args.region} --out {score_file} " 
                print(cmd)
                subprocess.call(shlex.split(cmd))
            else:print(f'{score_file} already generated')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description ="Calculation of PRS for summary stats.")

    parser.add_argument('--out', type=str, default = ".")
    parser.add_argument('--region',type = file_exists,help ='Path to list of regions to exclude')

    # required
    weight_list_parser = parser.add_mutually_exclusive_group(required = True)
    weight_list_parser.add_argument('--weight-list',type = file_exists,help ='List of any weight file formatted as FILE_ROOT.chrCHROM.weights.etc')
    weight_list_parser.add_argument('--weight',type =file_exists,help ='path to single weight file')
    
    # BED OR PGEN
    scores_input = parser.add_mutually_exclusive_group(required = True)
    scores_input.add_argument('--bed',type = file_exists,help ='Path to plink bed file')
    scores_input.add_argument('--pgen',type = file_exists,help ='Path to pgen file')


    args = parser.parse_args()
    args.suffix = ""

    make_sure_path_exists(args.out)
    scores(args)
