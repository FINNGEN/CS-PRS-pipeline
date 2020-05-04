#!/usr/bin/env python3

import argparse,os,subprocess,shlex
from utils import file_exists,make_sure_path_exists,tmp_bash,pretty_print,get_filepaths,basic_iterator,mem_mib,cpus

mem_mib = int(mem_mib)


def scores(args):

    scores_path = os.path.join(args.out,'scores')
    make_sure_path_exists(scores_path)

    # plink or pgen
    if args.bed: args.suffix += '.plink'
    elif args.pgen: args.suffix += '.pgen'


    # get final list of weights to run
    if args.weight:
        weight_iterator = [args.weight]
    else:
        # read list of files
        weight_iterator = basic_iterator(args.weight_list,columns= 0)
    
    weight_files = []
    for weight in weight_iterator:
        if os.path.isfile(weight):
            root_name = os.path.basename(weight).split('.weights')[0]
            print(root_name)
            score_file = os.path.join(scores_path,args.prefix + root_name + args.suffix)
            if os.path.isfile(score_file + '.sscore'):
                print(f'{score_file} already generated')
            weight_files.append((weight,score_file))
                    
    if not weight_files:
        print('no weights need to run, check inputs!')
        return
  
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
    for entry in weight_files:
        weight_file,score_file = entry
        pretty_print(weight_file)
        cmd = plink_cmd +  f" --memory {mem_mib} --score {weight_file} 2 4 6  header center --out {score_file}" 
        print(cmd)
        subprocess.call(shlex.split(cmd))
                

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description ="Calculation of PRS for summary stats.")

    parser.add_argument('--out', type=str, default = ".",help='Output Prefix')
   
    # requried
    weight_list_parser = parser.add_mutually_exclusive_group(required = True)
    weight_list_parser.add_argument('--weight-list',type = file_exists,help ='List of any weight file formatted as FILE_ROOT.chrCHROM.weights.etc')
    weight_list_parser.add_argument('--weight',type =file_exists,help ='path to single weight file')


    # BED OR PGEN
    scores_input = parser.add_mutually_exclusive_group(required = True)
    scores_input.add_argument('--bed',type = file_exists,help ='Path to plink bed file')
    scores_input.add_argument('--pgen',type = file_exists,help ='Path to pgen file')
        
    parser.add_argument('--prefix',type = str,help = "string to prepend to output",default = "")

    args = parser.parse_args()
    args.suffix = ".cs" 

    make_sure_path_exists(args.out)
    scores(args)
