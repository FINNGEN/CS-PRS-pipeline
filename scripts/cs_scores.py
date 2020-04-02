#!/usr/bin/env python3

import argparse,os,subprocess,shlex
from Tools.utils import file_exists,make_sure_path_exists,tmp_bash,pretty_print,get_filepaths,basic_iterator,mem_mib,cpus

mem_mib = int(mem_mib)



def scores(args):

    scores_path = os.path.join(args.out,'scores')
    make_sure_path_exists(scores_path)

    # get final list of weights to run
    if args.weight_path:
        weight_iterator = ( f for f in get_filepaths(args.weight_path) if f.endswith('.weights.txt'))
    else:
        # read list of files
        weight_iterator = basic_iterator(args.weight_list,columns= 0)

    if args.bed: args.suffix += '.plink'
    elif args.pgen: args.suffix += '.pgen'
    
    weight_files = []
    for weight in weight_iterator:
        if os.path.isfile(weight):
            root_name = os.path.basename(weight).split('.weights')[0]
            score_file = os.path.join(scores_path,root_name + args.suffix)
            if os.path.isfile(score_file + '.sscore'):
                print(f'{score_file} already generated')
            weight_files.append((weight,score_file))
                    
    if not weight_files:
        raise ValueError('weight list empty')

  
    # fix plink command for bed file
    if args.bed:
        file_root = args.bed.split('.bed')[0]
        plink_cmd = f"plink2 --bfile {file_root} "
    elif args.pgen:
        file_root = args.pgen.split('.pgen')[0]
        plink_cmd = f"plink2 --pfile {file_root} "            
    else:
        raise TypeError("wrong inputs")

    # check if freq file is provided
    if not args.freq_file:
        args.freq_file = file_root + '.afreq'
        #check if it's in the same path
        if not os.path.isfile(args.freq_file):
            print(f"{args.freq_file} missing")
            subprocess.call(shlex.split(f"{plink_cmd} --freq --out {file_root}"))
        else:
            print(f"{args.freq_file} found")

    plink_cmd += f' --read-freq {args.freq_file} '

    # now run for all weight_files
    for entry in weight_files:
        weight_file,score_file = entry
        pretty_print(weight_file)
        cmd = plink_cmd +  f" --memory {mem_mib} --score {weight_file} 2 4 6  header center --out {score_file}" 
        print(cmd)
        subprocess.call(shlex.split(cmd))
                

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description ="Calculation of PRS for summary stats.")

    parser.add_argument('--out', type=str, default = ".",
                    help='Output Prefix')
   
    # requried
    weight_list_parser = parser.add_mutually_exclusive_group()
    weight_list_parser.add_argument('--weight-list',type = file_exists,help ='List of any weight file formatted as FILE_ROOT.chrCHROM.weights.etc')
    weight_list_parser.add_argument('--weight-path',type = str,help ='path where weights are  located in format XXX.chrCHR.weightsXXX')

    parser.add_argument('--freq-file',help ='path to freq file',type = file_exists)
    parser.add_argument('--extract',help ='path to snps to extract',type = file_exists)

    # BED OR PGEN
    scores_input = parser.add_mutually_exclusive_group(required = True)
    scores_input.add_argument('--bed',type = file_exists,help ='Path to plink bed file')
    scores_input.add_argument('--pgen',type = file_exists,help ='Path to pgen file')
        
    parser.add_argument('--suffix',type = str,help = "string to append to output",default = "")

    args = parser.parse_args()
    args.suffix = ".cs." + args.suffix

    make_sure_path_exists(args.out)
    scores(args)
