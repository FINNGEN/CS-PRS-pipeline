#!/usr/bin/env python3

import argparse,os,gzip,itertools,subprocess,shlex,pickle,csv
import numpy as np
import os,gzip
from collections import defaultdict as dd
from utils import file_exists,return_open_func,make_sure_path_exists,identify_separator,basic_iterator,get_path_info,return_header,load_rsid_mapping


def check_inputs(args):
    """
    Sanitzes inputs
    """
    open_func = return_open_func(args.file)
    # read header no matter what for debug perposes
    separator = identify_separator(args.file)
    header = return_header(args.file)

    #check if metadata is integers
    int_meta_check = np.prod([elem.isdigit() for elem in args.metadata])

    # adjust columns if one wants to keep all file
    if args.columns == 'all':
        if int_meta_check:
            args.columns = list(map(str,range(len(header))))
        else:
            args.columns = header

    # check if output columns are integers
    int_col_check=  np.prod([elem.isdigit() for elem in args.columns])

    #both meta and out columns need to be consistent
    if int_meta_check ^ int_col_check:
        raise ValueError('both metadata and columns need to be of the same type')

    # check that index is within limit
    if int_meta_check and int_col_check:
        args.columns = [int(elem) for elem in args.columns]
        for index in args.columns:
            assert index < len(header),f'{index} column index out of bounds'
        out_index = args.columns
        args.metadata = [int(elem) for elem in args.metadata]
        for index in args.metadata:
            assert index < len(header),f'{index} metadata index out of bounds'
        meta_index = args.metadata

    # check that columns and metadat names are in header
    if (not int_meta_check) and (not int_col_check):
        for column in args.columns:
            assert column in header,f"{column} column not in header"
        for column in args.metadata:
            assert column in header,f"{column} metadata not in header"      

        # index of columns to write out
        out_index = [header.index(elem) for elem in args.columns]
        # index of columns to read
        meta_index = [header.index(elem) for elem in args.metadata]

    print('columns to keep',out_index)
    print('metadata columns',meta_index)
    return out_index,meta_index,header,separator

def parse_file(args):

    _,out_root,extension = get_path_info(args.file)

    out_index,meta_index,header,separator = check_inputs(args)
    write_func = gzip.open if args.gz else open

    out_file = os.path.join(args.out,out_root + ".CONVERT")
    if args.gz: out_file += '.gz'

    #if there is no header, don't skip first row when parsing file
    skip = 0 if args.no_header else 1
    iterator = basic_iterator(args.file,skiprows=skip)
    if args.to_rsid:
        rsid_dict = load_rsid_mapping(args.map,inverse = True)#chrompos --> rsid
        out_file = out_file.replace('CONVERT','rsid')
        print(f"saving to {out_file}")
        with write_func(out_file,'wt') as o:
            if not args.no_header:# write header,since it cannot be parsed
                o.write(separator.join([header[i] for i in out_index]) + '\n')
            for line in iterator:
                #extract chrom/pos from snp id (i.e. the first element of meta_index)
                snp = line[meta_index[0]]
                # return all integers in snp string. i assume the first 2 are chrom/pos
                #this handles mixed sumstats where rsids are present
                if 'rs' in snp:
                    o.write(separator.join([line[i] for i in out_index]) + '\n')
                    continue
                # this handles weird missing lines
                integers = ''.join((ch if ch.isdigit() else ' ') for ch in snp).split()
                if not integers:
                    continue
                    
                if 'X' in snp: # unless it's chrom X
                    chrom,pos = 'X',integers[0]
                else:
                    try:
                        chrom,pos,*_ = integers
                    except:
                        print(snp)
                        print(integers)
                    
                rsid = rsid_dict['_'.join([chrom,pos])]
                if rsid: # replace snpid with rsid
                    line[meta_index[0]] = rsid
                o.write(separator.join([line[i] for i in out_index]) + '\n')
            
    if args.to_chrompos:
        chrompos_dict = load_rsid_mapping(args.map) # risd -- >chrompos
        out_file = out_file.replace('CONVERT','chrompos')
        print(f"saving to {out_file}")
        with write_func(out_file,'wt') as o:
            if not args.no_header:
                o.write(separator.join([header[i] for i in out_index]) + '\n')
            for line in iterator:
                rsid,a1,a2 = [line[i] for i in meta_index]
                chrompos = chrompos_dict[rsid]
                if chrompos:
                    line[meta_index[0]] = f"chr{chrompos}_{a1}_{a2}"
                o.write(separator.join([line[i] for i in out_index]) + '\n')


if __name__ == '__main__':

    """
    Pipeline to convert any file from/to rsids to chrompos notation
    """
    parser = argparse.ArgumentParser(description ="Process a file")

    #Basic inputs
    parser.add_argument("-o",'--out',type = str, help = "folder in which to save the results", required = True)
    parser.add_argument("--file",'-f',type = file_exists, help = "File to map", required = True)             
    parser.add_argument("--map",type = file_exists, help = "Mapping file from rsid to chrompos", required = True)             
    parser.add_argument('--gz',action = 'store_true',help = ' Compress output file to gz',default = False)
    parser.add_argument('--no-header',action = 'store_true',help = 'Flag to use when no header is present',default = False)
    parser.add_argument('--metadata','-m',nargs = '*',required = True,help='columns required for parsing. Should be SNPID columns for to-rsid and SNPID,A1,A2 to go to chrompos')    
    parser.add_argument('--columns',help = 'column that need to be kept, either numerical integers or column names',action = 'store',type = str,nargs = '*', default = 'all')
    
    conv_type = parser.add_mutually_exclusive_group(required = True)
    conv_type.add_argument('--to-rsid',action = 'store_true')
    conv_type.add_argument('--to-chrompos',action = 'store_true')
        
    args = parser.parse_args()
    print(args)
    make_sure_path_exists(args.out)
    args.file = os.path.abspath(args.file)
    parse_file(args)
