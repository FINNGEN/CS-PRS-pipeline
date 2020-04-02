#!/usr/bin/env python3

import argparse,os,gzip,itertools,subprocess,shlex,pickle
import numpy as np
import os,gzip
from collections import defaultdict as dd
from Tools.utils import file_exists,return_open_func,make_sure_path_exists,identify_separator,basic_iterator,get_path_info


def check_inputs(args):
    """
    Sanitzes inputs
    """
    open_func = return_open_func(args.file)

    # read header no matter what for debug perposes
    with open_func(args.file) as i: header = i.readline().strip()
    separator = identify_separator(header)
    header = header.split(separator)

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

    print(out_index,meta_index)
    return out_index,meta_index,header,separator

def parse_file(args):

    _,out_root,extension = get_path_info(args.file)
    open_func = return_open_func(args.file)

    out_index,meta_index,header,separator = check_inputs(args)
    write_func = gzip.open if args.gz else open
   
    if args.to_rsid:
        rsid_dict = load_rsid_mapping(args.out,args.map)
        out_file = os.path.join(args.out,out_root + '.rsid' )
        print(f"saving to {out_file}")
        if args.gz: out_file += '.gz'
        iterator = basic_iterator(args.file,skiprows=1)
        with write_func(out_file,'wt') as o:
            if not args.no_header:
                o.write(separator.join([header[i] for i in out_index]) + '\n')
            for line in iterator:
                snp,a1,a2 = [line[i] for i in meta_index]
                integers = ''.join((ch if ch.isdigit() else ' ') for ch in snp).split()
                if 'X' in snp:
                    chrom,pos = 'X',integers[0]
                else:
                    chrom,pos,*_ = integers
                    
                rsid = rsid_dict['_'.join([chrom,pos])]
                if rsid:
                    line[meta_index[0]] = rsid
                o.write(separator.join([line[i] for i in out_index]) + '\n')
            
    if args.to_chrompos:
        chrompos_dict = load_chrompos_mapping(args.out,args.map)
        out_file = os.path.join(args.out,out_root + '.chrompos' )
        print(f"saving to {out_file}")
        if args.gz: out_file += '.gz'
        iterator = basic_iterator(args.file,skiprows=1)
        with write_func(out_file,'wt') as o:
            if not args.no_header:
                o.write(separator.join([header[i] for i in out_index]) + '\n')
            for line in iterator:
                rsid,a1,a2 = [line[i] for i in meta_index]
                chrompos = chrompos_dict[rsid]
                if chrompos:
                    line[meta_index[0]] = f"chr{chrompos}_{a1}_{a2}"
                o.write(separator.join([line[i] for i in out_index]) + '\n')
    
def load_chrompos_mapping(out_path,file_map):
    '''
    Loads the finngen rsid to chrom_pos mapping
    '''
    _,map_root,_ = get_path_info(file_map)
    dict_path  = os.path.join(out_path,map_root + '.chrompos.pickle')
    if os.path.isfile(dict_path):
        print('pickling chrompos dict...')
        with open(dict_path,'rb') as i: chrompos_dict = pickle.load(i)
    else:
        print('generating chrompos dict...')
        chrompos_dict = dd(str)
        iterator = basic_iterator(file_map)
        for entry in iterator:
            rsid,chrom_pos = entry
            chrompos_dict[rsid] = chrom_pos
        with open(dict_path,'wb') as o:
            pickle.dump(chrompos_dict,o,protocol = pickle.HIGHEST_PROTOCOL)
            
    print('done.')
    return chrompos_dict

def load_rsid_mapping(out_path,file_map):
    '''
    Loads the chrompos to rsid mapping
    '''

    _,map_root,_ = get_path_info(file_map)
    dict_path  = os.path.join(out_path,map_root + '.rsid.pickle')     
    if os.path.isfile(dict_path):
        print('pickling rsid dict...')
        with open(dict_path,'rb') as i: rsid_dict = pickle.load(i)
    else:
        print('generating rsid dict...')
        rsid_dict = dd(str)
        iterator = basic_iterator(file_map)
        for entry in iterator:
            rsid,chrom_pos = entry
            rsid_dict[chrom_pos] = rsid
        with open(dict_path,'wb') as o:
            pickle.dump(rsid_dict,o,protocol = pickle.HIGHEST_PROTOCOL)
    print('done.')
    return rsid_dict

if __name__ == '__main__':

    """
    Pipeline to convert any file from/to rsids to chrompos notation
    """
    parser = argparse.ArgumentParser(description ="Process a file")

    #Basic inputs
    parser.add_argument("-o",'--out',type = str, help = "folder in which to save the results", required = True)
    parser.add_argument("--file",'-f',type = file_exists, help = "File to map", required = True)             
    parser.add_argument("--map",type = file_exists, help = "Mapping file", required = True)             
    parser.add_argument('--gz',action = 'store_true',help = ' Compress output file to gz',default = False)
    parser.add_argument('--no-header',action = 'store_true',help = 'Flag to use when no header is present',default = False)
    parser.add_argument('--metadata','-m',nargs = 3,metavar = ("ID","A1","A2"),required = True)    
    parser.add_argument('--columns',action = 'store',type = str,nargs = '*', default = 'all')
    conv_type = parser.add_mutually_exclusive_group(required = True)
    conv_type.add_argument('--to-rsid',action = 'store_true')
    conv_type.add_argument('--to-chrompos',action = 'store_true')
        
    args = parser.parse_args()
    print(args)
    make_sure_path_exists(args.out)
    parse_file(args)
