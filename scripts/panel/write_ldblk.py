#!/usr/bin/python


"""
Read LD blocks and write as hdf5
"""

import numpy as np
import h5py,os,sys

print(sys.argv[1])
BLK_DIR = sys.argv[1] 
LABEL = "1kg"

with open(os.path.join(BLK_DIR,'blk_chr')) as ff:
    blk_chr = [int(line.strip()) for line in ff]

with open(os.path.join(BLK_DIR,'blk_size')) as ff:
    blk_size = [int(line.strip()) for line in ff]

n_blk = len(blk_chr)
n_chr = max(blk_chr)

print(n_blk,n_chr)

for chrom in range(1,n_chr+1):

    #out file
    out_dir =  os.path.join(BLK_DIR,f'ldblk_{LABEL}_chr')
    os.makedirs(out_dir,exist_ok=True)
    chr_name = os.path.join(out_dir,f'ldblk_{LABEL}_chr' + str(chrom) + '.hdf5')
    print(f'... parse chomosome {chrom} {chr_name} ...')


    hdf_chr = h5py.File(chr_name, 'w')
    blk_cnt = 0
    # loop over blocks and 
    for blk in range(n_blk):
        if blk_chr[blk] == chrom: #use only block that belongs to the chrom
            ld = []; snplist = []
            if blk_size[blk] > 0: # check that block is not empty
                blk_root = os.path.join(BLK_DIR, 'ldblk','ldblk_' + str(blk+1) + f'_{LABEL}')
                #read actual ld matrix
                with open(blk_root + '.ld') as ff:
                    ld = [[float(val) for val in (line.strip()).split()] for line in ff]
                print(f'blk {blk+1} {np.shape(ld)}')

                #get blocksnplist  
                with open(blk_root + '.snplist') as ff:
                    snplist = [str(line.strip()) for line in ff]
            else:
                print(f'blk {blk+1} {blk_size[blk]}')

         
            blk_cnt += 1
            hdf_blk = hdf_chr.create_group('blk_%d' % blk_cnt)
            hdf_blk.create_dataset('ldblk', data=np.array(ld), compression="gzip", compression_opts=9)
#            hdf_blk.create_dataset('snplist', data=np.array(snplist,dtype=h5py.special_dtype(vlen=str)), compression="gzip", compression_opts=9)
            data = [elem.encode() for elem in snplist]
            hdf_blk.create_dataset('snplist', data=data, compression="gzip", compression_opts=9)
            
