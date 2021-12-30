#!/bin/bash

set -e

blk=$1

dir_blk="/data/tge/Tian/UKBB_full/projects/genome_pred/software/ref_ld"
dir_1kg="/data/tge/Tian/UKBB_full/projects/genome_pred/g1000_eur/eur_qc"

plink \
    --bfile ${dir_1kg}/g1000_eur_qc \
    --keep-allele-order \
    --extract ${dir_blk}/snplist_ldblk/snplist_blk${blk} \
    --r square \
    --out ${dir_blk}/ldblk/ldblk${blk}_1kg \
    --memory 2000
