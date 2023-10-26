#!/bin/bash

DIR="/mnt/disks/data/panel/prs/hm3/" #path where to output/work
PANEL_DIR="/mnt/disks/data/panel/prs/sisuv4_panel_hm3" #path to root of bed file of panel
LABEL="1kg" #output label
BLOCKS="/home/pete/Dropbox/Projects/CS-PRS-pipeline/scripts/panel/blocks.bed" # bed file with regions

cwd=$(pwd)

cd $DIR
rm -rf ldblk && mkdir ldblk && rm -f blk_chr blk_size snpinfo_$LABEL
# this loop goes trough all the regions of the $BLOCKS file and for each of them calculates ld and creates ld matrix/snplist and freq file
echo -e "CHR\tSNP\tBP\tA1\tA2\tMAF" > snpinfo_$LABEL
while read -r -a arr
do
    # calculate ld,frq,snplist and bim files
    plink --bfile $PANEL_DIR --freq --make-just-bim --snps-only --write-snplist --keep-allele-order --chr ${arr[1]} --from-bp ${arr[2]} --to-bp  ${arr[3]}  --r square --out ./ldblk/ldblk_${arr[0]}"_"$LABEL
   # write down block's chr and #of snps
    echo ${arr[1]} |sed 's/chr//g' >> blk_chr
    # create bim file if missing and 
    touch ./ldblk/ldblk_${arr[0]}"_"$LABEL".bim" && wc -l < ./ldblk/ldblk_${arr[0]}"_"$LABEL".bim" >> blk_size
    # produce snpinfo data from bim and frq files and append to final file
    paste <(cut -f 1,2,4,5,6 ./ldblk/ldblk_${arr[0]}"_"$LABEL".bim") <( sed -E 1d ./ldblk/ldblk_${arr[0]}"_"$LABEL".frq"  |awk '{$1=$1;print}' | cut -d " " -f5) >> snpinfo_$LABEL

done < <( cat $BLOCKS |  sed -E 1d |   awk '{print NR,$0}' )

python3 $cwd/write_ldblk.py $DIR $LABEL
mv snpinfo_$LABEL "./ldblk_"$LABEL"_chr/snpinfo_"$LABEL"_hm3"


