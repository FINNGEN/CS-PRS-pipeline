#!/bin/bash

while read line; do arr=($line) &&   python3 ./munge.py --ss /mnt/disks/r5/PRS/gwas_raw/${arr[0]} -o /mnt/disks/r5/PRS/munge/ --effect_type ${arr[8]} --variant ${arr[9]} --chrom ${arr[10]} --pos ${arr[11]} --ref ${arr[12]} --alt ${arr[13]} --effect ${arr[14]} --pval ${arr[15]}  --rsid-map /mnt/disks/r5/PRS/variant_mapping/finngen.rsid.map.tsv --chrompos-map /mnt/disks/r5/PRS/variant_mapping/finngen.variants.tsv --chainfile /mnt/disks/r5/Data/lift/hg19ToHg38.over.chain.gz  |& tee -a /mnt/disks/r5/PRS/munge/munge.log; done < <(cat ../data/PRS_data.txt| sed -E 1d  )
