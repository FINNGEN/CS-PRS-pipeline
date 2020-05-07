#!/bin/bash

while read line; do arr=($line) &&   python3 ./munge.py --ss /mnt/disks/r5/PRS/gwas_raw/${arr[0]} -o /mnt/disks/r5/PRS/docker_test/munge/ --effect_type ${arr[8]} --variant ${arr[9]}  --chrom ${arr[10]} --pos ${arr[11]}  --ref ${arr[12]} --alt ${arr[13]} --effect ${arr[14]} --pval ${arr[15]}  --rsid-map /mnt/disks/r5/PRS/variant_mapping/finngen.rsid.map.tsv --chrompos-map /mnt/disks/r5/PRS/variant_mapping/finngen.variants.tsv --chainfile /mnt/disks/r5/Data/lift/hg19ToHg38.over.chain.gz --prefix finngen_R5  ; done < <(cat ../data/PRS_data.txt | sed -E 1d | grep MTAG_CP)



python3 ./cs_wrapper.py --bim-file /mnt/disks/r5/plink/R5.hm3.rsid.bim --ref-file /mnt/disks/r5/PRS/ldblk_1kg_eur/snpinfo_1kg_hm3 --map /mnt/disks/r5/PRS/variant_mapping/finngen.rsid.map.tsv --out /mnt/disks/r5/PRS/docker_test/ --N 257828 --sum-stats /mnt/disks/r5/PRS/docker_test/munge/finngen_R5_MTAG_CP.to10K.txt.munged.gz --parallel 2


python3 ./cs_scores.py --weight /mnt/disks/r5/PRS/docker_test/finngen_R5_MTAG_CP.to10K.txt.weights.txt --bed /mnt/disks/r5/PRS/plink/R5.hm3.test.bed --out /mnt/disks/r5/PRS/docker_test/ 
