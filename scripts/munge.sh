#!/bin/bash
DATA_PATH="/mnt/disks/r8/prs/test"
mkdir -p $DATA_PATH
echo $DATA_PATH


tmpfile=$(mktemp /tmp/abc-script.XXXXXX)
cat ../data/PRS_data.txt | sed -E 1d | grep -wf <(ls $DATA_PATH"/sumstats/") | head  -n 1 > $tmpfile


while read line;\
 do arr=($line) && echo ${arr[0]} && \
 python3 ./munge.py --ss $DATA_PATH/sumstats/${arr[0]}  \
  --effect_type ${arr[8]} --variant ${arr[9]} --chrom ${arr[10]} --pos ${arr[11]} --ref ${arr[12]} --alt ${arr[13]} --effect ${arr[14]} --pval ${arr[15]}  \
  -o $DATA_PATH/ \
  --rsid-map /mnt/disks/r8/data/variant_mapping/finngen.rsid.map.tsv \
  --chrompos-map /mnt/disks/r8/data/variant_mapping/finngen.variants.tsv \
  --chainfile /mnt/disks/r8/data/lift/hg19ToHg38.over.chain.gz \
  --prefix finngen_R8 \
  |& tee -a $DATA_PATH/munge.log  ; done < $tmpfile
