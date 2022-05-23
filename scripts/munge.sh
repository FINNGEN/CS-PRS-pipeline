#!/bin/bash
DATA_PATH="/mnt/disks/r9/prs/test/"
echo $DATA_PATH

PRS_DATA=$1

tmpfile=$(mktemp /tmp/abc-script.XXXXXX)
cat $PRS_DATA | sed -E 1d | grep -wf <(ls $DATA_PATH"/sumstats/")  > $tmpfile

echo "MATCHING FILES"
cat $tmpfile | cut -f1 
echo "MUNGING"
while read line;\
 do arr=($line) && echo ${arr[0]} && \
 python3 ./munge.py --ss $DATA_PATH/sumstats/${arr[0]}  \
  --effect_type ${arr[8]} --variant ${arr[9]} --chrom ${arr[10]} --pos ${arr[11]} --ref ${arr[12]} --alt ${arr[13]} --effect ${arr[14]} --pval ${arr[15]}  \
  -o $DATA_PATH/ \
  --rsid-map /mnt/disks/r9/data/variant_mapping/finngen.rsid.map.tsv \
  --chrompos-map /mnt/disks/r9/data/variant_mapping/finngen.variants.tsv \
  --chainfile /mnt/disks/r9/data/lift/hg19ToHg38.over.chain.gz \
  --prefix finngen_R9  \
  |& tee -a $DATA_PATH/munge.log  ; done < $tmpfile
