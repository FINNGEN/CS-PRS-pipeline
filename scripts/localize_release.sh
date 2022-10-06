#!/bin/bash


CROMWELL_NAME="fg-cromwell_fresh"
CROMWELL_ID=$1
OUT_PATH=$2

echo $CROMWELL_ID $OUT_PATH

mkdir -p $OUT_PATH/release/data/
mkdir -p $OUT_PATH/release/documentation/
mkdir -p $OUT_PATH/munged_rsid/

# SCORES & LOGS
echo "LOGS"
wc -l <  <(gsutil ls gs://$CROMWELL_NAME/prs_cs/$CROMWELL_ID/call-scores/**/finngen*log)
gsutil -m cp  gs://$CROMWELL_NAME/prs_cs/$CROMWELL_ID/call-scores/**/finngen*log $OUT_PATH/release/documentation/ > /dev/null 2>&1

echo "SCORES"
wc -l <  <( gsutil ls gs://$CROMWELL_NAME/prs_cs/$CROMWELL_ID/call-scores/**/finngen*sscore)
gsutil -m cp         gs://$CROMWELL_NAME/prs_cs/$CROMWELL_ID/call-scores/**/finngen*sscore $OUT_PATH/release/data/ > /dev/null 2>&1

# WEIGHTS
echo "WEIGHTS"
wc -l < <(gsutil ls  gs://$CROMWELL_NAME/prs_cs/$CROMWELL_ID/call-weights/**/finngen*weights.txt)
gsutil -m cp gs://$CROMWELL_NAME/prs_cs/$CROMWELL_ID/call-weights/**/finngen*weights.txt $OUT_PATH/release/data/ > /dev/null 2>&1

# MUNGED
echo "MUNGED SUMSTATS"
wc -l < <(gsutil ls gs://$CROMWELL_NAME/prs_cs/$CROMWELL_ID/call-weights/**/munge/finngen*rsid)
gsutil -m cp gs://$CROMWELL_NAME/prs_cs/$CROMWELL_ID/call-weights/**/munge/finngen*rsid $OUT_PATH/munged_rsid/ > /dev/null 2>&1




