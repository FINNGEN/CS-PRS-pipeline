#!/usr/bin/env bash
##
##     @script.name [option] ARGUMENTS...
##
## Options:
##     -h, --help              All client scripts have this, it can be omitted.
##     --var=VALUE        index for variant column chr:pos:ref:alt. This or next four must be specified
##     --chr=VALUE        Columnn indexes
##     --pos=VALUE
##     --ref=VALUE
##     --alt=VALUE
##     --sep=VALUE       column separator in data. if you want space give 'space' as easyoptions does not accept space as cmdline parameter



export LCTYPE=C
export LANG=C
my_dir="$(dirname "$0")"
source "${my_dir}/easyoptions.sh" || exit

inputfile=$1
liftedfile=$2

cat_cmd="cat"
file $inputfile | grep gzip > /dev/null 2>&1
if [[ $? -eq 0 ]];
then
    cat_cmd='zcat'
fi

cols=$( $cat_cmd < "$inputfile" | head -n 1 | awk 'BEGIN {FS="\t"}{ print NF}')


if [[ ! -n "$sep"  ]]
then
    ## tried to add optional sep but could not figure out how to get proper tab character to variable to be used below
    ## this does not work
    colsep=$'-t\t'
else
    sep=${sep/space/ /}
    cat_cmd=$cat_cmd" |  tr "$sep" '\t'"
fi

if [[ ! -n "$var"  ]]
then
    if [[ ! -n "$chr"  ]] || [[ ! -n "$pos"  ]] || [[ ! -n "$ref"  ]]  || [[ ! -n "$alt"  ]]
    then
        echo "chr pos ref alt columns must be specified if var not given"
        exit 1
    fi
    cols=$((cols+1))
    cat <( $cat_cmd < "$inputfile" | awk -v var=$var 'BEGIN {FS="\t"}  NR==1{ printf "#variant"; for(i=1;i<=NF; i++) if(i!=var) printf "\t"$i; printf "\tlift_chr\tlift_pos\tREF\tALT" ;printf "\n"; exit 0 }'  )  \
        <( join -1 1 -2 3 -t$'\t' -a 1 <(  $cat_cmd < "$inputfile" | awk -v chr=$chr -v pos=$pos -v ref=$ref -v alt=$alt ' BEGIN{OFS="\t"; FS="\t"} NR>1{ print $chr":"$pos":"$ref":"$alt,$0 }' | sort -b -k 1,1 -t $'\t'  ) <( awk 'BEGIN{ OFS="\t"; } { print $1,$2+1,$4}'  $liftedfile | sort -t$'\t' -b -k 3,3 ) \
| awk -v n_exp_col=$((cols+2)) 'BEGIN{ FS="\t"; OFS="\t"} { split($1,a,":"); printf $0;  if (NF<n_exp_col) { printf "\tNA\t0"; for(i=NF+2;i<n_exp_col;i++) {printf "\tNA" } };  printf "\t"a[3]"\t"a[4]"\n"; } ' | sort -t$'\t' -V -k $((cols+1)),$((cols+1))  -k $((cols+2)),$((cols+2)) ) |  bgzip > $(basename $inputfile)".lifted.gz"
else

cat <( $cat_cmd < "$inputfile" | awk -v var=$var 'BEGIN {FS="\t"}  NR==1{ printf "#variant"; for(i=1;i<=NF; i++) if(i!=var) printf "\t"$i; printf "\tlift_chr\tlift_pos\tREF\tALT" ;printf "\n"; exit 0 }'  )  \
    <(join -1 1 -2 3 -t$'\t' -a 1 <(  $cat_cmd < "$inputfile" | awk -v var=$var 'BEGIN {FS="\t"}  NR>1{ printf $var; for(i=1;i<=NF; i++) if(i!=var) printf "\t"$i; printf "\n"; }' | sort -b -k 1,1 -t $'\t'  ) \
        <( awk 'BEGIN{ OFS="\t";}{ print $1,$2+1,$4}'  $liftedfile | sort -t$'\t' -b -k 3,3 )  \
        | awk -v n_exp_col=$((cols+3)) 'BEGIN{ FS="\t"; OFS="\t"} { split($1,a,":"); printf $0;  if (NF<n_exp_col) { printf "\tNA\t0"; for(i=NF+2;i<=exp_col;i++) printf "\tNA" }  printf "\t"a[3]"\t"a[4]"\n" } ' | sort -t$'\t' -V -k $((cols+1)),$((cols+1))  -k $((cols+2)),$((cols+2)) ) |  bgzip > $(basename $inputfile)".lifted.gz"
fi

tabix -s $((cols+1)) -b $((cols+2)) -e $((cols+2)) $(basename $inputfile)".lifted.gz"
