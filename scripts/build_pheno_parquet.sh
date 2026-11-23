#!/bin/bash
# Builds a compact Parquet phenotype+covariate matrix from a large gz pheno file, for reuse
# across multiple corr.py runs instead of every run re-decompressing/re-parsing the same file.
#
# Phenotype columns are cast to a plain TINYINT (duckdb's 1-byte int) with NA encoded as -1 --
# real values are 0/1 only, so -1 is a safe, unused sentinel. This avoids ever materializing a
# nullable/masked representation (pandas' nullable Int8 costs double: values + a boolean mask,
# both 1 byte/cell) for a ~2800-column x ~500k-row matrix. corr.py's parquet fast-path converts
# -1 back to NaN on the small per-fit slice, not the whole matrix.
#
# Caller is responsible for only passing binary (0/1) phenotypes in --pheno-list -- a
# quantitative phenotype would silently get rounded/clamped into this TINYINT cast.
#
# Usage:
#   build_pheno_parquet.sh --pheno-file <gz> --pheno-list <txt> [--cov <comma-list>] --out <parquet>
set -euo pipefail

COV="SEX_IMPUTED,AGE_AT_DEATH_OR_END_OF_FOLLOWUP,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"

usage() {
    echo "Usage: $0 --pheno-file <gz> --pheno-list <txt> [--cov <comma-list>] --out <parquet>" >&2
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pheno-file) PHENO_FILE="$2"; shift 2 ;;
        --pheno-list) PHENO_LIST="$2"; shift 2 ;;
        --cov)        COV="$2";        shift 2 ;;
        --out)        OUT="$2";        shift 2 ;;
        *) echo "Unknown argument: $1" >&2; usage ;;
    esac
done

[[ -z "${PHENO_FILE:-}" || -z "${PHENO_LIST:-}" || -z "${OUT:-}" ]] && usage
[[ -f "$PHENO_FILE" ]] || { echo "no such file: $PHENO_FILE" >&2; exit 1; }
[[ -f "$PHENO_LIST" ]] || { echo "no such file: $PHENO_LIST" >&2; exit 1; }
command -v duckdb >/dev/null || { echo "duckdb not found on PATH" >&2; exit 1; }

echo "$(date +%T) reading header of $PHENO_FILE..." >&2
# `|| true`: head closing the pipe after one line sends zcat a SIGPIPE, which under pipefail
# would otherwise abort the whole script even though $HEADER is captured correctly regardless
HEADER=$(zcat -f "$PHENO_FILE" | head -1) || true
IFS=$'\t' read -r -a HEADER_COLS <<< "$HEADER"

# id column: FINNGENID if present, else IID -- same convention as corr.py's return_header check
ID_COL="IID"
for c in "${HEADER_COLS[@]}"; do
    [[ "$c" == "FINNGENID" ]] && ID_COL="FINNGENID"
done

# phenotypes present in BOTH the pheno file header and the requested pheno list -- same
# intersection corr.py's parallel() computes per-shard, done here once instead
mapfile -t SHARED_PHENOS < <(comm -12 \
    <(printf '%s\n' "${HEADER_COLS[@]}" | LC_ALL=C sort -u) \
    <(LC_ALL=C sort -u "$PHENO_LIST"))

echo "$(date +%T) ${#SHARED_PHENOS[@]} phenotypes shared between list and pheno file" >&2
[[ ${#SHARED_PHENOS[@]} -gt 0 ]] || { echo "no phenotypes in common -- check --pheno-list" >&2; exit 1; }

IFS=',' read -r -a COV_COLS <<< "$COV"

# TRY_CAST returns NULL on anything that isn't a valid tinyint (an empty field or a literal "NA"
# both fail to parse and become NULL), so COALESCE(...,-1) catches either missing-value
# convention without needing to know which one the file actually uses.
PHENO_EXPRS=""
for p in "${SHARED_PHENOS[@]}"; do
    PHENO_EXPRS="${PHENO_EXPRS}, COALESCE(TRY_CAST(\"$p\" AS TINYINT), -1) AS \"$p\""
done

COV_EXPRS=""
for c in "${COV_COLS[@]}"; do
    COV_EXPRS="${COV_EXPRS}, TRY_CAST(\"$c\" AS DOUBLE) AS \"$c\""
done

mkdir -p "$(dirname "$OUT")"

echo "$(date +%T) building $OUT..." >&2
# with ~2800 phenotype columns, the generated SELECT is too long to pass as a single -c argument
# (hits the OS's execve() argument-length limit) -- write it to a temp file and have duckdb read
# that instead, which has no such limit
TMPSQL=$(mktemp --suffix=.sql)
trap 'rm -f "$TMPSQL"' EXIT

# all_varchar=true: read every source column as text and do our own explicit casts above,
# instead of letting duckdb's CSV reader sniff types over a sample (which risks a wrong guess
# from a column whose only non-null values appear late in a ~500k-row file)
cat > "$TMPSQL" <<EOF
COPY (
    SELECT "$ID_COL" AS FINNGENID ${PHENO_EXPRS} ${COV_EXPRS}
    FROM read_csv('$PHENO_FILE', delim='\t', header=true, all_varchar=true)
) TO '$OUT' (FORMAT PARQUET);
EOF

duckdb -c ".read $TMPSQL"

echo "$(date +%T) done: $OUT" >&2
