#!/usr/bin/env bash
#
# Computes plink2 PRS scores from one or more PRS-CS weight files.
#
# Replaces cs_scores.py. Same plink2 --score call, same column indices (2 4 6: SNP A1 BETA --
# matches PRScs' own CHR SNP BP A1 A2 BETA output, unchanged since munge/cs_wrapper.sh moved to
# rsid-keyed weights). bim/fam are derived from --bed's own root name by default (so passing just
# --bed keeps behaving exactly like plink2's --bfile); pass --bim (and/or --fam) explicitly to
# override just that one file -- e.g. scoring against the rsid-keyed bim while --bed/--fam stay
# the original FinnGen genotype files.

set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") (--weight FILE | --weight-list FILE) (--bed FILE | --pgen FILE)
                         [--bim FILE] [--fam FILE] [--freq FILE] [--region FILE] [--out DIR]
                         [--memory INT]

  --weight       Path to a single PRS-CS weight file
  --weight-list  Path to a file listing weight file paths, one per line (first tab-separated column)
  --bed          Path to a plink1 bed file
  --pgen         Path to a plink2 pgen file (mutually exclusive with --bed)
  --bim          Override the bim file (default: --bed's own path with .bed -> .bim)
  --fam          Override the fam file (default: --bed's own path with .bed -> .fam)
  --freq         Override the allele-frequency (--read-freq) file (default: --bed/--pgen's own
                 path with .bed/.pgen -> .afreq). Needs to match whatever --bim's variant IDs
                 actually are -- e.g. if --bim is rsid-keyed but the default .afreq next to --bed
                 is still chrompos-keyed, none of its rows will match by ID and plink2 silently
                 falls back to computing frequencies from the loaded sample instead (correct at
                 large sample sizes, but expensive/slow, and skips whatever the .afreq was meant
                 to supply) -- pass a --freq file with matching (e.g. also rsid-keyed) IDs instead.
  --region       Path to a list of regions to additionally score while excluding them
  --out, -o      Output directory (default: .)
  --memory       Memory in MiB to pass to plink2 (default: total system RAM in MiB)

Each weight's score file is skipped if its .sscore already exists (rerun-safe).
EOF
}

WEIGHT=""
WEIGHT_LIST=""
BED=""
PGEN=""
BIM=""
FAM=""
FREQ=""
REGION=""
OUT="."
MEMORY=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --weight) WEIGHT="$2"; shift 2;;
        --weight-list) WEIGHT_LIST="$2"; shift 2;;
        --bed) BED="$2"; shift 2;;
        --pgen) PGEN="$2"; shift 2;;
        --bim) BIM="$2"; shift 2;;
        --fam) FAM="$2"; shift 2;;
        --freq) FREQ="$2"; shift 2;;
        --region) REGION="$2"; shift 2;;
        --out|-o) OUT="$2"; shift 2;;
        --memory) MEMORY="$2"; shift 2;;
        -h|--help) usage; exit 0;;
        *) echo "Unknown argument: $1" >&2; usage; exit 1;;
    esac
done

if [[ -z "$WEIGHT" && -z "$WEIGHT_LIST" ]]; then echo "one of --weight or --weight-list is required" >&2; usage; exit 1; fi
if [[ -n "$WEIGHT" && -n "$WEIGHT_LIST" ]]; then echo "--weight and --weight-list are mutually exclusive" >&2; usage; exit 1; fi
if [[ -z "$BED" && -z "$PGEN" ]]; then echo "one of --bed or --pgen is required" >&2; usage; exit 1; fi
if [[ -n "$BED" && -n "$PGEN" ]]; then echo "--bed and --pgen are mutually exclusive" >&2; usage; exit 1; fi

if [[ -z "$MEMORY" ]]; then
    MEMORY=$(( $(grep -i '^MemTotal:' /proc/meminfo | awk '{print $2}') / 1024 ))
fi

################################################################################################
# ---- CORE (WDL-portable) BLOCK -------------------------------------------------------------
# Every bash variable the rest of this script needs, spelled out here as one flat assignment
# per line -- this is the literal template for the WDL port: declare each as a task input, then
# inside `command <<< >>>` replace each right-hand side below with its `~{...}` expansion (e.g.
# WEIGHT=~{weight}) and paste everything from here down verbatim. Nothing below this line reads
# a variable that isn't assigned in this block.
################################################################################################
WEIGHT="$WEIGHT"
WEIGHT_LIST="$WEIGHT_LIST"
BED="$BED"
PGEN="$PGEN"
BIM="$BIM"
FAM="$FAM"
FREQ="$FREQ"
REGION="$REGION"
OUT="$OUT"
MEMORY="$MEMORY"
################################################################################################

SCORES_PATH="${OUT}/scores"
mkdir -p "$SCORES_PATH"

WEIGHT_FILES=()
if [[ -n "$WEIGHT" ]]; then
    WEIGHT_FILES=("$WEIGHT")
else
    while IFS=$'\t' read -r w _; do
        [[ -n "$w" ]] && WEIGHT_FILES+=("$w")
    done < "$WEIGHT_LIST"
fi

if [[ -n "$PGEN" ]]; then
    PLINK_ROOT="${PGEN%.pgen}"
    PLINK_INPUT_ARGS=(--pfile "$PLINK_ROOT")
else
    BED_ROOT="${BED%.bed}"
    BIM_FILE="${BIM:-${BED_ROOT}.bim}"
    FAM_FILE="${FAM:-${BED_ROOT}.fam}"
    PLINK_ROOT="$BED_ROOT"
    PLINK_INPUT_ARGS=(--bed "$BED" --bim "$BIM_FILE" --fam "$FAM_FILE")
fi

FREQ_ARGS=()
FREQ_FILE="${FREQ:-${PLINK_ROOT}.afreq}"
if [[ -f "$FREQ_FILE" ]]; then
    echo "freq file present: ${FREQ_FILE}"
    FREQ_ARGS=(--read-freq "$FREQ_FILE")
fi

for weight_file in "${WEIGHT_FILES[@]}"; do
    [[ -f "$weight_file" ]] || continue
    basename_w=$(basename "$weight_file")
    root_name="${basename_w%.weights*}"
    score_file="${SCORES_PATH}/${root_name}"

    if [[ ! -f "${score_file}.sscore" ]]; then
        # || true: one weight file's plink2 failure (e.g. no overlapping variants) shouldn't
        # abort the rest of the batch -- matches cs_scores.py's subprocess.call(), which never
        # checked the return code either
        plink2 "${PLINK_INPUT_ARGS[@]}" "${FREQ_ARGS[@]}" --out "$score_file" \
            --score "$weight_file" 2 4 6 header center list-variants ignore-dup-ids \
            --memory "$MEMORY" || echo "WARNING: plink2 failed for ${weight_file}"
    else
        echo "${score_file} already generated"
    fi

    if [[ -n "$REGION" && -s "$REGION" ]]; then
        region_root=$(basename "$REGION")
        region_root="${region_root%.*}"
        no_region_file="${score_file}.no_${region_root}"
        if [[ ! -f "${no_region_file}.sscore" ]]; then
            plink2 "${PLINK_INPUT_ARGS[@]}" "${FREQ_ARGS[@]}" --memory "$MEMORY" \
                --score "$weight_file" 2 4 6 center list-variants \
                --exclude range "$REGION" --out "$no_region_file" || echo "WARNING: plink2 (region exclusion) failed for ${weight_file}"
        else
            echo "${no_region_file} already generated"
        fi
    fi
done
