#!/usr/bin/env bash
#
# Runs PRS-CS per chromosome and concatenates the results into one weights file.
#
# Replaces cs_wrapper.py now that munge.py/munge_new.py/munge_fast.py emit rsid-keyed,
# PRScs-ready sumstats directly (SNP A1 A2 BETA STAT) -- there is no more rsid<->chrompos
# round-trip to do here, so the whole thing is just: run PRScs per chromosome, concatenate
# its own per-chromosome output.

set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") --ref-file FILE --bim-file FILE --sum-stats FILE --N INT --out DIR
                         [--kwargs "extra PRScs args"] [--chrom "1 2 3 ..."] [--prefix STR]
                         [--force] [--test]

  --ref-file    Path to the PRScs LD reference panel snpinfo file (e.g. snpinfo_1kg_hm3)
  --bim-file    Path to the target plink bim file (rsid-keyed)
  --sum-stats   Path to the rsid-keyed, PRScs-ready sumstats (SNP A1 A2 BETA STAT)
  --N           GWAS sample size
  --out, -o     Output directory
  --kwargs      Extra arguments passed through to PRScs.py, as one quoted string
  --chrom       Space-separated chromosome list to run (default: 1..22)
  --prefix      String to prepend to output filenames
  --force       Re-run chromosomes even if their weight file already exists
  --test        Run PRScs with --n_iter=100 for a quick smoke test
  --venv-dir    Where to create/reuse the local numpy<2.4 venv PRScs runs under
                (default: /mnt/disks/data/prs/venv). Must be on local disk, not a network/synced
                mount (e.g. Dropbox) -- venvs create many small files/symlinks in quick
                succession, which FUSE-backed network filesystems handle unreliably.

Chromosomes run one at a time, in order, each in its own PRScs process. This is deliberately
sequential (no concurrency): PRScs never parallelizes across chromosomes on its own, and running
one at a time means a failure on any one chromosome (e.g. an OOM kill) never loses previously
completed ones -- re-running just picks up where it left off.
EOF
}

KWARGS=""
PREFIX=""
FORCE=0
TEST=0
CHROM=""
VENV_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref-file) REF_FILE="$2"; shift 2;;
        --bim-file) BIM_FILE="$2"; shift 2;;
        --sum-stats) SUM_STATS="$2"; shift 2;;
        --N) N="$2"; shift 2;;
        --out|-o) OUT="$2"; shift 2;;
        --kwargs) KWARGS="$2"; shift 2;;
        --prefix) PREFIX="$2"; shift 2;;
        --force) FORCE=1; shift;;
        --test) TEST=1; shift;;
        --chrom) CHROM="$2"; shift 2;;
        --venv-dir) VENV_DIR="$2"; shift 2;;
        -h|--help) usage; exit 0;;
        *) echo "Unknown argument: $1" >&2; usage; exit 1;;
    esac
done

: "${REF_FILE:?--ref-file is required}"
: "${BIM_FILE:?--bim-file is required}"
: "${SUM_STATS:?--sum-stats is required}"
: "${N:?--N is required}"
: "${OUT:?--out is required}"

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PRSCS="${SCRIPT_DIR}/../PRScs/PRScs.py"

# ---- local-only: numpy>=2.4 breaks PRScs (see docker/requirements.txt for why) -- create/reuse
# a venv pinned from that same requirements file so PRScs always runs against a numpy that works.
# Deliberately kept ABOVE the CORE block: it has no WDL equivalent and should be deleted from
# here once the docker image itself pins numpy<2.4 -- inside docker there's no venv, PRSCS_PYTHON
# just becomes plain "python3" (the image's own interpreter, already on the right numpy).
# Defaults to a fixed local-disk path rather than next to this script: this script lives under
# ~/Dropbox/..., a FUSE network mount, and venv creation (many rapid small-file/symlink writes)
# is unreliable there. Override with --venv-dir.
PRSCS_VENV="${VENV_DIR:-/mnt/disks/data/prs/venv}"
if [[ ! -x "${PRSCS_VENV}/bin/python3" ]]; then
    echo "creating PRScs venv at ${PRSCS_VENV} (numpy<2.4, per docker/requirements.txt)"
    python3 -m venv "$PRSCS_VENV"
    "${PRSCS_VENV}/bin/pip" install -q -r "${SCRIPT_DIR}/../docker/requirements.txt"
fi
PRSCS_PYTHON="${PRSCS_VENV}/bin/python3"

################################################################################################
# ---- CORE (WDL-portable) BLOCK -------------------------------------------------------------
# Every bash variable the rest of this script needs, spelled out here as one flat assignment
# per line -- this is the literal template for the WDL port: declare each as a task input, then
# inside `command <<< >>>` replace each right-hand side below with its `~{...}` expansion (e.g.
# REF_FILE=~{ref_file}) and paste everything from here down verbatim. Nothing below this line
# reads a variable that isn't assigned in this block.
################################################################################################
REF_FILE="$REF_FILE"
BIM_FILE="$BIM_FILE"
SUM_STATS="$SUM_STATS"
N="$N"
OUT="$OUT"
KWARGS="$KWARGS"
PREFIX="$PREFIX"
CHROM="$CHROM"
FORCE="$FORCE"
TEST="$TEST"
PRSCS="$PRSCS"                 # WDL: absolute path to PRScs.py inside the docker image, e.g. /PRScs/PRScs.py
PRSCS_PYTHON="$PRSCS_PYTHON"   # WDL: just "python3" -- the docker image's numpy is already pinned <2.4
################################################################################################

[[ -n "$PREFIX" ]] && PREFIX="${PREFIX}_"

REF_DIR=$(dirname "$REF_FILE")
BIM_DIR=$(dirname "$BIM_FILE")
BIM_ROOT=$(basename "$BIM_FILE" .bim)
BIM_PREFIX="${BIM_DIR}/${BIM_ROOT}"

# strip a trailing .munged/.munged.gz (or plain .gz) the same way cs_wrapper.py did, so
# output filenames match what the rest of the pipeline (WDL glob patterns) expects
SS_BASENAME=$(basename "$SUM_STATS")
SS_ROOT="${PREFIX}${SS_BASENAME%.munged*}"
SS_ROOT="${SS_ROOT%.gz}"

WEIGHTS_PATH="${OUT}/weights"
LOG_PATH="${OUT}/logs"
mkdir -p "$WEIGHTS_PATH" "$LOG_PATH"

if [[ -z "$CHROM" ]]; then
    CHROM_LIST=($(seq 1 22))
else
    CHROM_LIST=($CHROM)
fi
echo "requested chrom list: ${CHROM_LIST[*]}"

if [[ "$TEST" -eq 1 ]]; then
    KWARGS="${KWARGS} --n_iter=100"
fi

TO_RUN=()
if [[ "$FORCE" -eq 1 ]]; then
    TO_RUN=("${CHROM_LIST[@]}")
else
    for c in "${CHROM_LIST[@]}"; do
        shopt -s nullglob
        existing=("${WEIGHTS_PATH}/${SS_ROOT}"*"chr${c}.txt")
        shopt -u nullglob
        [[ ${#existing[@]} -eq 0 ]] && TO_RUN+=("$c")
    done
fi

if [[ ${#TO_RUN[@]} -eq 0 ]]; then
    echo "All chromosomes ran"
else
    echo "final chrom list: ${TO_RUN[*]}"

    # PRScs' own file reads (parse_ref/parse_bim/parse_sumstats) use plain open(), no gzip
    # awareness -- unpack once into the same folder and reuse the unpacked copy after that
    if [[ "$SUM_STATS" == *.gz ]]; then
        UNZIPPED="${SUM_STATS%.gz}"
        if [[ -f "$UNZIPPED" ]]; then
            echo "WARNING: $SUM_STATS is gzipped; PRScs can't read gzip directly, reusing existing unpacked copy at $UNZIPPED"
        else
            echo "WARNING: $SUM_STATS is gzipped; PRScs can't read gzip directly, unpacking to $UNZIPPED"
            gunzip -k "$SUM_STATS"
        fi
        SUM_STATS="$UNZIPPED"
    fi

    OUT_DIR_PREFIX="${WEIGHTS_PATH}/${SS_ROOT}"

    for c in "${TO_RUN[@]}"; do
        echo "chromosome ${c}"
        LOG_FILE="${LOG_PATH}/${SS_ROOT}.${c}.weights.log"
        # shellcheck disable=SC2086
        "$PRSCS_PYTHON" -u "$PRSCS" --ref_dir "$REF_DIR" --bim_prefix "$BIM_PREFIX" --sst_file "$SUM_STATS" \
            --n_gwas "$N" --out_dir "$OUT_DIR_PREFIX" $KWARGS --chrom "$c" > "$LOG_FILE"
    done

    cat "${LOG_PATH}/${SS_ROOT}".*.weights.log > "${OUT}/${SS_ROOT}.weights.log" 2>/dev/null || true
fi

# ---------- merge: PRScs' own per-chromosome output is already rsid-keyed, no conversion needed ----------
OUT_FILE="${OUT}/${SS_ROOT}.weights.txt"
: > "$OUT_FILE"
for c in "${CHROM_LIST[@]}"; do
    shopt -s nullglob
    files=("${WEIGHTS_PATH}/${SS_ROOT}"*"chr${c}.txt")
    shopt -u nullglob
    for f in "${files[@]}"; do cat "$f" >> "$OUT_FILE"; done
done

echo "wrote ${OUT_FILE}"
