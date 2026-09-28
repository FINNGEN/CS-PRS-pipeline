# CS-PRS-pipeline

Calculates PRS weights (and, outside the sandbox, scores) from GWAS summary stats using
[PRScs](https://github.com/getian107/PRScs). Runs entirely in **rsid space**: sumstats, weights
and scoring all key on rsid, with no chrom_pos_ref_alt round-trip or allele-permutation step.

This repo is no longer responsible for FinnGen's full production trait list — real PRS
calculation for users happens via the `sandbox-unmodifiable-pipelines` port of `prs_weights.wdl`.
This repo is for developing and validating the pipeline itself, against a small, deliberately
diverse debug set of studies (`data/PRS_data.txt`), not the full ~140-study production list.

## Pipeline

1. **munge** (`scripts/munge.py`) — parses a raw sumstat, resolves each variant to a FinnGen
   rsid (via a chrompos↔rsid map, with liftover to hg38 for non-hg38 studies), and writes a
   5-column PRScs-ready file (`SNP A1 A2 BETA STAT`).
2. **weights** (`PRScs/PRScs.py`, invoked per chromosome) — MCMC posterior effect-size estimation
   against a reference LD panel.
3. **scores** (`scripts/cs_scores.sh`, plink2 `--score`) — applies weights to genotypes. Not part
   of the sandbox-bound weights-only WDL: a `.sscore` file is individual-level data.

## WDLs (`wdl/`)

- **`prs.wdl`** — full pipeline (munge → weights → scores), scattered over (study, chrom) since
  one study's 22 chromosomes take ~1.5-2 days run sequentially. Dev/test use only.
- **`prs_weights.wdl`** — weights-only, same scatter design, no scores task. This is what gets
  ported into `sandbox-unmodifiable-pipelines`.

Both take `data/PRS_data.txt` (18 columns: filename, phenotype metadata, and per-study column
names/types for `munge.py`, including `statistic`/`statistic_type` for SE-vs-pval handling).

## Docker

`docker/build_docker.py --image cs-prs --version X --registry {refinery,sandbox}` builds/pushes
to either registry. `docker/Dockerfile` itself is still the pre-migration version and needs
rewriting (numpy pin, `PRScs`/`scripts` layout) before either WDL can actually run against it.

## Legacy

`wdl/finngen_weights.wdl` and `wdl/sandbox/` predate the rsid migration (a separate FG-specific
LD panel workflow) and haven't been touched by it.
