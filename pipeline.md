# PRS pipeline (`prs_sandbox.wdl`) — inputs

## `gwas_meta`
Tab-separated table with a header line and one row per GWAS study. All 18 columns must be present
in this exact order (they are selected by position):

| # | Column | Used | Description |
|---|---|---|---|
| 1 | `filename` | yes | Sumstats file name inside `gwas_data_path` (gzipped). Names every per-study output, so it must be unique. |
| 2 | `pheno` | yes | Free-text study label; used as the correlation phenotype when `finngen_phenocode` is `NA`. |
| 3 | `n_total` | yes | GWAS sample size, passed to PRScs as `--n_gwas`. |
| 4 | `n_cases` | no | Informational. |
| 5 | `n_ctrls` | no | Informational. |
| 6 | `publication` | no | Informational (free text, may contain `;` and URLs). |
| 7 | `ancestry` | no | Informational. |
| 8 | `finngen_phenocode` | yes | Matching FinnGen endpoint (e.g. `G6_ALZHEIMER`), or `NA`. Used as the survival endpoint (`DEATH` if `NA`), as the correlation phenotype, and to look up `regions`. |
| 9 | `effect_type` | yes | `BETA` or `OR` (case-insensitive). OR is converted to log-OR. |
| 10 | `variant` | yes | Name of the variant-ID column in the sumstats (rsid or chrom/pos-like), or `NA` if there is none (then `chrom` and `pos` are required). |
| 11 | `chrom` | yes | Name of the chromosome column, or `NA` to parse it from `variant`. |
| 12 | `pos` | yes | Name of the position column, or `NA` to parse it from `variant`. |
| 13 | `effect_allele` | yes | Name of the effect-allele column. |
| 14 | `other_allele` | yes | Name of the other-allele column. |
| 15 | `effect` | yes | Name of the effect column (beta or OR, matching `effect_type`). |
| 16 | `statistic` | yes | Name of the standard-error or p-value column. |
| 17 | `statistic_type` | yes | `SE` or `PVAL`: whether `statistic` is a standard error or a p-value. |
| 18 | `build` | yes | Genome build of the sumstats: `hg18`, `hg19`, `hg37` or `hg38`. Non-hg38 studies are lifted over to hg38. |

Columns 4–7 are dropped before validation; every other value must match `^[a-zA-Z0-9_.:#+-]*$`
(no whitespace, quotes, `$`, backticks, `;`, `|`, `&` or `/`).

Example row:

```
AD_sumstats_Jansenetal.txt.gz	Alzheimers_Jansen_2019	455258	71880	383378	https://www.nature.com/articles/s41588-018-0311-9	European	G6_ALZHEIMER	BETA	SNP	CHR	BP	A1	A2	BETA	SE	SE	hg19
```

## `gwas_data_path`
Folder containing the sumstats files listed in `gwas_meta`. Must end with `/`, since each study's
file is `gwas_data_path + filename`.

## `prefix`
Short tag prepended to the per-study outputs and used to name the merged figures and the correlation table.

## `bed_file`
PLINK `.bed` of the dataset to score. The `.bim`, `.fam` and `.afreq` must sit next to it with the same
prefix. Variant IDs in the `.bim` can be rsids or `chr{CHR}_{POS}_{REF}_{ALT}`; the notation is detected
automatically.

## `regions`
Tab-separated file with a header, listing genomic regions to exclude for a given FinnGen endpoint:

```
pheno	regions
G6_ALZHEIMER	19_35000000_70000000
```

Regions are `CHROM_START_END`, separated by `;`. Studies whose `finngen_phenocode` is listed are scored
a second time with these regions removed (`.no_regions`).

## `pheno_list_file`
Text file with one FinnGen binary endpoint per line; each score is tested for association with every
listed endpoint. Set to `""` to skip the correlation step.
