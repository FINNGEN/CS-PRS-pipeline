version 1.0

# Not a sandbox-unmodifiable-pipelines WDL -- no injection-safety/input-datastore constraints
# apply here, unlike prs.wdl/prs_weights.wdl.
#
# Default docker/cpu/memory/preemptible/zones come from prs_survival_options.json's
# default_runtime_attributes (the generic bioinformatics image, 2 cpu, 4 GB) -- tasks that don't
# need anything else just declare `disks`. Only survival/correlate_pheno/merge_figs need the
# PRS-specific docker (corr.py/survival_analysis.py + poppler-utils), passed in via the inputs
# json as the `docker` workflow input.

workflow prs_analysis {
  input {
    File prs_metadata
    File regions
    File phenos_file
    File pheno_list_file
    File age_onset
    File score_file_list
    String prefix
    String docker
    String covars = "SEX_IMPUTED,AGE_AT_DEATH_OR_END_OF_FOLLOWUP,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
    Int cpus = 8
    Boolean test = false
  }

  # one gs:// score path per line, matching the read_lines(File) pattern used for ref_dir_list
  # elsewhere in this repo -- kept as a plain path list, not Array[File], so listing them doesn't
  # force Cromwell to localize every score file (matching in sumstats only needs the filenames)
  Array[String] score_files = read_lines(score_file_list)

  call prepare_pheno_matrix {
    input:
    phenos_file = phenos_file,
    pheno_list_file = pheno_list_file,
    covars = covars,
    docker = docker,
  }

  call sumstats {
    input:
    gwas_meta = prs_metadata,
    regions = regions,
    score_files = score_files,
    test = test,
  }
  Array[Array[String]] prs_data = read_tsv(sumstats.sstats)

  scatter (data in prs_data) {
    Int score_idx = data[3]
    String score_file = score_files[score_idx]

    String survival_pheno = if data[2] != "NA" then data[2] else "DEATH"
    call survival {
      input:
      score_file = score_file,
      study = data[0],
      pheno = survival_pheno,
      age_onset = age_onset,
      docker = docker,
    }

    String corr_pheno = if data[2] != "NA" then data[2] else data[1]
    call correlate_pheno {
      input:
      score_file = score_file,
      pheno = corr_pheno,
      phenos_file = prepare_pheno_matrix.pheno_parquet,
      pheno_list_file = pheno_list_file,
      covars = covars,
      cpus = cpus,
      docker = docker,
    }
  }

  call merge_figs {
    input:
    risk_figs = survival.risk_fig,
    survival_figs = survival.survival_fig,
    onset_figs = survival.onset_fig,
    AUC_figs = survival.auc_fig,
    AUC_logs = survival.auc_log,
    prefix = prefix,
    docker = docker,
  }

  call sort_pheno {
    input:
    prefix = prefix,
    corr_files = correlate_pheno.corr_file,
    log_files = correlate_pheno.log,
  }

  output {
    File survival_fig = merge_figs.survival_fig
    File onset_fig = merge_figs.onset_fig
    File risk_fig = merge_figs.risk_fig
    File auc_fig = merge_figs.AUC_fig
    File auc_log = merge_figs.AUC_log
    File sorted_pvals = sort_pheno.sorted_pvals
    File corr_logs = sort_pheno.corr_logs
  }
}


task prepare_pheno_matrix {
  input {
    File phenos_file
    File pheno_list_file
    String covars
    String docker
  }

  Int disk_size = ceil(size(phenos_file, "GB")) * 2 + 10

  command <<<
  set -euo pipefail
  H=$(zcat -f ~{phenos_file} | head -1) || true
  IFS=$'\t' read -r -a COLS <<< "$H"
  ID="IID"; for c in "${COLS[@]}"; do [[ "$c" == "FINNGENID" ]] && ID="FINNGENID"; done

  mapfile -t PHENOS < <(comm -12 <(printf '%s\n' "${COLS[@]}" | LC_ALL=C sort -u) <(LC_ALL=C sort -u ~{pheno_list_file}))
  echo "${#PHENOS[@]} phenotypes shared" >&2

  PSEL=""; for p in "${PHENOS[@]}"; do PSEL+=", COALESCE(TRY_CAST(\"$p\" AS TINYINT), -1) AS \"$p\""; done
  CSEL=""; IFS=',' read -r -a COVS <<< "~{covars}"; for c in "${COVS[@]}"; do CSEL+=", TRY_CAST(\"$c\" AS DOUBLE) AS \"$c\""; done

  # too many phenotype columns to pass as a single duckdb -c argument (hits the OS argv-length
  # limit at real R14 scale, ~2800 columns) -- write the query to a file instead
  cat > q.sql <<SQL
  COPY (SELECT "$ID" AS FINNGENID ${PSEL} ${CSEL}
        FROM read_csv('~{phenos_file}', delim='\t', header=true, all_varchar=true))
  TO 'pheno_matrix.parquet' (FORMAT PARQUET);
  SQL
  duckdb -c ".read q.sql"
  >>>

  output {
    File pheno_parquet = "pheno_matrix.parquet"
  }

  runtime {
    docker: "~{docker}"
    cpu: 4
    memory: "8 GB"
    disks: "local-disk ~{disk_size} HDD"
    zones: "europe-west1-b"
    preemptible: 1
  }
}


task sumstats {
  input {
    File gwas_meta
    File regions
    Array[String] score_files
    Boolean test
  }

  command <<<
  set -euo pipefail
  INPUT=$(sed -E 1d ~{gwas_meta})
  ~{true='INPUT=$(echo "$INPUT" | tail -n 2)' false='' test}

  echo "$INPUT" | cut -f 1,2,8 | sed -e 's/.gz//g' > sumstats_base.txt
  echo "$INPUT" | cut -f 1,2,8 | sed -e 's/.gz/.no_regions/g' | grep -wf <(cut -f1 ~{regions}) >> sumstats_base.txt

  : > score_index.txt
  i=0
  for f in ~{sep=" " score_files}; do
      echo -e "$(basename "$f")\t$i" >> score_index.txt
      i=$((i+1))
  done

  : > sumstats.txt
  while IFS=$'\t' read -r filename pheno phenocode; do
      # anchored suffix match, not substring: a plain filename (e.g. "AD_sumstats_Jansenetal.txt")
      # is a literal substring of its own ".no_regions" score file's basename
      # ("finngen_AD_sumstats_Jansenetal.txt.no_regions.sscore"), so a plain grep -F here would
      # match both and silently bind them to the same file via head -1. Requiring the score
      # basename to end in exactly "<filename>.sscore" disambiguates the two.
      idx=""
      while IFS=$'\t' read -r basename bidx; do
          case "$basename" in
              *"${filename}.sscore") idx="$bidx"; break ;;
          esac
      done < score_index.txt
      if [[ -z "$idx" ]]; then
          echo "ERROR: no matching score file found for study '$filename' in score_files" >&2
          exit 1
      fi
      echo -e "${filename}\t${pheno}\t${phenocode}\t${idx}" >> sumstats.txt
  done < sumstats_base.txt
  >>>

  output {
    File sstats = "sumstats.txt"
  }

  runtime {
    disks: "local-disk 5 HDD"
  }
}


task survival {
  input {
    String pheno
    String study
    File score_file
    File age_onset
    String docker
  }

  Int disk_size = ceil(size(age_onset, "GB") + size(score_file, "GB")) * 2 + 2

  command <<<
  set -euo pipefail
  mkdir -p survival
  python3 /scripts/survival_analysis.py \
      --scores ~{score_file} --age_file ~{age_onset} --pheno ~{pheno} --tag ~{study} --out survival/
  >>>

  output {
    File survival_fig = "survival/~{pheno}_~{study}_survival.pdf"
    File onset_fig = "survival/~{pheno}_~{study}_age_onset.pdf"
    File risk_fig = "survival/~{pheno}_~{study}_risk.pdf"
    File auc_fig = "survival/~{pheno}_~{study}_AUC.pdf"
    File auc_log = "survival/~{pheno}_~{study}_AUC.log"
  }

  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    disks: "local-disk ~{disk_size} HDD"
  }
}


task merge_figs {
  input {
    Array[File] survival_figs
    Array[File] onset_figs
    Array[File] risk_figs
    Array[File] AUC_figs
    Array[File] AUC_logs
    String prefix
    String docker
  }

  command <<<
  set -euo pipefail
  pdfunite ~{sep=" " survival_figs} ~{prefix}_survival.pdf
  pdfunite ~{sep=" " onset_figs} ~{prefix}_onset.pdf
  pdfunite ~{sep=" " AUC_figs} ~{prefix}_AUC.pdf
  pdfunite ~{sep=" " risk_figs} ~{prefix}_risk.pdf
  cat ~{sep=" " AUC_logs} > ~{prefix}_AUC.log
  >>>

  output {
    File survival_fig = "~{prefix}_survival.pdf"
    File onset_fig = "~{prefix}_onset.pdf"
    File risk_fig = "~{prefix}_risk.pdf"
    File AUC_fig = "~{prefix}_AUC.pdf"
    File AUC_log = "~{prefix}_AUC.log"
  }

  runtime {
    docker: "~{docker}"
    disks: "local-disk 5 HDD"
  }
}


task correlate_pheno {
  input {
    File phenos_file
    File pheno_list_file
    File score_file
    String pheno
    String covars
    Int cpus
    String docker
  }

  # corr.py runs one independent OS process per phenotype (GNU parallel -j cpus --memfree), not a
  # shared-memory pool -- per-worker cost varies by phenotype (sample size after dropna,
  # covariate patterns), measured up to ~3.9GB for one worker vs a ~1-1.5GB typical average. Not
  # over-provisioning the task's own memory for that worst case: parallel's --memfree throttles
  # job starts and requeues the youngest job if free memory runs low, and --retries backstops
  # whatever that doesn't catch in time, so an occasional OOM is an expected, self-healing,
  # cheap-to-absorb case rather than something to size the task around. --memfree is derived from
  # the same per-cpu multiplier as mem, not a separately hardcoded value, so they can't drift
  # out of sync if the multiplier changes.
  Int mem_multiplier = 2
  Int mem = cpus * mem_multiplier
  String out_file = basename(score_file, ".sscore") + "_corr.txt"
  String log_file = basename(score_file, ".sscore") + "_corr.log"
  Int disk_size = ceil(size(phenos_file, "GB") + size(score_file, "GB")) * 2 + 2

  command <<<
  set -euo pipefail
  python3 /scripts/corr.py --pheno-file ~{phenos_file} --pheno-list ~{pheno_list_file} \
      --pheno ~{pheno} --scores ~{score_file} --cov ~{covars} --cpus ~{cpus} \
      --memfree ~{mem_multiplier}G
  >>>

  output {
    File log = log_file
    File corr_file = out_file
  }

  runtime {
    docker: "~{docker}"
    cpu: cpus
    memory: "~{mem} GB"
    disks: "local-disk ~{disk_size} HDD"
  }
}


task sort_pheno {
  input {
    Array[File] corr_files
    Array[File] log_files
    String prefix
  }

  String out_file = prefix + "_prs_pheno_corr.tsv"
  String out_log = prefix + "_prs_pheno_corr.log"

  command <<<
  set -euo pipefail
  head -n1 ~{corr_files[0]} | awk '$3="log(pval)"' > tmp.txt

  bodies=()
  i=0
  for f in ~{sep=" " corr_files}; do
      out="body_${i}.tsv"
      tail -n +2 "$f" > "$out"
      bodies+=("$out")
      i=$((i+1))
  done

  sort -m -t $'\t' -k 3,3 -g "${bodies[@]}" \
      | awk -F'\t' 'BEGIN{OFS="\t"} $3!="NA"{if ($3==0) $3="inf"; else $3=-log($3)/log(10); print}' \
      | awk -F'\t' '$3>4' >> tmp.txt

  column -t tmp.txt > ~{out_file}
  cat ~{sep=" " log_files} > ~{out_log}
  >>>

  output {
    File sorted_pvals = out_file
    File corr_logs = out_log
  }

  runtime {
    disks: "local-disk 5 HDD"
  }
}
