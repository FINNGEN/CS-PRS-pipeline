version 1.0

workflow prs_cs {
  input {
    File gwas_meta
    String gwas_data_path
    String prefix

    String bed_file
    File bim_file
    File fam_file
    String ref_dir
    File regions

    # test mode: only chroms 20/21, only the first 2 studies in gwas_meta
    Boolean test = false
  }

  Array[Int] chrom_list = if test then [20, 21] else [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]

  # FinnGen's target build is hg38 (confirmed via tests/test.sh) -- hg38 studies need no liftover.
  Map[String, File] build_chains = {
    "hg18": "gs://finngen-production-library-green/prs/chains/hg18ToHg38.over.chain.gz",
    "hg19": "gs://finngen-production-library-green/prs/chains/hg19ToHg38.over.chain.gz",
    "hg37": "gs://finngen-production-library-green/prs/chains/hg19ToHg38.over.chain.gz",
    "hg38": "gs://finngen-production-library-green/prs/chains/empty_file.txt"
  }

  File rsid_map = "gs://finngen-production-library-green/prs/rsid_mapping/finngen.rsid.map.tsv.pickle.chrompos"

  String general_docker = "eu.gcr.io/finngen-sandbox-v3-containers/bioinformatics:0.7"
  String prs_docker = "eu.gcr.io/finngen-sandbox-v3-containers/cs-prs:r14-se.2"

  call build_dedup_panel {
    input:
    bed_file = bed_file,
    bim_file = bim_file,
    fam_file = fam_file,
    docker = general_docker,
  }

  call validate_inputs {
    input:
    gwas_meta = gwas_meta,
    prefix = prefix,
    test = test,
    docker = general_docker,
  }

  Array[Array[String]] gwas_traits = read_tsv(validate_inputs.sstats)

  scatter (gwas in gwas_traits) {
    String file_name        = gwas[0]
    String pheno             = gwas[1]
    String n_total           = gwas[2]
    String finngen_phenocode = gwas[3]
    String effect_type       = gwas[4]
    String variant           = gwas[5]
    String chrom_col         = gwas[6]
    String pos_col           = gwas[7]
    String effect_allele     = gwas[8]
    String other_allele      = gwas[9]
    String effect_col        = gwas[10]
    String statistic         = gwas[11]
    String statistic_type    = gwas[12]
    String build             = gwas[13]

    call munge {
      input:
      gwas_data_path = gwas_data_path,
      file_name = file_name,
      prefix = prefix,
      effect_type = effect_type,
      variant = variant,
      chrom = chrom_col,
      pos = pos_col,
      ref = effect_allele,
      alt = other_allele,
      effect = effect_col,
      statistic = statistic,
      statistic_type = statistic_type,
      rsid_map = rsid_map,
      chainfile = build_chains[build],
      ref_dir = ref_dir,
      docker = prs_docker,
    }

  }

  call munge_summary {
    input:
    snp_summaries = munge.snp_summary,
    docker = general_docker,
  }

  # single flat scatter over (study, chrom) instead of a nested scatter -- this Cromwell setup
  # doesn't track nested scatters properly (shows up as an untracked synthetic "ScatterAt..."
  # node). cross(study_pairs, chrom_list) builds the full (study, chrom) index space up front.
  Array[Pair[File, String]] study_pairs = zip(munge.munged_file_hm3, n_total)
  Array[Pair[Pair[File, String], Int]] weight_jobs = cross(study_pairs, chrom_list)

  scatter (job in weight_jobs) {
    call weights {
      input:
      munged_gwas = job.left.left,
      N = job.left.right,
      bim_file = bim_file,
      ref_dir = ref_dir,
      chrom = job.right,
      docker = prs_docker,
    }
  }

  scatter (i in range(length(gwas_traits))) {
    call weights_gather {
      input:
      root_name = basename(munge.munged_file_hm3[i], ".munged.hm3.gz"),
      all_weights = weights.weights,
      all_logs = weights.log,
      docker = general_docker,
    }

    call scores {
      input:
      weights = weights_gather.weights,
      bed_file = build_dedup_panel.bed,
      bim_file = build_dedup_panel.bim,
      fam_file = build_dedup_panel.fam,
      freq_file = build_dedup_panel.freq,
      root_name = basename(munge.munged_file_hm3[i], ".munged.hm3.gz"),
      regions = regions,
      pheno = finngen_phenocode[i],
      docker = general_docker,
    }
  }

  output {
    Array[File] weights_out = weights_gather.weights
    Array[File] weights_logs = weights_gather.log
    Array[File] munged_files = munge.munged_file
    Array[File] munged_files_hm3 = munge.munged_file_hm3
    File munge_summary_table = munge_summary.summary
    Array[File] scores_out = flatten(scores.scores)
    Array[File] scores_logs = flatten(scores.log)
  }
}



task validate_inputs {
  input {
    File gwas_meta
    String prefix
    Boolean test
    String docker
  }

  # cols 1,2,3,8-18 of PRS_data.txt; n_cases/n_ctrls/publication/ancestry dropped (never
  # interpolated into a command, would fail the charset check for no reason).
  command <<<
  set -euo pipefail

  cut -f 1,2,3,8,9,10,11,12,13,14,15,16,17,18 <(sed -E 1d ~{gwas_meta}) ~{true="| head -n2" false="" test} > sumstats.txt

  # injection-safe charset: alnum/_/./:/#/+/- only (covers every real value in PRS_data.txt,
  # rejects quotes/$/backticks/;/|/&/whitespace).
  SAFE='^[a-zA-Z0-9_.:#+-]*$'

  : > bad.txt
  grep -vP "$SAFE" <(tr '\t' '\n' < sumstats.txt) >> bad.txt || true
  grep -vP "$SAFE" <(printf '%s\n' "~{prefix}") >> bad.txt || true

  if [[ -s bad.txt ]]; then
      echo "Irregular/unsafe value(s) found in gwas_meta or prefix:" >&2
      cat bad.txt >&2
      exit 1
  fi
  >>>

  runtime {
    docker: "~{docker}"
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    zones: "europe-west1-b"
    preemptible: 1
  }

  output {
    File sstats = "sumstats.txt"
  }
}


task munge {
  input {
    String gwas_data_path
    String file_name
    String prefix

    String effect_type
    String variant
    String chrom
    String pos
    String ref
    String alt
    String effect
    String statistic
    String statistic_type

    File rsid_map
    File chainfile
    String ref_dir

    String docker
    Int disk_factor = 4
  }

  File ss = gwas_data_path + file_name
  String out_root = prefix + "_" + sub(file_name, ".gz", ".munged.gz")
  String out_root_hm3 = sub(out_root, ".munged.gz$", ".munged.hm3.gz")
  String stat_col = if statistic_type == "SE" then "SE" else "P"
  # tolerate a trailing slash in ref_dir, same as weights task
  String clean_ref_dir = sub(ref_dir, "/$", "")
  File snpinfo_file = "~{clean_ref_dir}/snpinfo_1kg_hm3"
  Int disk_size = ceil(size(chainfile, "GB")) + ceil(size(rsid_map, "GB")) + ceil(size(ss, "GB")) * disk_factor + 10

  command <<<
  set -euo pipefail
  python3 /scripts/munge.py -o . --ss ~{ss} \
      --effect_type "~{effect_type}" --variant "~{variant}" --chrom "~{chrom}" --pos "~{pos}" \
      --ref "~{ref}" --alt "~{alt}" --effect "~{effect}" \
      --statistic "~{statistic}" --statistic-type "~{statistic_type}" \
      --prefix "~{prefix}" --rsid-map ~{rsid_map} --chainfile ~{chainfile}

  mkdir -p rejected_variants_out
  cp -r tmp_parse/rejected_variants/. rejected_variants_out/ 2>/dev/null || true
  tar -czf rejected_variants.tar.gz -C rejected_variants_out .

  # PRScs only ever uses the HM3 SNPs in snpinfo_file, so pre-filter to that set for weights.
  awk 'NR>1{print $2}' ~{snpinfo_file} | LC_ALL=C sort -k1,1 > hm3_snps.sorted.txt
  zcat ~{out_root} | tail -n +2 | LC_ALL=C sort -k1,1 -t$'\t' > munged.sorted.txt
  LC_ALL=C join -t $'\t' -1 1 -2 1 munged.sorted.txt hm3_snps.sorted.txt > munged.hm3.txt
  { printf 'SNP\tA1\tA2\tBETA\t%s\n' "~{stat_col}"; cat munged.hm3.txt; } | gzip > ~{out_root_hm3}

  N_SNPS=$(( $(zcat ~{ss} | wc -l) - 1 ))
  MUNGED_SNPS=$(( $(zcat ~{out_root} | wc -l) - 1 ))
  HM3_SNPS=$(( $(zcat ~{out_root_hm3} | wc -l) - 1 ))
  echo "input SNPs: $N_SNPS"
  echo "munged SNPs: $MUNGED_SNPS"
  echo "HM3 SNPs: $HM3_SNPS"
  printf '%s\t%s\t%s\t%s\n' "~{file_name}" "$N_SNPS" "$MUNGED_SNPS" "$HM3_SNPS" > snp_summary.txt
  >>>

  output {
    File munged_file = "~{out_root}"
    File munged_file_hm3 = "~{out_root_hm3}"
    File rejected_variants_tar = "rejected_variants.tar.gz"
    File snp_summary = "snp_summary.txt"
  }

  runtime {
    docker: "~{docker}"
    cpu: 4
    memory: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    zones: "europe-west1-b"
    preemptible: 1
  }
}


task munge_summary {
  input {
    Array[File] snp_summaries
    String docker
  }

  command <<<
  set -euo pipefail
  {
    printf 'file_name\tn_snps\tmunged_snps\tpct_munged_snps\thm3_snps\tpct_hm3_snps\n'
    for f in ~{sep=" " snp_summaries}; do
        awk -F'\t' '{
            pct_munged = ($2>0) ? 100*$3/$2 : 0
            pct_hm3 = ($2>0) ? 100*$4/$2 : 0
            printf "%s\t%s\t%s\t%.2f\t%s\t%.2f\n", $1, $2, $3, pct_munged, $4, pct_hm3
        }' "$f"
    done
  } > munge_summary.txt
  >>>

  output {
    File summary = "munge_summary.txt"
  }

  runtime {
    docker: "~{docker}"
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    zones: "europe-west1-b"
    preemptible: 1
  }
}


task weights {
  input {
    File munged_gwas
    File bim_file
    String ref_dir
    String N
    Int chrom
    String docker
  }

  String root_name = basename(munged_gwas, ".munged.hm3.gz")
  # tolerate a trailing slash in ref_dir (e.g. ".../ldblk_1kg_eur/") -- left as-is it would build
  # a double-slash gs:// path, which is not a valid object key and would fail to localize
  String clean_ref_dir = sub(ref_dir, "/$", "")
  # only this chromosome's LD block + the (always-needed) snpinfo file are localized -- not the
  # other 21 chromosomes' multi-GB hdf5 files, which PRScs never opens for a single-chrom run
  File ldblk_file = "~{clean_ref_dir}/ldblk_1kg_chr~{chrom}.hdf5"
  File snpinfo_file = "~{clean_ref_dir}/snpinfo_1kg_hm3"
  Int disk_size = ceil(size(munged_gwas, "GB")) * 2 + 10

  command <<<
  set -euo pipefail

  BIM_PREFIX="~{sub(bim_file, '.bim$', '')}"

  # reference snpinfo_file to force its localization -- ldblk_file is already referenced below
  : ~{snpinfo_file}

  SUM_STATS="~{munged_gwas}"
  if [[ "$SUM_STATS" == *.gz ]]; then
      gunzip -k "$SUM_STATS"
      SUM_STATS="${SUM_STATS%.gz}"
  fi

  python3 -u /PRScs/PRScs.py \
      --ref_dir "$(dirname "~{ldblk_file}")" --bim_prefix "$BIM_PREFIX" --sst_file "$SUM_STATS" \
      --n_gwas ~{N} --out_dir ~{root_name} --chrom ~{chrom} | tee ~{root_name}.chr~{chrom}.weights.log

  mv "~{root_name}_pst_eff_a1_b0.5_phiauto_chr~{chrom}.txt" ~{root_name}.chr~{chrom}.weights.txt
  >>>

  output {
    File weights = "~{root_name}.chr~{chrom}.weights.txt"
    File log = "~{root_name}.chr~{chrom}.weights.log"
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


task weights_gather {
  input {
    Array[File] all_weights
    Array[File] all_logs
    String root_name
    String docker
  }

  # all_weights/all_logs is the FULL flat (study, chrom) array, not just this study's slice --
  # each weights.txt/log is named "<root_name>.chr<N>.weights.*", so filtering by root_name
  # prefix here picks out just this study's shards regardless of array order.
  command <<<
  set -euo pipefail

  weights_match=()
  for f in ~{sep=" " all_weights}; do
      [[ "$(basename "$f")" == "~{root_name}".chr*.weights.txt ]] && weights_match+=("$f")
  done
  logs_match=()
  for f in ~{sep=" " all_logs}; do
      [[ "$(basename "$f")" == "~{root_name}".chr*.weights.log ]] && logs_match+=("$f")
  done

  cat "${weights_match[@]}" > ~{root_name}.weights.txt
  cat "${logs_match[@]}" > ~{root_name}.weights.log
  >>>

  output {
    File weights = "~{root_name}.weights.txt"
    File log = "~{root_name}.weights.log"
  }

  runtime {
    docker: "~{docker}"
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 20 HDD"
    zones: "europe-west1-b"
    preemptible: 1
  }
}


task scores {
  input {
    File weights
    String bed_file
    String bim_file
    String fam_file
    String freq_file
    String root_name
    File regions
    String pheno
    String docker
  }

  Int disk_size = 20

  command <<<
  set -euo pipefail
  fuse_bed=$(echo "~{bed_file}" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|')
  fuse_bim=$(echo "~{bim_file}" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|')
  fuse_fam=$(echo "~{fam_file}" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|')
  fuse_freq=$(echo "~{freq_file}" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|')

  plink2 --bed "$fuse_bed" --bim "$fuse_bim" --fam "$fuse_fam" --read-freq "$fuse_freq" \
      --score ~{weights} 2 4 6 header center list-variants ignore-dup-ids \
      --out ~{root_name}

  # PRS_regions.txt: one row per phenocode that needs a region-excluded PRS variant too, with a
  # ';'-separated list of chrom_start_end regions (1-based, fully-closed -- matches plink2's
  # '--exclude bed1' convention). Most phenocodes have no row here, so regions.txt stays empty
  # and only the plain sscore above is produced.
  awk -F'\t' -v p="~{pheno}" '$1==p{print $2}' ~{regions} | tr ';' '\n' | tr '_' '\t' | sed '/^$/d' > regions.txt

  if [[ -s regions.txt ]]; then
      # ".no_regions" (dot, not underscore) matches prs_survival.wdl's sumstats task, which
      # derives this exact suffix from gwas_meta's filename via 's/.gz/.no_regions/g'
      plink2 --bed "$fuse_bed" --bim "$fuse_bim" --fam "$fuse_fam" --read-freq "$fuse_freq" \
          --exclude bed1 regions.txt \
          --score ~{weights} 2 4 6 header center list-variants ignore-dup-ids \
          --out ~{root_name}.no_regions
  fi
  >>>

  output {
    Array[File] scores = glob("~{root_name}*.sscore")
    Array[File] log = glob("~{root_name}*.log")
  }

  runtime {
    docker: "~{docker}"
    cpu: 16
    memory: "8 GB"
    disks: "local-disk ~{disk_size} HDD"
    zones: "europe-west1-b"
    preemptible: 1
  }
}

task build_dedup_panel {
  input {
    String bed_file
    File bim_file
    File fam_file
    String docker
  }

  command <<<
  set -euo pipefail
  fuse_bed=$(echo "~{bed_file}" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|')

  plink2 --bed "$fuse_bed" --bim ~{bim_file} --fam ~{fam_file} \
      --rm-dup force-first --make-bed --out dedup

  plink2 --bfile dedup --freq --out dedup
  >>>

  output {
    File bed = "dedup.bed"
    File bim = "dedup.bim"
    File fam = "dedup.fam"
    File freq = "dedup.afreq"
  }

  runtime {
    docker: "~{docker}"
    cpu: 8
    memory: "16 GB"
    disks: "local-disk 500 HDD"
    zones: "europe-west1-b"
    preemptible: 1
  }
}
