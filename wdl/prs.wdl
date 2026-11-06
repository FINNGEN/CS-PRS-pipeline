version 1.0

workflow prs_cs {
  input {
    File gwas_meta
    String gwas_data_path
    String prefix

    # target dataset: .bim/.fam/.afreq must sit next to it with the same prefix. Read through FUSE
    # for scoring; only its .bim is localized (rsid_bim task). Variant IDs may be rsid or chrom_pos_ref_alt.
    String bed_file
    String ref_dir
    File regions

    File phenos_file
    File pheno_list_file
    File age_onset
    String covars = "SEX_IMPUTED,AGE_AT_DEATH_OR_END_OF_FOLLOWUP,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
    Int corr_cpus = 8

    String prs_docker
    String survival_docker

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



  call validate_inputs {
    input:
    gwas_meta = gwas_meta,
    prefix = prefix,
    test = test,
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
  }

  # detects the variant ID notation of the input bim and builds its row-aligned companion in the
  # other notation: PRScs needs the rsid bim, scores need weights in the dataset's own notation
  call rsid_bim {
    input:
    bed_file = bed_file,
    rsid_map = rsid_map,
    docker = prs_docker,
  }

  # scatter over all (study, chrom) pairs to compute weights for each chromosome of each study
  Int n_chroms = length(chrom_list)
  Int n_weight_jobs = length(gwas_traits) * n_chroms
  scatter (job_idx in range(n_weight_jobs)) {
    Int study_idx = job_idx / n_chroms
    Int chrom_idx = job_idx % n_chroms
    call weights {
      input:
      munged_gwas = munge.munged_file_hm3[study_idx],
      N = n_total[study_idx],
      bim_file = rsid_bim.rsid_bim,
      ref_dir = ref_dir,
      chrom = chrom_list[chrom_idx],
      docker = prs_docker,
    }
  }

  # regroup the flat (study, chrom) shard paths into one row per study, so each gather shard only
  # localizes its own study's files (nested scatters are not available to slice them in WDL)
  call group_shards {
    input:
    weights = weights.weights,
    logs = weights.log,
    n = n_chroms,
  }

  # merge chrom weights into a single study-level weights file, and merge logs into a single study-level log
  scatter (i in range(length(gwas_traits))) {
    call weights_gather {
      input:
      root_name = basename(munge.munged_file_hm3[i], ".munged.hm3.gz"),
      all_weights = group_shards.weights_grouped[i],
      all_logs = group_shards.logs_grouped[i],
      rsid_bim = rsid_bim.rsid_bim,
      chrompos_bim = rsid_bim.chrompos_bim,
      docker = prs_docker,
    }
    # weights in the input dataset's own variant ID notation, used for scoring
    File score_weights = if rsid_bim.is_rsid then weights_gather.weights else weights_gather.weights_chrompos
  }

  # EVERYTHING BELOW HERE IS FOR SCORES AND SURVIVAL ANALYSIS
  # builds new table to loop over, with one row per (study, pheno, finngen_phenocode) and a boolean column indicating whether the study is in the regions file
  call gather_regions {
    input:
    sstats = validate_inputs.sstats,
    regions = regions,
  }
  Array[Array[String]] expanded_traits = read_tsv(gather_regions.expanded)

  # builds a phenotypes matrix in parquet format, with one row per FINNGENID and one column per phenotype in pheno_list_file, plus covariates for corr
  call prepare_pheno_matrix {
    input:
    phenos_file = phenos_file,
    pheno_list_file = pheno_list_file,
    covars = covars,
    docker = survival_docker,
  }

  # main action, loop over each study including with/without regions
  scatter (row in expanded_traits) {
    Int orig_idx = row[0]
    String row_pheno = row[2]
    String row_finngen_phenocode = row[4]
    Boolean row_is_no_regions = row[15] == "true"
    String row_root_name = basename(munge.munged_file_hm3[orig_idx], ".munged.hm3.gz")

    call scores {
      input:
      weights = score_weights[orig_idx],
      bed_file = bed_file,
      root_name = row_root_name,
      regions = regions,
      pheno = row_finngen_phenocode,
      is_no_regions = row_is_no_regions,
    }

    String survival_pheno = if row_finngen_phenocode != "NA" then row_finngen_phenocode else "DEATH"
    call survival {
      input:
      score_file = scores.scores,
      study = basename(scores.scores, ".sscore"),
      pheno = survival_pheno,
      age_onset = age_onset,
      docker = survival_docker,
    }

    String corr_pheno = if row_finngen_phenocode != "NA" then row_finngen_phenocode else row_pheno
    call correlate_pheno {
      input:
      score_file = scores.scores,
      pheno = corr_pheno,
      phenos_file = prepare_pheno_matrix.pheno_parquet,
      pheno_list_file = pheno_list_file,
      covars = covars,
      cpus = corr_cpus,
      docker = survival_docker,
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
    docker = survival_docker,
  }

  call sort_pheno {
    input:
    prefix = prefix,
    corr_files = correlate_pheno.corr_file,
    log_files = correlate_pheno.log,
  }

  output {
    Array[File] weights_out = weights_gather.weights
    Array[File] weights_chrompos_out = weights_gather.weights_chrompos
    Array[File] weights_logs = weights_gather.log
    Array[File] munged_files = munge.munged_file
    Array[File] munged_files_hm3 = munge.munged_file_hm3
    File munge_summary_table = munge_summary.summary
    Array[File] scores_out = scores.scores
    Array[File] scores_logs = scores.log
    File survival_fig = merge_figs.survival_fig
    File onset_fig = merge_figs.onset_fig
    File risk_fig = merge_figs.risk_fig
    File auc_fig = merge_figs.AUC_fig
    File auc_log = merge_figs.AUC_log
    File sorted_pvals = sort_pheno.sorted_pvals
    File corr_logs = sort_pheno.corr_logs
  }
}



task validate_inputs {
  input {
    File gwas_meta
    String prefix
    Boolean test
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
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 HDD"
  }

  output {
    File sstats = "sumstats.txt"
  }
}


task gather_regions {
  input {
    File sstats
    File regions
  }

  command <<<
  set -euo pipefail
  awk -F'\t' 'BEGIN{OFS="\t"}
      NR==FNR{if(FNR>1) r[$1]=1; next}
      {idx=FNR-1; print idx, $0, "false"; if ($4 in r) print idx, $0, "true"}
  ' ~{regions} ~{sstats} > expanded.tsv
  >>>

  output {
    File expanded = "expanded.tsv"
  }

  runtime {
    cpu: 1
    memory: "2 GB"
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
  }
}


task munge_summary {
  input {
    Array[File] snp_summaries
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
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 HDD"
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
  }
}


task weights_gather {
  input {
    Array[File] all_weights
    Array[File] all_logs
    String root_name
    File rsid_bim
    File chrompos_bim
    String docker
  }

  # all_weights/all_logs are this study's shards (grouped by group_shards), in chrom order. Each is
  # named "<root_name>.chr<N>.weights.*", so the root_name filter is kept as a guard against misgrouping.
  # The merged rsid weights are then translated back to chrom_pos_ref_alt IDs: the two bims are
  # row-aligned, and multiallelic sites sharing an rsid are resolved by the weight's A1/A2 pair
  # (PRScs always emits the exact bim alleles). Any unmappable or ambiguous weight fails the task.
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
  if [[ ${#weights_match[@]} -ne ~{length(all_weights)} || ${#logs_match[@]} -ne ~{length(all_logs)} ]]; then
      echo "shards not all from ~{root_name}: ${#weights_match[@]}/~{length(all_weights)} weights, ${#logs_match[@]}/~{length(all_logs)} logs" >&2
      exit 1
  fi

  cat "${weights_match[@]}" > ~{root_name}.weights.txt
  cat "${logs_match[@]}" > ~{root_name}.weights.log

  python3 - ~{chrompos_bim} ~{rsid_bim} ~{root_name}.weights.txt ~{root_name}.weights.chrompos.txt <<'EOF'
  import sys
  from collections import defaultdict
  bim, rbim, weights, out = sys.argv[1:]
  lookup = defaultdict(list)
  with open(bim) as b, open(rbim) as r:
      for n, (x, y) in enumerate(zip(b, r, strict=True), 1):
          x, y = x.split(), y.split()
          if x[0] != y[0] or x[2:] != y[2:]: sys.exit(f'bims not row-aligned at line {n}')
          lookup[y[1]].append((x[1], {x[4], x[5]}))
  rows, bad, seen = [], [], set()
  for line in open(weights):
      f = line.rstrip('\n').split('\t')
      hits = [vid for vid, alleles in lookup.get(f[1], []) if alleles == {f[3], f[4]}]
      if len(hits) != 1 or hits[0] in seen: bad.append(line); continue
      seen.add(hits[0]); f[1] = hits[0]; rows.append('\t'.join(f) + '\n')
  if bad or not rows: sys.exit(f'{len(bad)} unmappable/ambiguous weights out of {len(bad) + len(rows)}, e.g.:\n' + ''.join(bad[:5]))
  open(out, 'w').writelines(rows)
  print(f'{len(rows)} weights translated to chrom_pos_ref_alt')
  EOF
  >>>

  output {
    File weights = "~{root_name}.weights.txt"
    File weights_chrompos = "~{root_name}.weights.chrompos.txt"
    File log = "~{root_name}.weights.log"
  }

  runtime {
    docker: "~{docker}"
    cpu: 1
    disks: "local-disk 20 HDD"
  }
}


task group_shards {
  input {
    # String, not File: shard paths are only regrouped here, never localized
    Array[String] weights
    Array[String] logs
    Int n
  }

  # flat (study, chrom) shard paths -> one row of n paths per study (scatter index = study * n + chrom)
  command <<<
  set -euo pipefail
  group() { awk -v n=~{n} '{printf "%s%s", $0, (NR % n ? "\t" : "\n")} END{if (NR % n) {print "path count " NR " not a multiple of " n > "/dev/stderr"; exit 1}}' "$1"; }
  group ~{write_lines(weights)} > weights_grouped.tsv
  group ~{write_lines(logs)} > logs_grouped.tsv
  >>>

  output {
    Array[Array[String]] weights_grouped = read_tsv("weights_grouped.tsv")
    Array[Array[String]] logs_grouped = read_tsv("logs_grouped.tsv")
  }

  runtime {
    cpu: 1
    memory: "1 GB"
  }
}


task scores {
  input {
    File weights
    String bed_file
    String root_name
    File regions
    String pheno
    Boolean is_no_regions
  }

  # bim/fam/afreq share bed_file's prefix (the workflow requires them next to the input bed)
  String bim_file = sub(bed_file, ".bed$", ".bim")
  String fam_file = sub(bed_file, ".bed$", ".fam")
  String freq_file = sub(bed_file, ".bed$", ".afreq")

  Int disk_size = 20
  String out_root = if is_no_regions then root_name + ".no_regions" else root_name

  command <<<
  set -euo pipefail
  fuse_bed=$(echo "~{bed_file}" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|')
  fuse_bim=$(echo "~{bim_file}" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|')
  fuse_fam=$(echo "~{fam_file}" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|')
  fuse_freq=$(echo "~{freq_file}" | sed 's|gs://[^/]*/|/mnt/disks/gcs/|')

  EXCLUDE_ARGS=()
  if [[ "~{is_no_regions}" == "true" ]]; then
      awk -F'\t' -v p="~{pheno}" '$1==p{print $2}' ~{regions} | tr ';' '\n' | tr '_' '\t' | sed '/^$/d' > regions.txt
      EXCLUDE_ARGS=(--exclude bed1 regions.txt)
  fi

  plink2 --bed "$fuse_bed" --bim "$fuse_bim" --fam "$fuse_fam" --read-freq "$fuse_freq" \
      "${EXCLUDE_ARGS[@]}" \
      --score ~{weights} 2 4 6 center list-variants ignore-dup-ids \
      --out ~{out_root}
  >>>

  output {
    File scores = "~{out_root}.sscore"
    File log = "~{out_root}.log"
  }

  runtime {
    cpu: 16
    memory: "8 GB"
    disks: "local-disk ~{disk_size} HDD"
  }
}


task rsid_bim {
  input {
    String bed_file
    File rsid_map
    String docker
  }

  File bim = sub(bed_file, "\\.bed$", ".bim")
  String root = basename(bim, ".bim")
  Int disk_size = ceil(size(rsid_map, "GB")) + ceil(size(bim, "GB")) * 3 + 5

  # notation is detected from the first 10 IDs (rsid if they start with "rs"). The input bim is
  # copied to its own notation's output and the other one is built row-aligned (only col 2 changes):
  # chrompos -> rsid via the chrom_pos -> rsid map (unmapped positions keep their ID),
  # rsid -> chrompos as chr{CHR}_{POS}_{A2}_{A1} from the bim columns.
  command <<<
  set -euo pipefail
  python3 - ~{bim} ~{rsid_map} ~{root} <<'EOF'
  import itertools, shutil, sys
  sys.path.insert(0, '/scripts')
  from utils import load_rsid_mapping
  bim, rsid_map, root = sys.argv[1:]
  with open(bim) as f: ids = [line.split()[1] for line in itertools.islice(f, 10)]
  n_rs = sum(i.lower().startswith('rs') for i in ids)
  if not ids or 0 < n_rs < len(ids): sys.exit(f'empty bim or mixed variant ID notation in its first rows: {ids}')
  is_rsid = n_rs > 0
  if is_rsid:
      new_id, out = lambda f: f'chr{f[0]}_{f[3]}_{f[5]}_{f[4]}', f'{root}.chrompos.bim'
  else:
      to_rsid = load_rsid_mapping(rsid_map, inverse=True)
      new_id, out = lambda f: to_rsid.get(f'{f[0]}_{f[3]}', f[1]), f'{root}.rsid.bim'
  n = n_new = 0
  with open(bim) as i, open(out, 'w') as o:
      for line in i:
          f = line.split(); vid = new_id(f); n += 1; n_new += vid != f[1]; f[1] = vid; o.write('\t'.join(f) + '\n')
  shutil.copy(bim, f'{root}.rsid.bim' if is_rsid else f'{root}.chrompos.bim')
  open('is_rsid.txt', 'w').write(str(is_rsid).lower())
  print(f'{bim}: {"rsid" if is_rsid else "chrompos"} notation, built {out}')
  print(f'{n} rows, {n_new} mapped to {"chrompos" if is_rsid else "rsid"}, {n - n_new} kept original ID')
  EOF
  >>>

  output {
    File rsid_bim = "~{root}.rsid.bim"
    File chrompos_bim = "~{root}.chrompos.bim"
    Boolean is_rsid = read_boolean("is_rsid.txt")
  }

  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    disks: "local-disk ~{disk_size} HDD"
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
  Int disk_size = ceil(size(corr_files, "GB")) * 2

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
    disks: "local-disk ~{disk_size} HDD"
  }
}
