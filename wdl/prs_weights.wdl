version 1.0

# rsid-space PRS-CS weights-only pipeline: validate_inputs -> [munge -> [weights x chrom] -> weights_gather] x study
#
# Weights-only copy of prs.wdl for the sandbox-unmodifiable-pipelines port: scoring against
# individual-level genotypes (the scores task in prs.wdl) is not something an unmodifiable
# pipeline is allowed to do, since its output (.sscore, one row per sample) is individual-level
# data and this repo's WDLs there auto-export their outputs -- see
# sandbox-unmodifiable-pipelines/CLAUDE.md's two-expert sign-off policy. This file produces only
# per-study weight files (variant IDs + effect sizes -- no individual-level content) and stops
# there; scoring stays a separate, human-run step outside of any unmodifiable pipeline.
#
# gwas_meta is the 18-column data/PRS_data.txt (see that file's header). validate_inputs enforces
# an injection-safe charset on every cell that gets interpolated into a task command, per
# sandbox-unmodifiable-pipelines/CLAUDE.md's input-parameter policy.
#
# weights is scattered over (study, chrom) rather than one task looping all 22 chromosomes per
# study: a real 22-chromosome run for one study takes ~1.5-2 days sequential wall-clock (measured
# directly), which is a bad use of Cromwell/cloud resources when each chromosome is an independent
# unit of work. A chromosome with zero overlapping variants for a given study is not a special
# case: PRScs.py's own MCMC code path was traced directly (parse_genet.py / mcmc_gtb.py) and
# confirmed to produce a valid, merely empty, per-chromosome effects file rather than erroring
# when p=0 -- so every study always scatters over the fixed chrom_list below, no dynamic
# chromosome-list computation needed.

workflow prs_cs_weights {
  input {
    File gwas_meta
    String gwas_data_path
    String prefix

    File bim_file

    File ref_dir_list

    Array[Int] chrom_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
  }

  # FinnGen's target build is hg38 (not hg19) -- confirmed via tests/test.sh's own
  # chainfile-selection logic. A study declared hg38 needs no liftover (munge.py's chrompos
  # branch skips lift() when handed an empty chainfile); hg19/hg37/hg18 studies get lifted up to
  # hg38 before their chrompos rows can be matched against the (hg38) rsid map.
  Map[String, File] build_chains = {
    "hg18": "gs://finngen-production-library-green/prs/chains/hg18ToHg38.over.chain.gz",
    "hg19": "gs://finngen-production-library-green/prs/chains/hg19ToHg38.over.chain.gz",
    "hg37": "gs://finngen-production-library-green/prs/chains/hg19ToHg38.over.chain.gz",
    "hg38": "gs://finngen-production-library-green/prs/chains/empty_file.txt"
  }

  File rsid_map = "gs://finngen-production-library-green/prs/rsid_mapping/finngen.rsid.map.tsv.pickle.chrompos"

  # general_docker: plain plink2/bash tasks with no PRScs/munge.py dependency (validate_inputs,
  # weights_gather). prs_docker: tasks that need munge.py and/or PRScs itself.
  String general_docker = "eu.gcr.io/finngen-sandbox-v3-containers/bioinformatics:0.7"
  String prs_docker = "eu.gcr.io/finngen-sandbox-v3-containers/cs-prs:r14-se"

  call validate_inputs {
    input:
    gwas_meta = gwas_meta,
    prefix = prefix,
    docker = general_docker,
  }

  Array[Array[String]] gwas_traits = read_tsv(validate_inputs.sstats)

  scatter (gwas in gwas_traits) {
    # column layout after validate_inputs' cut -- see that task's header comment
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
      docker = prs_docker,
    }

    scatter (chrom in chrom_list) {
      call weights {
        input:
        munged_gwas = munge.munged_file,
        bim_file = bim_file,
        ref_dir_list = ref_dir_list,
        N = n_total,
        chrom = chrom,
        docker = prs_docker,
      }
    }

    call weights_gather {
      input:
      chrom_weights = weights.weights,
      chrom_logs = weights.log,
      root_name = basename(munge.munged_file, ".munged.gz"),
      docker = general_docker,
    }
  }

  output {
    Array[File] weights_out = weights_gather.weights
    Array[File] weights_logs = weights_gather.log
    Array[File] rejected_variants = munge.rejected_variants_tar
  }
}


task validate_inputs {
  input {
    File gwas_meta
    String prefix
    String docker
  }

  # column layout produced here (1-indexed source columns from data/PRS_data.txt's 18-column
  # header): filename(1) pheno(2) n_total(3) finngen_phenocode(8) effect_type(9) variant(10)
  # chrom(11) pos(12) effect_allele(13) other_allele(14) effect(15) statistic(16)
  # statistic_type(17) build(18). n_cases/n_ctrls/publication/ancestry are dropped here: they're
  # never interpolated into any task command (publication URLs in particular contain characters
  # -- ':','/','?','=' -- that would fail the safe-charset check below for no reason, since
  # nothing downstream ever reads them).
  command <<<
    set -euo pipefail

    cut -f 1,2,3,8,9,10,11,12,13,14,15,16,17,18 <(sed -E 1d ~{gwas_meta}) > sumstats.txt

    # injection-safe charset check (sandbox-unmodifiable-pipelines/CLAUDE.md input-parameter
    # policy, rules 5-7): alnum, underscore, dot, colon, hash, plus, hyphen only -- covers every
    # real GWAS column-name/filename value seen in data/PRS_data.txt (e.g. '#chrom',
    # 'Chr:Position', 'p-value', a '+' in one filename) while rejecting quotes, '$', backticks,
    # ';', '|', '&', whitespace, slashes, parens -- anything that could break out of the
    # double-quoted WDL interpolations these values feed into downstream.
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

    String docker
    Int disk_factor = 4
  }

  File ss = gwas_data_path + file_name
  String out_root = prefix + "_" + sub(file_name, ".gz", ".munged.gz")
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
  >>>

  output {
    File munged_file = "~{out_root}"
    File rejected_variants_tar = "rejected_variants.tar.gz"
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


task weights {
  input {
    File munged_gwas
    File bim_file
    File ref_dir_list
    String N
    Int chrom
    String docker
  }

  String root_name = basename(munged_gwas, ".munged.gz")
  Array[String] ref_dirs = read_lines(ref_dir_list)
  Int disk_size = ceil(size(munged_gwas, "GB")) * 2 + 10

  # cs_wrapper.sh reduced to its two essential steps for a single (study, chrom) shard: unzip,
  # then run PRScs directly. Everything else in the standalone script -- CHROM_LIST/TO_RUN
  # resumability via file-existence globbing, PREFIX, KWARGS, FORCE, TEST -- existed to let one
  # invocation safely cover many chromosomes over hours on one machine; here each (study, chrom)
  # pair is already its own isolated Cromwell call, so Cromwell's own call-caching is the
  # resumability mechanism, and none of the rest was ever used by this call anyway (KWARGS/PREFIX
  # empty, FORCE/TEST off).
  command <<<
    set -euo pipefail

    BIM_PREFIX="~{sub(bim_file, '.bim$', '')}"

    SUM_STATS="~{munged_gwas}"
    if [[ "$SUM_STATS" == *.gz ]]; then
        gunzip -k "$SUM_STATS"
        SUM_STATS="${SUM_STATS%.gz}"
    fi

    python3 -u /PRScs/PRScs.py \
        --ref_dir "$(dirname "~{ref_dirs[0]}")" --bim_prefix "$BIM_PREFIX" --sst_file "$SUM_STATS" \
        --n_gwas ~{N} --out_dir ~{root_name} --chrom ~{chrom} > ~{root_name}.weights.log

    # PRScs' own per-chromosome output file (already rsid-keyed) adopted under the plain name the
    # rest of the pipeline expects -- possibly empty when this study had no variants on this
    # chromosome, which is a valid PRScs output, not an error (see workflow-level comment).
    shopt -s nullglob
    eff_files=(~{root_name}*chr~{chrom}.txt)
    shopt -u nullglob
    : > ~{root_name}.weights.txt
    for f in "${eff_files[@]}"; do cat "$f" >> ~{root_name}.weights.txt; done
  >>>

  output {
    File weights = "~{root_name}.weights.txt"
    File log = "~{root_name}.weights.log"
  }

  runtime {
    docker: "~{docker}"
    cpu: 2
    memory: "8 GB"
    disks: "local-disk ~{disk_size} HDD"
    zones: "europe-west1-b"
    preemptible: 1
  }
}


task weights_gather {
  input {
    Array[File] chrom_weights
    Array[File] chrom_logs
    String root_name
    String docker
  }

  # each chrom_weights[i] is one chromosome's already rsid-keyed weight rows (cs_wrapper.sh's own
  # per-invocation merge step already produced a well-formed, single-chromosome weights.txt for
  # that shard -- possibly empty when that study had no variants on that chromosome, which is a
  # valid PRScs output, not an error, see the workflow-level comment above). Plain concatenation
  # in chrom_list order is the entire gather step. chrom_logs is gathered the same way purely for
  # debuggability -- one place to look for what each chromosome's PRScs run actually did.
  command <<<
    set -euo pipefail
    cat ~{sep=" " chrom_weights} > ~{root_name}.weights.txt
    cat ~{sep=" " chrom_logs} > ~{root_name}.weights.log
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
