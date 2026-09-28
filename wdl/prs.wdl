version 1.0

# rsid-space PRS-CS pipeline: rsid_map -> validate_inputs -> [munge -> [weights x chrom] -> weights_gather -> scores] x study
#
# gwas_meta is the 18-column data/PRS_data.txt (see that file's header). validate_inputs enforces
# an injection-safe charset on every cell that gets interpolated into a task command, per
# sandbox-unmodifiable-pipelines/CLAUDE.md's input-parameter policy (this repo isn't that one,
# but this WDL is the direct basis for the Phase D sandbox port, so it's built to the same bar
# from the start rather than retrofitted later).
#
# weights is scattered over (study, chrom) rather than one task looping all 22 chromosomes per
# study: a real 22-chromosome run for one study takes ~1.5-2 days sequential wall-clock (measured
# directly), which is a bad use of Cromwell/cloud resources when each chromosome is an independent
# unit of work. A chromosome with zero overlapping variants for a given study is not a special
# case: PRScs.py's own MCMC code path was traced directly (parse_genet.py / mcmc_gtb.py) and
# confirmed to produce a valid, merely empty, per-chromosome effects file rather than erroring
# when p=0 -- so every study always scatters over the fixed chrom_list below, no dynamic
# chromosome-list computation needed.

workflow prs_cs {
  input {
    File gwas_meta
    String gwas_data_path
    String prefix

    File bed_file
    File bim_file
    File fam_file
    File hm3_rsids
    File vcf_gz
    File ref_dir_list

    File regions_file

    Map[String, File] build_chains

    # Overrides rsid_map.rsid_bim wherever an rsid-keyed bim is needed (weights' --bim-file and
    # scores' --bim override) -- lets a manually-built rsid bim be supplied directly, independent
    # of the rsid_map task's own convert_rsids.py step (currently blocked -- see rsid_map task's
    # header comment). Leave unset to use rsid_map's own output as before.
    File? rsid_bim_override

    Array[Int] chrom_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
  }

  call validate_inputs {
    input:
    gwas_meta = gwas_meta,
    prefix = prefix,
  }

  call rsid_map {
    input:
    vcf_gz = vcf_gz,
    hm3_rsids = hm3_rsids,
    bed_file = bed_file,
    bim_file = bim_file,
    fam_file = fam_file,
  }

  File rsid_bim = select_first([rsid_bim_override, rsid_map.rsid_bim])

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
      rsid_map = rsid_map.rsid,
      chainfile = build_chains[build],
    }

    scatter (chrom in chrom_list) {
      call weights {
        input:
        munged_gwas = munge.munged_file,
        rsid_map = rsid_map.rsid,
        bim_file = rsid_bim,
        ref_dir_list = ref_dir_list,
        N = n_total,
        chrom = chrom,
      }
    }

    call weights_gather {
      input:
      chrom_weights = weights.weights,
      root_name = basename(munge.munged_file, ".munged.gz"),
    }

    call scores {
      input:
      weights = weights_gather.weights,
      bed_file = bed_file,
      fam_file = fam_file,
      rsid_bim = rsid_bim,
      rsid_afreq = rsid_map.rsid_afreq,
      finngen_phenocode = finngen_phenocode,
      regions_file = regions_file,
    }
  }

  output {
    Array[File] weights_out = weights_gather.weights
    Array[File] scores_out = scores.scores
    Array[Array[File]] scores_logs = scores.logs
    Array[File] rejected_variants = munge.rejected_variants_tar
  }
}


task validate_inputs {
  input {
    File gwas_meta
    String prefix
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
    # double-quoted ~{...} interpolations these values feed into downstream.
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
    docker: "eu.gcr.io/finngen-sandbox-v3-containers/cs-prs:r14-se"
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


task rsid_map {
  input {
    File vcf_gz
    File hm3_rsids
    File bed_file
    File bim_file
    File fam_file
  }

  String bim_root = basename(bim_file, ".bim")

  command <<<
    set -euo pipefail
    mkdir -p variant_mapping
    mv ~{vcf_gz} variant_mapping/~{basename(vcf_gz)}

    python3 /scripts/rsid_map.py -o . --bim ~{bim_file} --rsids ~{hm3_rsids} --prefix hm3
    python3 /scripts/convert_rsids.py -o . --file ~{bim_file} --no-header --to-rsid \
        --map variant_mapping/finngen.rsid.map.tsv --metadata 1
    mv ~{bim_root}.rsid ~{bim_root}.rsid.bim

    # rsid-relabeled, deduped allele-frequency file for --score/--read-freq: the chrompos->rsid
    # map is position-only, so real multi-allelic sites (~0.4% of HM3 variants, confirmed) get the
    # same rsid on more than one bim row -- --read-freq hard-errors on any duplicate ID, so those
    # rows are dropped from the .afreq (not the bim/bed: .bed is positionally bound to the exact
    # .bim row order/count, trimming the bim without a matching --make-bed would silently corrupt
    # genotype-to-variant correspondence). --score's own 'ignore-dup-ids' already handles the same
    # ambiguity safely on the scoring side (drops the ambiguous ID rather than ever misapplying a
    # weight to the wrong allele -- confirmed empirically). --bed/--bim/--fam given explicitly
    # (not --bfile) since bed_file/bim_file/fam_file are three independently-localized Cromwell
    # inputs, not guaranteed to share a directory the way a plain --bfile prefix would assume.
    plink2 --bed ~{bed_file} --bim ~{bim_root}.rsid.bim --fam ~{fam_file} \
        --rm-dup exclude-mismatch --freq \
        --out ~{bim_root}.rsid
  >>>

  runtime {
    docker: "eu.gcr.io/finngen-sandbox-v3-containers/cs-prs:r14-se"
    cpu: 4
    memory: "16 GB"
    disks: "local-disk 100 HDD"
    zones: "europe-west1-b"
    preemptible: 1
  }

  output {
    File rsid = "variant_mapping/finngen.rsid.map.tsv"
    File hm3_snplist = "variant_mapping/hm3.snplist"
    File rsid_bim = "~{bim_root}.rsid.bim"
    File rsid_afreq = "~{bim_root}.rsid.afreq"
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
    docker: "eu.gcr.io/finngen-sandbox-v3-containers/cs-prs:r14-se"
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
    File rsid_map
    File bim_file
    File ref_dir_list
    String N
    Int chrom
  }

  String root_name = basename(munged_gwas, ".munged.gz")
  Array[File] ref_dirs = read_lines(ref_dir_list)
  Int disk_size = ceil(size(munged_gwas, "GB")) * 2 + 10

  command <<<
    set -euo pipefail
    # force localization of every reference-panel file, not just ref_dirs[0] -- Cromwell only
    # localizes File values actually referenced in the command, and PRScs reads sibling
    # ldblk_1kg_chrN.hdf5/snpinfo_1kg_hm3 files from ref_dirs[0]'s own directory
    : ~{sep=" " ref_dirs}

    # test mode hardcoded off: cs_wrapper.sh's --test flag (--n_iter=100 smoke test) is never
    # passed here.
    /scripts/cs_wrapper.sh \
        --ref-file ~{ref_dirs[0]} --bim-file ~{bim_file} --sum-stats ~{munged_gwas} \
        --N ~{N} --out . --chrom "~{chrom}"
  >>>

  output {
    File weights = "~{root_name}.weights.txt"
    File log = "~{root_name}.weights.log"
  }

  runtime {
    docker: "eu.gcr.io/finngen-sandbox-v3-containers/cs-prs:r14-se"
    cpu: 2
    memory: "16 GB"
    disks: "local-disk ~{disk_size} HDD"
    zones: "europe-west1-b"
    preemptible: 1
  }
}


task weights_gather {
  input {
    Array[File] chrom_weights
    String root_name
  }

  # each chrom_weights[i] is one chromosome's already rsid-keyed weight rows (cs_wrapper.sh's own
  # per-invocation merge step already produced a well-formed, single-chromosome weights.txt for
  # that shard -- possibly empty when that study had no variants on that chromosome, which is a
  # valid PRScs output, not an error, see the workflow-level comment above). Plain concatenation
  # in chrom_list order is the entire gather step.
  command <<<
    set -euo pipefail
    cat ~{sep=" " chrom_weights} > ~{root_name}.weights.txt
  >>>

  output {
    File weights = "~{root_name}.weights.txt"
  }

  runtime {
    docker: "eu.gcr.io/finngen-sandbox-v3-containers/cs-prs:r14-se"
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
    File bed_file
    File fam_file
    File rsid_bim
    File rsid_afreq
    String finngen_phenocode
    File regions_file
  }

  String file_root = basename(weights, ".weights.txt")
  Int disk_size = ceil(size(bed_file, "GB")) + 10

  command <<<
    set -euo pipefail

    # finngen_phenocode was already charset-validated by validate_inputs; regions_file is a
    # workflow File input, never user-composed text -- this grep never interpolates raw,
    # unvalidated input into a shell context.
    if [[ "~{finngen_phenocode}" != "NA" ]]; then
        grep -P "^~{finngen_phenocode}\t" ~{regions_file} | cut -f 2 | tr ';' '\n' | sed 's/_/\t/g' > regions.txt || true
    else
        : > regions.txt
    fi

    /scripts/cs_scores.sh --weight ~{weights} --bed ~{bed_file} \
        --bim ~{rsid_bim} --fam ~{fam_file} --freq ~{rsid_afreq} \
        --region regions.txt --out .
  >>>

  output {
    Array[File] logs = glob("scores/~{file_root}*log")
    Array[File] scores = glob("scores/~{file_root}*sscore")
  }

  runtime {
    docker: "eu.gcr.io/finngen-sandbox-v3-containers/bioinformatics:0.7"
    cpu: 2
    memory: "8 GB"
    disks: "local-disk ~{disk_size} HDD"
    zones: "europe-west1-b"
    preemptible: 1
  }
}
