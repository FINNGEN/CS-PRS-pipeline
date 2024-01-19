version 1.0

workflow sandbox_prs_weights{
  input {
    File ss_meta
    String docker
    String prefix
    File rsid_map
  }

  call sumstats{
    input :
    ss_meta = ss_meta,
    docker = docker,
  }

  Array[Array[String]] ss_data = read_tsv(sumstats.sstats)
   scatter( data  in ss_data) {
    call munge {
      input :
      prefix=prefix,
      docker = docker,
      rsid_map = rsid_map,
      file_name = data[0],
      effect_type = data[3],
      variant = data[4],
      chrom = data[5],
      pos = data[6],
      ref = data[7],
      alt = data[8],
      effect = data[9],
      pval = data[10],
      build = data[11]
    }

    call weights {
      input :
      munged_ss = munge.munged_file,
      docker = docker,
      rsid_map=rsid_map,
      N = data[1],
    }
  }

  output {
    Array[File] rsid_weights = weights.weights_rsid
    Array[File] chrompos_weights = weights.weights
    Array[File] logs = weights.log
  }
}


task weights {
  input {
    String docker

    File munged_ss
    String N
    File rsid_map
    File bim_file
    File file_list 
  }

  String root_name = basename(munged_ss,'.munged.gz')
  Array[File] ref_files = read_lines(file_list)
  Int disk_size = ceil(size(munged_ss,'GB'))*2+10
  String rsid_weights = root_name + ".weights.rsid.txt"
  
  command <<<
  python3 /scripts/cs_wrapper.py --bim-file ~{bim_file} --ref-file ~{ref_files[0]}  --map ~{rsid_map} --out . --N ~{N} --sum-stats ~{munged_ss}     --parallel 1
  touch ~{root_name}.weights.log ~{root_name}.weights.txt

  >>>
  
  output {
    File munged_rsid = "/cromwell_root/munge/~{root_name}.munged.rsid"
    File log = "/cromwell_root/~{root_name}.weights.log"
    File weights = "/cromwell_root/~{root_name}.weights.txt"
    File weights_rsid =  "/cromwell_root/~{root_name}.weights.rsid.txt"
  }

  runtime {
    docker: "~{docker}"
    cpu: 8
    memory: "16 GB"
    disks: "local-disk ~{disk_size} HDD"
    zones: "europe-west1-b"
    preemptible: 1
  }
}


task munge {

  input {
    String ss_data_path
    String build
    Map[String,File] build_chains
    String file_name
    String prefix

    String effect_type
    String variant
    String chrom
    String pos
    String ref
    String alt
    String effect
    String pval
    File rsid_map
    File chrompos_map

    String docker
  }
  File chainfile = build_chains[build]
  File ss = ss_data_path + file_name
  String out_root =  prefix + "_" + sub(file_name,'.gz','.munged.gz')
  Int disk_size = ceil(size(chainfile,'GB')) + ceil(size(rsid_map,'GB')) + ceil(size(chrompos_map,'GB')) + ceil(size(ss,'GB'))*4+10
  
  command <<<
  python3 /scripts/munge.py  -o .  --ss ~{ss} --effect_type "~{effect_type}"  --variant "~{variant}"  --chrom "~{chrom}"  --pos "~{pos}"  --ref "~{ref}"   --alt "~{alt}"  --effect "~{effect}"  --pval "~{pval}"  --prefix "~{prefix}"  --rsid-map ~{rsid_map}  --chrompos-map ~{chrompos_map}  --chainfile ~{chainfile}
  ls *
  >>>

  output {
    File munged_file = "/cromwell_root/~{out_root}"
    Array[File] rejected_variants = glob("/cromwell_root/tmp_parse/rejected_variants/*")
  }
    
  runtime {
    docker: "~{docker}"
    cpu: "4"
    memory: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    zones: "europe-west1-b"
    preemptible: 1
  }
}

task sumstats {

  input {
    File ss_meta
    String docker
  }
  command <<<
   cat ~{ss_meta} | sed -E 1d | cut -f 1,3,8,9,10,11,12,13,14,15,16,17 > sumstats.txt
  >>>

  output {
    File sstats = "./sumstats.txt"
  }

  runtime {
    docker: "~{docker}"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk 2 HDD"
    zones: "europe-west1-b"
    preemptible: 1
  }
}
