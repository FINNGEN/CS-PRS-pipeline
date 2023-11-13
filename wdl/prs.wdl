version 1.0

workflow prs_cs{
  input {
    String gwas_data_path
    String docker
    Boolean test
    Map[String,File] build_chains
    String prefix
    File gwas_meta
  }
  
  call rsid_map {}

  call sumstats {
    input:
    gwas_meta = gwas_meta,
    docker = docker,
  }

  Array[Array[String]] gwas_traits = read_tsv(sumstats.sstats)

  scatter( gwas in gwas_traits) {
    String build = gwas[11]
    call munge {
      input :
      chainfile = build_chains[build],
      prefix=prefix,
      docker = docker,
      gwas_data_path = gwas_data_path,
      file_name = gwas[0],
      effect_type = gwas[3],
      variant = gwas[4],
      chrom = gwas[5],
      pos = gwas[6],
      ref = gwas[7],
      alt = gwas[8],
      effect = gwas[9],
      pval = gwas[10],
      rsid_map = rsid_map.rsid,
      chrompos_map = rsid_map.chrompos,
    }
    call weights {
      input:
      munged_gwas = munge.munged_file,
      rsid_map = rsid_map.rsid,
      docker = docker,
      N = gwas[1],
      test = test,
      bim_file = rsid_map.rsid_bim,
    }
    call scores {
      input:
      weights = weights.weights,
      docker = docker,
      pheno = gwas[2],
    }
  }
}

task scores {
  input {
    File weights
    File bed_file
    File regions
    String pheno
    String docker
    String? scores_docker
    Int cpu
    Int mem   
  }
  
  Int disk_size = ceil(size(bed_file,'GB')) + 10
  String? final_docker = if defined(scores_docker) then scores_docker else docker
  String file_root = basename(weights,'.weights.txt')
  File bim_file = sub(bed_file,'.bed','.bim')
  File fam_file = sub(bed_file,'.bed','.fam')
  File frq_file = sub(bed_file,'.bed','.afreq')
  
  command <<<

  cat ~{regions} | grep ~{pheno}  |  cut -f 2 | sed 's/;\t*/\n/g' | sed 's/_/\t/g' > regions.txt
  cat regions.txt

  python3 /scripts/cs_scores.py   --weight ~{weights}     --bed ~{bed_file}     --region regions.txt     --out .
  >>>

  output {
    Array[File] logs = glob("/cromwell_root/scores/~{file_root}*log")
    Array[File] scores = glob("/cromwell_root/scores/~{file_root}*sscore")
  }
  
  runtime {
    docker: "~{final_docker}"
    cpu: "~{cpu}"
    memory: "~{mem} GB"
    disks: "local-disk ~{disk_size} HDD"
    zones: "europe-west1-b"
    preemptible: 1
  }
}


task weights {

  input {
    File munged_gwas
    String N
    File rsid_map
    File bim_file
    Boolean test
    File file_list 
    String docker
    String? weights_docker
    Int cpu
    Int mem
  }
  String root_name = basename(munged_gwas,'.munged.gz')
  Array[File] ref_files = read_lines(file_list)
  Int disk_size = ceil(size(munged_gwas,'GB'))*2+10
  String? final_docker = if defined(weights_docker) then weights_docker else docker

  command <<<
  python3 /scripts/cs_wrapper.py --bim-file ~{bim_file} --ref-file ~{ref_files[0]}  --map ~{rsid_map} --out . --N ~{N} --sum-stats ~{munged_gwas}   ~{true="--test" false="" test}    --parallel 1
  
  touch ~{root_name}.weights.log ~{root_name}.weights.txt

  >>>
  
  output {
    File munged_rsid = "/cromwell_root/munge/~{root_name}.munged.rsid"
    File log = "/cromwell_root/~{root_name}.weights.log"
    File weights = "/cromwell_root/~{root_name}.weights.txt"
  }

  runtime {
    docker: "~{final_docker}"
    cpu: "~{cpu}"
    memory: "~{mem} GB"
    disks: "local-disk ~{disk_size} HDD"
    zones: "europe-west1-b"
    preemptible: 1
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
    String pval

    File rsid_map
    File chrompos_map

    File chainfile

    String docker
    String? munge_docker

    Int disk_factor

  }
  File ss = gwas_data_path + file_name
  String out_root =  prefix + "_" + sub(file_name,'.gz','.munged.gz')
  String out_cpra = prefix + "_" + sub(file_name,'.gz','.munged.cpra')
  
  Int disk_size = ceil(size(chainfile,'GB')) + ceil(size(rsid_map,'GB')) + ceil(size(chrompos_map,'GB')) + ceil(size(ss,'GB'))*disk_factor+10
  String? final_docker = if defined(munge_docker) then munge_docker else docker
  command <<<
  echo ~{disk_size} ~{disk_factor}
  
  python3 /scripts/munge.py  -o .  --ss ~{ss} --effect_type "~{effect_type}"  --variant "~{variant}"  --chrom "~{chrom}"  --pos "~{pos}"  --ref "~{ref}"   --alt "~{alt}"  --effect "~{effect}"  --pval "~{pval}"  --prefix "~{prefix}"  --rsid-map ~{rsid_map}  --chrompos-map ~{chrompos_map}  --chainfile ~{chainfile} 
  
  >>>

  output {
    File munged_file = "/cromwell_root/~{out_root}"
    File cpra_file   = "/cromwell_root/~{out_cpra}"

    Array[File] rejected_variants = glob("/cromwell_root/tmp_parse/rejected_variants/*")
  }
    
  runtime {
    docker: "~{final_docker}"
    cpu: "4"
    memory: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    zones: "europe-west1-b"
    preemptible: 1
  }
}

task rsid_map {

  input{
    File vcf_gz
    String rsid_docker 
    File hm3_rsids
    File bim_file
  }
  
  command <<<
  mkdir ./variant_mapping/
  mv ~{vcf_gz} ./variant_mapping/
    
  python3 /scripts/rsid_map.py  -o . --bim ~{bim_file}  --rsids ~{hm3_rsids}  --prefix hm3
  python3 /scripts/convert_rsids.py -o . --file ~{bim_file} --no-header --to-rsid --map ./variant_mapping/finngen.rsid.map.tsv  --metadata 1

  ls /cromwell_root/
  mv ~{sub(basename(bim_file),'.bim','.rsid')} ~{sub(basename(bim_file),'.bim','.rsid.bim')}
  >>>

  runtime {
    docker: "~{rsid_docker}"
    cpu: 4
    memory: "16 GB"
    disks: "local-disk 100 HDD"
    zones: "europe-west1-b"
    preemptible: 1
  }
  output {
    File rsid = "./variant_mapping/finngen.rsid.map.tsv"
    File chrompos = "./variant_mapping/finngen.variants.tsv"
    File hm3_snplist = "./variant_mapping/hm3.snplist"
    File rsid_bim = "/cromwell_root/" + sub(basename(bim_file),'.bim','.rsid.bim')
  }
}


task sumstats {

  input {
    File gwas_meta
    Int? last_sumstats
    String docker
  }
  command <<<

   cat ~{gwas_meta} | sed -E 1d | cut -f 1,3,8,9,10,11,12,13,14,15,16,17  ~{if defined(last_sumstats) then " | tail -n ~{last_sumstats}" else ""} > sumstats.txt
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
