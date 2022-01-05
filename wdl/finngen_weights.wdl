version 1.0

workflow finngen_weights{
  input {
    String docker
    String plink_root
    File pheno_list
    Boolean test
    String prefix
    Array[Array[String]] pheno_data = read_tsv(pheno_list)
    }

    scatter (entry in pheno_data) {

      call munge_sumstats {
      	input:
      	test = test,
      	prefix = prefix,
      	pheno = entry[0],
      	bim_file = plink_root +".bim",
      	docker = docker,
      }
      call weights {
      	input:
      	test = test,
      	munged_gwas= munge_sumstats.munged_sumstat,
      	N = entry[1],
      	pheno = entry[0],
      	bim_file = plink_root + ".bim",
      	docker = docker,
      }

    }

    call weights_chunk{ input: docker=docker,weights_list = weights.weights,pheno_list = weights.phenos }
    scatter (chunk in weights_chunk.chunk_list) {
      call scores{
      	input :
      	test = test,docker = docker,
      	weights_chunk = chunk,
      	plink_root = plink_root
      }
    }
}

task scores {

  input {
    Boolean test
    String plink_root
    File weights_chunk

    String docker
    String? scores_docker
    Int cpu
    Int mem
  }

  #String test_root = "gs://finngen-imputation-panel/sisu4/plink/sisuv4_panel_hm3"
  String test_root = "gs://r8_data/panel/sisuv4_panel_hm3"
  String root = if test then test_root else plink_root

  # i need to manipulat the input chunk to localize the file
  Array[Array[String]] by_type = transpose(read_tsv(weights_chunk))
  Array[String] phenos = by_type[1]
  Array[File] fnames = by_type[0]

  File bed_file = root + ".bed"
  File bim_file = root + ".bim"
  File fam_file = root + ".fam"
  File frq_file  = root + ".afreq"

  String? final_docker = if defined(scores_docker) then scores_docker else docker
  Int disk_size = ceil(size(bed_file,'GB')) + 10

  command <<<
    python3 /scripts/cs_scores.py \
    --weight-list ~{write_lines(fnames)} \
    --bed ~{bed_file} \
    --out .

  >>>

   output {
     Array[File] logs = glob("/cromwell_root/scores/*log")
     Array[File] scores = glob("/cromwell_root/scores/*sscore")
   }

   runtime {
     docker: "${final_docker}"
     cpu: "${cpu}"
     memory: "${mem} GB"
     disks: "local-disk ${disk_size} HDD"
     zones: "europe-west1-b europe-west1-c europe-west1-d"
     preemptible: 2
   }
 }


  # split outputs in chunks so that i pass the chunk weight to each machine
  task weights_chunk {
    input {
      Array[String] weights_list
      Array[String] pheno_list
      Int chunks
      String docker
      }

      command <<<
	paste ~{write_lines(weights_list)} ~{write_lines(pheno_list)} > tmp.txt
	split -n r/~{chunks} -d --additional-suffix=.txt  tmp.txt weight_chunk
      >>>

      output {Array[File] chunk_list =  glob("./weight_chunk*")}
      runtime {
	docker: "${docker}"
	cpu: 1
	memory: "4 GB"
	disks: "local-disk 10 HDD"
	zones: "europe-west1-b"
	preemptible: 2
	noAddress: true
    }

}

task weights {
  input {
    Boolean test

    File bim_file
    File ref_list
    File munged_gwas
    String N
    String pheno

    String docker
    String? weights_docker
    Int cpu
    Int mem
    }

    Array[File] ref_files = read_lines(ref_list)
    String? final_docker = if defined(weights_docker) then weights_docker else docker
    Int disk_size = ceil(size(munged_gwas,'GB'))*2+10
    String root_name = basename(munged_gwas,'.munged')
    Int pre =  if test then 2 else 0

    command <<<
      python3 /scripts/cs_wrapper.py --bim-file ~{bim_file} --ref-file ~{ref_files[0]}  \
      --out . \
      --N ~{N} \
      --sum-stats ~{munged_gwas} \
      ~{true="--test" false="" test} \
      --parallel 1
    >>>

    output {
      File weights = "/cromwell_root/${root_name}.weights.txt"
      File log = "/cromwell_root/${root_name}.weights.log"
      String phenos = pheno
      }

      runtime {
        docker: "${final_docker}"
	cpu: "${cpu}"
	memory: "${mem} GB"
        disks: "local-disk ~{disk_size} HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: "${pre}"
    }
  }

# munge sumstats to match prs format
task munge_sumstats {
  input{
    String prefix
    String columns
    Boolean test

    String pheno
    File bim_file
    String file_root

    String docker
    String? munge_docker
  }

  File pheno_file = sub(file_root,"PHENO",pheno)
  Int disk_size = ceil(size(pheno_file,'GB'))*2 + 10
  String? final_docker = if defined(munge_docker) then munge_docker else docker

  command <<<
    req_cols=~{columns}
    zcat ~{pheno_file}  | awk -F'\t' -v OFS="\t" -v cols=$req_cols \
      ' BEGIN{ split(cols,printcols,","); }
        NR==1{ for (i=1;i<=NF;i++) { h[$i]=i; }
               for(k in printcols) {
                 if (!(printcols[k] in h)) { print "ALL necessary columns not in sumstats. Columns needed:",cols> "/dev/stderr"; exit 1}
               }
             }
        { print "chr"$h[printcols[1]]"_"$h[printcols[2]]"_"$h[printcols[3]]"_"$h[printcols[4]],$h[printcols[4]],$h[printcols[3]],$h[printcols[5]],$h[printcols[6]]}'|\
    { awk -v OFS="\t" '{ if(NR==1) {$1="SNP"} print $1,$2,$3,$4,$5 }' ;} > ~{prefix}_~{pheno}.munged

  >>>

  output {
    File munged_sumstat = "${prefix}_${pheno}.munged"
  }

  runtime {
    docker: "${final_docker}"
    cpu : "4"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    mem : "4 GB"
    disks : "local-disk ${disk_size} HDD"
    preemptible: 2

  }

}
