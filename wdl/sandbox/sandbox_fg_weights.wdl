version 1.0

workflow finngen_weights{
  input {
    String docker
    File bim_file
    File pheno_list
    Boolean test
    String prefix
    Array[Array[String]] pheno_data = read_tsv(pheno_list)

    Boolean run_scores
    Boolean rsid_weights

  }
  
  scatter (entry in pheno_data) {
    call munge_sumstats {
      input:
      test = test,
      prefix = prefix,
      pheno = entry[0],
      bim_file = bim_file,
      docker = docker,
    }
    
    call weights {
      input:
      test = test,
      munged_gwas= munge_sumstats.munged_sumstat,
      N = entry[1],
      pheno = entry[0],
      bim_file = bim_file,
      docker = docker,
      }   
  }
  
  Array[File] final_weights = if (rsid_weights) then weights.rsid_weights else weights.weights
  
  if (run_scores) {
    call scores {
      input :
      docker = docker,
      weights = final_weights,
    }
  }
}
  
task scores {
  input {
    String plink_root
    Array[File] weights
    String docker
    Int cpu
    Int mem
  }
    
  File bed_file = plink_root + ".bed"
  File bim_file = plink_root + ".bim"
  File fam_file = plink_root + ".fam" 
  
  command <<<
    python3 /scripts/cs_scores.py \
    --weight-list ~{write_lines(weights)} \
    --bed ~{bed_file} \
    --out .
    
  >>>
  
  output {
     Array[File] logs = glob("./scores/*log")
    Array[File] scores = glob("./scores/*sscore")
  }
  
  runtime {
    docker: "~{docker}"
    cpu: "~{cpu}"
     memory: "~{mem} GB"
    disks: "local-disk ~{ceil(size(bed_file,'GB')) + 10} HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 2
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
    Int cpu
    Int mem
    File rsid_map
  }
  
  Array[File] ref_files = read_lines(ref_list)
  Int disk_size = ceil(size(munged_gwas,'GB'))*2+10
  String root_name = basename(munged_gwas,'.munged')
  String weights_root = "/cromwell_root/${root_name}.weights"
  Int pre =  if test then 2 else 0
  
  command <<<
    
    python3 /scripts/cs_wrapper.py --bim-file ~{bim_file} --ref-file ~{ref_files[0]}  \
    --out . \
     --N ~{N} \
    --sum-stats ~{munged_gwas} \
    ~{true="--test" false="" test} \
    --parallel 1
    
    python3 /scripts/convert_rsids.py \
    -o . \
    --file ~{weights_root}.txt \
    --no-header \
    --to-rsid \
    --map ~{rsid_map} \
    --metadata 1   
   >>>
   
   output {
     File weights = "${weights_root}.txt"
     File rsid_weights = "${weights_root}.rsid"
     File log = "${weights_root}.log"
     String phenos = pheno
   }
   
   runtime {
     docker: "${docker}"
     cpu: "${cpu}"
     memory: "${mem} GB"
     disks: "local-disk ~{disk_size} HDD"
     zones: "europe-west1-b europe-west1-c europe-west1-d"
     preemptible: "${pre}"
   }
 } 
 
 task munge_sumstats {
   input{
     String prefix
     String columns
     Boolean test
     
     String pheno
     File bim_file
     String file_root
     
     String docker
   }
   
   File pheno_file = sub(file_root,"PHENO",pheno)
   Int disk_size = ceil(size(pheno_file,'GB'))*2 + 10
   
   command <<<
     set -euxo pipefail
     req_cols=~{columns}
     zcat ~{pheno_file}  | awk -F'\t' -v OFS="\t" -v cols=$req_cols \
     ' BEGIN{ split(cols,printcols,","); }
     NR==1{ for (i=1;i<=NF;i++) { h[$i]=i; }
     for(k in printcols) {
       if (!(printcols[k] in h)) { print "ALL necessary columns not in sumstats. Columns needed:",cols> "/dev/stderr"; exit 1}
     }
   }
   { print "chr"$h[printcols[1]]"_"$h[printcols[2]]"_"$h[printcols[3]]"_"$h[printcols[4]],$h[printcols[4]],$h[printcols[3]],$h[printcols[5]],$h[printcols[6]]}' \
   | grep -wf <( cut -f2 ~{bim_file} ~{true=" | shuf | head -n 10000 " false="" test})  | awk -v OFS="\t" 'BEGIN {print "SNP","A1","A2","BETA","PVAL" } { print $0 }' > ~{prefix}_~{pheno}.munged
 >>>
 
 output {
   File munged_sumstat = "${prefix}_${pheno}.munged"
 }
 
 runtime {
   docker: "${docker}"
   cpu : "4"
   zones: "europe-west1-b europe-west1-c europe-west1-d"
   mem : "4 GB"
   disks : "local-disk ${disk_size} HDD"
   preemptible: 2
 }
}


