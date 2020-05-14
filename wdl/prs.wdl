workflow prs_cs{

    String gwas_data_path
    String docker
    Boolean test
    String prefix
    call rsid_map {
        input:
        docker = docker
        }

    File gwas_meta
    call sumstats {
        input:
        gwas_meta = gwas_meta,
        docker = docker,
        test = test
        }
    
    Array[Array[String]] gwas_traits = read_tsv(sumstats.sstats)
    scatter( gwas in gwas_traits) {
        call munge {
            input :
            gwas_data_path = gwas_data_path,
            file_name = gwas[0],
            effect_type = gwas[2],
            variant = gwas[3],
            chrom = gwas[4],
            pos = gwas[5],
            ref = gwas[6],
            alt = gwas[7],
            effect = gwas[8],
            pval = gwas[9],
            rsid_map = rsid_map.rsid,
            chrompos_map = rsid_map.chrompos,
            docker = docker,
            prefix = prefix
        }
        call weights {
            input:
            munged_gwas = munge.munged_file,
            rsid_map = rsid_map.rsid,
            docker = docker,
            N = gwas[1]

        }
        call scores {
            input:
            weights = weights.weights,
            docker = docker,
            test = test
            }
    }
}

task scores {

    File weights
    String file_root = basename(weights,'.weights.txt')

    Boolean test
    String bed_string
    File bed_file = if test then sub(bed_string,'.bed','.test.bed') else bed_string
    File bim_file = sub(bed_file,'.bed','.bim')
    File fam_file = sub(bed_file,'.bed','.fam')
    File frq_file = sub(bed_file,'.bed','.afreq')
    
    String docker
    String? scores_docker
    String? final_docker = if defined(scores_docker) then scores_docker else docker
    
    Int cpu
    Int mem
    Int disk_size = ceil(size(bed_file,'GB')) + 10
    
    command <<<
    python3 /scripts/cs_scores.py \
    --weight ${weights} \
    --bed ${bed_file} \
    --out .
    >>>

    output {
        File log = "/cromwell_root/scores/${file_root}.log"
        File scores = "/cromwell_root/scores/${file_root}.sscore"
        }
    
    runtime {
        docker: "${final_docker}"
        cpu: "${cpu}"
	memory: "${mem} GB"
        disks: "local-disk ${disk_size} HDD"
        zones: "europe-west1-b"
        preemptible: 1
    }
}


task weights {

    File munged_gwas
    String root_name = basename(munged_gwas,'.munged.gz')
    
    String N
    File rsid_map
    File bim_file
    
    File file_list
    Array[File] ref_files = read_lines(file_list)

    String docker
    String? weights_docker
    String? final_docker = if defined(weights_docker) then weights_docker else docker
    Int cpu
    Int mem
    
    command <<<
    python3 /scripts/cs_wrapper.py --bim-file ${bim_file} --ref-file ${ref_files[0]}  \
    --map ${rsid_map} --out . \
    --N ${N} \
    --sum-stats ${munged_gwas} \
    --parallel 2 
    >>>

    output {
        File munged_rsid = "/cromwell_root/munge/${root_name}.munged.rsid"
        File weights = "/cromwell_root/${root_name}.weights.txt"
    }
    
    runtime {
        docker: "${final_docker}"
        cpu: "${cpu}"
	memory: "${mem} GB"
        disks: "local-disk 15 HDD"
        zones: "europe-west1-b"
        preemptible: 1
    }

    
}

task munge {

    String gwas_data_path
    String file_name
    File ss = gwas_data_path + file_name
    String prefix
    String out_root = prefix + "_" +  sub(file_name,'.gz','.munged.gz')
    
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
    String? final_docker = if defined(munge_docker) then munge_docker else docker

    command <<<
    python3 /scripts/munge.py \
    -o . \
    --ss ${ss} \
    --effect_type "${effect_type}" \
    --variant "${variant}" \
    --chrom "${chrom}" \
    --pos "${pos}" \
    --ref "${ref}" \
    --alt "${alt}" \
    --effect "${effect}" \
    --pval "${pval}" \
    --rsid-map ${rsid_map} \
    --chrompos-map ${chrompos_map} \
    --chainfile ${chainfile} \
    --prefix ${prefix}
            
    >>>

    output {
        File munged_file = "/cromwell_root/${out_root}"
        Array[File] rejected_variants = glob("/cromwell_root/tmp_parse/rejected_variants/*")
    }

    runtime {
        docker: "${final_docker}"
        cpu: "4"
	memory: "16 GB"
        disks: "local-disk 10 HDD"
        zones: "europe-west1-b"
        preemptible: 1
    }
  

}

task rsid_map {

    String docker
    String? rsid_docker
    String? final_docker = if defined(rsid_docker) then rsid_docker else docker

    File hm3_rsids    
    File bim_file
    
    command <<<
    python3 /scripts/rsid_map.py \
    -o . \
    --bim ${bim_file} \
    --rsids ${hm3_rsids} \
    --prefix hm3 
    >>>

    runtime {
        docker: "${final_docker}"
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
        }
    }


task sumstats {

    File gwas_meta
    Boolean test
    String grep = if test then " | grep MTAG" else ""

    String docker
    command <<<

    cat ${gwas_meta} | sed -E 1d ${grep}| cut -f 1,3,9,10,11,12,13,14,15,16  > sumstats.txt
    >>>

    output {
        File sstats = "./sumstats.txt"
        }
    
    runtime {
        docker: "${docker}"
        cpu: "1"
	memory: "1 GB"
        disks: "local-disk 2 HDD"
        zones: "europe-west1-b"
        preemptible: 1
        }

    }
