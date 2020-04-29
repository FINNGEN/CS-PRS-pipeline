workflow prs_cs{

    String gwas_data_path
    String docker
    Boolean test

    call rsid_map {
        input:
        docker = docker
        }
}


task rsid_map {

    String docker
    String? rsid_docker
    String? final_docker = if defined(rsid_docker) then rsid_docker else docker
    Int mem
    Int cpu
    
    File bim_file

    command <<<
    python3 /scripts/rsid_map.py \
    -o . \
    --bim ${bim_file} 
    >>>

    runtime {
        docker: "${final_docker}"
        cpu: "${cpu}"
	memory: "${mem} GB"
        disks: "local-disk 100 HDD"
        zones: "europe-west1-b"
        preemptible: 1
        }
    output {
        File rsid = "./variant_mapping/rsid_map.tsv"
        File chrompos = "./variant_mapping/finngen.rsid.map.tsv"
        }
    }
