workflow prs_cs{

    File gwas_traits_table
    Array[Array[String]] gwas_traits = read_tsv(gwas_traits_table)
    String gwas_data_path
    String docker
    
    scatter (gwas in gwas_traits){
        call weights{
            input:
            sum_stats_name = gwas[0],
            docker = docker,
            gwas_data_path = gwas_data_path,
            N = gwas[2]
            }
        }
}


task weights {


    # INPUTS
    File file_list
    Array[File] ref_files = read_lines(file_list)
    
    String sum_stats_name        
    String gwas_data_path
    File sum_stats = gwas_data_path + sum_stats_name
    String basename = basename(sum_stats,'.munged.gz')

    Boolean test
    String test_flag = if test then " --test " else ""
    Array[String] chrom_list
    Array[String] final_chrom_list = if test then ["21","22"] else chrom_list
    
    File bim_file   
    File map_file
    Int N

    # RUNTIME OPTIONS
    String docker
    String? weights_docker
    String? final_docker = if defined(weights_docker) then weights_docker else docker

    Int disk_size = ceil(size(ref_files[0],'GB')) * 20 + 20
    Int mem
    Int cpu

    command {

        python3 -u /scripts/cs_wrapper.py \
        --ref-file ${ref_files[0] } \
        --bim-file ${bim_file} \
        --sum-stats ${sum_stats} \
        --N ${N} \
        --out /cromwell_root/results/ \
        --map ${map_file} \
        --chrom ${sep=" " final_chrom_list} \
        ${test_flag}
    }
    
    runtime {
        docker: "${final_docker}"
        cpu: "${cpu}"
	memory: "${mem} GB"
        disks: "local-disk " + "${disk_size}" + " HDD"
        zones: "europe-west1-b"
        preemptible: 1
    }
    output {
        File weights = "/cromwell_root/results/${basename}.weights.txt"
        File logs = "/cromwell_root/results/${basename}.weights.log"
        }
}
