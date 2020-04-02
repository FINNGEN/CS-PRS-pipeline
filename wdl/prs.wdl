workflow prs_cs{

    String gwas_data_path
    String docker
    Boolean test

    File gwas_traits_table
    File gwas_table = if test then sub(gwas_traits_table,'.csv','_test.csv') else gwas_traits_table
    Array[Array[String]] gwas_traits = read_tsv(gwas_table)
    scatter (gwas in gwas_traits){
        call weights{
            input:
            test = test,
            sum_stats_name = gwas[0],
            docker = docker,
            gwas_data_path = gwas_data_path,
            N = gwas[2]
            }
        }
    call scores{
        input:
        weight_files = weights.weights,
        docker = docker,
        test = test
        }
}


task scores{

    # WEIGHTS FROM PREVIOUS STEPS
    Array[File] weight_files   

    # SCORES INPUTS
    Boolean pgen
    Boolean hm3 
    Boolean test
    
    ## fix output file root
    String suffix = if hm3 then " --suffix hm3"  else "--suffix all" 
        
    ## fix input files
    String input_root
    String input_test = sub(input_root,'INPUTFILE','test') 
    String tmp_list = if pgen then sub(input_root,"INPUTFILE",'pgen') else sub(input_root,'INPUTFILE','plink')
    String input_list = if hm3 then  sub(tmp_list,'txt','hm3.txt') else tmp_list

    File file_list = if test then input_test else input_list
    Array[File] all_files = read_lines(file_list)

    String cmd = if pgen then " --pgen " else " --bed "

    # RUNTIME PARAMS
    String docker
    String? scores_docker
    String? final_docker = if defined(scores_docker) then scores_docker else docker
    Int cpu
    Int mem 
    Int disk_size = ceil(size(all_files[0],'GB')) + 100
    

    command {
        python3 /scripts/cs_scores.py \
        --out /cromwell_root/results/ \
        --weight-list ${write_lines(weight_files)} \
        ${cmd} ${all_files[0]}  ${suffix}
    }

    output {
        Array[File] scores = glob("/cromwell_root/results/scores/*sscore")
        Array[File] log_files = glob("/cromwell_root/results/scores/*log")
        
        }
    runtime {
        docker: "${final_docker}"
        cpu: "${cpu}"
	memory: "${mem} GB"
        disks: "local-disk " + "${disk_size}" + " HDD"       
        zones: "europe-west1-b"
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
    Int cpu
    Int mem = 6*cpu

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
