#!/usr/bin/env python3.6
from utils import make_sure_path_exists,file_exists,tmp_bash,mapcount,pretty_print
from tempfile import NamedTemporaryFile
import subprocess,shlex,argparse
import os.path

tab = '\"\\t\"'

def main(args):

    out_path = os.path.join(args.out,'variant_mapping')
    make_sure_path_exists(out_path)
    #download vcf files with snp info
    vcf_path ="ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/"
    vcf_name = "00-All.vcf.gz"
    
    if not os.path.isfile(os.path.join(out_path,vcf_name)): 
        cmd = f"wget -P {out_path} {os.path.join(vcf_path,vcf_name)}"
        subprocess.call(shlex.split(cmd))
                          
    else:
        print(f"{vcf_name} already downloaded")


    # now we proceed to extract info about the variants.
    # CHROM_POS --> RSID
    rsid_map = os.path.join(out_path,'rsid_map.tsv')
    if not os.path.isfile(rsid_map):
        print('extracting rsid to chrompos mapping from vcf')
        cmd = f"""zcat {os.path.join(out_path,vcf_name)} | grep -v "^#" |  awk '{{print $1"_"$2{tab}$3}}' | sort -k1   > {rsid_map}"""
        print(cmd)
        tmp_bash(cmd)
    else:
        print(f"{rsid_map} already generated")


    # check how many variants we are dealing with
    variant_count= os.path.join(out_path,'variant_count.txt')
    if not os.path.isfile(variant_count):
        lines = mapcount(rsid_map)
        with open(variant_count,'wt') as o: o.write(str(lines))
    else:
        with open(variant_count,'rt') as i: lines = i.readline().strip()

    print(f"{lines} present in the mapping file")


    # now we create a snp specific chrom_pos to ref/alt mapping from the bim file
    # CHROMPOS --> REF --> ALT from input .bim
    args.finngen_variants = os.path.join(out_path,"finngen.variants.tsv")
    if not os.path.isfile(args.finngen_variants):          
        #cmd = f"""cat {args.bim} |  awk '{{print $1"_"$4{tab}$5{tab}$6}}' | sort -k1 > {args.finngen_variants}"""
        cmd = f"""cat {args.bim} | cut -f2| cut -c 4-|  awk '{{gsub("_",{tab},$0); print;}}'|  awk '{{print $1"_"$2{tab}$3{tab}$4}}' | sort -k1 > {args.finngen_variants}"""
        tmp_bash(cmd)
    else:
        print(f"{args.finngen_variants} already generated")

    print(f"{mapcount(args.finngen_variants)} variants present in the finngen bim file")
     

    # now we join the two files based on the chrom_pos in order to build a rsid to chrompos mapping for finngen variants
    # RSID --> CHROMPOS
    args.finngen_rsids = os.path.join(out_path,"finngen.rsid.map.tsv")
    if not os.path.isfile(args.finngen_rsids):
        cmd = f"join {rsid_map} {args.finngen_variants}  | awk '{{print $2{tab}$1}}' | sort -k1  > {args.finngen_rsids}"
        print(cmd)
        tmp_bash(cmd)
    else:
        print(f"{args.finngen_rsids} already generated")

    print(f"{mapcount(args.finngen_rsids)} finngen variants have an rsid mapping")

    pretty_print("RSID FILTERING")
    if args.rsids:
        # sort rsids for joining purposes
        sorted_rsids = os.path.join(out_path,args.prefix + ".rsids")
        if not os.path.isfile(sorted_rsids):
            cmd = f' sort {args.rsids} > {sorted_rsids}'
            print(cmd)
            tmp_bash(cmd)
        else:
            print(f"input rsids already sorted {sorted_rsids}")
        print(f"{mapcount(sorted_rsids)} rsids provided")

        # join sorted rsids with original rsid mapping in order to have a new rsid --> chrom_pos mapping
        filtered_rsid_map = os.path.join(out_path,args.prefix + ".rsid.map.tsv")
        if not os.path.isfile(filtered_rsid_map):
            cmd = f"join {sorted_rsids}  {args.finngen_rsids} |  sort -k2> {filtered_rsid_map}"
            print(cmd)
            tmp_bash(cmd)
        else:
            print(f"subsetting of rsid map already generated")
        print(f"{mapcount(filtered_rsid_map)} rsids shared")

        # using the new chrom_pos mapping do the same in order to filter down the original list of variants.
        filtered_variants = os.path.join(out_path,args.prefix + '.variants.tsv')
        if not os.path.isfile(filtered_variants):
            cmd = f"join -2 2  {args.finngen_variants} {filtered_rsid_map} |  awk '{{print $1{tab}$2{tab}$3}}'> {filtered_variants}"
            print(cmd)
            tmp_bash(cmd)
        else:
            print(f"subsetting of variants already generated")

        print(f"{mapcount(filtered_variants)} variants in risd file.")
        
        snplist = os.path.join(out_path,args.prefix +'.snplist')
        if not os.path.isfile(snplist):
            cmd = f"""cat {filtered_variants} | awk '{{print "chr"$1"_"$2"_"$3}}' > {snplist} """
            print(cmd)
            tmp_bash(cmd)
        else:
            print(f"{snplist} already generated")

        print(f"{mapcount(snplist)} variants in snplist")
        
if __name__ == '__main__':

    """
    Creates a complete rsid to chrom_pos mapping and produces the specific mapping for a bim file (Fingen variants). Also supports rsid filtering.
    """
    
    parser = argparse.ArgumentParser(description ="Process a file")


    #Basic inputs
    parser.add_argument("-o",'--out',type = str, help = "folder in which to save the results", required = True)
    parser.add_argument("--bim",type = file_exists, help = "bim file", required = True)
    parser.add_argument("--rsids",type = file_exists, help = "optional list of rsids", default = False)
    parser.add_argument("--prefix",type = str, help = "prefix of output files for rsid filtering",default ='test')

    args = parser.parse_args()


    main(args)
