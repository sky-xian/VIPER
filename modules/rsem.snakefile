#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab

#---------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Aug, 09, 2016
#---------------------------

def getFastq(wildcards):
    return config["samples"][wildcards.sample]

rule rsem_align:
    input:
        getFastq
    output:
        rsem_transcript_out = protected("analysis/RSEM/{sample}/{sample}.isoforms.results"),
        rsem_genes_out = protected("analysis/RSEM/{sample}/{sample}.genes.results")
    threads: 8
    message: "Running RSEM on {wildcards.sample}"
    params:
        sample_name = lambda wildcards: wildcards.sample,
        stranded = "--strand-specific" if config["stranded"] else "",
        paired_end = "--paired-end" if len(config["samples"][config["ordered_sample_list"][0]]) == 2 else ""
    shell:
        "{config[rsem_path]}/rsem-calculate-expression {params.stranded} --num-threads {threads} --star --star-gzipped-read-file "
        "{params.paired_end} {input} "
        "{config[rsem_ref]} analysis/RSEM/{params.sample_name}/{params.sample_name} " 
       
rule rsem_iso_matrix:
    input:
        rsem_iso_files = expand( "analysis/RSEM/{sample}/{sample}.isoforms.results", sample=config["ordered_sample_list"] )
    output:
        rsem_iso_matrix = "analysis/" + config["token"] + "/RSEM/tpm_iso_matrix.csv"
    message: "Running RSEM matrix generation rule for isoforms"
    run:
        args = " -f ".join( input.rsem_iso_files )
        shell("perl viper/modules/scripts/raw_and_fpkm_count_matrix.pl --column 5 --header -f {args} 1>{output.rsem_iso_matrix}")


rule rsem_gene_matrix:
    input:
        rsem_gene_files = expand( "analysis/RSEM/{sample}/{sample}.genes.results", sample=config["ordered_sample_list"] )
    output:
        rsem_gene_matrix = "analysis/" + config["token"] + "/RSEM/tpm_gene_matrix.csv"
    message: "Running RSEM matrix generation rule for genes"
    run:
        args = " -f ".join( input.rsem_gene_files )
        shell( "perl viper/modules/scripts/raw_and_fpkm_count_matrix.pl --column 5 --header -f {args} 1>{output.rsem_gene_matrix}" )


