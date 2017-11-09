#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab

#---------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Aug, 09, 2016
#---------------------------

_logfile = "analysis/logs/rsem.txt"

rule rsem:
    input:
        bam="analysis/STAR/{sample}/{sample}.Aligned.toTranscriptome.out.bam"
    output:
        rsem_transcript_out = protected("analysis/RSEM/{sample}/{sample}.isoforms.results"),
        rsem_genes_out = protected("analysis/RSEM/{sample}/{sample}.genes.results")
    threads: 8
    message: "Running RSEM on {wildcards.sample}"
    benchmark:
        "benchmarks/{sample}/{sample}.rsem_align.txt"
    log: _logfile
    params:
        sample_name = lambda wildcards: wildcards.sample,
        stranded = "--strand-specific" if config["stranded"] else "",
        paired_end = "--paired-end" if len(config["samples"][config["ordered_sample_list"][0]]) == 2 else ""
    shell:
        "rsem-calculate-expression -p {threads} {params.stranded} {params.paired_end} --bam --no-bam-output --estimate-rspd --append-names {input} {config[rsem_ref]} analysis/RSEM/{params.sample_name}/{params.sample_name} > {log}"
       
rule rsem_iso_matrix:
    input:
        rsem_iso_files = expand( "analysis/RSEM/{sample}/{sample}.isoforms.results", sample=config["ordered_sample_list"] ),
        metasheet = config['metasheet']
    output:
        rsem_iso_matrix = "analysis/" + config["token"] + "/RSEM/tpm_iso_matrix.csv"
    message: "Running RSEM matrix generation rule for isoforms"
    benchmark:
        "benchmarks/" + config["token"] + "/rsem_iso_matrix.txt"
    run:
        args = " -f ".join( input.rsem_iso_files )
        shell("perl viper/modules/scripts/raw_and_fpkm_count_matrix.pl --column 5 --metasheet {input.metasheet} --header -f {args} 1>{output.rsem_iso_matrix}")

rule rsem_gene_matrix:
    input:
        rsem_gene_files = expand( "analysis/RSEM/{sample}/{sample}.genes.results", sample=config["ordered_sample_list"] ),
        metasheet = config["metasheet"]
    output:
        rsem_gene_matrix = "analysis/" + config["token"] + "/RSEM/tpm_gene_matrix.csv"
    message: "Running RSEM matrix generation rule for genes"
    benchmark:
        "benchmarks/" + config["token"] + "/rsem_gene_matrix.txt"
    run:
        args = " -f ".join( input.rsem_gene_files )
        shell( "perl viper/modules/scripts/raw_and_fpkm_count_matrix.pl --column 5 --metasheet {input.metasheet} --header -f {args} 1>{output.rsem_gene_matrix}" )


