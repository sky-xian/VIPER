#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab

#---------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Aug, 09, 2016
#---------------------------

from scripts.utils import _getGCTfile

_logfile = "analysis/logs/rsem.txt"

rule rsem:
    input:
        bam="analysis/STAR/{sample}/{sample}.Aligned.toTranscriptome.out.bam"
    output:
        rsem_transcript_out = protected("analysis/rsem/{sample}/{sample}.isoforms.results"),
        rsem_genes_out = protected("analysis/rsem/{sample}/{sample}.genes.results")
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
        "rsem-calculate-expression -p {threads} {params.stranded} {params.paired_end} --bam --no-bam-output --estimate-rspd --append-names {input} {config[rsem_ref]} analysis/rsem/{params.sample_name}/{params.sample_name} > {log}"
       
rule rsem_iso_tpm_matrix:
    input:
        rsem_iso_files = expand( "analysis/rsem/{sample}/{sample}.isoforms.results", sample=config["ordered_sample_list"] ),
        metasheet = config['metasheet']
    output:
        rsem_iso_matrix = "analysis/" + config["token"] + "/rsem/iso_tpm_matrix.csv"
    message: "Running RSEM matrix generation rule for isoforms"
    benchmark:
        "benchmarks/" + config["token"] + "/rsem_iso_tpm_matrix.txt"
    run:
        args = " -f ".join( input.rsem_iso_files )
        shell("perl viper/modules/scripts/raw_and_fpkm_count_matrix.pl --column 5 --metasheet {input.metasheet} --header -f {args} 1>{output.rsem_iso_matrix}")

rule rsem_gene_tpm_matrix:
    input:
        rsem_gene_files = expand( "analysis/rsem/{sample}/{sample}.genes.results", sample=config["ordered_sample_list"] ),
        metasheet = config["metasheet"]
    output:
        rsem_gene_matrix = "analysis/" + config["token"] + "/rsem/gene_tpm_matrix.csv"
    message: "Running RSEM matrix generation rule for gene tpms"
    benchmark:
        "benchmarks/" + config["token"] + "/rsem_gene_tpm_matrix.txt"
    run:
        args = " -f ".join( input.rsem_gene_files )
        shell( "perl viper/modules/scripts/raw_and_fpkm_count_matrix.pl --column 5 --metasheet {input.metasheet} --header -f {args} 1>{output.rsem_gene_matrix}" )

rule rsem_gene_process:
    """RSEM gene.results produces duplicated names in the gene ids.  
    This rule removes that"""
    input:
        "analysis/rsem/{sample}/{sample}.genes.results"
    output:
        "analysis/rsem/{sample}/{sample}.genes.processed.txt"
    message: "Processing the RSEM genes.results output"
    benchmark: 
        "benchmarks/" + config["token"] + "/rsem_gene_process.txt"
    run:
        shell("viper/modules/scripts/rsem_process_genes.py -f {input} 1> {output}")

rule rsem_gene_ct_matrix:
    """Generate a matrix of (expected) gene counts from RSEM outputs"""
    input:
        rsem_gene_files = expand( "analysis/rsem/{sample}/{sample}.genes.processed.txt", sample=config["ordered_sample_list"] ),
        metasheet = config["metasheet"]
    output:
        rsem_gene_matrix = "analysis/" + config["token"] + "/rsem/rsem_gene_ct_matrix.csv"
    message: "Running RSEM matrix generation rule for gene counts"
    benchmark:
        "benchmarks/" + config["token"] + "/rsem_gene_ct_matrix.txt"
    run:
        args = " -f ".join( input.rsem_gene_files )
        shell( "perl viper/modules/scripts/raw_and_fpkm_count_matrix.pl --column 4 --metasheet {input.metasheet} --header -f {args} 1>{output.rsem_gene_matrix}" )

rule rsem_filter_gene_ct_matrix:
    """filters the rsem gene count.
    REPLACES: preprocess.snakefile- filter_Cuff_matrix"""
    input:
        #TODO: handle batch_effect correction
        tpmFile = "analysis/" + config["token"] + "/rsem/rsem_gene_ct_matrix.csv",
        annotFile=config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        filtered_tpm = "analysis/" + config["token"] + "/rsem/rsem_gene_ct_matrix.filtered.csv"
    params:
        sample_names = " ".join(config["ordered_sample_list"])
    message: "Generating Pre-processed RSEM TPM matrix file"
    benchmark:
        "benchmarks/" + config["token"] + "/rsem_filter_gene_ct_matrix.txt"
    shell:
        "Rscript viper/modules/scripts/rsem_filter_gene_ct_matrix.R "
        "--tpm_file {input.tpmFile} "
        "--min_samples {config[min_num_samples_expressing_at_threshold]} "
        "--TPM_cutoff {config[TPM_threshold]} "
        "--filter_miRNA {config[filter_mirna]} "
        "--numgenes {config[numgenes_plots]} "
        "--out_file {output.filtered_tpm} "
        "--sample_names {params.sample_names} "
        "Rscript viper/modules/scripts/rsem_filter_gene_ct_matrix.R "
        "--tpm_file {input.tpmFile} "
        "--min_samples {config[min_num_samples_expressing_at_threshold]} "
        "--TPM_cutoff {config[TPM_threshold]} "
        "--filter_miRNA {config[filter_mirna]} "
        "--numgenes {config[numgenes_plots]} "
        "--out_file {output.filtered_tpm} "
        "--sample_names {params.sample_names} "

rule batch_effect_removal_rsem:
    input:
        gene_cts = "analysis/" + config["token"] + "/rsem/rsem_gene_ct_matrix.csv",
        annotFile = config["metasheet"]
    output:
        rsemcsvoutput="analysis/" + config["token"] + "/rsem/batch_corrected_rsem_gene_ct_matrix.csv",
        rsempdfoutput="analysis/" + config["token"] + "/rsem/rsem_combat_qc.pdf"
    params:
        batch_column="batch",
        datatype = "rsem"
    message: "Removing batch effect from Gene Count matrix, if errors, check metasheet for batches, refer to README for specifics"
    #priority: 2
    benchmark:
        "benchmarks/" + config["token"] + "/batch_effect_removal_rsem.txt"
    shell:
        "Rscript viper/modules/scripts/batch_effect_removal.R {input.gene_cts} {input.annotFile} {params.batch_column} "
        "{params.datatype} {output.rsemcsvoutput} {output.rsempdfoutput} "
        " && mv {input.gene_cts} analysis/{config[token]}/rsemlinks/without_batch_correction_rsem_gene_ct_matrix.csv "

rule tpm_plot:
    input:
        tpm_mat = _getGCTfile(config),
        annotFile = config["metasheet"]
    output:
        tpm_png = "analysis/" + config["token"] + "/plots/gene_counts.tpm.png"
    message: "Plot gene counts at various tpm cutoffs"
    benchmark:
        "benchmarks/" + config["token"] + "/rsem_tpm_plot.txt"
    shell:
        "Rscript viper/modules/scripts/fpkm_plot.R {input.tpm_mat} {output.tpm_png}"
