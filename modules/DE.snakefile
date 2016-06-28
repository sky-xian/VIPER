#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.utils import _getSTARcounts
import pandas as pd

metadata = pd.read_table(config['metasheet'], index_col=0, sep=',')

def _getColumn(comparison):
    return metadata["comp_{}".format(comparison)]

def _getComparison(name, group):
    comp = _getColumn(name)
    return metadata[comp == group].index

def _getSamples(wildcards):
    comp = _getColumn(wildcards.comparison)
    return comp.dropna().index


rule limma_and_deseq:
    input:
        counts = _getSTARcounts(config)
    output:
        limma = "analysis/diffexp/{comparison}/{comparison}.limma.csv",
        deseq = "analysis/diffexp/{comparison}/{comparison}.deseq.csv",
        deseqSum = "analysis/diffexp/{comparison}/{comparison}.deseq.sum.csv",
        #annotations
        limma_annot = "analysis/diffexp/{comparison}/{comparison}.limma.annot.csv",
        deseq_annot = "analysis/diffexp/{comparison}/{comparison}.deseq.annot.csv",
    params:
        s1=lambda wildcards: ",".join(_getComparison(wildcards.comparison, 1)),
        s2=lambda wildcards: ",".join(_getComparison(wildcards.comparison, 2)),
        gene_annotation = config['gene_annotation']
    message: "Running differential expression analysis using limma and deseq for {wildcards.comparison}"
    shell:
        """ Rscript viper/modules/scripts/DEseq.R \"{input.counts}\" \"{params.s1}\" \"{params.s2}\" 
            {output.limma} {output.deseq} {output.limma_annot} {output.deseq_annot} 
            {output.deseqSum} {params.gene_annotation} """

rule deseq_limma_fc_plot:
    input:
        deseq = "analysis/diffexp/{comparison}/{comparison}.deseq.csv",
        limma = "analysis/diffexp/{comparison}/{comparison}.limma.csv"
    output:
        out_csv = "analysis/diffexp/{comparison}/deseq_limma_fc_corr.csv",
        out_png = "analysis/diffexp/{comparison}/deseq_limma_fc_corr.png"
    shell:
        "Rscript viper/modules/scripts/deseq_limma_fc_corr.R {input.deseq} {input.limma} {output.out_csv} {output.out_png}"


rule fetch_DE_gene_list:
    input:
        deseq_file_list=expand("analysis/diffexp/{comparison}/{comparison}.deseq.csv",comparison=config["comparisons"]),
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        csv="analysis/diffexp/de_summary.csv",
        png="analysis/diffexp/de_summary.png"
    message: "Creating Differential Expression summary"
    run:
        deseq_file_string = ' -f '.join(input.deseq_file_list)
        shell("perl viper/modules/scripts/get_de_summary_table.pl -f {deseq_file_string} 1>{output.csv}")
        shell("Rscript viper/modules/scripts/de_summary.R {output.csv} {output.png}")


#Generate volcano plots for each comparison
rule volcano_plot:
    input:
        deseq = "analysis/diffexp/{comparison}/{comparison}.deseq.csv",
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        plot = "analysis/diffexp/{comparison}/{comparison}_volcano.pdf",
        png = "analysis/plots/images/{comparison}_volcano.png"
    message: "Creating volcano plots for Differential Expressions for {wildcards.comparison}"
    shell:
        "Rscript viper/modules/scripts/volcano_plot.R {input.deseq} {output.plot} {output.png}"


