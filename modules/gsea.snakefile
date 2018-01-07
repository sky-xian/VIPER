#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule gsea:
    """Use clusterProfile enchricher fn to perform gene set enrichment analysis
    using MSigDB (default) or user defined db. Uses pvalAdj cut-off of 0.05 to 
    determine gene list.
    """
    input:
        deseq = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.deseq.csv",
    output:
        gene_list = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.gene_list.txt",
        gsea = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.gene_set.enrichment.txt",
        barplot = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.gene_set.enrichment.barplot.png",
        dotplot = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.gene_set.enrichment.dotplot.png",
    params:
        db = config["gsea_db"],
        out_path = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}",
        title = "{comparison}"
    message: "Running Gene Set Enrichment Analysis on {wildcards.comparison}"
    shell:
        "Rscript viper/modules/scripts/gsea.R {input.deseq} {params.db} \"{params.title}\" {params.out_path}"


