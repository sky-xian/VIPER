#!/usr/binmenv python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#-----------------------------------
# @authors: Tosh, Mahesh Vangala
# @emails: , vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#-----------------------------------

rule goterm_analysis:
    input:
        deseq = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.deseq.csv",
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        out_file = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.goterm.done"
    params:
        csv = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.goterm.csv",
        plot = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.goterm.pdf",
        png = "analysis/" + config["token"] + "/plots/images/{comparison}_goterm.png",
        gotermadjpvalcutoff = config["goterm_adjpval_cutoff"],
        numgoterms = config["numgoterms"],
        reference = config["reference"]
    message: "Creating Goterm Analysis plots for Differential Expressions for {wildcards.comparison}"
    shell:
        "Rscript viper/modules/scripts/goterm_analysis.R {input.deseq} {params.gotermadjpvalcutoff} "
        "{params.numgoterms} {params.reference} {params.csv} {params.plot} {params.png} && "
        " touch {output.out_file} "

rule kegg_analysis:
    input:
        deseq = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.deseq.csv",
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        out_file = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.kegg.done"
    params:
        keggpvalcutoff = config["kegg_pval_cutoff"],
        numkeggpathways = config["numkeggpathways"],
        kegg_table_up = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.kegg.up.csv",
        kegg_table_down = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.kegg.down.csv",
        keggsummary_pdf = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.keggsummary.pdf",
        keggsummary_png = "analysis/" + config["token"] + "/plots/images/{comparison}.keggsummary.png",
        gsea_table = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.gsea.csv",
        gsea_pdf = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.gsea.pdf",
        kegg_dir = "analysis/" + config["token"] + "/diffexp/{comparison}/kegg_pathways/",
        reference = config["reference"],
        temp_dir = "analysis/" + config["token"] + "/diffexp/{comparison}/temp/"
    message: "Creating Kegg Pathway Analysis for Differential Expressions for {wildcards.comparison}"
    shell:
        "mkdir {params.temp_dir} && "
        "Rscript viper/modules/scripts/kegg_pathway.R {input.deseq} {params.keggpvalcutoff} "
        "{params.numkeggpathways} {params.kegg_dir} {params.reference} {params.temp_dir} "
        "{params.kegg_table_up} {params.kegg_table_down} {params.keggsummary_pdf} " 
        "{params.keggsummary_png} {params.gsea_table} {params.gsea_pdf} && "
        "touch {output.out_file} && "
        "rm -rf {params.temp_dir} "


