#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#----------------------------------------
# @authors: Tosh, Mahesh Vangala
# @emails: , vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#----------------------------------------

from scripts.utils import _getCuffCounts

rule pca_plot:
    input:
        rpkmFile = _getCuffCounts(config),
        annotFile=config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        expand("analysis/plots/images/pca_plot_{metacol}.png", metacol=config["metacols"]),
        pca_plot_out="analysis/plots/pca_plot.pdf"
    params:
        RPKM_threshold = config["RPKM_threshold"],
        min_num_samples_expressing_at_threshold = config["min_num_samples_expressing_at_threshold"],
        filter_mirna = config["filter_mirna"],
        SSnumgenes = config["SSnumgenes"]
    message: "Generating PCA plots"
    shell:
        "Rscript viper/modules/scripts/pca_plot_new.R {input.rpkmFile} {input.annotFile} {params.RPKM_threshold} "
        "{params.min_num_samples_expressing_at_threshold} {params.filter_mirna} {params.SSnumgenes} {output.pca_plot_out} "

rule heatmapSS_plot:
    input:
        rpkmFile = _getCuffCounts(config),
        annotFile=config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        ss_plot_out="analysis/plots/heatmapSS_plot.pdf",
        ss_txt_out="analysis/plots/heatmapSS.txt"
    params:
        RPKM_threshold = config["RPKM_threshold"],
        min_num_samples_expressing_at_threshold = config["min_num_samples_expressing_at_threshold"],
        filter_mirna = config["filter_mirna"],
        SSnumgenes = config["SSnumgenes"]
    message: "Generating Sample-Sample Heatmap"
    shell:
        "mkdir -p analysis/plots/images && Rscript viper/modules/scripts/heatmapSS_plot.R {input.rpkmFile} "
        "{input.annotFile} {params.RPKM_threshold} {params.min_num_samples_expressing_at_threshold} "
        "{params.filter_mirna} {params.SSnumgenes} {output.ss_plot_out} {output.ss_txt_out} "


rule heatmapSF_plot:
    input:
        rpkmFile = _getCuffCounts(config),
        annotFile=config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        sf_plot_out="analysis/plots/heatmapSF_plot.pdf",
        sf_txt_out="analysis/plots/heatmapSF.txt"
    params:
        RPKM_threshold = config["RPKM_threshold"],
        min_num_samples_expressing_at_threshold = config["min_num_samples_expressing_at_threshold"],
        filter_mirna = config["filter_mirna"],
        SFnumgenes = config["SFnumgenes"],
        num_kmeans_clust = config["num_kmeans_clust"]
    message: "Generating Sample-Feature heatmap"
    shell:
        "mkdir -p analysis/plots/images && Rscript viper/modules/scripts/heatmapSF_plot.R {input.rpkmFile} "
        "{input.annotFile} {params.RPKM_threshold} {params.min_num_samples_expressing_at_threshold} " 
        "{params.filter_mirna} {params.SFnumgenes} {params.num_kmeans_clust} {output.sf_plot_out} {output.sf_txt_out} "


