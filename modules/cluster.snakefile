#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#----------------------------------------
# @authors: Tosh, Mahesh Vangala
# @emails: , vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#----------------------------------------

from scripts.utils import _getProcessedCuffCounts

rule pca_plot:
    input:
        rpkmFile = _getProcessedCuffCounts(config),
        annotFile = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        expand("analysis/" + config["token"] + "/plots/images/pca_plot_{metacol}.png", metacol=config["metacols"]),
        pca_plot_out="analysis/" + config["token"] + "/plots/pca_plot.pdf",
        pca_out_dir = "analysis/" + config["token"] + "/plots/"
    message: "Generating PCA plots"
    shell:
        "Rscript viper/modules/scripts/pca_plot.R {input.rpkmFile} {input.annotFile} {output.pca_out_dir} "


rule heatmapSS_plot:
    input:
        rpkmFile = _getProcessedCuffCounts(config),
        annotFile=config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        ss_plot_out="analysis/" + config["token"] + "/plots/heatmapSS_plot.pdf",
        ss_txt_out="analysis/" + config["token"] + "/plots/heatmapSS.txt",
        ss_out_dir = "analysis/" + config["token"] + "/plots/"
    message: "Generating Sample-Sample Heatmap"
    shell:
        "mkdir -p analysis/plots/images && Rscript viper/modules/scripts/heatmapSS_plot.R {input.rpkmFile} "
        "{input.annotFile} {output.ss_out_dir} "


rule heatmapSF_plot:
    input:
        rpkmFile = _getProcessedCuffCounts(config),
        annotFile=config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        sf_plot_out="analysis/" + config["token"] + "/plots/heatmapSF_plot.pdf",
        sf_txt_out="analysis/" + config["token"] + "/plots/heatmapSF.txt",
        sf_out_dir = "analysis/" + config["token"] + "/plots/"
    params:
        num_kmeans_clust = config["num_kmeans_clust"]
    message: "Generating Sample-Feature heatmap"
    shell:
        "mkdir -p analysis/plots/images && Rscript viper/modules/scripts/heatmapSF_plot.R {input.rpkmFile} "
        "{input.annotFile} {params.num_kmeans_clust} {output.sf_out_dir} "
