#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#-----------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#-----------------------------------

def getTargetInfo(config):
    targetFiles = ["analysis/STAR/STAR_Align_Report.csv",
        "analysis/STAR/STAR_Align_Report.png"]
    targetFiles.append(_getSTARcounts(config))
    targetFiles.extend(["analysis/STAR/star_combat_qc.pdf", 
        "analysis/cufflinks/cuff_combat_qc.pdf"] if config["batch_effect_removal"] == "true" else[])
    targetFiles.extend([_getCuffCounts(config), 
                        _fusionOutput(config), 
                        _insertSizeOutput(config), 
                        _rRNAmetrics(config), 
                        _readQC(config), 
                        _bw(config),
                        _SNP(config),
                        _DE(config),
                        _cluster(config),
                        _pathway(config)])
    return targetFiles

## Returns proper count files for with and without batch effect correction
def _getSTARcounts(config):
    if config["batch_effect_removal"] == "true":
        return "analysis/STAR/batch_corrected_STAR_Gene_Counts.csv"
    else:
        return "analysis/STAR/STAR_Gene_Counts.csv"

def _getCuffCounts(config):
    cuff_files = ["analysis/plots/gene_counts.fpkm.png"]
    if config["batch_effect_removal"] == "true":
        cuff_files.append("analysis/cufflinks/batch_corrected_Cuff_Gene_Counts.csv")
    else:
        cuff_files.append("analysis/cufflinks/Cuff_Gene_Counts.csv")
    return cuff_files

def _getProcessedCuffCounts(config):
    return "analysis/cufflinks/Cuff_Gene_Counts.filtered.csv"

def _fusionOutput(config):
    fusion_out = []
    if len(config["samples"][config["ordered_sample_list"][0]]) == 2:
        fusion_out.append("analysis/STAR_Fusion/STAR_Fusion_Report.png")
    return fusion_out

def _insertSizeOutput(config):
    insert_size_out_files = []
    if len(config["samples"][config["ordered_sample_list"][0]]) == 2:
        for sample in config["ordered_sample_list"]:
            insert_size_out_files.append( "analysis/RSeQC/insert_size/" + sample + "/" + sample + ".histogram.pdf" )
    return insert_size_out_files

def _rRNAmetrics(config):
    if config["star_rRNA_index"] is not None:
        return "analysis/STAR_rRNA/STAR_rRNA_Align_Report.csv"
    else:
        return []

def _cluster(config):
    cluster_files = ["analysis/plots/pca_plot.pdf",
                    "analysis/plots/heatmapSS_plot.pdf",
                    "analysis/plots/heatmapSF_plot.pdf"]
    return cluster_files

def _DE(config):
    de_list = []
    if config["comparisons"]:
        de_list.append("analysis/diffexp/de_summary.png")
        de_list.extend([["analysis/diffexp/" + comp + "/" + comp + "_volcano.pdf",
                        "analysis/diffexp/" + comp + "/deseq_limma_fc_corr.png"]
            for comp in config["comparisons"]])
    return de_list

def _SNP(config):
    snp_files = ["analysis/plots/sampleSNPcorr_plot.hla.png"]
    if ('snp_scan_genome' in config and config['snp_scan_genome'].upper() == 'TRUE'):
        snp_files.extend([["analysis/snp/" + sample + "/" + sample + ".snp.genome.vcf", 
            "analysis/snp/" + sample + "/" + sample + ".snpEff.annot.vcf"] for sample in config["ordered_sample_list"]])
    return snp_files

def _readQC(config):
    qc_files = []
    qc_files.append("analysis/RSeQC/read_distrib/read_distrib.png")
    qc_files.append("analysis/RSeQC/gene_body_cvg/geneBodyCoverage.heatMap.png")
    qc_files.extend(["analysis/RSeQC/junction_saturation/" + sample + "/" + sample + ".junctionSaturation_plot.pdf" 
        for sample in config["ordered_sample_list"]])
    return qc_files


def _bw(config):
    bw_files = []
    bw_files.extend(["analysis/bam2bw/" + sample + "/" + sample + ".bw" for sample in config["ordered_sample_list"]])
    return bw_files


def _pathway(config):
    path_files = []
    path_files.extend([["analysis/diffexp/" + comp + "/" + comp + ".goterm.done",
                        "analysis/diffexp/" + comp + "/" + comp + ".kegg.done"] 
        for comp in config["comparisons"]])
    return path_files






