#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

def getTargetInfo(config):
    targetFiles = ["analysis/STAR/STAR_Align_Report.csv",
        "analysis/STAR/STAR_Align_Report.png"]
    targetFiles.append(_getSTARcounts(config))
    targetFiles.append(["analysis/STAR/star_combat_qc.pdf", 
        "analysis/cufflinks/cuff_combat_qc.pdf"] if config["batch_effect_removal"] == "true" else[])
    return targetFiles

## Returns proper count files for with and without batch effect correction
def _getSTARcounts(config, normalized=False):
    if config["batch_effect_removal"] == "true":
        return "analysis/STAR/batch_corrected_STAR_Gene_Counts.csv"
    else:
        return "analysis/STAR/STAR_Gene_Counts.csv"

def _getCuffCounts(config, normalized=False):
    if config["batch_effect_removal"] == "true":
        return "analysis/cufflinks/batch_corrected_Cuff_Gene_Counts.csv"
    else:
        return "analysis/cufflinks/Cuff_Gene_Counts.csv"

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

def _DEsummaryOutPNG(config):
    file_list = []
    if config["comparisons"]:
        file_list.append("analysis/diffexp/de_summary.png")
    return file_list

def _runSNPgenome(config):
    ls = []
    if ('snp_scan_genome' in config) and (config['snp_scan_genome'].upper() == 'TRUE'):
        for sample in config["ordered_sample_list"]:
            #NOTE: LINE BELOW IS VERY ugly, but it's the only way it will work!
            ls.append("analysis/snp/"+sample+"/"+sample+".snp.genome.vcf")
            ls.append("analysis/snp/"+sample+"/"+sample+".snpEff.annot.vcf")
    return ls

