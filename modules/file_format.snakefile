#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#------------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#------------------------------------

rule get_chrom_size:
    output:
        "analysis/bam2bw/" + config['reference'] + ".Chromsizes.txt"
    params:
        config['reference']
    message: "Fetching chromosome sizes"
    shell:
        "fetchChromSizes {params} 1>{output}"
        " && if [ -e /zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/ERCC/input/ERCC92.chromInfo ]; then"
        " cat /zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/ERCC/input/ERCC92.chromInfo 1>>{output}; fi"


rule bam_to_bigwig:
    input:
        bam="analysis/STAR/{sample}/{sample}.sorted.bam",
        chrom_size="analysis/bam2bw/" + config['reference'] + ".Chromsizes.txt"
    output:
        protected("analysis/bam2bw/{sample}/{sample}.bw")
    params:
        "analysis/bam2bw/{sample}/{sample}"
    message: "Converting {wildcards.sample} bam to bigwig"
    shell:
        "bedtools genomecov -bg -split -ibam {input.bam} -g {input.chrom_size} 1> {params}.bg"
        " && bedSort {params}.bg {params}.sorted.bg"
        " && bedGraphToBigWig {params}.sorted.bg {input.chrom_size} {output}"


