#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

_HLA_regions = {'hg19':"chr6:28477797-33448354", 'mm9':'chr17:34111604-36221194'}

#call snps from the samples
#NOTE: lots of duplicated code below!--ONE SET for chr6 (default) and another
#for genome-wide
#------------------------------------------------------------------------------
# snp calling for chr6 (default)
#------------------------------------------------------------------------------
rule call_snps_hla:
    input:
        bam="analysis/STAR/{sample}/{sample}.sorted.bam",
        ref_fa=config["ref_fasta"],
    output:
        protected("analysis/snp/{sample}/{sample}.snp.hla.txt")
    params:
        varscan_path = config["varscan_path"],
        region = _HLA_regions[config['reference']]
    message: "Running varscan for snp analysis for ch6 fingerprint region"
    shell:
        "samtools mpileup -r \"{params.region}\" -f {input.ref_fa} {input.bam} | awk \'$4 != 0\' | "
        "{params.varscan_path} pileup2snp - --min-coverage 20 --min-reads2 4 > {output}"

#calculate sample snps correlation using all samples
rule sample_snps_corr_hla:
    input:
        snps = expand("analysis/snp/{sample}/{sample}.snp.hla.txt", sample=config["ordered_sample_list"]),
        metasheet = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        snp_matrix="analysis/snp/snp_corr.hla.txt",
        snp_png="analysis/plots/sampleSNPcorr_plot.hla.png",
        snp_pdf="analysis/plots/sampleSNPcorr_plot.hla.pdf"
    message: "Running snp correlations for HLA fingerprint region"
    shell:
        "{config[python2]} viper/modules/scripts/sampleSNPcorr.py {snps}> {output.snp_matrix} && "
        "Rscript viper/modules/scripts/sampleSNPcorr_plot.R {output.snp_matrix} {input.metasheet} {output.snp_png} {output.snp_pdf}"


#------------------------------------------------------------------------------
# snp calling GENOME wide (hidden config.yaml flag- 'snp_scan_genome:True'
#------------------------------------------------------------------------------

rule call_snps_genome:
    input:
        bam="analysis/STAR/{sample}/{sample}.sorted.bam",
        ref_fa=config["ref_fasta"],
    output:
        protected("analysis/snp/{sample}/{sample}.snp.genome.vcf")
    params:
        varscan_path=config["varscan_path"]
    message: "Running varscan for snp analysis genome wide"
    shell:
        "samtools mpileup -f {input.ref_fa} {input.bam} | awk \'$4 != 0\' | "
        "{params.varscan_path} mpileup2snp - --min-coverage 20 --min-reads2 4 --output-vcf > {output}"


rule snpEff_annot:
    input:
        vcf="analysis/snp/{sample}/{sample}.snp.genome.vcf"
    output:
        protected("analysis/snp/{sample}/{sample}.snpEff.annot.vcf")
    params:
        snpEff_conf=config["snpEff_conf"],
        snpEff_db=config['snpEff_db']
    message: "Running varscan for snpEff annotation analysis"
    shell:
        "snpEff -Xmx2G -c {params.snpEff_conf} {params.snpEff_db} {input.vcf} > {output}"
                                                                                                                                                                                           787,1         Bot


