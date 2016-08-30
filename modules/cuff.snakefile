#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#--------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#--------------------------------

from scripts.utils import _getCuffCounts

cuff_command=""

if( config["stranded"] ):
    cuff_command="--library-type " + config["library_type"]


rule run_cufflinks:
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        protected("analysis/cufflinks/{sample}/{sample}.genes.fpkm_tracking")
    threads: 4
    message: "Running Cufflinks on {wildcards.sample}"
    params:
        library_command=cuff_command
    shell:
        "cufflinks -o analysis/cufflinks/{wildcards.sample} -p {threads} -G {config[gtf_file]} {params.library_command} {input}"
        " && mv analysis/cufflinks/{wildcards.sample}/genes.fpkm_tracking {output}"
        " && mv analysis/cufflinks/{wildcards.sample}/isoforms.fpkm_tracking analysis/cufflinks/{wildcards.sample}/{wildcards.sample}.isoforms.fpkm_tracking"

rule generate_cuff_matrix:
    input:
        cuff_gene_fpkms=expand( "analysis/cufflinks/{sample}/{sample}.genes.fpkm_tracking", sample=config["ordered_sample_list"] ),
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        "analysis/" + config["token"] + "/cufflinks/Cuff_Gene_Counts.csv"
    message: "Generating expression matrix using cufflinks counts"
    priority: 3
    run:
        fpkm_files= " -f ".join( input.cuff_gene_fpkms )
        shell( "perl viper/modules/scripts/raw_and_fpkm_count_matrix.pl -c -d -f {fpkm_files} 1>{output}" )


rule batch_effect_removal_cufflinks:
    input:
        cuffmat = "analysis/" + config["token"] + "/cufflinks/Cuff_Gene_Counts.csv",
        annotFile = config["metasheet"]
    output:
        cuffcsvoutput="analysis/" + config["token"] + "/cufflinks/batch_corrected_Cuff_Gene_Counts.csv",
        cuffpdfoutput="analysis/" + config["token"] + "/cufflinks/cuff_combat_qc.pdf"
    params:
        batch_column="batch",
        datatype = "cufflinks"
    message: "Removing batch effect from Cufflinks Gene Count matrix, if errors, check metasheet for batches, refer to README for specifics"
    priority: 2
    shell:
        "Rscript viper/modules/scripts/batch_effect_removal.R {input.cuffmat} {input.annotFile} {params.batch_column} "
        "{params.datatype} {output.cuffcsvoutput} {output.cuffpdfoutput} "
        " && mv {input.cuffmat} analysis/{config[token]}/cufflinks/without_batch_correction_Cuff_Gene_Counts.csv "


rule fpkm_plot:
    input:
        cuffmat = _getCuffCounts(config)[1],
        annotFile = config["metasheet"]
    output:
        fpkm_png = "analysis/" + config["token"] + "/plots/gene_counts.fpkm.png"
    message: "Plot gene counts at various fpkm cutoffs"
    shell:
        "Rscript viper/modules/scripts/fpkm_plot.R {input.cuffmat} {output.fpkm_png}"
