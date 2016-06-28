#!/usr/bin/env python
# vim: syntax=python tabstop=4 expandtab

def getFastq(wildcards):
    return config["samples"][wildcards.sample]

strand_command=""
rRNA_strand_command=""

if( config["stranded"] ):
    strand_command="--outFilterIntronMotifs RemoveNoncanonical"
    rRNA_strand_command="--outFilterIntronMotifs RemoveNoncanonical"
else:
    strand_command="--outSAMstrandField intronMotif"
    rRNA_strand_command="--outSAMstrandField intronMotif"

run_fusion= True if len(config["samples"][config["ordered_sample_list"][0]]) == 2 else False
gz_command="--readFilesCommand zcat" if config["samples"][config["ordered_sample_list"][0]][0][-3:] == '.gz' else ""

if( run_fusion ):
    if( config["stranded"] ):
        strand_command = """ --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outReadsUnmapped None 
                             --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 
                             --alignMatesGapMax 200000 --alignIntronMax 200000 """
    else:
        strand_command = """ --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outReadsUnmapped None 
                             --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10
                             --alignMatesGapMax 200000 --alignIntronMax 200000 --outSAMstrandField intronMotif """


rule run_STAR:
    input:
        getFastq
    output:
        bam=protected("analysis/STAR/{sample}/{sample}.sorted.bam"),
        counts="analysis/STAR/{sample}/{sample}.counts.tab",
        log_file="analysis/STAR/{sample}/{sample}.Log.final.out"
    params:
        stranded=strand_command,
        gz_support=gz_command,
        prefix=lambda wildcards: "analysis/STAR/{sample}/{sample}".format(sample=wildcards.sample),
        readgroup=lambda wildcards: "ID:{sample} PL:illumina LB:{sample} SM:{sample}".format(sample=wildcards.sample)
    threads: 8
    message: "Running STAR Alignment on {wildcards.sample}"
    shell:
        "STAR --runMode alignReads --runThreadN {threads} --genomeDir {config[star_index]}"
        " --sjdbGTFfile {config[gtf_file]}"
        " --readFilesIn {input} {params.gz_support} --outFileNamePrefix {params.prefix}."
        "  --outSAMstrandField intronMotif"
        "  --outSAMmode Full --outSAMattributes All {params.stranded} --outSAMattrRGline {params.readgroup} --outSAMtype BAM SortedByCoordinate"
        "  --limitBAMsortRAM 45000000000 --quantMode GeneCounts"
        " && mv {params.prefix}.Aligned.sortedByCoord.out.bam {output.bam}"
        " && mv {params.prefix}.ReadsPerGene.out.tab {output.counts}"
        " && samtools index {output.bam}"


rule generate_STAR_report:
    input:
        star_log_files=expand( "analysis/STAR/{sample}/{sample}.Log.final.out", sample=config["ordered_sample_list"] ),
        star_gene_count_files=expand( "analysis/STAR/{sample}/{sample}.counts.tab", sample=config["ordered_sample_list"] ),
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        csv="analysis/STAR/STAR_Align_Report.csv",
        png="analysis/STAR/STAR_Align_Report.png",
        gene_counts="analysis/STAR/STAR_Gene_Counts.csv"
    message: "Generating STAR report"
    priority: 3
    run:
        log_files = " -l ".join( input.star_log_files )
        count_files = " -f ".join( input.star_gene_count_files )
        shell( "perl viper/modules/scripts/STAR_reports.pl -l {log_files} 1>{output.csv}" )
        shell( "Rscript viper/modules/scripts/map_stats.R {output.csv} {output.png}" )
        shell( "perl viper/modules/scripts/raw_and_fpkm_count_matrix.pl -f {count_files} 1>{output.gene_counts}" )


rule batch_effect_removal_star:
    input:
        starmat = "analysis/STAR/STAR_Gene_Counts.csv",
        annotFile = config["metasheet"]
    output:
        starcsvoutput="analysis/STAR/batch_corrected_STAR_Gene_Counts.csv",
        starpdfoutput="analysis/STAR/star_combat_qc.pdf"
    params:
        batch_column="batch",
        datatype = "star"
    message: "Removing batch effect from STAR Gene Count matrix, if errors, check metasheet for batches, refer to README for specifics"
    priority: 2
    shell:
        "Rscript viper/modules/scripts/batch_effect_removal.R {input.starmat} {input.annotFile} "
        "{params.batch_column} {params.datatype} {output.starcsvoutput} {output.starpdfoutput} "
        " && mv {input.starmat} analysis/STAR/without_batch_correction_STAR_Gene_Counts.csv"


rule run_STAR_fusion:
    input:
        bam="analysis/STAR/{sample}/{sample}.sorted.bam" #just to make sure STAR output is available before STAR_Fusion
    output:
        protected("analysis/STAR_Fusion/{sample}/{sample}.fusion_candidates.final")
    log:
        "analysis/STAR_Fusion/{sample}/{sample}.star_fusion.log"
    message: "Running STAR fusion on {wildcards.sample}"
    shell:
        "STAR-Fusion --chimeric_junction analysis/STAR/{wildcards.sample}/{wildcards.sample}.Chimeric.out.junction "
        "--genome_lib_dir {config[genome_lib_dir]} --output_dir analysis/STAR_Fusion/{wildcards.sample} >& {log}"
        " && mv analysis/STAR_Fusion/{wildcards.sample}/star-fusion.fusion_candidates.final {output}"
        " && mv analysis/STAR_Fusion/{wildcards.sample}/star-fusion.fusion_candidates.final.abridged"
        " analysis/STAR_Fusion/{wildcards.sample}/{wildcards.sample}.fusion_candidates.final.abridged"
        " && touch {output}" # For some sample, final.abridged is created but not .final file; temp hack before further investigate into this


rule run_STAR_fusion_report:
    input:
        sf_list = expand("analysis/STAR_Fusion/{sample}/{sample}.fusion_candidates.final.abridged", sample=config["ordered_sample_list"]),
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        csv="analysis/STAR_Fusion/STAR_Fusion_Report.csv",
        png="analysis/STAR_Fusion/STAR_Fusion_Report.png"
    message: "Generating STAR fusion report"
    shell:
        "python viper/modules/scripts/STAR_Fusion_report.py -f {input.sf_list} 1>{output.csv} "
        "&& Rscript viper/modules/scripts/STAR_Fusion_report.R {output.csv} {output.png}"


rule run_rRNA_STAR:
    input:
        getFastq
    output:
        bam=protected("analysis/STAR_rRNA/{sample}/{sample}.sorted.bam"),
        log_file="analysis/STAR_rRNA/{sample}/{sample}.Log.final.out"
    params:
        stranded=rRNA_strand_command,
        prefix=lambda wildcards: "analysis/STAR_rRNA/{sample}/{sample}".format(sample=wildcards.sample),
        readgroup=lambda wildcards: "ID:{sample} PL:illumina LB:{sample} SM:{sample}".format(sample=wildcards.sample)
    threads: 8
    message: "Running rRNA STAR for {wildcards.sample}"
    shell:
        "STAR --runMode alignReads --runThreadN {threads} --genomeDir {config[star_rRNA_index]}"
        " --readFilesIn {input} --readFilesCommand zcat --outFileNamePrefix {params.prefix}."
        "  --outSAMmode Full --outSAMattributes All {params.stranded} --outSAMattrRGline {params.readgroup} --outSAMtype BAM SortedByCoordinate"
        "  --limitBAMsortRAM 45000000000"
        " && mv {params.prefix}.Aligned.sortedByCoord.out.bam {output.bam}"
        " && samtools index {output.bam}"


rule generate_rRNA_STAR_report:
    input:
        star_log_files=expand( "analysis/STAR_rRNA/{sample}/{sample}.Log.final.out", sample=config["ordered_sample_list"] ),
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        csv="analysis/STAR_rRNA/STAR_rRNA_Align_Report.csv",
        png="analysis/STAR_rRNA/STAR_rRNA_Align_Report.png"
    message: "Generating STAR rRNA report"
    run:
        log_files = " -l ".join( input.star_log_files )
        shell( "perl viper/modules/scripts/STAR_reports.pl -l {log_files} 1>{output.csv}" )
        shell( "Rscript viper/modules/scripts/map_stats_rRNA.R {output.csv} {output.png}" )

