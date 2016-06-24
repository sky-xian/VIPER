#!/usr/bin/env python
# vim: syntax=python tabstop=4 expandtab

rule run_STAR:
    input:
        get_fastq
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
        " && /usr/bin/samtools index {output.bam}"
