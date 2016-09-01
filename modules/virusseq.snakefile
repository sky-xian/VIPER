#Stub file
rule virusseq_all:
    input:
        ["analysis/virusseq/"+sample+"/"+sample+".virusseq.transcripts.gtf" for sample in config['ordered_sample_list']],
        ["analysis/virusseq/"+sample+"/"+sample+".virusseq.filtered.gtf" for sample in config['ordered_sample_list']],
        "analysis/virusseq/virusseq_table.csv",
        "analysis/virusseq/virusseq_summary.txt",

def getUnmappedReads(wildcards):
    ls = ["analysis/STAR/%s/%s.Unmapped.out.mate1" % (wildcards.sample, wildcards.sample)]
    if len(config["samples"][config["ordered_sample_list"][0]]) == 2:
        ls.append("analysis/STAR/%s/%s.Unmapped.out.mate2" % (wildcards.sample, wildcards.sample))
    return ls

rule virusseq_map:
    input:
        getUnmappedReads
    output:
        "analysis/virusseq/{sample}/STAR/{sample}.virus.Aligned.sortedByCoord.out.bam",
        "analysis/virusseq/{sample}/STAR/{sample}.virus.ReadsPerGene.out.tab"
    params:
        prefix=lambda wildcards: "analysis/virusseq/{sample}/STAR/{sample}.virus.".format(sample=wildcards.sample),
        readgroup=lambda wildcards: "ID:{sample} PL:illumina LB:{sample} SM:{sample}".format(sample=wildcards.sample)
    message: "Mapping unmapped reads to hg19Virus"
    threads: 8
    shell:
        "STAR --runMode alignReads --runThreadN {threads} --genomeDir {config[virusseq_index]}"
        " --sjdbGTFfile {config[virusseq_gtf_file]}"
        "  --readFilesIn {input} --outFileNamePrefix {params.prefix}"
        "  --outSAMstrandField intronMotif"
        # STRANDED information DROPPED
        #"  --outSAMmode Full --outSAMattributes All {params.stranded} --outSAMattrRGline {params.readgroup} --outSAMtype BAM SortedByCoordinate"
        "  --outSAMmode Full --outSAMattributes All --outSAMattrRGline {params.readgroup} --outSAMtype BAM SortedByCoordinate"
        "  --limitBAMsortRAM 45000000000 --quantMode GeneCounts"

if( config["stranded"] ):
    cuff_command="--library-type " + config["library_type"]

rule virusseq_cuff:
    input:
        "analysis/virusseq/{sample}/STAR/{sample}.virus.Aligned.sortedByCoord.out.bam"
    output:
        "analysis/virusseq/{sample}/{sample}.virusseq.transcripts.gtf"
    threads: 4
    params:
        library_command=cuff_command
    message: "Running Cufflinks on viral-mapped reads"
    shell:
        "cufflinks -o analysis/virusseq/{wildcards.sample} -p {threads} -G {config[virusseq_gtf_file]} {params.library_command} {input} && mv analysis/virusseq/{wildcards.sample}/transcripts.gtf analysis/virusseq/{wildcards.sample}/{wildcards.sample}.virusseq.transcripts.gtf"

rule virusseq_filterTranscripts:
    input:
        "analysis/virusseq/{sample}/{sample}.virusseq.transcripts.gtf"
    output:
        "analysis/virusseq/{sample}/{sample}.virusseq.filtered.gtf"
    shell:
        "viper/modules/scripts/virusseq_filterTranscripts.py -f {input} > {output}"

rule virusseq_table:
    input:
        filteredFPKMs = expand("analysis/virusseq/{sample}/{sample}.virusseq.filtered.gtf", sample=config["ordered_sample_list"]),
        readCounts = expand("analysis/virusseq/{sample}/STAR/{sample}.virus.ReadsPerGene.out.tab", sample=config["ordered_sample_list"])
    output:
        table="analysis/virusseq/virusseq_table.csv",
    message: "Generating virusseq output table"
    run:
        fpkms = " -f ".join(input.filteredFPKMs)
        counts = " -c ".join(input.readCounts)
        shell("viper/modules/scripts/virusseq_table.py -f {fpkms} -c {counts} > {output}")

rule virusseq_summarize:
    input:
        "analysis/virusseq/virusseq_table.csv",
    output:
        "analysis/virusseq/virusseq_summary.txt",
    message: "Summarizing virusseq output"
    shell:
        "viper/modules/scripts/virusseq_summarize.py -f {input} > {output}"
