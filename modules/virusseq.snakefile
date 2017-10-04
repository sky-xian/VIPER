#Stub file
rule virusseq_all:
    input:
        ["analysis/virusseq/"+sample+"/"+sample+".virusseq.transcripts.gtf" for sample in config['ordered_sample_list']],
        ["analysis/virusseq/"+sample+"/"+sample+".virusseq.filtered.gtf" for sample in config['ordered_sample_list']],
        #GENERATE .bw for each of the alignments
        ["analysis/virusseq/"+sample+"/STAR/"+sample+".virus.Aligned.sortedByCoord.out.bw" for sample in config['ordered_sample_list']],
        ["analysis/virusseq/"+sample+"/STAR/"+sample+".virus.junctions.bed" for sample in config['ordered_sample_list']],
        "analysis/" + config["token"] + "/virusseq/virusseq_table.csv",
        "analysis/" + config["token"] + "/virusseq/virusseq_summary.csv",
        "analysis/" + config["token"] + "/virusseq/virusseq_Cuff_Isoform_Counts.csv",


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
        "analysis/virusseq/{sample}/STAR/{sample}.virus.ReadsPerGene.out.tab",
        "analysis/virusseq/{sample}/STAR/{sample}.virus.SJ.out.tab"
    params:
        prefix=lambda wildcards: "analysis/virusseq/{sample}/STAR/{sample}.virus.".format(sample=wildcards.sample),
        readgroup=lambda wildcards: "ID:{sample} PL:illumina LB:{sample} SM:{sample}".format(sample=wildcards.sample)
    message: "Mapping unmapped reads to hg19Virus"
    benchmark:
        "benchmarks/{sample}/{sample}.virusseq_map.txt"
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
        "analysis/virusseq/{sample}/{sample}.virusseq.transcripts.gtf",
        "analysis/virusseq/{sample}/isoforms.fpkm_tracking"
    threads: 4
    params:
        library_command=cuff_command
    message: "Running Cufflinks on viral-mapped reads"
    benchmark:
        "benchmarks/{sample}/{sample}.virusseq_cuff.txt"
    shell:
        "cufflinks -o analysis/virusseq/{wildcards.sample} -p {threads} -G {config[virusseq_gtf_file]} {params.library_command} {input} && mv analysis/virusseq/{wildcards.sample}/transcripts.gtf analysis/virusseq/{wildcards.sample}/{wildcards.sample}.virusseq.transcripts.gtf"

rule virusseq_filterTranscripts:
    input:
        "analysis/virusseq/{sample}/{sample}.virusseq.transcripts.gtf"
    output:
        "analysis/virusseq/{sample}/{sample}.virusseq.filtered.gtf"
    benchmark:
        "benchmarks/{sample}/{sample}.virusseq_filterTranscripts.txt"
    shell:
        "viper/modules/scripts/virusseq_filterTranscripts.py -f {input} > {output}"

rule virusseq_table:
    input:
        filteredFPKMs = expand("analysis/virusseq/{sample}/{sample}.virusseq.filtered.gtf", sample=config["ordered_sample_list"]),
        readCounts = expand("analysis/virusseq/{sample}/STAR/{sample}.virus.ReadsPerGene.out.tab", sample=config["ordered_sample_list"])
    output:
        table="analysis/" + config["token"] + "/virusseq/virusseq_table.csv",
    message: "Generating virusseq output table"
    benchmark:
        "benchmarks/" + config["token"] + "/virusseq_table.txt"
    run:
        fpkms = " -f ".join(input.filteredFPKMs)
        counts = " -c ".join(input.readCounts)
        shell("viper/modules/scripts/virusseq_table.py -f {fpkms} -c {counts} > {output}")

rule virusseq_summarize:
    input:
        "analysis/" + config["token"] + "/virusseq/virusseq_table.csv",
    output:
        "analysis/" + config["token"] + "/virusseq/virusseq_summary.csv",
    message: "Summarizing virusseq output"
    benchmark:
        "benchmarks/" + config["token"] + "/virusseq_summarize.txt"
    shell:
        "viper/modules/scripts/virusseq_summarize.py -f {input} > {output}"

rule virusseq_bamToBdg:
    """Convert bam to bedGraph"""
    input:
        "analysis/virusseq/{sample}/STAR/{sample}.virus.Aligned.sortedByCoord.out.bam",
    output:
        "analysis/virusseq/{sample}/STAR/{sample}.virus.Aligned.sortedByCoord.out.bg"
    benchmark:
        "benchmarks/{sample}/{sample}.virusseq_bamToBdg.txt"
    shell:
        "bedtools genomecov -bg -split -ibam {input} -g {config[virusseq_chrom_len]} > {output}"

rule virusseq_sortBdg:
    """sort bedGraph"""
    input:
        "analysis/virusseq/{sample}/STAR/{sample}.virus.Aligned.sortedByCoord.out.bg"
    output:
        "analysis/virusseq/{sample}/STAR/{sample}.virus.Aligned.sortedByCoord.out.sorted.bg"
    benchmark:
        "benchmarks/{sample}/{sample}.virusseq_sortBdg.txt"
    shell:
        "bedSort {input} {output}"


rule virusseq_bdgToBw:
    """Convert bedGraph to bigwig"""
    input:
        "analysis/virusseq/{sample}/STAR/{sample}.virus.Aligned.sortedByCoord.out.sorted.bg"
    output:
        "analysis/virusseq/{sample}/STAR/{sample}.virus.Aligned.sortedByCoord.out.bw"
    benchmark:
        "benchmarks/{sample}/{sample}.virusseq_bdgToBw.txt"
    shell:
        "bedGraphToBigWig {input} {config[virusseq_chrom_len]} {output}"

rule virusseq_SJtab2JunctionsBed:
    """Convert STAR's SJ.out.tab to (tophat) junctions.bed BED12 format"""
    input:
        "analysis/virusseq/{sample}/STAR/{sample}.virus.SJ.out.tab"
    output:
        "analysis/virusseq/{sample}/STAR/{sample}.virus.junctions.bed"
    benchmark:
        "benchmarks/{sample}/{sample}.virusseq_SJtab2JunctionsBed.txt"
    shell:
        "viper/modules/scripts/STAR_SJtab2JunctionsBed.py -f {input} > {output}"

rule virusseq_gen_cuff_isoform_matrix:
    """Collect all of the virusseq isoform fpkms"""
    input:
        cuff_gene_fpkms=expand( "analysis/virusseq/{sample}/isoforms.fpkm_tracking", sample=config["ordered_sample_list"] ),
    output:
        "analysis/" + config["token"] + "/virusseq/virusseq_Cuff_Isoform_Counts.csv",
    message: "Generating expression matrix using cufflinks isoform counts"
    benchmark:
        "benchmarks/" + config["token"] + "/virusseq_gen_cuff_isoform_matrix.txt"
    priority: 3
    params:
        #What to call our col 0
        iid="Transcript_ID"
    run:
        fpkm_files= " -f ".join(input.cuff_gene_fpkms)
        shell("viper/modules/scripts/cuff_collect_fpkm.py -n {params.iid} -f {fpkm_files} > {output}")
