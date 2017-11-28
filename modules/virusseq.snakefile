#Virusseq module
_logfile = "analysis/logs/virusseq.txt"

rule virusseq_all:
    input:
        ["analysis/virusseq/"+sample+"/rsem/"+sample+".virusseq.genes.results" for sample in config['ordered_sample_list']],
        ["analysis/virusseq/"+sample+"/"+sample+".virusseq.filtered.genes" for sample in config['ordered_sample_list']],
        #GENERATE .bw for each of the alignments
        ["analysis/virusseq/"+sample+"/STAR/"+sample+".virus.Aligned.sortedByCoord.out.bw" for sample in config['ordered_sample_list']],
        ["analysis/virusseq/"+sample+"/STAR/"+sample+".virus.junctions.bed" for sample in config['ordered_sample_list']],
        "analysis/" + config["token"] + "/virusseq/virusseq_table.csv",
        "analysis/" + config["token"] + "/virusseq/virusseq_summary.csv",
        "analysis/" + config["token"] + "/virusseq/virusseq_Isoform_Counts.csv",


def getUnmappedReads(wildcards):
    ls = ["analysis/STAR/%s/%s.Unmapped.out.mate1" % (wildcards.sample, wildcards.sample)]
    if len(config["samples"][config["ordered_sample_list"][0]]) == 2:
        ls.append("analysis/STAR/%s/%s.Unmapped.out.mate2" % (wildcards.sample, wildcards.sample))
    return ls

rule virusseq_map:
    input:
        getUnmappedReads
    output:
        bam="analysis/virusseq/{sample}/STAR/{sample}.virus.Aligned.sortedByCoord.out.bam",
        counts="analysis/virusseq/{sample}/STAR/{sample}.virus.Aligned.toTranscriptome.out.bam",
        sjtab="analysis/virusseq/{sample}/STAR/{sample}.virus.SJ.out.tab"
    params:
        prefix=lambda wildcards: "analysis/virusseq/{sample}/STAR/{sample}.virus.".format(sample=wildcards.sample),
        readgroup=lambda wildcards: "ID:{sample} PL:illumina LB:{sample} SM:{sample}".format(sample=wildcards.sample)
    message: "Mapping unmapped reads to hg19Virus"
    benchmark:
        "benchmarks/{sample}/{sample}.virusseq_map.txt"
    threads: 8
    shell:
        "STAR --runMode alignReads --runThreadN {threads} --genomeDir {config[virusseq_index]}"
        "  --readFilesIn {input} --outFileNamePrefix {params.prefix}"
        "  --outSAMstrandField intronMotif"
        "  --outSAMmode Full --outSAMattributes All --outSAMattrRGline {params.readgroup} --outSAMtype BAM SortedByCoordinate"
        "  --limitBAMsortRAM 45000000000 --quantMode TranscriptomeSAM"

rule virusseq_rsem:
    """Quantify virusseq transcripts using RSEM"""
    input:
        bam="analysis/virusseq/{sample}/STAR/{sample}.virus.Aligned.toTranscriptome.out.bam"
    output:
        rsem_transcript_out = protected("analysis/virusseq/{sample}/rsem/{sample}.virusseq.isoforms.results"),
        rsem_genes_out = protected("analysis/virusseq/{sample}/rsem/{sample}.virusseq.genes.results")
    threads: 8
    message: "Running virusseq RSEM on {wildcards.sample}"
    benchmark:
        "benchmarks/{sample}/{sample}.virusseq_rsem.txt"
    log: _logfile
    params:
        sample_name = lambda wildcards: wildcards.sample,
        stranded = "--strand-specific" if config["stranded"] else "",
        paired_end = "--paired-end" if len(config["samples"][config["ordered_sample_list"][0]]) == 2 else ""
    shell:
        "rsem-calculate-expression -p {threads} {params.stranded} {params.paired_end} --bam --no-bam-output --estimate-rspd --append-names {input} {config[virusseq_rsem_ref]} analysis/virusseq/{params.sample_name}/rsem/{params.sample_name}.virusseq > {log}"

rule virusseq_processRsem:
    """Remove duplicated gene_ids from rsem results"""
    input: 
        "analysis/virusseq/{sample}/rsem/{sample}.virusseq.genes.results"
    output:
        "analysis/virusseq/{sample}/rsem/{sample}.virusseq.genes.processed.txt"
    message: "Processing the RSEM gene.results output"
    benchmark: 
        "benchmarks/{sample}/{sample}.virusseq_processRsem.txt"
    run:
        shell("viper/modules/scripts/rsem_process_genes.py -f {input} 1> {output}")

rule virusseq_annoteRsem:
    """Add chr, start, end to rsem results (additional cols)"""
    input: 
        "analysis/virusseq/{sample}/rsem/{sample}.virusseq.genes.processed.txt"
    output:
        "analysis/virusseq/{sample}/rsem/{sample}.virusseq.genes.annot.txt"
    message: "Annotating the RSEM gene.results output"
    benchmark: 
        "benchmarks/{sample}/{sample}.virusseq_annotateRsem.txt"
    params:
        gtf=config['virusseq_ucsc_gtf']
    run:
        shell("viper/modules/scripts/virusseq_annotate_rsem.py -b {input} -g {params.gtf} 1> {output}")

rule virusseq_filterTranscripts:
    """Filter for potential viral transcripts, i.e. on chrM and TPM > 0.0"""
    input:
        "analysis/virusseq/{sample}/rsem/{sample}.virusseq.genes.annot.txt"
    output:
        "analysis/virusseq/{sample}/{sample}.virusseq.filtered.genes"
    benchmark:
        "benchmarks/{sample}/{sample}.virusseq_filterTranscripts.txt"
    shell:
        "viper/modules/scripts/virusseq_filterTranscripts.py -f {input} > {output}"

rule virusseq_table:
    input:
        filteredTPMs = expand("analysis/virusseq/{sample}/{sample}.virusseq.filtered.genes", sample=config["ordered_sample_list"]),
    output:
        table="analysis/" + config["token"] + "/virusseq/virusseq_table.csv",
    message: "Generating virusseq output table"
    benchmark:
        "benchmarks/" + config["token"] + "/virusseq_table.txt"
    run:
        tpms = " -f ".join(input.filteredTPMs)
        shell("viper/modules/scripts/virusseq_table.py -f {tpms} > {output}")

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
        "viper/modules/scripts/STAR_SJtab2JunctionsBed.py -f {input} >{output}"

rule virusseq_rsem_isoform_matrix:
    """Collect all of the virusseq isoform tpms"""
    input:
        rsem_iso=expand( "analysis/virusseq/{sample}/rsem/{sample}.virusseq.isoforms.results", sample=config["ordered_sample_list"] ),
    output:
        "analysis/" + config["token"] + "/virusseq/virusseq_Isoform_Counts.csv"
    message: "Generating expression matrix using isoform counts"
    benchmark:
        "benchmarks/" + config["token"] + "/virusseq_rsem_isoform_matrix.txt"
    run:
        rsem_files= " -f ".join(input.rsem_iso)
        shell("viper/modules/scripts/virusseq_rsem_collect_tpm.py -f {rsem_files} > {output}")
    
