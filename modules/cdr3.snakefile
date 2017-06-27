#CDR3 VIPER Module- performs CDR3 analysis on the samples
#Q: does this really need to be Tokenized?
#A: I don't think so b/c it's intimately tied to the SAMPLE/mapping
rule cdr3_all:
    input:
        ["analysis/cdr3/"+sample+"/"+sample+".sorted.bam.fa" for sample in config['ordered_sample_list']],
        ["analysis/cdr3/CPK.csv"],
        ["analysis/cdr3/CPK.png"],

rule index_sortedBam:
    """index the bam"""
    input:
        bam="analysis/STAR/{sample}/{sample}.sorted.bam",
    output:
        bai="analysis/STAR/{sample}/{sample}.sorted.bam.bai",
    message:
        "CDR3: indexing bam file"
    shell:
        "samtools index {input}"

rule CDR3_TRUST:
    """perform CDR3 analysis using TRUST2"""
    input:
        bam="analysis/STAR/{sample}/{sample}.sorted.bam",
        bai="analysis/STAR/{sample}/{sample}.sorted.bam.bai",
    output:
        fa="analysis/cdr3/{sample}/{sample}.sorted.bam.fa",
        log="analysis/cdr3/{sample}/{sample}_TRUST.log",
    params:
        outdir="analysis/cdr3/{sample}/",
        pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
        genome=config["reference"],
    message:
        "CDR3: Performing CDR3 analysis using TRUST"
    shell:
        "{params.pypath} {config[python2]} {config[trust_path]} -f {input.bam} -o {params.outdir} -a -g {params.genome} > {output.log}"

rule calc_CPK:
    """RULE to calculate the CPK of each trust cdr .fa output
    NOTE: CPK = (# lines /2) *1000 / (first number of est_lib_size)
    """
    input:
        #NOTE: cdr3_calc_cpk.py has been modified to take directories, but
        #because snakemake inputs need links to outputs, we keep the longer
        #files-as-input call
        cdr_files = expand("analysis/cdr3/{sample}/{sample}.sorted.bam.fa", sample=config["ordered_sample_list"]),
    output:
        "analysis/cdr3/CPK.csv"
    message:
        """CDR3: calculating CPKs"""
    run:
        files = " -f ".join(input.cdr_files)
        shell("viper/modules/scripts/cdr3_calc_cpk.py -f {files} > {output}")

rule CPK_boxplot:
    """RULE to generate CPK boxplot image based on cdr3/CPK.csv"""
    input:
        "analysis/cdr3/CPK.csv"
    output:
        "analysis/cdr3/CPK.png"
    message:
        """CDR3: generating CPK boxplot"""
    shell:
        "Rscript viper/modules/scripts/cdr3_cpk_plot.R {input} {output}"
