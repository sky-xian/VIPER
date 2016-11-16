#CDR3 VIPER Module- performs CDR3 analysis on the samples
#Q: does this really need to be Tokenized?
#A: I don't think so b/c it's intimately tied to the SAMPLE/mapping
rule cdr3_all:
    input:
        ["analysis/cdr3/"+sample+"/"+sample+".sorted.bam.fa" for sample in config['ordered_sample_list']],
        ["analysis/cdr3/"+sample+"/"+sample+".sorted.bam-Locs.bam" for sample in config['ordered_sample_list']],
        ["analysis/cdr3/"+sample+"/"+sample+".sorted.bam-unmapped.bam" for sample in config['ordered_sample_list']],


rule CDR3_TRUST:
    """perform CDR3 analysis using TRUST2"""
    input:
        bam="analysis/STAR/{sample}/{sample}.sorted.bam",
    output:
        fa="analysis/cdr3/{sample}/{sample}.sorted.bam.fa",
        loc="analysis/cdr3/{sample}/{sample}.sorted.bam-Locs.bam",
        umb="analysis/cdr3/{sample}/{sample}.sorted.bam-unmapped.bam",
        log="analysis/cdr3/{sample}/{sample}_TRUST.log",
    params:
        outdir="analysis/cdr3/{sample}/",
        pypath="PYTHONPATH=%s" % config["python2_pythonpath"]
    message:
        "CDR3: Performing CDR3 analysis using TRUST2"
    shell:
        "{params.pypath} {config[python2]} viper/modules/scripts/TRUST2.py -f {input.bam} -o {params.outdir} -a > {output.log}"
