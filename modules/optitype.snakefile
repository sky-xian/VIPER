#!/usr/bin/env python

rule optitype_get_chr6_reads:
    input:
        bam="analysis/STAR/{sample}/{sample}.sorted.bam",
        bai="analysis/STAR/{sample}/{sample}.sorted.bam.bai"
    output:
        "analysis/STAR/{sample}/{sample}.sorted.chr6.bam"
    message: "OPTITYPE: obtaining chr6 reads"
    shell:
        "samtools view -b -h {input.bam} chr6 > {output}"

rule optitype_convert_bam2fq:
    input:
        bam="analysis/STAR/{sample}/{sample}.sorted.chr6.bam"
    output:
        #USING _mates from align.snakefile
        mates = expand( "analysis/STAR/{{sample}}/{{sample}}.chr6.{mate}.fq", mate=_mates)
    message: "OPTITYPE: converting chr6 bam to fastq"
    run:
        out_files = " ".join(["-%s %s" % (i+1,m) for (i,m) in enumerate(output.mates)])
        #out_files= [m for (i,m) in enumerate(output.mates)]
        shell("samtools bam2fq {out_files} {input}")

rule optitype:
    input:
        mates = expand( "analysis/STAR/{{sample}}/{{sample}}.chr6.{mate}.fq", mate=_mates)
    output:
        results="analysis/optitype/{sample}/{sample}_result.tsv",
        coverage="analysis/optitype/{sample}/{sample}_coverage_plot.pdf"
    message: "OPTITYPE: running optitype"
    params:
        pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
        outpath="analysis/optitype/{sample}",
        name="{sample}"
    shell:
        "{params.pypath} {config[python2]} {config[optitype_path]}/OptiTypePipeline.py -i {input} -r -o {params.outpath} -c viper/static/optitype/config.ini -p {params.name}"
