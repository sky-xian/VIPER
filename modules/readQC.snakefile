#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#----------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#----------------------------------

rule read_distrib_qc:
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        protected("analysis/RSeQC/read_distrib/{sample}.txt")
    message: "Running RseQC read distribution on {wildcards.sample}"
    params: pypath="PYTHONPATH=%s" % config["python2_pythonpath"]
    shell:
        "{params.pypath} {config[python2]} {config[rseqc_path]}/read_distribution.py"
        " --input-file={input}"
        " --refgene={config[bed_file]} 1>{output}"

rule read_distrib_qc_matrix:
    input:
        read_distrib_files=expand( "analysis/RSeQC/read_distrib/{sample}.txt", 
            sample=config["ordered_sample_list"] ),
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        matrix="analysis/RSeQC/read_distrib/read_distrib.matrix.tab",
        png="analysis/RSeQC/read_distrib/read_distrib.png"
    message: "Creating RseQC read distribution matrix"
    run:
        file_list_with_flag = " -f ".join( input.read_distrib_files )
        shell( "perl viper/modules/scripts/read_distrib_matrix.pl -f {file_list_with_flag} 1>{output.matrix}" )
        shell( "Rscript viper/modules/scripts/read_distrib.R {output.matrix} {output.png}" )


rule gene_body_cvg_qc:
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        protected("analysis/RSeQC/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.curves.png"),
        protected("analysis/RSeQC/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.r")
    message: "Creating gene body coverage curves"
    params: pypath="PYTHONPATH=%s" % config["python2_pythonpath"]
    shell:
        "{params.pypath} {config[python2]} {config[rseqc_path]}/geneBody_coverage.py -i {input} -r {config[bed_file]}"
        " -f png -o analysis/RSeQC/gene_body_cvg/{wildcards.sample}/{wildcards.sample}"


rule plot_gene_body_cvg:
    input:
        samples_list=expand("analysis/RSeQC/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.r", 
            sample=config["ordered_sample_list"] ),
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        rscript="analysis/RSeQC/gene_body_cvg/geneBodyCoverage.r",
        png="analysis/RSeQC/gene_body_cvg/geneBodyCoverage.heatMap.png",
        png_curves="analysis/RSeQC/gene_body_cvg/geneBodyCoverage.curves.png"
    message: "Plotting gene body coverage"
    shell:
        "perl viper/modules/scripts/plot_gene_body_cvg.pl --rfile {output.rscript} --png {output.png} --curves_png {output.png_curves}"
        " {input.samples_list} && Rscript {output.rscript}"

rule junction_saturation:
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        protected("analysis/RSeQC/junction_saturation/{sample}/{sample}.junctionSaturation_plot.pdf")
    message: "Determining junction saturation for {wildcards.sample}"
    params: pypath="PYTHONPATH=%s" % config["python2_pythonpath"]
    shell:
        "{params.pypath} {config[python2]} {config[rseqc_path]}/junction_saturation.py -i {input} -r {config[bed_file]}"
        " -o analysis/RSeQC/junction_saturation/{wildcards.sample}/{wildcards.sample}"


rule collect_insert_size:
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        protected("analysis/RSeQC/insert_size/{sample}/{sample}.histogram.pdf")
    message: "Collecting insert size for {wildcards.sample}"
    shell:
        "{config[picard_path]} CollectInsertSizeMetrics"
        " H={output} I={input} O=analysis/RSeQC/insert_size/{wildcards.sample}/{wildcards.sample} R={config[ref_fasta]}"



