rule immunology_all:
    input:
        "analysis/"+config["token"]+"/immunology/relative_abundance.txt",
        "analysis/"+config["token"]+"/immunology/output.pdf",
        "analysis/"+config["token"]+"/immunology/TIMER_results.pdf",

rule estimate_immune_abundance:
    input:
        fpkm_collected="analysis/"+config["token"]+"/cufflinks/Cuff_Gene_Counts.csv"
    output:
        "analysis/"+config["token"]+"/immunology/relative_abundance.txt",
        "analysis/"+config["token"]+"/immunology/output.pdf",
        "analysis/"+config["token"]+"/immunology/TIMER_results.pdf",
    message: "Estimating immune cell abundance output"
    run:
        shell( "Rscript viper/modules/scripts/immunology.R {input.fpkm_collected} {cancer_type} --staticdir=viper/static/immunology --outdir=`pwd`/analysis/immunology/", cancer_type=config["cancer_type"])
