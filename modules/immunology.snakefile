rule immunology_all:
    input:
        "analysis/immunology/relative_abundance.txt",
        "analysis/immunology/output.pdf",
        "analysis/immunology/TIMER_results.pdf",

rule estimate_immune_abundance:
    input:
        #fpkm_collected="analysis/cufflinks/Cuff_FPKM_Collected.tsv"
        fpkm_collected="analysis/default/cufflinks/Cuff_Gene_Counts.csv"
    output:
        "analysis/immunology/relative_abundance.txt",
        "analysis/immunology/output.pdf",
        "analysis/immunology/TIMER_results.pdf",
    message: "Estimating immune cell abundance output"
    run:
        shell( "Rscript viper/modules/scripts/immunology.R {input.fpkm_collected} {cancer_type} --staticdir=viper/static/immunology --outdir=`pwd`/analysis/immunology/", cancer_type=config["cancer_type"])
