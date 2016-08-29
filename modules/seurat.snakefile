#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab

#-------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Aug, 25, 2016
#-------------------------

rule run_seurat:
    input:
        gene_matrix = "analysis/RSEM/tpm_gene_matrix.csv"
    output:
        seurat_out = "analysis/seurat/seurat.done"
    message: "Running Seurat - tSNE"
    shell:
        "{config[seurat_path]}/Rscript --default-packages=methods,utils viper/modules/scripts/run_seurat.R {input.gene_matrix} "
        "&& touch {output.seurat_out}"
