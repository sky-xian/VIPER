#!/usr/bin/env Rscript
# vim: syntax=r tabstop=4 expandtab

#---------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Aug, 25, 2016
#---------------------------------

#suppressMessages(library(argparse))
if( is.element("Seurat", installed.packages())){
    suppressMessages(library(Seurat))
} else {
    suppressMessages(require("devtools"))
    source('http://bioconductor.org/biocLite.R')
    install_github("Storeylab/lfa")
    install_github("satijalab/seurat")
    suppressMessages(require(Seurat))
}

options(error = function() traceback(2))

runSeurat <- function( matrix_file ){
    sc.data <- read.csv( matrix_file, header=TRUE, row.names=1 )
    sc.data <- log( sc.data + 1 )
    sc <- new( "seurat", raw.data=sc.data )
    sc <- Setup( sc, project="Viper Single Cell Analysis", min.cells = 3, names.field = 1, names.delim = "_", min.genes = 1000, is.expr=1, )
    #png( "analysis/seurat/mean_var_plot.png", width = 8, height = 8, unit="in",res=300 )
    sc <- MeanVarPlot( sc, y.cutoff = 2, x.low.cutoff = 1, fxn.x = expMean, fxn.y = logVarDivMean )
    #dev.off()
    #png( "analysis/seurat/linear_PCA.png", width = 8, height = 8, unit="in",res=300)
    sc <- PCA( sc, do.print=FALSE )
    PCAPlot( sc, 1, 2, pt.size = 2 )
    #dev.off()
    
    #png( "analysis/seurat/non_linear_tSNE_PCA.png", width = 8, height = 8, unit="in",res=300 )
    sc <- RunTSNE( sc, dims.use = 1:2, max_iter=2000 )
    TSNEPlot( sc )
    #dev.off()
}

parse_args <- function() {
    parser <- argparse::ArgumentParser(description="Takes TPM matrix file and runs seurat tSNE methods")
    parser$add_argument('-m', '--matrix_file', help="TPM matrix file", type="character", nargs=1)
    args <- parser$parse_args()
    return (args)
}

#args <- parse_args()
args <- commandArgs( trailingOnly = TRUE )
runSeurat( args[1] )
