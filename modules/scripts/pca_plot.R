#!/usr/bin/env Rscript
#-------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: May, 23, 2016
#--------------------

### Added these scripts to hardcoded supp_fns to avoid online communication issues
#if( is.element("ggbiplot", installed.packages())){
#      #suppressMessages(library(ggbiplot))
#    library(ggbiplot)
#} else {
#      #suppressMessages(require("devtools"))
#    require("devtools")
#      install_github("vangalamaheshh/ggbiplot")
#      #suppressMessages(require(ggbiplot))
#    require(ggbiplot)
#}

options(error = function() traceback(2))

library(ggplot2)
library(plyr)
library(scales)
library(grid)
library(ggrepel)
suppressMessages(source('viper/modules/scripts/supp_fns.R'))

#source("/mnt/cfce-stor1/home/mgc31/code/viperproject/viper/modules/scripts/supp_fns.R")
#rpkmFile = "/mnt/cfce-stor1/home/mgc31/code/viperproject/analysis/test2_all/cufflinks/Cuff_Gene_Counts.filtered.csv"
#metaFile = "/mnt/cfce-stor1/home/mgc31/code/viperproject/weinstock_test_metasheet.csv"
#pca_out_dir = "/mnt/cfce-stor1/home/mgc31/code/viperproject/analysis/test2_all/plots/"

pca_plot <- function(rpkmTable, annot, pca_out_dir) {
  rpkm.pca <- prcomp(t(rpkmTable), center = TRUE, scale. = TRUE)
  plot.var <- ggscreeplot(rpkm.pca)
  suppressMessages(ggsave(paste(pca_out_dir,"images/pca_plot_scree.png", sep="")))
  all_plots <- list()
  for (ann in colnames(annot)){
    g <- ggbiplot(rpkm.pca, groups = as.character(annot[,ann]), scale = 1, var.scale = 1, obs.scale = 1,
                labels=colnames(rpkmTable), choices = 1:2,
                ellipse=FALSE, circle = TRUE, var.axes = FALSE)
    g <- g + scale_color_discrete(name = ann)
    g <- g + theme(legend.direction = 'horizontal',
                   legend.position = 'top',
                   legend.title = element_text(face="bold"))
    all_plots <- c(all_plots, list(g))
    suppressMessages(ggsave(paste(pca_out_dir, "images/pca_plot_", ann, ".png", sep="")))
  }

  pdf(paste(pca_out_dir, "pca_plot.pdf", sep=""))
  capture.output(print(c(all_plots,list(plot.var))))
  junk <- dev.off()
}


args <- commandArgs( trailingOnly = TRUE )
rpkmFile <- args[1]
metaFile <- args[2]
pca_out_dir <- args[3]

rpkmTable <- read.csv(rpkmFile, header=T, check.names=F,
                        row.names=1, stringsAsFactors=FALSE, dec='.')
annot <- read.csv(metaFile, sep=",", header=T, row.names=1,
                      stringsAsFactors=FALSE, check.names=F, comment.char='#')
if(any(grepl("comp_*", colnames(annot)))) {
  annot <- annot[, !grepl('comp_*', colnames(annot)), drop = F]
}
pca_plot(rpkmTable, annot, pca_out_dir)
