#!/usr/bin/env Rscript
#-------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: May, 23, 2016
#--------------------

library(dplyr)
if( is.element("ggbiplot", installed.packages())){
  library(ggbiplot)
} else {
  require("devtools")
  install_github("vqv/ggbiplot")
  require(ggbiplot)
}

options(error = function() traceback(2))

preprocess <- function(rpkm_file, metasheet, filter_miRNA=TRUE, 
                       min_genes=250, min_samples=4, rpkm_cutoff=2.0) {
  rpkmTable <- read.csv(rpkm_file, header=T, check.names=F, 
                        row.names=1, stringsAsFactors=FALSE, dec='.')
  for (n in names(rpkmTable)) {
    rpkmTable[n] <- apply(rpkmTable[n], 1, as.numeric)
  }
  rpkmTable <- na.omit(rpkmTable)
  tmp_ann <- read.csv(metasheet, sep=",", header=T, row.names=1, 
                      stringsAsFactors=FALSE, check.names=F)
  tmp_ann <- dplyr::select(tmp_ann, -(starts_with("comp_")))
  df <- dplyr::select_(rpkmTable, .dots=rownames(tmp_ann))
  sub_df <- df[apply(df, 1, function(x) length(x[x>=rpkm_cutoff])>min_samples),]
  sub_df <- log2(sub_df + 1)
  if (filter_miRNA == TRUE) {
    sub_df <- sub_df[ !grepl("MIR|SNO",rownames(sub_df)), ]
  }
  min_genes = min(min_genes, nrow(sub_df))
  ## Calculate CVs for all genes (rows)
  mean_rpkm <- apply(sub_df,1,mean)
  var_rpkm <- apply(sub_df,1,var)
  cv_rpkm <- abs(var_rpkm/mean_rpkm)
  ## Select out the most highly variable genes into the dataframe 'Exp_data'
  exp_data <- sub_df[order(cv_rpkm,decreasing=T)[1:min_genes],]
  return (list(exp_data=exp_data, tmp_ann=tmp_ann))
}

pca_plot <- function(rpkmTable, annot, pca_plot_out) {
  rpkm.pca <- prcomp(t(rpkmTable), center = TRUE, scale. = TRUE)
  pc_var <- signif(100.0 * summary(rpkm.pca)[[6]][2,1:3], digits = 3)
  pc_var <- data.frame(PCA=names(pc_var), Variance=pc_var)
  plot.var <- ggplot(pc_var, aes(x=PCA,y=Variance))
  plot.var <- plot.var + geom_bar(stat="identity") + theme_bw() 
  plot.var <- plot.var + ylab("% Variance") + xlab("PCA")
  ggsave("analysis/plots/images/pca_plot_scree.png")
  all_plots <- list()
  for (ann in colnames(annot)){
    g <- ggbiplot(rpkm.pca, groups = annot[,ann], scale = 0, var.scale = 0,
                labels=colnames(rpkmTable), ellipse = TRUE,
                labels.size=3, circle = TRUE, var.axes = FALSE)
    g <- g + scale_color_discrete(name = '')
    g <- g + theme(legend.direction = 'horizontal',
                 legend.position = 'top')
    all_plots <- c(all_plots, list(g))
    ggsave(paste("analysis/plots/images/pca_plot_",ann,".png",sep=""))
  }
  pdf(pca_plot_out)
  print(c(all_plots,list(plot.var)))
  dev.off()
}


args <- commandArgs( trailingOnly = TRUE )
rpkmFile <- args[1]
annotFile <- args[2]
rpkm_cutoff <- args[3]
min_samples <- args[4]
filter_miRNA <- args[5]
min_genes <- args[6]
pca_plot_out <- args[7]


info <- preprocess(rpkmFile, annotFile, filter_miRNA, min_genes, 
                       min_samples, rpkm_cutoff)
pca_plot(info$exp_data, info$tmp_ann, pca_plot_out)
