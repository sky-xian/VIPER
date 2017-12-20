## Load required packages
suppressMessages(library("gplots"))
suppressWarnings(suppressMessages(library("ComplexHeatmap")))
suppressMessages(library("circlize"))
suppressMessages(library("viridis"))
suppressMessages(library('dplyr'))
suppressMessages(source('viper/modules/scripts/supp_fns.R'))

## Enable stack trace
#options(error = function() traceback(2))

heatmapSS_plot <- function(rpkmTable,annot, ss_out_dir) {
    
    ## Read in and Log Transform Data
    Exp_data <- log2(rpkmTable+1)
    #CHECK: DROP cols that are all 0s
    Exp_data <- Exp_data[, -(which(colSums(Exp_data) == 0))]

    ## Calc. spearman correlation
    cordata <- cor(Exp_data, method="spearman")

    # NOTES on clustering, not used for now
    # Distance options: euclidean (default), maximum, canberra, binary, minkowski, manhattan
    # Cluster options: complete (default), single, average, mcquitty, median, centroid, ward
    rowdistance = dist(as.matrix(cordata), method = "euclidean")
    rowcluster = hclust(rowdistance, method = "ward.D2")
    coldistance = dist(t(as.matrix(cordata)), method = "euclidean")
    colcluster = hclust(coldistance, method = "ward.D2")
    
    ## make SS (sample-sample) heatmap
    ma_nolym <- max(cordata)
    mi_nolym <- min(cordata)
    my.breaks_nolym<-c(mi_nolym,seq(mi_nolym + 0.01, ma_nolym - 0.01,length.out=99),ma_nolym)
    
    pdf(file = paste(ss_out_dir, "heatmapSS_plot.pdf", sep=""))

    ha1 <- make_complexHeatmap_annotation(annot)

    mapplot <-Heatmap(t(as.matrix(cordata)),
                     col = colorRamp2(my.breaks_nolym,  bluered(101), transparency = 0),
                     #column_dend_height = unit(2, "cm"),
                     #heatmap_legend_param = list(title = "exp. level"),
                     column_title = "Sample-Sample Correlation",
                     #row_title = "Samples",
                     show_row_names = TRUE, show_column_names = TRUE,
                     #row_names_max_width = unit(3, "mm"),
                     row_names_gp = gpar(fontsize = 6),
                     column_names_gp = gpar(fontsize = 6),
                     #cluster_rows = TRUE,
                     #cluster_columns=TRUE,
                     cluster_rows = rowcluster,
                     cluster_columns = colcluster,
                     show_heatmap_legend = TRUE,
                     heatmap_legend_param=list(title="corr", title_gp=gpar(fontsize=8), labels_gp=gpar(fontsize=8)),
                     #row_dend_width = unit(5, "mm"),
                     #width=unit(60,"cm"),
                     top_annotation=ha1,
                     )
    draw(mapplot)
    for(an in colnames(annot[1:ncol(annot)])) {
        decorate_annotation(an, {
            grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp=gpar(fontsize=6), check=TRUE)
            grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right", gp=gpar(fontsize=6), check=TRUE)
        })
    }
    dev.off()
    
    png(file=paste(ss_out_dir, "images/heatmapSS_plot.png", sep=""), width = 8, height = 8, unit="in",res=300)
    draw(mapplot)
    for(an in colnames(annot[1:ncol(annot)])) {
        decorate_annotation(an, {
            grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp=gpar(fontsize=6), check=TRUE)
            grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right", gp=gpar(fontsize=6), check=TRUE)
        })
    }
    dev.off()

    #WRITE output to file
    output<-as.matrix(cordata)
    output<-output[rowcluster$order, colcluster$order]
    write.table(output, file=paste(ss_out_dir, "heatmapSS.txt",sep=""), quote=F, col.names = NA, sep="\t")
    
}


args <- commandArgs( trailingOnly = TRUE )
rpkmFile=args[1]
annotFile=args[2]
ss_out_dir=args[3]

rpkmTable <- read.csv(rpkmFile, header=T, check.names=F, row.names=1, stringsAsFactors=FALSE, dec='.')

annot <- read.csv(annotFile, sep=",", header=T, row.names=1, stringsAsFactors=FALSE, check.names=F, comment.char='#')
if(any(grepl("comp_*", colnames(annot)))) {
  annot <- annot[, !grepl('Pair', colnames(annot)), drop = F]
  annot <- annot[, !grepl('comp_*', colnames(annot)), drop = F]
}

heatmapSS_plot(rpkmTable,annot, ss_out_dir)
