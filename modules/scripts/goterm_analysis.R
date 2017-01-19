suppressMessages(library("AnnotationDbi"))
suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library("org.Mm.eg.db"))
suppressMessages(library("GOstats"))
suppressMessages(library("edgeR"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))
suppressMessages(library("ggalt"))
suppressMessages(library("mutoss"))

## The traceback is actually necessary to not break pipe at the stop step, so leave on
options(error = function() traceback(2))

## Example inputs
#deseq_file = "/mnt/cfce-stor1/home/mgc31/code/viperproject/analysis/test_all/diffexp/BCellvsTCell/BCellvsTCell.deseq.csv"
#adjpvalcutoff = 0.1
#numgoterms = 5
#reference = "hg19"
#up_goterm_csv = "/mnt/cfce-stor1/home/mgc31/code/viperproject/uptest.csv"
#down_goterm_csv = "/mnt/cfce-stor1/home/mgc31/code/viperproject/downtest.csv"
#goterm_pdf = "/mnt/cfce-stor1/home/mgc31/code/viperproject/test.pdf"
#up_goterm_png = "/mnt/cfce-stor1/home/mgc31/code/viperproject/up_test.png"
#down_goterm_png = "/mnt/cfce-stor1/home/mgc31/code/viperproject/down_test.png"

goterm_analysis_f <- function(deseq_file, adjpvalcutoff,numgoterms,reference, up_goterm_csv,down_goterm_csv,goterm_pdf,up_goterm_png,down_goterm_png) {

    adjpvalcutoff = as.numeric(adjpvalcutoff)
    numgoterms = as.numeric(numgoterms)

    ## Read in detable
    detable = read.table(deseq_file, header=TRUE, sep=",", fill=TRUE)
    rownames(detable) <- detable[,1]

    ## Append ENTREZ IDs from loaded in db
    if (reference == "hg19") {IDdb = org.Hs.eg.db}
    if (reference == "mm9") {IDdb = org.Mm.eg.db}
    detable$entrez <- mapIds(IDdb,
                         keys=rownames(detable),
                         column="ENTREZID",
                         keytype="SYMBOL",
                         multiVals="first")

    ## Select genes that pass the adjPval cutoff and select those entrez IDs as pop, set rest as universe.
    uptopgenes <- subset(detable, detable$padj < adjpvalcutoff & detable$log2FoldChange > 1)
    downtopgenes <- subset(detable, detable$padj < adjpvalcutoff & detable$log2FoldChange < -1)

    ## Failsafes to quit program if there are no differentially expressed genes that pass log2fc and adjpval
    if(nrow(uptopgenes) < 10) {stop(paste("Not enough significant genes for GOterm analysis, Need at least 10 and there are only ", nrow(uptopgenes)," upregulated genes at the current adjpval of ", adjpvalcutoff, sep=""))}
    if(nrow(downtopgenes) < 10) {stop(paste("Not enough significant genes for GOterm analysis, Need at least 10 and there are only ", nrow(downtopgenes)," downregulated genes at the current adjpval of ", adjpvalcutoff, sep=""))}
    
    upselectedIDs = na.omit(uptopgenes$entrez)
    downselectedIDs = na.omit(downtopgenes$entrez)
    universeIDs = na.omit(detable$entrez)

    ## Select database
    if (reference == "hg19") {annotDB = "org.Hs.eg.db"}
    if (reference == "mm9") {annotDB = "org.Mm.eg.db"}

    ## Run GOstats for upgenes
    upgoParams <- new("GOHyperGParams", geneIds = upselectedIDs, universeGeneIds = universeIDs, annotation = annotDB, ontology = "BP", pvalueCutoff = 1, conditional = FALSE, testDirection = "over")
    upgoResults <- hyperGTest(upgoParams)

    ## Summary table has columns: GOBPID, Pvalue, OddsRatio, ExpCount, Count, Size, Term.
    updf = summary(upgoResults)
    upfdrvalues = p.adjust(updf$Pvalue, method="fdr", n=length(updf$Pvalue))
    updf$adjPvalue = upfdrvalues
    updf$logadjpval = -log(updf$adjPvalue)

    ## Run GOstats for upgenes
    downgoParams <- new("GOHyperGParams", geneIds = downselectedIDs, universeGeneIds = universeIDs, annotation = annotDB, ontology = "BP", pvalueCutoff = 1, conditional = FALSE, testDirection = "over")
    downgoResults <- hyperGTest(downgoParams)

    ## Summary table has columns: GOBPID, Pvalue, OddsRatio, ExpCount, Count, Size, Term.
    downdf = summary(downgoResults)
    downfdrvalues = p.adjust(downdf$Pvalue, method="fdr", n=length(downdf$Pvalue))
    downdf$adjPvalue = downfdrvalues
    downdf$logadjpval = -log(downdf$adjPvalue)

    ## Write out Results for up
    xx = gsub(",","", as.matrix(updf[,7]))
    updf[,7] = xx
    write.table(updf, file = up_goterm_csv, col.names=T, row.names=F, quote=F, sep=",")

     ## Write out Results for down
    xx = gsub(",","", as.matrix(downdf[,7]))
    downdf[,7] = xx
    write.table(downdf, file = down_goterm_csv, col.names=T, row.names=F, quote=F, sep=",")
    
    ## Create title for plot
    temptitle = tail(unlist(strsplit(goterm_pdf, split="/")), n=1)
    temptitle = head(unlist(strsplit(temptitle, split="[.]")), n=1)
    uptitle = paste(temptitle, "_Top_", numgoterms, "_GOterms_UP", sep="")
    downtitle = paste(temptitle, "_Top_", numgoterms, "_GOterms_DOWN", sep="")

    ####### NEW
    plot_list <- list()
    ## Create Up Plot
    go.term.width = 20
    upgo.df <- as.data.frame(updf[numgoterms:1,], stringsAsFactors=FALSE)
    to_fold <- which(lapply(as.character((upgo.df$Term)), nchar) > go.term.width)

    foldednames <- sapply(upgo.df$Term[to_fold],function(x) paste(strwrap(x,width = go.term.width),collapse = "\n"))
    levels(upgo.df$Term) <- c(levels(upgo.df$Term), foldednames)
    upgo.df[to_fold,1] <- foldednames

    ggup <- ggplot(upgo.df, aes(y=reorder(Term, adjpval), x=logadjpval)) +
        geom_segment(aes(y=reorder(Term, logadjpval), yend=reorder(Term, logadjpval), x=0, xend=logadjpval)) +
        geom_point(aes(y=reorder(Term, logadjpval), x=logadjpval), color = "steelblue", size=5) +
        labs(x="- Log(Q-value)", y=NULL, title=uptitle) +
        theme(panel.grid.major.y=element_blank()) +
        theme(panel.grid.minor=element_blank()) +
        theme(axis.text.y=element_text(size=10)) +
        theme(plot.title = element_text(size=8, margin = margin(10, 0, 20, 0)))
    plot_list[[1]] = ggup


    ## Create Down Plot
    go.term.width = 20
    downgo.df <- as.data.frame(downdf[numgoterms:1,], stringsAsFactors=FALSE)
    to_fold <- which(lapply(as.character((downgo.df$Term)), nchar) > go.term.width)

    foldednames <- sapply(downgo.df$Term[to_fold],function(x) paste(strwrap(x,width = go.term.width),collapse = "\n"))
    levels(downgo.df$Term) <- c(levels(downgo.df$Term), foldednames)
    downgo.df[to_fold,1] <- foldednames

    ggdown <- ggplot(downgo.df, aes(y=reorder(Term, logadjpval), x=adjpval)) +
        geom_segment(aes(y=reorder(Term, logadjpval), yend=reorder(Term, logadjpval), x=0, xend=logadjpval)) +
        geom_point(aes(y=reorder(Term, logadjpval), x=logadjpval), color = "steelblue", size=5) +
        labs(x="- Log(Q-value)", y=NULL, title=downtitle) +
        theme(panel.grid.major.y=element_blank()) +
        theme(panel.grid.minor=element_blank()) +
        theme(axis.text.y=element_text(size=10)) +
        theme(plot.title = element_text(size=8, margin = margin(10, 0, 20, 0)))
    plot_list[[2]] = ggdown



    
    ## Plot out results
    #plot_list = list()
    #upgobarplot <- ggplot(updf[numgoterms:1,], aes(factor(Term, levels=unique(Term)), logadjpval)) +
    #  ylab("-log(adjPvalue)") + xlab("Go term") +
    #geom_bar(stat = "identity", fill="palegreen3") + theme_bw(base_size = 12) +
    #  scale_x_discrete(labels = function(x) str_wrap(x, width = 40, indent = 2),"\n") + coord_flip() + ggtitle(uptitle)
    #plot_list[[1]] = upgobarplot

    #downgobarplot <- ggplot(downdf[numgoterms:1,], aes(factor(Term, levels=unique(Term)), logadjpval)) +
    #  ylab("-log(adjPvalue)") + xlab("Go term") +
    #geom_bar(stat = "identity", fill="palegreen3") + theme_bw(base_size = 12) +
    #  scale_x_discrete(labels = function(x) str_wrap(x, width = 40, indent = 2),"\n") + coord_flip() + ggtitle(downtitle)
    #plot_list[[2]] = downgobarplot
    
    pdf(goterm_pdf)
    for (i in 1:2) { print(plot_list[[i]]) }
    junk <- dev.off()

    png(up_goterm_png, width = 8, height = 8, unit="in",res=300)
    print(plot_list[[1]])
    junk <- dev.off()

    png(down_goterm_png, width = 8, height = 8, unit="in",res=300)
    print(plot_list[[2]])
    junk <- dev.off()

}

args <- commandArgs( trailingOnly = TRUE )
deseq_file = args[1]
adjpvalcutoff = args[2]
numgoterms = args[3]
reference = args[4]
up_goterm_csv = args[5]
down_goterm_csv = args[6]
goterm_pdf = args[7]
up_goterm_png = args[8]
down_goterm_png = args[9]

goterm_analysis_f(deseq_file, adjpvalcutoff,numgoterms,reference, up_goterm_csv,down_goterm_csv,goterm_pdf,up_goterm_png,down_goterm_png)
