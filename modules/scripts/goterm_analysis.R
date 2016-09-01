suppressMessages(library("AnnotationDbi"))
suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library("org.Mm.eg.db"))
suppressMessages(library("GOstats"))
suppressMessages(library("edgeR"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))
suppressMessages(library("ggalt"))
#suppressMessages(library("scales"))
suppressMessages(library("mutoss"))

## The traceback is actually necessary to not break pipe at the stop step, so leave on  
options(error = function() traceback(2))

goterm_analysis_f <- function(deseq_file, adjpvalcutoff,numgoterms,reference, goterm_csv,goterm_pdf,goterm_png) {

    adjpvalcutoff = as.numeric(adjpvalcutoff)
    numgoterms = as.numeric(numgoterms)
    
    ## Read in detable
    detable = read.table(deseq_file, header=TRUE, sep=",", fill=TRUE)
    rownames(detable) <- detable[,1]
    
    ## Append ENSEMBL and ENTREZ IDs from loaded in db
    if (reference == "hg19") {IDdb = org.Hs.eg.db}
    if (reference == "mm9") {IDdb = org.Mm.eg.db}
    
    detable$ensembl <- mapIds(IDdb,
                         keys=rownames(detable),
                         column="ENSEMBL",
                         keytype="SYMBOL",
                         multiVals="first")
    detable$entrez <- mapIds(IDdb,
                         keys=rownames(detable),
                         column="ENTREZID",
                         keytype="SYMBOL",
                         multiVals="first")

    ## Select genes that pass the adjPval cutoff and select those entrez IDs as pop, set rest as universe.
    topgenes <- subset(detable, detable$padj < adjpvalcutoff)

    if(nrow(topgenes) < numgoterms) {stop(paste("Not enough significant genes, there are only ", nrow(topgenes)," genes at the current pval of ", adjpvalcutoff, sep=""))}
    
    selectedIDs = topgenes$entrez
    universeIDs = detable$entrez

    if (reference == "hg19") {annotDB = "org.Hs.eg.db"}
    if (reference == "mm9") {annotDB = "org.Mm.eg.db"}
    
    ## Run GOstats
    goParams <- new("GOHyperGParams",
                    geneIds = selectedIDs,
                    universeGeneIds = universeIDs,
                    annotation = annotDB,
                    ontology = "BP", ## CC, BP, MF
                    pvalueCutoff = adjpvalcutoff,
                    conditional = TRUE,
                    testDirection = "over")

    goResults <- hyperGTest(goParams)

    ## Summary table has columns: GOBPID, Pvalue, OddsRatio, ExpCount, Count, Size, Term.
    df = summary(goResults)

    by <- BY(df$Pvalue, 0.05)
    df$adjPvalue <- by[["adjPValues"]]
    df$logpval = -log(df$Pvalue)

    ## Write out Results
    xx = gsub(",","", as.matrix(df[,7]))
    df[,7] = xx    
    write.table(df, file = goterm_csv, col.names=T, row.names=F, quote=F, sep=",")
    
    ## Create title for plot
    temptitle = tail(unlist(strsplit(goterm_pdf, split="/")), n=1)
    temptitle = head(unlist(strsplit(temptitle, split="[.]")), n=1)
    title = paste(temptitle, "_Top_", numgoterms, "_GOterms", sep="")

    go_bar_plot <- ggplot(df[numgoterms:1,], aes(factor(Term, levels=unique(Term)), logpval)) + 
      ylab("-log(Pvalue)") + xlab("Go term") +
      geom_bar(stat = "identity", fill="palegreen3") + theme_bw(base_size = 12) +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 40, indent = 2),"\n") + coord_flip() + ggtitle(title)

    ggsave(goterm_pdf, width=11, height=8.5, unit="in")
    ggsave(goterm_png, width=10, height=8, unit="in")

}

args <- commandArgs( trailingOnly = TRUE )
deseq_file = args[1]
adjpvalcutoff = args[2]
numgoterms = args[3]
reference = args[4]
goterm_csv = args[5]
goterm_pdf = args[6]
goterm_png = args[7]

goterm_analysis_f(deseq_file, adjpvalcutoff,numgoterms,reference, goterm_csv,goterm_pdf,goterm_png)



#goterm_anlysis_f(
#    snakemake@input[[""]],
#    snakemake@params[[""]],
#    snakemake@output[[""]]
#    )
