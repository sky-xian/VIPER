suppressMessages(library("dplyr"))
suppressMessages(library("AnnotationDbi"))
suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library("gage"))
suppressMessages(library("gageData"))
suppressMessages(library("pathview"))
suppressMessages(library("clusterProfiler"))
suppressMessages(library("XML"))
#data(kegg.sets.hs)
#data(sigmet.idx.hs)

## "cleans up the gene sets"
#kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
kegg.sets.hs = kegg.gsets()
kegg.sets.hs = kegg.sets.hs$kg.sets

kegg.set = kegg.gsets()
ks = kegg.set$kg.sets
kss = kegg.set$kg.sets[kegg.set$sigmet.idx]
names(kss) = gsub("/","",names(kss))

## Removing pathways that I know don't load properly... no idea why
kss[which(names(kss) %in% c("hsa01200 Carbon metabolism"))] <- NULL

kegg_pathway_f<- function(deseq_file, keggpvalcutoff,numkeggpathways,kegg_dir,reference,temp_dir, kegg_table_up,kegg_table_down,gsea_table,gsea_pdf) {

    ## These are here until we update snakemake
    numkeggpathways = as.numeric(numkeggpathways)
    keggpvalcutoff = as.numeric(keggpvalcutoff)

    ## Will need this path stuff for later as kegg output is very messy
    mainDir = substr(kegg_dir, 1, nchar(kegg_dir)-14)
    dir.create(file.path(mainDir, "kegg_pathways/"), showWarnings = FALSE)

    ## Read in deseq table
    detable = read.table(deseq_file, header=TRUE, sep=",", fill=TRUE)
    rownames(detable) <- detable[,1]

    ## Append ENSEMBL and ENTREZ IDs from loaded in db
    if (reference == "hg19") {IDdb = org.Hs.eg.db}
    #if (reference == "mm9") {IDdb = ##########}
        
    detable$ensembl = mapIds(IDdb,
                        keys=row.names(detable),
                        column="ENSEMBL",
                        keytype="SYMBOL",
                        multiVals="first")
    detable$entrez = mapIds(IDdb,
                        keys=row.names(detable),
                        column="ENTREZID",
                        keytype="SYMBOL",
                        multiVals="first")
    detable$description = mapIds(IDdb,
                        keys=row.names(detable),
                        column="GENENAME",
                        keytype="SYMBOL",
                        multiVals="first")

    ## Couple failsafes
    detable = na.omit(detable)
    detable = detable[is.finite(detable$log2FoldChange),]
    
    ## Setting up gage input, needs the log2fc with the entrez id
    gageinput = detable$log2FoldChange
    names(gageinput) = detable$entrez

    ## Run gage
    keggres = gage(gageinput, gsets = kss, same.dir=TRUE)

    kegg_up = keggres$greater
    kegg_up = cbind(rownames(kegg_up), kegg_up)
    colnames(kegg_up)[1] = "Kegg_pathway"
    xx = gsub(",","", as.matrix(kegg_up[,1]))
    kegg_up[,1] = xx
    write.table(kegg_up, file = kegg_table_up, quote=F, col.names=TRUE, row.names=FALSE, sep=",")

    kegg_down= keggres$less
    kegg_down= cbind(rownames(kegg_down), kegg_down)
    colnames(kegg_down)[1] = "Kegg_pathway"
    xx = gsub(",","", as.matrix(kegg_down[,1]))
    kegg_down[,1] = xx
    write.table(kegg_down, file = kegg_table_down, quote=F, col.names=TRUE, row.names=FALSE, sep=",")

    
    ## Stop run if the params don't match for output
    kegg_output_filter = subset(kegg_up, kegg_up[,4] < keggpvalcutoff)
    if(nrow(kegg_output_filter) < numkeggpathways) {stop(paste("Only ",nrow(kegg_output_filter), " pathways pass the current keggpvalcutoff of ", keggpvalcutoff, ", please run again with increased pval. Check comp.kegg.txt for details", sep="")) }

    
    ## Get the pathways
    keggrespathways = keggres$stats
    keggrespathways = keggrespathways[order(-abs(keggrespathways[,1])),]
    keggrespathways = rownames(keggrespathways)[1:numkeggpathways]
    keggresids = substr(keggrespathways, start=1, stop=8)

    ## Plot using pathview
    
    normwd = getwd()
    setwd(temp_dir)
    
    for ( i in 1:numkeggpathways) {
        pvout <- pathview(gene.data=gageinput,              ## Gene list
                          pathway.id=keggresids[i],         ## Which pathway
                          species = "hsa",                  ## Species
                          #limit = list(gene=max(abs(gageinput)),cpd=1),
                          #kegg.dir = temp_dir               ## Save directory
                          low = list(gene = "blue", cpd = "purple"),
                          mid = list(gene = "gray", spd = "gray"),
                          high = list(gene = "red", cpd = "gray")  ## Color scale
                     )
    }
    setwd(normwd)

    ## Renaming files
    # Create variable with keggrespathways sorted and pull out name
    sortkeggrespathways = sort(keggrespathways)
    newnames = substr(sortkeggrespathways, 10, nchar(sortkeggrespathways))
    newnames = gsub(" ", "_", newnames)
    newnames = paste0(newnames, "_", match(sortkeggrespathways,keggrespathways))
    
    # Read in the list of made png files
    png_files <- list.files(temp_dir, pattern=glob2rx("*.pathview.png"))
    file.rename(paste0(temp_dir,png_files), paste0(kegg_dir, newnames, ".png"))

    # Repeat for xml files
    xml_files <- list.files(temp_dir, pattern=glob2rx("*.xml"))
    file.rename(paste0(temp_dir,xml_files), paste0(kegg_dir, newnames, ".xml"))

    
    ## GSEA Analysis
    gseainput = sort(gageinput, decreasing=TRUE)
    gseainput = gseainput[is.finite(gseainput)]

    fullgsea <- gseKEGG(geneList = gseainput,
                        organism     = "human",
                        nPerm        = 100,
                        minGSSize    = 1,
                        pvalueCutoff = 0.99,
                        verbose      = FALSE,
                        use_internal_data = FALSE)

    gsea_data = summary(fullgsea)
    gsea_data = gsea_data[order(-abs(gsea_data$NES)),]
    xx = gsub(",","", as.matrix(gsea_data[,2]))
    gsea_data[,2] = xx
    write.table(gsea_data, file = gsea_table, quote=FALSE, sep= ",", row.names=FALSE, col.names=TRUE)

    pdf(gsea_pdf)
    plot.new()
    mtext("GSEA_plots")
    for ( i in 1:10 ) {
        gseaplot(fullgsea, geneSetID = gsea_data[i,1])
        mtext(gsea_data[i,2])
    }
    dev.off()
    
    }

args <- commandArgs( trailingOnly = TRUE )
deseq_file = args[1]
keggpvalcutoff = args[2]
numkeggpathways = args[3]
kegg_dir = args[4]
reference = args[5]
temp_dir = args[6]
kegg_table_up = args[7]
kegg_table_down = args[8]
gsea_table = args[9]
gsea_pdf = args[10]


kegg_pathway_f(deseq_file, keggpvalcutoff,numkeggpathways,kegg_dir,reference,temp_dir, kegg_table_up,kegg_table_down,gsea_table,gsea_pdf)




#kegg_pathway_f(
#    snakemake@input[[""]],
#    snakemake@params[[""]],
#    snakemake@output[[""]]
#    )
