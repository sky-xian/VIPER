suppressMessages(library("limma"))
suppressWarnings(suppressMessages(library("DESeq2")))
suppressMessages(library("edgeR"))
suppressMessages(library(tximport))
suppressMessages(library(readr))

limma_and_deseq_f <- function(arg_rsem, arg_names, arg_s1, arg_s2, limma, deseq, limma_annot, deseq_annot, deseqSum_out, gene_annotation) {
    #READ in gene_annotation table--even though gene descriptions are quoted
    #in the annotations, we have to quote it again!
    if( grepl(".bz2$",gene_annotation) ) {
        gene_annot <- read.table(bzfile(gene_annotation), header=TRUE, sep=",", quote="", fill=TRUE, row.names=NULL)
    }
    else {
        gene_annot <- read.table(gene_annotation, header=TRUE, sep=",", quote="", fill=TRUE, row.names=NULL)
    }

    gene_annot[,'Gene.Description'] <- sapply(gene_annot[,'Gene.Description'],
                                              function(x){paste0("\"",x,"\"")},
                                              simplify="vector")
    
    #Read in lists for comparison, read in rsem files using tximport,
    #and associate sample names
    treatlist = strsplit(arg_s2,',')[[1]]
    ctrllist = strsplit(arg_s1,',')[[1]]
    rsem_files = strsplit(arg_rsem,",")[[1]]
    names = strsplit(arg_names,",")[[1]]
    names(rsem_files) <- names

    #HERE we select out the sample/cols that are related to the analysis
    rsem_files <- rsem_files[names(rsem_files) %in% c(treatlist,ctrllist)]
    
    #how to use tximport
    #ref: http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
    txi.deseq <- tximport(rsem_files, type = "rsem")

    #MAKE the condition data frame
    ntreat <- length(treatlist)
    nctrl <- length(ctrllist)
    #condition is VECTOR of ntreat 'treat' followed by nctrl 'control'
    #exmple: TREAT TREAT TREAT CONTROL CONTROL CONTROL
    condition = c(rep('treat',ntreat),rep('control',nctrl))
    #now make a table that pairs each sampleName with it's condition
    colData <- as.data.frame(cbind(c(treatlist,ctrllist),condition))

    preparationD <- function (txi, sampleTable, condition){
        #IF there are transcript lengths == 0, set it to 1
        #ref: ref: https://support.bioconductor.org/p/84304/
        txi$length[txi$length == 0] <- 1
        
        dds <- DESeqDataSetFromTximport(txi, sampleTable, design= ~ condition)
        dds <- dds[rowSums(counts(dds)) > 0, ]
        dds <- DESeq(dds)
        res <- results(dds)
        
        #Get summary stats
        summary <- c(sum(res$padj<0.05 & res$log2FoldChange>0.0, na.rm=TRUE),
                     sum(res$padj<0.05 & res$log2FoldChange>0.5, na.rm=TRUE),
                     sum(res$padj<0.05 & res$log2FoldChange>1.0, na.rm=TRUE),
                     sum(res$padj<0.01 & res$log2FoldChange>0.0, na.rm=TRUE),
                     sum(res$padj<0.01 & res$log2FoldChange>0.5, na.rm=TRUE),
                     sum(res$padj<0.01 & res$log2FoldChange>1.0, na.rm=TRUE))
        sumTable <- matrix(summary, nrow=3, ncol=2)
        rownames(sumTable)<-c('log2FC > 0.0','log2FC > 0.5', 'log2FC > 1.0')
        colnames(sumTable)<-c('padj < 0.05','padj < 0.01')
        #Write/output summary stats
        write.table(sumTable, deseqSum_out, quote=FALSE, sep=",")
        return (res)
    }
    #CALL
    deseq_result <- preparationD(txi.deseq, colData, condition)
    
    ## Setting the first column name to 'id'
    deseq_result <- cbind(id=rownames(deseq_result), as.matrix(deseq_result))
    ## Sort by padj. and remove padj = NA
    deseq_result <-deseq_result[order(as.numeric(deseq_result[,"padj"])),]
    ## Write output deseq
    write.table(deseq_result,deseq,sep=',',col.names=T,row.names=F,quote=F)

    ## ANNOTATE w/ local .csv biomart annotation file
    deseq_annotations <- merge(deseq_result, gene_annot, by= "id", all.x = TRUE, sort=FALSE)
    ## Sort by padj. and remove padj = NA
    write.table(deseq_annotations,deseq_annot,sep=',',col.names=T,row.names=F,quote=F)
    
    if (ntreat>1){
        #STAR limma calls-
        #prepare data for limma import
        #ref: http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
        txi.limma <- tximport(rsem_files, type = "rsem",
                              countsFromAbundance = "lengthScaledTPM")
        #NORMALIZING FACTOR
        nf <- calcNormFactors(txi.limma$counts)
        design <- model.matrix(~condition)
        v <- voom(txi.limma$counts, design,
                  lib.size=colSums(txi.limma$counts)*nf)
        fit <- lmFit(v,design)
        fit <- eBayes(fit)
        #END limma calls
        
        limma_result=topTable(fit,coef=2,n=length(txi.limma$counts[,1]),sort="p")
        limma_result <- cbind(id=rownames(limma_result), limma_result)
        ## ANNOTATE w/ local .csv biomart annotation file
        limma_annotations <- merge(limma_result, gene_annot, by= "id", all.x = TRUE, sort=FALSE)
        ## MOVING to just csv files
        write.table(limma_result,limma_out,sep=',',col.names=T,row.names=F,quote=F)
        ## WRITE annotation table--limma.annot.csv
        write.table(limma_annotations,limma_annot,sep=',',col.names=T,row.names=F,quote=F)
    }
}

args <- commandArgs( trailingOnly = TRUE )
arg_rsem = args[1]
arg_names = args[2]
arg_s1 = args[3]
arg_s2 = args[4]
limma_out=args[5]
deseq_out=args[6]
limma_annot = args[7]
deseq_annot = args[8]
deseqSum_out=args[9]
gene_annotation=args[10]

limma_and_deseq_f(arg_rsem, arg_names, arg_s1, arg_s2, limma_out, deseq_out, limma_annot, deseq_annot, deseqSum_out, gene_annotation)

        
