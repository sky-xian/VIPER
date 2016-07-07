suppressMessages(library(limma))
suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))

limma_and_deseq_f <- function(arg_counts, meta, arg_s1, arg_s2, gene_annotation, method, limma_out, deseq_out, deseqSum_out, limma_annot, deseq_annot) {
    
    #READ in gene_annotation table--even though gene descriptions are quoted
    #in the annotations, we have to quote it again!
    if( grepl(".bz2$",gene_annotation) ) {
        gene_annot <- read.table(bzfile(gene_annotation), header=TRUE, sep=",", quote="", fill=TRUE, row.names=NULL)
    } else {
        gene_annot <- read.table(gene_annotation, header=TRUE, sep=",", quote="", fill=TRUE, row.names=NULL)
    }
    #Quote the gene_descriptions--NEED " instead of '
    gene_annot[,'Gene.Description'] <- sapply(gene_annot[,'Gene.Description'], function(x){paste0("\"",x,"\"")}, simplify="vector")

    ## Read in count matrix and control/treat lists
    treatlist = strsplit(arg_s2,',')[[1]]
    ctrllist = strsplit(arg_s1,',')[[1]]
    countmat <- read.table(arg_counts, header=TRUE, sep=",", row.names=1, check.names=FALSE)
    countmat = round(countmat)

    ## Pick out matrixes for the control and for the treat, and count the amounts
    ctrllist = as.data.frame(countmat[ ,colnames(countmat) %in% ctrllist])
    treatlist = as.data.frame(countmat[ ,colnames(countmat) %in% treatlist])
    ntreat = ncol(treatlist)
    nctrl = ncol(ctrllist)

    ## Create matrix with desired columns, in order or all treat, then all control
    data = cbind(treatlist,ctrllist)
    condition = c(rep('treat',ntreat),rep('control',nctrl))

    ## Read in metadata, will need if performing paired test
    metadata = read.table(meta, header=T, sep=",", row.names=1)
    if (method == "paired") {
        Pair = metadata[,"Pair",drop=FALSE]
        Pair = as.data.frame(Pair[rownames(Pair) %in% colnames(data),,drop=FALSE])
        Pair = unique(Pair)
        Pair = rbind(Pair,Pair)
    }

    ## Create colData, with sample names, and then condition, PLUS ANY NEW PARAMS
    colData <- as.data.frame(cbind(colnames(data),condition))
    if (method == "paired") {colData <- cbind(colData, Pair)}
    
    preparationD <- function (data, colData){

        DEdesign = ~ condition
        if (method == "paired") DEdesign = ~ Pair + condition
        
        dds <- DESeqDataSetFromMatrix(countData=data, colData=colData, design= DEdesign)

        dds <- dds[rowSums(counts(dds)) > 0, ]
        dds <- DESeq(dds)
        res <- results(dds)
        
        #LEN: get summary stats
        summary <- c(sum(res$padj<0.05 & res$log2FoldChange>0.0, na.rm=TRUE),
                     sum(res$padj<0.05 & res$log2FoldChange>0.5, na.rm=TRUE),
                     sum(res$padj<0.05 & res$log2FoldChange>1.0, na.rm=TRUE),
                     sum(res$padj<0.01 & res$log2FoldChange>0.0, na.rm=TRUE),
                     sum(res$padj<0.01 & res$log2FoldChange>0.5, na.rm=TRUE),
                     sum(res$padj<0.01 & res$log2FoldChange>1.0, na.rm=TRUE))
        sumTable <- matrix(summary, nrow=3, ncol=2)
        
        rownames(sumTable)<-c('log2FC > 0.0','log2FC > 0.5', 'log2FC > 1.0')
        colnames(sumTable)<-c('padj < 0.05','padj < 0.01')
        write.table(sumTable, deseqSum_out, quote=FALSE, sep=",")

        return (res)
    }
    deseq_result= preparationD(data,colData)
    deseq_result <- cbind(id=rownames(deseq_result), as.matrix(deseq_result))
    deseq_result <- deseq_result[order(as.numeric(deseq_result[,"padj"])),]

    write.table(deseq_result,deseq_out,sep=',',col.names=T,row.names=F,quote=F)

    ## Annotate w/ local .csv biomart annotation file
    deseq_annotations <- merge(deseq_result, gene_annot, by= "id", all.x = TRUE, sort=FALSE)
    write.table(deseq_annotations,deseq_annot,sep=',',col.names=T,row.names=F,quote=F)
    
    
    preparationL <- function(data,ntreat,nctrl){

        ## Define design
        Ldesign <- model.matrix(~condition)
        if (method == "paired") {
            Pair = factor(Pair[[1]])
            Ldesign = model.matrix(~ Pair + condition)
        }

        nf <- calcNormFactors(data)
        y <- voom (data, Ldesign, lib.size = colSums(data)*nf)
        fit <- lmFit(y,Ldesign)        
        fit <- eBayes(fit)
        
        res = topTable(fit,coef="conditiontreat",n=length(data[,1]),sort="p")
        
        return(res)
    }
    limma_result = preparationL(data,ntreat,nctrl)
    limma_result <- cbind(id=rownames(limma_result), limma_result)
    limma_result <- limma_result[order(as.numeric(limma_result[,5])),]
        
    write.table(limma_result,limma_out,sep=',',col.names=T,row.names=F,quote=F)

    ## ANNOTATE w/ local .csv biomart annotation file 
    limma_annotations <- merge(limma_result, gene_annot, by= "id", all.x = TRUE, sort=FALSE)
    write.table(limma_annotations,limma_annot,sep=',',col.names=T,row.names=F,quote=F)

}

limma_and_deseq_f(
    snakemake@input[["counts"]],
    snakemake@input[["meta"]],
    snakemake@params[["s1"]],
    snakemake@params[["s2"]],
    snakemake@params[["gene_annotation"]],
    snakemake@params[["method"]],
    snakemake@output[["limma_out"]],
    snakemake@output[["deseq_out"]],
    snakemake@output[["deseqSum_out"]],
    snakemake@output[["limma_annot"]],
    snakemake@output[["deseq_annot"]]
    )
