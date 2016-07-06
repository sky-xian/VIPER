suppressMessages(library(limma))
suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))

#limma_and_deseq_f <- function(arg_counts, meta, arg_s1, arg_s2, method, limma, deseq, limma_annot, deseq_annot, deseqSum_out, gene_annotation) {
limma_and_deseq_f <- function(arg_counts, meta, arg_s1, arg_s2, gene_annotation, method, limma_out, deseq_out, deseqSum_out, limma_annot, deseq_annot) {
    
    #READ in gene_annotation table--even though gene descriptions are quoted
    #in the annotations, we have to quote it again!
    if( grepl(".bz2$",gene_annotation) ) {
        gene_annot <- read.table(bzfile(gene_annotation), header=TRUE, sep=",", quote="", fill=TRUE, row.names=NULL)
    }
    else {
        gene_annot <- read.table(gene_annotation, header=TRUE, sep=",", quote="", fill=TRUE, row.names=NULL)
    }
    #DROP--move to the annotation files
    #Quote the gene_descriptions--VERSION 1
    #gene_annot[,'Gene.Description'] <- sapply(gene_annot[,'Gene.Description'],
    #                                          dQuote, simplify="vector")
    #Quote the gene_descriptions--NEED " instead of '
    gene_annot[,'Gene.Description'] <- sapply(gene_annot[,'Gene.Description'],
                                              function(x){paste0("\"",x,"\"")},
                                              simplify="vector")
    
    ## Read in lists for comparison, read in count matrix, and do "rounding" failsafe
    treatlist = strsplit(arg_s2,',')[[1]]
    ctrllist = strsplit(arg_s1,',')[[1]]
    countmat <- read.table(arg_counts, header=TRUE, sep=",", row.names=1, check.names=FALSE)
    countmat = round(countmat)
    
    ctrllist = as.data.frame(countmat[ ,colnames(countmat) %in% ctrllist])
    treatlist = as.data.frame(countmat[ ,colnames(countmat) %in% treatlist])
    
    ntreat = ncol(treatlist)
    nctrl = ncol(ctrllist)

    data = cbind(treatlist,ctrllist)

    condition = c(rep('treat',ntreat),rep('control',nctrl))

    colData <- as.data.frame(cbind(colnames(data),condition))
    
    preparationD <- function (countData, colData){
       
	dds <- DESeqDataSetFromMatrix(countData=data, colData=colData, design= ~ condition)

        dds <- dds[rowSums(counts(dds)) > 0, ]
        dds <- DESeq(dds)
        res <- results(dds)
        
        #ORDER results by XXX
        #res <- results[order(res$padj)]
        #summarize data
        #sum(res$padj < 0.1, na.rm=TRUE)
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
        #LEN: write/output summary stats
        write.table(sumTable, deseqSum_out, quote=FALSE, sep=",")

        #nbinomTest() default:
        #nbinomTest(cds, condA, condB, pvals_only = FALSE)
        return (res)
    }
    deseq_result= preparationD(data,colData)
    
    preparationL <- function(counts,ntreat,nctrl){
        #counts <- read.delim(table, header = FALSE, row.names = 1)
        mType <- c(rep('treat',ntreat),rep('control',nctrl))
        ###In case of needing to filtration uncomment two below commands:
        #isexpr <- rowSums(cpm(counts)>1)>= each
        #counts <- counts[isexpr,]
        nf <- calcNormFactors(counts)
        #calcNormFactors() default:
        #calcNormFactors(object, method=c("TMM","RLE","upperquartile","none"), refColumn = NULL,
        #logratioTrim = .3, sumTrim = 0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75)
        design <- model.matrix(~mType)
        y <- voom (counts, design, lib.size = colSums(counts)*nf)
        #voom() default:
        #voom(counts, design = NULL, lib.size = NULL, normalize.method = "none", plot = FALSE, span=0.5, ...)
        fit <- lmFit(y,design)        
        #lmFit() default:
        #lmFit(object,design=NULL,ndups=1,spacing=1,block=NULL,correlation,weights=NULL,method=c("ls","robust")) 
        fit <- eBayes(fit)
        
        res = topTable(fit,coef=2,n=length(counts[,1]),sort="p")
        
        #eBayes() default:
        #eBayes(fit, proportion=0.01, stdev.coef.lim=c(0.1,4), trend=FALSE, robust=FALSE, winsor.tail.p=c(0.05,0.1))
        return(res)
    }
    if (ntreat>1){
        limma_result = preparationL(data,ntreat,nctrl)
        limma_result <- cbind(id=rownames(limma_result), limma_result)
        #ANNOTATE w/ local .csv biomart annotation file
        limma_annotations <- merge(limma_result, gene_annot, by= "id", all.x = TRUE, sort=FALSE)
        print("test1")

        #MOVING to just csv files
        #write.table(limma_result,limma,sep='\t',col.names=T,row.names=F,quote=F)
        write.table(limma_result,limma_out,sep=',',col.names=T,row.names=F,quote=F)

        print("test2")
        #WRITE annotation table--limma.annot.csv
        write.table(limma_annotations,limma_annot,sep=',',col.names=T,row.names=F,quote=F)

        print("test3")
    }
    #LEN: setting the first column name to 'id'
    deseq_result <- cbind(id=rownames(deseq_result), as.matrix(deseq_result))

    print("test4")
    #LEN:  sort by padj. and remove padj = NA
    #deseq_result <-deseq_result[order(as.numeric(deseq_result[,"padj"]), na.last = NA),]
    deseq_result <-deseq_result[order(as.numeric(deseq_result[,"padj"])),]

    print("test5")
    #WRITE output deseq
    #MOVING to just csv files
    #write.table(deseq_result,deseq,sep='\t',col.names=T,row.names=F,quote=F)
    write.table(deseq_result,deseq_out,sep=',',col.names=T,row.names=F,quote=F)

    print("test6")

    #ANNOTATE w/ local .csv biomart annotation file
    deseq_annotations <- merge(deseq_result, gene_annot, by= "id", all.x = TRUE, sort=FALSE)

    print("test7")
    #LEN:  sort by padj. and remove padj = NA
    #gene_annot <- gene_annot[order(as.numeric(gene_annot[,"padj"]), na.last = NA),]
    #MOVING to just csv files
    #write.table(gene_annot,paste(deseq, "annot", sep="."),sep='\t',col.names=T,row.names=F,quote=F)
    write.table(deseq_annotations,deseq_annot,sep=',',col.names=T,row.names=F,quote=F)
    print("test8")

}

args <- commandArgs( trailingOnly = TRUE )
#arg_counts = strsplit(args[1]," ")[[1]]
#arg_s1 = strsplit(args[2], ",")
#arg_s2 = strsplit(args[3], ",")
arg_counts = args[1]
meta = args[2]
arg_s1 = args[3]
arg_s2 = args[4]
method = args[5]
limma_out=args[6]
deseq_out=args[7]
limma_annot = args[8]
deseq_annot = args[9]
deseqSum_out=args[10]
gene_annotation=args[11]

#limma_and_deseq_f(arg_counts, meta, arg_s1, arg_s2, method, limma_out, deseq_out, limma_annot, deseq_annot, deseqSum_out, gene_annotation)

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
