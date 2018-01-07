#Script to run enricher fn on most significantly diffexp genes

suppressMessages(library("clusterProfiler"))

gsea <- function(deseqTable, gsea_db, comp_title, out_path) {
   #get set of significant diffexp genes--threshold 0.05
   genes <- deseqTable[deseqTable$padj < 0.05,]

   #sort by padj --DROP b/c this is already sorted by padj
   #genes <- genes[order(genes$padj),]

   #get list of gene names
   gene_list <- rownames(genes)
   
   #WRITE this as out_gene_list
   #write.list(gene_list, "gene_list.txt", quote=F, col.names = NA, sep="\t")
   write(gene_list, file=paste0(out_path, ".gene_list.txt"), ncolumns=1)

   #CHECK if the egmt output already exists
   out_egmt = paste0(out_path, ".egmt.Rda")
   if (file.exists(out_egmt)) {
      #print("Exists")
      #LOAD
      load(out_egmt)
   } else {
      #print("New")
      #START from scratch
      #running enrichment on gene symbol
      #READ in gmt
      gmt <- read.gmt(gsea_db)
      egmt <- enricher(gene_list, TERM2GENE=gmt)
      save(egmt, file=out_egmt)
   }
   #WRITE out summary
   write.csv(summary(egmt), file=paste0(out_path, ".gene_set.enrichment.txt"))

   #generate barplot
   png(paste0(out_path, ".gene_set.enrichment.barplot.png"), width = 8, height = 8, unit="in",res=300)
   #NOTE: if I just did barplot(...), the png doesn't get saved.  
   #I need the print(p)
   p1 <- barplot(egmt, font.size=6, title=comp_title)
   print(p1)
   junk <- dev.off()

   #generate dotplot
   png(paste0(out_path, ".gene_set.enrichment.dotplot.png"), width = 8, height = 8, unit="in",res=300)
   p2<-dotplot(egmt, font.size=6, title=comp_title)
   print(p2)
   junk <- dev.off()
   
   #TRY GSEA plot
   #top_hit<-summary(egmt)[1,]
   #gseaplot(egmt, TERM2GENE=gmt)#geneSetID=top_hit$ID)
   #y<-GSEA(gene_list, TERM2GENE=gmt)
   #head(summary(y))
}

## Read in arguments
args <- commandArgs( trailingOnly = TRUE )
deseqFile=args[1]
gsea_db=args[2]
title=args[3]
out_path=args[4]

deseqTable <- read.csv(deseqFile, header=T, check.names=F, row.names=1, stringsAsFactors=FALSE, dec='.')
gsea(deseqTable, gsea_db, title, out_path) 
