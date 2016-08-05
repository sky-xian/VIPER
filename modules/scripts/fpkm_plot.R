#!/usr/bin/env Rscript

# vim: syntax=r tabstop=4 expandtab

#---------------------------------
# @authors: Zach Herbert, Mahesh Vangala
# @emails: , vangalamaheshh@gmail.com
# @date: Aug, 05, 2016
#---------------------------------

library(reshape2)
library(ggplot2)

options(error = function() traceback(2))
args <- commandArgs( trailingOnly = TRUE )

cuff.csv <- args[1]
plot.pdf <- args[2]

cdat <- read.csv(cuff.csv, header=T, row.names=1)

fpkm_0.1 <- colSums(apply(cdat, 2,function(x) ifelse(x>0.1, 1, 0)))
fpkm_1 <- colSums(apply(cdat, 2,function(x) ifelse(x>1, 1, 0)))

fpkm <- cbind(fpkm_0.1,fpkm_1)
mfpkm <- melt(fpkm)

yl <- 1.2*(max(mfpkm$value))

pdf(plot.pdf)
ggplot(mfpkm,aes(Var1,value,fill=Var2))+
  ggtitle("\nGenes Detected\n")+
  geom_bar(stat="identity",position="dodge",width=.6)+
  geom_text(aes(label=value),vjust=-0.2, position=position_dodge(width=0.5),hjust=-0.25, size=2)+
  coord_cartesian(ylim=c(0,yl))+
  facet_grid(Var2~.)+
  theme(text = element_text(size=9), axis.text.x = element_text(angle=0)) +
  coord_flip()
dev.off()
