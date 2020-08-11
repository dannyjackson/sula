#!/usr/bin/env Rscript
#Fst.R

#a file that is accessed by fst_slidingwindow.sh to create a pdf of an fst sliding window scan (just like it sounds dummy)

args = commandArgs()

outDir <- args[7]
name <- args[8]

library(qqman)
fst<-read.table(paste0("/data5/sulidae/final/fst/bfbo_pebo/fst.txt.windowed.weir.fst"), header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = paste0(outDir,"/",name,".pdf"), width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="BIN_START",p="WEIGHTED_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst",cex = 0.2))
dev.off()
