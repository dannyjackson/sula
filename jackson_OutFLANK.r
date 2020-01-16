#!/usr/bin/env Rscript
#compliments jackson_OutFLANK.sh

args = commandArgs()
popNames = substr(args[6], 0, 100)

library("OutFLANK")
library("vcfR")

SNPmat <- read.table("outflank_SNPmat.txt")

locusNames <- read.table("outflank_loci_names.txt")
write.table(popNames, file="test.txt")

#Eliminate rows with NA
SNPmat <- SNPmat[complete.cases(SNPmat), ]


FstDataFrame <- MakeDiploidFSTMat(SNPmat,locusNames,popNames)

write.table(FstDataFrame, file = "outflank_FstDataFrame.txt")

pdf(file = "outflank_plot1.pdf", width = 20, height = 20, useDingbats=FALSE)
plot(FstDataFrame$FST, FstDataFrame$FSTNoCorr, xlim = c(-0.01,0.3),
     ylim = c(-0.01, 0.3), pch = 20)
abline(0, 1) # Checking the effect of sample size on Fst since FSTNoCorr will be used in the follow
  dev.off()


pdf(file = "outflank_plot2.pdf", width = 20, height = 20, useDingbats=FALSE)
hist(FstDataFrame$FSTNoCorr)
  dev.off()

object <- OutFLANK(FstDataFrame, NumberOfSamples=12, qthreshold = 0.1,
               RightTrimFraction = 0.4)

pdf(file = "outflank_plot3.pdf", width = 20, height = 20, useDingbats=FALSE)
OutFLANKResultsPlotter(object, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.01,
                      binwidth = 0.005, Zoom = FALSE, RightZoomFraction = 0.5,
                      titletext = NULL)
  dev.off()


# Identify Outlier
outliers <- object$results$LocusName[object$results$OutlierFlag == TRUE]

write.table(outliers, file="outflank_outliers.txt")
