#!/usr/bin/env Rscript
#Plots the fd values output by the Simon Martin ABBA BABA sliding windows

#a file that is accessed by abba_complete.sh

args = commandArgs()

project_name = substr(args[grep("project_name_", args)],14,100000)
outputdirectory = substr(args[grep("outputdirectory_", args)],17,100000)

data <- read.csv(paste0(outputdirectory,"/",project_name,"_slidingwindows.csv.gz"))

pdf(file = paste0(outputdirectory,"/",project_name,".slidingwindow_fdstats_plot.pdf"), width = 20, height = 7, useDingbats=FALSE)
hist(data$fd, xlim=c(0,1), breaks =   c(-1000,0,0.05,0.1,0.15,0.20,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1000), ylim=c(0,3))
  dev.off()
