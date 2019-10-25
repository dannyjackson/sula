#!/usr/bin/env Rscript
#Plots the fd values output by the Simon Martin ABBA BABA sliding windows

#a file that is accessed by abba_complete.sh

args = commandArgs()

project_name = substr(args[grep("project_name_", args)],14,100000)
outputdirectory = substr(args[grep("outputdirectory_", args)],17,100000)

data <- read.csv(paste0(outputdirectory,"/",project_name,"_slidingwindows.csv.gz"))

#convert all fd values to 0 at sites where D is negative
for (x in 1:length(data)){
data[[x]]$fd = ifelse(data[[x]]$D < 0, 0, data[[x]]$fd)
    }

pdf(file = paste0(outputdirectory,"/",project_name,".slidingwindow_fdstats_plot.pdf"), width = 20, height = 7, useDingbats=FALSE)
hist(data$fd, xlim=c(0,1), breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1000), ylim=c(0,1000))
  dev.off()
