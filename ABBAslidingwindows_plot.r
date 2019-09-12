#!/usr/bin/env Rscript
#Plots the output of the Simon Martin script ABBABABAwindows.py

#a file that is accessed by abba_complete.sh

args = commandArgs()

project_name = substr(args[grep("project_name_", args)],11,100000)
outputdirectory = substr(args[grep("outputdirectory_", args)],17,100000)

AB_files <- c(paste0($output_directory,"/",project_name,"_slidingwindows.csv.gz"))
AB_tables = lapply(AB_files, read.csv)
head(AB_tables[[1]])

#convert all fd values to 0 at sites where D is negative
for (x in 1:length(AB_tables)){
AB_tables[[x]]$fd = ifelse(AB_tables[[x]]$D < 0, 0, AB_tables[[x]]$fd)
    }

par(mfrow=c(length(AB_tables), 1), mar = c(4,4,1,1))

pdf(file = paste0(outputdirectory,"/",project_name,"slidingwindow_fdstats_plot.pdf"), width = 20, height = 7, useDingbats=FALSE)
for (x in 1:length(AB_tables)){
    plot(AB_tables[[x]]$fd,
    type = "l", xlim=c(0,17e6),ylim=c(0,1),ylab="Admixture Proportion",xlab="Position")
    rect(1000000,0,1250000,1, col = rgb(0,0,0,0.2), border=NA)
    }
  dev.off()
