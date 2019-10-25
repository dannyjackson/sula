#!/usr/bin/env Rscript
#Plots the fd values output by the Simon Martin ABBA BABA sliding windows

#a file that is accessed by abba_complete.sh

args = commandArgs()

project_name = substr(args[grep("project_name_", args)],14,100000)
outputdirectory = substr(args[grep("outputdirectory_", args)],17,100000)

data <- read.csv(paste0(outputdirectory,"/",project_name,"_slidingwindows.csv.gz"))

data_fd <- c(data$fd)

lesszero <- data_fd[data_fd<0]

cat(paste("Number of fd sites < 0 = ", length(lesszero)),file=paste0(outputdirectory,"/",project_name,".abbawholegenome.stats.txt"),sep="\n",append=TRUE)

greaterone <- data_fd[data_fd>1]

cat(paste("Number of fd sites > 1 = ", length(greaterone)),file=paste0(outputdirectory,"/",project_name,".abbawholegenome.stats.txt"),sep="\n",append=TRUE)
cat(paste0("List of fd values > 1 :"),file=paste0(outputdirectory,"/",project_name,".abbawholegenome.stats.txt"),sep="\n",append=TRUE)

cat(paste("Length of total values before NA filtering : ", length(data_fd[!is.na(data_fd)])),file=paste0(outputdirectory,"/",project_name,".abbawholegenome.stats.txt"),sep="\n",append=TRUE)

data_fd[data_fd<0] <- NA
data_fd[data_fd>1] <- NA

cat(paste("Length of total values after NA filtering : ", length(data_fd[!is.na(data_fd)])),file=paste0(outputdirectory,"/",project_name,".abbawholegenome.stats.txt"),sep="\n",append=TRUE)

pdf(file = paste0(outputdirectory,"/",project_name,".slidingwindow_fdstats_plot.pdf"), width = 20, height = 7, useDingbats=FALSE)
hist(data_fd, xlim=c(0,1))
  dev.off()
