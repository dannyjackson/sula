#!/usr/bin/env Rscript
#Genome wide ABBA BABA test with block jackknife procedure

#a file that is accessed by abba_complete.sh

args = commandArgs()

inputfile = substr(args[grep("inputfile_", args)],11,100000)
outputdirectory = substr(args[grep("outputdirectory_", args)],17,100000)
simonhmartin_directory = substr(args[grep("simonhmartin_directory_", args)],24,100000)
P1 = paste0(substr(args[grep("population1_", args)],13,100000))
P2 = paste0(substr(args[grep("population2_", args)],13,100000))
P3 = paste0(substr(args[grep("population3_", args)],13,100000))

print("this1")
print(simonhmartin_directory)
D.stat <- function(p1, p2, p3) {
    ABBA <- (1 - p1) * p2 * p3
    BABA <- p1 * (1 - p2) * p3
    (sum(ABBA) - sum(BABA)) / (sum(ABBA) + sum(BABA))
    }

freq_table = read.table(paste0(inputfile), header=T, as.is=T)

#check for buried NAs
for (i in 1:length(freq_table[1,])){
  print(anyNA(freq_table[,i]))
}

#removing these values
freq_table2<-freq_table
freq_table<-na.omit(freq_table2)

D <- D.stat(freq_table[,P1], freq_table[,P2], freq_table[,P3])

#print_this_1 <- print(paste("D =", round(D,4)))
#write(print_this_1, file = paste0(outputdirectory,".abbawholegenome.stats.txt") append = FALSE)
cat(paste("D =", round(D,4)),file=paste0(outputdirectory,".abbawholegenome.stats.txt",sep="\n",append=TRUE))
print("this?")
source(paste0(simonhmartin_directory,"/jackknife.R"))

block_indices <- get_block_indices(block_size=1e6,
                                   positions=freq_table$position,
                                   chromosomes=freq_table$scaffold)

n_blocks <- length(block_indices)

#print_this_2 <- print(paste("Genome divided into", n_blocks, "blocks."))
#write(print_this_2, file = paste0(outputdirectory,".abbawholegenome.stats.txt") append = TRUE)
cat(paste("Genome divided into", n_blocks, "blocks."),file=paste0(outputdirectory,".abbawholegenome.stats.txt",sep="\n",append=TRUE)

D_sd <- get_jackknife_sd(block_indices=block_indices,
                         FUN=D.stat,
                         freq_table[,P1], freq_table[,P2], freq_table[,P3])

#print_this_3 <- print(paste("D standard deviation = ", round(D_sd,4)))
#write(print_this_3, file = paste0(outputdirectory,".abbawholegenome.stats.txt") append = TRUE)
cat(paste("D standard deviation = ", round(D_sd,4)),file=paste0(outputdirectory,".abbawholegenome.stats.txt",sep="\n",append=TRUE)

D_err <- D_sd/sqrt(n_blocks)
D_Z <- D / D_err

#print_this_4 <- print(paste("D Z score = ", round(D_Z,3)))
#write(print_this_4, file = paste0(outputdirectory,".abbawholegenome.stats.txt", append = TRUE)

cat(paste("D Z score = ", round(D_Z,3)),file=paste0(outputdirectory,".abbawholegenome.stats.txt",sep="\n",append=TRUE)
