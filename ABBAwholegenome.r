#!/usr/bin/env Rscript
#Genome wide ABBA BABA test with block jackknife procedue

#a file that is accessed by abba_complete.sh

args = commandArgs()

inputfile = substr(args[grep("inputfile_", args)],11,100000)
outputdirectory = substr(args[grep("outputdirectory_", args)],17,100000)
simonhmartin_directory = substr(args[grep("simonhmartin_directory", args)],17,100000)
P1 = paste0(substr(args[grep("population1_", args)],13,100000))
P2 = paste0(substr(args[grep("population2_", args)],13,100000))
P3 = paste0(substr(args[grep("population3_", args)],13,100000))


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

print(paste("D =", round(D,4)))

source(paste0(simonhmartin_directory,"jackknife.R"))

block_indices <- get_block_indices(block_size=1e6,
                                   positions=freq_table$position,
                                   chromosomes=freq_table$scaffold)

n_blocks <- length(block_indices)

print(paste("Genome divided into", n_blocks, "blocks."))

D_sd <- get_jackknife_sd(block_indices=block_indices,
                         FUN=D.stat,
                         freq_table[,P1], freq_table[,P2], freq_table[,P3])

print(paste("D standard deviation = ", round(D_sd,4)))

D_err <- D_sd/sqrt(n_blocks)
D_Z <- D / D_err

print(paste("D Z score = ", round(D_Z,3)))
