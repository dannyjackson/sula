#!/usr/bin/env Rscript
#Genome wide ABBA BABA test with block jackknife procedure

args = commandArgs()

inputfile = "brbo_all.geno.tsv"
outputdirectory = "/data5/sulidae/final/abba/brbo/"
simonhmartin_directory = "~/genomics_general-master/"
project_name = "brbo"
inputpops = "brbo_pops_1.txt"

D.stat <- function(p1, p2, p3) {
    ABBA <- (1 - p1) * p2 * p3
    BABA <- p1 * (1 - p2) * p3
    (sum(ABBA) - sum(BABA)) / (sum(ABBA) + sum(BABA))
    }

freq_table = read.table(paste0(inputfile), header=T, as.is=T)

populations = read.table("brbo_pops_1.txt", header=T, as.is=T)

#removing NA values
freq_table2<-freq_table
freq_table<-na.omit(freq_table2)

source(paste0(simonhmartin_directory,"/jackknife.R"))

block_indices <- get_block_indices(block_size=1e6,
                                   positions=freq_table$position,
                                   chromosomes=freq_table$scaffold)

n_blocks <- length(block_indices)

cat(paste("Genome divided into", n_blocks, "blocks."),file=paste0(outputdirectory,"/",project_name,".abbawholegenome.stats.txt"),sep="\n",append=TRUE)

for (row in 1:nrow(populations)){
  P1 <- populations[row, "P1"]
  P2 <- populations[row, "P2"]
  P3 <- populations[row, "P3"]

  D <- D.stat(freq_table[,P1], freq_table[,P2], freq_table[,P3])

  D_sd <- get_jackknife_sd(block_indices=block_indices,
                           FUN=D.stat,
                           freq_table[,P1], freq_table[,P2], freq_table[,P3])
  D_err <- D_sd/sqrt(n_blocks)
  D_Z <- D / D_err
  D_p <- 2*pnorm(-abs(D_Z))

  cat(paste(P1,P2,P3,round(D,4), round(D_sd,4), round(D_Z,3), round(D_p,3), sep="\t"), file=paste0(outputdirectory,"/",project_name,".abbawholegenome.stats.txt"),sep="\n",append=TRUE)
}
