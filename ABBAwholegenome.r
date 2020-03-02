#!/usr/bin/env Rscript
#Genome wide ABBA BABA test with block jackknife procedure

#a file that is accessed by abba_complete.sh

args = commandArgs()

inputfile = substr(args[grep("inputfile_", args)],11,100000)
inputadmx = substr(args[grep("inputadmx_", args)],11,100000)
outputdirectory = substr(args[grep("outputdirectory_", args)],17,100000)
simonhmartin_directory = substr(args[grep("simonhmartin_directory_", args)],24,100000)
project_name = substr(args[grep("project_name_", args)],14,100000)

D.stat <- function(p1, p2, p3) {
    ABBA <- (1 - p1) * p2 * p3
    BABA <- p1 * (1 - p2) * p3
    (sum(ABBA) - sum(BABA)) / (sum(ABBA) + sum(BABA))
    }

freq_table = read.table(paste0(inputfile), header=T, as.is=T)

#removing NA values
freq_table2<-freq_table
freq_table<-na.omit(freq_table2)

D <- D.stat(freq_table[,"P1"], freq_table[,"P2"], freq_table[,"P3"])

cat(paste("D =", round(D,4)),file=paste0(outputdirectory,"/",project_name,".abbawholegenome.stats.txt"),sep="\n",append=TRUE)

source(paste0(simonhmartin_directory,"/jackknife.R"))

block_indices <- get_block_indices(block_size=1e6,
                                   positions=freq_table$position,
                                   chromosomes=freq_table$scaffold)

n_blocks <- length(block_indices)

cat(paste("Genome divided into", n_blocks, "blocks."),file=paste0(outputdirectory,"/",project_name,".abbawholegenome.stats.txt"),sep="\n",append=TRUE)

D_sd <- get_jackknife_sd(block_indices=block_indices,
                         FUN=D.stat,
                         freq_table[,P1], freq_table[,P2], freq_table[,P3])

cat(paste("D standard deviation = ", round(D_sd,4)),file=paste0(outputdirectory,"/",project_name,".abbawholegenome.stats.txt"),sep="\n",append=TRUE)

D_err <- D_sd/sqrt(n_blocks)
D_Z <- D / D_err
D_p <- 2*pnorm(-abs(D_Z))

cat(paste("D Z score = ", round(D_Z,3)),file=paste0(outputdirectory,"/",project_name,".abbawholegenome.stats.txt"),sep="\n",append=TRUE)

cat(paste("D p value = ", round(D_p,3)),file=paste0(outputdirectory,"/",project_name,".abbawholegenome.stats.txt"),sep="\n",append=TRUE)



# Estimate Admixture Proportion

admx_table = read.table(paste0(inputadmx), header=T, as.is=T)

#removing NA values
admx_table2<-admx_table
admx_table<-na.omit(admx_table2)

ABBA_1_2_3a = abba(admx_table[,"P1"], admx_table[,"P2"], admx_table[,"P3a"])
BABA_1_2_3a = baba(admx_table[,"P1"], admx_table[,"P2"], admx_table[,"P3a"])

ABBA_1_3b_3a = abba(admx_table[,"P1"], admx_table[,"P3b"], admx_table[,"P3a"])
BABA_1_3b_3a = baba(admx_table[,"P1"], admx_table[,"P3b"], admx_table[,"P3a"])

f = (sum(ABBA_1_2_3a) - sum(BABA_1_2_3a))/
     (sum(ABBA_1_3b_3a) - sum(BABA_1_3b_3a))

cat(paste("Admixture proportion (f value) = ", round(f,3)),file=paste0(outputdirectory,"/",project_name,".abbawholegenome.stats.txt"),sep="\n",append=TRUE)
