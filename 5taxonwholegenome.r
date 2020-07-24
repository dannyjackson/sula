#!/usr/bin/env Rscript

Dfo.stat <- function(p1, p2, p3, p4) {

    BABAA <- p1 * (1 - p2) * p3 * (1 - p4)
    BBBAA <- p1 * p2 * p3 * (1 - p4)
    ABABA <- (1 - p1) * p2 * (1 - p3) * p4
    AAABA <- (1 - p1) * (1 - p2) * (1 - p3) * p4

    BAABA <- p1 * (1 - p2) * (1 - p3) * p4
    BBABA <- p1 * p2 * (1 - p3) * p4
    ABBAA <- (1 - p1) * p2 * p3 * (1 - p4)
    AABAA <- (1 - p1) * (1 - p2) * p3 * (1 - p4)

    first_fo <- (sum(BABAA) + sum(BBBAA) + sum(ABABA) + sum(AAABA))
    second_fo <- (sum(BAABA) + sum(BBABA) + sum(ABBAA) + sum(AABAA))

    (first_fo - second_fo) / (first_fo + second_fo)
    }

Dil.stat <- function(p1, p2, p3, p4) {

    ABBAA <- (1 - p1) * p2 * p3 * (1 - p4)
    BBBAA <- p1 * p2 * p3 * (1 - p4)
    BAABA <- p1 * (1 - p2) * (1 - p3) * p4
    AAABA <- (1 - p1) * (1 - p2) * (1 - p3) * p4

    ABABA <- (1 - p1) * p2 * (1 - p3) * p4
    BBABA <- p1 * p2 * (1 - p3) * p4
    BABAA <- p1 * (1 - p2) * p3 * (1 - p4)
    AABAA <- (1 - p1) * (1 - p2) * p3 * (1 - p4)

    first_il <- (sum(ABBAA) + sum(BBBAA) + sum(BAABA) + sum(AAABA))
    second_il <- (sum(ABABA) + sum(BBABA) + sum(BABAA) + sum(AABAA))

    (first_il - second_il) / (first_il + second_il)
    }


Dfi.stat <- function(p1, p2, p3, p4) {

    BABAA <- p1 * (1 - p2) * p3 * (1 - p4)
    BABBA <- p1 * (1 - p2) * p3 * p4
    ABABA <- (1 - p1) * p2 * (1 - p3) * p4
    ABAAA <- (1 - p1) * p2 * (1 - p3) * (1 - p4)

    ABBAA <- (1 - p1) * p2 * p3 * (1 - p4)
    ABBBA <- (1 - p1) * p2 * p3 * p4
    BAABA <- p1 * (1 - p2) * (1 - p3) * p4
    BAAAA <- p1 * (1 - p2) * (1 - p3) * (1 - p4)

    first_fi <- (sum(BABAA) + sum(BABBA) + sum(ABABA) + sum(ABAAA))
    second_fi <- (sum(ABBAA) + sum(ABBBA) + sum(BAABA) + sum(BAAAA))

    (first_fi - second_fi) / (first_fi + second_fi)
    }

Dol.stat <- function(p1, p2, p3, p4) {

    BAABA <- p1 * (1 - p2) * (1 - p3) * p4
    BABBA <- p1 * (1 - p2) * p3 * p4
    ABBAA <- (1 - p1) * p2 * p3 * (1 - p4)
    ABAAA <- (1 - p1) * p2 * (1 - p3) * (1 - p4)

    ABABA <- (1 - p1) * p2 * (1 - p3) * p4
    ABBBA <- (1 - p1) * p2 * p3 * p4
    BABAA <- p1 * (1 - p2) * p3 * (1 - p4)
    BAAAA <- p1 * (1 - p2) * (1 - p3) * (1 - p4)

    first_ol <- (sum(BAABA) + sum(BABBA) + sum(ABBAA) + sum(ABAAA))
    second_ol <- (sum(ABABA) + sum(ABBBA) + sum(BABAA) + sum(BAAAA))

    (first_ol - second_ol) / (first_ol + second_ol)
    }


freq_table = read.table("mabo.geno.tsv", header=T, as.is=T)

#removing NA values
freq_table2<-freq_table
freq_table<-na.omit(freq_table2)

Dfo <- Dfo.stat(freq_table[,"P1"], freq_table[,"P2"], freq_table[,"P3"], freq_table[,"P4"])
Dil <- Dil.stat(freq_table[,"P1"], freq_table[,"P2"], freq_table[,"P3"], freq_table[,"P4"])
Dfi <- Dfi.stat(freq_table[,"P1"], freq_table[,"P2"], freq_table[,"P3"], freq_table[,"P4"])
Dol <- Dol.stat(freq_table[,"P1"], freq_table[,"P2"], freq_table[,"P3"], freq_table[,"P4"])


Dfo
Dil
Dfi
Dol

cat(paste("Dfo =", round(Dfo,4)),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)
cat(paste("Dil =", round(Dil,4)),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)
cat(paste("Dfi =", round(Dfi,4)),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)
cat(paste("Dol =", round(Dol,4)),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)


source("~/genomics_general-master/jackknife.R")

block_indices <- get_block_indices(block_size=1e6,
                                   positions=freq_table$position,
                                   chromosomes=freq_table$scaffold)

n_blocks <- length(block_indices)

cat(paste("Genome divided into", n_blocks, "blocks."),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)


Dfo_sd <- get_jackknife_sd(block_indices=block_indices,
                         FUN=Dfo.stat,
                         freq_table[,"P1"], freq_table[,"P2"], freq_table[,"P3"], freq_table[,"P4"])


Dfo_err <- Dfo_sd/sqrt(n_blocks)
Dfo_Z <- Dfo / Dfo_err
Dfo_p <- 2*pnorm(-abs(Dfo_Z))

cat(paste("Dfo_sd =", round(Dfo_sd,4)),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)
cat(paste("Dil_err =", round(Dil_err,4)),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)
cat(paste("Dfo_Z =", round(Dfo_Z,4)),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)
cat(paste("Dfo_p =", round(Dfo_p,4)),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)




Dil_sd <- get_jackknife_sd(block_indices=block_indices,
                         FUN=Dil.stat,
                         freq_table[,"P1"], freq_table[,"P2"], freq_table[,"P3"], freq_table[,"P4"])


Dil_err <- Dil_sd/sqrt(n_blocks)
Dil_Z <- Dil / Dil_err
Dil_p <- 2*pnorm(-abs(Dil_Z))

cat(paste("Dil_sd =", round(Dil_sd,4)),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)
cat(paste("Dil_err =", round(Dil_err,4)),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)
cat(paste("Dil_Z =", round(Dil_Z,4)),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)
cat(paste("Dil_p =", round(Dil_p,4)),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)



Dfi_sd <- get_jackknife_sd(block_indices=block_indices,
                         FUN=Dfi.stat,
                         freq_table[,"P1"], freq_table[,"P2"], freq_table[,"P3"], freq_table[,"P4"])



Dfi_err <- Dfi_sd/sqrt(n_blocks)
Dfi_Z <- Dfi / Dfi_err
Dfi_p <- 2*pnorm(-abs(Dfi_Z))

cat(paste("Dfi_sd =", round(Dfi_sd,4)),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)
cat(paste("Dfi_err =", round(Dfi_err,4)),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)
cat(paste("Dfi_Z =", round(Dfi_Z,4)),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)
cat(paste("Dfi_p =", round(Dfi_p,4)),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)



Dol_sd <- get_jackknife_sd(block_indices=block_indices,
                         FUN=Dol.stat,
                         freq_table[,"P1"], freq_table[,"P2"], freq_table[,"P3"], freq_table[,"P4"])


Dol_err <- Dol_sd/sqrt(n_blocks)
Dol_Z <- Dol / Dol_err
Dol_p <- 2*pnorm(-abs(Dol_Z))

cat(paste("Dol_sd =", round(Dol_sd,4)),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)
cat(paste("Dol_err =", round(Dol_err,4)),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)
cat(paste("Dol_Z =", round(Dol_Z,4)),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)
cat(paste("Dol_p =", round(Dol_p,4)),file=paste0("5taxon.abbawholegenome.stats.txt"),sep="\n",append=TRUE)

























# Estimate Admixture Proportion

admx_table = read.table(paste0(inputadmx), header=T, as.is=T)

#removing NA values
admx_table2<-admx_table
admx_table<-na.omit(admx_table2)


block_indices <- get_block_indices(block_size=1e6,
                                   positions=admx_table$position,
                                   chromosomes=admx_table$scaffold)

n_blocks <- length(block_indices)

cat(paste("Genome divided into", n_blocks, "blocks."),file=paste0(outputdirectory,"/",project_name,".abbawholegenome.stats.txt"),sep="\n",append=TRUE)


D_err <- D_sd/sqrt(n_blocks)
D_Z <- D / D_err
D_p <- 2*pnorm(-abs(D_Z))

cat(paste("D Z score = ", round(D_Z,3)),file=paste0(outputdirectory,"/",project_name,".abbawholegenome.stats.txt"),sep="\n",append=TRUE)

cat(paste("D p value = ", round(D_p,3)),file=paste0(outputdirectory,"/",project_name,".abbawholegenome.stats.txt"),sep="\n",append=TRUE)








f.stat <- function(p1, p2, p3a, p3b) {
    ABBA_numerator <- (1 - p1) * p2 * p3a
    BABA_numerator <- p1 * (1 - p2) * p3a

    ABBA_denominator <- (1 - p1) * p3b * p3a
    BABA_denominator <- p1 * (1 - p3b) * p3a

    (sum(ABBA_numerator) - sum(BABA_numerator)) /
    (sum(ABBA_denominator) - sum(BABA_denominator))
    }

f <- f.stat(admx_table[,"P1"], admx_table[,"P2"], admx_table[,"P3a"], admx_table[,"P3b"])


cat(paste("Admixture proportion (f value) = ", round(f,3)),file=paste0(outputdirectory,"/",project_name,".abbawholegenome.stats.txt"),sep="\n",append=TRUE)


f_sd <- get_jackknife_sd(block_indices=block_indices,
                         FUN=f.stat,
                         admx_table[,"P1"], admx_table[,"P2"], admx_table[,"P3a"], admx_table[,"P3b"])

cat(paste("Admixture proportion (f value) = ", round(f,3)),file=paste0(outputdirectory,"/",project_name,".abbawholegenome.stats.txt"),sep="\n",append=TRUE)


f_err <- f_sd/sqrt(n_blocks)

f_CI_lower <- f - 1.96*f_err
f_CI_upper <- f + 1.96*f_err

cat(paste("95% confidence interval of f =", round(f_CI_lower,4), round(f_CI_upper,4)), file=paste0(outputdirectory,"/",project_name,".abbawholegenome.stats.txt"),sep="\n",append=TRUE)

cat(paste("f err =", round(f_err,4)), file=paste0(outputdirectory,"/",project_name,".abbawholegenome.stats.txt"),sep="\n",append=TRUE)
