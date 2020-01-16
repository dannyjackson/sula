#!/bin/bash
#Runs Outflank given any vcf as input

if [ $# -lt 1 ]
  then
    echo "Runs OutFLANK given any vcf as input.
    [-v] vcf file
    [-o] Output directory for files
    [-p] Populations as comma separated string, e.g. bfbo,bfbo,pebo,pebo"
  else
    while getopts v:o: option
    do
    case "${option}"
    in
    v) dataset=${OPTARG};;
    o) outDir=${OPTARG};;
    p) pops=${OPTARG};;

    esac
    done

vcftools --vcf $dataset --mac 1 --non-ref-af 0.1 --recode --out outflank_invariant


# Convert vcf file to 012 file
vcftools --vcf outflank_invariant.recode.vcf --012 --out outflank_matrix.txt
# Replace the -1 by 9 in the 012 file
sed -i 's/-1/9/g' outflank_matrix.txt.012

# Remove the first column
awk '{$1=""; print substr($0,1)}' outflank_matrix.txt.012 > outflank_SNPmat.txt

# Take the loci name
cut -f 2 outflank_matrix.txt.012.pos > loci_names.txt

# Call R script

Rscript ~/sula/jackson_OutFLANK.R $pops

fi
