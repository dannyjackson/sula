#!/bin/bash
#Runs Outflank given any vcf as input

if [ $# -lt 1 ]
  then
    echo "Runs OutFLANK given any vcf as input. Requires file in working directory titled pops.txt
    [-v] vcf file, [-p] number of individuals in population 1, [-q] number of individuals in population 2"
  else
    while getopts v:p:q: option
    do
    case "${option}"
    in
    v) dataset=${OPTARG};;
    p) p1_indv=${OPTARG};;
    q) p2_indv=${OPTARG};;
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
cut -f 2 outflank_matrix.txt.012.pos > outflank_loci_names.txt

# Remove sites that are invariant within each population
python ~/sula/filter_missingPop.py $p1_indv $p2_indv

# Call R script

Rscript ~/sula/jackson_OutFLANK.r

fi
