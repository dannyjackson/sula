#!/bin/bash

#shell script to run sliding window ABBA BABA tests, and then to export regions of the genome that show an f statistic higher than a determined value, and then to take those regions and pull gene names from gff files of the chicken genome.
#input is a vcf file. Also requires a file produced by satsuma termes "satsuma_summary.chained.out "
#requires chicken gff files to be in a folder with the names "chicken_1.gff" "chicken_2.gff" etc with the number corresponding to the chromosome
#requires Simon Marton's scripts as well

if [ $# -lt 1 ]
  then
    echo "Pulls list of lines from chicken gff file that Satsuma aligned to regions of high introgression, as identified with sliding ABBA BABA windows
    [-z] Name of project; name that all output files will be prefixed with
    [-a] Complete path to input vcf file
    [-b] Output directory for files. Must be empty.
    [-c] Path to satsuma_summary_chained.out file
    [-d] Path to directory containing gff files
    [-e] Path to directory containing all supporing files from github
    [-f] Path to Simon Martin's genomics_general-master directory
    [-g] File with two tab separated columns. Column 1 has sample names as they appear in the vcf, and column 2 has the population names to which each sample will be attributed. Populations should be listed together, not with some individuals from Pop1 followed by Pop2 followed by more Pop1. IMPORTANTLY the last population will be treated as your outgroup in the ABBA BABA analysis, or your P4, so be sure to list it last."

  else
    while getopts a:b:c:d:e:f: option
    do
    case "${option}"
    in
    z) name=${OPTARG};;
    a) vcf=${OPTARG};;
    b) outDir=${OPTARG};;
    c) chain=${OPTARG};;
    d) ref=${OPTARG};;
    e) path=${OPTARG};;
    f) simon=${OPTARG};;
    g) pops=${OPTARG};;

    esac
    done

mkdir $outDir/temp/
cat -n $pops | sort -uk3 | awk '{print $3}' > $outDir/temp/simplepops.txt

python $simon/VCF_processing/parseVCF.py -i $vcf -o $outDir/$name.geno.gz

#all of this is a work in progress. Trying to use the pop code file to set -p variables in simon martin's script.
#while read simplepops;
#do
#  echo '-p' $simplepops >> $outDir/temp/simplepops_variables.txt
#  tr '\n' ' ' $outDir/temp/simplepops_variables.txt >

#cat -n temp.txt | sort -uk3 | awk '{print $3}' | while read temp;
#do echo '-p' $temp
#done | tr '\n' ' '

python $simon/freq.py -g $outDir/$name.geno.gz \
-p BRBO -p MABO -p NABO -p BFBO -p PEBO -p RFBO --popsFile $pops -f phased --target derived \
-o A_1_nauto.geno.derFreq.tsv
done < $outDir/temp/simplepops.txt

Rscript $path/ABBAsliding.R $dataset $outDir $pops $name;

rm $outDir/temp/*
rm -r $outDir/temp/
