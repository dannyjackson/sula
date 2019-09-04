#!/bin/bash

#shell script to pull list of gene names from chicken gff file
#requires chicken gff files to be in a folder with the names "chicken_chr1.gff" "chicken_chr2.gff" etc
#input is a csv file of a list of regions output by Simon Martin's ABBA BABA sliding window scans

if [ $# -lt 1 ]
  then
    echo "Pulls list of lines from chicken gff file that Satsuma aligned to regions of high introgression, as identified with sliding ABBA BABA windows
    [-i] Path to input csv file
    [-o] Output directory for files
    [-c] Path to satsuma_summary_chained.out file
    [-r] Path to directory containing gff files"

  else
    while getopts i:o:c:r: option
    do
    case "${option}"
    in
    i) dataset=${OPTARG};;
    o) outDir=${OPTARG};;
    c) chain=${OPTARG};;
    r) ref=${OPTARG};;

    esac
    done

    esac
    done

sed 's|,|\t|g' $dataset.csv > $dataset.tsv

awk '{if ($11 > 0.3) print $0;}' $dataset.tsv > $dataset_0pt3.tsv

while read -r file;
do
scaffold=$(awk '{print $1}' <<< "$file")
scaffoldstart=$(awk '{print $2}' <<< "$file")
scaffoldend=$(awk '{print $3}' <<< "$file")

awk '{if ($4 == "$scaffold") print $0;}' /data5/satsuma_output/maskedbooby_chicken_output/satsuma_summary.chained.out | awk '{if ($5 > $scaffoldstart) print $0;}' | awk '{if ($5 < $scaffoldend) print $0;}' >> $dataset_satsuma.tsv
done

#determine the range of the chromosome to which this region mapped on the chicken
awk '{print $2}' $dataset_satsuma.tsv > working.txt
awk '{print $3}' $dataset_satsuma.tsv >> working.txt
start=$(sort working.txt | head -n 1)
end=$(sort working.txt | tail -n 1)

#gives a range of 178084359 178105995
#So I want to awk the chicken chromosome 1 gff file for any genes in this region
#gff $4 is start; $5 is end
#I imaging that there might be a rare occasion where a gene straddles my boundary. In that case, one of these likely would grab it with the others would leave it off. Or not. Who knows. Regardless, you only need one of them.
awk '{if ($4 > $start && $4 < $end || $5 > $start && $5 < $end) print $0}' ~/reference_datasets/Gallus_gallus.GRCg6a.97.chromosome.1.gff3 > $dataset_gff.tsv

sed -i 's/*\_chromosome_//g'
sed -i 's/:.*//g' phala_scaffoldnumbers.txt
