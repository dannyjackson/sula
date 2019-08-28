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
