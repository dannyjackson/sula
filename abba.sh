#!/bin/bash

#A shell script to perform sliding window ABBA BABA tests on whole genome SNP data

if [ $# -lt 1 ]
  then
    echo "An automated bash interface for calling python and R scripts involved in Simon Martin's 4 taxon D-statistic"

  else
    while getopts a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p: option
    do
    case "${option}"
    in
    a) project_name=${OPTARG};;
    b) output_directory=${OPTARG};;
    c) github_directory=${OPTARG};;
    d) simonhmartin_directory=${OPTARG};;
    e) path_to_geno_file=${OPTARG};;
    f) path_to_populations_file=${OPTARG};;
    g) path_to_admx_populations_file=${OPTARG};;

    esac
    done

    pops=$(python $github_directory/parsepops.py $path_to_populations_file)

    python $simonhmartin_directory/freq.py -g /data5/sulidae/final/abba/sula.geno.gz \
    $pops --popsFile $path_to_populations_file -f phased --target derived \
    -o $output_directory/$project_name.geno.tsv

    pops_admx=$(python $github_directory/parsepops.py $path_to_admx_populations_file)

    python $simonhmartin_directory/freq.py -g /data5/sulidae/final/abba/sula.geno.gz \
    $pops_admx --popsFile $path_to_admx_populations_file -f phased --target derived \
    -o $output_directory/$project_name.admx.geno.tsv

    Rscript ${github_directory}/ABBAwholegenome.r inputfile_${output_directory}/${project_name}.geno.tsv inputadmx_${output_directory}/${project_name}.admx.geno.tsv outputdirectory_$output_directory simonhmartin_directory_$simonhmartin_directory project_name_$project_name
fi 
