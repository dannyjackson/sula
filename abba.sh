#!/bin/bash

#A shell script to perform sliding window ABBA BABA tests on whole genome SNP data

if [ $# -lt 1 ]
  then
    echo "I still have to write an introduction to this script."

  else
    while getopts a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p: option
    do
    case "${option}"
    in
    a) project_name=${OPTARG};;
    b) output_directory=${OPTARG};;
    c) github_directory=${OPTARG};;
    d) simonhmartin_directory=${OPTARG};;
    e) path_to_vcf_file=${OPTARG};;
    f) path_to_populations_file=${OPTARG};;
    g) path_to_admx_populations_file=${OPTARG};;
    j) populationoutgroup=${OPTARG};;
    k) proceed_through_stage_2=${OPTARG};;
    l) windowsize=${OPTARG};;
    m) minimumsnps=${OPTARG};;
    n) threads=${OPTARG};;
    o) fstat_threshold=${OPTARG};;


    esac
    done

    #this converts the vcf file into a format that Simon Martin refers to as a ".geno", see his website for more information https://github.com/simonhmartin/genomics_general/tree/master/VCF_processing

    if [ ${path_to_vcf_file: -3} = "vcf" ]
      then
        echo "file read as vcf"
        python $simonhmartin_directory/VCF_processing/parseVCF.py -i $path_to_vcf_file -o $output_directory/$project_name.geno.gz

        pops=$(python $github_directory/parsepops.py $path_to_populations_file)

        python $simonhmartin_directory/freq.py -g $output_directory/$project_name.geno.gz \
        $pops --popsFile $path_to_populations_file -f phased --target derived \
        -o $output_directory/$project_name.geno.tsv
    fi

    if [ ${path_to_vcf_file: -6} = "vcf.gz" ]
      then
        echo "file read as gzipped vcf"
        python $simonhmartin_directory/VCF_processing/parseVCF.py -i $path_to_vcf_file -o $output_directory/$project_name.geno.gz

        pops=$(python $github_directory/parsepops.py $path_to_populations_file)

        python $simonhmartin_directory/freq.py -g $output_directory/$project_name.geno.gz \
        $pops --popsFile $path_to_populations_file -f phased --target derived \
        -o $output_directory/$project_name.geno.tsv
    fi

    if [ ${path_to_vcf_file: -4} =  "geno" ]
      then
        echo "file read as geno"
        pops=$(python $github_directory/parsepops.py $path_to_populations_file)

        python $simonhmartin_directory/freq.py -g $path_to_vcf_file \
        $pops --popsFile $path_to_populations_file -f phased --target derived \
        -o $output_directory/$project_name.geno.tsv
    fi

    if [ ${path_to_vcf_file: -7} = "geno.gz" ]
      then
        echo "file read as gzipped geno"
        pops=$(python $github_directory/parsepops.py $path_to_populations_file)

        python $simonhmartin_directory/freq.py -g $path_to_vcf_file \
        $pops --popsFile $path_to_populations_file -f phased --target derived \
        -o $output_directory/$project_name.geno.tsv
    fi


#it works up to here!

# Quanitfy admixture proportion -- make a new tsv using 3a and 3b pops

    pops_admx=$(python $github_directory/parsepops.py $path_to_admx_populations_file)

    python $simonhmartin_directory/freq.py -g $output_directory/$project_name.geno.gz \
    $pops_admx --popsFile $path_to_admx_populations_file -f phased --target derived \
    -o $output_directory/$project_name.admx.geno.tsv

    Rscript ${github_directory}/ABBAwholegenome.r inputfile_${output_directory}/${project_name}.geno.tsv inputadmx_${output_directory}/${project_name}.admx.geno.tsv outputdirectory_$output_directory simonhmartin_directory_$simonhmartin_directory project_name_$project_name



    if [ ${path_to_vcf_file: -3} = "vcf" ] || [ ${path_to_vcf_file: -6} = "vcf.gz" ]
      then
        python $simonhmartin_directory/ABBABABAwindows.py \
        -g ${output_directory}/${project_name}.geno.gz -f phased \
        -o ${output_directory}/${project_name}_slidingwindows.csv.gz \
        -P1 ${population1} -P2 ${population2} -P3 ${population3} -O ${populationoutgroup} \
        --popsFile ${path_to_populations_file} -w ${windowsize} -m ${minimumsnps} --T ${threads}
    fi

    if [ ${path_to_vcf_file: -4} =  "geno" ] || [ ${path_to_vcf_file: -7} = "geno.gz" ]
      then
        python $simonhmartin_directory/ABBABABAwindows.py \
        -g $path_to_vcf_file -f phased \
        -o ${output_directory}/${project_name}_slidingwindows.csv.gz \
        -P1 ${population1} -P2 ${population2} -P3 ${population3} -O ${populationoutgroup} \
        --popsFile ${path_to_populations_file} -w ${windowsize} -m ${minimumsnps} --T ${threads}
    fi


    python $simonhmartin_directory/ABBABABAwindows.py -g data/${output_directory}/${project_name}.geno.gz -f phased -o -o ${output_directory}/${project_name}.w50m1s10.csv.gz -P1 flo -P2 txn -P3 ama -O slv --popsFile data/bar92.pop.txt -w 50000 -m 1000 -s 20000 --minData 0.5 --T 2
fi
