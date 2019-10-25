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
    g) population1=${OPTARG};;
    h) population2=${OPTARG};;
    i) population3=${OPTARG};;
    j) populationoutgroup=${OPTARG};;
    k) proceed_through_stage_2=${OPTARG};;
    l) windowsize=${OPTARG};;
    m) minimumsnps=${OPTARG};;
    n) threads=${OPTARG};;
    o) fstat_threshold=${OPTARG};;
    p) path_to_satsuma_summary_chained_file=${OPTARG};;
    q) path_to_directory_containing_gff_files=${OPTARG};;
    r) path_to_referencegenome_fasta=${OPTARG};;
    s) path_to_referencegenome_bed=${OPTARG};;


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


    Rscript ${github_directory}/ABBAwholegenome.r inputfile_${output_directory}/${project_name}.geno.tsv outputdirectory_$output_directory simonhmartin_directory_$simonhmartin_directory project_name_$project_name population1_$population1 population2_$population2 population3_$population3

    #it works up to here!

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

    Rscript $github_directory/ABBAslidingwindows_plot.r project_name_$project_name outputdirectory_$output_directory

    python $github_directory/subset_ABBABABAwindows_output.py ${output_directory}/${project_name} $fstat_threshold

    python $github_directory/subset_satsumachain_by_ABBA.py ${output_directory}/${project_name} ${path_to_satsuma_summary_chained_file}
fi
