#!/bin/bash

#A shell script to perform sliding window ABBA BABA tests on whole genome SNP data, with an option to output a list of candidate genes from high regions of introgression. Uses scripts developed by Simon Martin to perform sliding window analyses.

#Program stages are as follows. The first stage is required, but each stage after that is optional. This means that stage 1 and 2 can be performed, or stage 1 and 3, or simply stage 1.

#0) Parses the vcf file into a geno file. Can be skipped by providing the geno file instead of the vcf file.
#1)Performs a genome wide ABBA BABA test. Outputs a block jackknife procedure to estimate standard deviation for D statistic. Outputs the D statistic, the D standard deviation, and the Z score in a text file with the suffix _genomewide_abbababa_stats.txt
#2) Perform sliding window ABBA / BABA analyses that look for regions of introgression between population 2 and population 3 (most strongly sensitive to introgression from P3 -> P2)
#3) Plot as a pdf the points at which introgression has been detected (note that if the input data is in scaffolds of relatively small size, this feature is not very meaningful as the scaffolds are not ordered in a biologically relevant way)
#4) Requires that the reference genome used to generate the VCF file (hereafter refered to as the primary reference genome) does not have its own gff file but instead has been aligned to another genome using the SatsumaSynteny program (hereafter referred to as the secondary reference genome). This stage first identifies windows in which high introgression has been identified, and then identifies regions of the secondary reference genome where the regions of high introgression from the primary reference genome have aligned to. Finally, it looks at the gff file and pulls the names of genes that are within the regions of the secondary genome which have been identified as candidate regions of introgression. This list can easily be input into a gene ontology program, such as DAVID.

#Inputs required for each stage:

#0) Administrative inputs. A name for the project, a location of the directory for the output files of this program, the location of the directory into which the supporting executables have been put (just those from the dannyjackson github), and the location of the directory into which the supporting executables of Simon Martin's github have been put (https://github.com/simonhmartin/genomics_general).

#1) A VCF file, a file assigning individuals (column 1) to populations (column 2), and a defined outgroup

#2) The only thing necessary here is the go ahead from the user to proceed to this step.

#3) The only thing necessary here is the go ahead from the user to proceed to this step.

#4) Here we need two files, one is an output from Satsuma titled satsuma_summary.chained.out and the other is a path to the directory containing the gff files of the genome to which Satsuma aligned the primary reference genome. These should each be labelled as "secondaryreference_chromosome1.gff" "secondaryreference_chromosome2.gff" etc. using chromosome numbers that correspond with those used in the satsuma_summary_chained.out file.

if [ $# -lt 1 ]
  then
    echo "Pulls list of lines from chicken gff file that Satsuma aligned to regions of high introgression, as identified with sliding ABBA BABA windows

    #administrative inputs
    [-a] Name of project; name that all output files will be prefixed with
    [-b] Output directory for files. Must be empty.
    [-c] Path to directory containing all supporing files from github
    [-d] Path to Simon Martin's genomics_general-master directory

    #stage 1
    [-1a] Complete path to input vcf or vcf.gz file, or to .geno or geno.gz file.
    [-1b] File with two tab separated columns. Column 1 has sample names as they appear in the vcf, and column 2 has the population names to which each sample will be attributed. Populations should be listed together, not with some individuals from Pop1 followed by Pop2 followed by more Pop1. IMPORTANTLY the last population will be treated as your outgroup in the ABBA BABA analysis, or your P4, so be sure to list it last.

    #stage 2

    #stage 3
    [-2a] write "y" if you want a plot, write "n" or leave empty if you don't

    #stage 4
    [-3a] Path to satsuma_summary_chained.out file
    [-3b] Path to directory containing gff files"

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
        python $simonhmartin_directory/VCF_processing/parseVCF.py -i $path_to_vcf_file -o $output_directory/$project_name.geno.gz

        pops=$(python $github_directory/parsepops.py $path_to_populations_file)

        python $simonhmartin_directory/freq.py -g $output_directory/$project_name.geno.gz \
        $pops --popsFile $path_to_populations_file -f phased --target derived \
        -o $output_directory/$project_name.geno.tsv
    fi

    if [ ${path_to_vcf_file: -5} = "vcf.gz" ]
      then
        python $simonhmartin_directory/VCF_processing/parseVCF.py -i $path_to_vcf_file -o $output_directory/$project_name.geno.gz

        pops=$(python $github_directory/parsepopsfile.py $path_to_populations_file)

        python $simonhmartin_directory/freq.py -g $output_directory/$project_name.geno.gz \
        $pops --popsFile $path_to_populations_file -f phased --target derived \
        -o $output_directory/$project_name.geno.tsv
    fi

    if [ ${path_to_vcf_file: -4} =  "geno" ]
      then
        pops=$(python $github_directory/parsepopsfile.py $path_to_populations_file)

        python $simonhmartin_directory/freq.py -g $path_to_vcf_file \
        $pops --popsFile $path_to_populations_file -f phased --target derived \
        -o $output_directory/$project_name.geno.tsv
    fi

    if [ ${path_to_vcf_file: -6} = "geno.gz" ]
      then
        pops=$(python $github_directory/parsepopsfile.py $path_to_populations_file)

        python $simonhmartin_directory/freq.py -g $path_to_vcf_file \
        $pops --popsFile $path_to_populations_file -f phased --target derived \
        -o $output_directory/$project_name.geno.tsv
    fi

    echo "a"
    Rscript $github_directory/ABBAwholegenome.r inputfile_$output_directory/$project_name.geno.tsv outputdirectory_$output_directory simonhmartin_directory_$simonhmartin_directory population1_$population1 population2_$population2 population3_$population3 > $output_directory/$project_name_wholegenomestats.txt
    echo "b"
    python $simonhmartin_directory/ABBABABAwindows.py \
    -g ${output_directory}/${project_name}.geno.gz -f phased \
    -o ${output_directory}/${project_name}_slidingwindows.csv.gz \
    -P1 $population1 -P2 $population2 -P3 $population3 -O $populationoutgroup \
    --popsFile $path_to_populations_file -w $windowsize -m $minimumsnps --T $threads
    echo "c"
    Rscript $github_directory/ABBAslidingwindows_plot.r project_name_$project_name outputdirectory_$output_directory

    #I think this next one is writing its output to a weird place. Double check it when the test run finshes.
    python $github_directory/subset_ABBABABAwindows_output.py $output_directory/${project_name}_slidingwindows.csv.gz $output_directory/$project_name $fstat_threshold

    #edit $IDK_OUTPUTOFABOVE to reflect however we output from the subsetting step.
    awk 'BEGIN {FS="\t"}; {print $1 FS $2 FS $3}' $IDK_OUTPUTOFABOVE > $output_directory/$project_name_slidingwindows.bed

    #this isn't quite right... we want to read through each line of the bed file and use that line as the single range to extract...one at a time so we can trace it back to the abba test
    #bedtools getfasta [OPTIONS] -fi $path_to_referencegenome_fasta -bed $output_directory/$project_name_slidingwindows.bed

    while -r file,
    do
      echo $file > $output_directory/$project_name_temp.fsa

      bedtools getfasta -fi $path_to_referencegenome_fasta -bed $output_directory/$project_name_temp.fsa

      blastn -query $path_to_referencegenome_fasta -db "nt" -out $output_directory/$project_name_tempblast.txt

      echo $file >> $output_directory/$project_name_allblast.txt

      cat $output_directory/$project_name_tempblast.txt >> $output_directory/$project_name_allblast.txt

    done < $output_directory/$project_name_slidingwindows.bed



fi
