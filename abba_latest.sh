#!/bin/bash

#A shell script to perform sliding window ABBA BABA tests on whole genome SNP data, with an option to output a list of candidate genes from high regions of introgression. Uses scripts developed by Simon Martin to perform sliding window analyses.

#Program stages are as follows. The first stage is required, but each stage after that is optional. This means that stage 1 and 2 can be performed, or stage 1 and 3, or simply stage 1.

#1) Perform sliding window ABBA / BABA analyses that look for regions of introgression between population 2 and population 3 (most strongly sensitive to introgression from P3 -> P2)
#2) Plot as a pdf the points at which introgression has been detected (note that if the input data is in scaffolds of relatively small size, this feature is not very meaningful as the scaffolds are not ordered in a biologically relevant way)
#3) Requires that the reference genome used to generate the VCF file (hereafter refered to as the primary reference genome) does not have its own gff file but instead has been aligned to another genome using the SatsumaSynteny program (hereafter referred to as the secondary reference genome). This stage first identifies windows in which high introgression has been identified, and then identifies regions of the secondary reference genome where the regions of high introgression from the primary reference genome have aligned to. Finally, it looks at the gff file and pulls the names of genes that are within the regions of the secondary genome which have been identified as candidate regions of introgression. This list can easily be input into a gene ontology program, such as DAVID.

#Inputs required for each stage:

#0) Administrative inputs. A name for the project, a location of the directory for the output files of this program, the location of the directory into which the supporting executables have been put (just those from the dannyjackson github), and the location of the directory into which the supporting executables of Simon Martin's github have been put (https://github.com/simonhmartin/genomics_general).

#1) A VCF file, a file assigning individuals (column 1) to populations (column 2), and a defined outgroup

#2) The only thing necessary here is the go ahead from the user to proceed to this step.

#3) Here we need two files, one is an output from Satsuma titled satsuma_summary.chained.out and the other is a path to the directory containing the gff files of the genome to which Satsuma aligned the primary reference genome. These should each be labelled as "secondaryreference_chromosome1.gff" "secondaryreference_chromosome2.gff" etc. using chromosome numbers that correspond with those used in the satsuma_summary_chained.out file.

if [ $# -lt 1 ]
  then
    echo "Pulls list of lines from chicken gff file that Satsuma aligned to regions of high introgression, as identified with sliding ABBA BABA windows

    #administrative inputs
    [-a] Name of project; name that all output files will be prefixed with
    [-b] Output directory for files. Must be empty.
    [-c] Path to directory containing all supporing files from github
    [-d] Path to Simon Martin's genomics_general-master directory

    #stage 1
    [-1a] Complete path to input vcf file
    [-1b] File with two tab separated columns. Column 1 has sample names as they appear in the vcf, and column 2 has the population names to which each sample will be attributed. Populations should be listed together, not with some individuals from Pop1 followed by Pop2 followed by more Pop1. IMPORTANTLY the last population will be treated as your outgroup in the ABBA BABA analysis, or your P4, so be sure to list it last.
    [-1c] The outgroup population

    #stage 2
    [-2a] write "y" if you want a plot, write "n" or leave empty if you don't

    #stage 3
    [-3a] Path to satsuma_summary_chained.out file
    [-3b] Path to directory containing gff files"

  else
    while getopts a:b:c:d:1a:1b:1c:2a:3a:3b option
    do
    case "${option}"
    in
    a) project_name=${OPTARG};;
    b) output_directory=${OPTARG};;
    c) github_directory=${OPTARG};;
    d) simonhmartin_directory=${OPTARG};;
    1a) path_to_vcf_file=${OPTARG};;
    1b) path_to_populations_file=${OPTARG};;
    1c) outgroup=${OPTARG};;
    2a) proceed_through_stage_2=${OPTARG};;
    3a) path_to_satsuma_summary_chained_file=${OPTARG};;
    3b) path_to_directory_containing_gff_files=${OPTARG};;

    esac
    done

    #this converts the vcf file into a format that Simon Martin refers to as a ".geno", see his website for more information https://github.com/simonhmartin/genomics_general/tree/master/VCF_processing

    python $simonhmartin_directory/VCF_processing/parseVCF.py -i $path_to_vcf_file -o $output_directory/$project_name.geno.gz
