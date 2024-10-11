# Command line scripts for all analyses in the paper, "Evidence from whole genomes for interspecific introgression and pantropical gene flow in boobies, a genus of tropical seabirds"
=
The files trim-and-QC.sh and align-and-sort.sh can be found here: https://github.com/erikrfunk/whole_genome_bioinformatics/ (thanks Erik for the useful scripts!)

    ~/whole_genome_bioinformatics/trim-and-QC.sh -i /data5/sulafilenames.txt -p /data5/ -f R1_001.fastq.gz -r R2_001.fastq.gz -a /data5/NexteraPE-PE.fa -t 12

Important settings within the script: ILLUMINACLIP:$adapters:1:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:90

    align-and-sort.sh -i /data5/sulidae/reference_lists/sulafilenames.txt -r ~/reference_datasets/flightless/ncbi-genomes-2020-08-14/GCA_002173475.1_Pharrisi_ref_V1/GCA_002173475.1_Pharrisi_ref_V1_genomic.fna -t 12 -p /data5/sulidae/my_datasets/trimming_step/fastqs/ -b bam_files_flightless -s sorted_bam_files_flightless

    ref="~/reference_datasets/flightless/ncbi-genomes-2020-08-14/GCA_002173475.1_Pharrisi_ref_V1/GCA_002173475.1_Pharrisi_ref_V1_genomic.fna"
    bamdir="/data5/sulidae/my_datasets/sorted_bam_files/"
    ID="sula_flightless"
    /usr/bin/bcftools mpileup -Ou -f "$ref" -a FORMAT/AD,DP,INFO/AD,SP "$bamdir"*_sorted_RGadded_dupmarked.bam | /usr/bin/bcftools call -mv -V indels > "$ID"_snps_multiallelic.vcf

    awk '{print $1"_"$2}' sula_flightless.vcf | grep -v '#' > locinames.txt
    grep -v '#' sula_flightless.vcf > temp.vcf
    awk -F "\t" 'FNR==NR{a[NR]=$1;next}{$3=a[FNR]}1' locinames.txt temp.vcf > newfile.vcf
    tr ' ' '\t' < newfile.vcf 1<> newfile.vcf
    grep '#' sula_flightless.vcf > sula_flightless_ID.vcf
    cat newfile.vcf >> sula_flightless_ID.vcf

# generate summary statistics

    ~/plink --vcf /data5/sulidae/my_datasets/sula_flightless_ID.vcf --allow-extra-chr --missing --cluster-missing --freq

    ~/plink --vcf /data5/sulidae/my_datasets/sula_flightless_ID.vcf --allow-extra-chr --missing --cluster missing --within cluster_species.txt --freq

    ~/plink --vcf /data5/sulidae/my_datasets/sula_flightless_ID.vcf --allow-extra-chr --missing --cluster missing --within cluster_pop.txt --freq

# filter the vcf

    #filters by quality
    bcftools view -i 'QUAL>100' sula_flightless_ID.vcf > sula_flightless_qualitysort.vcf

    #filters by depth and removes indels
    vcftools --vcf sula_flightless_qualitysort.vcf --min-meanDP 2 --max-meanDP 8 --remove-indels --recode --out sula_flightless_filtered

    ~/plink --vcf /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.02 --mind 0.2 --maf 0.01 --recode vcf-iid --out sula_flightless_flightless_mind2

    ~/plink --vcf /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.02 --mind 0.1 --maf 0.01 --recode vcf-iid --out sula_flightless_cleaned


    cp sula_flightless_cleaned.vcf sula_flightless_cleaned_zip.vcf
    bgzip sula_flightless_cleaned_zip.vcf

    bcftools index sula_flightless_cleaned_zip.vcf.gz


# RAxML

~/vcf2phylip/vcf2phylip.py -i sula_flightless_flightless_mind2.vcf

mkdir pruned
cd pruned
python ~/sula/filter_invariants_all.py ../sula_flightless_flightless_mind2.min4.phy
mv variantsites.phy ../variantsites_mind2.phy
mv variantsites_kept.txt ../variantsites_mind2_kept.txt
cd ..
rm -r pruned


raxmlHPC -m ASC_GTRCAT --asc-corr felsenstein -f d -d -k -n sula_flightless_b -q /data5/sulidae/final/to_flightless/raxml/partition_file.txt -s /data5/sulidae/my_datasets/variantsites_mind2.phy -T 6 -p 12345 -N 10 ­-b 12345 -V




# Examine structure with PCA

# All

echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BRBO" >> pops.txt
echo -e "BRBO" >> pops.txt
echo -e "BRBO" >> pops.txt
echo -e "BRBO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "RFBO" >> pops.txt
echo -e "RFBO" >> pops.txt
echo -e "RFBO" >> pops.txt
echo -e "RFBO" >> pops.txt
echo -e "RFBO" >> pops.txt
echo -e "RFBO" >> pops.txt


~/genomics/PCA_r.sh -v /data5/sulidae/my_datasets/sula_flightless_flightless_mind2.vcf -o /data5/sulidae/final/to_flightless/PCA/all/ -p /data5/sulidae/final/to_flightless/PCA/all/pops.txt -n all -s y

#four
bcftools view -s MABO301,MABO302,MABO304,MABO305,MABO306,NABO401,NABO402,NABO403,NABO404,NABO405,NABO406,BFBO501,BFBO502,BFBO503,BFBO504,BFBO505,BFBO506,PEBO601,PEBO602,PEBO603,PEBO604,PEBO605,PEBO606 /data5/sulidae/my_datasets/sula_flightless_flightless_mind2.vcf --force-samples > /data5/sulidae/final/to_flightless/four_pca.vcf


echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt


~/genomics/PCA_r.sh -v /data5/sulidae/final/to_flightless/four_pca.vcf -o /data5/sulidae/final/to_flightless/PCA/four/ -p /data5/sulidae/final/to_flightless/PCA/four/pops.txt -n four -s y

#sister_mn

bcftools view -s MABO301,MABO302,MABO304,MABO305,MABO306,NABO401,NABO402,NABO403,NABO404,NABO405,NABO406 /data5/sulidae/my_datasets/sula_flightless_flightless_mind2.vcf --force-samples > /data5/sulidae/final/to_flightless/mana_pca.vcf

echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt

echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt



~/genomics/PCA_r.sh -v /data5/sulidae/final/to_flightless/mana_pca.vcf -o /data5/sulidae/final/to_flightless/PCA/sister_mn/ -p /data5/sulidae/final/to_flightless/PCA/sister_mn/pops.txt -n sister_mn -s y


#sister_bp

bcftools view -s BFBO501,BFBO502,BFBO503,BFBO504,BFBO505,BFBO506,PEBO601,PEBO602,PEBO603,PEBO604,PEBO605,PEBO606 /data5/sulidae/my_datasets/sula_flightless_flightless_mind2.vcf --force-samples > /data5/sulidae/final/to_flightless/PCA/sister_bp/bfpe_pca.vcf

echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt

echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt

~/genomics/PCA_r.sh -v /data5/sulidae/final/to_flightless/PCA/sister_bp/bfpe_pca.vcf -o /data5/sulidae/final/to_flightless/PCA/sister_bp/ -p /data5/sulidae/final/to_flightless/PCA/sister_bp/pops.txt -n sister_bp -s y


#species
mkdir bfbo pebo mabo nabo brbo rfbo

#bfbo

bcftools view -s BFBO501,BFBO502,BFBO503,BFBO504,BFBO505,BFBO506 /data5/sulidae/my_datasets/sula_flightless_flightless_mind2.vcf --force-samples > /data5/sulidae/final/to_flightless/PCA/bfbo/bfbo_pca.vcf

echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt


~/genomics/PCA_r.sh -v /data5/sulidae/final/to_flightless/PCA/bfbo/bfbo_pca.vcf -o /data5/sulidae/final/to_flightless/PCA/bfbo/ -p /data5/sulidae/final/to_flightless/PCA/bfbo/pops.txt -n bfbo -s y

#pebo

bcftools view -s PEBO601,PEBO602,PEBO603,PEBO604,PEBO605,PEBO606 /data5/sulidae/my_datasets/sula_flightless_flightless_mind2.vcf --force-samples > /data5/sulidae/final/to_flightless/PCA/pebo/pebo_pca.vcf

echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt


~/genomics/PCA_r.sh -v /data5/sulidae/final/to_flightless/PCA/pebo/pebo_pca.vcf -o /data5/sulidae/final/to_flightless/PCA/pebo/ -p /data5/sulidae/final/to_flightless/PCA/pebo/pops.txt -n pebo -s y


#### resume ####
#mabo

bcftools view -s MABO301,MABO302,MABO304,MABO305,MABO306 /data5/sulidae/my_datasets/sula_flightless_flightless_mind2.vcf --force-samples > /data5/sulidae/final/to_flightless/PCA/mabo/mabo_pca.vcf

echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt


~/genomics/PCA_r.sh -v /data5/sulidae/final/to_flightless/PCA/mabo/mabo_pca.vcf -o /data5/sulidae/final/to_flightless/PCA/mabo/ -p /data5/sulidae/final/to_flightless/PCA/mabo/pops.txt -n mabo -s y

#nabo

bcftools view -s 	NABO401,NABO402,NABO403,NABO404,NABO405,NABO406 /data5/sulidae/my_datasets/sula_flightless_flightless_mind2.vcf --force-samples > /data5/sulidae/final/to_flightless/PCA/nabo/nabo_pca.vcf

echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt


~/genomics/PCA_r.sh -v /data5/sulidae/final/to_flightless/PCA/nabo/nabo_pca.vcf -o /data5/sulidae/final/to_flightless/PCA/nabo/ -p /data5/sulidae/final/to_flightless/PCA/nabo/pops.txt -n nabo -s y

#brbo

bcftools view -s BRBO201,BRBO202,BRBO203,BRBO204,BRBO205 /data5/sulidae/my_datasets/sula_flightless_flightless_mind2.vcf --force-samples > /data5/sulidae/final/to_flightless/PCA/brbo/brbo_pca.vcf

echo -e "BRBO" >> pops.txt
echo -e "BRBO" >> pops.txt
echo -e "BRBO" >> pops.txt
echo -e "BRBO" >> pops.txt
echo -e "BRBO" >> pops.txt


~/genomics/PCA_r.sh -v /data5/sulidae/final/to_flightless/PCA/brbo/brbo_pca.vcf -o /data5/sulidae/final/to_flightless/PCA/brbo/ -p /data5/sulidae/final/to_flightless/PCA/brbo/pops.txt -n brbo -s y

#rfbo

bcftools view -s RFBO102,RFBO103,RFBO104,RFBO105,RFBO106,RFBO101 /data5/sulidae/my_datasets/sula_flightless_flightless_mind2.vcf --force-samples > /data5/sulidae/final/to_flightless/PCA/rfbo/rfbo_pca.vcf

echo -e "RFBO" >> pops.txt
echo -e "RFBO" >> pops.txt
echo -e "RFBO" >> pops.txt
echo -e "RFBO" >> pops.txt
echo -e "RFBO" >> pops.txt
echo -e "RFBO" >> pops.txt


~/genomics/PCA_r.sh -v /data5/sulidae/final/to_flightless/PCA/rfbo/rfbo_pca.vcf -o /data5/sulidae/final/to_flightless/PCA/rfbo/ -p /data5/sulidae/final/to_flightless/PCA/rfbo/pops.txt -n rfbo -s y










# Pruned individuals PCA

mkdir all  bfbo  brbo  four  mabo  nabo  pebo  rfbo  sister_bp  sister_mn
# All

bcftools view -s MABO301,MABO302,MABO304,MABO305,MABO306,NABO401,NABO402,NABO403,NABO404,NABO405,NABO406,BFBO501,BFBO502,BFBO503,BFBO504,BFBO505,BFBO506,BRBO201,BRBO202,BRBO203,BRBO204,BRBO205,PEBO601,PEBO602,PEBO603,PEBO604,PEBO605,PEBO606,RFBO102,RFBO103,RFBO104,RFBO105,RFBO106,RFBO101 /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf > /data5/sulidae/final/all_pruned_pca.vcf

echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BRBO" >> pops.txt
echo -e "BRBO" >> pops.txt
echo -e "BRBO" >> pops.txt
echo -e "BRBO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "RFBO" >> pops.txt
echo -e "RFBO" >> pops.txt
echo -e "RFBO" >> pops.txt
echo -e "RFBO" >> pops.txt
echo -e "RFBO" >> pops.txt
echo -e "RFBO" >> pops.txt


~/genomics/PCA_r.sh -v /data5/sulidae/final/all_pruned_pca.vcf -o /data5/sulidae/final/PCA/all/ -p /data5/sulidae/final/PCA/all/pops.txt -n all -s y

#four
bcftools view -s MABO301,MABO302,MABO304,MABO305,MABO306,NABO402,NABO403,NABO404,NABO405,NABO406,BFBO501,BFBO502,BFBO503,BFBO504,BFBO505,PEBO601,PEBO603,PEBO604,PEBO605,PEBO606 /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf > /data5/sulidae/final/four_pruned_pca.vcf


echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt


~/genomics/PCA_r.sh -v /data5/sulidae/final/four_pruned_pca.vcf -o /data5/sulidae/final/PCA/four/ -p /data5/sulidae/final/PCA/four/pops.txt -n four -s y

#sister_mn

bcftools view -s MABO301,MABO302,MABO304,MABO305,MABO306,NABO402,NABO403,NABO404,NABO405,NABO406 /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf > /data5/sulidae/final/mana_pruned_pca.vcf

echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt

echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt



~/genomics/PCA_r.sh -v /data5/sulidae/final/mana_pruned_pca.vcf -o /data5/sulidae/final/PCA/sister_mn/ -p /data5/sulidae/final/PCA/sister_mn/pops.txt -n sister_mn -s y


#sister_bp

bcftools view -s BFBO501,BFBO502,BFBO503,BFBO504,BFBO505,PEBO601,PEBO603,PEBO604,PEBO605,PEBO606 /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf > /data5/sulidae/final/PCA/sister_bp/bfpe_pruned_pca.vcf

echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt

echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt

~/genomics/PCA_r.sh -v /data5/sulidae/final/PCA/sister_bp/bfpe_pruned_pca.vcf -o /data5/sulidae/final/PCA/sister_bp/ -p /data5/sulidae/final/PCA/sister_bp/pops.txt -n sister_bp -s y


#species
mkdir bfbo pebo mabo nabo brbo rfbo

#bfbo

bcftools view -s BFBO501,BFBO502,BFBO503,BFBO504,BFBO505 /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf > /data5/sulidae/final/PCA/bfbo/bfbo_pruned_pca.vcf

echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt
echo -e "BFBO" >> pops.txt


~/genomics/PCA_r.sh -v /data5/sulidae/final/PCA/bfbo/bfbo_pruned_pca.vcf -o /data5/sulidae/final/PCA/bfbo/ -p /data5/sulidae/final/PCA/bfbo/pops.txt -n bfbo -s y

#pebo

bcftools view -s PEBO601,PEBO603,PEBO604,PEBO605,PEBO606 /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf > /data5/sulidae/final/PCA/pebo/pebo_pruned_pca.vcf

echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt
echo -e "PEBO" >> pops.txt


~/genomics/PCA_r.sh -v /data5/sulidae/final/PCA/pebo/pebo_pruned_pca.vcf -o /data5/sulidae/final/PCA/pebo/ -p /data5/sulidae/final/PCA/pebo/pops.txt -n pebo -s y

#mabo

bcftools view -s MABO301,MABO302,MABO304,MABO305,MABO306 /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf > /data5/sulidae/final/PCA/mabo/mabo_pca.vcf

echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt
echo -e "MABO" >> pops.txt


~/genomics/PCA_r.sh -v /data5/sulidae/final/PCA/mabo/mabo_pca.vcf -o /data5/sulidae/final/PCA/mabo/ -p /data5/sulidae/final/PCA/mabo/pops.txt -n mabo -s y

#nabo

bcftools view -s 	NABO402,NABO403,NABO404,NABO405,NABO406 /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf > /data5/sulidae/final/PCA/nabo/nabo_pruned_pca.vcf

echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt
echo -e "NABO" >> pops.txt


~/genomics/PCA_r.sh -v /data5/sulidae/final/PCA/nabo/nabo_pruned_pca.vcf -o /data5/sulidae/final/PCA/nabo/ -p /data5/sulidae/final/PCA/nabo/pops.txt -n nabo -s y

#brbo

bcftools view -s BRBO201,BRBO202,BRBO203,BRBO205 /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf > /data5/sulidae/final/PCA/brbo/brbo_pruned_pca.vcf

echo -e "BRBO" >> pops.txt
echo -e "BRBO" >> pops.txt
echo -e "BRBO" >> pops.txt
echo -e "BRBO" >> pops.txt


~/genomics/PCA_r.sh -v /data5/sulidae/final/PCA/brbo/brbo_pruned_pca.vcf -o /data5/sulidae/final/PCA/brbo/ -p /data5/sulidae/final/PCA/brbo/pops.txt -n brbo -s y

#rfbo

bcftools view -s RFBO102,RFBO103,RFBO104,RFBO105,RFBO106,RFBO101 /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf > /data5/sulidae/final/PCA/rfbo/rfbo_pca.vcf

echo -e "RFBO" >> pops.txt
echo -e "RFBO" >> pops.txt
echo -e "RFBO" >> pops.txt
echo -e "RFBO" >> pops.txt
echo -e "RFBO" >> pops.txt
echo -e "RFBO" >> pops.txt


~/genomics/PCA_r.sh -v /data5/sulidae/final/PCA/rfbo/rfbo_pca.vcf -o /data5/sulidae/final/PCA/rfbo/ -p /data5/sulidae/final/PCA/rfbo/pops.txt -n rfbo -s y


# ABBA BABA

# Tests of introgression between species

# Convert vcf to geno
python ~/genomics_general-master/VCF_processing/parseVCF.py -i /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf -o /data5/sulidae/final/to_flightless/abba/sula.geno.gz

~/sula/abba.sh -a [NAME OF RUN] -b [PATH TO OUTPUT] -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/abba/sula.geno.gz -f /data5/sulidae/final/abba/NAMEOFRUN/pops.txt -g P1 -h P2 -i P3 -j P4 -k y  -n 2

# BF_GCA    BF_LOBOS      PE  RF

echo -e "BFBO504\tP1" > pops.txt
echo -e "BFBO502\tP2" >> pops.txt
echo -e "BFBO503\tP2" >> pops.txt
echo -e "PEBO601\tP3" >> pops.txt
echo -e "PEBO603\tP3" >> pops.txt
echo -e "PEBO604\tP3" >> pops.txt
echo -e "PEBO605\tP3" >> pops.txt
echo -e "PEBO606\tP3" >> pops.txt
echo -e "RFBO101\tP4" >> pops.txt
echo -e "RFBO102\tP4" >> pops.txt
echo -e "RFBO103\tP4" >> pops.txt
echo -e "RFBO104\tP4" >> pops.txt
echo -e "RFBO105\tP4" >> pops.txt
echo -e "RFBO106\tP4" >> pops.txt


~/sula/abba.sh -a BfBfPeRf_500kb -b /data5/sulidae/final/to_flightless/abba/BfBfPeRf/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/to_flightless/abba/sula.geno.gz -f /data5/sulidae/final/to_flightless/abba/BfBfPeRf/pops.txt -k y -n 2


# MA_ATLCAR MA_INDOPA     NA  RF

echo -e "MABO304\tP1" > pops.txt
echo -e "MABO305\tP1" >> pops.txt
echo -e "MABO302\tP2" >> pops.txt
echo -e "MABO306\tP2" >> pops.txt
echo -e "NABO402\tP3" >> pops.txt
echo -e "NABO403\tP3" >> pops.txt
echo -e "NABO404\tP3" >> pops.txt
echo -e "NABO405\tP3" >> pops.txt
echo -e "NABO406\tP3" >> pops.txt
echo -e "RFBO101\tP4" >> pops.txt
echo -e "RFBO102\tP4" >> pops.txt
echo -e "RFBO103\tP4" >> pops.txt
echo -e "RFBO104\tP4" >> pops.txt
echo -e "RFBO105\tP4" >> pops.txt
echo -e "RFBO106\tP4" >> pops.txt

~/sula/abba.sh -a MaMaNaRf_500kb -b /data5/sulidae/final/to_flightless/abba/MaMaNaRf/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/to_flightless/abba/sula.geno.gz -f /data5/sulidae/final/to_flightless/abba/MaMaNaRf/pops.txt -k y -n 2

# PE  BF  NA  RF

echo -e "PEBO601\tP1" > pops.txt
echo -e "PEBO603\tP1" >> pops.txt
echo -e "PEBO604\tP1" >> pops.txt
echo -e "PEBO605\tP1" >> pops.txt
echo -e "PEBO606\tP1" >> pops.txt
echo -e "BFBO501\tP2" >> pops.txt
echo -e "BFBO502\tP2" >> pops.txt
echo -e "BFBO503\tP2" >> pops.txt
echo -e "BFBO504\tP2" >> pops.txt
echo -e "BFBO505\tP2" >> pops.txt

echo -e "NABO401\tP3" >> pops.txt
echo -e "NABO402\tP3" >> pops.txt
echo -e "NABO403\tP3" >> pops.txt
echo -e "NABO404\tP3" >> pops.txt
echo -e "NABO405\tP3" >> pops.txt
echo -e "NABO406\tP3" >> pops.txt

echo -e "RFBO101\tP4" >> pops.txt
echo -e "RFBO102\tP4" >> pops.txt
echo -e "RFBO103\tP4" >> pops.txt
echo -e "RFBO104\tP4" >> pops.txt
echo -e "RFBO105\tP4" >> pops.txt
echo -e "RFBO106\tP4" >> pops.txt



~/sula/abba.sh -a PeBfNaRf -b /data5/sulidae/final/to_flightless/abba/PeBfNaRf/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/to_flightless/abba/sula.geno.gz -f /data5/sulidae/final/to_flightless/abba/PeBfNaRf/pops.txt -k y -n 2


../

# PE BF BR RF
mkdir PeBfBrRf

cd PeBfBrRf

echo -e "PEBO601\tP1" > pops.txt
echo -e "PEBO602\tP1" >> pops.txt
echo -e "PEBO603\tP1" >> pops.txt
echo -e "PEBO604\tP1" >> pops.txt
echo -e "PEBO605\tP1" >> pops.txt
echo -e "PEBO606\tP1" >> pops.txt
echo -e "BFBO501\tP2" >> pops.txt
echo -e "BFBO502\tP2" >> pops.txt
echo -e "BFBO503\tP2" >> pops.txt
echo -e "BFBO504\tP2" >> pops.txt
echo -e "BFBO505\tP2" >> pops.txt
echo -e "BFBO506\tP2" >> pops.txt
echo -e "BRBO201\tP3" >> pops.txt
echo -e "BRBO202\tP3" >> pops.txt
echo -e "BRBO205\tP3" >> pops.txt
echo -e "RFBO101\tP4" >> pops.txt
echo -e "RFBO102\tP4" >> pops.txt
echo -e "RFBO103\tP4" >> pops.txt
echo -e "RFBO104\tP4" >> pops.txt
echo -e "RFBO105\tP4" >> pops.txt
echo -e "RFBO106\tP4" >> pops.txt



~/sula/abba.sh -a PeBfBrRf -b /data5/sulidae/final/to_flightless/abba/PeBfBrRf/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/to_flightless/abba/sula.geno.gz -f /data5/sulidae/final/to_flightless/abba/PeBfBrRf/pops.txt -k y -n 2

cd ../

# BF_GCA    BF_GALAPAGOS  NA  RF
mkdir BfBfNaRf
cd BfBfNaRf

echo -e "BFBO504\tP1" > pops.txt
echo -e "BFBO502\tP2" >> pops.txt
echo -e "BFBO503\tP2" >> pops.txt
echo -e "NABO402\tP3" >> pops.txt
echo -e "NABO403\tP3" >> pops.txt
echo -e "NABO404\tP3" >> pops.txt
echo -e "NABO405\tP3" >> pops.txt
echo -e "NABO406\tP3" >> pops.txt
echo -e "RFBO101\tP4" >> pops.txt
echo -e "RFBO102\tP4" >> pops.txt
echo -e "RFBO103\tP4" >> pops.txt
echo -e "RFBO104\tP4" >> pops.txt
echo -e "RFBO105\tP4" >> pops.txt
echo -e "RFBO106\tP4" >> pops.txt


~/sula/abba.sh -a BfBfNaRf_500kb -b /data5/sulidae/final/to_flightless/abba/BfBfNaRf/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/to_flightless/abba/sula.geno.gz -f /data5/sulidae/final/to_flightless/abba/BfBfNaRf/pops.txt -k y -n 2

~/sula/abba.sh -a BfBfNaRf_10kb -b /data5/sulidae/final/abba/BfBfNaRf/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/abba/sula.geno.gz -f /data5/sulidae/final/abba/BfBfNaRf/pops.txt -k y -n 2

~/sula/abba.sh -a BfBfNaRf_5kb -b /data5/sulidae/final/abba/BfBfNaRf/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/abba/sula.geno.gz -f /data5/sulidae/final/abba/BfBfNaRf/pops.txt -k y -n 2



# Tests of introgression within species


# rfbo

echo -e "RFBO101\tR1" >> pops.txt
echo -e "RFBO102\tR2" >> pops.txt
echo -e "RFBO103\tR3" >> pops.txt
echo -e "RFBO104\tR4" >> pops.txt
echo -e "RFBO105\tR5" >> pops.txt
echo -e "RFBO106\tR6" >> pops.txt

echo -e "MABO304\tO" >> pops.txt
echo -e "MABO305\tO" >> pops.txt
echo -e "MABO302\tO" >> pops.txt
echo -e "MABO306\tO" >> pops.txt

echo -e "P1\tP2\tP3" >> abbapops.txt
echo -e "R1\tR2\tR3" >> abbapops.txt
echo -e "R1\tR2\tR5" >> abbapops.txt
echo -e "R1\tR2\tR4" >> abbapops.txt
echo -e "R1\tR2\tR6" >> abbapops.txt
echo -e "R1\tR3\tR6" >> abbapops.txt
echo -e "R1\tR3\tR5" >> abbapops.txt
echo -e "R1\tR3\tR4" >> abbapops.txt
echo -e "R2\tR3\tR6" >> abbapops.txt
echo -e "R2\tR3\tR5" >> abbapops.txt
echo -e "R2\tR3\tR4" >> abbapops.txt
echo -e "R2\tR6\tR5" >> abbapops.txt
echo -e "R2\tR6\tR4" >> abbapops.txt
echo -e "R3\tR5\tR4" >> abbapops.txt
echo -e "R3\tR5\tR6" >> abbapops.txt
echo -e "R1\tR5\tR4" >> abbapops.txt
echo -e "R2\tR5\tR4" >> abbapops.txt
echo -e "R3\tR5\tR4" >> abbapops.txt
echo -e "R1\tR4\tR6" >> abbapops.txt
echo -e "R2\tR4\tR6" >> abbapops.txt
echo -e "R3\tR4\tR6" >> abbapops.txt
echo -e "R5\tR4\tR6" >> abbapops.txt

pops=$(python ~/sula/parsepops.py /data5/sulidae/final/to_flightless/abba/rfbo/pops.txt)

python ~/genomics_general-master/freq.py -g /data5/sulidae/final/to_flightless/abba/sula.geno.gz \
$pops --popsFile /data5/sulidae/final/to_flightless/abba/rfbo/pops.txt -f phased --target derived \
-o /data5/sulidae/final/to_flightless/abba/rfbo/rfbo.geno.tsv

Rscript ~/sula/ABBA_RFBO.r
Rscript ~/sula/ABBA_RFBO_2.r
Rscript ~/sula/ABBA_RFBO_3.r



# mabo

echo -e "MABO301\tM1" >> pops.txt
echo -e "MABO302\tM2" >> pops.txt
echo -e "MABO304\tM4" >> pops.txt
echo -e "MABO305\tM5" >> pops.txt
echo -e "MABO306\tM6" >> pops.txt

echo -e "RFBO101\tO" >> pops.txt
echo -e "RFBO102\tO" >> pops.txt
echo -e "RFBO103\tO" >> pops.txt
echo -e "RFBO104\tO" >> pops.txt
echo -e "RFBO105\tO" >> pops.txt
echo -e "RFBO106\tO" >> pops.txt



echo -e "MABO301\tM1" >> pops_45.txt
echo -e "MABO302\tM2" >> pops_45.txt
echo -e "MABO306\tM6" >> pops_45.txt
echo -e "MABO304\tM45" >> pops_45.txt
echo -e "MABO305\tM45" >> pops_45.txt



echo -e "RFBO101\tO" >> pops_45.txt
echo -e "RFBO102\tO" >> pops_45.txt
echo -e "RFBO103\tO" >> pops_45.txt
echo -e "RFBO104\tO" >> pops_45.txt
echo -e "RFBO105\tO" >> pops_45.txt
echo -e "RFBO106\tO" >> pops_45.txt


pops=$(python ~/sula/parsepops.py /data5/sulidae/final/to_flightless/abba/mabo/pops.txt)

python ~/genomics_general-master/freq.py -g /data5/sulidae/final/to_flightless/abba/sula.geno.gz \
$pops --popsFile /data5/sulidae/final/to_flightless/abba/mabo/pops.txt -f phased --target derived \
-o /data5/sulidae/final/to_flightless/abba/mabo/mabo.geno.tsv

pops=$(python ~/sula/parsepops.py /data5/sulidae/final/to_flightless/abba/mabo/pops_45.txt)

python ~/genomics_general-master/freq.py -g /data5/sulidae/final/to_flightless/abba/sula.geno.gz \
$pops --popsFile /data5/sulidae/final/to_flightless/abba/mabo/pops_45.txt -f phased --target derived \
-o /data5/sulidae/final/to_flightless/abba/mabo/mabo_45.geno.tsv

Rscript ~/sula/ABBA_MABO_1.r
Rscript ~/sula/ABBA_MABO_2.r
Rscript ~/sula/ABBA_MABO_3.r
Rscript ~/sula/ABBA_MABO_4.r



M4  M5  M162

echo -e "MABO304\tP1" > pops.txt
echo -e "MABO305\tP2" >> pops.txt
echo -e "MABO301\tP3" >> pops.txt
echo -e "MABO302\tP3" >> pops.txt
echo -e "MABO306\tP3" >> pops.txt
echo -e "RFBO101\tP4" >> pops.txt
echo -e "RFBO102\tP4" >> pops.txt
echo -e "RFBO103\tP4" >> pops.txt
echo -e "RFBO104\tP4" >> pops.txt
echo -e "RFBO105\tP4" >> pops.txt
echo -e "RFBO106\tP4" >> pops.txt



~/sula/abba.sh -a mabo_4_5_162 -b /data5/sulidae/final/to_flightless/abba/mabo/mabo_4_5_162/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/to_flightless/abba/sula.geno.gz -f /data5/sulidae/final/to_flightless/abba/mabo/mabo_4_5_162/pops.txt -k y -n 2

M1,6  M2  M45

echo -e "MABO301\tP1" > pops.txt
echo -e "MABO306\tP1" >> pops.txt
echo -e "MABO302\tP2" >> pops.txt
echo -e "MABO304\tP3" >> pops.txt
echo -e "MABO305\tP3" >> pops.txt
echo -e "RFBO101\tP4" >> pops.txt
echo -e "RFBO102\tP4" >> pops.txt
echo -e "RFBO103\tP4" >> pops.txt
echo -e "RFBO104\tP4" >> pops.txt
echo -e "RFBO105\tP4" >> pops.txt
echo -e "RFBO106\tP4" >> pops.txt



~/sula/abba.sh -a mabo_16_2_45 -b /data5/sulidae/final/to_flightless/abba/mabo/mabo_16_2_45/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/abba/sula.geno.gz -f /data5/sulidae/final/to_flightless/abba/mabo/mabo_16_2_45/pops.txt -k y -n 2


# brbo

# all

echo -e "BRBO201\tB1" > pops.txt
echo -e "BRBO202\tB2" >> pops.txt
echo -e "BRBO203\tB3" >> pops.txt
echo -e "BRBO204\tB4" >> pops.txt
echo -e "BRBO205\tB5" >> pops.txt

echo -e "RFBO101\tO" >> pops.txt
echo -e "RFBO102\tO" >> pops.txt
echo -e "RFBO103\tO" >> pops.txt
echo -e "RFBO104\tO" >> pops.txt
echo -e "RFBO105\tO" >> pops.txt
echo -e "RFBO106\tO" >> pops.txt



pops=$(python ~/sula/parsepops.py /data5/sulidae/final/to_flightless/abba/brbo/pops.txt)

python ~/genomics_general-master/freq.py -g /data5/sulidae/final/to_flightless/abba/sula.geno.gz \
$pops --popsFile /data5/sulidae/final/to_flightless/abba/brbo/pops.txt -f phased --target derived \
-o /data5/sulidae/final/to_flightless/abba/brbo/brbo_all.geno.tsv

Rscript ~/sula/ABBA_BRBO_1.r

# 35
echo 'brbo 35' >> brbo.abbawholegenome.stats.txt

echo -e "BRBO201\tB1" >> pops_35.txt
echo -e "BRBO202\tB2" >> pops_35.txt
echo -e "BRBO203\tB3" >> pops_35.txt
echo -e "BRBO205\tB3" >> pops_35.txt
echo -e "BRBO204\tB4" >> pops_35.txt

echo -e "RFBO101\tO" >> pops_35.txt
echo -e "RFBO102\tO" >> pops_35.txt
echo -e "RFBO103\tO" >> pops_35.txt
echo -e "RFBO104\tO" >> pops_35.txt
echo -e "RFBO105\tO" >> pops_35.txt
echo -e "RFBO106\tO" >> pops_35.txt



pops=$(python ~/sula/parsepops.py /data5/sulidae/final/to_flightless/abba/brbo/pops_35.txt)

python ~/genomics_general-master/freq.py -g /data5/sulidae/final/to_flightless/abba/sula.geno.gz \
$pops --popsFile /data5/sulidae/final/to_flightless/abba/brbo/pops_35.txt -f phased --target derived \
-o /data5/sulidae/final/to_flightless/abba/brbo/brbo_35.geno.tsv

Rscript ~/sula/ABBA_BRBO_35.r



B3  B5  B12
echo 'brbo 3 5 12' >> brbo.abbawholegenome.stats.txt

echo -e "BRBO203\tP1" > pops.txt
echo -e "BRBO205\tP2" >> pops.txt
echo -e "BRBO201\tP3" >> pops.txt
echo -e "BRBO202\tP3" >> pops.txt
echo -e "RFBO101\tP4" >> pops.txt
echo -e "RFBO102\tP4" >> pops.txt
echo -e "RFBO103\tP4" >> pops.txt
echo -e "RFBO104\tP4" >> pops.txt
echo -e "RFBO105\tP4" >> pops.txt
echo -e "RFBO106\tP4" >> pops.txt



~/sula/abba.sh -a brbo_3_5_12 -b /data5/sulidae/final/to_flightless/abba/brbo/brbo_3_5_12/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/abba/to_flightless/sula.geno.gz -f /data5/sulidae/final/to_flightless/abba/brbo/brbo_3_5_12/pops.txt -k y -n 2

B12 B4 B35
echo -e "BRBO201\tP1" > pops.txt
echo -e "BRBO202\tP1" >> pops.txt
echo -e "BRBO204\tP2" >> pops.txt
echo -e "BRBO203\tP3" >> pops.txt
echo -e "BRBO205\tP3" >> pops.txt
echo -e "RFBO101\tP4" >> pops.txt
echo -e "RFBO102\tP4" >> pops.txt
echo -e "RFBO103\tP4" >> pops.txt
echo -e "RFBO104\tP4" >> pops.txt
echo -e "RFBO105\tP4" >> pops.txt
echo -e "RFBO106\tP4" >> pops.txt



~/sula/abba.sh -a brbo_12_4_35 -b /data5/sulidae/final/to_flightless/abba/brbo/brbo_12_4_35/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/to_flightless/abba/sula.geno.gz -f /data5/sulidae/final/to_flightless/abba/brbo/brbo_12_4_35/pops.txt -k y -n 2


# bfbo

echo -e "BFBO501\tP1" >> pops.txt
echo -e "BFBO502\tP1" >> pops.txt
echo -e "BFBO503\tP1" >> pops.txt
echo -e "BFBO505\tP2" >> pops.txt
echo -e "BFBO506\tP2" >> pops.txt
echo -e "BFBO504\tP3" >> pops.txt

echo -e "RFBO101\tP4" >> pops.txt
echo -e "RFBO102\tP4" >> pops.txt
echo -e "RFBO103\tP4" >> pops.txt
echo -e "RFBO104\tP4" >> pops.txt
echo -e "RFBO105\tP4" >> pops.txt
echo -e "RFBO106\tP4" >> pops.txt


~/sula/abba.sh -a bfbo -b /data5/sulidae/final/to_flightless/abba/bfbo/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/to_flightless/abba/sula.geno.gz -f /data5/sulidae/final/to_flightless/abba/bfbo/pops.txt -k y -n 2
