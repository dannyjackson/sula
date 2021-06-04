Sulidae Master Script

~/whole_genome_bioinformatics/trim-and-QC.sh -i /data5/sulafilenames.txt -p /data5/ -f R1_001.fastq.gz -r R2_001.fastq.gz -a /data5/NexteraPE-PE.fa -t 12
#important settings: ILLUMINACLIP:$adapters:1:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:90

align-and-sort.sh -i /data5/sulidae/reference_lists/sulafilenames.txt -r ~/reference_datasets/flightless/ncbi-genomes-2020-08-14/GCA_002173475.1_Pharrisi_ref_V1/GCA_002173475.1_Pharrisi_ref_V1_genomic.fna -t 12 -p /data5/sulidae/my_datasets/trimming_step/fastqs/ -b bam_files_flightless -s sorted_bam_files_flightless

ref="~/reference_datasets/flightless/ncbi-genomes-2020-08-14/GCA_002173475.1_Pharrisi_ref_V1/GCA_002173475.1_Pharrisi_ref_V1_genomic.fna"
bamdir="/data5/sulidae/my_datasets/sorted_bam_files/"
ID="sula_flightless"
/usr/bin/bcftools mpileup -Ou -f "$ref" -a FORMAT/AD,DP,INFO/AD,SP "$bamdir"*_sorted_RGadded_dupmarked.bam | /usr/bin/bcftools call -mv -V indels > "$ID"_snps_multiallelic.vcf

awk '{gsub("RFBO_PAL_63", "RFBO101"); print}' sula_flightless_snps_multiallelic.vcf > sula_flightless.vcf

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

#Note that the vcf is significantly reduced by the maf 0.01 command. Consider leaving this off for any analyses that don't require it (I think only RAxML requires it?)


~/plink --vcf /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.02 --mind 0.1 --maf 0.01 --recode vcf-iid --out sula_flightless_cleaned










# redoing the phylogeny with all the samples
~/plink --vcf /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno --maf 0.01 --recode vcf-iid --out sula_flightless_cleaned_all

cp sula_flightless_cleaned.vcf sula_flightless_cleaned_zip.vcf
bgzip sula_flightless_cleaned_zip.vcf

bcftools index sula_flightless_cleaned_zip.vcf.gz



# ... the geno and maf tags reduce this a whole lot... see if leaving them off works with the invariant filter script ?

/data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf --

cp sula_flightless_cleaned.vcf sula_flightless_cleaned_zip.vcf
bgzip sula_flightless_cleaned_zip.vcf

bcftools index sula_flightless_cleaned_zip.vcf.gz




# generate population statistics

echo -e "MABO301" >> mabo.txt
echo -e "MABO302" >> mabo.txt
echo -e "MABO304" >> mabo.txt
echo -e "MABO305" >> mabo.txt
echo -e "MABO306" >> mabo.txt

echo -e "NABO402" >> nabo.txt
echo -e "NABO403" >> nabo.txt
echo -e "NABO404" >> nabo.txt
echo -e "NABO405" >> nabo.txt
echo -e "NABO406" >> nabo.txt


vcftools --vcf /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf --weir-fst-pop mabo.txt --weir-fst-pop nabo.txt --out mabo_vs_nabo


echo -e "MABO302" >> mabo_indopac.txt
echo -e "MABO306" >> mabo_indopac.txt

echo -e "MABO304" >> mabo_atlcar.txt
echo -e "MABO305" >> mabo_atlcar.txt


vcftools --vcf /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf --weir-fst-pop mabo_atlcar.txt --weir-fst-pop mabo_indopac.txt --out mabo_vs_mabo



echo -e "BFBO501" >> bfbo.txt
echo -e "BFBO502" >> bfbo.txt
echo -e "BFBO503" >> bfbo.txt
echo -e "BFBO504" >> bfbo.txt
echo -e "BFBO505" >> bfbo.txt


echo -e "PEBO601" >> pebo.txt
echo -e "PEBO603" >> pebo.txt
echo -e "PEBO604" >> pebo.txt
echo -e "PEBO605" >> pebo.txt
echo -e "PEBO606" >> pebo.txt


vcftools --vcf /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf --weir-fst-pop bfbo.txt --weir-fst-pop pebo.txt --out bfbo_vs_pebo


echo -e "BRBO201" >> brbo_p.txt
echo -e "BRBO202" >> brbo_p.txt


echo -e "BRBO203" >> brbo_a.txt
echo -e "BRBO205" >> brbo_a.txt



vcftools --vcf /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf --weir-fst-pop brbo_p.txt --weir-fst-pop brbo_a.txt --out brboa_brbop



# get just fixed sites from vcftools output

grep -E "(^|\s)1(\s|$)" bfbo_vs_pebo.weir.fst > bfbo_vs_pebo.fixed

grep -E "(^|\s)1(\s|$)" mabo_vs_nabo.weir.fst > mabo_vs_nabo.fixed

grep -E "(^|\s)1(\s|$)" mabo_vs_mabo.weir.fst > mabo_vs_mabo.fixed


# run RAxML


# all individuals included

~/vcf2phylip/vcf2phylip.py -i sula_flightless_cleaned_all.vcf
python ~/sula/filter_invariants_.py sula_flightless_cleaned_all.min4.phy



# pruned samples

~/vcf2phylip/vcf2phylip.py -i sula_flightless_flightless_mind2.vcf

mkdir pruned
cd pruned
python ~/sula/filter_invariants_all.py ../sula_flightless_flightless_mind2.min4.phy
mv variantsites.phy ../variantsites_mind2.phy
mv variantsites_kept.txt ../variantsites_mind2_kept.txt
cd ..
rm -r pruned


raxmlHPC -m ASC_GTRCAT --asc-corr felsenstein -f d -d -k -n sula_flightless_b -q /data5/sulidae/final/to_flightless/raxml/partition_file.txt -s /data5/sulidae/my_datasets/variantsites_mind2.phy -T 6 -p 12345 -N 10 ­-b 12345 -V

# I've run it through here... not sure what the notes below are doing. doublecheck the "trying to bootstrap" section



# trying to bootstrap
# added in -f d
raxmlHPC ­-f a -x 12345 -p 12345 -# 1000 -n sula_flightless -m ASC_GTRCAT --asc-corr felsenstein -s /data5/sulidae/my_datasets/variantsites_mind2.phy -q /data5/sulidae/final/to_flightless/raxml/partition_file.txt -T 6

raxmlHPC -f b -m PROTGAMMAILG -n sula.tre -t RAxML_bestTree.sula_flightless_b -z  RAxML_bootstrap.sula_flightless

  # ASC specifies that it needs to correct for ascertainment bias (SNPs/invariant sites) and the GTRGAMMA specifies that I want a General Time Reversible model with optimization of substitution rates under the GAMMA model of rate heterogeneity
  # -d starts ML optimization from a random starting tree, instead of the default randomized stepwise addition Maximum Parsimony starting tree,which the manual states "sometimes yields better - more diverse - starting trees for the analysis of broad phylogenomic alignments that have a very strong phylogenetic signal"

  #-f k addresses the issue of very long branch lengths in partitioned datasets with missing data by stealing branch lengths from those partitions that have data on both sides of the branch under consideration. If there are several other partitions that have data it computes a weighted average for the stolen branch length based on the site counts in each partition from which it stole a branch length. The nice property of this algorithm is that by changing the branch lengths in this way, the likelihood of the tree is not changed.

  # -J STRICT computes a strict consensus tree

  # -k specifies that bootstrapped trees should be printed with branch lengths

  # -n specifies the name of the output file

  # -q path to partition file. Specifies regions of alignment for which an individual model of nucleotide substitution should be estimated.


/data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf


~/vcf2phylip/vcf2phylip.py -i /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf

python ~/sula/filter_invariants.py sula_flightless_cleaned_all.min4.phy

raxmlHPC -m ASC_GTRCAT --asc-corr felsenstein -f d -d -k -n sula_flightless_b -q /data5/sulidae/my_datasets/partition_file.txt -s /data5/sulidae/my_datasets/variantsites.phy -T 6 -p 12345 -N 10 ­-b 12345 -V







# SVDQuartets

# all individuals
# Convert vcf to nexus file
git clone https://github.com/ODiogoSilva/ElConcatenero.git

python ElConcatenero/ElConcatenero.py -if phylip -of nexus -in variantsites_mind2.phy -o variantsites_mind2.nex


sed 's/datatype=mixed ()/datatype=dna/g' /data5/sulidae/my_datasets/variantsites_mind2.nex > test.nex
mv test.nex /data5/sulidae/my_datasets/variantsites_mind2.nex


~/programs/whole_genome_bioinformatics/paup4a166_ubuntu64

exe /data5/sulidae/my_datasets/variantsites_mind2.nex
SVDQuartets nquartets=500000 nreps=200 bootstrap=standard treeFile=sula.tre
contree all /strict=yes saveSupport=Both treefile=sulaconsensus_strict.tre;
contree / strict=no majrule=yes treefile=sulamajrule_con.tree;

SVDQuartets nquartets=500000 nreps=10000 bootstrap=standard treeFile=sula_10kreps.tre
contree all /strict=yes saveSupport=Both treefile=sulaconsensus_strict_10kreps.tre;
contree / strict=no majrule=yes treefile=sulamajrule_con_10kreps.tree;

# pruned individuals

# Convert vcf to nexus file
python ElConcatenero/ElConcatenero.py -if phylip -of nexus -in ../../variantsites_mind2.phy -o variantsites_mind2

sed -r 's/mixed [(][)]/DNA/g' variantsites_mind2.nex > new.nex
mv new.nex variantsites_mind2.nex


~/programs/whole_genome_bioinformatics/paup4a166_ubuntu64

exe variantsites_mind2.nex
SVDQuartets nquartets=500000 nreps=200 bootstrap=standard treeFile=sula.tre
contree all /strict=yes saveSupport=Both treefile=sulaconsensus_strict_mind2.tre;
contree / strict=no majrule=yes treefile=sulamajrule_con_mind2.tree;






# ABBA BABA

# Tests of introgression between species

# Convert vcf to geno
python ~/genomics_general-master/VCF_processing/parseVCF.py -i /data5/sulidae/my_datasets/sula_flightless_filtered.recode.vcf -o /data5/sulidae/final/to_flightless/abba/sula.geno.gz

~/sula/abba.sh -a [NAME OF RUN] -b [PATH TO OUTPUT] -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/abba/sula.geno.gz -f /data5/sulidae/final/abba/NAMEOFRUN/pops.txt -g P1 -h P2 -i P3 -j P4 -k y -l 15000 -m 100 -n 2 -o 0.9 -p ~/my_datasets/satsuma_summary.chained.out

# BF_GCA    BF_LOBOS      PE  RF

echo -e "BFBO504\tP1" > pops.txt
echo -e "BFBO502\tP2" >> pops.txt
echo -e "BFBO503\tP2" >> pops.txt
echo -e "PEBO601\tP3" >> pops.txt
echo -e "PEBO602\tP3" >> pops.txt
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

echo -e "BFBO504\tP1" > pops_admx.txt
echo -e "BFBO502\tP2" >> pops_admx.txt
echo -e "BFBO503\tP2" >> pops_admx.txt
echo -e "PEBO603\tP3a" >> pops_admx.txt
echo -e "PEBO604\tP3b" >> pops_admx.txt
echo -e "PEBO602\tP3b" >> pops_admx.txt
echo -e "PEBO606\tP3b" >> pops_admx.txt
echo -e "RFBO101\tP4" >> pops_admx.txt
echo -e "RFBO102\tP4" >> pops_admx.txt
echo -e "RFBO103\tP4" >> pops_admx.txt
echo -e "RFBO104\tP4" >> pops_admx.txt
echo -e "RFBO105\tP4" >> pops_admx.txt
echo -e "RFBO106\tP4" >> pops_admx.txt

~/sula/abba.sh -a BfBfPeRf_500kb -b /data5/sulidae/final/to_flightless/abba/BfBfPeRf/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/to_flightless/abba/sula.geno.gz -f /data5/sulidae/final/to_flightless/abba/BfBfPeRf/pops.txt -g /data5/sulidae/final/to_flightless/abba/BfBfPeRf/pops_admx.txt -k y -l 500000 -m 100 -n 2

~/sula/abba.sh -a BfBfPeRf_10kb -b /data5/sulidae/final/abba/BfBfPeRf/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/abba/sula.geno.gz -f /data5/sulidae/final/abba/BfBfPeRf/pops.txt -g /data5/sulidae/final/abba/BfBfPeRf/pops_admx.txt -k y -l 10000 -m 25 -n 2

~/sula/abba.sh -a BfBfPeRf_5kb -b /data5/sulidae/final/abba/BfBfPeRf/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/abba/sula.geno.gz -f /data5/sulidae/final/abba/BfBfPeRf/pops.txt -g /data5/sulidae/final/abba/BfBfPeRf/pops_admx.txt -k y -l 5000 -m 3 -n 2

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

echo -e "MABO304\tP1" > pops_admx.txt
echo -e "MABO305\tP1" >> pops_admx.txt
echo -e "MABO302\tP2" >> pops_admx.txt
echo -e "MABO306\tP2" >> pops_admx.txt
echo -e "NABO402\tP3" >> pops_admx.txt
echo -e "NABO403\tP3a" >> pops_admx.txt
echo -e "NABO405\tP3a" >> pops_admx.txt
echo -e "NABO404\tP3b" >> pops_admx.txt
echo -e "NABO406\tP3b" >> pops_admx.txt
echo -e "RFBO101\tP4" >> pops_admx.txt
echo -e "RFBO102\tP4" >> pops_admx.txt
echo -e "RFBO103\tP4" >> pops_admx.txt
echo -e "RFBO104\tP4" >> pops_admx.txt
echo -e "RFBO105\tP4" >> pops_admx.txt
echo -e "RFBO106\tP4" >> pops_admx.txt

~/sula/abba.sh -a MaMaNaRf_500kb -b /data5/sulidae/final/to_flightless/abba/MaMaNaRf/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/to_flightless/abba/sula.geno.gz -f /data5/sulidae/final/to_flightless/abba/MaMaNaRf/pops.txt -g /data5/sulidae/final/to_flightless/abba/MaMaNaRf/pops_admx.txt -k y -l 500000 -m 100 -n 2

~/sula/abba.sh -a MaMaNaRf_10kb -b /data5/sulidae/final/abba/MaMaNaRf/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/abba/sula.geno.gz -f /data5/sulidae/final/abba/MaMaNaRf/pops.txt -g /data5/sulidae/final/abba/MaMaNaRf/pops_admx.txt -k y -l 10000 -m 25 -n 2

~/sula/abba.sh -a MaMaNaRf_5kb -b /data5/sulidae/final/abba/MaMaNaRf/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/abba/sula.geno.gz -f /data5/sulidae/final/abba/MaMaNaRf/pops.txt -g /data5/sulidae/final/abba/MaMaNaRf/pops_admx.txt -k y -l 5000 -m 3 -n 2

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



~/sula/abba.sh -a PeBfBrRf -b /data5/sulidae/final/to_flightless/abba/PeBfBrRf/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/to_flightless/abba/sula.geno.gz -f /data5/sulidae/final/to_flightless/abba/PeBfBrRf/pops.txt -k y -l 500000 -m 100 -n 2

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


echo -e "BFBO504\tP1" > pops_admx.txt
echo -e "BFBO502\tP2" >> pops_admx.txt
echo -e "BFBO503\tP2" >> pops_admx.txt
echo -e "NABO402\tP3a" >> pops_admx.txt
echo -e "NABO403\tP3a" >> pops_admx.txt
echo -e "NABO404\tP3b" >> pops_admx.txt
echo -e "NABO405\tP3b" >> pops_admx.txt
echo -e "NABO406\tP3b" >> pops_admx.txt
echo -e "RFBO101\tP4" >> pops_admx.txt
echo -e "RFBO102\tP4" >> pops_admx.txt
echo -e "RFBO103\tP4" >> pops_admx.txt
echo -e "RFBO104\tP4" >> pops_admx.txt
echo -e "RFBO105\tP4" >> pops_admx.txt
echo -e "RFBO106\tP4" >> pops_admx.txt

~/sula/abba.sh -a BfBfNaRf_500kb -b /data5/sulidae/final/to_flightless/abba/BfBfNaRf/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/to_flightless/abba/sula.geno.gz -f /data5/sulidae/final/to_flightless/abba/BfBfNaRf/pops.txt -g /data5/sulidae/final/to_flightless/abba/BfBfNaRf/pops_admx.txt -k y -l 500000 -m 100 -n 2

~/sula/abba.sh -a BfBfNaRf_10kb -b /data5/sulidae/final/abba/BfBfNaRf/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/abba/sula.geno.gz -f /data5/sulidae/final/abba/BfBfNaRf/pops.txt -g /data5/sulidae/final/abba/BfBfNaRf/pops_admx.txt -k y -l 500000 -m 25 -n 2

~/sula/abba.sh -a BfBfNaRf_5kb -b /data5/sulidae/final/abba/BfBfNaRf/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/abba/sula.geno.gz -f /data5/sulidae/final/abba/BfBfNaRf/pops.txt -g /data5/sulidae/final/abba/BfBfNaRf/pops_admx.txt -k y -l 5000 -m 3 -n 2



# PE  PE  BF  BF  BR  RF

echo -e "PEBO601\tP1" > pops.txt

echo -e "PEBO603\tP2" >> pops.txt
echo -e "PEBO605\tP2" >> pops.txt

echo -e "BFBO502\tP3" >> pops.txt
echo -e "BFBO503\tP3" >> pops.txt

echo -e "BFBO504\tP4" >> pops.txt

echo -e "RFBO101\tP5" >> pops.txt
echo -e "RFBO102\tP5" >> pops.txt
echo -e "RFBO103\tP5" >> pops.txt
echo -e "RFBO104\tP5" >> pops.txt
echo -e "RFBO105\tP5" >> pops.txt
echo -e "RFBO106\tP5" >> pops.txt

~/sula/abba.sh -a BfBfNaRf_500kb -b /data5/sulidae/final/to_flightless/abba/BfBfNaRf/ -c ~/sula/ -d ~/genomics_general-master/ -e /data5/sulidae/final/to_flightless/abba/sula.geno.gz -f /data5/sulidae/final/to_flightless/abba/BfBfNaRf/pops.txt -g /data5/sulidae/final/to_flightless/abba/BfBfNaRf/pops_admx.txt -k y -l 500000 -m 100 -n 2



python ~/genomics_general-master/VCF_processing/parseVCF.py -i /data5/sulidae/final/abba/sula.geno.gz -o /data5/sulidae/final/abba/5_PePeBfBfRf/PePeBfBfRf.geno.gz

pops=$(python ~/sula/parsepops.py /data5/sulidae/final/abba/5_PePeBfBfRf/pops.txt)

python ~/genomics_general-master/freq.py -g /data5/sulidae/final/abba/5_PePeBfBfRf/PePeBfBfRf.geno.gz \
$pops --popsFile /data5/sulidae/final/abba/5_PePeBfBfRf/pops.txt -f phased --target derived \
-o /data5/sulidae/final/abba/5_PePeBfBfRf/PePeBfBfRf.geno.tsv



# BF  BF  PE  PE  BR  RF


echo -e "BFBO504\tP1" > pops.txt

echo -e "BFBO502\tP2" >> pops.txt
echo -e "BFBO503\tP2" >> pops.txt

echo -e "PEBO603\tP3" >> pops.txt
echo -e "PEBO605\tP3" >> pops.txt

echo -e "PEBO601\tP4" >> pops.txt

echo -e "RFBO101\tP5" >> pops.txt
echo -e "RFBO102\tP5" >> pops.txt
echo -e "RFBO103\tP5" >> pops.txt
echo -e "RFBO104\tP5" >> pops.txt
echo -e "RFBO105\tP5" >> pops.txt
echo -e "RFBO106\tP5" >> pops.txt

screen -r 68776

python ~/genomics_general-master/VCF_processing/parseVCF.py -i /data5/sulidae/final/abba/sula.geno.gz -o /data5/sulidae/final/abba/5_BfBfPePeRf/BfBfPePeRf.geno.gz

pops=$(python ~/sula/parsepops.py /data5/sulidae/final/abba/5_BfBfPePeRf/pops.txt)

python ~/genomics_general-master/freq.py -g /data5/sulidae/final/abba/5_BfBfPePeRf/BfBfPePeRf.geno.gz \
$pops --popsFile /data5/sulidae/final/abba/5_BfBfPePeRf/pops.txt -f phased --target derived \
-o /data5/sulidae/final/abba/5_BfBfPePeRf/BfBfPePeRf.geno.tsv







# Subset ABBA BABA sliding windows by fd value, and pull those scaffolds from the gff

python ~/sula/subset_ABBABABAwindows_output.py MaMaNaRf_500kb 0.9
python ~/sula/subset_ABBABABAwindows_output.py MaMaNaRf_10kb 0.9
python ~/sula/subset_ABBABABAwindows_output.py MaMaNaRf_5kb 0.9

zcat MaMaNaRf_10kb_slidingwindows.csv.gz | gawk --posix -f "," '$12 >= 0.9'

# Subset gff file to genes found in sliding window output.

file = "MaMaNaRf_10kb_slidingwindows.subsetfd.txt"
gff_file = "/home/daja5529/reference_datasets/gigadb/Phalacrocorax_carbo.gff"
project_name = "test"

python ~/sula/subset_gff_by_slidingwindow.py MaMaNaRf_500kb /home/daja5529/reference_datasets/gigadb/Phalacrocorax_carbo.gff

python ~/sula/subset_gff_by_slidingwindow.py MaMaNaRf_10kb /home/daja5529/reference_datasets/gigadb/Phalacrocorax_carbo.gff

python ~/sula/subset_gff_by_slidingwindow.py MaMaNaRf_5kb /home/daja5529/reference_datasets/gigadb/Phalacrocorax_carbo.gff



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

# I need to rethink this one... include MABO301

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

#hmmm
I still need to compute:
M4  M5  M162
M1,6  M2  M45
Maybe best to do by hand

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


### resume ###
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
