# ANGSD Notes
# Trim for quality and adapters (Trimmomatic)
# Align, sort, markdups
# Clip overlapping read pairs

clipping () {
echo "clipping" "$@" >> /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/clippingstats.txt 

/home/daja5529/programs/bamutil_local/bam clipOverlap --in /data5/sulidae/my_datasets/trimmedfiles/sorted_bam_files/"$@"_sorted_RGadded_dupmarked.bam --out /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/"$@".sorted.marked.clipped.bam --stats --params


echo "done " $@ >>/data5/sulidae/my_datasets/trimmedfiles/clipoverlap/clippingstats.txt 
}


export -f clipping 

parallel -j 12 clipping :::: /data5/sulidae/reference_lists/sulafilenames.txt  


/home/daja5529/programs/bamutil_local/bam clipOverlap --in /data5/sulidae/my_datasets/trimmedfiles/sorted_bam_files/BFBO501_sorted_RGadded_dupmarked.bam --out /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO501_sorted.marked.clipped.bam --stats --params --poolSize 10000000
BFBO501_sorted_RGadded_dupmarked.bam --poolSize 1000000


while read -r bird
do 
/home/daja5529/programs/bamutil_local/bam clipOverlap --in /data5/sulidae/my_datasets/trimmedfiles/sorted_bam_files/"$bird"_sorted_RGadded_dupmarked.bam --out /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/"$bird"_sorted.marked.clipped.bam --stats --params --poolSize 10000000

done < /data5/sulidae/reference_lists/sulafilenames.txt  


# Indel realignment



while read -r bird;
do 
  samtools index "$bird".sorted.marked.clipped.bam
done < /data5/sulidae/reference_lists/sulafilenames.txt 



cd /data5/sulidae/my_datasets/trimmedfiles/indelmaps/

java -jar picard.jar Create


while read -r bird;
do 
~/programs/jdk1.8.0_411/bin/java -jar ~/programs/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R genome.fa \
    -I /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/"$bird".sorted.marked.clipped.bam \
    -o /data5/sulidae/my_datasets/trimmedfiles/indelmaps/"$bird".intervals

done < /data5/sulidae/reference_lists/sulafilenames.txt 
SequenceDictionary -R genome.fa -O genome.dict

~/programs/jdk1.8.0_411/bin/java -jar ~/programs/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R genome.fa \
    -I /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO101.sorted.marked.clipped.bam \
    -o /data5/sulidae/my_datasets/trimmedfiles/indelmaps/RFBO101.intervals

# waiting for 34 files
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/*bam | wc -l

while read -r bird;
do 
  ~/programs/jdk1.8.0_411/bin/java -jar ~/programs/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  -R genome.fa \
  --consensusDeterminationModel USE_READS \
  -I /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/"$bird".sorted.marked.clipped.bam \
  --targetIntervals /data5/sulidae/my_datasets/trimmedfiles/indelmaps/"$bird".intervals \
  -o /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/"$bird".realigned.bam

done < /data5/sulidae/reference_lists/sulafilenames.txt 



# comupute depth stats
while read -r bird;
do 
~/programs/angsd/angsd -I /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/"$bird".realigned.bam -doDepth 1 -nThreads 8 -out "$bird".depth -doCounts 1

done < /data5/sulidae/reference_lists/sulafilenames.txt 

ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/*bam > bamlist.txt
~/programs/angsd/angsd -bam bamlist.txt -doDepth 1 -nThreads 8 -out all.depth -doCounts 1

# compute quality scores

while read -r bird;
do 
~/programs/angsd/angsd -I /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/"$bird".realigned.bam -doQsDist 1 -nThreads 8 -out "$bird".quality -doCounts 1

done < /data5/sulidae/reference_lists/sulafilenames.txt 

# compute mapping quality
# python script for summarizing depth stats from ANGSD
python
import numpy as np
import pandas as pd 

df = pd.read_csv('all.depth.depthSample', sep='\t', header=None)
cols = df.shape[1]
multiplier = list(range(1, cols+1))

df_adj = df.mul(multiplier, axis = 1)

counts = df.sum(axis = 1)
totaldepth = df_adj.sum(axis = 1)
avgdepth = totaldepth / counts

# 6.68432

df_names = pd.read_csv('/xdisk/mcnew/dannyjackson/finches/reference_lists/sample_species_treatment.txt', sep='\t')

df_final = pd.concat([df_names, avgdepth], axis=1)

df_final = df_final.rename({0: 'avgdepth'}, axis=1)

df_final.groupby(['species', 'treatment']).mean()

df_final.sort_values(by=['avgdepth']).round(decimals=2)

df_final[df_final["species"] == 'CRA'].sort_values(by=['avgdepth']).round(decimals=2)

df_final[df_final["species"] == 'FOR'].sort_values(by=['avgdepth']).round(decimals=2)

df_final[df_final["species"] == 'PAR'].sort_values(by=['avgdepth']).round(decimals=2)

print(df_final.groupby(['species', 'treatment']).size())

df_final.groupby(['species', 'treatment']).mean('MEAN_DEPTH')



# MA_ATLCAR MA_INDOPA     NA  RF
cd /data5/sulidae/angsd/MANA


ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/MABO304*bam > angsd.MANA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/MABO305*bam >> angsd.MANA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/MABO302*bam >> angsd.MANA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/MABO306*bam >> angsd.MANA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/NABO402*bam >> angsd.MANA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/NABO403*bam >> angsd.MANA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/NABO404*bam >> angsd.MANA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/NABO405*bam >> angsd.MANA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/NABO406*bam >> angsd.MANA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO101*bam >> angsd.MANA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO102*bam >> angsd.MANA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO103*bam >> angsd.MANA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO104*bam >> angsd.MANA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO105*bam >> angsd.MANA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO106*bam >> angsd.MANA.pops

echo "2" > angsd.MANA.abba
echo "2" >> angsd.MANA.abba
echo "5" >> angsd.MANA.abba
echo "6" >> angsd.MANA.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.MANA.pops -sizeFile angsd.MANA.abba -doCounts 1 -out bam -anc ~/reference_datasets/flightless/ncbi-genomes-2020-08-14/GCA_002173475.1_Pharrisi_ref_V1/GCA_002173475.1_Pharrisi_ref_V1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 500000 -enhance 1 

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.MANA.abba 




# Test 2a	Brown (3,5)	Brown (1,2)	Blue-footed (All)	Red-footed (all)
mkdir testa testb testc testd teste testf
cd testa 

ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO203*bam > angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO205*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO201*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO202*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO501*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO502*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO503*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO504*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO505*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO506*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO101*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO102*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO103*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO104*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO105*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO106*bam >> angsd.BRBF.pops

echo "2" > angsd.BRBF.abba
echo "2" >> angsd.BRBF.abba
echo "6" >> angsd.BRBF.abba
echo "6" >> angsd.BRBF.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BRBF.pops -sizeFile angsd.BRBF.abba -doCounts 1 -out bam -anc ~/reference_datasets/flightless/ncbi-genomes-2020-08-14/GCA_002173475.1_Pharrisi_ref_V1/GCA_002173475.1_Pharrisi_ref_V1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 500000 -enhance 1 

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam.Angsd" out="result" sizeFile=angsd.BRBF.abba 


# Test 2b Blue-footed (1,2,3,5)	Blue-footed (4)	Brown (1)	Red-footed (all)
cd ../testb


ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO501*bam > angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO502*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO503*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO505*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO504*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO201*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO101*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO102*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO103*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO104*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO105*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO106*bam >> angsd.BRBF.pops

echo "4" > angsd.BRBF.abba
echo "1" >> angsd.BRBF.abba
echo "1" >> angsd.BRBF.abba
echo "6" >> angsd.BRBF.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BRBF.pops -sizeFile angsd.BRBF.abba -doCounts 1 -out bam -anc ~/reference_datasets/flightless/ncbi-genomes-2020-08-14/GCA_002173475.1_Pharrisi_ref_V1/GCA_002173475.1_Pharrisi_ref_V1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 500000 -enhance 1 

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam.Angsd" out="result" sizeFile=angsd.BRBF.abba 


# Test 2c	Brown (3,5)	Brown (1,2)	Peruvian (All)	Red-footed (all)
cd ../testc


ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO203*bam > angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO205*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO201*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO202*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/PEBO601*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/PEBO602*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/PEBO603*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/PEBO604*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/PEBO605*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/PEBO606*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO101*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO102*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO103*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO104*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO105*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO106*bam >> angsd.BRBF.pops

echo "2" > angsd.BRBF.abba
echo "2" >> angsd.BRBF.abba
echo "6" >> angsd.BRBF.abba
echo "6" >> angsd.BRBF.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BRBF.pops -sizeFile angsd.BRBF.abba -doCounts 1 -out bam -anc ~/reference_datasets/flightless/ncbi-genomes-2020-08-14/GCA_002173475.1_Pharrisi_ref_V1/GCA_002173475.1_Pharrisi_ref_V1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 500000 -enhance 1 


Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.BRBF.abba 


# Test 2d	Brown (3,5)	Brown (1,2)	Masked (All)	Red-footed (all)
cd ../testd


ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO203*bam > angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO205*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO201*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO202*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/MABO301*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/MABO302*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/MABO304*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/MABO305*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/MABO306*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO101*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO102*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO103*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO104*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO105*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO106*bam >> angsd.BRBF.pops

echo "2" > angsd.BRBF.abba
echo "2" >> angsd.BRBF.abba
echo "5" >> angsd.BRBF.abba
echo "6" >> angsd.BRBF.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BRBF.pops -sizeFile angsd.BRBF.abba -doCounts 1 -out bam -anc ~/reference_datasets/flightless/ncbi-genomes-2020-08-14/GCA_002173475.1_Pharrisi_ref_V1/GCA_002173475.1_Pharrisi_ref_V1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 500000 -enhance 1 

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.BRBF.abba 


# Test 2e	Brown (3,5)	Brown (1,2)	Nazca (All)	Red-footed (all)
cd ../teste


ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO203*bam > angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO205*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO201*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO202*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/NABO401*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/NABO402*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/NABO403*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/NABO404*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/NABO405*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/NABO406*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO101*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO102*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO103*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO104*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO105*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO106*bam >> angsd.BRBF.pops

echo "2" > angsd.BRBF.abba
echo "2" >> angsd.BRBF.abba
echo "6" >> angsd.BRBF.abba
echo "6" >> angsd.BRBF.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BRBF.pops -sizeFile angsd.BRBF.abba -doCounts 1 -out bam -anc ~/reference_datasets/flightless/ncbi-genomes-2020-08-14/GCA_002173475.1_Pharrisi_ref_V1/GCA_002173475.1_Pharrisi_ref_V1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 500000 -enhance 1 


Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam.Angsd" out="result" sizeFile=angsd.BRBF.abba 



# Test 2f	Peruvian (all)	Blue-footed (all)	Brown (1,2)	Red-footed (all)
cd ../testf


ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/PEBO601*bam > angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/PEBO602*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/PEBO603*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/PEBO604*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/PEBO605*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/PEBO606*bam >> angsd.BRBF.pops

ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO501*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO502*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO503*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO504*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO505*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO506*bam >> angsd.BRBF.pops

ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO201*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO202*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO203*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO204*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO205*bam >> angsd.BRBF.pops

ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO101*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO102*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO103*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO104*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO105*bam >> angsd.BRBF.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO106*bam >> angsd.BRBF.pops

echo "6" > angsd.BRBF.abba
echo "6" >> angsd.BRBF.abba
echo "5" >> angsd.BRBF.abba
echo "6" >> angsd.BRBF.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BRBF.pops -sizeFile angsd.BRBF.abba -doCounts 1 -out bam -anc ~/reference_datasets/flightless/ncbi-genomes-2020-08-14/GCA_002173475.1_Pharrisi_ref_V1/GCA_002173475.1_Pharrisi_ref_V1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 500000 -enhance 1 


Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam.Angsd" out="result" sizeFile=angsd.BRBF.abba 



# Test 3	Blue-footed (4)	Blue-footed (2,3)	Peruvian (all)	Red-footed (all)
cd /data5/sulidae/angsd/BFPE

ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO504*bam > angsd.BFPE.pops

ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO502*bam >> angsd.BFPE.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BFBO503*bam >> angsd.BFPE.pops

ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/PEBO601*bam >> angsd.BFPE.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/PEBO603*bam >> angsd.BFPE.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/PEBO604*bam >> angsd.BFPE.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/PEBO605*bam >> angsd.BFPE.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/PEBO606*bam >> angsd.BFPE.pops



ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO101*bam >> angsd.BFPE.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO102*bam >> angsd.BFPE.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO103*bam >> angsd.BFPE.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO104*bam >> angsd.BFPE.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO105*bam >> angsd.BFPE.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO106*bam >> angsd.BFPE.pops

echo "1" > angsd.BFPE.abba
echo "2" >> angsd.BFPE.abba
echo "5" >> angsd.BFPE.abba
echo "6" >> angsd.BFPE.abba

echo "BFBO_CA" > angsd.BFPE.abba.name
echo "BFBO_peru" >> angsd.BFPE.abba.name
echo "PEBO" >> angsd.BFPE.abba.name
echo "RFBO" >> angsd.BFPE.abba.name

~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BFPE.pops -sizeFile angsd.BFPE.abba -doCounts 1 -out bam.500k -anc ~/reference_datasets/flightless/ncbi-genomes-2020-08-14/GCA_002173475.1_Pharrisi_ref_V1/GCA_002173475.1_Pharrisi_ref_V1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 500000 -enhance 1 

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam.500k" out="result" sizeFile=angsd.BFPE.abba



# Test 4a	Peruvian (all)	Blue-footed (all)	Nazca (all)	Red-footed (all)
cd /data5/sulidae/angsd/BFNA/testa

ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/PEBO601*bam > angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/PEBO603*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/PEBO604*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/PEBO605*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/PEBO606*bam >> angsd.BFNA.pops

ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/BFBO501*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/BFBO502*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/BFBO503*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/BFBO504*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/BFBO505*bam >> angsd.BFNA.pops

ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/NABO402*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/NABO403*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/NABO404*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/NABO405*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/NABO406*bam >> angsd.BFNA.pops



ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/RFBO101*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/RFBO102*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/RFBO103*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/RFBO104*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/RFBO105*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/RFBO106*bam >> angsd.BFNA.pops

echo "5" > angsd.BFNA.abba
echo "5" >> angsd.BFNA.abba
echo "5" >> angsd.BFNA.abba
echo "6" >> angsd.BFNA.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BFNA.pops -sizeFile angsd.BFNA.abba -doCounts 1 -out bam -anc ~/reference_datasets/flightless/ncbi-genomes-2020-08-14/GCA_002173475.1_Pharrisi_ref_V1/GCA_002173475.1_Pharrisi_ref_V1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 500000 -enhance 1 


Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam.Angsd" out="result" sizeFile=angsd.BFNA.abba 



# Test 4b	Blue-footed (4)	Blue-footed (5)	Nazca (all)	Red-footed (all)
cd /data5/sulidae/angsd/BFNA/testb

ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/BFBO504*bam > angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/BFBO505*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/NABO402*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/NABO403*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/NABO404*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/NABO405*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/NABO406*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/RFBO101*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/RFBO102*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/RFBO103*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/RFBO104*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/RFBO105*bam >> angsd.BFNA.pops
ls /data5/sulidae/my_datasets/trimmedfiles/indelrealignment/RFBO106*bam >> angsd.BFNA.pops

echo "1" > angsd.BFNA.abba
echo "1" >> angsd.BFNA.abba
echo "5" >> angsd.BFNA.abba
echo "6" >> angsd.BFNA.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BFNA.pops -sizeFile angsd.BFNA.abba -doCounts 1 -out bam -anc ~/reference_datasets/flightless/ncbi-genomes-2020-08-14/GCA_002173475.1_Pharrisi_ref_V1/GCA_002173475.1_Pharrisi_ref_V1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 500000 -enhance 1 


Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam.Angsd" out="result" sizeFile=angsd.BFNA.abba 


# Test 5	Masked (1, 2)	Masked (4,5)	Brown (3,5)	Red-footed (all)
cd /data5/sulidae/angsd/MABR

ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/MABO301*bam > angsd.MABR.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/MABO302*bam >> angsd.MABR.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/MABO304*bam >> angsd.MABR.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/MABO305*bam >> angsd.MABR.pops

ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO203*bam >> angsd.MABR.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/BRBO205*bam >> angsd.MABR.pops

ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO101*bam >> angsd.MABR.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO102*bam >> angsd.MABR.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO103*bam >> angsd.MABR.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO104*bam >> angsd.MABR.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO105*bam >> angsd.MABR.pops
ls /data5/sulidae/my_datasets/trimmedfiles/clipoverlap/RFBO106*bam >> angsd.MABR.pops

echo "2" > angsd.MABR.abba
echo "2" >> angsd.MABR.abba
echo "2" >> angsd.MABR.abba
echo "6" >> angsd.MABR.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.MABR.pops -sizeFile angsd.MABR.abba -doCounts 1 -out bam -anc ~/reference_datasets/flightless/ncbi-genomes-2020-08-14/GCA_002173475.1_Pharrisi_ref_V1/GCA_002173475.1_Pharrisi_ref_V1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 500000 -enhance 1 

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam.Angsd" out="result" sizeFile=angsd.MABR.abba 
