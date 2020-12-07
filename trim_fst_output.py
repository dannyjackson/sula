# A script to drop fst data from scaffolds below 5k:

fai_file = "/home/daja5529/reference_datasets/flightless/ncbi-genomes-2020-08-14/GCA_002173475.1_Pharrisi_ref_V1/GCA_002173475.1_Pharrisi_ref_V1_genomic.fna.fai"
fst_file = "fst_10k.windowed.weir.fst"

fai_subset = []
final = []


with open(fai_file) as f:
  for column in f:
    if float(column.split("\t")[1]) > 100000 :
      fai_subset.append((column.split("\t")[0].split("NEVG0")[1]))

with open (fst_file) as f:
  for item in f:
    if item.split("\t")[0] in fai_subset:
      final.append(item)

with open("trimmed_25.fst.txt","w") as file:
    file.writelines("%s" %item for item in final )



# A script to subset fst file to just fixed sites:

fst_file = "fst_snps.weir.fst"

no_nan = []
final = []
exclude = ['-nan', 'WEIR_AND_COCKERHAM_FST']
with open(fst_file) as f:
  for row in f:
    if row.split("\t")[2].split("\n")[0] not in exclude:
        no_nan.append(row)

for item in no_nan:
    if float(item.split("\t")[2].split("\n")[0]) > 0.99 :
      final.append(item)


with open("fixed.fst.txt","w") as file:
    file.writelines("%s" %item for item in final )




# Identify genes with fixed Fst between species

import pandas as pd

fai_file = pd.read_csv("/data5/sulidae/reference_datasets/Phalacrocorax_noC.fa.fai",sep="\t")

fixed_file = pd.read_csv("fixed.fst.txt",sep="\t", header=None)

fixed_file = fixed_file.drop([2],axis=1)

if fixed_file[]





# Looking for rows in the gff with 1. matching column 1 and 2. $2 of fixed is $4 < $5 in gff


gff_file = open("/home/daja5529/reference_datasets/gigadb/Phalacrocorax_carbo.gff", "r")
chrom = []
subset = []
all = []

for item in gff_file.readlines():
    if item[0] is 'S':
        subset.append(item[8:])
    if item[0] is 's':
        subset.append(item[8:])

for item in subset:
    if item.split("\t")[0].split("caffold")[1] in fixed_dict:
        chrom.append(item)





item = gff_file.readlines()

for row in item:
    if row[0] is 'C':
        next(row)
    if row.split("\t")[0].split("caffold")[1] in fixed_dict:
        subset.append(row)
    row = item.readlines()

gff_file.close()
