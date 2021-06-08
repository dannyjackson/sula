import numpy as np
import pandas as pd
import sys

file = "/home/daja5529/reference_datasets/flightless/ncbi-genomes-2020-08-14/GCA_002173475.1_Pharrisi_ref_V1/GCA_002173475.1_Pharrisi_ref_V1_genomic_sortedbyscaffoldsize.fna.fai"

with open(file) as fin:
     rows = ( line.split('\t') for line in fin )
     d = { row[0]:row[1] for row in rows }

inv_d = {v: k for k, v in d.items()}



df = pd.read_csv(sys.argv[1], sep='\t')

df['CHROM'] = df['CHROM'].map(inv_d)

df.to_csv(sys.argv[1],index=False)
