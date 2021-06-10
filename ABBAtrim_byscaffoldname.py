import sys
import numpy as np
import pandas as pd

mask = pd.read_csv(sys.argv[1], sep='\t', names="a")

df = pd.read_csv(sys.argv[2], sep='\t')

filtered = df[df['scaffold'].isin(mask['a'])]

filtered.to_csv('trimmed.geno.tsv',sep='\t', na_rep='NA', index=False)
