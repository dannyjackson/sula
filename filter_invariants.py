# Python script for removing invariant sites from a phylip file
import sys, subprocess

phylip = sys.argv[1]

code = []
samples = []
with open(phylip) as f:
  next(f)
  for columns in f:
    code.append(columns[10:len(columns)])
    samples.append(columns[0:9])

sequence = list(range(0,len(code[0])))

variant_sites = set()
unique_alleles = set('N')
for i in sequence:
  allele_column = [l[i] for l in code]
  unique_alleles = set('N')
  for site in allele_column:
    if site not in unique_alleles:
      unique_alleles.add(site)
    if len(unique_alleles) > 2:
      variant_sites.add(i)

x = range(0,len(variant_sites))
y = range(0,len(code))

temp = []
for v in variant_sites:
  guess = [s[v] for s in code]
  temp.append(guess)

x = range(0,len(temp))
y = range(0,len(temp[0]))

header = []
header.append(str(len(samples)))
header.append(" ")
header.append(str(len(temp)))

trimmed_matrix = []
trimmed_matrix.append(header)
trimmed_matrix.append("\n")
for w in y:
  trimmed_matrix.append(samples[w])
  trimmed_matrix.append("\t")
  for z in x:
    trimmed_matrix.append(temp[z][w])
  trimmed_matrix.append("\n")

with open("variantsites_kept.txt","w") as file:
    file.writelines('%s\n' % item for item in variant_sites)

with open("variantsites.phy","w") as file:
  file.writelines( "%s" % "".join(item) for item in trimmed_matrix )
