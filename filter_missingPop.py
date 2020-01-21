# python

# 1. Iterate through each column of the 012 file
# 2. Evaluate each site to identify if it is missing for all BFBO individuals
# 3. Evaluate each site to identify if it is missing for all PEBO individuals
# 4. If either 2 or 3 is TRUE, remove the site from both the 012 file and the locus_names file.
import sys, subprocess

p1_indv = sys.argv[1]
p2_indv = sys.argv[2]

print(sys.argv[1])
print(sys.argv[2])
code = []
with open("outflank_SNPmat.txt") as f:
    for column in f:
        code.append(column.split())

locinames = []
with open("outflank_loci_names.txt") as loci:
    for thing in loci:
        locinames.append(thing.split())

sequence = list(range(0,len(code[0])))

p1_missing = set()
p2_missing = set()
pop1_generator = list(range(0,int(p1_indv)))
pop2_generator = list(range(int(p1_indv),int(p2_indv)))
removesites = []

for i in sequence:
    allele_column = [l[i] for l in code]
    for site in pop1_generator:
        if int(allele_column[int(site)]) is 9:
            p1_missing.add(site)
        if len(p1_missing) is 6:
            removesites.append(i)
    for site in pop2_generator:
        if int(allele_column[int(site)]) is 9:
            p2_missing.add(site)
        if len(p2_missing) is 6:
            removesites.append(i)
    p1_missing = set()
    p2_missing = set()

numbersites = range(len(removesites))
print(len(numbersites))

for rm in numbersites:
    sequence.remove(removesites[rm])


x = range(len(sequence))
y = range(len(code))

trimmed_matrix = []
for w in y:
    for z in x:
        trimmed_matrix.append(code[w][sequence[z]])
        trimmed_matrix.append(" ")
    trimmed_matrix.append("\n")

trimmed_loci = []
for z in x:
    trimmed_loci.append(locinames[sequence[z]-1])
    trimmed_loci.append("\n")


with open("trimmed.012.txt","w") as file:
    file.writelines( "%s" % "".join(item) for item in trimmed_matrix )

with open("trimmed.locinames.txt","w") as file:
  file.writelines( "%s" % "".join(item) for item in trimmed_loci )
