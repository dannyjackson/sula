# python script for trimming abba sliding window file

import sys

file = "bfpe_slidingwindows.subsetfd.txt"

dataset = []

with open(file) as f:
    for column in f:
        dataset.append(column.split(","))

sequence = list(range(0,len(dataset[0])))

keep = []
for item in sequence:
  if dataset[item][9] > 0.9:
    keep.append(dataset[item])

s = ","
s.join(str(keep))

with open("trimmed.abba.txt","w") as file:
    file.writelines( "%s" % "".join(str(item)) for item in keep )
