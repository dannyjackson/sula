import sys

file = sys.argv[1] + "_slidingwindows.subsetfd.txt"
project_name = sys.argv[1] + "_satsumachained.abbasubset.txt"
reference = sys.argv[2]

subsetSatsuma = open(project_name, "w+")



with open(file) as f:
    f_lines = f.readlines()

with open(reference) as r:
    r_lines = r.readlines()

f_list = [i.split(",")[0].rstrip() for i in f_lines]
r_list = [j.split()[3].rstrip() for j in r_lines]

for i in range(len(r_list)):
    if r_list[i] in f_list:
        subsetSatsuma.write(r_lines[i])

subsetSatsuma.close()
