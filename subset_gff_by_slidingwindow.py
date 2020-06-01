import sys

file = sys.argv[1] + "_slidingwindows.subsetfd.txt"
gff_file = sys.argv[2]
project_name = sys.argv[1]

outputList = open(project_name, "a")

relevant_scaffolds_start = {}
relevant_scaffolds_end = {}
relevant_scaffolds = {}
gffsubset = []

with open(file) as f:
    line = f.readline()
    while line:
        scaffold = line.split(",")[0]
        start = line.split(",")[1]
        end = line.split(",")[2]
        both = start + "," + end
        relevant_scaffolds[scaffold] = both
        line = f.readline()

with open(gff_file) as g:
  gline = g.readline()
  gstart = gline.split()[3]
  gend = gline.split()[4]
  while gline:
      gscaffold = gline.split()[0]
      if not gline.split()[8].startswith("ID"):
            gline = g.readline()
            continue
      if gscaffold in relevant_scaffolds:
        window = relevant_scaffolds[gscaffold]
        wstart = window.split(",")[0]
        wend = window.split(",")[1]
        if  wstart < gstart < wend or wstart < gend < wend:
            gffsubset.append(gline.split()[8].split("=")[3][1:][:-2])
      gline = g.readline()

with open(project_name + "genelist.txt","w") as file:
    file.writelines('%s\n' % item for item in gffsubset)
