import sys, subprocess

file = sys.argv[1]

uniquePops = set()
uniquePopsList = []
arguments = ""
with open(file) as f:
    line = f.readline()
    while line:
        pop = line.split()[1]
        if pop not in uniquePops:
            uniquePopsList.append(pop)
            arguments = arguments + " -p " + pop
        uniquePops.add(pop)
        line = f.readline()

print arguments
