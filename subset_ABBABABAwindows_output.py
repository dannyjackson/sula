import sys

file = sys.argv[1]

outputList = []

with open(file) as f:
    line = f.readline()
    while line:
        fdstat = line.split(",")[10]
        if fdstat > 0.3:
            outputList.append(line)
            outputList.write("\n")
        line = f.readline()

print outputList
