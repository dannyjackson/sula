import sys

file = sys.argv[1]

outputList = []

with open(file) as f:
    line = f.readline()
    while line:
        fdstat = line.split(",")[10]
        if fdstat > 0.3:
            outputList.append(line)
        line = f.readline()

write("\n".join(outputList))
