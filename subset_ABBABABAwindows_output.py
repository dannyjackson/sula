import sys

file = sys.argv[1]
project_name = sys.argv[2] + ".txt"

outputList = open(project_name, "a")

with open(file) as f:
    line = f.readline()
    while line:
        fdstat = line.split(",")[10]
        if fdstat > 0.3:
            outputList.write(line)
        line = f.readline()

outputList.close()
