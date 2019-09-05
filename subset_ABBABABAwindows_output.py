import sys
import csv

file = sys.argv[1]

outputList = csv.writer()

with open(file) as f:
    line = f.readline()
    while line:
        fdstat = line.split(",")[10]
        if fdstat > 0.3:
            outputList.writerow(line)
        line = f.readline()
