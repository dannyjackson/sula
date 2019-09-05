import sys
import csv

file = sys.argv[1]

outputList = []

with open(file) as csv_file:
    line = f.readline()
    while line:
        fdstat = line.split(",")[10]
        if fdstat > 0.3:
            outputList.writerow(line)
        line = f.readline()
