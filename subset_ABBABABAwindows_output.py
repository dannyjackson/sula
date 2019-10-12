import sys

file = sys.argv[1]
project_name = sys.argv[2] + ".subsetabbababa.txt"
threshold = sys.argv[3]

outputList = open(project_name, "a")

with open(file) as f:
    line = f.readline()
    while line:
        fdstat = line.split(",")[9]
        if fdstat > threshold:
            outputList.write(line)
        line = f.readline()

with open(file) as f:
    line = f.readline()
    while line:
        fdstat = line.split(",")[9]
        if ( fdstat == "fd" ):
            outputList.write(line)
        if ( fdstat != "fd" ):
          if ( float(fdstat) > float(threshold) ):
              outputList.write(line)
        line = f.readline()

outputList.close()
