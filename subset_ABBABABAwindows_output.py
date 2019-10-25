import sys, gzip

file = sys.argv[1] + "_slidingwindows.csv"
project_name = sys.argv[1] + "_slidingwindows.subsetfd.txt"
threshold = sys.argv[2]

outputList = open(project_name, "a")

with gzip.open(file) as f:
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
