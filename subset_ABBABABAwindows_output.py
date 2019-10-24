import sys

file = sys.argv[0]
project_name = sys.argv[1] + ".subsetabbababa.txt"
threshold = sys.argv[2]

outputList = open(project_name, "a")

    with open(file) as f:
        line = f.readline()
        while line:
            fdstat = line.split(",")[0]
            if fdstat > threshold:
                outputList.write(line)
            line = f.readline()

outputList.close()
