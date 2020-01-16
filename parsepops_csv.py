import sys, subprocess

file = sys.argv[1]

arguments = ""
with open(file) as f:
    line = f.readline()
    while line:
        pop = line.split()
        arguments = arguments + pop + ","
        arguments = arguments[:-1]
        line = f.readline()

print arguments
