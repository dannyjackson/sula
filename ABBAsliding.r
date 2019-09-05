#!/usr/bin/env Rscript
#Genome wide ABBA BABA test with block jackknife procedue

#a file that is accessed by abba_complete.sh

args = commandArgs()

outputdirectory = substr(args[grep("outputdirectory_", args)],17,100000)
population1 = substr(args[grep("population1_", args)],13,100000)
population2 = substr(args[grep("population2_", args)],13,100000)
population3 = substr(args[grep("population3_", args)],13,100000)
