from utils import *

#PROBLEM 1
fileName = 'genome-assembly.txt'
reads = readDataFromFile(fileName)

#PROBLEM 2
print(meanLength(fileName))

#PROBLEM 3 & 4
overlaps = getAllOverlaps(reads)

#PROBLEM 5
prettyPrint(overlaps)

#PROBLEM 6
name = findFirstRead(overlaps)

#PROBLEM 7
order = findOrder(name, overlaps)

#PROBLEM 8
genome = assembleGenome(order, reads, overlaps)
print(genome)