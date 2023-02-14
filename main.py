from utils import *

fileName = 'genome-assembly.txt'
reads = readDataFromFile(fileName)
print(meanLength(fileName))
overlaps = getAllOverlaps(reads)
prettyPrint(overlaps)
name = findFirstRead(overlaps)
order = findOrder(name, overlaps)
genome = assembleGenome(order, reads, overlaps)
print(genome)