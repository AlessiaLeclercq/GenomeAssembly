import os 
import numpy as np

# PROBLEM 1
def readDataFromFile(fileName: str) -> dict:
    '''Given a file name read the data which are in the format <sequence name>:<sequence>'''
    reads = dict()
    current_dir = os.getcwd()
    current_file = os.path.join(current_dir, fileName)
    
    with open(current_file, "r") as f:
        for line in f:
            read_name = line.split(":")[0].rstrip()
            sequence = line.split(":")[1].rstrip()
            reads[read_name] = sequence
    f.close()
    return reads


# PROBLEM 2
def meanLength(fileName:str) -> float:
    '''Given a file compute the mean length of the sequences'''
    reads = readDataFromFile(fileName)
    return sum([len(sequence) for sequence in reads.values()])/len(reads)



# PROBLEM 3
def getOverlap(left:str, right:str) -> str:
    '''Computes the overlap between two sequences'''
    max_overlap = min(len(left), len(right))
    for cut_index in range(max_overlap, 0, -1):
        #get the laft and right cuts
        left_cut = left[-cut_index:]
        right_cut = right[:cut_index]
        #compare to check if they overlap 
        if left_cut == right_cut:
            return left_cut
    return ""


# PROBLEM 4
def getAllOverlaps(reads:dict) -> dict:
    '''Find the overlaps between all sequences'''
    overlaps_dict = dict()
    #iterating over sequences
    keys = list(reads.keys())
    #initialize internal dictionaries
    for key in keys:
        overlaps_dict[key] = dict()
        
    #search for all overlaps
    for out_index in range(len(keys)):
        read1 = keys[out_index]
        #compare with the other sequences
        for in_index in range(out_index+1, len(keys)):
            read2 = keys[in_index]
            #using the string as left search
            left_search = len(getOverlap(reads[read1], reads[read2]))
            overlaps_dict[read1][read2] = left_search
            #perfrom the opposite
            right_search = len(getOverlap(reads[read2], reads[read1]))
            overlaps_dict[read2][read1] = right_search

    return overlaps_dict


# PROBLEM 5
def prettyPrint(overlaps: dict) -> None:
    '''Prints the overlaps in matrix form'''
    key_list = sorted(overlaps.keys())
    #display the header 
    print("%3s" %" " , end = ' ')
    for key in key_list:
        print("%3s" %key , end = ' ')
    print()
    
    #display the content
    for key in key_list:
        print("%3s" %key , end = ' ')
        for match in key_list:
            if match == key:
                print("%3s" %"-" , end = ' ')
            else: 
                print("%3d" %overlaps[key][match], end = ' ')
        print()


# PROBLEM 6
def is_significant_left(overlaps: dict, sequence_name:str)-> bool:
    '''Finds the number of significant overlaps (>2) when positioned to the left'''
    return sum([1 for match, ovrlp in overlaps[sequence_name].items() if ovrlp>2])

    
def is_significant_right(overlaps: dict, sequence_name: str)-> bool:
    '''Finds whether a sequence has a significant overlap (>2) when positioned to the right'''
    significant = False
    for match in overlaps.keys():
        if match != sequence_name:
            if overlaps[match][sequence_name] > 2:
                significant = True
    return significant 

def findFirstRead(overlaps: dict) -> str:
    '''Finds the first (leftmost) read, hence the one that only has a significant (>2) overlap to its right end, i.e., 
    it only has a good overlaps when positioned to the left of other reads'''
    for sequence_name in overlaps.keys():
        if is_significant_left(overlaps, sequence_name)==1 and not is_significant_right(overlaps, sequence_name):
            return sequence_name
    return None 


# PROBLEM 7
def findKey(d:dict) -> str:
    '''Finds the read that has the largest overlap to the right end of the current read'''
    return max(list(d.items()), key = lambda item: item[1])[0]


def findOrder(name: str, overlaps: dict) -> list:
    '''Gets the first read and the dictionary of all overlaps and returns a list of read names in the order in 
    which they represent the genomic sequence'''
    readOrder = [name]
    readOrder = recursiveOrder(name, overlaps, readOrder)
    return readOrder


def recursiveOrder(current_read: str, overlaps: dict, readOrder: list) -> list:
    '''Recursive construction of the list of read names representing the genomic sequence'''
    #find the next read (overlapping) and add it to the sequence
    next_read = findKey(overlaps[current_read])
    readOrder.append(next_read)
    #check if next_read has significant overlapping and srecursive call
    if is_significant_left(overlaps, next_read):
        recursiveOrder(next_read, overlaps, readOrder)
    
    return readOrder


# PROBLEM 8
def assembleGenome(readOrder: list, reads: dict, overlaps: dict)-> str:
    '''Given the order, the sequences and the overlaps assemble the genome'''
    genome = reads[readOrder[0]]
    previous = readOrder[0]
    for current in readOrder[1:]:
        overlapping_length = overlaps[previous][current]
        genome = genome + reads[current][overlapping_length : ]
        previous = current
    return genome