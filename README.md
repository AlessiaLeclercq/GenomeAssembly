# GenomeAssembly
An approach to genome assembly is presented alligning short reads (sequences) into long ones (chromosomes). \
The tasks consists in the following: 
  - Read the sequences from a file in the format sequence_name:sequence and store it into a dictionary. 
  - Compute the mean length of the reads.
  - Find the maximal length overlaps between reads assuming that no reads are completely nested one into another. Given a left and a right reads, in order to overlap the the 3’ end of the left read must match the 5’ end of the right read. Furthermore we assume that there are no sequencing errors, therefore the sequence match is perfect. 
  - Find the maximal length overlaps between the reads in both orientations and store them in a dictionary of dictionaries in the format $d[left][right] = length$.
  - Print in matrix format the content of the previous dictionary.
  - Identify the first leftmost read. It is the one that only has a significant (>2) overlap to its right, i.e., it only has a good overlap when positioned to the left of other reads and has no good overlaps when positioned to the right. 
  - Recursively find the correct order of the reads. Starting from the read identified in the previous point, search for the largest ovelapping read among the available ones and keep on going until a read that has no significant ovelaps (when positioned to the left) is found. 
  - Reconstruct the genomic sequence given the order of the reads, the overlapping sequences and the reads.
  
The code has been divided into [functions](https://github.com/AlessiaLeclercq/GenomeAssembly/blob/main/utils.py) implementing the steps and a [main](https://github.com/AlessiaLeclercq/GenomeAssembly/blob/main/main.py) simulating the result. 
