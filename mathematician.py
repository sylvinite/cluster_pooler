import re
import pdb
from collections import defaultdict
from ecdsa.util import string_to_number
from ecdsa import ecdsa

class Tree(defaultdict):
    def __init__(self, value=None):
        super(Tree, self).__init__(Tree)
        self.value = value

def idealBins(readIndexPair):
    indexer = defaultdict()
    readIt = 0
    binSize = 320000000
    for f in readIndexPair:
        pdb.set_trace()
        readIt += f
        indexer[readIndexPair[f]] == 1
    minBins = readIt/binSize
    return minBins

fileName = open('/Users/isaksylvin/Documents/P2652.txt', 'r')
counter = 0
sampleReadPair = {}
readIndexPair = {}
pattern=re.compile("(P[0-9]{3,5}_[0-9]{3,5}) ([0-9]+) ([0-9]+\.[0-9]+) ([A-Z]+)")
#pattern=re.compile("(P[0-9]{3,5}_[0-9]{3,5})")
while fileName.readline():
    line = fileName.readline()
    if pattern.search(line):
        sample = pattern.search(line).group(1)
        reads = int(pattern.search(line).group(2))
        lanePerc = float(pattern.search(line).group(3))
        index = pattern.search(line).group(4)
        sampleReadPair[sample] = reads
        readIndexPair[reads] = index
        
        counter+=1
print idealBins(readIndexPair)

#Attempt 1: No internal order, no threshold, no initial placement of many index clusters etc.
pants = list()
#for f in readIndexPair:
    
    
        


        
        
    