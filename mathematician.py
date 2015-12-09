import re
import pdb
from collections import defaultdict
from collections import Counter

def minBins_readspace(reads):
    return sum(reads)/binMax

def minBins_indexIt(index):
    count = Counter()
    for tag in index:
        count[tag] += 1
    return count

""" Attempt 1: No internal order, no threshold, no initial placement of many index clusters etc.  
"""  
def crude_Binner(sample, reads, index):
    if minBins_indexIt(index) > minBins_readspace(reads):
        initialBins = minBins_indexIt(index)
    else:
        initialBins = minBins_readspace(reads)
    
    binList = list()
    bin = list()
    n = 0
    while n < initialBins:
        binList.append(bin)
    


fileName = open('./P2652.txt', 'r')
counter = 0
binMax = 320000000
sample = list()
reads = list()
lanePerc = list() 
index = list()
pattern=re.compile("(P[0-9]{3,5}_[0-9]{3,5}) ([0-9]+) ([0-9]+\.[0-9]+) ([A-Z]+)")
for line in fileName:
    if pattern.search(line):
        sample.append( pattern.search(line).group(1) )
        reads.append(binMax - int(pattern.search(line).group(2)) )
        lanePerc.append( float(pattern.search(line).group(3)) )
        index.append( pattern.search(line).group(4) )
        counter+=1

print minBins_readspace(reads)
print minBins_indexIt(index)

    
        


        
        
    