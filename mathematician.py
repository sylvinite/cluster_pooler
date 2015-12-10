import re
import pdb
from collections import defaultdict

#Attempt one, just pulls randomly to new bins
def crudeBinning(struct):
    bins = defaultdict(list)
    binCounter = 1
    for k, v in struct.items():   
        newBin = True
        for num in bins:
            if reads < spaceLeft(bins[num]):
                sample = {}
                sample[k] = v
                bins[num].append(sample)
                newBin = False
                break
        if newBin:
            sample = {}
            sample[k] = v
            bins[binCounter].append(sample)
            binCounter += 1
    return bins
        
def spaceLeft(struct): 
    binMax = 320000000
    total = 0
    if not (struct == []):
        for ceg in struct:
           total  += ceg.values()[0].keys()[0]
    return binMax - total

def printer(struct):
    for f in struct.items():
        print "Bin #" , f[0]
        pdb.set_trace()
        print struct.items()

counter = 0
pattern=re.compile("(P[0-9]{3,5}_[0-9]{3,5}) ([0-9]+) ([0-9]+\.[0-9]+) ([A-Z]+)")
binMax = 320000000
fileName = open('./P2652.txt', 'r')
init = dict()

#Initializes data
for line in fileName:
    if pattern.search(line):
        sample = pattern.search(line).group(1) 
        reads = (binMax - int(pattern.search(line).group(2)) )
        lanePerc =  float(pattern.search(line).group(3)) 
        index =  pattern.search(line).group(4) 
        counter+=1
        init[sample] = dict()
        init[sample][reads] = index

xy = crudeBinning(init)
printer(xy)
testa = list()
        
        
    
        


        
        
    