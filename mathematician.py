import re
from collections import defaultdict, Counter
import pdb

binMax = 320000000

""" Initializes data from file
    Caution: File composition is hardcoded
"""
def readdata(file):
    initialBin = defaultdict(set)
    pattern=re.compile("(P[0-9]{3,5}_[0-9]{3,5}) ([0-9]+) ([0-9]+\.[0-9]+) ([A-Z]+)") 

    for line in file:
        if pattern.search(line):
            sample = pattern.search(line).group(1) 
            reads = (binMax - int(pattern.search(line).group(2)) )
            lanePerc =  float(pattern.search(line).group(3)) 
            index =  pattern.search(line).group(4) 
            initialBin[sample] = [reads, index]

    return initialBin

""" Pools data into bins in the most mundane way
    Caution: File composition is hardcoded
"""
def simpleBinner(init):
    bin = defaultdict(set)
    
    for k, v in init.items():
        for num in xrange(1, 200):
            if bin[num] == set([]):
                bin[num] = dict()
            if v[0] < 0:
                break
            if v[0] < avail_bin_space(bin[num]) and not tag_present(v[1], bin[num]):
                bin[num][k] = v
                break
    return bin

""" Checks if specific tag is present in a bin
    Caution: value[1] is hardcoded to represent tag
"""
def tag_present(tag, bin):
    for k, v in bin.items():
        if v[1] == tag:
            return True
    return False
""" Calculates available space left in a bin
    Caution: value[0] is hardcoded to represent (binMax - clusterSize)
"""    
def avail_bin_space(bin):
    sum = 0
    if not bin == {}:
        for k, v in bin.items():
            if k == 'P2652_157':
                pdb.set_trace()
            sum += v[0]
        return binMax - sum
    else:
        return binMax

""" Prints info for each sample in each bin
    Caution: Most value positions are hardcodedd
"""    
def bin_printer(bin):
    for num in bin:
        
        print '***Bin #' , num ,'***'
        for item in bin[num].items():
            print item[0] ,"{0:.2f}".format((float)(item[1][0])/1000000),'M', item[1][1]
        print
def bin_stats(bin):
    tag_its = Counter()
    emptySpace = dict()
    totalBins = bin.keys()[-1]
    
    for num in bin.keys():
        emptySpace[num] = (avail_bin_space(bin[num])/(binMax/1000000)) 
    spaceUsed = (totalBins*binMax) - sum(emptySpace.values())
    #Odd calculation to avoid truncation
    minBins_size = (float)(spaceUsed/(binMax/100))/100
    for binNo in bin.values():
        for tag in binNo.values():
            barcode = tag[1]
            tag_its[barcode] += 1
        
    print 'Total bins: ', totalBins
    print 'Minimal amount of bins (size): ', minBins_size
    print 'Most common tags: ' , tag_its.most_common()
    print 'Empty space per bin (ppm): ', emptySpace.values()
    print


        
    
file = open('.\P2652.txt', 'r')
bin = readdata(file)
output = simpleBinner(bin)
bin_stats(output)
bin_printer(output)

