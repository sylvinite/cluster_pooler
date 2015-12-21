import re
from collections import defaultdict, OrderedDict, Counter
import time
import pdb


minClust = 320000000
clust_per_lane = 320000000

""" Initializes data from file
    Caution: File composition is hardcoded
"""
def readdata(file):
    initialBin = defaultdict(set)
    pattern=re.compile("(P[0-9]{3,5}_[0-9]{3,5}) ([0-9]+) ([0-9]+\.[0-9]+) ([A-Z]+)") 

    for line in file:
        if pattern.search(line):
            sample = pattern.search(line).group(1) 
            remain = (minClust - int(pattern.search(line).group(2)) )
            lanePerc =  float(pattern.search(line).group(3)) 
            index =  pattern.search(line).group(4)
            if (remain > 0): 
                initialBin[sample] = [remain, index]

    return initialBin

""" Pools data into bins in the most mundane way
    Uses small first (to min splits and bins)
    Caution: File composition is hardcoded
"""
def simpleBinner(init):
    #Small first
    init = OrderedDict(sorted(init.iteritems(),key=lambda (k,v): v[0],reverse=False))
    bin = defaultdict(set)
    maxBins = 200
    
    for k, v in init.items():
        for num in xrange(1, maxBins):
            #If bin is empty
            if num > len(bin):
                bin[num] = dict()
            if not tag_present(v[1], bin[num]) and suf_reads(v[0], bin[num]):
                bin[num][k] = v
                break
    return bin

""" Pools data into bins by minimizing overexpression
    Finds ideal bin for any sample; never splits; goes through elements by least OE first
    Caution: File composition is hardcoded
"""
def greedyBinnerDeep(init):
    bin = defaultdict(set)
    #Big first, naturally has least overexpression
    maxBins = 200
    #bin index, overexpression
    while not init == {}:
        least_oe = clust_per_lane
        least_oe_index = 0
        for k, v in init.items():
            for num in xrange(1, maxBins):
                #If we reached empty bin, add number and break loop
                if num > len(bin):
                    emptyBinIndex = num
                    if clust_per_lane - v[0] < least_oe:
                        oe_key = k 
                        oe_value = v 
                        least_oe = clust_per_lane - v[0]  
                    break
                #Check if item fits in bin
                elif not tag_present(v[1], bin[num]) and suf_reads(v[0], bin[num]):
                        #Check if adding the item here gives the least OE, and add info in case
                        if clust_per_lane/len(bin[num].keys())+1 - v[0] < least_oe:
                            oe_key = k 
                            oe_value = v 
                            least_oe_index = num
                            least_oe = clust_per_lane/len(bin[num].keys())+1 - v[0]  
        #Adds item
        if least_oe_index == 0:
            bin[emptyBinIndex] = dict()
            if not tag_present(v[1], bin[num]) and suf_reads(v[0], bin[num]):
                    bin[emptyBinIndex][oe_key] = oe_value
                    del init[oe_key]
            else:
                print "Error, algo should absolutely never reach here"
                pdb.set_trace()
        else:
            bin[least_oe_index][oe_key] = oe_value
            del init[oe_key]
            
    return bin


""" Takes elements in descending order and 
    adds to bin until overexpression for first sample hits negative
    Reinserts the remainders into new bins as sub-elements in order of size
    Nicknamed confident since it assumes ending OE won't overshadow the rest
    Caution: File composition is hardcoded
"""
def confidentBinner(init, threshold):
    #Big first (OD is actually backwards, but pop works correctly
    init = OrderedDict(sorted(init.iteritems(),key=lambda (k,v): v[0],reverse=False))
    bin = defaultdict(set)
    maxBins = 200
    
    while not init == {}:
        init = OrderedDict(sorted(init.iteritems(),key=lambda (k,v): v[0],reverse=False))
        k, v = init.popitem()
        for num in xrange(1, maxBins):
            #If bin is empty
            if num > len(bin):
                bin[num] = dict()
            if not tag_present(v[1], bin[num]) and suf_reads(v[0], bin[num]):
                bin[num][k] = v
                expressed = clust_per_lane/len(bin[num].keys())
                if expressed < v[0]:
                        init[k] = [v[0] - expressed, v[1]]
                break
            #This if statement is muy importante to be right
            elif not tag_present(v[1], bin[num]) and decreases_oe(v[0], bin[num], threshold):
                bin[num][k] = v
                #Retain all values that weren't pushed into OE, but subtracts the reads entered
                expressed = clust_per_lane/len(bin[num].keys())
                for item in bin[num].items():
                    if expressed < item[1][0]:
                        init[item[0]] = [item[1][0] - expressed, item[1][1]]
                break
                    
                    
    return bin

""" Tallies OE for a lane/bin
    :param threshold Percentage OE of clust_per_lane that this function is ok with
    :return Whether adding the element to the lane/bin reduces lane/bin OE
"""
def decreases_oe(reads, bin, threshold):
    existing_samples = len(bin.keys())
    init_oe = 0
    post_oe = 0
    if bin == {}:
        return False
    else:
        for k, v in bin.items():
            if (clust_per_lane/existing_samples - v[0] > 0):
                init_oe += clust_per_lane/existing_samples - v[0]
        for k, v in bin.items():
            if (clust_per_lane/(existing_samples + 1) - v[0] > 0):
                post_oe += clust_per_lane/(existing_samples + 1) - v[0]
        if (clust_per_lane/(existing_samples + 1) - reads > 0):
            post_oe += clust_per_lane/(existing_samples + 1) - reads
        
        if (post_oe < init_oe - threshold*clust_per_lane):
            return True
        else:
            return False

""" Checks if specific tag is present in a bin
    Caution: value[1] is hardcoded to represent tag
"""
def tag_present(tag, bin):
    
    for k, v in bin.items():
        if v[1] == tag:
            return True
    return False
    
""" Calculates how much each sample is overflowing and adds it to the print structure
"""
def read_overflow(bin):
    for num in xrange(1, len(bin.keys())+1):
        samples = len(bin[num].keys())
        clust_per_sample = clust_per_lane/samples
        
        for k, v in bin[num].items():
            if isinstance(v, tuple):
                pdb.set_trace()
            #Pos = Overexpressed, negative = Needs more to pass
            if clust_per_sample - v[0] > 0:
                v.append(clust_per_sample - v[0])
            else:
                v.append(0)
    return bin

""" Calculates if a sample reaches its threshold in the current bin
"""
def suf_reads(thres, bin):
    if bin == {}:
        if thres < clust_per_lane:
            return True
        else:
            print "Sample needs more than a single lane; unhandled exception!"
            pdb.set_trace()
            return False
    else:
        #If space for own sample read requirements
        if thres < clust_per_lane/(len(bin.keys())+1):
            for k, v in bin.items():
                #If adding doesnt break an existing sample read requirements
                if v[0] > clust_per_lane/(len(bin.keys())+1):
                    return False
            return True
    return False

""" Prints info for each sample in each bin
    Now hardlocked for CONFIDENT
    Caution: Most value positions are hardcoded
""" 
def gen_stats(bin):
    tag_its = Counter()
    clusters = 0

    
    for tag in bin.values():
        barcode = tag[1]
        tag_its[barcode] += 1
        clusters += tag[0]
        

    print 'Reads per lane:', clust_per_lane/1000000, 'M'
    print 'Sample reads requirement:', minClust/1000000, 'M'
    print 'Most common tags (descending):' , tag_its.most_common()
    print 'Lowest possible bins:', clusters/clust_per_lane
    print
    
""" Prints info for each sample in each bin
    Now hardlocked for CONFIDENT
    Caution: Most value positions are hardcoded
"""    
def bin_printer(bin):
    for num in bin:
        print '***Bin #' , num ,'***'
        print 'Sample Reads_Prior TAG Expressed'
        for item in bin[num].items():
            print item[0] ,"{0:.2f}".format((float)(item[1][0])/1000000),'M', item[1][1], clust_per_lane/(len(bin[num].keys())*1000000), 'M'
        print

""" Prints stats after an algorithm has been ran. Probably pretty prone to break
"""
def algo_stats(init, bin):
    tag_its = Counter()
    clusters = 0

    
    for tag in init.values():
        barcode = tag[1]
        tag_its[barcode] += 1
        clusters += tag[0]

    bins = bin.keys()[-1]
    algo_clusters = (bins*clust_per_lane)/1000000
    print "Bins used:", bins, '; Ideal:', clusters/clust_per_lane
    print "Overexpression:", algo_clusters, 'M', "OR", "{0:.2f}".format(float(bins*100)/(clusters/clust_per_lane)-100), "%"

file = open('.\P2652.txt', 'r')
bin = readdata(file)
gen_stats(bin)

def confidentCalibration(bin):
    best_thres = 1000000
    worst_thres = 0
    least_bins = 1000000
    most_bins = 0
    
    for thres in xrange(0, 8000, 50):
        thres = float(thres)/1000
        output = confidentBinner(bin, thres)
        if output.keys()[-1] < least_bins:
            best_thres = thres
            least_bins = output.keys()[-1]
        elif output.keys()[-1] > most_bins:
            worst_thres = thres
            most_bins = output.keys()[-1]
        
    #return {best_thres, least_bins, worst_thres, most_bins}
    return best_thres
        
    
"""
print "SIMPLE"
simple_start = time.time()
simple_output = simpleBinner(bin)
simple_end = time.time()
simple_output = read_overflow(simple_output)
bin_printer(simple_output)
print("Simple time:", "%s seconds" % (simple_end - simple_start))
"""
"""
#Ironically works almost identical to simple. Without splitting its impossible to get good output
print "GREEDY"
greedy_start = time.time()
greedy_output = greedyBinnerDeep(bin)
greedy_end = time.time()
greedy_output = read_overflow(greedy_output)
bin_printer(greedy_output)
print("Greedy time:", "%s seconds" % (greedy_end - greedy_start))
"""

print "CONFIDENT"
conf_start = time.time()
threshold = confidentCalibration(bin)
conf_output = confidentBinner(bin, threshold)
conf_end = time.time()
bin_printer(conf_output)
print("Confident time:", "%s seconds" % (conf_end - conf_start))
algo_stats(bin, conf_output)


"""
bin_stats(simple_output)
bin_stats(greedy_output)


"""

