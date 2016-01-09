import re
from collections import defaultdict, OrderedDict, Counter
import time
import pdb
import sys

sys.setrecursionlimit(10000) # Overrides recursion limit

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
    if(minClust > clust_per_lane):
        print "Sample needs more than a single lane, use a more complex algo"
        pdb.set_trace()
        return init
    
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
    if(minClust > clust_per_lane):
        print "Sample needs more than a single lane, use a more complex algo"
        pdb.set_trace()
        return init

    bin = defaultdict(set)
    #Big first, naturally has least overexpression
    maxBins = 2000
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
    if(minClust > clust_per_lane):
        print "Sample needs more than a single lane, use a more complex algo"
        pdb.set_trace()
        return init

    bin = defaultdict(set)
    maxBins = 200
    magic_number = 4
    tag_its = Counter()
    
    for tag in init.values():
        tag_its[tag[1]] += 1
    unique = len(tag_its)
    
    while not init == {}:
        #Big first (OD is actually backwards but pop gathers biggest)
        init = OrderedDict(sorted(init.iteritems(),key=lambda (k,v): v[0],reverse=False))
        k, v = init.popitem()
        for num in xrange(1, maxBins):
            #If bin is empty
            if num > len(bin):
                bin[num] = dict()
            if not tag_present(v[1], bin[num]) and suf_reads(v[0], bin[num]):
                bin[num][k] = v
                expressed = clust_per_lane/float(len(bin[num].keys()))
                break
            #broken_samples works better than decreases_oe
            elif not tag_present(v[1], bin[num]) and broken_samples(v[0], bin[num]) < magic_number:
                bin[num][k] = v
                #Retain all values that weren't pushed into OE, but subtracts the reads entered
                expressed = clust_per_lane/float(len(bin[num].keys()))
                for item in bin[num].items():
                    if expressed < item[1][0]:
                        init[item[0]] = [item[1][0] - expressed, item[1][1]]
                break
                    
                    
    return bin

""" Scans a bin for unique tags
"""
def unique_tags(bin):
    tag_its = Counter()
    for tag in bin.values():
        tag_its[tag[1]] += 1
    unique = len(tag_its)
    return unique

""" Creates an amount of bins equal to unique barcodes and place reads into sets that fits the sample.
    If difference is greater than offset, attempts to split the sample into smaller bins to resolve it.
    If a bin cant be filled with samples they are reevaluated
"""
def divide_n_conquer(init, offset):
    bin = defaultdict(set)
    output = defaultdict(set)
    #Scans for unique tags
    unique = unique_tags(init)
    
    #Creates all necessary bins
    #xrange never actually hits maximum
    for num in xrange(1,unique+1):
        bin[num] = dict()
    
    while not init == {}:
        init = OrderedDict(sorted(init.iteritems(),key=lambda (k,v): v[0],reverse=False))
        k, v = init.popitem()
        dividers = searchspace(v[0], offset, unique)
        if not dividers:
            for minBin in xrange(1, unique+1):
                if clust_per_lane/minBin < v[0]:
                    minBin -= 1
                    break
            bin[minBin][k] = v
        else:
            z = 0
            for n in dividers:
                bin[n][k] = [v[0] - float(z), v[1]]
                z += clust_per_lane/float(n)
                
    #Tally into actual bins
    currentBin = 1
    for n in bin:
        bin[n] = OrderedDict(sorted(bin[n].iteritems(),key=lambda (k,v): v[0],reverse=True))
        while len(bin[n]) >= n and unique_tags(bin[n]) >= n:
            counter = 0
            output[currentBin] = dict()
            while counter < n:
                for k,v in bin[n].items():
                    if not tag_present(v[1], output[currentBin]):
                        output[currentBin][k] = v
                        del bin[n][k]
                        counter += 1
                        break
            currentBin +=1
            
    #Tally remains into bins
    #Could be improved with new offset and unique param
    for boxNum in bin:
        status = True
        while not bin[boxNum] == OrderedDict():
            output[currentBin] = dict() #May create an unnecessary dict. Fix later in that case
            #If output can be created from a box (or subset), push to output
            if unique_tags(bin[boxNum]) >= boxNum:
                for k, v in bin[boxNum].items():
                    if not tag_present(v[1], output[currentBin]) and len(output[currentBin]) < boxNum:
                        #Push to output
                        output[currentBin][k] = v
                        del bin[boxNum][k]
                currentBin+=1
                
            #If all elements can be split
            #Check that operations work
            elif status == True:
                for k, v in bin[boxNum].items():
                    if boxNum + 3 < 13:
                        q = searchspace_inner(clust_per_lane/float(boxNum), offset, 2, list(xrange(boxNum+1, boxNum+3)), unique)
                        if q == False:
                            status = False
                        else:
                            #Checks that a proposed spot isnt occupied by the same sample
                            for numbers in q:
                                if k in bin[numbers].keys():
                                    status = False
                                    break
                    else: 
                        status = False
                #Actually split the items down        
                if status == True:
                    while not bin[boxNum] == OrderedDict():
                        k, v = bin[boxNum].popitem()
                        divs = searchspace_inner(clust_per_lane/float(boxNum), offset, 2, list(xrange(boxNum+1, boxNum+3)), unique)
                        z = 0
                        for num in divs:
                            #Push to bin
                            bin[num][k] = [v[0] - z, v[1]] #Copy to bigger div, right hand is to copy corret value-part
                            z += clust_per_lane/float(num)
            
            #Else, recruit a single element from smaller boxes
            else:
                plusOne = boxNum
                still_running = True
                while plusOne < unique and still_running:
                    plusOne += 1
                    for k, v in bin[plusOne].items():
                        if not tag_present(v[1], bin[boxNum]):
                            #Push to bin
                            bin[boxNum][k] = v
                            del bin[plusOne][k]
                            still_running = False
                            break
                        
                #If no elements can be recruited then still_running is True!
                if still_running:
                    counter = 0
                    #Nothing can be done. Partion (sub)box to an output bin
                    for k,v in bin[boxNum].items():
                        if not tag_present(v[1], output[currentBin]):
                            #Push to output
                            output[currentBin][k] = v 
                            del bin[boxNum][k]
                            counter += 1
                        if counter == boxNum:
                            break
                    currentBin+=1
    return output

""" Fits a read length into a sum of integer divisions, where every integer is unique and
    the maximum integer is equal to the amount of unique barcodes for the dataset
    :param numToPlace Read to fit into the sum
    :param offset Maximum padding around each read (partial overexpression)
    :param unique_barcodes Unique barcodes
    :return A list of integers to divide the clust_per_lane by
"""
def searchspace(numToPlace, offset, unique_barcodes):
    chunks = 1
    positions = list(xrange(1, chunks+1))
    pos = searchspace_inner(numToPlace, offset, chunks, positions, unique_barcodes)
    return pos

""" Inner function for searchspace
    :param numToPlace Read to fit into the sum
    :param offset Maximum padding around each read (partial overexpression)
    :param chunks Minimum dividers to use
    :param posits Starting dividers
    :param max_chunks Unique barcodes
    :return A list of integers to divide the clust_per_lane by
"""
#remove chunks parameter, can be replaced by len(posits) or similar
def searchspace_inner(numToPlace, offset, chunks, posits, max_chunks):
    #Should really search in a better way than [1]->[12]->[1,2] to find minimal
    #Checks if goal was reached
    sum = 0
    for value in posits:
        sum += clust_per_lane/float(value)
    if numToPlace < sum and numToPlace*(1+offset) >= sum:
        return posits
    elif numToPlace < clust_per_lane/float(max_chunks):
            posits[-1] = max_chunks
            return posits
    else:
        for holes in xrange(1, chunks+1):
            #Moves last box
            if holes == 1 and posits[-1*holes] < max_chunks:
                posits[-1*holes] += 1
                return searchspace_inner(numToPlace, offset, chunks, posits, max_chunks)
            #Moves not-last box and resets rest
            elif len(posits) > 1 and posits[-1*holes] < posits[-1*(holes-1)] -1:
                posits[-1*holes] += 1
                
                for spots in xrange(holes-1, 0, -1):
                    #Should only increase if its possible to increase the spot
                    if posits[-1*spots -1] < posits[-1*spots]:
                        posits[-1*spots] = posits[-1*holes] + 1
                return searchspace_inner(numToPlace, offset, chunks, posits, max_chunks)
        #We reached the end and want to split
        if chunks < max_chunks:
            if numToPlace < sum:
                # Error, want to split but split increases size"
                return False
            return searchspace_inner(numToPlace, offset, chunks+1, list(xrange(1, chunks+2)), max_chunks)
        #We reached the complete end
        else:
            print "Error, no solution found"
            return False

def bin_filled(bin):
    bin_filled = 0
    for k,v in bin.items():
        bin_filled += v[0]
    return bin_filled

""" Checks how many samples in a bin a new insert breaks
    :param read new insert
    :return broken samples (including new one)
"""
def broken_samples(read, bin):
    broken_samples = 0
    new_size = clust_per_lane/(len(bin.keys())+1)
    for k, v in bin.items():
        if new_size < v[0]:
            broken_samples += 1
    if read < new_size:
            broken_samples += 1
    return broken_samples

""" Checks to see if amount of unique tags changed from data that previously has been multi-binned
"""       
def remaining_splits(bin):
    rem_pool = defaultdict(set)
    for n in xrange(1,13):
        while not bin[n] == {}:
            k,v = bin[n].popitem()
            rem_pool[k] = v
    rem_tags = unique_tags(rem_pool)  
    rem_pool = OrderedDict(sorted(rem_pool.iteritems(),key=lambda (k,v): v[0],reverse=False))
    return rem_tags

""" Tallies OE for a lane/bin
    :param threshold Percentage OE of clust_per_lane that this function is ok with
    :return Whether adding the element to the lane/bin reduces lane/bin OE
"""
def decreases_oe(reads, bin, threshold, unique):
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
            #OE can actually never be less than clust_per_lane/unique tags
            elif v[0] - clust_per_lane/(existing_samples + 1)  < clust_per_lane/unique:
                post_oe += clust_per_lane/unique - (v[0] - clust_per_lane/(existing_samples + 1))
        if (clust_per_lane/(existing_samples + 1) - reads > 0):
            post_oe += clust_per_lane/(existing_samples + 1) - reads
        
        if (post_oe < init_oe - threshold*clust_per_lane):
            return True
        else:
            return False

""" Checks if specific tag is present in a bin
    Caution: value[1] is hardcoded to represent tag
    :param tag barcode to look for
    :param bin bin where the tag may exist
    :return presense of barcode
"""
def tag_present(tag, bin):
    if bin == set([]):
        print "Error, input bin is just a set. Assign to dict"
        pdb.set_trace()
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
    print 'Unique tags:', len(tag_its)
    print 'Lowest possible bins (size alone):', clusters/clust_per_lane
    print
    
""" Prints info for each sample in each bin
    Caution: Most value positions are hardcoded
"""    
def bin_printer(bin):
    for num in bin:
        used_space = 0
        print '***Bin #' , num ,'***'
        print 'Sample Reads_Rem TAG Expressed'
        for item in bin[num].items():
            print item[0] ,"{0:.2f}".format((float)(item[1][0])/1000000),'M', item[1][1], "{0:.2f}".format(float(clust_per_lane/1000000)/(len(bin[num].keys()))), 'M'
            if item[1][0] < clust_per_lane/(len(bin[num].keys())):
               used_space += item[1][0]
            else:
                used_space += clust_per_lane/(len(bin[num].keys()))
        print '--- --- --- --- --- ---'
        print 'Items:', len(bin[num].keys()), '"Unused" space:', "{0:.2f}".format(1-(float)(used_space)/clust_per_lane)
        print

""" Prints stats after an algorithm has been ran
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
    print "Overexpression:", algo_clusters, 'M', "OR", "{0:.2f}".format(float(bins)/(clusters/clust_per_lane)), "OPT"
    
def confidentCalibration(bin):
    best_thres = 1000000
    worst_thres = 0
    least_bins = 1000000
    most_bins = 0
    
    print "Calibrating Confident..."
    
    for thres in xrange(0, 8000, 50):
        thres = float(thres)/1000
        output = confidentBinner(bin, thres)
        if output.keys()[-1] < least_bins:
            best_thres = thres
            least_bins = output.keys()[-1]
        elif output.keys()[-1] > most_bins:
            worst_thres = thres
            most_bins = output.keys()[-1]
        
    print "Calibration complete!"
    #return {best_thres, least_bins, worst_thres, most_bins}
    return best_thres

def dnc_calibration(bin):
    best_thres = 1000000
    worst_thres = 0
    least_bins = 1000000
    most_bins = 0
    
    print "Calibrating DnC..."
    
    for thres in xrange(1, 1000, 2):
        thres = float(thres)/100
        print thres
        output = divide_n_conquer(bin, thres)
        if output.keys()[-1] < least_bins:
            best_thres = thres
            least_bins = output.keys()[-1]
        elif output.keys()[-1] > most_bins:
            worst_thres = thres
            most_bins = output.keys()[-1]
        
    print "Calibration complete!"
    #return {best_thres, least_bins, worst_thres, most_bins}
    return best_thres

file = open('.\P2652.txt', 'r')
minClust = 320000000
clust_per_lane = 320000000

bin = readdata(file)
gen_stats(bin)

print "FIRST FIT ALGO"
simple_start = time.time()
simple_output = simpleBinner(bin)
simple_end = time.time()
simple_output = read_overflow(simple_output)
bin_printer(simple_output)
print("Simple time:", "%s seconds" % (simple_end - simple_start))

"""
#Ironically works almost identical to simple. Without splitting its impossible to get good output
#Also has some kind of destruction issue
print "GREEDY"
greedy_start = time.time()
greedy_output = greedyBinnerDeep(bin)
simple_output = read_overflow(greedy_output)
greedy_end = time.time()
bin_printer(greedy_output)
print("Greedy time:", "%s seconds" % (greedy_end - greedy_start))
pdb.set_trace()
"""
print "CONFIDENT"
conf_start = time.time()
#threshold = confidentCalibration(bin)
threshold = 0.05
conf_output = confidentBinner(bin, threshold)
conf_end = time.time()
bin_printer(conf_output)
print("Confident time:", "%s seconds" % (conf_end - conf_start))
print("Threshold:", threshold)

print "DIVIDE AND CONQUER"
dnc_start = time.time()
#threshold = dnc_calibration(bin)
threshold = 0.1
dnc_output = divide_n_conquer(bin, threshold)
dnc_end = time.time()
bin_printer(dnc_output)
print("DnC time:", "%s seconds" % (dnc_end - dnc_start))
print("Threshold:", threshold)


print "SIMPLE STATS"
algo_stats(bin, simple_output)
print "CONF STATS"
algo_stats(bin, conf_output)
print "DNC STATS"
algo_stats(bin, dnc_output)

