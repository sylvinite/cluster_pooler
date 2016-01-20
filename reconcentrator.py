import pdb
import json
import couchdb
import re
import math
from collections import defaultdict
import unicodedata
import csv
import copy

#Assumes ind. sample conc measurements have failed. As such it relies on changing relative volume on already normalized samples and structure
#Structure are retained as conc measurements failure means there's no way to know conc. delta between samples from seperate poolss

def connection():
    couch = couchdb.Server('http://isak:Purpleplant89@tools-dev.scilifelab.se:5984')
    return couch

#Fetches the structure of a project
def proj_struct(couch, project):
    db = couch['x_flowcells']
    view = db.view('names/project_ids_list')
    fc_track = defaultdict(set)
    
    #Adds flowcells to ALL projects. Due to intractions its easier to just get FCs for ALL projects
    for rec in view.rows:
        fc = ''.join(rec.key)
        fc = unicodedata.normalize('NFKD', fc).encode('ascii','ignore')
        id = ''.join(rec.id)
        id = unicodedata.normalize('NFKD', id).encode('ascii','ignore')
        for projs in rec.value:
            projs = ''.join(projs)
            projs = unicodedata.normalize('NFKD', projs).encode('ascii','ignore')
            if fc_track[projs] == set([]):
                fc_track[projs] = dict()
            fc_track[projs][fc] = id
            
    #Adds lanes and samples to flowcells, includes samples from other projects if they share lane
    for fc, id in fc_track[project].items():
        entry = db[id]['illumina']['Demultiplex_Stats']['Barcode_lane_statistics']
        for index in xrange(0, len(entry)):
            lane = entry[index]['Lane']
            sample = entry[index]['Sample']
            if 'Clusters' in entry[index]:
                clusters = entry[index]['Clusters']
            else: 
                clusters = entry[index]['PF Clusters']
            clusters = int(re.sub(r",", "", clusters))
            
            
            if not isinstance(fc_track[project][fc], dict):
                fc_track[project][fc] = dict()
            if not lane in fc_track[project][fc]:
                fc_track[project][fc][lane] = dict()
            #Only counts samples for the given project, other samples are "auto-filled"
            if project in sample:
                fc_track[project][fc][lane][sample] = clusters
            else:
                fc_track[project][fc][lane][sample] = target_clusters
    #Removes any lanes that don't have any part project samples
    for fc, lanes in fc_track[project].items():
        for lane,sample in lanes.items():
            if not any(project in s for s in sample.keys()):
                   del fc_track[project][fc][lane]
    return fc_track[project]

def aggregator(struct, target_clusters):
    clusters_rem = dict()
    clusters_expr = dict()
    lane_maps = dict()
    ideal_ratios = dict()
    acc_ratios = dict()
    req_lanes = dict()
    
    #Calculate clusters read per sample
    for fc, lanes in struct.items():
        for lane, samples in lanes.items():
            for sample, value in samples.items():
                if not sample in clusters_rem:
                    clusters_rem[sample] = target_clusters
                    clusters_expr[sample] = 0
                clusters_rem[sample] -= value
                clusters_expr[sample] += value
        
    #Concat structure into unique pools
    counter = 1
    for fc, lanes in struct.items():
        for lane, samples in lanes.items():
            mapping = sorted(samples.keys(), reverse=True)
            if not mapping in lane_maps.values():
                lane_maps[counter] = mapping 
                counter +=1

    #Gives how many percent of the lane should give clusters for a specific sample
    for index in lane_maps:
        summ = 0
        for entry in lane_maps[index]:
            if clusters_rem[entry] > 0:
                summ += clusters_rem[entry]
        for entry in lane_maps[index]:
            if not index in ideal_ratios:
                ideal_ratios[index] = list()
            if clusters_rem[entry] > 0:
                ideal_ratios[index].append(clusters_rem[entry]/float(summ))
            else: 
                ideal_ratios[index].append(0.0)
        #Minimal number of required lanes per pool
        req_lanes[index] = summ/float(target_clusters)
        
    #Have to be rounded up, rounding down when only using duplicates makes no sense
    total_lanes = map(math.ceil, req_lanes.values())
    
    
    #Since some samples are strong and some weaksauce, 10% in ideal_ratios does not mean 10% of lane volume
    #As such, ideal_ratios need to be divided by actual_reads/expected_reads
    #Compares against sequenced lane and ignores Undetermined clusters
    #No conc diff, puts all changes into volume alteration
    for ind in xrange(1, len(lane_maps.keys())+1):
        #Bases w/o sample are not expected
        exp = 1/float(len(lane_maps[ind])-1)
        laneTypeExpr = 0
        counter = 0
        for sample in lane_maps[ind]:
            if not sample == 'Undetermined':
                laneTypeExpr += clusters_expr[sample]
        for sample in lane_maps[ind]:
            act = clusters_expr[sample]/float(laneTypeExpr)
            ideal_ratios[ind][counter] = ideal_ratios[ind][counter]*(exp/act)
            counter += 1
                   
    #normalizes
    for index in xrange(1, len(ideal_ratios.keys())+1):
        curSum = sum(ideal_ratios[index])    
        for sample in xrange(0, len(ideal_ratios[index])):
            ideal_ratios[index][sample] = (ideal_ratios[index][sample]/curSum)*100
            
    
    #Iteratively rounds to whole percent (min pipette for volume) to reach 100%
    acc_ratios = copy.deepcopy(ideal_ratios)
    for index in xrange(1, len(ideal_ratios.keys())+1):
        for sample in xrange(0, len(ideal_ratios[index])):
            acc_ratios[index][sample] = math.ceil(ideal_ratios[index][sample])
        if sum(acc_ratios[index]) == 100:
            break
        else:
            while sum(acc_ratios[index]) > 100:
                stuck = True
                for sample in xrange(1, len(ideal_ratios[index])):
                    need = ideal_ratios[index][sample]*req_lanes.values()[index-1]
                    cur = (acc_ratios[index][sample] - 1)*total_lanes[index-1]
                    if sum(acc_ratios[index]) > 100 and cur >= need:
                        acc_ratios[index][sample] -= 1
                        stuck = False
                    if sum(acc_ratios[index])== 100:
                        break
                if(stuck):
                    total_lanes[index-1] += 1
                 
                    
    # ratio * req_lanes.values() = needed
    # ratio * total_lanes = current
    # means a sample can take any whole number between the two
    # Opt for ratio * total_lanes and then go down by 1perc (one sample at a time, that dont pass ratio * req_lanes.values() = needed

    
    #Print output including duplicates
    bin = 0
    for index in xrange(1, len(lane_maps)+1):
        for dupes in xrange(1, int(total_lanes[index-1])+1):
            bin  += 1
            print ""
            print "Bin", bin
            for counter in xrange(1, len(lane_maps[index])):
                print lane_maps[index][counter], acc_ratios[index][counter], "%"
    OPT = sum(total_lanes)/sum(req_lanes.values())
    print "Min lanes:", sum(req_lanes.values()) , "Total lanes:", sum(total_lanes), "OPT:", OPT
    
    #Creates csv   
    #add timestamp and proj name
    name = "_repool.csv" 
    with open(name, 'w') as csvfile:
        writer = csv.writer(csvfile)
        #<source plate ID>,<source well>,<volume>,<destination plate ID>,<destination well>
        #writer.writerow()
    
    #CHECK FOR ODD CASES (2 PROJECTS ONE LANE); NON-DUPLICATE LANES WITH SHARED SAMPLE
    


target_clusters = 320*1000000
clusters_per_lane = 380*1000000

couch = connection()
structure = proj_struct(couch, 'P2652')
aggregator(structure, target_clusters)
pdb.set_trace()