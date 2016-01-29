#!/usr/bin/env python2.7

import couchdb
import re
import math
from collections import defaultdict, Counter
import unicodedata
import csv
import copy
import click

from time import time
from datetime import datetime

from genologics.config import BASEURI, USERNAME, PASSWORD
from genologics.lims import Lims
from genologics.entities import Process

#Assumes ind. sample conc measurements have failed. As such it relies on changing relative volume on already normalized samples and structure
#Structure are retained as conc measurements failure means there's no way to know conc. delta between samples from seperate poolss
def connection():
    couch = couchdb.Server('http://isak:Purpleplant89@tools-dev.scilifelab.se:5984')
    return couch

#Fetches the structure of a project
def proj_struct(couch, project, target_clusters):
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

def aggregator(struct, target_clusters, project, destid):
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
    
    #Crude way to check that no samples are in different TYPES of lanes
    tempList = list()
    for k, v in lane_maps.items():
        for index in xrange(1,len(v)):
            if not v[index] == 'Undetermined':
                tempList.append(v[index])
    counter = Counter(tempList)
    for values in counter.itervalues():
        if values > 1: 
            raise Exception('Error: This app does NOT handle situations where a sample is present in lanes/well with differing structure!')

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
    # ideal_ratio * req_lanes.values() = needed
    # acc_ratio * total_lanes = current
    # means a sample can take any whole number between the two
    # Go down from current to needed 1% at a time until one hits sum(lane) == 100
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
     
                    
    #Gathers the container id and well name for all samples in project
    #Cred to Denis for providing a base epp
    location = dict()
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    allProjects = lims.get_projects()
    for proj in allProjects:
        if proj.id == project:
            projName = proj.name 
            break

    #All normalization processes for project
    norms=['Library Normalization (MiSeq) 4.0', 'Library Normalization (Illumina SBS) 4.0','Library Normalization (HiSeq X) 1.0']
    pros=lims.get_processes(type=norms, projectname=projName)
    #For all processes
    for p in pros:
        #For all artifacts in process
        for o in p.all_outputs():
            #If artifact is analyte type and has project name in sample
            if o.type=="Analyte" and project in o.name:
                location[o.name.split()[0]] = list()
                location[o.name.split()[0]].append(o.location[0].id)
                location[o.name.split()[0]].append(o.location[1])
    
    
    #Print stats including duplicates
    timestamp = datetime.fromtimestamp(time()).strftime('%Y-%m-%d_%H:%M')
    sumName = projName,  "_summary_", timestamp,".txt"
    sumName = ''.join(sumName)
    with open(sumName, "w") as summary:
        OPT = sum(total_lanes)/sum(req_lanes.values())
        output = "Ideal lanes (same schema): ", str(sum(req_lanes.values())) , ", Total lanes: ", str(sum(total_lanes)), ", OPT: ", str(round(OPT,3)),'\n'
        output = ''.join(output)
        summary.write( output )
        output = "Unique pools: ", str(len(total_lanes)), ", Average pool duplication: ", str(sum(total_lanes)/float(len(total_lanes))) ,'\n'
        output = ''.join(output)
        summary.write( output )
        
        bin = 0
        for index in xrange(1, len(lane_maps)+1):
            bin  += 1
            summary.write('\n')
            output = "Wells ", str(bin) , '-' , str(bin+int(total_lanes[index-1])-1),':','\n'
            output = ''.join(output)
            summary.write( output )
            bin += int(total_lanes[index-1]-1)
            for counter in xrange(1, len(lane_maps[index])):
                output = str(lane_maps[index][counter]),' ', str(acc_ratios[index][counter]), "%",'\n'
                output = ''.join(output)
                summary.write( output )

    
    #Creates csv   
    
    name = projName,"_repool_",timestamp,".csv"
    name = ''.join(name)
    wells = ['Empty','A','B','C','D','E','F','G','H']
    wellIndex = [1, 1]
    destNo = 0
    
    with open(name, 'w') as csvfile:
        writer = csv.writer(csvfile)
        for index in xrange(1, len(lane_maps)+1):
            for dupes in xrange(1, int(total_lanes[index-1])+1):
                for counter in xrange(1, len(lane_maps[index])):
                    #<source plate ID>,<source well>,<volume>,<destination plate ID>,<destination well>
                    #Destination well 200 microL, minimum pipette 2 microL; acc_ratios multiplied by 2.
                    sample = lane_maps[index][counter]
                    position = wells[wellIndex[1]],':',str(wellIndex[0])
                    position = ''.join(position)
                    output = location[sample][0],location[sample][1],str(int(acc_ratios[index][counter]*2)),str(destid[destNo]),position
                    if not acc_ratios[index][counter] == 0:
                        writer.writerow(output)
                
                #Increment wellsindex
                if not acc_ratios[index][counter] == 0:
                    if not wellIndex[1] >= 8:
                        wellIndex[1] += 1
                    else:
                        wellIndex[1] = 1
                        if not wellIndex[0] >= 8:
                            wellIndex[0] += 1
                        else:
                            wellIndex[0] = 1
                            destNo += 1
                            try:
                                destid[destNo]
                            except IndexError:
                                print "Critical error; not enough destination plates provided"
                            
@click.command()
@click.option('--project_id', required=True,help='REQUIRED: ID of project to repool. Examples:P2652, P1312 etc.')
@click.option('--dest_plate_list', default=['dp_1','dp_2','dp_3','dp_4','dp_5'], 
              help='List of destination plates for the robot\'s csv file. Include too many rather than too few; excess will be unused Default:[dp_1,dp_2,dp_3,dp_4,dp_5]') 
@click.option('--target_clusters', default=320*1000000, help='Threshold of clusters per sample. \nDefault:320*1000000')
@click.option('--clusters_per_lane', default=380*1000000, help='Expected clusters generated by a single lane/well. \nDefault:380*1000000')            

def main(target_clusters, clusters_per_lane, project_id, dest_plate_list):
    """Application that calculates samples under threshold for a project, then calculate the optimal composition for reaching the threshold
    without altering concentrations nor the structure of the pools. Outputs both a summary as well as a functional csv file."""
    couch = connection()
    structure = proj_struct(couch, project_id, target_clusters)
    aggregator(structure, target_clusters, project_id, dest_plate_list)

if __name__ == '__main__':
    main()