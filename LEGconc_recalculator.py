import pdb
import json
import couchdb
import re
from collections import defaultdict
import unicodedata

def initial_data():
    couch = couchdb.Server('http://isak:Purpleplant89@tools-dev.scilifelab.se:5984')
    db = couch['x_flowcells']
    view = db.view('names/project_ids_list')
    #Creates a proj -> FCs structure
    fc_track = defaultdict(set)
    
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
    return db, fc_track

def actual_app(db, fc_track, project):
    #Collects data in structure [sample][fc][lane] 
    sample_output = dict()
    #Reads expressed per sample
    sample_total = dict()
    #Samples for a given lane
    samples_in_lane = dict()
    
    #Create another hash sample -> FC that gives reads expressed
    for p in fc_track:
        for fc in fc_track[p]:
            samples_in_lane[fc] = dict()
            if 'illumina' in db[fc_track[p][fc]]:
                xflow_dem = db[fc_track[p][fc]]['illumina']['Demultiplex_Stats']['Barcode_lane_statistics']
                for index in xrange(0, len(xflow_dem)):
                    sample_name = xflow_dem[index]['Sample']
                    lane = xflow_dem[index]['Lane']
                    if 'Clusters' in xflow_dem[index]:
                        clusters = xflow_dem[index]['Clusters']
                    else:
                        clusters = xflow_dem[index]['PF Clusters']
                    #Difference between clusters and raw clusters?
                    if not sample_name == "Undetermined" or not sample_name == "unknown":
                        if not lane in samples_in_lane[fc]:
                            samples_in_lane[fc][lane] = 0
                        samples_in_lane[fc][lane] += 1
                    #Removes commas
                    clusters = int(re.sub(r",", "", clusters))
                    if not sample_name in sample_output.keys():
                        sample_output[sample_name] = dict()
                    if not fc in sample_output[sample_name].keys():
                        sample_output[sample_name][fc] = dict()
                    sample_output[sample_name][fc][lane] = clusters
                    #Also sets sample total
                    if not sample_name in sample_total.keys():
                        sample_total[sample_name] = 0
                    sample_total[sample_name] += clusters

                
    #Check that all samples of THE project are in pools of identical fractions. Might have to check actual sample names later; or omit this
    for sample in sample_output:
        if project in sample:
            a_fc = sample_output[sample].keys()[0]
            a_lane = sample_output[sample][a_fc].keys()[0]
            that_lane = samples_in_lane[a_fc][a_lane]
            for fc in sample_output[sample]:
                for lane in sample_output[sample][fc]:
                    this_lane = samples_in_lane[fc][lane]
                    if this_lane != that_lane:
                        print "Error! Sample appeared in differing pools"
                        pdb.set_trace()
    
    #Gives optimal sample ratio for each lane
    unexpr_clust = 0
    for sample in sample_total:
        if project in sample:
            if target_clusters > sample_total[sample]:
                unexpr_clust += (target_clusters - sample_total[sample])
    pdb.set_trace()
    
    
    #Converts n into c or v difference
    #Calculates the least amount of necessary lane duplications
    #Outputs the mapping (print format)
    #Outputs the mapping (csv format)
    
    #REMEMBER, ASSUMES LANES ARE 100% DUPLICATE

target_clusters = 320*1000000
clusters_per_lane = 380*1000000

[db, fc_track] = initial_data()
actual_app(db, fc_track, 'P2652')
