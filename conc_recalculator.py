import re
from collections import defaultdict, OrderedDict, Counter
import time
import pdb
import sys

import socket
import os
import couchdb
import glob
import re
import logging

sys.setrecursionlimit(10000) # Overrides recursion limit

def connection():
    couch = couchdb.Server('http://isak:Purpleplant89@tools-dev.scilifelab.se:5984')
    db = couch['x_flowcells']
    view = db.view('names/project_ids_list')
    pdb.set_trace()

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
    #pdb.set_trace()
    mini = 0
    for k,v in initialBin.items():
        if v[0] > mini:
            mini = v[0]
    print mini
    
    return initialBin

target_clusters = 320*1000000
clusters_per_lane = 380*1000000
connection()