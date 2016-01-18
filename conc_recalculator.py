import pdb
import json
import couchdb
from collections import defaultdict
import unicodedata

def connection():
    couch = couchdb.Server('http://isak:Purpleplant89@tools-dev.scilifelab.se:5984')
    db = couch['x_flowcells']
    view = db.view('names/project_ids_list')

    fc_track = defaultdict(set)
    
    #Creates a proj -> FCs structure
    for rec in view.rows:
        fc = ''.join(rec.key)
        fc = unicodedata.normalize('NFKD', fc).encode('ascii','ignore')
        proj = ''.join(rec.value)
        proj = unicodedata.normalize('NFKD', proj).encode('ascii','ignore')
        id = ''.join(rec.id)
        id = unicodedata.normalize('NFKD', id).encode('ascii','ignore')
        if fc_track[proj] == set([]):
            fc_track[proj] = dict()
        fc_track[proj][fc] = id
    
    #Creates FC -> sample_info structure
    map_fun = """function(doc) {
      emit(doc['name'], doc['illumina']);
    }"""
    
    for red in (db.query(map_fun)).rows:
        pdb.set_trace()
        
    pdb.set_trace()

target_clusters = 320*1000000
clusters_per_lane = 380*1000000

view = connection()
