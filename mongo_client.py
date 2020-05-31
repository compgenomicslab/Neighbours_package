## connection to mongo client ###

from pymongo import MongoClient

def mongo_connect():
        global coll_unigenes
        global coll_clusters
        db = None
        if not db:
                client = MongoClient()
                client = MongoClient('fat01', 27017, maxPoolSize=10)
                db = client.gmgc_unigenes
                coll_unigene = db.neighbour
                coll_cluster = db.emapper_v2
                coll_e5 = db.eggnog_v5

        return [coll_unigene,coll_cluster,coll_e5]

## changed variable of collections