from elasticsearch import Elasticsearch
import pandas as pd
import numpy as np

es = Elasticsearch(['http://esnode-khadija.cbrc.kaust.edu.sa:9400'])


def essearch(qterm, index='drugbank-5.1.8', **kwargs):
    qc = {
        "size": 1,
        "query": {
            "query_string": {
                "query": qterm
            }}}
    cr = es.search(index=index, body=qc, **kwargs)
    return cr['hits']['hits'], cr['hits']['total']


def test_synonyms():
    for cancer in ["Bladder","Bone", "Breast", "Colon", "Kidney", "Liver",
                   "Luk", "Lung", "NonHodgkinLymph", "Ovarian", "Pancreatic",
                   "PlasmaCell", "Rectal", "Thyroid"]:
        dr_syns = []
        dr_product = []
        header = ('DrugBank ID','NCI Drug Name','DrugBank Drug Name','All Synonyms')
        dr_syns.append(header)
        print('\n\n**** %s *******\n' % cancer)
        file = open(f"Datasets/Approved_cancer_drugs/{cancer}_"
                    f"approved_drugs.txt").read()
        for ilac in file.splitlines(keepends=False):
            records, nrecords = essearch(ilac, size=1)
            print(f"{ilac};  {nrecords} DrugBank records matched")
            if nrecords == 0:
                print("*********** no drugbank records matched ************\n")
                continue
            for rec in records:
                _id = rec['_id']
                rec = rec['_source']
                if rec['synonyms'] is None:
                    s = ""
                elif not isinstance(rec['synonyms']['synonym'], list):
                    s = rec['synonyms']['synonym']['#text']
                
                else:
                    s = (",".join(s["#text"] for s in rec['synonyms']['synonym'] if s['language']=='english'))
                
                print("%s, %s, %s" % (_id, rec['name'], s))
                DR = _id, ilac, rec['name'], s
                dr_syns.append(DR)
                
                # if 'products' in rec:
                #     products = (",".join({pr['name'] for pr in rec['Brand Names']}))

                # else:
                #     products = "no product name listed"
                #dr_product.append(products)
        dr_syns_df = pd.DataFrame(dr_syns[1:], columns=dr_syns[0])
        print("\n****^DRUG INFO^****")
        print(dr_syns_df)
        file_name = 'Datasets/cancer_dr_syn/'+str(cancer)+'_DB_syn.csv'
        dr_syns_df.to_csv(file_name, header=True,  index_label=None)

        print()

if __name__ == "__main__":
    test_synonyms()
