import json
import os

from elasticsearch import Elasticsearch

from cdi.dbqueries import dbqueries

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


cancertypes = ["new_anti_cancer",
               "Bladder",
               "Bone",
               "Breast", "Colon", "Kidney", "Liver",
               "Luk", "Lung", "NonHodgkinLymph", "Ovarian", "Pancreatic",
               "PlasmaCell", "Rectal",
               "Thyroid"
               ]

def test_readaminoacidseqs():
    targets = json.load(open("targets.json"))
    for cancertype in cancertypes:
        for drug in targets[cancertype]:
            print("Drug name: %s, id: %s" % (drug["name"], drug["drugbankid"]))
            for target in drug["targets"]:
                print("Target gene name: %s" % target["gene"])
                print(target["amino-acid-sequence"][:200])


def read_targets(folder="approved_cancer_drugs"):
    r = dict()  # final obj to store drug target data for each cancer type
    missingtargetgenes = set()
    ttd = folder[:4] == "ttd_"
    for cancer in cancertypes:
        print('\n\n**** %s *******\n' % cancer)
        infile = f"../data/{folder}/"
        if ttd:
            infile += f"DT_{cancer}.csv"
        else:
            infile += f"{cancer}_approved_drugs.txt"
        if not os.path.exists(infile):
            r[cancer] = []
            continue
        drugnames = open(infile).read().splitlines(keepends=False)
        if ttd:
            drugnames = drugnames[1:]
        drugs = list()
        for ilac in drugnames:
            if ttd:
                genes = {g for g in ilac.split(",")[-1].split("; ")
                         if len(g) > 0}
                ilac = ilac.split(",")[0]
                print("TTD genes:", genes)

            records, nrecords = essearch(ilac, size=1)
            if nrecords == 0:
                print("*********** no Drugbank records matched ************\n")
                continue
            for rec in records:
                _id = rec['_id']
                rec = rec['_source']
                if rec['synonyms'] is None:
                    s = []
                elif not isinstance(rec['synonyms']['synonym'], list):
                    s = [rec['synonyms']['synonym']['#text']]
                else:
                    s = [s["#text"]
                         for s in rec['synonyms']['synonym']
                         if s['language'] == 'english']

                if 'products' in rec:
                    products = list({pr['name'] for pr in rec['products']})
                else:
                    products = []
                trgts = list()
                if 'targets' in rec:
                    for trgt in rec['targets']:
                        if 'polypeptide' in trgt:
                            pp = trgt['polypeptide']
                            if isinstance(pp, list):
                                for pp_ in pp:
                                    gene = pp_['gene-name']
                                    if ttd and gene in genes:
                                        genes.remove(gene)
                                    trgts.append({
                                        'gene': gene,
                                        'id': pp_['id'],
                                        "amino-acid-sequence":
                                            pp_["amino-acid-sequence"]["#text"],
                                        # "gene-sequence":
                                        #     pp_["gene-sequence"]["#text"]
                                    })
                            else:
                                if "#text" not in pp["gene-sequence"]:
                                    # 2 polypeptides match here:
                                    # A9X444, P54756
                                    print("polypeptide %s's gene sequence"
                                          " is missing" % pp['id'])
                                    # geneseq = ""
                                # else:
                                # geneseq = pp["gene-sequence"]["#text"]
                                gene = pp['gene-name']
                                if ttd and gene in genes:
                                    genes.remove(gene)
                                trgts.append({
                                    'gene': gene,
                                    'id': pp['id'],
                                    "amino-acid-sequence":
                                        pp["amino-acid-sequence"]["#text"],
                                    # "gene-sequence": geneseq
                                })
                drugs.append({
                    "queryname": ilac,
                    "name": rec['name'],
                    "drugbankid": _id,
                    "numberofrecordsmatched": nrecords,
                    "synonyms": s,
                    "products": products,
                    "targets": trgts
                })
                if ttd and len(genes) > 0:
                    print("*** Target genes not in drugbank record: ***")
                    print(genes)
                    drugs[-1]["missingtargetgenes"] = list(genes)
                    missingtargetgenes.update(genes)
        r[cancer] = drugs
    json.dump(r, open(f"{folder}.json", "w"), indent=4)
    return r, missingtargetgenes

def read_missinggeneseqs(folder="ttd_approved_cancer_drugs"):
    r = dict()
    infile = f"../data/{folder}/aaseqs.fasta"
    from Bio import SeqIO
    for record in SeqIO.parse(infile, "fasta"):
        _id = record.id
        if _id in ['Candi-TMP1', 'HSV-UL30']:
            _id = _id.replace('-', ' ')
        r[_id] = ">%s\n%s" % (_id, str(record.seq))
    return r


def test_read_missinggeneseqs():
    read_missinggeneseqs()

def test_merge_targets():
    # Merge NIC/TTD drugs/targets

    nic = read_targets("approved_cancer_drugs")[0]
    ttd, missingtargetgenes = read_targets("ttd_approved_cancer_drugs")

    r = dbqueries.get_target_genes_aaseqs(list(missingtargetgenes), 18000)
    c = read_missinggeneseqs()
    for j in missingtargetgenes - {gene for gene in r}:
        assert j in c
        r[j] = c[j]
    mergeddts = dict()
    for cancertype in cancertypes:
        nicdrugs = {drug["name"]:drug for drug in nic[cancertype]}

        for drug in ttd[cancertype]:
            print([i['gene'] for i in drug["targets"]])
            if "missingtargetgenes" in drug:
                for target in drug["missingtargetgenes"]:
                    drug["targets"].append(
                        {
                            "gene": target,
                            "amino-acid-sequence": r[target]
                        })
            if drug["name"] not in nicdrugs:
                print("*** Drug %s not in NIC list ***" % drug["name"])
                print([i['gene'] for i in drug["targets"]])
                nicdrugs[drug["name"]] = drug
            else:
                # Merge targets
                nicdrug = nicdrugs[drug["name"]]
                nictargets = {i['gene'] for i in nicdrug["targets"]}
                for target_ in drug["targets"]:
                    if target_["gene"] not in nictargets:
                        n = len(nicdrug["targets"])
                        nicdrug["targets"].append(target_)
                        assert len(nicdrug["targets"]) == n+1
                nicdrugs[drug["name"]] = nicdrug

        # Save nicdrugs to mergeddts
        mergeddts[cancertype] = []
        for drug in nicdrugs.values():
            mergeddts[cancertype].append(drug)

        json.dump(mergeddts, open("targets_.json", "w"), indent=4)


if __name__ == "__main__":
    test_merge_targets()