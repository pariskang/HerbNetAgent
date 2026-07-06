#!/usr/bin/env python3
"""
Step 2: Map each bioactive constituent to human protein targets using ChEMBL
*measured* bioactivity data (reproducible, curated).

Method:
  compound -> PubChem InChIKey -> ChEMBL molecule (exact structure match)
           -> ChEMBL activities (Homo sapiens, pChEMBL >= 5, i.e. <=10 uM)
           -> target_chembl_id -> gene symbol (SINGLE PROTEIN targets)
pChEMBL >= 5 is the conventional 'active' cutoff; >=6 flagged as high-confidence.
"""
import requests, csv, time, json, sys, os
UA={'User-Agent':'Mozilla/5.0 (research/network-pharmacology)'}
CHEMBL="https://www.ebi.ac.uk/chembl/api/data"
PCHEMBL_MIN=5.0

def get(url, params=None, tries=4):
    for i in range(tries):
        try:
            r=requests.get(url, params=params, headers=UA, timeout=40)
            if r.status_code==200: return r.json()
            if r.status_code==404: return None
        except Exception: pass
        time.sleep(1.5*(i+1))
    return None

def inchikey_for_cid(cid):
    d=get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/InChIKey/JSON")
    try: return d["PropertyTable"]["Properties"][0]["InChIKey"]
    except Exception: return None

def chembl_id_for_inchikey(ik):
    d=get(f"{CHEMBL}/molecule.json",{'molecule_structures__standard_inchi_key':ik,'limit':1})
    try:
        ms=d.get("molecules",[])
        return ms[0]["molecule_chembl_id"] if ms else None
    except Exception: return None

def activities_targets(chembl_id):
    """Return dict target_chembl_id -> max pchembl (human activities)."""
    targets={}
    url=f"{CHEMBL}/activity.json"
    params={'molecule_chembl_id':chembl_id,'pchembl_value__gte':PCHEMBL_MIN,
            'target_organism':'Homo sapiens','limit':1000}
    d=get(url,params)
    while d:
        for a in d.get("activities",[]):
            t=a.get("target_chembl_id"); pv=a.get("pchembl_value")
            if not t or pv is None: continue
            pv=float(pv)
            if t not in targets or pv>targets[t]: targets[t]=pv
        nxt=d.get("page_meta",{}).get("next")
        if not nxt: break
        d=get("https://www.ebi.ac.uk"+nxt)
    return targets

TARGET_CACHE={}
def target_gene(tid):
    if tid in TARGET_CACHE: return TARGET_CACHE[tid]
    d=get(f"{CHEMBL}/target/{tid}.json")
    gene=None; ttype=None
    if d:
        ttype=d.get("target_type")
        if ttype=="SINGLE PROTEIN":
            for comp in d.get("target_components",[]):
                for syn in comp.get("target_component_synonyms",[]):
                    if syn.get("syn_type")=="GENE_SYMBOL":
                        gene=syn.get("component_synonym"); break
                if gene: break
    TARGET_CACHE[tid]=(gene,ttype)
    return gene,ttype

# ---- run ----
comps=list(csv.DictReader(open("data/processed/compounds.csv")))
edges=[]  # compound, gene, pchembl
for c in comps:
    cid=c["pubchem_cid"]
    if not cid: continue
    ik=inchikey_for_cid(cid); time.sleep(0.05)
    if not ik:
        print(f"  {c['compound']:26s} no InChIKey"); continue
    mid=chembl_id_for_inchikey(ik); time.sleep(0.05)
    if not mid:
        print(f"  {c['compound']:26s} no ChEMBL match"); continue
    tmap=activities_targets(mid); time.sleep(0.05)
    genes={}
    for tid,pv in tmap.items():
        g,tt=target_gene(tid)
        if g:
            if g not in genes or pv>genes[g]: genes[g]=pv
    for g,pv in genes.items():
        edges.append({"compound":c["compound"],"gene":g,"max_pchembl":round(pv,2)})
    print(f"  {c['compound']:26s} ChEMBL={mid:14s} targets={len(genes)}")

with open("data/processed/compound_targets.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["compound","gene","max_pchembl"]); w.writeheader(); w.writerows(edges)

drug_targets=sorted(set(e["gene"] for e in edges))
json.dump(drug_targets,open("data/processed/drug_targets.json","w"),indent=1)
print(f"\n[compound_targets] {len(edges)} compound-target edges; {len(drug_targets)} unique drug targets")
print("wrote data/processed/compound_targets.csv and drug_targets.json")
