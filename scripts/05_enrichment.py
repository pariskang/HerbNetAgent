#!/usr/bin/env python3
"""
Step 5: GO (BP/MF/CC) + KEGG pathway + Reactome enrichment of the mechanism
target set via the Enrichr API (Fisher exact + Benjamini-Hochberg adjusted p).
"""
import requests, json, csv, time, os
UA={'User-Agent':'Mozilla/5.0 (research/network-pharmacology)'}
BASE="https://maayanlab.cloud/Enrichr"

genes=json.load(open("data/processed/network_genes.json"))
print(f"enriching {len(genes)} mechanism target genes")

def add_list(genes):
    payload={'list':(None,'\n'.join(genes)),'description':(None,'DHZCW_venous_thrombosis')}
    r=requests.post(BASE+"/addList",files=payload,headers=UA,timeout=60)
    return r.json()['userListId']

def enrich(uid, lib):
    for i in range(4):
        try:
            r=requests.get(BASE+"/enrich",params={'userListId':uid,'backgroundType':lib},headers=UA,timeout=60)
            if r.status_code==200: return r.json().get(lib,[])
        except Exception: pass
        time.sleep(2)
    return []

uid=add_list(genes); time.sleep(1)
LIBS={
 "GO_Biological_Process_2023":"GO_BP",
 "GO_Molecular_Function_2023":"GO_MF",
 "GO_Cellular_Component_2023":"GO_CC",
 "KEGG_2021_Human":"KEGG",
 "Reactome_2022":"Reactome",
}
allrows=[]
for lib,tag in LIBS.items():
    res=enrich(uid,lib); time.sleep(0.8)
    # Enrichr row: [rank, term, pval, zscore, combined, genes, adj_pval, ...]
    sig=[r for r in res if r[6]<0.05]
    print(f"{tag:10s}: {len(res)} terms, {len(sig)} with adj-p<0.05")
    for r in res[:60]:
        allrows.append({"library":tag,"term":r[1],"pval":r[2],"adj_pval":r[6],
                        "combined_score":round(r[4],2),"n_genes":len(r[5]),
                        "genes":";".join(r[5])})

with open("results/tables/enrichment.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["library","term","pval","adj_pval","combined_score","n_genes","genes"])
    w.writeheader(); w.writerows(allrows)

print("\n== Top KEGG pathways (adj-p<0.05) ==")
for r in [x for x in allrows if x["library"]=="KEGG" and x["adj_pval"]<0.05][:15]:
    print(f"  adjP={r['adj_pval']:.2e}  n={r['n_genes']:2d}  {r['term']}")
print("\n== Top GO-BP (adj-p<0.05) ==")
for r in [x for x in allrows if x["library"]=="GO_BP" and x["adj_pval"]<0.05][:15]:
    print(f"  adjP={r['adj_pval']:.2e}  n={r['n_genes']:2d}  {r['term']}")
print("\nwrote results/tables/enrichment.csv")
