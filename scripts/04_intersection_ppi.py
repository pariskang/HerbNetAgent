#!/usr/bin/env python3
"""
Step 4: Intersect formula targets (ChEMBL, measured) with venous-thrombosis
disease targets (Open Targets); build the STRING PPI network on the intersection;
rank hub genes by degree & betweenness centrality.
"""
import csv, json, requests, os
import networkx as nx
UA={'User-Agent':'Mozilla/5.0 (research/network-pharmacology)'}

DIS_SCORE_MIN=0.1  # Open Targets meaningful-association cutoff for primary analysis

drug=set(json.load(open("data/processed/drug_targets.json")))
disease_rows=list(csv.DictReader(open("data/processed/disease_targets.csv")))
dscore={r["gene"]:float(r["max_assoc_score"]) for r in disease_rows}
dtop={r["gene"]:r["top_disease"] for r in disease_rows}
disease=set(g for g,s in dscore.items() if s>=DIS_SCORE_MIN)

inter=sorted(drug & disease)
inter_all=sorted(drug & set(dscore))  # any disease evidence
print(f"drug targets: {len(drug)} | disease targets(score>={DIS_SCORE_MIN}): {len(disease)} | any-evidence: {len(dscore)}")
print(f"intersection (score>={DIS_SCORE_MIN}): {len(inter)} genes")
print(f"intersection (any evidence): {len(inter_all)} genes")

# per-compound pChEMBL to annotate intersection
ct=list(csv.DictReader(open("data/processed/compound_targets.csv")))
gene_best_pchembl={}
gene_compounds={}
for e in ct:
    g=e["gene"]; pv=float(e["max_pchembl"])
    gene_best_pchembl[g]=max(gene_best_pchembl.get(g,0),pv)
    gene_compounds.setdefault(g,set()).add(e["compound"])

# STRING PPI network on intersection (use any-evidence set to keep network connected, tag score)
genes=inter_all
def string_network(genes, required_score=400):
    r=requests.post("https://string-db.org/api/tsv/network",
        data={'identifiers':'\r'.join(genes),'species':9606,'required_score':required_score},
        headers=UA,timeout=90)
    edges=[]
    for i,line in enumerate(r.text.strip().splitlines()):
        if i==0: continue
        p=line.split('\t')
        if len(p)>=6: edges.append((p[2],p[3],float(p[5])))  # preferredName_A, B, score
    return edges

edges=string_network(genes, required_score=400)
G=nx.Graph()
G.add_nodes_from(genes)
for a,b,w in edges: G.add_edge(a,b,weight=w)
print(f"STRING PPI: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges (conf>=0.4)")

deg=dict(G.degree())
btw=nx.betweenness_centrality(G)
clo=nx.closeness_centrality(G)

rows=[]
for g in genes:
    rows.append({
        "gene":g,
        "disease_score":round(dscore.get(g,0),4),
        "top_disease":dtop.get(g,""),
        "drug_max_pchembl":round(gene_best_pchembl.get(g,0),2),
        "n_compounds":len(gene_compounds.get(g,[])),
        "degree":deg.get(g,0),
        "betweenness":round(btw.get(g,0),4),
        "closeness":round(clo.get(g,0),4),
        "meaningful_assoc":int(dscore.get(g,0)>=DIS_SCORE_MIN),
    })
rows.sort(key=lambda x:(-x["degree"],-x["betweenness"]))
with open("results/tables/intersection_targets.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

# hub genes = top by degree within the largest connected component
if G.number_of_edges()>0:
    lcc=max(nx.connected_components(G),key=len)
    print(f"largest connected component: {len(lcc)} genes")
hub=rows[:20]
print("\nTop hub genes (degree | betweenness | disease_score | pChEMBL):")
for h in hub[:20]:
    print(f"  {h['gene']:10s} deg={h['degree']:3d}  btw={h['betweenness']:.3f}  "
          f"dis={h['disease_score']:.3f}  pIC50={h['drug_max_pchembl']:.1f}  n_cpd={h['n_compounds']}")

json.dump([h['gene'] for h in hub], open("data/processed/hub_genes.json","w"))
json.dump(genes, open("data/processed/network_genes.json","w"))
# save edges for figure
with open("results/tables/ppi_edges.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["source","target","string_score"])
    for a,b,wt in edges: w.writerow([a,b,round(wt,3)])
print("\nwrote results/tables/intersection_targets.csv, ppi_edges.csv, hub_genes.json")
