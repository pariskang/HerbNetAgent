#!/usr/bin/env python3
"""Step 6: publication-style figures for the network-pharmacology analysis."""
import csv, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import networkx as nx
import numpy as np
os.makedirs("results/figures",exist_ok=True)
plt.rcParams.update({'font.size':10})
C={'herb':'#2E7D5B','cpd':'#E08A3C','tgt':'#3B6FA0','hub':'#C0392B','dis':'#7D5BA6'}

# ---------- Fig 1: Venn (drug vs disease) ----------
drug=set(json.load(open("data/processed/drug_targets.json")))
drows=list(csv.DictReader(open("data/processed/disease_targets.csv")))
dscore={r["gene"]:float(r["max_assoc_score"]) for r in drows}
disease=set(g for g,s in dscore.items() if s>=0.1)
only_d=len(drug-disease); only_dis=len(disease-drug); both=len(drug&disease)
fig,ax=plt.subplots(figsize=(6,4.2))
ax.add_patch(Circle((0.38,0.5),0.30,alpha=0.45,color=C['cpd']))
ax.add_patch(Circle((0.62,0.5),0.30,alpha=0.45,color=C['dis']))
ax.text(0.24,0.5,str(only_d),ha='center',va='center',fontsize=14,weight='bold')
ax.text(0.76,0.5,str(only_dis),ha='center',va='center',fontsize=14,weight='bold')
ax.text(0.50,0.5,str(both),ha='center',va='center',fontsize=16,weight='bold',color='white')
ax.text(0.30,0.86,f"Dahuang Zhechong Wan\ntargets (ChEMBL)\nn={len(drug)}",ha='center',fontsize=9)
ax.text(0.72,0.86,f"Venous thrombosis\ntargets (Open Targets)\nn={len(disease)}",ha='center',fontsize=9)
ax.set_xlim(0,1);ax.set_ylim(0.1,1);ax.axis('off')
ax.set_title("Shared targets: formula ∩ venous thrombosis",weight='bold')
plt.tight_layout();plt.savefig("results/figures/fig1_venn.png",dpi=150);plt.close()

# ---------- Fig 2: PPI network, hub-highlighted ----------
edges=list(csv.DictReader(open("results/tables/ppi_edges.csv")))
tg=list(csv.DictReader(open("results/tables/intersection_targets.csv")))
deg={r["gene"]:int(r["degree"]) for r in tg}
hub=set(json.load(open("data/processed/hub_genes.json")))
G=nx.Graph()
for r in tg: G.add_node(r["gene"])
for e in edges: G.add_edge(e["source"],e["target"],weight=float(e["string_score"]))
G.remove_nodes_from([n for n in list(G.nodes()) if G.degree(n)==0])
# disease-relevant core targets (score>=0.1) to highlight
dcore=set(json.load(open("data/processed/drug_targets.json"))) & disease
if G.number_of_nodes()>0:
    pos=nx.spring_layout(G,k=1.1,seed=42,iterations=140)
    fig,ax=plt.subplots(figsize=(13,10))
    sizes=[30+deg.get(n,0)*6 for n in G.nodes()]
    def ncol(n):
        if n in dcore: return C['hub']      # red = disease-relevant core target
        if n in hub:   return C['cpd']       # orange = topology hub only
        return C['tgt']                        # blue = other
    cols=[ncol(n) for n in G.nodes()]
    nx.draw_networkx_edges(G,pos,alpha=0.08,width=0.4,ax=ax)
    nx.draw_networkx_nodes(G,pos,node_size=sizes,node_color=cols,alpha=0.92,linewidths=0.4,edgecolors='white',ax=ax)
    lbl={n:n for n in G.nodes() if n in dcore or deg.get(n,0)>=np.percentile(list(deg.values()),88)}
    nx.draw_networkx_labels(G,pos,labels=lbl,font_size=8,font_weight='bold',ax=ax)
    from matplotlib.lines import Line2D
    leg=[Line2D([0],[0],marker='o',color='w',markerfacecolor=C['hub'],markersize=11,label='disease-relevant core target (assoc≥0.1)'),
         Line2D([0],[0],marker='o',color='w',markerfacecolor=C['cpd'],markersize=11,label='topology hub (generic)'),
         Line2D([0],[0],marker='o',color='w',markerfacecolor=C['tgt'],markersize=9,label='other shared target')]
    ax.legend(handles=leg,loc='lower left',fontsize=9,frameon=True)
    ax.set_title(f"PPI network of formula × venous-thrombosis targets (STRING conf≥0.4)\n{G.number_of_nodes()} nodes, {G.number_of_edges()} edges — node size ∝ degree",weight='bold')
    ax.axis('off');plt.tight_layout();plt.savefig("results/figures/fig2_ppi_network.png",dpi=150);plt.close()

# ---------- Fig 3: hub gene centrality bar ----------
top=sorted(tg,key=lambda r:-int(r["degree"]))[:20]
fig,ax=plt.subplots(figsize=(7,6))
y=range(len(top))
ax.barh([t["gene"] for t in top][::-1],[int(t["degree"]) for t in top][::-1],color=C['hub'],alpha=0.85)
ax.set_xlabel("Degree (PPI connectivity)");ax.set_title("Top-20 hub targets by network degree",weight='bold')
plt.tight_layout();plt.savefig("results/figures/fig3_hubs.png",dpi=150);plt.close()

# ---------- Fig 4: thrombosis-relevant enrichment ----------
cur=list(csv.DictReader(open("results/tables/enrichment_thrombosis_relevant.csv")))
cur.sort(key=lambda r:float(r["adj_pval"]))
cur=cur[:16]
fig,ax=plt.subplots(figsize=(9.5,6.5))
labels=[r["term"].split(" R-HSA")[0].split(" (GO:")[0][:46]+f"  [{r['library']}]" for r in cur][::-1]
vals=[-np.log10(float(r["adj_pval"])) for r in cur][::-1]
palette={'KEGG':'#2E7D5B','Reactome':'#7D5BA6','GO_BP':'#3B6FA0','GO_MF':'#E08A3C','GO_CC':'#C0392B'}
cols=[palette.get(r["library"],'#555') for r in cur][::-1]
ax.barh(labels,vals,color=cols,alpha=0.88)
ax.set_xlabel("-log10(adjusted p)")
ax.set_title("Thrombosis / cardiovascular / inflammation pathways\nenriched among formula × venous-thrombosis targets",weight='bold',fontsize=11)
plt.tight_layout();plt.savefig("results/figures/fig4_enrichment.png",dpi=150);plt.close()

# ---------- Fig 5: targets per herb ----------
ct=list(csv.DictReader(open("data/processed/compound_targets.csv")))
comp2herb={}
for r in csv.DictReader(open("data/processed/compounds.csv")):
    comp2herb[r["compound"]]=r["herbs_zh"].split(";")
herb_tgts={}
inter=drug&disease
for e in ct:
    if e["gene"] in inter:
        for h in comp2herb.get(e["compound"],[]):
            herb_tgts.setdefault(h,set()).add(e["gene"])
PINYIN={"大黄":"Dahuang","黄芩":"Huangqin","甘草":"Gancao","桃仁":"Taoren","苦杏仁":"Kuxingren",
        "白芍":"Baishao","地黄":"Dihuang","干漆":"Ganqi","水蛭":"Shuizhi","土鳖虫":"Tubiechong",
        "虻虫":"Mengchong","蛴螬":"Qicao"}
items=sorted(herb_tgts.items(),key=lambda x:-len(x[1]))
fig,ax=plt.subplots(figsize=(8,5))
ax.bar([PINYIN.get(h,h) for h,_ in items],[len(s) for _,s in items],color=C['herb'],alpha=0.85)
ax.set_ylabel("# core mechanism targets (of 16)")
ax.set_title("Herb contribution to venous-thrombosis mechanism targets",weight='bold')
plt.xticks(rotation=40,ha='right');plt.tight_layout();plt.savefig("results/figures/fig5_herb_contribution.png",dpi=150);plt.close()
print("figures written to results/figures/")
