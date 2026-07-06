#!/usr/bin/env python3
"""Curate thrombosis/cardiovascular/inflammation-relevant enriched pathways."""
import csv
enr=list(csv.DictReader(open("results/tables/enrichment.csv")))
KW=["platelet","coagul","hemostasis","haemostasis","fibrin","complement","thromb","atheroscler",
    "hif","hypox","inflammat","vegf","nitric oxide","arachidon","leukocyte","wound","fluid shear",
    "angiogen","interleukin"," tnf","nf-kappa","nlrp","lipoxygen","prostagland","vascular","endothel",
    "chemokine","cytokine","integrin","oxidative","eicosan"]
def rel(t):
    tl=t.lower(); return any(k in tl for k in KW)
rows=[r for r in enr if rel(r["term"]) and float(r["adj_pval"])<0.05]
rows.sort(key=lambda r:float(r["adj_pval"]))
with open("results/tables/enrichment_thrombosis_relevant.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=enr[0].keys()); w.writeheader(); w.writerows(rows)
print(f"curated {len(rows)} thrombosis-relevant enriched terms")
