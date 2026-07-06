#!/usr/bin/env python3
"""
Step 14: Mechanism-complementarity ("增效"/synergy rationale) analysis of
rivaroxaban + Dahuang Zhechong Wan across the pathophysiological axes of venous
thrombosis. Rivaroxaban = pure FXa (F10) inhibitor (ChEMBL mechanism confirmed);
the formula adds platelet, inflammation, resolution and endothelial coverage.

We classify each agent's targets into the 6 axes of venous-thrombosis biology
and quantify single-agent vs combination coverage (this supports COMPLEMENTARY
multi-axis coverage — a rationale for synergy, NOT proof of synergy).
"""
import json, csv
from collections import defaultdict

drug=set(json.load(open("data/processed/drug_targets.json")))        # 208 formula targets
core=set(json.load(open("data/processed/hub_genes.json")))           # top hubs
inter=[r["gene"] for r in csv.DictReader(open("results/tables/intersection_targets.csv"))]
riva={"F10","F2"}   # rivaroxaban: potent FXa + weak prothrombin (ChEMBL)
riva_potent={"F10"}

# --- 6 pathophysiological axes of venous thrombosis (Virchow + resolution) ---
AXES={
 "凝血级联 Coagulation":      ["F2","F10","F5","F7","F3","F9","SERPINC1","PROC","PROS1","THBD","TFPI"],
 "纤溶 Fibrinolysis":         ["PLG","PLAT","PLAU","SERPINE1","SERPINB2"],
 "血小板活化 Platelet":        ["PTGS1","PTGS2","ALOX12","ALOX15","TBXA2R","P2RY12","P2RY1","ITGB3","GP6","GP1BA","GRK6"],
 "血栓炎症 Thrombo-inflammation":["NFKB1","NLRP3","JAK2","STAT6","IL2","IL6","IL1B","TNF","PIK3R1","TLR4","CXCL8"],
 "血栓消散/重塑 Resolution":    ["MMP2","MMP9","MMP12","MMP14","TIMP1","PLAU"],
 "内皮/血管 Endothelium":      ["KDR","VEGFA","NOS3","HIF1A","VWF","SELE","ICAM1","VCAM1"],
}

def covered(genes, axis_list):
    return sorted(set(genes) & set(axis_list))

rows=[]
print(f"{'Axis':30s} {'Rivaroxaban':16s} {'Formula':28s} Combination")
combo_axes=0; riva_axes=0; formula_axes=0
for axis,glist in AXES.items():
    r=covered(riva, glist)
    f=covered(drug, glist)
    both=sorted(set(r)|set(f))
    riva_axes+= 1 if r else 0
    formula_axes+= 1 if f else 0
    combo_axes+= 1 if both else 0
    rows.append({"axis":axis,"rivaroxaban":";".join(r),"formula":";".join(f),
                 "combination":";".join(both),
                 "riva_n":len(r),"formula_n":len(f),"combo_n":len(both)})
    print(f"{axis:30s} {(','.join(r) or '—'):16s} {(','.join(f) or '—'):28s} {','.join(both) or '—'}")

n=len(AXES)
print(f"\nAxis coverage: rivaroxaban {riva_axes}/{n} | formula {formula_axes}/{n} | combination {combo_axes}/{n}")

# complementarity: overlap of the two target sets (low overlap = high complementarity)
overlap=riva & drug
jac=len(overlap)/len(riva|drug)
print(f"Target overlap rivaroxaban∩formula: {sorted(overlap)}  (Jaccard={jac:.3f} -> highly complementary)")

with open("results/tables/mechanism_coverage.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=["axis","riva_n","formula_n","combo_n","rivaroxaban","formula","combination"])
    w.writeheader(); w.writerows(rows)

summary={"axes":n,"riva_axes":riva_axes,"formula_axes":formula_axes,"combo_axes":combo_axes,
         "overlap":sorted(overlap),"jaccard":round(jac,3),
         "riva_only_axes":[a["axis"] for a in rows if a["riva_n"] and not a["formula_n"]],
         "formula_adds_axes":[a["axis"] for a in rows if a["formula_n"] and not a["riva_n"]]}
json.dump(summary,open("data/processed/synergy_summary.json","w"),ensure_ascii=False,indent=1)
print("\nFormula adds these axes rivaroxaban lacks:", summary["formula_adds_axes"])
print("wrote results/tables/mechanism_coverage.csv and synergy_summary.json")
