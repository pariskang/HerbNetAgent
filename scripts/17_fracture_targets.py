#!/usr/bin/env python3
"""
Step 17: Bring bone-fracture-healing biology into the analysis. Objectively assess
how Dahuang Zhechong Wan (and rivaroxaban) intersect fracture-healing targets,
and — critically — the DIRECTION (pro- vs anti-healing) of those actions, mapped
onto the healing phases (the patient is ~6 weeks = hard-callus/remodeling phase).
"""
import requests, json, csv, time
UA={'User-Agent':'Mozilla/5.0 (research)'}
API="https://api.platform.opentargets.org/api/v4/graphql"

BONE={"MONDO_0005315":"bone fracture","MONDO_0005298":"osteoporosis","EFO_0009707":"fracture nonunion"}
def assoc(efo,size=1500):
    q="""query($id:String!,$sz:Int!){disease(efoId:$id){name associatedTargets(page:{index:0,size:$sz}){count rows{target{approvedSymbol} score}}}}"""
    r=requests.post(API,json={'query':q,'variables':{'id':efo,'sz':size}},headers=UA,timeout=60).json()
    d=r['data']['disease']
    return {row['target']['approvedSymbol']:row['score'] for row in d['associatedTargets']['rows']} if d else {}

bone_score={}
for efo,name in BONE.items():
    m=assoc(efo); 
    for g,s in m.items(): bone_score[g]=max(bone_score.get(g,0),s)
    print(f"{name}: {len(m)} targets")
print(f"union bone-healing associated genes: {len(bone_score)}")

drug=set(json.load(open("data/processed/drug_targets.json")))
riva={"F10","F2"}
inter=sorted(drug & set(g for g,s in bone_score.items() if s>=0.05), key=lambda g:-bone_score[g])
print(f"\nformula ∩ bone (score>=0.05): {len(inter)} genes")

# ---- Curated fracture-healing modules + DIRECTION of formula action ----
# phase: I=inflammation/hematoma(d0-7), A=angiogenesis(throughout), C=chondro/soft callus(wk1-3),
#        O=osteogenesis/hard callus(wk3-6+), R=remodeling(wk6+), K=coagulation scaffold(early)
# effect_on_healing of the formula's action on that target (from known pharmacology/literature):
MODULES=[
 # gene, phase, formula_action, net_effect_on_fracture_healing, note
 ("PTGS2","I/O","抑制(quercetin/apigenin)","⚠️ 不利","COX-2 是骨折愈合必需;NSAID样抑制延迟骨愈合"),
 ("PTGS1","I","抑制(oleic acid)","中性/轻不利","COX-1"),
 ("KDR","A","抑制(quercetin,VEGFR2)","⚠️ 潜在不利","血管新生对骨痂血运关键;抗VEGFR2或碍血运"),
 ("HIF1A","A","调节","🟡 双向","缺氧诱导血管新生;HIF激活利于骨愈合"),
 ("MMP9","R","抑制(isoliquiritigenin)","⚠️ 潜在不利","MMP9缺失延迟骨痂血管化/重塑"),
 ("MMP2","R","抑制(isoliquiritigenin/luteolin)","⚠️ 潜在不利","参与骨痂重塑"),
 ("MMP13","R","(未直接命中)","—","软骨/骨重塑主力胶原酶"),
 ("NFKB1","I","抑制(quercetin)","🟡 双向","过度炎症有害,但早期炎症必需"),
 ("NLRP3","I","抑制(isoliquiritigenin)","🟡 偏有利","抑制过度炎症小体或利后期愈合"),
 ("TNF","I","下调","🟡 双向","早期需要,慢性过量有害"),
 ("IL6","I","下调","🟡 双向","同上"),
 ("ESR1","O","激动样(植物雌激素)","🟢 潜在有利","雌激素信号促成骨/抗骨丢失"),
 ("CTNNB1","O","命中(Wnt/β-catenin)","🟢 潜在有利","Wnt通路促成骨"),
 ("RUNX2","O","(黄酮文献促成骨)","🟢 潜在有利","成骨主控转录因子"),
 ("BMP2","O","(间接)","🟢 有利","诱导成骨"),
 ("ALPL","O","(黄酮促ALP)","🟢 潜在有利","成骨标志"),
 ("F2","K","抑制(quercetin+水蛭素)","⚠️ 早期不利/现阶段影响小","凝血酶参与早期血肿支架;6周已过血肿期"),
 ("F10","K","利伐沙班抑制","⚠️ 早期不利/现阶段影响小","抗凝影响早期血肿;6周影响较小"),
]

rows=[]
for gene,phase,action,effect,note in MODULES:
    in_formula = gene in drug
    in_riva = gene in riva
    bs = round(bone_score.get(gene,0),3)
    rows.append({"gene":gene,"phase":phase,"in_formula":int(in_formula),"in_rivaroxaban":int(in_riva),
                 "bone_assoc_score":bs,"formula_action":action,"effect_on_healing":effect,"note":note})
with open("results/tables/fracture_targets.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["gene","phase","in_formula","in_rivaroxaban","bone_assoc_score","formula_action","effect_on_healing","note"])
    w.writeheader(); w.writerows(rows)

# tally directions among targets the formula actually hits
hit=[r for r in rows if r["in_formula"] or r["in_rivaroxaban"]]
adverse=[r for r in hit if "不利" in r["effect_on_healing"]]
favor=[r for r in hit if "有利" in r["effect_on_healing"]]
mixed=[r for r in hit if "双向" in r["effect_on_healing"]]
print("\n=== formula/rivaroxaban actions on fracture-healing targets ===")
for r in hit:
    tag="F" if r["in_formula"] else ""; tag+="R" if r["in_rivaroxaban"] else ""
    print(f"  {r['gene']:7s}[{r['phase']:3s}] {r['effect_on_healing']:14s} {r['formula_action']:22s} ({tag})")
print(f"\nADVERSE(潜在不利): {[r['gene'] for r in adverse]}")
print(f"FAVORABLE(潜在有利): {[r['gene'] for r in favor]}")
print(f"MIXED(双向): {[r['gene'] for r in mixed]}")

summary={"bone_union_genes":len(bone_score),"formula_bone_intersection":inter[:25],
         "n_intersection":len(inter),
         "adverse":[r['gene'] for r in adverse],"favorable":[r['gene'] for r in favor],
         "mixed":[r['gene'] for r in mixed]}
json.dump(summary,open("data/processed/fracture_summary.json","w"),ensure_ascii=False,indent=1)
print("\nwrote results/tables/fracture_targets.csv and fracture_summary.json")
