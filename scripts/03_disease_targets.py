#!/usr/bin/env python3
"""
Step 3: Collect disease-associated targets for venous thrombosis from the
Open Targets Platform (aggregates genetics, known drugs, literature, RNA, etc.).

Lower-limb intermuscular (muscular calf vein) thrombosis is a distal deep-vein
thrombosis / venous thromboembolism. We union the core venous-thrombosis disease
nodes and keep, per gene, the maximum overall association score + contributing disease.
"""
import requests, csv, json, time, os
UA={'User-Agent':'Mozilla/5.0 (research/network-pharmacology)'}
API="https://api.platform.opentargets.org/api/v4/graphql"

DISEASES={
 "EFO_0003907":"deep vein thrombosis",
 "MONDO_0005399":"venous thromboembolism",
 "MONDO_0002305":"thrombophilia",
 "MONDO_0000831":"thrombotic disease",
 "MONDO_0005279":"pulmonary embolism",
}

def assoc_targets(efo, page_size=1000):
    q="""query($id:String!,$idx:Int!,$sz:Int!){
      disease(efoId:$id){ name
        associatedTargets(page:{index:$idx,size:$sz}){
          count rows{ target{ id approvedSymbol } score } } } }"""
    out=[]; idx=0
    while True:
        r=requests.post(API,json={'query':q,'variables':{'id':efo,'idx':idx,'sz':page_size}},
                        headers=UA,timeout=60)
        d=r.json()['data']['disease']
        if not d: break
        rows=d['associatedTargets']['rows']; cnt=d['associatedTargets']['count']
        for row in rows:
            sym=row['target']['approvedSymbol']; sc=row['score']
            out.append((sym,sc))
        idx+=1
        if idx*page_size>=cnt or not rows: break
        time.sleep(0.3)
    return out

best={}  # gene -> (max_score, disease)
for efo,name in DISEASES.items():
    rows=assoc_targets(efo)
    print(f"{name:24s} ({efo}): {len(rows)} targets")
    for sym,sc in rows:
        if sym not in best or sc>best[sym][0]:
            best[sym]=(sc,name)

with open("data/processed/disease_targets.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["gene","max_assoc_score","top_disease"])
    for g,(sc,dis) in sorted(best.items(),key=lambda x:-x[1][0]):
        w.writerow([g,round(sc,5),dis])

# save gene sets at a few score thresholds for transparency
for thr in [0.0,0.05,0.1,0.2]:
    genes=sorted(g for g,(sc,_) in best.items() if sc>=thr)
    json.dump(genes,open(f"data/processed/disease_targets_thr{thr}.json","w"))
    print(f"  score>={thr}: {len(genes)} genes")
print(f"\n[disease_targets] {len(best)} unique disease-associated genes (union) written")
