#!/usr/bin/env python3
"""
Step 8: Rivaroxaban (利伐沙班) co-administration analysis for Dahuang Zhechong Wan.

Two interaction axes, both data-driven from ChEMBL measured bioactivity:

(A) Pharmacokinetic (PK): rivaroxaban is cleared by CYP3A4/CYP2J2 and is a
    substrate of P-glycoprotein (ABCB1) and BCRP (ABCG2). We test whether the
    46 formula constituents *inhibit* these enzymes/transporters (which would
    raise rivaroxaban exposure -> bleeding risk).
(B) Pharmacodynamic (PD): rivaroxaban is a direct Factor Xa (F10) inhibitor.
    We check whether the formula engages the SAME coagulation/platelet axis
    (F2 thrombin, PTGS1/2, ALOX12/15, F10) -> additive anticoagulation.
"""
import requests, csv, json, time, os
from collections import defaultdict
UA={'User-Agent':'Mozilla/5.0 (research/ddi)'}
CHEMBL="https://www.ebi.ac.uk/chembl/api/data"

ADME_PD={
 "CHEMBL340":  ("CYP3A4",  "PK: 主要代谢酶 (rivaroxaban ~51% CYP)"),
 "CHEMBL3491": ("CYP2J2",  "PK: 次要代谢酶"),
 "CHEMBL4302": ("ABCB1/P-gp","PK: 外排转运体 (rivaroxaban 底物)"),
 "CHEMBL5393": ("ABCG2/BCRP","PK: 外排转运体 (rivaroxaban 底物)"),
 "CHEMBL244":  ("F10/FXa", "PD: 利伐沙班的直接作用靶点"),
}

def get(url, params=None, tries=4):
    for i in range(tries):
        try:
            r=requests.get(url, params=params, headers=UA, timeout=40)
            if r.status_code==200: return r.json()
            if r.status_code==404: return None
        except Exception: pass
        time.sleep(1.2*(i+1))
    return None

# --- map constituents -> ChEMBL molecule id (cache) ---
cache_path="data/processed/compound_chembl.csv"
comp2chembl={}
if os.path.exists(cache_path):
    for r in csv.DictReader(open(cache_path)):
        comp2chembl[r["compound"]]=r["chembl_id"]
else:
    comps=list(csv.DictReader(open("data/processed/compounds.csv")))
    rows=[]
    for c in comps:
        cid=c["pubchem_cid"]
        if not cid: continue
        d=get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/InChIKey/JSON")
        try: ik=d["PropertyTable"]["Properties"][0]["InChIKey"]
        except Exception: ik=None
        time.sleep(0.05)
        mid=None
        if ik:
            d2=get(f"{CHEMBL}/molecule.json",{'molecule_structures__standard_inchi_key':ik,'limit':1})
            try:
                ms=d2.get("molecules",[]); mid=ms[0]["molecule_chembl_id"] if ms else None
            except Exception: pass
        time.sleep(0.05)
        if mid: comp2chembl[c["compound"]]=mid; rows.append({"compound":c["compound"],"chembl_id":mid})
        print(f"  {c['compound']:24s} -> {mid}")
    with open(cache_path,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["compound","chembl_id"]); w.writeheader(); w.writerows(rows)

chembl2comp={v:k for k,v in comp2chembl.items()}
mol_ids=list(comp2chembl.values())
print(f"\nmapped {len(mol_ids)} constituents to ChEMBL molecules")

# --- (A) PK: query activities of constituents against each ADME/PD target ---
def activities_for_target(tid, mol_ids):
    """Return list of (molecule_chembl_id, pchembl, std_type, std_value, units)."""
    out=[]
    # chunk molecule ids to keep URL length sane
    for i in range(0,len(mol_ids),20):
        chunk=",".join(mol_ids[i:i+20])
        d=get(f"{CHEMBL}/activity.json",
              {'target_chembl_id':tid,'molecule_chembl_id__in':chunk,'limit':1000})
        while d:
            for a in d.get("activities",[]):
                out.append((a.get("molecule_chembl_id"),a.get("pchembl_value"),
                            a.get("standard_type"),a.get("standard_value"),a.get("standard_units")))
            nxt=d.get("page_meta",{}).get("next")
            d=get("https://www.ebi.ac.uk"+nxt) if nxt else None
    return out

ddi_rows=[]
summary=defaultdict(list)
for tid,(name,role) in ADME_PD.items():
    acts=activities_for_target(tid, mol_ids)
    # best pchembl per molecule
    best={}
    for mid,pv,stype,sval,units in acts:
        if pv is None: continue
        pv=float(pv)
        if mid not in best or pv>best[mid][0]: best[mid]=(pv,stype,sval,units)
    print(f"\n### {name} ({tid}) — {role}")
    for mid,(pv,stype,sval,units) in sorted(best.items(),key=lambda x:-x[1][0]):
        cpd=chembl2comp.get(mid,mid)
        pot = "strong(<1µM)" if pv>=6 else ("moderate(1-10µM)" if pv>=5 else "weak(>10µM)")
        print(f"   {cpd:22s} pChEMBL={pv:.2f} [{pot}]  {stype}={sval}{units or ''}")
        ddi_rows.append({"axis":"PK" if tid!="CHEMBL244" else "PD","target":name,"target_chembl":tid,
                         "role":role,"compound":cpd,"pchembl":round(pv,2),"potency":pot,
                         "assay_type":stype,"value":sval,"units":units or ""})
        summary[name].append((cpd,pv))
    if not best: print("   (no measured activity among constituents)")

with open("results/tables/rivaroxaban_ddi.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["axis","target","target_chembl","role","compound","pchembl","potency","assay_type","value","units"])
    w.writeheader(); w.writerows(ddi_rows)

# --- (B) PD: formula's coagulation/platelet axis vs rivaroxaban (F10) ---
drug=set(json.load(open("data/processed/drug_targets.json")))
COAG_PLT={"F2":"凝血酶(thrombin)","F10":"Xa因子(rivaroxaban 靶点)","PTGS1":"COX-1(血小板)",
          "PTGS2":"COX-2(血小板/炎症)","ALOX12":"12-脂氧合酶(血小板)","ALOX15":"15-脂氧合酶",
          "PLG":"纤溶酶原","PLAT":"tPA","SERPINE1":"PAI-1"}
pd_hits={g:d for g,d in COAG_PLT.items() if g in drug}
json.dump({"pk_modulation":{k:[(c,round(p,2)) for c,p in v] for k,v in summary.items()},
           "pd_coag_platelet_overlap":pd_hits},
          open("data/processed/rivaroxaban_ddi_summary.json","w"),ensure_ascii=False,indent=1)

print("\n=== PD overlap (formula 命中的凝血/血小板靶点，与 rivaroxaban 同一止血系统) ===")
for g,d in pd_hits.items(): print(f"   {g:8s} {d}")
print("\nwrote results/tables/rivaroxaban_ddi.csv and summary json")
