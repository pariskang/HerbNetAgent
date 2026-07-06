#!/usr/bin/env python3
"""
Step 10: Quantitative static mechanistic DDI model for a SPECIFIC regimen:
    rivaroxaban 10 mg qd  +  half a Tongrentang honey pill (1.5 g pill) of
    Dahuang Zhechong Wan.

Question: does halving/lowering the doses meaningfully reduce the interaction?

Approach (FDA 2020 in-vitro DDI static model, [I]/Ki screening):
  For each key inhibitory constituent, estimate the delivered dose in 1.5 g pill,
  the intestinal concentration [I]gut = molar dose / 250 mL, and the systemic
  unbound Cmax; then compute R = [I]/Ki against CYP3A4 / P-gp / BCRP (the
  proteins that clear rivaroxaban), using ChEMBL-measured Ki/IC50 values.
  FDA screening cut-offs: intestinal [I]gut/Ki >= 10, or hepatic 1+[I]u/Ki >= 1.02
  -> interaction cannot be excluded.

ALL constituent-content values are literature/Pharmacopoeia order-of-magnitude
ESTIMATES for free aglycone and carry large uncertainty; results are screening-
level, not a validated PBPK prediction.
"""
import csv, json, math

MW=json.load(open("/tmp/mw.json"))

# ---- Formula composition: Chinese Pharmacopoeia 2020, crude-drug grams per batch ----
FORMULA={  # herb_zh : grams
 "大黄":300,"黄芩":60,"甘草":90,"桃仁":120,"苦杏仁":120,"白芍":120,"地黄":300,
 "干漆":30,"虻虫":45,"水蛭":60,"蛴螬":45,"土鳖虫":30}
TOTAL=sum(FORMULA.values())  # 1320 g crude drug per batch
frac={h:g/TOTAL for h,g in FORMULA.items()}

# ---- Dose being evaluated ----
PILL_MASS=1.5        # g, half of a 3 g Tongrentang 大蜜丸 (honey pill)
HONEY_RATIO=1.1      # 大蜜丸 crude-drug : honey ~ 1 : 1.1  => crude drug fraction
CRUDE_IN_PILL=PILL_MASS/(1+HONEY_RATIO)   # ~0.71 g crude drug equivalent
print(f"Pill mass evaluated: {PILL_MASS} g honey-pill  ->  ~{CRUDE_IN_PILL:.2f} g crude drug equivalent")
print(f"(reference full adult dose ~6 g pill/day = {6/(1+HONEY_RATIO):.2f} g crude drug; ratio {6/PILL_MASS:.0f}x)\n")

# ---- Constituent content: free-aglycone % of the SOURCE herb (ESTIMATE, high uncertainty) ----
# name : (source herb_zh, content_fraction_of_herb, note)
CONTENT={
 "emodin":          ("大黄",0.0010,"游离蒽醌~0.1%"),
 "baicalein":       ("黄芩",0.0050,"苷元~0.5%(苷 baicalin 更高但非该抑制剂)"),
 "quercetin":       ("黄芩",0.0002,"游离槲皮素极低,多为糖苷"),
 "apigenin":        ("黄芩",0.0002,"游离芹菜素低"),
 "kaempferol":      ("甘草",0.0002,"游离低"),
 "luteolin":        ("黄芩",0.0002,"游离低"),
 "isoliquiritigenin":("甘草",0.0003,"查耳酮,低"),
 "licochalcone A":  ("甘草",0.0002,"仅胀果甘草含量高,乌拉尔甘草极低-高度不确定"),
 "naringenin":      ("甘草",0.0003,"低"),
 "glabridin":       ("甘草",0.0010,"光甘草定,光果甘草-不确定"),
}

# ---- Ki / IC50 (nM) from ChEMBL DDI table ----
ddi=list(csv.DictReader(open("results/tables/rivaroxaban_ddi.csv")))
Ki={}  # (compound,target)->nM
for r in ddi:
    try: v=float(r["value"])
    except: continue
    key=(r["compound"],r["target"])
    if key not in Ki or v<Ki[key]: Ki[key]=v   # most potent (lowest) value
TARGETS=["CYP3A4","ABCB1/P-gp","ABCG2/BCRP"]

GUT_VOL=0.25  # L, FDA intestinal model
# crude systemic assumptions for oral flavonoids (order of magnitude)
F_ABS=0.10          # ~10% fraction reaching portal (low BA flavonoids)
VD=100.0            # L, large apparent Vd
FU=0.05             # unbound fraction ~5%

def model_at(crude_g, label):
    rows=[]
    print(f"===== {label}  (crude drug {crude_g:.2f} g) =====")
    for name,(herb,cfrac,note) in CONTENT.items():
        if name not in MW: continue
        mw=MW[name]
        herb_mass_g = crude_g*frac.get(herb,0)           # g of source herb in the dose
        dose_mg = herb_mass_g*1000*cfrac                  # mg constituent delivered
        dose_umol = dose_mg/mw*1000                        # µmol
        Igut_uM = dose_umol/GUT_VOL/1000                   # µmol/L*(1/1000?) -> compute in µM
        # dose_umol/0.25 L = µmol/L = µM  (careful): µmol / L = µM
        Igut_uM = dose_umol/GUT_VOL                         # this is µmol/L = µM
        # systemic unbound Cmax estimate
        Cmax_uM = (F_ABS*dose_umol)/VD                     # µmol/L = µM (total)
        Cmax_u_uM = Cmax_uM*FU
        rec={"constituent":name,"herb":herb,"dose_mg":round(dose_mg,3),"dose_umol":round(dose_umol,3),
             "Igut_uM":round(Igut_uM,3),"Cmax_u_uM":round(Cmax_u_uM,5),"note":note}
        for t in TARGETS:
            ki_nM=Ki.get((name,t))
            if ki_nM:
                ki_uM=ki_nM/1000.0
                rec[f"{t}_Ki_uM"]=round(ki_uM,4)
                rec[f"{t}_Igut/Ki"]=round(Igut_uM/ki_uM,2)
                rec[f"{t}_R1sys"]=round(1+Cmax_u_uM/ki_uM,4)
            else:
                rec[f"{t}_Ki_uM"]="";rec[f"{t}_Igut/Ki"]="";rec[f"{t}_R1sys"]=""
        rows.append(rec)
    # print flagged
    for r in rows:
        flags=[]
        for t in TARGETS:
            g=r.get(f"{t}_Igut/Ki")
            if isinstance(g,(int,float)) and g>=10: flags.append(f"{t} Igut/Ki={g} ⚠️>10")
        fl=" | ".join(flags) if flags else "systemic & intestinal below FDA threshold"
        print(f"  {r['constituent']:18s} dose={r['dose_mg']:.3f}mg  Igut={r['Igut_uM']:.2f}µM  -> {fl}")
    print()
    return rows

rows_half=model_at(CRUDE_IN_PILL, "评估方案: 1.5 g 蜜丸 (半粒)")
rows_full=model_at(6/(1+HONEY_RATIO), "对照: 6 g 蜜丸/日 (常规成人量)")

with open("results/tables/dose_ddi_model.csv","w",newline="") as f:
    fn=list(rows_half[0].keys())
    w=csv.DictWriter(f,fieldnames=fn); w.writeheader()
    for r in rows_half: w.writerow(r)

# save for report
json.dump({"pill_mass":PILL_MASS,"crude_in_pill":round(CRUDE_IN_PILL,3),
           "half":rows_half,"full":rows_full},
          open("data/processed/dose_ddi_model.json","w"),ensure_ascii=False,indent=1)
print("wrote results/tables/dose_ddi_model.csv and json")
