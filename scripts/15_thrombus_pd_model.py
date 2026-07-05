#!/usr/bin/env python3
"""
Step 15: CONCEPTUAL pharmacodynamic model of muscular-vein thrombus resolution,
to formalise the "synergy / faster healing" hypothesis and generate testable
predictions with an explicit bleeding trade-off.

*** ILLUSTRATIVE MODEL — rate constants are assumed, NOT fitted to clinical data.
    It shows the STRUCTURE of the hypothesis and what to measure, not proof. ***

Thrombus burden T(t) (1.0 = initial):
  dT/dt = G * (fxa*thr*plt) * T*(Tmax-T)      # propagation (coagulation+platelet)
          - K * res_eff * T                    # resolution (fibrinolysis+MMP+anti-inflam)
  res_eff = res * (1 + 0.6*(1-infl))           # lower inflammation -> faster resolution
Node activities (0-1) per scenario are derived from the mechanism-coverage map.
Bleeding proxy = weighted systemic anticoagulant+antiplatelet intensity.
"""
import numpy as np
from scipy.integrate import solve_ivp
import json, csv

# node activities (fraction of normal) under each regimen — illustrative, from coverage map
SC={
 "Control":       dict(fxa=1.00, thr=1.00, plt=1.00, infl=1.00, res=1.00),
 "Rivaroxaban 10mg": dict(fxa=0.15, thr=0.45, plt=1.00, infl=0.90, res=1.00),
 "Formula 1.5g":  dict(fxa=1.00, thr=0.60, plt=0.72, infl=0.60, res=1.45),
 "Combination":   dict(fxa=0.15, thr=0.32, plt=0.72, infl=0.55, res=1.45),
}
G=0.20; K=0.045; Tmax=1.3; T0=1.0; DAYS=35

def run(p):
    def rhs(t,T):
        T=T[0]
        prop=G*(p["fxa"]*p["thr"]*p["plt"])*T*(Tmax-T)
        res_eff=p["res"]*(1+0.6*(1-p["infl"]))
        reso=K*res_eff*T
        return [prop-reso]
    sol=solve_ivp(rhs,[0,DAYS],[T0],max_step=0.2,dense_output=True,rtol=1e-8,atol=1e-10)
    tt=np.linspace(0,DAYS,DAYS*4+1); Tt=sol.sol(tt)[0]
    # time to 50% resolution
    t50=None
    for i in range(len(tt)):
        if Tt[i]<=0.5*T0: t50=tt[i]; break
    resid=Tt[-1]
    return tt,Tt,t50,resid

def bleeding_proxy(p):
    return round(0.5*(1-p["fxa"])+0.3*(1-p["thr"])+0.2*(1-p["plt"]),3)

profiles={}; rows=[]
print(f"{'Scenario':18s} {'t50(d)':>7s} {'residual@35d':>12s} {'bleed_proxy':>11s} {'benefit/bleed':>13s}")
for name,p in SC.items():
    tt,Tt,t50,resid=run(p)
    bp=bleeding_proxy(p)
    speed=(1/t50) if t50 else 0
    ratio=round(speed/bp,3) if bp>0 else float('inf')
    profiles[name]={"t":tt.tolist(),"T":Tt.tolist()}
    rows.append({"scenario":name,"t50_days":round(t50,1) if t50 else None,
                 "residual_day35":round(resid,3),"bleeding_proxy":bp,
                 "resolution_per_bleed":ratio})
    print(f"{name:18s} {(round(t50,1) if t50 else '—'):>7} {resid:>12.3f} {bp:>11.3f} {str(ratio):>13}")

with open("results/tables/thrombus_pd_model.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["scenario","t50_days","residual_day35","bleeding_proxy","resolution_per_bleed"])
    w.writeheader(); w.writerows(rows)
json.dump({"scenarios":SC,"rows":rows,"profiles":profiles},
          open("data/processed/thrombus_pd_model.json","w"),ensure_ascii=False)
print("\nwrote results/tables/thrombus_pd_model.csv and json")
print("NOTE: illustrative/conceptual model — predictions to be tested, not proof.")
