#!/usr/bin/env python3
"""
Step 12: Dynamic (reduced-)PBPK model of rivaroxaban ± Dahuang Zhechong Wan
constituents, to quantify the fold-change in rivaroxaban exposure (AUCR).

This upgrades the earlier static [I]/Ki screen into a TIME-COURSE model:
a perpetrator (quercetin, real human PK) is dosed to steady state; its
time-varying plasma (hepatic/renal) and gut-lumen concentrations dynamically
modulate rivaroxaban clearance and gut availability through:
  - hepatic CYP3A4 reversible inhibition        (fm,CYP3A4 of rivaroxaban CL)
  - renal P-gp/BCRP active-secretion inhibition (fraction of renal CL)
  - intestinal CYP3A4 inhibition                (raises gut availability Fg)
  - intestinal BCRP/P-gp efflux inhibition      (raises absorbed fraction Fa)
Rivaroxaban AUC over the steady-state interval is compared with vs without the
perpetrator: AUCR = AUC_inhibited / AUC_control.

Model is calibrated against the known ketoconazole interaction (clinical AUCR
~2.6). Human PK parameters are literature values (see PARAMS + report).
NOTE: reduced/semi-mechanistic PBPK, not a full whole-body organ model.
"""
import numpy as np
from scipy.integrate import solve_ivp
import json, csv

# ---------------- Rivaroxaban (victim) human PK ----------------
RIV = dict(
    MW=435.9, dose_mg=10.0,
    ka=1.24,          # /h absorption
    V=50.0,           # L central volume
    fu=0.075,         # unbound plasma fraction (~92.5% bound)
    thalf=6.0,        # h terminal (used to derive CL)
    Fa0=0.95, Fg0=0.95,   # baseline absorbed fraction / gut-wall availability
    # fractional contributions to TOTAL clearance (calibrated to ketoconazole AUCR~2.6)
    f_CYP3A4=0.30, f_CYP2J2=0.14, f_hydrolysis=0.10,
    f_renal_secretion=0.30,   # active tubular secretion via P-gp/BCRP
    f_other=0.16,             # passive renal + biliary (uninhibited)
)
RIV["CL"]=0.693*RIV["V"]/RIV["thalf"]       # L/h total
RIV["dose_umol"]=RIV["dose_mg"]/RIV["MW"]*1000

# ---------------- Quercetin (perpetrator) human PK ----------------
# Free aglycone drives inhibition; oral bioavailability is very low.
QUER = dict(
    MW=302.23,
    ka=0.6, V=200.0, F=0.02, fu=0.015,   # low BA, high binding
    CL=70.0,                              # L/h (apparent, gives short aglycone t1/2)
    Ki_CYP3A4=1.0,    # µM (IC50~2.1 -> Ki~1)
    Ki_pgp=7.08,      # µM
    Ki_bcrp=0.03,     # µM (30 nM, potent)
)
Vgut=0.30   # L, gut-lumen distribution volume for luminal transporter/enzyme inhibition
Qh=90.0     # L/h hepatic blood flow (for hepatic-inlet estimate)

def simulate_AUCR(quercetin_mg, days=20, strong_inhibitor=None):
    """Return rivaroxaban AUC(0-24) at steady state with the given daily quercetin
    dose co-administered. strong_inhibitor=dict overrides A-terms for validation."""
    q_dose_umol = quercetin_mg/QUER["MW"]*1000 if quercetin_mg else 0.0

    def rhs(t, y):
        Aqg, Aqc, Arg, Arc = y
        # perpetrator concentrations
        Cq = Aqc/QUER["V"]                       # µM total plasma quercetin
        Iq_sys_u = QUER["fu"]*Cq                 # unbound systemic (renal transporters)
        # hepatic inlet unbound (systemic + portal absorption term)
        portal = (QUER["F"]*QUER["ka"]*Aqg)/Qh
        Iq_h_u = QUER["fu"]*(Cq + portal)
        Iq_gut = Aqg/Vgut                        # µM luminal (gut enzymes/transporters)

        if strong_inhibitor:                     # validation override (e.g. ketoconazole)
            A_h3A4  = strong_inhibitor["A_h3A4"]
            A_renal = strong_inhibitor["A_renal"]
            A_g3A4  = strong_inhibitor["A_g3A4"]
            A_geff  = strong_inhibitor["A_geff"]
        else:
            A_h3A4  = 1.0/(1.0 + Iq_h_u/QUER["Ki_CYP3A4"])
            A_renal = 1.0/(1.0 + Iq_sys_u/QUER["Ki_pgp"] + Iq_sys_u/QUER["Ki_bcrp"])
            A_g3A4  = 1.0/(1.0 + Iq_gut/QUER["Ki_CYP3A4"])
            A_geff  = 1.0/(1.0 + Iq_gut/QUER["Ki_bcrp"] + Iq_gut/QUER["Ki_pgp"])

        # rivaroxaban time-varying clearance & availability
        CL = RIV["CL"]*(RIV["f_CYP3A4"]*A_h3A4 + RIV["f_CYP2J2"] + RIV["f_hydrolysis"]
                        + RIV["f_renal_secretion"]*A_renal + RIV["f_other"])
        Fg = RIV["Fg0"] + (1-RIV["Fg0"])*(1-A_g3A4)
        Fa = RIV["Fa0"] + (1-RIV["Fa0"])*(1-A_geff)
        Feff = Fa*Fg

        dAqg = -QUER["ka"]*Aqg
        dAqc = QUER["ka"]*Aqg*QUER["F"] - (QUER["CL"]/QUER["V"])*Aqc
        dArg = -RIV["ka"]*Arg
        dArc = RIV["ka"]*Arg*Feff - (CL/RIV["V"])*Arc
        return [dAqg, dAqc, dArg, dArc]

    y=[0,0,0,0]
    last_auc=None
    for d in range(days):
        y[0]+=q_dose_umol      # dose perpetrator into gut
        y[2]+=RIV["dose_umol"] # dose rivaroxaban into gut
        sol=solve_ivp(rhs,[0,24],y,max_step=0.25,rtol=1e-7,atol=1e-9,dense_output=True)
        y=[sol.y[i,-1] for i in range(4)]
        if d==days-1:
            tt=np.linspace(0,24,481)
            Cr=sol.sol(tt)[3]/RIV["V"]     # µM
            last_auc=np.trapezoid(Cr,tt)       # µM·h
    return last_auc

# ---- baseline (no inhibitor) ----
auc0=simulate_AUCR(0)
print(f"Rivaroxaban baseline AUC(0-24,ss) = {auc0:.4f} µM·h  (CL={RIV['CL']:.2f} L/h)\n")

# ---- validation: strong CYP3A4 + P-gp inhibitor (ketoconazole-like) ----
strong=dict(A_h3A4=0.03, A_renal=0.12, A_g3A4=0.05, A_geff=0.05)  # near-complete inhibition
auc_keto=simulate_AUCR(0, strong_inhibitor=strong)
print(f"[validation] strong CYP3A4+P-gp inhibitor  -> AUCR = {auc_keto/auc0:.2f}x  (clinical ketoconazole ~2.6x)\n")

# ---- formula-delivered quercetin dose (from 1.5 g pill dose model) ----
M=json.load(open("data/processed/dose_ddi_model.json"))
q_half=[r for r in M["half"] if r["constituent"]=="quercetin"][0]["dose_mg"]
q_full=[r for r in M["full"] if r["constituent"]=="quercetin"][0]["dose_mg"]
print(f"quercetin delivered: half pill(1.5g)={q_half:.4f} mg | full day(6g)={q_full:.4f} mg\n")

# ---- dose-response sweep ----
sweep=[q_half, q_full, 1.0, 10.0, 50.0, 100.0, 250.0, 500.0, 1000.0]
labels={q_half:"formula 1.5g pill", q_full:"formula 6g/day", 500.0:"typical quercetin supplement"}
results=[]
print("quercetin_dose_mg   AUCR   note")
for dose in sorted(set(sweep)):
    r=simulate_AUCR(dose)/auc0
    note=labels.get(dose,"")
    results.append({"quercetin_mg":dose,"AUCR":round(r,4),"pct_change":round((r-1)*100,2),"note":note})
    print(f"  {dose:10.4f}       {r:.4f}   {note}")

with open("results/tables/pbpk_aucr.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["quercetin_mg","AUCR","pct_change","note"]); w.writeheader(); w.writerows(results)

# threshold: quercetin dose to reach AUCR 1.25 (clinically notable)
def auc_dose(dose): return simulate_AUCR(dose)/auc0
lo,hi=100.0,3000.0
for _ in range(30):
    mid=(lo+hi)/2
    if auc_dose(mid)<1.25: lo=mid
    else: hi=mid
thr125=(lo+hi)/2
print(f"\nquercetin dose to reach rivaroxaban AUCR=1.25 : ~{thr125:.0f} mg (~{thr125/max(q_half,1e-9):.0f}x the 1.5g-pill content)")

summary=dict(baseline_auc=round(auc0,4),
             validation_strong_AUCR=round(auc_keto/auc0,2),
             quercetin_half_mg=q_half, quercetin_full_mg=q_full,
             AUCR_half=round(simulate_AUCR(q_half)/auc0,4),
             AUCR_full=round(simulate_AUCR(q_full)/auc0,4),
             AUCR_500mg_supplement=round(simulate_AUCR(500)/auc0,3),
             quercetin_mg_for_AUCR1p25=round(thr125,0))
json.dump(summary,open("data/processed/pbpk_summary.json","w"),indent=1)

# ---- save steady-state rivaroxaban concentration-time profiles for plotting ----
def profile(quercetin_mg=0, strong=None, days=20):
    q=quercetin_mg/QUER["MW"]*1000 if quercetin_mg else 0.0
    def rhs(t,y):  # reuse global model via closure by re-calling simulate internals
        return None
    # re-run simulate but capture final-day profile
    y=[0,0,0,0]
    for d in range(days):
        y[0]+=q; y[2]+=RIV["dose_umol"]
        # inline rhs identical to simulate_AUCR
        def f(t,yy):
            Aqg,Aqc,Arg,Arc=yy
            Cq=Aqc/QUER["V"]; Iq_sys_u=QUER["fu"]*Cq
            portal=(QUER["F"]*QUER["ka"]*Aqg)/Qh; Iq_h_u=QUER["fu"]*(Cq+portal)
            Iq_gut=Aqg/Vgut
            if strong:
                A_h3A4,A_renal,A_g3A4,A_geff=strong["A_h3A4"],strong["A_renal"],strong["A_g3A4"],strong["A_geff"]
            else:
                A_h3A4=1/(1+Iq_h_u/QUER["Ki_CYP3A4"]); A_renal=1/(1+Iq_sys_u/QUER["Ki_pgp"]+Iq_sys_u/QUER["Ki_bcrp"])
                A_g3A4=1/(1+Iq_gut/QUER["Ki_CYP3A4"]); A_geff=1/(1+Iq_gut/QUER["Ki_bcrp"]+Iq_gut/QUER["Ki_pgp"])
            CL=RIV["CL"]*(RIV["f_CYP3A4"]*A_h3A4+RIV["f_CYP2J2"]+RIV["f_hydrolysis"]+RIV["f_renal_secretion"]*A_renal+RIV["f_other"])
            Fg=RIV["Fg0"]+(1-RIV["Fg0"])*(1-A_g3A4); Fa=RIV["Fa0"]+(1-RIV["Fa0"])*(1-A_geff)
            return [-QUER["ka"]*Aqg, QUER["ka"]*Aqg*QUER["F"]-(QUER["CL"]/QUER["V"])*Aqc,
                    -RIV["ka"]*Arg, RIV["ka"]*Arg*Fa*Fg-(CL/RIV["V"])*Arc]
        sol=solve_ivp(f,[0,24],y,max_step=0.25,rtol=1e-7,atol=1e-9,dense_output=True)
        y=[sol.y[i,-1] for i in range(4)]
        if d==days-1:
            tt=np.linspace(0,24,241)
            ngml=sol.sol(tt)[3]/RIV["V"]*RIV["MW"]   # µM -> ng/mL (µmol/L * g/mol = mg/L = µg/mL? -> *1000 ng)
            return tt.tolist(),(sol.sol(tt)[3]/RIV["V"]*RIV["MW"]).tolist()  # µg/L = ng/mL? µM*MW= µg/L
prof={}
tt,base=profile(0); prof["t"]=tt; prof["baseline"]=base
_,ph=profile(q_half); prof["formula_1p5g"]=ph
_,ps=profile(0,strong=strong); prof["strong_inhibitor"]=ps
json.dump(prof,open("data/processed/pbpk_profiles.json","w"))
print("\nSUMMARY:",json.dumps(summary,indent=1))
print("wrote results/tables/pbpk_aucr.csv, pbpk_summary.json, pbpk_profiles.json")
