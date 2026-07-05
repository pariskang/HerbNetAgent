#!/usr/bin/env python3
"""Step 19: architecture diagram + machine-readable module spec for the
next-generation quantitative TCM intelligence agent (HerbNetAgent 2.0)."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

LAYERS=[
 ("L0","Orchestration & Reasoning (LLM agent core)",
  "planner · tool-router · memory · uncertainty & provenance tracker  [ChemCrow/Coscientist-style]","#4A4A6A","steps 01-18 driver"),
 ("L1","Knowledge & Data Layer  (unified knowledge graph)",
  "TCMSP · HERB · BATMAN-TCM 2.0 · SymMap | PubChem · ChEMBL | Open Targets · STRING · KEGG/Reactome | FAERS/TCM-drug PV","#2E6B8A","steps 01-03"),
 ("L2","Molecular Representation & Structure",
  "compound QC (RDKit) · quantitative content (UPLC-MS/MS) · protein structure (AlphaFold3 / ESMFold)","#2E7D5B","steps 01,10"),
 ("L3","Binding & Target-Engagement Engine",
  "co-folding + affinity (Boltz-2 / AlphaFold3 / Chai-1) · docking (DiffDock/Vina) · FEP · target pred (SwissTP/SEA)","#3E8E7E","step 02 (ChEMBL)"),
 ("L4","Network & EFFICACY (yaoxiao) Engine",
  "disease module · network proximity + Complementary Exposure (Cheng/Barabasi) · enrichment · combo synergy","#C08A3C","steps 04,05,14,17"),
 ("L5","ADMET / PK (yaodai) Engine",
  "ADMET-AI · ADMETlab 3.0 (CYP/transporter/tox) -> whole-body PBPK/QSP (PK-Sim/MoBi) · virtual populations · DDI","#C0603C","steps 08,10,12"),
 ("L6","Systems PD & Outcome Simulation",
  "QSP ODE: exposure->target->pathway->phenotype (coagulation-thrombus-bleeding) · tissue models (fracture)","#A0443C","steps 12,15,18"),
 ("L7","Safety, Toxicology & Clinical Translation",
  "ProTox/hERG/DILI · DDI & contraindication · evidence-grading (compute vs clinical) · auto validation protocols","#8A3C5A","steps 09,11,16,briefing"),
]
CROSS=["Uncertainty Quantification","Provenance / Explainability","Human-in-the-loop clinical guardrails"]

fig,ax=plt.subplots(figsize=(14,10))
n=len(LAYERS); H=1.0; gap=0.18; y=n*(H+gap)
xL,wL=0.5,9.6
for i,(code,title,tools,color,proto) in enumerate(LAYERS):
    yy=y-(i+1)*(H+gap)
    ax.add_patch(FancyBboxPatch((xL,yy),wL,H,boxstyle="round,pad=0.02,rounding_size=0.06",
                 linewidth=1.2,edgecolor="white",facecolor=color,alpha=0.92,zorder=2))
    ax.text(xL+0.18,yy+H*0.66,f"{code}  {title}",fontsize=11,weight="bold",color="white",va="center",zorder=3)
    ax.text(xL+0.28,yy+H*0.26,tools,fontsize=8.0,color="white",va="center",zorder=3)
    # prototype tag on the right
    ax.text(xL+wL+0.15,yy+H*0.5,f"prototyped:\n{proto}",fontsize=7.5,color="#333",va="center")
    # downward data-flow arrow
    if i<n-1:
        ax.add_patch(FancyArrowPatch((xL+wL*0.5,yy),(xL+wL*0.5,yy-gap),
                     arrowstyle="-|>",mutation_scale=12,color="#666",lw=1.2,zorder=1))
# cross-cutting bar (left)
ax.add_patch(FancyBboxPatch((-0.15,y-n*(H+gap)+0.02),0.5,n*(H+gap)-0.1,boxstyle="round,pad=0.02",
             facecolor="#DDD",edgecolor="none",alpha=0.8,zorder=0))
ax.text(0.1,y-n*(H+gap)/2,"  |  ".join(CROSS),rotation=90,ha="center",va="center",fontsize=8.5,color="#333",weight="bold")

ax.text(xL+wL/2,y-0.05,"HerbNetAgent 2.0 — Quantitative Multi-Scale TCM Intelligence Agent",
        ha="center",fontsize=14,weight="bold")
ax.text(xL+wL/2,y-0.42,"molecule → structure → affinity → network → exposure(PBPK) → effect(QSP) → safety → evidence-graded report",
        ha="center",fontsize=9.5,style="italic",color="#444")
ax.set_xlim(-0.5,11.2); ax.set_ylim(-0.1,y+0.1); ax.axis("off")
plt.tight_layout(); plt.savefig("results/figures/fig11_agent_architecture.png",dpi=150,bbox_inches="tight"); plt.close()
print("wrote results/figures/fig11_agent_architecture.png")

# machine-readable spec
spec={"name":"HerbNetAgent 2.0","subtitle":"Quantitative Multi-Scale TCM Intelligence Agent",
 "pipeline":"molecule -> structure -> affinity -> network -> PBPK exposure -> QSP effect -> safety -> evidence-graded report",
 "layers":[{"id":c,"title":t,"methods":tools.split(" · "),"prototype":p} for c,t,tools,_,p in LAYERS],
 "cross_cutting":CROSS,
 "efficacy_engine":"L2->L3->L4 (structure->affinity->network proximity/Complementary Exposure)",
 "pk_engine":"L2->L5->L6 (ADMET-AI/ADMETlab -> PBPK PK-Sim -> QSP PD)"}
json.dump(spec,open("data/processed/agent_spec.json","w"),ensure_ascii=False,indent=1)
print("wrote data/processed/agent_spec.json")
