#!/usr/bin/env python3
"""
Step 20: structure-based binding affinity for the two headline targets
F2 (thrombin) and MMP9 (gelatinase B) via the L3 Boltz-2 backend.

Runs Boltz-2 for real when the `boltz` binary + a GPU are present (writes YAML,
co-folds protein+ligand, reads affinity_pred_value). Otherwise it reports that
the structure-based engine is unavailable and shows the measured ChEMBL value the
BindingEngine falls back to — so the pipeline is honest either way.

    # with a GPU + `pip install boltz`:
    python3 scripts/20_boltz2_affinity.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from herbnetagent.knowledge import Knowledge
from herbnetagent.structure import StructureProvider, UNIPROT
from herbnetagent.boltz2 import Boltz2Backend
from herbnetagent.binding import ChEMBLBackend, BindingEngine

PROBES = {
    "F2":   ["quercetin", "beta-sitosterol"],
    "MMP9": ["isoliquiritigenin", "quercetin", "luteolin"],
}

k = Knowledge()
sp = StructureProvider()
b2 = Boltz2Backend(knowledge=k, structure_provider=sp)
measured = ChEMBLBackend(k.compound_targets())
engine = BindingEngine([b2, measured])

print(f"Boltz-2 available: {b2.available()}  (binary on PATH + GPU)\n")

for gene, compounds in PROBES.items():
    acc = UNIPROT.get(gene)
    print(f"### {gene}  (UniProt {acc}, {os.path.basename(sp.cache)}/AF-{acc}.pdb)")
    for cpd in compounds:
        if b2.available():
            q = b2.affinity(cpd, gene)
            tag = "Boltz-2 predicted"
        else:
            q = None
            tag = "predicted N/A -> measured fallback"
        if q is None:
            q = measured.affinity(cpd, gene)
        line = q.graded() if q else "no data"
        print(f"   {cpd:20s} [{tag:32s}] {line}")
    print()

print("Note: with `boltz` installed the F2/MMP9 rows show predicted IC50 (nM) +"
      " P(binder); here they show the measured ChEMBL fallback.")
