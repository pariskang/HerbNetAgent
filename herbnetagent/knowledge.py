"""
herbnetagent.knowledge — L1 data layer.

Cache-first access to the primary databases. Anything already computed by the
01-19 pipeline lives in data/processed and is reused for a fast, deterministic
run; when a datum is missing we fall back to the live API (PubChem/ChEMBL/
Open Targets/STRING). This is the layer that would grow into a knowledge graph.
"""
from __future__ import annotations

import csv
import json
import os
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "data", "processed")
TABLES = os.path.join(ROOT, "results", "tables")

try:
    import requests
    _UA = {"User-Agent": "HerbNetAgent/2.0 (research)"}
except Exception:  # requests always present in this env, but stay defensive
    requests = None
    _UA = {}


def _read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def _read_json(path):
    with open(path) as f:
        return json.load(f)


class Knowledge:
    """Unified accessor. Prefers the on-disk cache, falls back to live APIs."""

    def __init__(self, proc: str = PROC, tables: str = TABLES):
        self.proc, self.tables = proc, tables

    # ---- formula / compounds -------------------------------------------------
    def compounds(self) -> list[dict]:
        return _read_csv(os.path.join(self.proc, "compounds.csv"))

    def herbs(self) -> list[dict]:
        return _read_csv(os.path.join(self.proc, "herbs.csv"))

    # ---- compound -> targets -------------------------------------------------
    def compound_targets(self) -> list[dict]:
        return _read_csv(os.path.join(self.proc, "compound_targets.csv"))

    def formula_targets(self) -> set[str]:
        return set(_read_json(os.path.join(self.proc, "drug_targets.json")))

    # ---- disease -> targets --------------------------------------------------
    def disease_targets(self, min_score: float = 0.0) -> dict[str, float]:
        rows = _read_csv(os.path.join(self.proc, "disease_targets.csv"))
        return {r["gene"]: float(r["max_assoc_score"])
                for r in rows if float(r["max_assoc_score"]) >= min_score}

    # ---- live fallbacks (used only when cache misses) ------------------------
    def chembl_targets_for_drug(self, chembl_id: str, pchembl_min: float = 6.0) -> tuple:
        """Measured human targets for an arbitrary drug (e.g. a co-administered
        Western drug not in the cached formula). Returns (gene, max_pchembl)."""
        if requests is None:
            return ()
        base = "https://www.ebi.ac.uk/chembl/api/data"
        try:
            r = requests.get(f"{base}/activity.json", headers=_UA, timeout=40,
                             params={"molecule_chembl_id": chembl_id,
                                     "pchembl_value__gte": pchembl_min,
                                     "target_organism": "Homo sapiens", "limit": 200})
            acts = r.json().get("activities", [])
            best: dict[str, float] = {}
            for a in acts:
                t, pv = a.get("target_chembl_id"), a.get("pchembl_value")
                if t and pv:
                    best[t] = max(best.get(t, 0), float(pv))
            out = []
            for t, pv in best.items():
                d = requests.get(f"{base}/target/{t}.json", headers=_UA, timeout=30).json()
                if d.get("target_type") != "SINGLE PROTEIN":
                    continue
                for c in d.get("target_components", []):
                    for s in c.get("target_component_synonyms", []):
                        if s.get("syn_type") == "GENE_SYMBOL":
                            out.append((s["component_synonym"], round(pv, 2)))
                            break
                time.sleep(0.03)
            return tuple(out)
        except Exception:
            return ()

    # ---- known co-drug target sets (curated, with ChEMBL id) -----------------
    KNOWN_CODRUGS = {
        "rivaroxaban": {"chembl": "CHEMBL198362",
                        "targets": {"F10": 9.4, "F2": 6.1},
                        "class": "direct FXa (F10) inhibitor",
                        "cleared_by": ["CYP3A4", "CYP2J2", "ABCB1/P-gp", "ABCG2/BCRP"]},
        "warfarin": {"chembl": "CHEMBL1464",
                     "targets": {"VKORC1": 7.0}, "class": "VKORC1 inhibitor",
                     "cleared_by": ["CYP2C9", "CYP3A4"]},
        "aspirin": {"chembl": "CHEMBL25",
                    "targets": {"PTGS1": 7.0, "PTGS2": 6.0}, "class": "COX inhibitor",
                    "cleared_by": ["esterase"]},
    }

    def codrug(self, name: str) -> dict | None:
        info = self.KNOWN_CODRUGS.get(name.lower())
        if info:
            return dict(info, name=name.lower())
        return None
