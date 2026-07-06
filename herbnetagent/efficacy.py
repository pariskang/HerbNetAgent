"""
herbnetagent.efficacy — L4 network & efficacy engine.

Implements the quantitative network-medicine metrics the architecture calls for:
  * network proximity (Menche/Guney/Barabasi): closeness of a drug's targets to a
    disease module on the interactome, with a degree-preserving null z-score;
  * Complementary Exposure (Cheng et al., Nat Commun 2019): a drug *pair* is a good
    combination when BOTH sets are proximal to the disease module AND the two sets
    are topologically SEPARATED from each other (complementary, not redundant);
  * mechanism-axis coverage: how many pathophysiology axes each agent + the
    combination touches.

The interactome here is a STRING-induced subgraph over the relevant genes
(cached to data/processed/interactome_edges.tsv). This is an approximation of the
full interactome and is labelled as such — the metric structure is the point.
"""
from __future__ import annotations

import os
import random
import statistics

import networkx as nx

from .core import zscore

_HERE = os.path.dirname(os.path.abspath(__file__))
_CACHE = os.path.join(os.path.dirname(_HERE), "data", "processed", "interactome_edges.tsv")

# venous-thrombosis pathophysiology axes (also used by step 14)
AXES = {
    "Coagulation":   ["F2", "F10", "F5", "F7", "F3", "F9", "SERPINC1", "PROC", "PROS1", "THBD", "TFPI"],
    "Fibrinolysis":  ["PLG", "PLAT", "PLAU", "SERPINE1", "SERPINB2"],
    "Platelet":      ["PTGS1", "PTGS2", "ALOX12", "ALOX15", "TBXA2R", "P2RY12",
                      "P2RY1", "ITGB3", "GP6", "GP1BA", "GRK6"],
    "Thrombo-inflammation": ["NFKB1", "NLRP3", "JAK2", "STAT6", "IL2", "IL6",
                             "IL1B", "TNF", "PIK3R1", "TLR4", "CXCL8"],
    "Resolution":    ["MMP2", "MMP9", "MMP12", "MMP14", "TIMP1"],
    "Endothelium":   ["KDR", "VEGFA", "NOS3", "HIF1A", "VWF", "SELE", "ICAM1", "VCAM1"],
}


def _fetch_string_edges(genes, required_score=400):
    import requests
    r = requests.post("https://string-db.org/api/tsv/network",
                      data={"identifiers": "\r".join(sorted(genes)), "species": 9606,
                            "required_score": required_score},
                      headers={"User-Agent": "HerbNetAgent/2.0"}, timeout=90)
    edges = []
    for i, line in enumerate(r.text.strip().splitlines()):
        if i == 0:
            continue
        p = line.split("\t")
        if len(p) >= 6:
            edges.append((p[2], p[3]))
    return edges


class InteractomeProximity:
    """Network-proximity + separation on a (cached) STRING-induced interactome."""

    def __init__(self, graph: nx.Graph):
        self.G = graph
        self._deg_bins = self._bin_by_degree()

    # --- construction ---------------------------------------------------------
    @classmethod
    def build(cls, genes, use_cache=True):
        edges = []
        if use_cache and os.path.exists(_CACHE):
            with open(_CACHE) as f:
                edges = [tuple(ln.split("\t")[:2]) for ln in f.read().strip().splitlines()]
        if not edges:
            try:
                edges = _fetch_string_edges(genes)
                with open(_CACHE, "w") as f:
                    f.write("\n".join(f"{a}\t{b}" for a, b in edges))
            except Exception:
                edges = []
        G = nx.Graph()
        G.add_nodes_from(genes)
        G.add_edges_from(edges)
        return cls(G)

    def _bin_by_degree(self):
        bins: dict[int, list[str]] = {}
        for n, d in self.G.degree():
            bins.setdefault(d, []).append(n)
        return bins

    # --- metrics --------------------------------------------------------------
    def _mean_min_dist(self, A, B):
        A = [g for g in A if g in self.G]
        B = [g for g in B if g in self.G]
        if not A or not B:
            return None
        ds = []
        for a in A:
            try:
                lengths = nx.single_source_shortest_path_length(self.G, a)
            except Exception:
                continue
            dd = [lengths[b] for b in B if b in lengths]
            if dd:
                ds.append(min(dd))
        return statistics.mean(ds) if ds else None

    def _degree_matched_sample(self, ref, rng):
        out = []
        for g in ref:
            d = self.G.degree(g) if g in self.G else None
            pool = self._deg_bins.get(d) or list(self.G.nodes())
            out.append(rng.choice(pool))
        return out

    def proximity(self, drug_targets, disease_module, n_rand=200, seed=42):
        """z-score of drug->disease closeness vs degree-preserving null.
        Negative z = closer than chance = proximal to the disease module."""
        rng = random.Random(seed)
        d = self._mean_min_dist(drug_targets, disease_module)
        if d is None:
            return None
        null = []
        for _ in range(n_rand):
            rd = self._degree_matched_sample(drug_targets, rng)
            v = self._mean_min_dist(rd, disease_module)
            if v is not None:
                null.append(v)
        if not null:
            return None
        z = zscore(d, statistics.mean(null), statistics.pstdev(null) or 1e-9)
        return {"d": round(d, 3), "z": round(z, 2), "proximal": z < -0.5}

    def separation(self, A, B):
        """s_AB = <d_AB> - (<d_AA>+<d_BB>)/2. >0 => topologically separated."""
        dAB = self._mean_min_dist(A, B)
        dAA = self._mean_min_dist(A, A)
        dBB = self._mean_min_dist(B, B)
        if None in (dAB, dAA, dBB):
            return None
        return round(dAB - (dAA + dBB) / 2, 3)


def axis_coverage(target_set):
    return {axis: sorted(set(target_set) & set(genes)) for axis, genes in AXES.items()}


class EfficacyEngine:
    def __init__(self, knowledge):
        self.k = knowledge

    def assess(self, formula_targets, disease_scores, codrug_targets=None,
               dis_min=0.1, use_network=True):
        disease_module = {g for g, s in disease_scores.items() if s >= dis_min}
        inter = sorted(formula_targets & disease_module)

        result = {
            "n_formula_targets": len(formula_targets),
            "n_disease_module": len(disease_module),
            "shared_targets": inter,
            "formula_axis_coverage": {a: gs for a, gs in axis_coverage(formula_targets).items() if gs},
        }

        if codrug_targets is not None:
            cd = set(codrug_targets)
            result["codrug_targets"] = sorted(cd)
            result["codrug_axis_coverage"] = {a: gs for a, gs in axis_coverage(cd).items() if gs}
            overlap = formula_targets & cd
            union = formula_targets | cd
            result["target_overlap"] = sorted(overlap)
            result["jaccard"] = round(len(overlap) / len(union), 4) if union else 0.0
            result["combined_axis_coverage"] = {
                a: gs for a, gs in axis_coverage(union).items() if gs}
            result["formula_adds_axes"] = [
                a for a in AXES if (set(formula_targets) & set(AXES[a]))
                and not (cd & set(AXES[a]))]

        if use_network:
            genes = set(formula_targets) | set(disease_module)
            if codrug_targets:
                genes |= set(codrug_targets)
            try:
                iap = InteractomeProximity.build(genes)
                prox_f = iap.proximity(formula_targets, disease_module)
                result["network_proximity_formula"] = prox_f
                result["network_note"] = ("proximity z uses a STRING-induced subgraph "
                    "(module-biased null) — interpret cautiously; the robust "
                    "complementarity signals are separation and module engagement.")
                if codrug_targets:
                    prox_c = iap.proximity(set(codrug_targets), disease_module)
                    sep = iap.separation(formula_targets, set(codrug_targets))
                    result["network_proximity_codrug"] = prox_c
                    result["separation_formula_codrug"] = sep
                    # Complementary Exposure (Cheng et al.): BOTH engage the disease
                    # module AND the two target sets are topologically SEPARATED.
                    # "engages module" = shares >=1 target with the module OR proximal
                    # (robust to the induced-subgraph proximity artifact).
                    f_engages = bool(formula_targets & disease_module) or (prox_f and prox_f["proximal"])
                    c_engages = bool(set(codrug_targets) & disease_module) or (prox_c and prox_c["proximal"])
                    separated = (sep is not None and sep > 0) or result.get("jaccard", 1) < 0.05
                    result["complementary_exposure"] = bool(f_engages and c_engages and separated)
                    result["complementary_exposure_basis"] = {
                        "formula_engages_module": bool(f_engages),
                        "codrug_engages_module": bool(c_engages),
                        "topologically_separated": bool(separated),
                        "separation": sep, "jaccard": result.get("jaccard")}
            except Exception as e:
                result["network_proximity_error"] = str(e)
        return result
