"""
herbnetagent.boltz2 — L3 structure-based binding backend (Boltz-2 / AlphaFold3).

A *real* implementation skeleton: it resolves compound SMILES + target sequence,
writes a valid Boltz-2 YAML (protein chain + ligand + affinity property), runs the
`boltz predict` CLI, and parses `affinity_pred_value` / `affinity_probability_binary`
from the output. It runs for real wherever the `boltz` binary + a GPU are present;
where they are not (as here) `available()` is False and `affinity()` returns None,
so BindingEngine transparently falls back to the measured ChEMBL backend.

Boltz-2 affinity convention (verified against docs, jwohlwend/boltz):
  affinity_pred_value = log10(IC50) with IC50 in µM   (LOWER = stronger binding)
  => IC50[nM] = 10**(value + 3);  pIC50 = 6 - value
  affinity_probability_binary in [0,1] = P(ligand is a binder)
NOTE: verify the numeric convention against your installed Boltz-2 version before
using values quantitatively — it has changed across releases.
"""
from __future__ import annotations

import glob
import json
import os
import shutil
import subprocess
import tempfile

from .core import Evidence, Quantity
from .structure import StructureProvider

_HERE = os.path.dirname(os.path.abspath(__file__))
RESULT_CACHE = os.path.join(os.path.dirname(_HERE), "data", "structures", "boltz_cache.json")


def boltz_value_to_nM(pred_value: float) -> float:
    """affinity_pred_value (log10 IC50 in µM) -> IC50 in nM."""
    return 10 ** (pred_value + 3)


def boltz_value_to_pic50(pred_value: float) -> float:
    return 6.0 - pred_value


def build_yaml(sequence: str, smiles: str) -> str:
    """Boltz-2 input YAML: one protein chain + one ligand, request affinity."""
    return (
        "version: 1\n"
        "sequences:\n"
        "  - protein:\n"
        "      id: A\n"
        f"      sequence: \"{sequence}\"\n"
        "  - ligand:\n"
        "      id: B\n"
        f"      smiles: \"{smiles}\"\n"
        "properties:\n"
        "  - affinity:\n"
        "      binder: B\n"
    )


def parse_affinity(out_dir: str) -> dict | None:
    """Find and parse the affinity_*.json Boltz-2 writes under the output dir."""
    hits = glob.glob(os.path.join(out_dir, "**", "affinity_*.json"), recursive=True)
    if not hits:
        return None
    with open(hits[0]) as f:
        d = json.load(f)
    val = d.get("affinity_pred_value")
    if val is None:
        return None
    return {
        "affinity_pred_value": float(val),
        "affinity_probability_binary": d.get("affinity_probability_binary"),
        "IC50_nM": round(boltz_value_to_nM(float(val)), 2),
        "pIC50": round(boltz_value_to_pic50(float(val)), 2),
    }


class Boltz2Backend:
    """Co-folding + affinity backend. Falls through (returns None) when unavailable."""
    name = "Boltz-2(predicted)"
    evidence = Evidence.COMPUTATIONAL

    def __init__(self, knowledge=None, structure_provider: StructureProvider | None = None,
                 boltz_bin: str = "boltz", use_msa_server: bool = True,
                 extra_args: tuple = (), cache_path: str = RESULT_CACHE):
        self.k = knowledge
        self.sp = structure_provider or StructureProvider()
        self.boltz_bin = boltz_bin
        self.use_msa_server = use_msa_server
        self.extra_args = tuple(extra_args)
        self.cache_path = cache_path
        self._cache = self._load_cache()

    # ---- availability --------------------------------------------------------
    def available(self) -> bool:
        return shutil.which(self.boltz_bin) is not None

    # ---- SMILES resolution ---------------------------------------------------
    def _smiles(self, compound: str) -> str | None:
        if self.k is not None:
            for c in self.k.compounds():
                if c["compound"] == compound and c.get("smiles"):
                    return c["smiles"]
        return None

    # ---- main entry ----------------------------------------------------------
    def affinity(self, compound: str, gene: str) -> Quantity | None:
        key = f"{compound}|{gene}"
        if key in self._cache:
            return self._to_quantity(self._cache[key], compound, gene)
        if not self.available():
            return None  # no boltz binary here -> let ChEMBL backend serve
        smiles = self._smiles(compound)
        seq = self.sp.sequence(gene)
        if not smiles or not seq:
            return None
        res = self._run(seq, smiles)
        if res is None:
            return None
        self._cache[key] = res
        self._save_cache()
        return self._to_quantity(res, compound, gene)

    def _run(self, sequence: str, smiles: str) -> dict | None:
        with tempfile.TemporaryDirectory() as td:
            yml = os.path.join(td, "job.yaml")
            out = os.path.join(td, "out")
            with open(yml, "w") as f:
                f.write(build_yaml(sequence, smiles))
            cmd = [self.boltz_bin, "predict", yml, "--out_dir", out]
            if self.use_msa_server:
                cmd.append("--use_msa_server")
            cmd += list(self.extra_args)
            try:
                subprocess.run(cmd, check=True, capture_output=True, timeout=3600)
            except Exception:
                return None
            return parse_affinity(out)

    def _to_quantity(self, res: dict, compound: str, gene: str) -> Quantity:
        return Quantity(
            value=res["IC50_nM"], unit="nM (predicted IC50)",
            evidence=Evidence.COMPUTATIONAL, source="Boltz-2",
            note=f"{compound}->{gene}; pIC50={res.get('pIC50')}, "
                 f"P(binder)={res.get('affinity_probability_binary')}")

    # ---- cache ---------------------------------------------------------------
    def _load_cache(self) -> dict:
        try:
            with open(self.cache_path) as f:
                return json.load(f)
        except Exception:
            return {}

    def _save_cache(self) -> None:
        try:
            os.makedirs(os.path.dirname(self.cache_path), exist_ok=True)
            with open(self.cache_path, "w") as f:
                json.dump(self._cache, f, indent=1)
        except Exception:
            pass


class AlphaFold3Backend(Boltz2Backend):
    """Same interface via the AlphaFold3 server / local weights.

    AF3's public server has usage limits and a distinct I/O format; a deployment
    overrides `_run` to submit the AF3 job and read its ranking/affinity outputs.
    Left as a thin subclass so the wiring is explicit.
    """
    name = "AlphaFold3(predicted)"

    def _run(self, sequence: str, smiles: str) -> dict | None:  # pragma: no cover
        raise NotImplementedError("wire up AF3 server/local weights here")
