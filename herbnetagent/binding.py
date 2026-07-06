"""
herbnetagent.binding — L3 binding & target-engagement engine.

Pluggable backends behind one interface. Today the ChEMBL backend serves
*measured* affinities (the strongest evidence we can get without a lab). The
Boltz2Backend is a stub with the exact call signature a GPU deployment would
fill in — so the rest of the agent is written against the interface, not the tool.
"""
from __future__ import annotations
from typing import Protocol
from .core import Quantity, Evidence, pchembl_to_nM


class BindingBackend(Protocol):
    name: str
    def affinity(self, compound: str, gene: str) -> Quantity | None: ...


class ChEMBLBackend:
    """Serve measured pChEMBL affinities from the cached compound->target table."""
    name = "ChEMBL(measured)"
    evidence = Evidence.IN_VITRO

    def __init__(self, compound_targets: list[dict]):
        self._best: dict[tuple[str, str], float] = {}
        for e in compound_targets:
            k = (e["compound"], e["gene"])
            pv = float(e["max_pchembl"])
            self._best[k] = max(self._best.get(k, 0.0), pv)

    def affinity(self, compound: str, gene: str) -> Quantity | None:
        pv = self._best.get((compound, gene))
        if pv is None:
            return None
        return Quantity(value=round(pchembl_to_nM(pv), 1), unit="nM (Ki/IC50-eq)",
                        evidence=self.evidence, source="ChEMBL",
                        note=f"pChEMBL={pv:.2f}")


# Real structure-based backend (Boltz-2 / AlphaFold3) lives in boltz2.py; it runs
# where the `boltz` binary + GPU exist and otherwise falls through to ChEMBL.
from .boltz2 import Boltz2Backend, AlphaFold3Backend  # noqa: E402,F401


class BindingEngine:
    """Resolve affinity through an ordered list of backends (best evidence first)."""

    def __init__(self, backends: list[BindingBackend]):
        self.backends = backends

    def affinity(self, compound: str, gene: str) -> Quantity | None:
        for b in self.backends:
            q = b.affinity(compound, gene)
            if q is not None:
                return q
        return None

    def target_profile(self, compound: str, genes) -> dict[str, Quantity]:
        out = {}
        for g in genes:
            q = self.affinity(compound, g)
            if q is not None:
                out[g] = q
        return out
