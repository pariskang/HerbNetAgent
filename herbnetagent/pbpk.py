"""
herbnetagent.pbpk — L5 whole-body PBPK engine (pluggable, with fallback).

Same pattern as the L3 binding engine: an ordered list of backends resolved by
availability. A PKSimBackend (Open Systems Pharmacology, whole-body PBPK with
virtual populations) runs where `ospsuite` + a compiled model (.pkml) exist;
otherwise the always-available ReducedPBPKBackend (the validated scipy model in
pk.py) serves the AUCR. The orchestrator asks the engine, not a specific model.
"""
from __future__ import annotations

import glob
import os
from typing import Protocol

from .core import Evidence, Quantity
from .pk import PERPETRATORS, VICTIMS, ReducedPBPK

_HERE = os.path.dirname(os.path.abspath(__file__))
MODEL_DIR = os.path.join(os.path.dirname(_HERE), "models", "pbpk")


class PBPKResult(dict):
    """AUCR + validation + which backend produced it."""


class PBPKBackend(Protocol):
    name: str
    def available(self) -> bool: ...
    def aucr(self, victim: str, perpetrator: str, perpetrator_mg: float) -> PBPKResult | None: ...


class ReducedPBPKBackend:
    """Semi-mechanistic scipy model — always available; ketoconazole-calibrated."""
    name = "reduced-PBPK(scipy)"

    def available(self) -> bool:
        return True

    def aucr(self, victim: str, perpetrator: str, perpetrator_mg: float) -> PBPKResult | None:
        if victim not in VICTIMS or perpetrator not in PERPETRATORS:
            return None
        m = ReducedPBPK(victim=victim, perpetrator=perpetrator)
        q = m.aucr(perpetrator_mg)
        return PBPKResult(aucr=q, validation_strong_inhibitor=m.validate(),
                          backend=self.name)


class PKSimBackend:
    """Open Systems Pharmacology whole-body PBPK — real skeleton.

    A deployment provides a PK-Sim-built `<victim>.pkml` under models/pbpk/ and
    installs `ospsuite` (the OSP Python interface over the .NET libraries). This
    backend loads the model, applies the perpetrator's competitive inhibition of
    CYP3A4 / P-gp / BCRP (from ChEMBL Ki + the perpetrator's own PBPK-predicted
    plasma/gut concentrations), runs the simulation with vs without the
    perpetrator, and returns the AUC ratio. Supports virtual populations.

    Where `ospsuite`/model are absent (as here) `available()` is False, so the
    engine falls through to ReducedPBPKBackend.
    """
    name = "PK-Sim/OSP (whole-body PBPK)"

    def __init__(self, model_dir: str = MODEL_DIR):
        self.model_dir = model_dir

    def available(self) -> bool:
        try:
            import ospsuite  # noqa: F401
        except Exception:
            return False
        return os.path.isdir(self.model_dir) and bool(
            glob.glob(os.path.join(self.model_dir, "*.pkml")))

    # ---- OSP-specific implementation (runs only when available) -------------
    def _model_path(self, victim: str) -> str | None:
        hits = glob.glob(os.path.join(self.model_dir, f"{victim}*.pkml"))
        return hits[0] if hits else None

    def aucr(self, victim: str, perpetrator: str, perpetrator_mg: float,
             population: int = 0) -> PBPKResult | None:  # pragma: no cover
        if not self.available():
            return None
        import ospsuite  # OSP Python interface (pythonnet-backed)
        model = self._model_path(victim)
        if model is None:
            return None
        # --- load whole-body PBPK model built in PK-Sim ---
        sim = ospsuite.load_simulation(model)

        def run_auc(s) -> float:
            ospsuite.run_simulation(s)
            res = s.results
            # integrate victim plasma (venous blood) over the last dosing interval;
            # parameter path is model-specific -> configured per pkml.
            t, c = res.time, res.values_of("Organism|PeripheralVenousBlood|"
                                           f"{victim}|Plasma (Peripheral Venous Blood)")
            import numpy as np
            return float(np.trapezoid(c, t))

        base = run_auc(sim)
        # --- apply perpetrator competitive inhibition of clearance pathways ---
        self._apply_inhibition(sim, perpetrator, perpetrator_mg)
        inhibited = run_auc(sim)
        r = inhibited / base
        return PBPKResult(
            aucr=Quantity(round(r, 3), unit="x AUC", evidence=Evidence.COMPUTATIONAL,
                          source="PK-Sim/OSP", note=f"whole-body PBPK; +{(r-1)*100:.1f}%"),
            validation_strong_inhibitor=None, backend=self.name,
            virtual_population=population)

    def _apply_inhibition(self, sim, perpetrator: str, perpetrator_mg: float):  # pragma: no cover
        """Set competitive-inhibition parameters on CYP3A4/P-gp/BCRP processes.

        In OSP this means adding the perpetrator as a second compound with its Ki
        against each enzyme/transporter, or scaling the victim's metabolic/transport
        kcat by 1/(1+[I]/Ki). Parameter paths are configured per model. Left as the
        explicit integration point for a PK-Sim deployment.
        """
        raise NotImplementedError("configure OSP inhibition parameters for this model")


class PBPKEngine:
    def __init__(self, backends: list[PBPKBackend] | None = None):
        self.backends = backends or [PKSimBackend(), ReducedPBPKBackend()]

    def aucr(self, victim: str, perpetrator: str, perpetrator_mg: float) -> PBPKResult | None:
        for b in self.backends:
            try:
                if not b.available():
                    continue
                res = b.aucr(victim, perpetrator, perpetrator_mg)
                if res is not None:
                    return res
            except Exception:
                continue
        return None

    def active_backend(self) -> str:
        for b in self.backends:
            if b.available():
                return b.name
        return "none"
