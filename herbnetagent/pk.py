"""
herbnetagent.pk — L5/L6 ADMET / PK / systems-PD engine.

A reduced (semi-mechanistic) dynamic PBPK model: a perpetrator (real human PK)
dosed to steady state dynamically modulates a victim drug's clearance and gut
availability through hepatic CYP3A4, renal P-gp/BCRP secretion and intestinal
CYP3A4 + efflux — yielding the AUC ratio (AUCR). Calibrated against the known
ketoconazole interaction (clinical AUCR ~2.6x). This is the same model validated
in scripts/12; here it is a reusable class the orchestrator calls for any
perpetrator/victim pair.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from .core import Evidence, Quantity

# victim library: fractional clearance routes + PK, for AUCR prediction
VICTIMS = {
    "rivaroxaban": dict(
        MW=435.9, dose_mg=10.0, ka=1.24, V=50.0, fu=0.075, thalf=6.0,
        Fa0=0.95, Fg0=0.95,
        f_CYP3A4=0.30, f_CYP2J2=0.14, f_hydrolysis=0.10,
        f_renal_secretion=0.30, f_other=0.16),
}

# perpetrator library: human PK + measured inhibition constants (uM)
PERPETRATORS = {
    "quercetin": dict(MW=302.23, ka=0.6, V=200.0, F=0.02, fu=0.015, CL=70.0,
                      Ki_CYP3A4=1.0, Ki_pgp=7.08, Ki_bcrp=0.03),
}
GUT_VOL = 0.25
QH = 90.0


class ReducedPBPK:
    """Dynamic AUCR predictor for a (victim, perpetrator) pair."""

    def __init__(self, victim: str = "rivaroxaban", perpetrator: str = "quercetin"):
        self.riv = dict(VICTIMS[victim])
        self.per = dict(PERPETRATORS[perpetrator])
        self.riv["CL"] = 0.693 * self.riv["V"] / self.riv["thalf"]
        self.riv["dose_umol"] = self.riv["dose_mg"] / self.riv["MW"] * 1000
        self.victim, self.perpetrator = victim, perpetrator

    def _simulate(self, perp_mg, days=20, strong=None):
        q = perp_mg / self.per["MW"] * 1000 if perp_mg else 0.0
        R, P = self.riv, self.per

        def rhs(t, y):
            Aqg, Aqc, Arg, Arc = y
            Cq = Aqc / P["V"]
            Iq_sys = P["fu"] * Cq
            Iq_h = P["fu"] * (Cq + (P["F"] * P["ka"] * Aqg) / QH)
            Iq_gut = Aqg / GUT_VOL
            if strong:
                A_h, A_r, A_g, A_e = strong
            else:
                A_h = 1 / (1 + Iq_h / P["Ki_CYP3A4"])
                A_r = 1 / (1 + Iq_sys / P["Ki_pgp"] + Iq_sys / P["Ki_bcrp"])
                A_g = 1 / (1 + Iq_gut / P["Ki_CYP3A4"])
                A_e = 1 / (1 + Iq_gut / P["Ki_bcrp"] + Iq_gut / P["Ki_pgp"])
            CL = R["CL"] * (R["f_CYP3A4"] * A_h + R["f_CYP2J2"] + R["f_hydrolysis"]
                            + R["f_renal_secretion"] * A_r + R["f_other"])
            Fg = R["Fg0"] + (1 - R["Fg0"]) * (1 - A_g)
            Fa = R["Fa0"] + (1 - R["Fa0"]) * (1 - A_e)
            return [-P["ka"] * Aqg, P["ka"] * Aqg * P["F"] - (P["CL"] / P["V"]) * Aqc,
                    -R["ka"] * Arg, R["ka"] * Arg * Fa * Fg - (CL / R["V"]) * Arc]

        y = [0, 0, 0, 0]
        auc = None
        for d in range(days):
            y[0] += q
            y[2] += R["dose_umol"]
            sol = solve_ivp(rhs, [0, 24], y, max_step=0.25, rtol=1e-7, atol=1e-9,
                            dense_output=True)
            y = [sol.y[i, -1] for i in range(4)]
            if d == days - 1:
                tt = np.linspace(0, 24, 241)
                auc = float(np.trapezoid(sol.sol(tt)[3] / R["V"], tt))
        return auc

    def validate(self):
        base = self._simulate(0)
        strong = self._simulate(0, strong=(0.03, 0.12, 0.05, 0.05))
        return round(strong / base, 2)

    def aucr(self, perpetrator_mg: float) -> Quantity:
        base = self._simulate(0)
        auc = self._simulate(perpetrator_mg)
        r = auc / base
        return Quantity(value=round(r, 3), unit="x AUC",
                        evidence=Evidence.COMPUTATIONAL,
                        source=f"reduced-PBPK({self.victim}/{self.perpetrator})",
                        note=f"perp {perpetrator_mg:.4g} mg; +{(r-1)*100:.1f}% exposure")
