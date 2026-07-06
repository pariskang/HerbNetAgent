"""
herbnetagent.safety — L7 safety, DDI & clinical-translation guardrails.

Turns the mechanistic findings into an evidence-graded risk read-out and, per the
project's honesty policy, never upgrades a computational prediction to a clinical
claim. Also flags pharmacodynamic (bleeding) additivity that PK numbers miss.
"""
from __future__ import annotations
from .core import Quantity, Evidence

# formula constituents with known bleeding-relevant pharmacology not captured by
# the small-molecule PK model (the peptide anticoagulant is the key blind spot)
PD_BLEEDING_FLAGS = {
    "hirudin (from 水蛭/leech)": "direct thrombin (F2) inhibitor — adds to any anticoagulant",
}


class SafetyEngine:
    def __init__(self, knowledge):
        self.k = knowledge

    def ddi(self, codrug: dict, aucr: Quantity | None, shared_pd_axes: list[str]):
        """Assemble a graded DDI verdict for formula + co-drug."""
        flags = []

        # PK axis
        pk_level = "negligible"
        if aucr is not None:
            r = aucr.value
            pk_level = ("negligible" if r < 1.25 else
                        "moderate" if r < 2 else "major")
            flags.append(Quantity(
                value=f"victim exposure {r:.2f}x ({pk_level})", unit="",
                evidence=Evidence.COMPUTATIONAL, source=aucr.source,
                note="PK interaction via CYP3A4/P-gp/BCRP"))

        # PD axis — bleeding additivity (the honest, dominant concern)
        pd_axes = [a for a in shared_pd_axes if a in ("Coagulation", "Platelet")]
        pd_concern = bool(pd_axes) or codrug.get("class", "").lower().find("inhibitor") >= 0
        if pd_concern:
            flags.append(Quantity(
                value="additive bleeding risk (PD)", unit="",
                evidence=Evidence.COMPUTATIONAL, source="mechanism-coverage",
                note=f"overlapping haemostatic axes: {pd_axes or ['anticoagulant co-drug']}; "
                     f"plus {', '.join(PD_BLEEDING_FLAGS)}"))

        verdict = {
            "pk_interaction_level": pk_level,
            "pd_bleeding_additivity": pd_concern,
            "flags": [f.graded() for f in flags],
            "bottom_line": self._bottom_line(pk_level, pd_concern),
            "evidence_ceiling": Evidence.COMPUTATIONAL.label,
        }
        return verdict

    @staticmethod
    def _bottom_line(pk_level, pd_concern):
        if pk_level == "negligible" and pd_concern:
            return ("PK 相互作用可忽略；主要残余风险为 PD 出血叠加 —— 需临床监护，"
                    "不作为自我用药依据（证据仅达计算级）")
        if pk_level in ("moderate", "major"):
            return f"PK 暴露升高({pk_level})叠加出血风险 —— 应避免或严密监护"
        return "机制层面未见显著相互作用信号；仍须临床个体化判断"

    def evidence_disclaimer(self):
        return ("本智能体所有输出的最高证据等级为 **计算预测**；"
                "临床用药决策须由医师依据体外/动物/RCT 证据个体化判断。")
