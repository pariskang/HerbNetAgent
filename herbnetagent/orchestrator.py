"""
herbnetagent.orchestrator — L0 driver.

Wires the layers into one end-to-end assessment for a (formula, disease, [co-drug],
[dose]) request, tracks provenance, and emits a structured result + graded report.
This is the programmatic form of the manual 01-19 pipeline.
"""
from __future__ import annotations
import json
from .core import Registry, ProvenanceLog, Evidence, Quantity
from .knowledge import Knowledge
from .binding import BindingEngine, ChEMBLBackend, Boltz2Backend
from .efficacy import EfficacyEngine, AXES
from .pk import ReducedPBPK, VICTIMS
from .safety import SafetyEngine


class HerbNetAgent:
    def __init__(self, knowledge: Knowledge | None = None):
        self.k = knowledge or Knowledge()
        self.prov = ProvenanceLog()
        self.reg = Registry()
        # L3 binding: measured first, predicted (Boltz-2) as future backend
        self.reg.register("binding", BindingEngine(
            [Boltz2Backend(), ChEMBLBackend(self.k.compound_targets())]), priority=1)
        self.reg.register("efficacy", EfficacyEngine(self.k), priority=1)
        self.reg.register("safety", SafetyEngine(self.k), priority=1)

    # ------------------------------------------------------------------ run
    def assess(self, disease: str = "venous thrombosis", codrug: str | None = None,
               codrug_perpetrator_mg: float | None = None,
               perpetrator: str = "quercetin", use_network: bool = True) -> dict:
        out: dict = {"request": {"disease": disease, "codrug": codrug,
                                 "perpetrator_mg": codrug_perpetrator_mg}}

        # L1/L2 knowledge
        formula_targets = self.k.formula_targets()
        disease_scores = self.k.disease_targets()
        self.prov.record("L1", "Knowledge", {"disease": disease},
                         f"{len(formula_targets)} formula targets, "
                         f"{len(disease_scores)} disease genes", Evidence.IN_VITRO)

        cd = self.k.codrug(codrug) if codrug else None
        cd_targets = set(cd["targets"]) if cd else None

        # L4 efficacy
        eff_engine = self.reg.get("efficacy")
        eff = eff_engine.assess(formula_targets, disease_scores, cd_targets,
                                use_network=use_network)
        out["efficacy"] = eff
        self.prov.record("L4", "EfficacyEngine",
                         {"codrug": codrug},
                         f"{len(eff['shared_targets'])} shared targets; "
                         f"complementary_exposure={eff.get('complementary_exposure')}",
                         Evidence.COMPUTATIONAL)

        # L3 binding — annotate the top shared targets with measured affinity
        binder = self.reg.get("binding")
        top = eff["shared_targets"][:8]
        aff = {}
        for g in top:
            # find best compound for this gene from the cached edges
            best = None
            for e in self.k.compound_targets():
                if e["gene"] == g:
                    q = binder.affinity(e["compound"], g)
                    if q and (best is None or q.value < best[1].value):
                        best = (e["compound"], q)
            if best:
                aff[g] = {"compound": best[0], "affinity": best[1].graded()}
        out["binding_highlights"] = aff
        self.prov.record("L3", "BindingEngine", {"targets": top},
                         f"annotated {len(aff)} targets (backend: {binder.backends[-1].name})",
                         Evidence.IN_VITRO)

        # L5/L6 PK — AUCR if a co-drug + perpetrator dose is supplied
        aucr = None
        if codrug and codrug_perpetrator_mg is not None and codrug in VICTIMS:
            pbpk = ReducedPBPK(victim=codrug, perpetrator=perpetrator)
            val_aucr = pbpk.validate()
            aucr = pbpk.aucr(codrug_perpetrator_mg)
            out["pk"] = {"aucr": aucr.to_dict(),
                         "model_validation_strong_inhibitor_aucr": val_aucr,
                         "graded": aucr.graded()}
            self.prov.record("L5", "ReducedPBPK",
                             {"victim": codrug, "perp_mg": codrug_perpetrator_mg},
                             f"AUCR={aucr.value}x (validated {val_aucr}x)",
                             Evidence.COMPUTATIONAL)

        # L7 safety
        safety = self.reg.get("safety")
        shared_axes = [a for a in AXES if (formula_targets & set(AXES[a]))
                       and (cd_targets and cd_targets & set(AXES[a]))] if cd_targets else []
        if cd:
            out["safety"] = safety.ddi(cd, aucr, shared_axes)
            self.prov.record("L7", "SafetyEngine", {"codrug": codrug},
                             out["safety"]["bottom_line"], Evidence.COMPUTATIONAL)
        out["disclaimer"] = safety.evidence_disclaimer()
        out["provenance"] = self.prov.to_dict()
        return out

    # --------------------------------------------------------------- report
    @staticmethod
    def render_markdown(res: dict) -> str:
        L = []
        w = L.append
        req = res["request"]
        w(f"# HerbNetAgent 2.0 — 评估结果")
        w("")
        w(f"> 疾病: **{req['disease']}**"
          + (f" · 联用: **{req['codrug']}**" if req["codrug"] else "")
          + (f" · 摄动剂量: {req['perpetrator_mg']} mg" if req["perpetrator_mg"] is not None else ""))
        w("")
        eff = res["efficacy"]
        w("## L4 药效 (Efficacy)")
        w(f"- 方剂靶点 {eff['n_formula_targets']}，疾病模块 {eff['n_disease_module']}，"
          f"共享机制靶点 **{len(eff['shared_targets'])}**：{', '.join(eff['shared_targets'][:15])}")
        w(f"- 方剂覆盖病理轴：{', '.join(a for a in eff['formula_axis_coverage'])}")
        if "complementary_exposure" in eff:
            w(f"- 与 {req['codrug']} 的靶点重叠 Jaccard={eff['jaccard']}；"
              f"方剂独家补充轴：{', '.join(eff.get('formula_adds_axes', []))}")
            basis = eff.get("complementary_exposure_basis", {})
            if basis:
                w(f"- 网络拓扑：分离度 s={basis.get('separation')}（>0 表示两靶点集拓扑分离）、"
                  f"Jaccard={basis.get('jaccard')}；双方均engage疾病模块="
                  f"{basis.get('formula_engages_module') and basis.get('codrug_engages_module')}")
            w(f"- **互补暴露 (Cheng/Barabási): {eff.get('complementary_exposure')}** "
              f"（判据：双方均engage疾病模块 且 两靶点集拓扑分离 → 有效药物组合的网络特征）")
            if eff.get("network_note"):
                w(f"  - _注：{eff['network_note']}_")
        if res.get("binding_highlights"):
            w("")
            w("### L3 结合亲和力（实测）")
            for g, d in res["binding_highlights"].items():
                w(f"- **{g}** ← {d['compound']}: {d['affinity']}")
        if "pk" in res:
            w("")
            w("## L5 药代 (PK / DDI)")
            w(f"- 动态 PBPK 预测 **AUCR = {res['pk']['aucr']['value']}× "
              f"({res['pk']['graded']})**")
            w(f"- 模型验证（强抑制剂）AUCR ≈ {res['pk']['model_validation_strong_inhibitor_aucr']}× "
              f"(临床酮康唑 ~2.6×)")
        if "safety" in res:
            s = res["safety"]
            w("")
            w("## L7 安全 (Safety / DDI)")
            w(f"- PK 相互作用：**{s['pk_interaction_level']}** · PD 出血叠加：**{s['pd_bleeding_additivity']}**")
            for f in s["flags"]:
                w(f"  - {f}")
            w(f"- **结论**：{s['bottom_line']}")
        w("")
        w(f"> ⚠️ {res['disclaimer']}")
        return "\n".join(L)
