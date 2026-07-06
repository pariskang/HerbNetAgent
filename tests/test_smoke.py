"""Smoke + regression tests: the packaged agent must reproduce the manual
01-19 pipeline's headline numbers (offline, from the cached knowledge)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from herbnetagent import HerbNetAgent, Knowledge, Evidence
from herbnetagent.pk import ReducedPBPK
from herbnetagent.binding import BindingEngine, ChEMBLBackend


def test_knowledge_loads():
    k = Knowledge()
    assert len(k.formula_targets()) == 208
    assert k.codrug("rivaroxaban")["targets"]["F10"] == 9.4


def test_efficacy_shared_targets():
    agent = HerbNetAgent()
    res = agent.assess(codrug="rivaroxaban", use_network=False)
    eff = res["efficacy"]
    assert len(eff["shared_targets"]) == 16
    assert "F2" in eff["shared_targets"]
    # rivaroxaban (FXa) is complementary: overlap only F2
    assert eff["jaccard"] < 0.02
    assert set(eff["formula_adds_axes"]) == {"Platelet", "Thrombo-inflammation",
                                             "Resolution", "Endothelium"}


def test_binding_affinity_f2_quercetin():
    k = Knowledge()
    be = BindingEngine([ChEMBLBackend(k.compound_targets())])
    q = be.affinity("quercetin", "F2")
    assert q is not None
    assert q.evidence == Evidence.IN_VITRO
    assert 30 < q.value < 50          # ~39 nM, matches step 02


def test_pbpk_aucr_and_validation():
    pbpk = ReducedPBPK("rivaroxaban", "quercetin")
    assert 2.2 < pbpk.validate() < 2.8            # ketoconazole calibration ~2.5x
    aucr = pbpk.aucr(0.006)                        # 1.5 g pill quercetin
    assert 1.0 < aucr.value < 1.10                 # ~1.03x, negligible


def test_safety_verdict():
    agent = HerbNetAgent()
    res = agent.assess(codrug="rivaroxaban", codrug_perpetrator_mg=0.006, use_network=False)
    s = res["safety"]
    assert s["pk_interaction_level"] == "negligible"
    assert s["pd_bleeding_additivity"] is True
    assert res["disclaimer"]


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn(); print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} tests passed")
