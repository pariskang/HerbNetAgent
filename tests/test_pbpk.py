"""L5 PBPK engine: PK-Sim backend is unavailable here (no ospsuite/.pkml), so the
engine must fall back to the validated reduced scipy model and reproduce ~1.03x."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from herbnetagent.pbpk import PBPKEngine, PKSimBackend, ReducedPBPKBackend


def test_pksim_unavailable_without_ospsuite():
    assert PKSimBackend().available() is False   # no ospsuite / model in CI


def test_reduced_backend_always_available():
    assert ReducedPBPKBackend().available() is True


def test_engine_falls_back_and_reproduces_aucr():
    eng = PBPKEngine()
    assert eng.active_backend() == "reduced-PBPK(scipy)"
    res = eng.aucr("rivaroxaban", "quercetin", 0.006)
    assert res is not None
    assert res["backend"] == "reduced-PBPK(scipy)"
    assert 1.0 < res["aucr"].value < 1.10                 # ~1.035x
    assert 2.2 < res["validation_strong_inhibitor"] < 2.8  # ketoconazole calibration


def test_engine_unknown_victim_returns_none():
    eng = PBPKEngine()
    assert eng.aucr("unknown-drug", "quercetin", 1.0) is None


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} tests passed")
