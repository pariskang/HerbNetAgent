"""Tests for the Boltz-2 structure-based backend that need NO GPU:
YAML generation, affinity-JSON parsing, unit conversion, and graceful
fall-through when the boltz binary is absent (the CI/default case)."""
import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from herbnetagent.binding import BindingEngine, ChEMBLBackend
from herbnetagent.boltz2 import (
    Boltz2Backend,
    boltz_value_to_nM,
    boltz_value_to_pic50,
    build_yaml,
    parse_affinity,
)
from herbnetagent.knowledge import Knowledge


def test_yaml_has_protein_ligand_affinity():
    y = build_yaml("MADEUPSEQ", "O=C(O)c1ccccc1")
    assert "protein:" in y and "ligand:" in y
    assert "smiles: \"O=C(O)c1ccccc1\"" in y
    assert "affinity:" in y and "binder: B" in y


def test_affinity_unit_conversion():
    # affinity_pred_value = log10(IC50 in uM); lower = stronger
    # value 0 -> IC50 1 uM = 1000 nM ; pIC50 6
    assert abs(boltz_value_to_nM(0.0) - 1000.0) < 1e-6
    assert abs(boltz_value_to_pic50(0.0) - 6.0) < 1e-6
    # value -2 -> IC50 0.01 uM = 10 nM ; pIC50 8
    assert abs(boltz_value_to_nM(-2.0) - 10.0) < 1e-6
    assert abs(boltz_value_to_pic50(-2.0) - 8.0) < 1e-6


def test_parse_affinity_from_fixture():
    with tempfile.TemporaryDirectory() as td:
        sub = os.path.join(td, "boltz_results_job", "predictions", "job")
        os.makedirs(sub)
        with open(os.path.join(sub, "affinity_job.json"), "w") as f:
            json.dump({"affinity_pred_value": -1.5,
                       "affinity_probability_binary": 0.87}, f)
        res = parse_affinity(td)
    assert res is not None
    assert res["pIC50"] == 7.5
    assert abs(res["IC50_nM"] - round(boltz_value_to_nM(-1.5), 2)) < 1e-6
    assert res["affinity_probability_binary"] == 0.87


def test_backend_falls_through_without_boltz():
    # No boltz binary in CI -> available() False -> affinity() returns None,
    # so BindingEngine uses the measured ChEMBL value instead.
    k = Knowledge()
    b2 = Boltz2Backend(knowledge=k, boltz_bin="boltz-does-not-exist")
    assert b2.available() is False
    assert b2.affinity("quercetin", "F2") is None
    engine = BindingEngine([b2, ChEMBLBackend(k.compound_targets())])
    q = engine.affinity("quercetin", "F2")
    assert q is not None and q.source == "ChEMBL"   # graceful fallback


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} tests passed")
