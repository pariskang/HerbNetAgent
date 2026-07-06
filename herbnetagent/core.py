"""
herbnetagent.core — L0 primitives shared by every layer.

Three ideas are first-class citizens here (the design principles from
docs/AGENT_ARCHITECTURE.md): evidence grading, uncertainty, and provenance.
Nothing in the agent returns a bare number: it returns a Quantity that knows
how confident it is and where it came from.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import IntEnum
from typing import Any


class Evidence(IntEnum):
    """Evidence strength ladder — computation is the weakest rung.

    The whole point of the project's honesty policy: a prediction is not proof.
    Reports must never state a COMPUTATIONAL result with the certainty of a TRIAL.
    """
    COMPUTATIONAL = 1   # network pharmacology / docking / PBPK prediction
    IN_VITRO = 2        # measured assay (e.g. ChEMBL pChEMBL)
    IN_VIVO = 3         # animal model
    HUMAN_PK = 4        # clinical PK / DDI study
    RCT = 5             # randomised controlled trial

    @property
    def label(self) -> str:
        return {1: "计算预测", 2: "体外实测", 3: "动物", 4: "人体PK", 5: "RCT"}[int(self)]


@dataclass
class Quantity:
    """A value that carries its own uncertainty, evidence level and provenance."""
    value: Any
    unit: str = ""
    ci: tuple[float, float] | None = None       # (low, high)
    evidence: Evidence = Evidence.COMPUTATIONAL
    source: str = ""                            # database / model + version
    note: str = ""

    def graded(self) -> str:
        v = f"{self.value}{(' ' + self.unit) if self.unit else ''}"
        if self.ci:
            v += f" [{self.ci[0]:.3g}, {self.ci[1]:.3g}]"
        return f"{v}  ⟨{self.evidence.label}·{self.source or '—'}⟩"

    def to_dict(self) -> dict:
        d = asdict(self)
        d["evidence"] = int(self.evidence)
        d["evidence_label"] = self.evidence.label
        return d


@dataclass
class ProvenanceLog:
    """Append-only record of every engine call — makes runs auditable/reproducible."""
    records: list[dict] = field(default_factory=list)

    def record(self, layer: str, tool: str, inputs: dict, summary: str,
               evidence: Evidence = Evidence.COMPUTATIONAL) -> None:
        self.records.append({
            "layer": layer, "tool": tool, "inputs": inputs,
            "summary": summary, "evidence": evidence.label,
        })

    def to_dict(self) -> list[dict]:
        return self.records


class Registry:
    """Minimal tool/engine registry (the L0 tool-router substrate).

    Engines register under a capability name; the orchestrator resolves the best
    available backend (e.g. binding -> ChEMBL now, Boltz-2 when a GPU appears).
    """
    def __init__(self) -> None:
        self._tools: dict[str, Any] = {}

    def register(self, name: str, impl: Any, priority: int = 0) -> None:
        cur = self._tools.get(name)
        if cur is None or priority >= cur[1]:
            self._tools[name] = (impl, priority)

    def get(self, name: str) -> Any:
        if name not in self._tools:
            raise KeyError(f"no engine registered for capability '{name}'")
        return self._tools[name][0]

    def has(self, name: str) -> bool:
        return name in self._tools

    def capabilities(self) -> list[str]:
        return sorted(self._tools)


def zscore(value: float, null_mean: float, null_sd: float) -> float:
    if null_sd == 0:
        return 0.0
    return (value - null_mean) / null_sd


def pchembl_to_nM(pchembl: float) -> float:
    """pChEMBL (-log10 M) -> nM."""
    return 10 ** (-pchembl) * 1e9
