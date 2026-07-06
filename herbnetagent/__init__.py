"""
HerbNetAgent 2.0 — Quantitative Multi-Scale TCM Intelligence Agent.

A runnable implementation of the 8-layer architecture in docs/AGENT_ARCHITECTURE.md:
molecule → structure → affinity → network → PBPK exposure → QSP effect → safety →
evidence-graded report. The manual 01-19 pipeline is refactored here into a
reusable, provenance-tracking, uncertainty-aware library.
"""
from .core import Evidence, ProvenanceLog, Quantity, Registry
from .knowledge import Knowledge
from .orchestrator import HerbNetAgent

__version__ = "2.0.0"
__all__ = ["HerbNetAgent", "Knowledge", "Quantity", "Evidence",
           "ProvenanceLog", "Registry"]
