"""
herbnetagent.structure — L2 structure provider.

Resolves a gene symbol to (UniProt accession, protein sequence, AlphaFold-DB
structure). Boltz-2/AlphaFold3 fold from sequence, so the sequence is the minimum
requirement; the AlphaFold PDB is optional (useful for docking/visualisation).

A curated accession registry covers the targets this project cares about (so the
skeleton is deterministic for F2/MMP9 etc.); anything else falls back to the
UniProt REST search. Everything is cached under data/structures/.
"""
from __future__ import annotations
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(os.path.dirname(_HERE), "data", "structures")
os.makedirs(CACHE, exist_ok=True)

try:
    import requests
    _UA = {"User-Agent": "HerbNetAgent/2.0 (research)"}
except Exception:
    requests = None
    _UA = {}

# curated gene -> reviewed human UniProt accession (Swiss-Prot)
UNIPROT = {
    "F2": "P00734", "F10": "P00742", "MMP9": "P14780", "MMP2": "P08253",
    "PTGS2": "P35354", "PTGS1": "P23219", "ALOX12": "P18054", "ALOX15": "P16050",
    "KDR": "P35968", "NLRP3": "Q96P20", "JAK2": "O60674", "NFKB1": "P19838",
    "ESR1": "P03372", "CTNNB1": "P35222", "PIK3R1": "P27986", "STAT6": "P42226",
}


class StructureProvider:
    def __init__(self, cache: str = CACHE):
        self.cache = cache

    # ---- gene -> UniProt accession ------------------------------------------
    def accession(self, gene: str) -> str | None:
        if gene in UNIPROT:
            return UNIPROT[gene]
        if requests is None:
            return None
        try:
            r = requests.get("https://rest.uniprot.org/uniprotkb/search", headers=_UA, timeout=30,
                             params={"query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
                                     "fields": "accession", "format": "tsv", "size": 1})
            lines = r.text.strip().splitlines()
            return lines[1] if len(lines) > 1 else None
        except Exception:
            return None

    # ---- accession -> sequence (cached FASTA) -------------------------------
    def sequence(self, gene: str) -> str | None:
        acc = self.accession(gene)
        if not acc:
            return None
        fp = os.path.join(self.cache, f"{acc}.fasta")
        if os.path.exists(fp):
            return self._read_fasta(fp)
        if requests is None:
            return None
        try:
            r = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", headers=_UA, timeout=30)
            if r.status_code == 200 and r.text.startswith(">"):
                with open(fp, "w") as f:
                    f.write(r.text)
                return self._read_fasta(fp)
        except Exception:
            pass
        return None

    # ---- accession -> AlphaFold-DB structure (optional) ---------------------
    def alphafold_pdb(self, gene: str) -> str | None:
        acc = self.accession(gene)
        if not acc:
            return None
        fp = os.path.join(self.cache, f"AF-{acc}.pdb")
        if os.path.exists(fp):
            return fp
        if requests is None:
            return None
        try:
            url = f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v4.pdb"
            r = requests.get(url, headers=_UA, timeout=60)
            if r.status_code == 200 and r.text.startswith(("HEADER", "ATOM", "MODEL")):
                with open(fp, "w") as f:
                    f.write(r.text)
                return fp
        except Exception:
            pass
        return None

    @staticmethod
    def _read_fasta(fp: str) -> str:
        with open(fp) as f:
            return "".join(l.strip() for l in f if not l.startswith(">"))
