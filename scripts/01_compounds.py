#!/usr/bin/env python3
"""
Step 1: Compile Dahuang Zhechong Wan (大黄蛰虫丸) formula composition and its
TCMSP-documented bioactive constituents. Verify each compound against PubChem
(programmatic name->CID + canonical SMILES retrieval).

Formula source: Chinese Pharmacopoeia (2020); classic formula from Jin Gui Yao Lue.
Bioactive constituent selection: TCMSP active ingredients (oral bioavailability
OB>=30% & drug-likeness DL>=0.18) plus principal pharmacological markers per herb.
"""
import requests, csv, time, json, sys, os

os.makedirs("data/processed", exist_ok=True)
UA = {'User-Agent': 'Mozilla/5.0 (research/network-pharmacology)'}

# 12-herb composition of Dahuang Zhechong Wan (Chinese Pharmacopoeia 2020)
HERBS = {
    "Dahuang (Rhei Radix et Rhizoma, Rhubarb)":        "大黄",
    "Tubiechong (Eupolyphaga/Steleophaga, Ground beetle)": "土鳖虫",
    "Shuizhi (Hirudo, Leech)":                          "水蛭",
    "Mengchong (Tabanus, Gadfly)":                      "虻虫",
    "Qicao (Holotrichia, Grub)":                        "蛴螬",
    "Ganqi (Toxicodendri Resina, Dried lacquer)":       "干漆",
    "Taoren (Persicae Semen, Peach kernel)":            "桃仁",
    "Kuxingren (Armeniacae Semen Amarum, Apricot kernel)": "苦杏仁",
    "Huangqin (Scutellariae Radix)":                    "黄芩",
    "Dihuang (Rehmanniae Radix)":                       "地黄",
    "Baishao (Paeoniae Radix Alba)":                    "白芍",
    "Gancao (Glycyrrhizae Radix, Licorice)":            "甘草",
}

# Bioactive constituents: (compound_name, [herb Chinese names it derives from])
# Compiled from TCMSP active-ingredient records for each constituent herb.
COMPOUNDS = [
    # 大黄 Rhubarb anthraquinones / others
    ("emodin", ["大黄"]),
    ("aloe-emodin", ["大黄"]),
    ("rhein", ["大黄"]),
    ("chrysophanol", ["大黄"]),
    ("physcion", ["大黄"]),
    ("emodin-8-glucoside", ["大黄"]),
    ("torachrysone", ["大黄"]),
    ("sennoside A", ["大黄"]),
    # 黄芩 Scutellaria flavones
    ("baicalein", ["黄芩"]),
    ("baicalin", ["黄芩"]),
    ("wogonin", ["黄芩"]),
    ("wogonoside", ["黄芩"]),
    ("oroxylin A", ["黄芩"]),
    ("scutellarein", ["黄芩"]),
    ("skullcapflavone II", ["黄芩"]),
    # 甘草 Licorice
    ("glycyrrhizic acid", ["甘草"]),
    ("liquiritigenin", ["甘草"]),
    ("isoliquiritigenin", ["甘草"]),
    ("formononetin", ["甘草"]),
    ("glabridin", ["甘草"]),
    ("licochalcone A", ["甘草"]),
    ("glycyrol", ["甘草"]),
    # 白芍 Paeonia
    ("paeoniflorin", ["白芍"]),
    ("albiflorin", ["白芍"]),
    ("paeonol", ["白芍", "大黄"]),
    ("(+)-catechin", ["白芍", "桃仁"]),
    ("benzoic acid", ["白芍"]),
    # 桃仁 & 苦杏仁 kernels
    ("amygdalin", ["桃仁", "苦杏仁"]),
    ("prunasin", ["桃仁", "苦杏仁"]),
    ("oleic acid", ["桃仁", "苦杏仁"]),
    # 地黄 Rehmannia
    ("catalpol", ["地黄"]),
    ("aucubin", ["地黄"]),
    ("acteoside", ["地黄"]),
    # 干漆 Lacquer flavonoids
    ("fisetin", ["干漆"]),
    ("butein", ["干漆"]),
    ("butin", ["干漆"]),
    ("sulfuretin", ["干漆"]),
    ("gallic acid", ["干漆", "大黄"]),
    # Shared / ubiquitous phytosterols & flavonols (multiple herbs)
    ("quercetin", ["甘草", "桃仁", "苦杏仁", "黄芩"]),
    ("kaempferol", ["甘草", "黄芩", "白芍"]),
    ("luteolin", ["黄芩", "甘草"]),
    ("beta-sitosterol", ["大黄", "桃仁", "苦杏仁", "白芍", "地黄", "甘草"]),
    ("stigmasterol", ["桃仁", "苦杏仁", "白芍", "地黄"]),
    ("naringenin", ["甘草"]),
    ("medicarpin", ["甘草"]),
    ("apigenin", ["黄芩", "甘草"]),
]

def pubchem_cid_smiles(name):
    """Resolve name -> (CID, canonical SMILES) via PubChem PUG-REST."""
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
    try:
        r = requests.get(base + requests.utils.quote(name) + "/property/CanonicalSMILES,IUPACName/JSON",
                         headers=UA, timeout=25)
        if r.status_code == 200:
            props = r.json()["PropertyTable"]["Properties"][0]
            return props.get("CID"), props.get("CanonicalSMILES", "")
    except Exception as e:
        pass
    return None, None

rows = []
seen = {}
for name, herbs in COMPOUNDS:
    cid, smiles = pubchem_cid_smiles(name)
    time.sleep(0.22)  # PubChem rate limit courtesy
    status = "ok" if cid else "NOT_FOUND"
    print(f"{status:10s} {name:28s} CID={cid}")
    rows.append({"compound": name, "herbs_zh": ";".join(herbs),
                 "n_herbs": len(herbs), "pubchem_cid": cid or "", "smiles": smiles or ""})

with open("data/processed/compounds.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["compound","herbs_zh","n_herbs","pubchem_cid","smiles"])
    w.writeheader(); w.writerows(rows)

with open("data/processed/herbs.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["herb_en","herb_zh"])
    for en, zh in HERBS.items(): w.writerow([en, zh])

n_ok = sum(1 for r in rows if r["pubchem_cid"])
print(f"\n[compounds] {n_ok}/{len(rows)} resolved to PubChem CIDs across {len(HERBS)} herbs")
print("wrote data/processed/compounds.csv and herbs.csv")
