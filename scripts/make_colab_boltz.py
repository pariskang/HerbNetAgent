#!/usr/bin/env python3
"""Generate notebooks/boltz2_F2_MMP9_validation.ipynb — a GPU Colab notebook that
runs Boltz-2 for real on F2/MMP9 x formula-compound pairs and validates predicted
affinity against the measured ChEMBL pChEMBL values used by HerbNetAgent's L3."""
import json, os

# measured ground-truth (ChEMBL, from data/processed/compound_targets.csv)
GROUND_TRUTH = [
    ("F2", "P00734", "quercetin",         "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O", 7.41),
    ("F2", "P00734", "beta-sitosterol",   "CC[C@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C)C(C)C", 6.57),
    ("MMP9", "P14780", "isoliquiritigenin","C1=CC(=CC=C1/C=C/C(=O)C2=C(C=C(C=C2)O)O)O", 8.00),
    ("MMP9", "P14780", "quercetin",        "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O", 6.76),
    ("MMP9", "P14780", "luteolin",         "C1=CC(=C(C=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O)O", 5.25),
]

def md(src): return {"cell_type": "markdown", "metadata": {}, "source": src.splitlines(keepends=True)}
def code(src): return {"cell_type": "code", "metadata": {}, "execution_count": None,
                       "outputs": [], "source": src.splitlines(keepends=True)}

cells = []

cells.append(md(
"""# HerbNetAgent · L3 结构基础亲和力验证 — Boltz-2 on F2 / MMP9

在 **GPU** 上用 **Boltz-2**（共折叠 + 结合亲和力）对大黄蛰虫丸关键成分与
**F2（凝血酶）/ MMP9（明胶酶B）** 做真实预测，并与 **ChEMBL 实测 pChEMBL** 对比校验。

> Runtime → Change runtime type → **T4 GPU**。整跑约 15–30 min（5 个复合物 × MSA + 扩散）。
> 产出 `boltz_predictions.json`，格式与 `herbnetagent` 的 `Boltz2Backend` 缓存兼容，可直接放回仓库
> `data/structures/boltz_cache.json` 让智能体用**预测值**替代实测回退。

**Boltz-2 亲和力约定（已核实）**：`affinity_pred_value = log10(IC50/µM)`，**越低越强**；
`pIC50 = 6 − affinity_pred_value`，`IC50[nM] = 10^(value+3)`；`affinity_probability_binary` = 结合概率。
"""))

cells.append(md("## 0. 环境检查 + 安装 Boltz-2"))
cells.append(code(
"""!nvidia-smi -L || echo 'NO GPU — set Runtime>Change runtime type>T4 GPU'
!pip -q install boltz requests pandas matplotlib scipy
import subprocess, shutil
print('boltz on PATH:', shutil.which('boltz'))
"""))

cells.append(md("## 1. 目标 + 探针（含 ChEMBL 实测真值）"))
cells.append(code(
"GROUND_TRUTH = " + json.dumps(GROUND_TRUTH, ensure_ascii=False, indent=4) + "\n"
"""# (gene, uniprot, compound, smiles, measured_pChEMBL)
import pandas as pd
gt = pd.DataFrame(GROUND_TRUTH, columns=['gene','uniprot','compound','smiles','measured_pChEMBL'])
gt
"""))

cells.append(md("## 2. 拉取靶点序列（UniProt）"))
cells.append(code(
"""import requests, functools
@functools.lru_cache
def uniprot_seq(acc):
    r = requests.get(f'https://rest.uniprot.org/uniprotkb/{acc}.fasta', timeout=60)
    return ''.join(l for l in r.text.splitlines() if not l.startswith('>'))
SEQ = {acc: uniprot_seq(acc) for acc in gt['uniprot'].unique()}
for a,s in SEQ.items(): print(a, len(s), 'aa')
"""))

cells.append(md("## 3. 生成 Boltz-2 输入 YAML（protein 链 + ligand + affinity）"))
cells.append(code(
"""import os, re
os.makedirs('jobs', exist_ok=True)
def job_name(gene, compound): return re.sub(r'[^A-Za-z0-9]+','_', f'{gene}_{compound}')
def write_yaml(gene, acc, compound, smiles):
    name = job_name(gene, compound)
    yml = (
        'version: 1\\n'
        'sequences:\\n'
        '  - protein:\\n'
        '      id: A\\n'
        f'      sequence: \"{SEQ[acc]}\"\\n'
        '  - ligand:\\n'
        '      id: B\\n'
        f'      smiles: \"{smiles}\"\\n'
        'properties:\\n'
        '  - affinity:\\n'
        '      binder: B\\n')
    p = f'jobs/{name}.yaml'
    open(p,'w').write(yml); return name, p
jobs = [write_yaml(r.gene, r.uniprot, r.compound, r.smiles) for r in gt.itertuples()]
jobs
"""))

cells.append(md(
"""## 4. 运行 Boltz-2（每个复合物）

`--use_msa_server` 用公共 MSA 服务器自动建 MSA。首个任务会下载权重（~数 GB）。
"""))
cells.append(code(
"""import subprocess, glob, json
os.makedirs('out', exist_ok=True)
def run_boltz(yaml_path):
    cmd = ['boltz','predict', yaml_path, '--use_msa_server', '--out_dir','out']
    print('>>', ' '.join(cmd))
    subprocess.run(cmd, check=False)
for name, path in jobs:
    run_boltz(path)
print('done')
"""))

cells.append(md("## 5. 解析预测亲和力 → pIC50 → 对比实测"))
cells.append(code(
"""def parse_affinity(name):
    hits = glob.glob(f'out/**/affinity_{name}.json', recursive=True) or \\
           glob.glob(f'out/**/predictions/**/*affinity*.json', recursive=True)
    if not hits: return None
    d = json.load(open(hits[0]))
    v = d.get('affinity_pred_value')
    if v is None: return None
    return dict(affinity_pred_value=float(v),
                pIC50_pred=round(6.0-float(v),2),
                IC50_nM_pred=round(10**(float(v)+3),1),
                P_binder=d.get('affinity_probability_binary'))
rows=[]
for r in gt.itertuples():
    name = job_name(r.gene, r.compound)
    res = parse_affinity(name) or {}
    rows.append(dict(gene=r.gene, compound=r.compound,
                     measured_pChEMBL=r.measured_pChEMBL, **res))
res_df = pd.DataFrame(rows)
res_df
"""))

cells.append(md("## 6. 校验：预测 pIC50 vs 实测 pChEMBL（相关性 + RMSE）"))
cells.append(code(
"""import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
v = res_df.dropna(subset=['pIC50_pred'])
if len(v) >= 2:
    r,_ = pearsonr(v['measured_pChEMBL'], v['pIC50_pred'])
    rmse = float(np.sqrt(np.mean((v['measured_pChEMBL']-v['pIC50_pred'])**2)))
    print(f'n={len(v)}  Pearson r={r:.2f}  RMSE={rmse:.2f} pIC50 units')
    fig,ax=plt.subplots(figsize=(5,5))
    ax.scatter(v['measured_pChEMBL'], v['pIC50_pred'], s=90, c='#2E7D5B')
    for _,row in v.iterrows():
        ax.annotate(f"{row['gene']}·{row['compound']}", (row['measured_pChEMBL'], row['pIC50_pred']),
                    fontsize=7, xytext=(4,4), textcoords='offset points')
    lo,hi=5,9; ax.plot([lo,hi],[lo,hi],'--',c='grey'); ax.set_xlim(lo,hi); ax.set_ylim(lo,hi)
    ax.set_xlabel('measured pChEMBL (ChEMBL)'); ax.set_ylabel('Boltz-2 predicted pIC50')
    ax.set_title('L3 structure-based affinity validation\\nBoltz-2 vs measured (F2 / MMP9)')
    plt.tight_layout(); plt.savefig('boltz2_validation.png', dpi=150); plt.show()
else:
    print('Not enough parsed predictions — check the boltz run logs above.')
"""))

cells.append(md(
"""## 7. 导出为 HerbNetAgent 缓存格式

把预测写成 `boltz_cache.json`（键 `compound|gene`），放回仓库 `data/structures/` 后，
`Boltz2Backend` 会优先返回**预测值**（证据级=计算预测），否则回退 ChEMBL 实测。
"""))
cells.append(code(
"""cache={}
for r in res_df.dropna(subset=['IC50_nM_pred']).itertuples():
    cache[f'{r.compound}|{r.gene}'] = dict(
        affinity_pred_value=float(6.0-r.pIC50_pred),
        affinity_probability_binary=r.P_binder,
        IC50_nM=float(r.IC50_nM_pred), pIC50=float(r.pIC50_pred))
json.dump(cache, open('boltz_cache.json','w'), indent=1)
print(json.dumps(cache, indent=1))
from google.colab import files  # optional download
# files.download('boltz_cache.json'); files.download('boltz2_validation.png')
"""))

cells.append(md(
"""---
> ⚠️ Boltz-2 亲和力为**计算预测**；此处 5 点为工作流演示（样本小）。
> 严谨校验应扩大到更多成分-靶点对并与实测/实验对照。**不构成医疗建议。**
"""))

nb = {"cells": cells,
      "metadata": {"accelerator": "GPU",
                   "colab": {"provenance": [], "gpuType": "T4"},
                   "kernelspec": {"display_name": "Python 3", "name": "python3"},
                   "language_info": {"name": "python"}},
      "nbformat": 4, "nbformat_minor": 0}

os.makedirs("notebooks", exist_ok=True)
out = "notebooks/boltz2_F2_MMP9_validation.ipynb"
json.dump(nb, open(out, "w"), ensure_ascii=False, indent=1)
print(f"wrote {out} ({len(cells)} cells)")
