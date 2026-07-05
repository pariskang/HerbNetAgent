# HerbNetAgent — 大黄蛰虫丸 × 下肢肌间静脉血栓 网络药理学评估

Reproducible network-pharmacology assessment of whether **Dahuang Zhechong Wan
(大黄蛰虫丸)** has a plausible mechanistic basis for treating **lower-limb
intermuscular (distal) venous thrombosis**.

> ⚠️ **免责声明 / Disclaimer:** 本项目为**计算性科研分析，仅供研究参考，不构成医疗建议**。
> 网络药理学提供的是**机制合理性假说**，**不能替代随机对照试验（RCT）**证明临床疗效。
> This is a computational research analysis for hypothesis generation only — **not medical advice**,
> and **not** a substitute for clinical trials.

## 📌 核心结论 / Key finding

大黄蛰虫丸在分子机制层面**具备明确、多层次的抗静脉血栓作用基础**（有强机制支撑的阳性假说），
但"机制成立 ≠ 临床确有疗效"。完整评估见 **[`docs/REPORT.md`](docs/REPORT.md)**。

- **12 味药 → 46 个 PubChem 核验成分 → 208 个 ChEMBL 实测靶点**
- 与 **静脉血栓疾病基因（Open Targets, 3286 个）** 交集 → **16 个高置信核心机制靶点**
- **头号发现：凝血酶 F2**（易栓症头号基因，关联评分 0.80）被黄酮 **quercetin 直接结合（≈39 nM）**，
  与水蛭素(hirudin)直接抑制凝血酶相互印证
- 富集通路：**止血 / 血小板活化 / 花生四烯酸代谢 / HIF-1(淤滞) / VEGF-内皮 / TNF-白介素炎症** —— 对应 Virchow 三要素

## 💊 联用分析 / Rivaroxaban co-administration

**大黄蛰虫丸能否与利伐沙班（Rivaroxaban）联用？** 见 **[`docs/RIVAROXABAN_DDI.md`](docs/RIVAROXABAN_DDI.md)**。

数据驱动结论：**不建议无专科监护下常规联用**，存在**双重出血风险放大**——
- **药效学(PD)**：利伐沙班抑制 FXa，方剂命中凝血酶 F2 + 血小板 COX/LOX + 含水蛭素 → 抗凝作用叠加；
- **药代(PK)**：方中 quercetin / licochalcone A / apigenin 等抑制 **CYP3A4 + P-gp + BCRP**（利伐沙班的代谢酶与外排转运体）→ 可能升高其血药浓度。quercetin 三条清除通路全抑制，licochalcone A 为强效 BCRP 抑制剂（IC50≈6 nM）。

## 🔬 方法与数据源 / Pipeline & data sources

全流程由公开一级数据库 API 驱动，**无人工挑数、可独立复现**：

| 环节 | 数据源 | 关键方法 |
|---|---|---|
| 成分核验 | **PubChem** PUG-REST | name → CID → SMILES |
| 成分–靶点 | **ChEMBL** REST | 人源单蛋白、**实测** pChEMBL≥5 |
| 疾病–靶点 | **Open Targets** GraphQL | 5 个静脉血栓节点并集 + 关联评分 |
| PPI 网络 | **STRING** v12 REST | 置信度≥0.4 + networkx 中心性 |
| 功能富集 | **Enrichr** | GO / KEGG / Reactome, BH adj-p<0.05 |

## 📂 目录结构 / Layout

```
scripts/     01_compounds → 07_report  (7-step pipeline, 可依次运行)
data/        raw/ + processed/  (成分、靶点、疾病基因等中间数据)
results/     tables/ (交集靶点、PPI边、富集表) + figures/ (5 张图)
docs/        REPORT.md  ← 完整中文评估报告（含图表）
```

## ▶️ 复现 / Reproduce

```bash
pip install -r requirements.txt
python3 scripts/01_compounds.py        # 方剂成分 → PubChem 核验
python3 scripts/02_compound_targets.py # 成分 → ChEMBL 实测靶点 (~数分钟)
python3 scripts/03_disease_targets.py  # 疾病 → Open Targets 靶点
python3 scripts/04_intersection_ppi.py # 交集 + STRING PPI + hub
python3 scripts/05_enrichment.py       # Enrichr 富集
python3 scripts/05b_curate_pathways.py # 筛选血栓相关通路
python3 scripts/06_figures.py          # 生成图表
python3 scripts/07_report.py           # 生成 docs/REPORT.md
python3 scripts/08_rivaroxaban_ddi.py  # 利伐沙班联用: CYP3A4/P-gp/BCRP + PD 分析
python3 scripts/09_ddi_report.py       # 生成 docs/RIVAROXABAN_DDI.md + 图
```

## 🖼️ 主要图表 / Figures

| 图 | 内容 |
|---|---|
| `fig1_venn.png` | 方剂靶点 ∩ 静脉血栓靶点 |
| `fig2_ppi_network.png` | PPI 网络（红=血栓核心靶点，橙=通用枢纽） |
| `fig3_hubs.png` | Top-20 网络度 hub |
| `fig4_enrichment.png` | 血栓/心血管/炎症相关通路富集 |
| `fig5_herb_contribution.png` | 各药材对核心靶点的贡献 |
