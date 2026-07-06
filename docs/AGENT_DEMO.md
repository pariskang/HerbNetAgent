# HerbNetAgent 2.0 — 评估结果

> 疾病: **venous thrombosis** · 联用: **rivaroxaban** · 摄动剂量: 0.006 mg

## L4 药效 (Efficacy)
- 方剂靶点 208，疾病模块 626，共享机制靶点 **16**：ALOX12, ALOX15, CDK6, DAPK1, F2, GRK6, IL2, JAK2, MAPT, NEK6, NFKB1, NLRP3, PIK3R1, PTGS1, PTGS2
- 方剂覆盖病理轴：Coagulation, Platelet, Thrombo-inflammation, Resolution, Endothelium
- 与 rivaroxaban 的靶点重叠 Jaccard=0.0048；方剂独家补充轴：Platelet, Thrombo-inflammation, Resolution, Endothelium
- 网络拓扑：分离度 s=2.137（>0 表示两靶点集拓扑分离）、Jaccard=0.0048；双方均engage疾病模块=True
- **互补暴露 (Cheng/Barabási): True** （判据：双方均engage疾病模块 且 两靶点集拓扑分离 → 有效药物组合的网络特征）
  - _注：proximity z uses a STRING-induced subgraph (module-biased null) — interpret cautiously; the robust complementarity signals are separation and module engagement._

### L3 结合亲和力（实测）
- **ALOX12** ← quercetin: 436.5 nM (Ki/IC50-eq)  ⟨体外实测·ChEMBL⟩
- **ALOX15** ← baicalein: 602.6 nM (Ki/IC50-eq)  ⟨体外实测·ChEMBL⟩
- **CDK6** ← fisetin: 851.1 nM (Ki/IC50-eq)  ⟨体外实测·ChEMBL⟩
- **DAPK1** ← isoliquiritigenin: 2454.7 nM (Ki/IC50-eq)  ⟨体外实测·ChEMBL⟩
- **F2** ← quercetin: 38.9 nM (Ki/IC50-eq)  ⟨体外实测·ChEMBL⟩
- **GRK6** ← scutellarein: 6606.9 nM (Ki/IC50-eq)  ⟨体外实测·ChEMBL⟩
- **IL2** ← formononetin: 5495.4 nM (Ki/IC50-eq)  ⟨体外实测·ChEMBL⟩
- **JAK2** ← sulfuretin: 9332.5 nM (Ki/IC50-eq)  ⟨体外实测·ChEMBL⟩

## L5 药代 (PK / DDI)
- 动态 PBPK 预测 **AUCR = 1.035× (1.035 x AUC  ⟨计算预测·reduced-PBPK(rivaroxaban/quercetin)⟩)**
- 模型验证（强抑制剂）AUCR ≈ 2.48× (临床酮康唑 ~2.6×)

## L7 安全 (Safety / DDI)
- PK 相互作用：**negligible** · PD 出血叠加：**True**
  - victim exposure 1.03x (negligible)  ⟨计算预测·reduced-PBPK(rivaroxaban/quercetin)⟩
  - additive bleeding risk (PD)  ⟨计算预测·mechanism-coverage⟩
- **结论**：PK 相互作用可忽略；主要残余风险为 PD 出血叠加 —— 需临床监护，不作为自我用药依据（证据仅达计算级）

> ⚠️ 本智能体所有输出的最高证据等级为 **计算预测**；临床用药决策须由医师依据体外/动物/RCT 证据个体化判断。