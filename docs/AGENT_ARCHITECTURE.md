# HerbNetAgent 2.0 — 面向中药定量计算的多尺度智能体架构

> **一个把中药"成分—结构—亲和力—网络—暴露—效应—安全"全链条定量打通、
> 客观评估药效(efficacy)与药代(PK/ADMET)的创新性 AI 智能体蓝图。**
> 设计锚定当前最顶级研究成果，并以本仓库 01–18 步流程作为**可运行原型**验证其可行性。

![架构总览](../results/figures/fig11_agent_architecture.png)

---

## 一、要解决的空白（为什么需要"新一代"智能体）

现有中药计算工具各自为战、且**多停留在定性网络药理学**：

| 现状 | 局限 |
|---|---|
| TCMSP / HERB / **BATMAN-TCM 2.0**（170 万+预测互作）| 是**数据库**，给"成分→靶点"清单，不做定量效应/暴露 |
| 经典网络药理学（靶点交集 + 富集）| **定性**："有关系"，不回答"多大剂量、多强效应、是否达到有效浓度" |
| ADMET-AI / ADMETlab 3.0（119 端点）| 强在**单分子**性质预测，不做整方/多成分系统暴露 |
| AlphaFold3 / Boltz-2（结构+亲和力）| 强在**单复合物**，未与中药整方网络/PBPK 打通 |
| ChemCrow / Coscientist（LLM 科学智能体）| 通用化学/自动实验，**无中药系统药理学专用编排** |

> **核心创新命题：** 把上述前沿"点能力"编排成一个**端到端定量闭环**——
> 从**分子结构**到**结合亲和力**到**靶点网络**到**全身 PBPK 暴露**到**系统药效(QSP)**到**安全/毒理**，
> 每一步都带**不确定度**与**证据分级**，并显式区分"计算合理性"与"临床证据"。
> 这正是本仓库对大黄蛰虫丸 × 利伐沙班案例做到的（步骤 01–18），本文将其升维为通用智能体。

---

## 二、设计原则（Design Principles）

1. **定量优先（Quantitative-first）**：一切结论落到数值 + 置信区间（亲和力 nM、暴露 AUC、效应 % 、AUCR 倍数），而非"相关/富集"。
2. **多尺度贯通（Multi-scale）**：molecule → protein structure → affinity → interactome → pathway → PBPK exposure → QSP phenotype → outcome。
3. **整方建模（Formula-level）**：原生处理"12 味药 / 数十成分"的混合物——按**真实含量 × 生物利用度**加权，建模**联合扰动**，而非单分子。
4. **不确定度贯穿（Uncertainty-aware）**：每个预测带置信度；低置信触发更高保真方法或实验建议。
5. **可证伪 & 自设计实验（Falsifiable）**：智能体为自己的每条结论**自动生成可验证/可证伪的实验方案**（呼应 Coscientist 自驱实验思路）。
6. **负责任 AI（Responsible）**：显式"计算 ≠ 临床"分级、DDI/禁忌/毒理护栏、人在环（human-in-the-loop）临床把关。
7. **可复现 & 溯源（Provenance）**：每个数值可追溯到一级数据库/模型版本（本仓库全程 API 驱动即为示范）。

---

## 三、八层架构（L0–L7）

### L0 — 编排与推理内核（LLM Agent Core）
- **职责**：任务分解、工具路由、记忆、不确定度与溯源追踪、报告生成。
- **顶级思路**：**ChemCrow**（用 LLM 作推理引擎调用 18 个化学工具，Nature Machine Intelligence 2024）、**Coscientist**（GPT-4 自主设计/执行实验，Nature 2023）、多智能体 DrugAgent。LLM **不靠内部知识作答**，而是**调用可信工具**并核验。
- **本仓库原型**：整个 01–18 步的编排即为一个手工版 L0。

### L1 — 知识与数据层（统一知识图谱）
- **职责**：把异构数据统一为知识图谱（herb–compound–target–pathway–disease–ADR）。
- **数据源**：中药 **TCMSP / HERB / BATMAN-TCM 2.0 / SymMap / TCMBank**；化学 **PubChem / ChEMBL**；靶点—疾病 **Open Targets / DisGeNET**；互作 **STRING**；通路 **KEGG / Reactome**；药物警戒 **FAERS / 中西药联用不良反应库**。图谱范式参考 **PrimeKG / Hetionet**。
- **本仓库原型**：步骤 01–03（PubChem 核验、ChEMBL 靶点、Open Targets 疾病基因）。

### L2 — 分子表征与结构预测
- **职责**：成分结构质控（RDKit）、**整方真实含量定量（UPLC-MS/MS）**、靶点蛋白结构（无晶体时用 **AlphaFold3 / ESMFold**）。
- **本仓库原型**：步骤 01（PubChem SMILES/CID）、步骤 10（含量估算——真实系统应换成 UPLC-MS/MS 实测）。

### L3 — 结合与靶点占据引擎（Binding & Target Engagement）
- **职责**：给出**定量结合亲和力**，取代/校正纯数据库靶点。
- **顶级方法**：**Boltz-2**（MIT/Recursion 2025，共折叠 + **结合亲和力**联合预测，approaching FEP 精度、远快于物理法）、**AlphaFold3**、**Chai-1 / Protenix**（共折叠）；**DiffDock**（ML 对接）、**AutoDock Vina**；**FEP+** 精修；靶点预测 **SwissTargetPrediction / SEA**。
- **本仓库原型**：步骤 02 用 ChEMBL **实测** pChEMBL 作为亲和力真值——真实系统对缺数据的成分—靶点用 Boltz-2 补齐。

### L4 — 网络与药效引擎（EFFICACY / 药效）
- **职责**：把"成分作用于哪些靶点"升级为"**对疾病模块的定量扰动 + 组合协同**"。
- **顶级方法**：**Barabási 网络医学**——疾病模块、**网络邻近度（network proximity）**、以及药物组合的 **互补暴露原理（Complementary Exposure, Cheng et al., Nature Communications 2019）**：只有对疾病模块"互补而不冗余"覆盖的组合才显著有效——**这正是本项目量化利伐沙班(FXa)+方剂(5/6 轴)互补性的理论基石**。叠加 GO/KEGG 富集、机制覆盖矩阵。
- **本仓库原型**：步骤 04（交集+PPI+hub）、05（富集）、14（机制互补/Complementary Exposure）、17（骨愈合靶点方向）。

### L5 — ADMET / 药代引擎（PK / 药代）
- **职责**：从"能不能结合"到"**口服后能否在靶部位达到有效浓度**"。
- **顶级方法**：**ADMET-AI**（Chemprop-RDKit，TDC 41 数据集）、**ADMETlab 3.0**（119 端点、DMPNN、**带不确定度**、API）预测吸收/CYP 代谢/转运体/毒性；再驱动**全身 PBPK/QSP**——**Open Systems Pharmacology 的 PK-Sim / MoBi**（开源、含人/多种动物生理库）+ **虚拟人群（virtual populations）** 刻画个体差异；DDI 用静态 [I]/Ki 筛查 + 动态 PBPK。
- **本仓库原型**：步骤 08（CYP/P-gp/BCRP 实测抑制）、10（静态 [I]/Ki 剂量模型）、12（**动态 PBPK**，校准酮康唑 AUCR≈2.5×，预测方剂对利伐沙班仅 +3%）。

### L6 — 系统药效与结局模拟（Systems PD & Outcome）
- **职责**：把暴露→靶点占据→通路→**表型/临床结局**用机制模型连起来，输出可量化终点。
- **顶级方法**：**定量系统药理学（QSP）** 机制 ODE 模型 + 虚拟人群；疾病/组织专用模型（凝血—血栓—出血；骨折愈合分期）。
- **本仓库原型**：步骤 12（利伐沙班 PK/PD）、15（血栓消散概念 PD 模型）、18（骨折愈合分期靶点映射）。

### L7 — 安全、毒理与临床转化（Safety & Translation）
- **职责**：毒性预测、DDI/禁忌、药物警戒信号、**证据分级**、**自动生成验证实验/临床方案**。
- **顶级方法**：**ProTox-III / hERG / DILI** 毒性预测；FAERS 真实世界信号；证据分级明确区分"计算级 vs 体外 vs 动物 vs RCT"；实验方案自动化（Chou-Talalay CI、动物血栓模型、RCT 设计）。
- **本仓库原型**：步骤 09/11（DDI 安全报告）、16（协同 + 实验路线图）、临床情景简报（骨折患者沟通材料）。

### 横切关注点（Cross-cutting）
- **不确定度量化**（贝叶斯/集成/conformal）、**溯源与可解释性**、**人在环临床护栏**——贯穿所有层。

---

## 四、两大客观评估引擎

### （A）药效引擎（Efficacy） = L2 → L3 → L4
```
成分结构 → (Boltz-2)定量亲和力 → 靶点占据谱 → 疾病模块网络邻近度
        → 通路富集 → [组合] 互补暴露协同评分 → 机制覆盖图 + 定量药效假设
```
**输出**：疾病模块扰动强度、hub 靶点优先级、组合协同/拮抗判定（网络邻近度）、机制覆盖率——**全部带置信度**。

### （B）药代引擎（PK/ADMET） = L2 → L5 → L6
```
成分结构+真实含量 → (ADMET-AI/ADMETlab)吸收·CYP·转运体·毒性
        → (PK-Sim PBPK + 虚拟人群) 血浆/组织暴露时程 → DDI(AUCR倍数)
        → (QSP ODE) 暴露→效应→结局 → 个体化预测 + 不确定区间
```
**输出**：Cmax/AUC、达峰、组织暴露、DDI 定量倍数（如本项目 AUCR≈1.03×）、有效浓度是否达标、出血/毒性风险——**全部带虚拟人群变异带**。

---

## 五、创新点（相较现有工作的"新"）

1. **端到端定量闭环**：首次把"结构→亲和力(Boltz-2)→网络→PBPK→QSP→安全"串成一条**带数值与不确定度**的中药评估流水线；现有 TCM 工具止步于定性富集。
2. **整方多成分定量**：以**真实含量 × 生物利用度**加权建模整方"联合扰动"，而非单分子——解决中药"复方"本质。
3. **网络互补暴露 → 中西药组合**：用 Cheng/Barabási 框架把"TCM+化药联用"从经验升级为**可计算的协同/冗余判定**（本项目已量化利伐沙班+方剂互补 5/6 轴）。
4. **可证伪与自设计实验**：智能体为每条结论自动产出**分级验证方案**（体外 CI → 动物模型 → RCT），把"假设"与"证据"闭环。
5. **不确定度 + 证据分级为一等公民**：显式输出"计算/体外/动物/临床"证据等级，杜绝"网络药理学=疗效"的过度解读。
6. **负责任临床转化层**：内置 DDI/禁忌/毒理护栏与人在环把关（本项目对"骨折患者是否联用"坚持交由临床即为示范）。

---

## 六、不确定度、验证与负责任 AI

- **不确定度**：亲和力（Boltz-2 置信/pLDDT-类）、ADMET（ADMETlab 3.0 自带不确定度）、PBPK（虚拟人群变异带）、网络（置换检验 p 值）逐层传播；低置信自动升级方法或转实验。
- **前瞻验证**：模型须以**已知真值**校准（如本项目用酮康唑 AUCR≈2.6× 校准 PBPK；用 ChEMBL 实测亲和力校准 L3）。
- **证据分级**：`计算预测 < 体外 < 动物 < 人体 PK < RCT`——报告显式标注每条结论所处等级。
- **护栏**：孕妇/出血/肝肾禁忌、DDI 红线、"不构成医疗建议"贯穿；**最终临床决策交由医师**。

---

## 七、实施路线图与技术栈

| 阶段 | 目标 | 关键技术 |
|---|---|---|
| **P1 (0–6 月)** 打通骨架 | L0–L4 定量化 | 现有仓库 + Boltz-2/AlphaFold3 接入 L3；知识图谱 (Neo4j/PrimeKG) |
| **P2 (6–12 月)** 药代闭环 | L5–L6 | ADMET-AI/ADMETlab API + **PK-Sim/MoBi** PBPK/QSP + 虚拟人群 |
| **P3 (12–18 月)** 智能编排 | L0 自主化 | LangGraph/多智能体；工具化每层；自动实验方案 (L7) |
| **P4 (18–30 月)** 验证 | 前瞻校准 | 回顾性 DDI/疗效基准；与湿实验/临床数据闭环 |

**技术栈**：Python（RDKit, scipy, networkx）· Boltz-2/AlphaFold3（GPU）· ADMET-AI/ADMETlab 3.0 · PK-Sim/MoBi(OSP) · Open Targets/ChEMBL/STRING API · Neo4j 知识图谱 · LangGraph LLM 编排 · conformal/贝叶斯不确定度。

---

## 八、本仓库 = 可运行原型（Proof of Concept）

| 智能体层 | 本仓库对应步骤 | 已验证 |
|---|---|---|
| L1 数据 | 01–03 | PubChem/ChEMBL/Open Targets 全 API 驱动 ✅ |
| L2 结构/表征 | 01,10 | 46 成分核验；含量估算 ✅ |
| L3 结合 | 02 | ChEMBL 实测亲和力（真系统换 Boltz-2）✅ |
| L4 药效 | 04,05,14,17 | 交集/PPI/富集/**互补暴露**/骨愈合方向 ✅ |
| L5 药代 | 08,10,12 | 静态 DDI + **动态 PBPK（校准酮康唑）** ✅ |
| L6 系统 PD | 12,15,18 | 血栓消散 + 骨折分期模型 ✅ |
| L7 安全/转化 | 09,11,16 | DDI 报告 + 实验路线图 + 临床简报 ✅ |

> 即：本项目对"大黄蛰虫丸 × 利伐沙班 × 骨折"案例，已**手工跑通** L1–L7 的一条纵切面（含定量 PBPK AUCR、网络互补、骨愈合方向、证据分级）——**HerbNetAgent 2.0 就是把这条纵切面工具化、自动化、并接入 Boltz-2/PK-Sim 等最强模块**。

---

## 九、参考（顶级研究成果）

- **Boltz-2**（共折叠 + 结合亲和力）: https://pmc.ncbi.nlm.nih.gov/articles/PMC12262699/
- **AlphaFold3 验证**: https://www.ebi.ac.uk/training/online/courses/alphafold/alphafold-3-and-alphafold-server/
- **共折叠开源对比 (ABCFold)**: https://academic.oup.com/bioinformaticsadvances/article/5/1/vbaf153/8176613
- **网络药物组合 / 互补暴露 (Cheng et al., Nature Comms 2019)**: https://www.nature.com/articles/s41467-019-09186-x
- **网络医学虚拟筛选综述**: https://www.mdpi.com/1424-8247/17/7/899
- **ADMETlab 3.0 (NAR 2024)**: https://academic.oup.com/nar/article/52/W1/W422/7640525
- **ADMET-AI (TDC)**: https://www.biorxiv.org/content/10.1101/2023.12.28.573531
- **ChemCrow (Nature Mach. Intell. 2024)**: https://www.nature.com/articles/s42256-024-00832-8
- **Coscientist (Nature 2023) / 综述**: https://pubs.acs.org/doi/10.1021/jacsau.6c00213
- **BATMAN-TCM 2.0 (NAR 2024)**: https://academic.oup.com/nar/article/52/D1/D1110/7334089
- **AI in TCM 药理综述 2025**: https://journals.sagepub.com/doi/10.1177/1934578X251405983
- **Open Systems Pharmacology (PK-Sim/MoBi)**: https://www.open-systems-pharmacology.org/
- **QSP 虚拟人群**: https://www.biorxiv.org/content/10.1101/2025.06.24.661262

---

> **免责声明：** 本文为**科研方法学与系统架构设计**，其临床用途须经严格实验与监管验证；
> 智能体的任何输出均为**计算假设**，不构成医疗建议。
