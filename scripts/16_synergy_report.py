#!/usr/bin/env python3
"""Step 16: synergy/complementarity figures + report (honest assessment + experiment roadmap)."""
import json, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

cov=list(csv.DictReader(open("results/tables/mechanism_coverage.csv")))
syn=json.load(open("data/processed/synergy_summary.json"))
pd=json.load(open("data/processed/thrombus_pd_model.json"))
pdrows={r["scenario"]:r for r in pd["rows"]}

# ---------- Figure 9: axis coverage + thrombus dynamics + trade-off ----------
fig=plt.figure(figsize=(15,5.2))
COL={"Control":"#7f8c8d","Rivaroxaban 10mg":"#3B6FA0","Formula 1.5g":"#2E7D5B","Combination":"#C0392B"}

# panel A: axis coverage (stacked count riva vs formula-adds)
axА=fig.add_subplot(1,3,1)
axes_labels=[c["axis"].split(" ")[-1] for c in cov]
riva_n=[int(c["riva_n"]) for c in cov]
formula_only=[max(int(c["combo_n"])-int(c["riva_n"]),0) for c in cov]
y=np.arange(len(cov))
axА.barh(y,riva_n,color="#3B6FA0",label="Rivaroxaban (FXa)")
axА.barh(y,formula_only,left=riva_n,color="#2E7D5B",label="Formula adds")
axА.set_yticks(y); axА.set_yticklabels(axes_labels,fontsize=9); axА.invert_yaxis()
axА.set_xlabel("# targets covered"); axА.set_title("Mechanism coverage across\nvenous-thrombosis axes",weight="bold",fontsize=11)
axА.legend(fontsize=8.5,loc="lower right")

# panel B: thrombus resolution curves
axB=fig.add_subplot(1,3,2)
for name,pr in pd["profiles"].items():
    axB.plot(pr["t"],pr["T"],lw=2.2,color=COL[name],label=name)
axB.axhline(0.5,ls=":",c="grey",lw=1); axB.text(0.5,0.52,"50% resolution",fontsize=8,color="grey")
axB.set_xlabel("days"); axB.set_ylabel("relative thrombus burden")
axB.set_title("Conceptual thrombus resolution\n(illustrative PD model)",weight="bold",fontsize=11)
axB.legend(fontsize=8.5); axB.set_ylim(0,1.15)

# panel C: efficacy-bleeding trade-off
axC=fig.add_subplot(1,3,3)
for name,r in pdrows.items():
    if r["t50_days"] is None: continue
    speed=1/ r["t50_days"]
    axC.scatter(r["bleeding_proxy"],speed,s=140,color=COL[name],zorder=5)
    axC.annotate(name,(r["bleeding_proxy"],speed),textcoords="offset points",xytext=(8,4),fontsize=8.5)
axC.set_xlabel("bleeding-risk proxy →"); axC.set_ylabel("resolution speed (1/t50) →")
axC.set_title("Faster healing vs bleeding trade-off",weight="bold",fontsize=11)
axC.grid(alpha=0.25)
plt.tight_layout(); plt.savefig("results/figures/fig9_synergy.png",dpi=150); plt.close()
print("wrote results/figures/fig9_synergy.png")

# ---------------- report ----------------
comb=pdrows["Combination"]; riv=pdrows["Rivaroxaban 10mg"]; fm=pdrows["Formula 1.5g"]
L=[]; w=L.append
w("# 附：协同增效 / 减毒 / 加速痊愈 —— 机制评估与实验验证路线图")
w("")
w("> 机制互补性(数据驱动) + 概念性 PD 模型(示意) + 严谨实验设计")
w("> **核心立场：以下是“有机制支撑的科学假说 + 如何验证”，不是“已被证实的疗效结论”。**")
w("")
w("---")
w("")
w("## 一、诚实的科学定性（先说结论）")
w("")
w("| 命题 | 机制层面 | 是否已证实 |")
w("|---|---|---|")
w("| **增效（互补覆盖）** | ✅ 合理：利伐沙班仅覆盖凝血 1 个轴，方剂补充血小板/炎症/消散/内皮 4 个轴 | ⚠️ 需实验(协同指数) |")
w("| **加速痊愈** | ✅ 有基础：方剂独家提供 MMP 血栓重塑 + 抗炎，促再通/机化清除 | ⚠️ 需动物+临床终点 |")
w("| **减毒（降低出血）** | ❌ 不成立：方剂**增加**而非降低抗栓负荷，出血方向叠加 | — |")
w("| **减毒（剂量/疗程节约）** | 🟡 可能：若加速痊愈，或可缩短疗程/维持低强度(10mg)→降低累积出血暴露 | ⚠️ 需临床验证 |")
w("| **减毒（非出血毒性）** | 🟡 可能：黄酮抗氧化/护内皮/护胃或缓解部分血管氧化损伤 | ⚠️ 推测 |")
w("")
w("> **一句话：** “协同增效、加速痊愈”是 **有强机制支撑的阳性假说**；但“减毒”若指“降低出血”则 **不成立**——"
  "合理的减毒只能是 **剂量/疗程节约** 与 **方剂自身的血管保护**。任何结论都须经下文实验验证。")
w("")
w("---")
w("")
w("## 二、增效的数据基础：机制互补性")
w("")
w(f"- 利伐沙班经 ChEMBL 确认为 **纯 FXa(F10) 抑制剂**（pChEMBL 9.4），仅覆盖 **凝血级联 1 个轴**。")
w(f"- 方剂覆盖 **{syn['formula_axes']}/{syn['axes']} 个轴**；二者靶点重叠仅 **F2**（Jaccard={syn['jaccard']}）→ **高度互补**。")
w(f"- **方剂独家补充利伐沙班完全缺失的 4 个轴**：{('、'.join(syn['formula_adds_axes']))}。")
w("")
w("| 血栓病理轴 | 利伐沙班 | 方剂补充的靶点 |")
w("|---|---|---|")
for c in cov:
    axis=c["axis"]; r=c["rivaroxaban"] or "—"; f=c["formula"] or "—"
    w(f"| {axis} | {r} | {f} |")
w("")
w("> **要点：** 利伐沙班“**防止血栓增大/播散**”（阻断 FXa），方剂“**促进已形成血栓的消散/再通 + 抗炎护内皮**”"
  "（MMP2/9/12 重塑、NLRP3/NF-κB 抗炎、KDR/HIF 内皮）。两者针对血栓“**从形成到清除**”的不同阶段，"
  "构成互补——这是“增效/加速痊愈”的机制根据。**注意：两者均未覆盖直接纤溶(PLG/tPA)，是共同空白。**")
w("")
w("![机制互补 + 血栓消散 + 权衡](../results/figures/fig9_synergy.png)")
w("")
w("---")
w("")
w("## 三、加速痊愈的概念性 PD 模型（示意，非证据）")
w("")
w("> ⚠️ 该模型速率常数为 **假设值**，用于 **形式化假说、生成可检验预测并展示出血权衡**，"
  "**不构成疗效证明**。")
w("")
w("| 方案 | 血栓消散半数时间 t50 | 第35天残余 | 出血代理 | 单位出血获益 |")
w("|---|---|---|---|---|")
for name in ["Control","Rivaroxaban 10mg","Formula 1.5g","Combination"]:
    r=pdrows[name]
    w(f"| {name} | {r['t50_days'] if r['t50_days'] else '未达50%'} d | {r['residual_day35']} | {r['bleeding_proxy']} | {r['resolution_per_bleed']} |")
w("")
w(f"**模型预测（待验证）：**")
w(f"- **联合方案血栓消散最快**（t50≈{comb['t50_days']}d vs 利伐沙班{riv['t50_days']}d vs 方剂{fm['t50_days']}d）、残余最少 → 支持“加速痊愈”。")
w(f"- **但联合方案出血代理最高**（{comb['bleeding_proxy']} vs 利伐沙班{riv['bleeding_proxy']}）→ **“减毒(降出血)”不成立**。")
w(f"- **单位出血的消散获益**：联合({comb['resolution_per_bleed']}) > 利伐沙班单用({riv['resolution_per_bleed']}) → "
  "提示联合可能以 **更短疗程/更低强度** 达到同等消散（**剂量-疗程节约式“减毒”**）；方剂单用最温和但最慢。")
w("")
w("---")
w("")
w("## 四、需要补充的科学实验（验证路线图）")
w("")
w("要把上述假说变成结论，须按 **计算 → 体外 → 体内 → 临床** four阶段验证，且每阶段都设 **出血安全性** 终点。")
w("")
w("### 阶段 0 — 计算/结构（1–2 月）")
w("- **分子对接 + MD 模拟**：quercetin/baicalein/isoliquiritigenin 对 F2、PTGS2、ALOX12、MMP9、NLRP3、KDR 的结合模式与亲和力（AutoDock Vina/Glide + GROMACS）。")
w("- **成分定量**：UPLC-MS/MS 实测同仁堂大蜜丸中 quercetin、baicalein、emodin、甘草酸、水蛭素活性等真实含量（校正本项目的含量估算）。")
w("- **完善 PBPK/PD**：将实测含量 + 水蛭素抗凝活性纳入定量出血叠加模型。")
w("")
w("### 阶段 1 — 体外（3–6 月）")
w("- **抗凝/抗栓功能试验**：PT、aPTT、**抗 FXa 活性**、**凝血酶生成试验(TGA/ETP)**、血栓弹力图(TEG)——比较利伐沙班、方剂提取物、联合。")
w("- **血小板聚集**（ADP/胶原/花生四烯酸诱导）、TXB2、12-HETE 测定（验证 COX/LOX 抑制）。")
w("- **纤溶/重塑**：D-二聚体、**MMP2/9 明胶酶谱(zymography)**、体外血栓块溶解试验。")
w("- **内皮/炎症**：HUVEC 管腔形成、NO 生成、NLRP3/IL-1β、NF-κB 报告基因。")
w("- **协同定量**：**Chou-Talalay 联合指数(CI)** 与 **等效线图(isobologram)**——判定协同(CI<1)/相加/拮抗。")
w("")
w("### 阶段 2 — 体内（6–12 月）")
w("- **静脉血栓动物模型**：大鼠/小鼠 **下腔静脉(IVC)狭窄或结扎淤滞模型**（最贴合静脉淤滞）、FeCl₃ 或光化学模型作补充。")
w("- **分组**：溶媒、利伐沙班、方剂(高/低剂量)、联合(析因设计)，含 **人体等效剂量换算**。")
w("- **疗效终点**：**血栓重量/长度、血栓消退率、血管再通率、组织学机化/再内皮化、炎症浸润**、时间曲线（验证“加速痊愈”）。")
w("- **安全性终点**：**尾出血时间、出血量、血红蛋白、内脏出血评分**（量化“减毒”真伪）。")
w("- **PK/PD**：联合下利伐沙班血浆浓度与抗 FXa（验证本项目 PBPK 的 AUCR≈1.03×预测）。")
w("- **机制确证**：核心靶点(F2、MMP9、NLRP3、KDR)的 WB/IHC/qPCR；必要时基因敲除/抑制剂反证。")
w("")
w("### 阶段 3 — 临床（1–3 年）")
w("- **设计**：远端 DVT/肌间静脉血栓患者的 **前瞻性随机对照(优先)或前瞻性队列**；利伐沙班 ± 大黄蛰虫丸。")
w("- **主要终点**：**血栓消退/再通时间（加速痊愈）**、血栓进展/复发率。")
w("- **次要/安全终点**：**大出血与临床相关非大出血(ISTH 标准)**、血红蛋白、肝肾功能、血钾血压（甘草）。")
w("- **样本量**：基于预试验效应量估算并预注册；设 **独立数据安全监察委员会(DSMB)**。")
w("- **药物警戒**：结合 **中西药联用不良反应监测/真实世界数据** 追踪出血信号。")
w("")
w("### 关键判定标准")
w("| 假说 | 支持性证据阈值 |")
w("|---|---|")
w("| 增效/协同 | 体外 CI<1 + 体内联合组血栓消退显著优于单药(析因交互 p<0.05) |")
w("| 加速痊愈 | 体内/临床 **再通时间显著缩短**、残余血栓更少 |")
w("| 减毒(剂量节约) | 联合达同等消退所需 **利伐沙班强度/疗程更低**，且大出血不劣于标准 |")
w("| 无 PK 风险 | 临床 PK 证实利伐沙班 AUC 无显著升高(呼应 AUCR≈1.03×) |")
w("")
w("---")
w("")
w("## 五、结论")
w("")
w("> **“利伐沙班 10 mg qd + 半粒 1.5 g 大黄蛰虫丸”具备“协同增效、加速痊愈”的坚实机制基础**"
  "（互补覆盖 5/6 病理轴、独家提供血栓消散与抗炎），**PK 相互作用可忽略**（AUCR≈1.03×）。")
w(">")
w("> **但：** (1) “减毒”若指降低出血则 **不成立**（抗栓作用叠加）；合理减毒只能是 **剂量/疗程节约** 与 **血管保护**；"
  "(2) “增效/加速痊愈”目前是 **计算与概念模型层面的阳性假说**，**必须经上文体外协同(CI)、动物血栓模型与临床 RCT 验证**方能成立。")
w(">")
w("> **临床落地前提：** 专科监护、出血监测、个体化风险评估——**不可据此自我联用**。")
w("")
w("> **免责声明：** 本节为计算性科研分析与实验设计建议，**不构成医疗建议**。")
w("")
open("docs/SYNERGY_AND_EXPERIMENTS.md","w").write("\n".join(L))
print(f"wrote docs/SYNERGY_AND_EXPERIMENTS.md ({len(L)} lines)")
