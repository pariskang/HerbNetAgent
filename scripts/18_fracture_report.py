#!/usr/bin/env python3
"""Step 18: fracture-healing timeline figure + objective assessment report."""
import json, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

rows=list(csv.DictReader(open("results/tables/fracture_targets.csv")))
S=json.load(open("data/processed/fracture_summary.json"))

# ---- Figure: healing-phase timeline with target actions ----
# phases (week ranges)
phases=[("Inflammation /\nhematoma",0,1,"#f5d5d5"),
        ("Soft callus\n(cartilage)",1,3,"#f5ecc9"),
        ("Hard callus\n(bone)",3,8,"#d5e8d5"),
        ("Remodeling",6,14,"#d5dce8")]
# target -> approx week to place (by phase code)
phase_week={"I":0.5,"K":0.5,"C":2,"A":5,"O":5,"R":9,"I/O":4}
colmap={"⚠️":"#C0392B","🟢":"#2E7D5B","🟡":"#E0A63C","中":"#888","—":"#bbb"}
def eff_color(e):
    for k,c in colmap.items():
        if e.startswith(k): return c
    return "#888"

fig,ax=plt.subplots(figsize=(13,6.2))
for name,x0,x1,c in phases:
    ax.add_patch(Rectangle((x0,0),x1-x0,10,color=c,alpha=0.55,zorder=0))
    ax.text((x0+x1)/2,9.4,name,ha="center",va="top",fontsize=9,weight="bold")
# patient marker at ~6 weeks
ax.axvline(6,color="k",ls="--",lw=1.8)
ax.text(6,10.15,"~6 wk (patient): hard callus -> remodeling",ha="center",fontsize=9,weight="bold")

hits=[r for r in rows if int(r["in_formula"]) or int(r["in_rivaroxaban"])]
# lay out labels vertically to avoid overlap, grouped by placed week
import collections
buckets=collections.defaultdict(list)
for r in hits:
    wk=phase_week.get(r["phase"],5); buckets[wk].append(r)
for wk,rs in buckets.items():
    for i,r in enumerate(rs):
        y=7.5-i*1.15
        col=eff_color(r["effect_on_healing"])
        who="F" if int(r["in_formula"]) else ""
        who+="R" if int(r["in_rivaroxaban"]) else ""
        ax.scatter([wk],[y],s=120,color=col,zorder=5,edgecolors="white")
        ax.text(wk+0.18,y,f"{r['gene']} ({who})",fontsize=8.5,va="center")
ax.set_xlim(0,14.5); ax.set_ylim(0,10.6); ax.set_yticks([])
ax.set_xlabel("weeks after fracture")
# legend
from matplotlib.lines import Line2D
leg=[Line2D([0],[0],marker='o',color='w',markerfacecolor="#C0392B",markersize=11,label="potentially ADVERSE to healing"),
     Line2D([0],[0],marker='o',color='w',markerfacecolor="#2E7D5B",markersize=11,label="potentially FAVORABLE"),
     Line2D([0],[0],marker='o',color='w',markerfacecolor="#E0A63C",markersize=11,label="mixed / bidirectional"),
     Line2D([0],[0],marker='o',color='w',markerfacecolor="#888",markersize=9,label="neutral/minor")]
ax.legend(handles=leg,loc="lower right",fontsize=8.5,framealpha=0.9)
ax.set_title("Formula & rivaroxaban actions mapped onto fracture-healing phases\n(F=formula, R=rivaroxaban; direction from known pharmacology)",weight="bold",fontsize=12)
plt.tight_layout(); plt.savefig("results/figures/fig10_fracture_timeline.png",dpi=150); plt.close()
print("wrote results/figures/fig10_fracture_timeline.png")

# ---- Report ----
adverse=S["adverse"]; favor=S["favorable"]; mixed=S["mixed"]
def rget(g): return next(r for r in rows if r["gene"]==g)
L=[]; w=L.append
w("# 附：结合骨折愈合靶点的客观评估（腓骨骨折 ~6 周）")
w("")
w("> 数据驱动、方向标注的骨折愈合维度分析 · **客观呈现利弊，非用药建议**")
w("> 本文回答的是“方剂/利伐沙班对**骨愈合**这一维度可能是帮助还是妨碍”，而非“您是否该用”。")
w("")
w("---")
w("")
w("## 一句话客观结论")
w("")
w(f"> **对“骨折愈合”这一维度，计算显示是一幅“双刃、以潜在不利信号更确凿”的图景**——"
  f"方剂与 **{S['n_intersection']} 个骨愈合相关基因**重叠；在可判定方向的关键靶点中，"
  f"**潜在不利 {len(adverse)} 个、潜在有利 {len(favor)} 个、双向 {len(mixed)} 个**。")
w(">")
w("> **关键：您 6 周处于“骨痂/重塑期”，此期最依赖 血管新生(VEGFR2)、基质重塑(MMP)、COX-2 信号——"
  "而方剂恰好抑制这三者**，机制上属较确凿的“可能妨碍此阶段愈合”的信号；其“促成骨”效应(植物雌激素/Wnt)相对较弱且多为体外证据。")
w(">")
w("> **因此，从骨愈合角度，本组合并不能得出“帮助骨折痊愈”的结论；反而存在若干需要正视的潜在不利机制。**"
  "（但见下方“剂量与不确定性”——实际影响幅度可能受低暴露限制，须临床判断。）")
w("")
w("![骨折愈合阶段与靶点作用](../results/figures/fig10_fracture_timeline.png)")
w("")
w("---")
w("")
w("## 一、方剂/利伐沙班作用于骨愈合靶点（含方向）")
w("")
w("| 靶点 | 愈合阶段 | 谁作用 | 对骨愈合方向 | 机制说明 |")
w("|---|---|---|---|---|")
who_map=lambda r: ("方剂" if int(r["in_formula"]) else "")+("+利伐沙班" if int(r["in_rivaroxaban"]) else "")
phase_zh={"I":"炎症/血肿","K":"凝血支架","C":"软骨痂","A":"血管新生","O":"成骨/硬骨痂","R":"重塑","I/O":"炎症→成骨"}
for r in rows:
    if not (int(r["in_formula"]) or int(r["in_rivaroxaban"])): continue
    w(f"| **{r['gene']}** | {phase_zh.get(r['phase'],r['phase'])} | {who_map(r) or '—'} | {r['effect_on_healing']} | {r['note']} |")
w("")
w("---")
w("")
w("## 二、潜在不利信号（值得正视，机制较确凿）")
w("")
w("1. **COX-2 抑制（PTGS2）——证据最强的顾虑。** COX-2 是骨折愈合所必需；**NSAID 类 COX 抑制剂已被大量研究证实可延迟/影响骨折愈合**。方中 quercetin/apigenin 抑制 PTGS2 → 机制上同向。")
w("2. **抗血管新生（KDR/VEGFR2）。** 骨痂再血管化对愈合关键；黄酮(quercetin)是 VEGFR2 抑制剂 → 理论上不利于此期血运。")
w("3. **抑制基质金属蛋白酶（MMP2/9）。** MMP 参与骨痂重塑与血管化，MMP9 缺失动物骨愈合延迟；方中 isoliquiritigenin 强效抑制 → 理论上不利于**您当前的重塑期**。")
w("4. **抗凝对早期血肿（F2/F10）。** 凝血酶/纤维蛋白血肿是早期愈合支架；利伐沙班+方剂(水蛭素)抗凝理论上不利，但**6 周已过血肿期，此项现阶段影响较小**。")
w("")
w("## 三、潜在有利信号（较弱/多为体外）")
w("")
w("1. **植物雌激素样(ESR1)** 与 **Wnt/β-catenin(CTNNB1)**：黄酮(quercetin/kaempferol/baicalein)有促成骨、抗骨丢失的体外/动物报道。")
w("2. **抑制过度炎症(NLRP3)**：后期减轻慢性炎症或利愈合。")
w("3. **中医传统**：活血化瘀在骨折中期常用于改善局部血运、消肿——但 **大黄蛰虫丸为峻猛破血之剂，非典型接骨方**，是否对证须中医辨证。")
w("")
w("---")
w("")
w("## 四、剂量与不确定性（务必平衡看待）")
w("")
w("- 本项目 PBPK 已显示：半粒蜜丸递送的黄酮为 **µg 级、系统暴露极低**（对利伐沙班 AUC 仅 +3%）。"
  "**同理，这些黄酮对骨的系统性抗血管/抗 MMP/抗 COX 作用，实际幅度也可能被低暴露所限**——即“不利信号”未必转化为明显的临床骨愈合损害。")
w("- 网络药理学靶点表 **不含剂量与体内方向确证**；上述“方向”来自已知药理，属 **理论机制层面**，需实验/临床验证。")
w("- 因此客观表述是：**“机制上存在若干对重塑期不利的信号，且不支持‘促进骨折愈合’；但实际影响幅度不确定、可能有限。”** 既不夸大有害，也不宜宣称有益。")
w("")
w("---")
w("")
w("## 五、客观小结")
w("")
w("| 维度 | 客观判断 |")
w("|---|---|")
w("| 对骨折愈合“有帮助”? | **证据不支持**；促成骨信号较弱且体外为主 |")
w("| 对骨折愈合“有妨碍”? | **存在较确凿的机制性顾虑**（COX-2、VEGFR2、MMP 抑制，且正落在您的重塑期）|")
w("| 实际幅度 | **不确定**，可能受 µg 级低暴露限制 |")
w("| 与出血/前述结论的关系 | 叠加“出血风险增加(PD)”一并考量 → **审慎** |")
w("")
w("> **客观结论：** 结合骨折愈合靶点后，**没有数据支持“该组合能促进腓骨骨折痊愈”**；相反，在您所处的 **6 周重塑期**，"
  "方剂对 **COX-2、血管新生、MMP 重塑** 的抑制构成 **机制上需正视的潜在不利**（幅度不确定）。这一维度**进一步支持“需骨科/临床医师个体化评估、不宜自行联用”**的结论。")
w("")
w("> **免责声明：** 本节为计算性科研分析，**不构成诊断或治疗建议**。是否用药、如何促进骨折愈合，请以骨科医生与中医师面诊评估为准。")
w("")
open("docs/FRACTURE_HEALING_ASSESSMENT.md","w").write("\n".join(L))
print(f"wrote docs/FRACTURE_HEALING_ASSESSMENT.md ({len(L)} lines)")
