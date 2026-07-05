#!/usr/bin/env python3
"""Step 9: rivaroxaban co-administration figure + report (Chinese)."""
import csv, json
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ddi=list(csv.DictReader(open("results/tables/rivaroxaban_ddi.csv")))
summ=json.load(open("data/processed/rivaroxaban_ddi_summary.json"))

# ---------- Figure: PK modulation heat-ish grouped bars ----------
PK_TARGETS=["CYP3A4","ABCB1/P-gp","ABCG2/BCRP"]
col={"CYP3A4":"#C0392B","ABCB1/P-gp":"#E08A3C","ABCG2/BCRP":"#7D5BA6"}
data=defaultdict(dict)   # compound -> target -> pchembl
for r in ddi:
    if r["target"] in PK_TARGETS:
        data[r["compound"]][r["target"]]=float(r["pchembl"])
# rank compounds by max pchembl across PK targets
comps=sorted(data.keys(), key=lambda c:-max(data[c].values()))
fig,ax=plt.subplots(figsize=(10,6.5))
y=np.arange(len(comps)); h=0.26
for i,t in enumerate(PK_TARGETS):
    vals=[data[c].get(t,0) for c in comps]
    ax.barh(y+(1-i)*h, vals, height=h, color=col[t], alpha=0.9, label=t)
ax.set_yticks(y); ax.set_yticklabels(comps)
ax.invert_yaxis()
ax.axvline(6,ls='--',c='grey',lw=1); ax.axvline(5,ls=':',c='grey',lw=1)
ax.text(6.02,len(comps)-0.5,'pChEMBL 6\n(strong,<1µM)',fontsize=8,color='grey')
ax.set_xlabel("pChEMBL (higher = stronger inhibition)")
ax.set_title("Formula constituents inhibiting rivaroxaban ADME proteins\n(CYP3A4 metabolism · P-gp/BCRP efflux transport)",weight='bold')
ax.legend(loc='lower right',fontsize=9)
plt.tight_layout(); plt.savefig("results/figures/fig6_rivaroxaban_pk.png",dpi=150); plt.close()
print("wrote results/figures/fig6_rivaroxaban_pk.png")

# ---------- Report ----------
def rows_for(t):
    xs=[r for r in ddi if r["target"]==t]
    xs.sort(key=lambda r:-float(r["pchembl"])); return xs
pd_overlap=summ["pd_coag_platelet_overlap"]

L=[]; w=L.append
w("# 附：大黄蛰虫丸 与 利伐沙班(Rivaroxaban) 联用分析")
w("")
w("> 药物相互作用(DDI)计算评估 · 数据源 ChEMBL 实测活性 + 利伐沙班已知 ADME 特性")
w("> **本节为科研机制分析，不构成用药建议；联用决策须由临床医师个体化判断。**")
w("")
w("---")
w("")
w("## 结论先行")
w("")
w("> **不建议在无专科医师监护下自行将大黄蛰虫丸与利伐沙班联用。** 计算证据显示二者存在 **双重相互作用风险**：")
w("> 1. **药效学(PD)——出血风险叠加**：利伐沙班抑制 Xa 因子，方剂命中 **凝血酶 F2、血小板 COX(PTGS1/2)、"
  "脂氧合酶 ALOX12/15**，并含 **水蛭素(直接抗凝血酶)**——二者作用于 **同一止血系统的不同环节**，抗凝作用叠加。")
w("> 2. **药代动力学(PK)——升高利伐沙班血药浓度**：利伐沙班经 **CYP3A4** 代谢、是 **P-gp/BCRP** 底物；"
  "方中 **quercetin、licochalcone A、apigenin、kaempferol、luteolin** 等对这些酶/转运体有实测抑制作用，"
  "尤其 **quercetin 同时抑制 CYP3A4 + P-gp + BCRP**——正是利伐沙班说明书警示会 **增加暴露量与出血风险** 的模式。")
w(">")
w("> **净效应：出血风险很可能升高。** 若临床确需联用，须严密监测出血、评估肝肾功能，并知悉 **目前无 RCT/临床药代研究证实其安全性**。")
w("")
w("---")
w("")
w("## 一、利伐沙班药理背景")
w("")
w("| 属性 | 内容 |")
w("|---|---|")
w("| 药理类别 | 口服直接 **Xa 因子(FXa/F10)** 抑制剂 (DOAC) |")
w("| 代谢 | ~1/2 经肝 **CYP3A4/3A5、CYP2J2** 及 CYP 非依赖水解；~1/3 经肾以原形排泄 |")
w("| 转运体 | **P-糖蛋白(P-gp/ABCB1)** 与 **BCRP(ABCG2)** 底物 |")
w("| 关键 DDI | **CYP3A4+P-gp 联合强抑制剂**(如唑类抗真菌、利托那韦)→暴露↑→出血；强诱导剂(利福平、圣约翰草)→暴露↓→血栓 |")
w("| 出血叠加 | 与其他抗凝/抗血小板/NSAID/影响止血的药物合用→出血风险↑ |")
w("")
w("---")
w("")
w("## 二、药效学(PD)相互作用：抗凝/抗血小板作用叠加")
w("")
w("利伐沙班的靶点是 **Xa 因子(F10)**。分析显示方中小分子 **未直接命中 F10**（非同一位点竞争），"
  "但命中凝血级联下游与血小板通路的多个靶点，与利伐沙班形成 **对止血系统的多点联合抑制**：")
w("")
w("| 方剂命中靶点 | 作用 | 与利伐沙班的关系 |")
w("|---|---|---|")
rolemap={"F2":"凝血酶——凝血级联终点","PTGS1":"COX-1——血小板聚集","PTGS2":"COX-2——血小板/炎症",
 "ALOX12":"12-脂氧合酶——血小板活化","ALOX15":"15-脂氧合酶——脂质炎症"}
for g in ["F2","PTGS1","PTGS2","ALOX12","ALOX15"]:
    if g in pd_overlap:
        w(f"| **{g}** | {rolemap.get(g,'')} | 利伐沙班抑制其上游 FXa；本方抑制下游/血小板 → **序贯叠加** |")
w("")
w("此外，方中 **水蛭(水蛭素 hirudin)** 是 **直接凝血酶抑制剂**，与利伐沙班(抗 FXa)构成 **双重抗凝**——"
  "类似“FXa 抑制剂 + 凝血酶抑制剂”联用，抗凝强度显著增强，**出血风险相应升高**。")
w("")
w("> **PD 小结：** 联用 = 利伐沙班(抗 FXa) + 本方(抗凝血酶/抗血小板/纤溶) → **抗凝-抗栓作用叠加，出血风险升高**。"
  "这与中医“破血逐瘀”方剂 + 抗凝西药联用需警惕出血的临床共识一致。")
w("")
w("---")
w("")
w("## 三、药代动力学(PK)相互作用：可能升高利伐沙班暴露量")
w("")
w("利伐沙班的清除依赖 **CYP3A4**(代谢) 与 **P-gp / BCRP**(外排转运)。下表为方中成分对这些蛋白的 **ChEMBL 实测抑制活性**"
  "（pChEMBL 越高抑制越强；≥6 强、5–6 中等、<5 弱）：")
w("")
for t,note in [("CYP3A4","主要代谢酶——抑制→利伐沙班代谢减慢→血药浓度↑"),
               ("ABCB1/P-gp","肠道/肾外排转运体——抑制→吸收↑/排泄↓→暴露↑"),
               ("ABCG2/BCRP","肠道/胆/肾外排转运体——抑制→暴露↑")]:
    rs=rows_for(t)
    w(f"**{t}**（{note}）")
    w("")
    if rs:
        w("| 成分 | pChEMBL | 强度 | 实测值 |")
        w("|---|---|---|---|")
        for r in rs[:8]:
            w(f"| {r['compound']} | {r['pchembl']} | {r['potency']} | {r['assay_type']}={r['value']}{r['units']} |")
    else:
        w("_(成分中无实测数据)_")
    w("")
w("**关键警示 —— 多通路联合抑制：**")
w("- **quercetin** 同时抑制 **CYP3A4(IC50≈2.1µM) + P-gp + BCRP(EC50≈30nM)** 三条清除通路；")
w("- **licochalcone A** 为 **强效 BCRP 抑制剂(IC50≈6nM)**；")
w("- apigenin/kaempferol/luteolin 亦同时触及 CYP3A4 与 P-gp/BCRP。")
w("")
w("即方剂整体呈现 **“CYP3A4 + P-gp + BCRP 多重抑制”** 特征——正是利伐沙班说明书列为 **增加暴露、需避免/慎用** 的相互作用类型。")
w("")
w("![方剂成分对利伐沙班 ADME 蛋白的抑制](../results/figures/fig6_rivaroxaban_pk.png)")
w("")
w("---")
w("")
w("## 四、其他安全性考量")
w("")
w("- **甘草(glycyrrhizin)**：长期/大量使用可致 **假性醛固酮增多症**（低血钾、高血压、水肿），并影响 CYP3A4/P-gp，"
  "与利伐沙班的合并用药与电解质监测需注意。")
w("- **大黄(蒽醌)**：具致泻作用，可改变胃肠转运与药物吸收，长期使用致电解质紊乱。")
w("- **孕妇禁用**：方剂破血逐瘀 + 利伐沙班均为妊娠禁忌/慎用。")
w("- **围手术期/高出血风险人群**（消化道溃疡、近期出血、严重肝肾功能不全、血小板减少）风险更高。")
w("")
w("---")
w("")
w("## 五、证据强度与局限")
w("")
w("- 本分析基于 **体外实测活性(ChEMBL)** 与利伐沙班 **已知 ADME 特性** 的 **机制推断**，"
  "**非临床药代/药效研究**，也 **未检索到该组合的 RCT 或人体 PK 相互作用数据**。")
w("- 体外抑制强度（尤其高通量 qHTS 的 Potency 值）**可能高估体内相关性**；多数黄酮 **口服生物利用度低**，"
  "实际体内 CYP3A4/P-gp 抑制程度可能弱于体外——但 **BCRP 强抑制(quercetin/licochalcone A)** 与 **PD 出血叠加** 仍值得警惕。")
w("- 复方为混合物，成分间相互作用、实际含量与暴露量未纳入定量 PBPK 模型。")
w("")
w("---")
w("")
w("## 六、结论与建议")
w("")
w("> **能否联用？——计算证据倾向于“不建议常规联用，须专科医师个体化评估与监护”。**")
w(">")
w("> - **机制上**：PD（抗凝叠加）+ PK（升高利伐沙班暴露）**双重放大出血风险**，方向一致、互相加强。")
w("> - **证据上**：缺乏人体安全性/药代数据，属 **高不确定性**。")
w("")
w("**若临床仍需联用，建议：**")
w("1. 由 **血液科/心内科/临床药师** 主导，个体化评估血栓与出血风险（如 HAS-BLED）；")
w("2. **密切监测出血征象**（牙龈/皮下瘀斑、黑便、血尿、血红蛋白），必要时查 **抗 FXa 活性**；")
w("3. 评估并监测 **肝肾功能**（利伐沙班清除依赖）与 **血钾/血压**（甘草）；")
w("4. 避免同时叠加 NSAID/抗血小板药；孕妇、高出血风险者 **不联用**；")
w("5. 记录、随访，遵循循证指南——**不作为自我药疗**。")
w("")
w("> **免责声明：** 本节为计算性科研分析，**不构成医疗建议**。抗凝治疗与中西药联用请务必咨询专业医师与临床药师。")
w("")
open("docs/RIVAROXABAN_DDI.md","w").write("\n".join(L))
print(f"wrote docs/RIVAROXABAN_DDI.md ({len(L)} lines)")
