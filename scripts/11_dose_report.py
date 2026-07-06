#!/usr/bin/env python3
"""Step 11: figure + report for the rivaroxaban-10mg + 1.5 g-pill regimen."""
import csv, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

M=json.load(open("data/processed/dose_ddi_model.json"))
half={r["constituent"]:r for r in M["half"]}
full={r["constituent"]:r for r in M["full"]}
order=[c for c in half if isinstance(half[c].get("ABCG2/BCRP_Igut/Ki"),(int,float))]
order=sorted(order,key=lambda c:-full[c]["ABCG2/BCRP_Igut/Ki"])

# ---------- Figure: intestinal BCRP Igut/Ki, half vs full, vs threshold ----------
fig,ax=plt.subplots(figsize=(9,5.5))
x=np.arange(len(order)); wdt=0.38
hv=[half[c]["ABCG2/BCRP_Igut/Ki"] for c in order]
fv=[full[c]["ABCG2/BCRP_Igut/Ki"] for c in order]
ax.bar(x-wdt/2,fv,wdt,label="6 g pill/day (full adult dose)",color="#C0392B",alpha=0.9)
ax.bar(x+wdt/2,hv,wdt,label="1.5 g pill (half pill — evaluated regimen)",color="#2E7D5B",alpha=0.9)
ax.axhline(10,ls="--",c="k",lw=1.2)
ax.text(len(order)-0.5,10.6,"FDA screen threshold = 10",ha="right",fontsize=9)
ax.set_yscale("log")
ax.set_xticks(x); ax.set_xticklabels(order,rotation=30,ha="right")
ax.set_ylabel("Intestinal [I]gut / Ki  (BCRP/ABCG2, log scale)")
ax.set_title("Rivaroxaban efflux-transporter (BCRP) interaction vs formula dose\nlower dose drops most constituents below the screening threshold",weight="bold",fontsize=11)
ax.legend(fontsize=9)
plt.tight_layout(); plt.savefig("results/figures/fig7_dose_bcrp.png",dpi=150); plt.close()
print("wrote results/figures/fig7_dose_bcrp.png")

# ---------- licorice / glycyrrhizic acid safety at this dose ----------
CRUDE=M["crude_in_pill"]
gancao_g=CRUDE*(90/1320); ga_mg=gancao_g*1000*0.03   # ~3% glycyrrhizic acid
gancao_full=(6/2.1)*(90/1320); ga_full=gancao_full*1000*0.03

# ---------- Report ----------
def flags(rowset):
    out=[]
    for c,r in rowset.items():
        for t in ["CYP3A4","ABCB1/P-gp","ABCG2/BCRP"]:
            g=r.get(f"{t}_Igut/Ki")
            if isinstance(g,(int,float)) and g>=10: out.append((c,t,g))
    return out
fh=flags(half); ff=flags(full)

L=[]; w=L.append
w("# 附：特定方案定量评估 —— 利伐沙班 10 mg qd + 半粒大黄蛰虫丸(1.5 g 蜜丸)")
w("")
w("> 静态机制性 DDI 剂量模型（FDA 2020 [I]/Ki 法）· ChEMBL 实测 Ki + 药典配比 + 文献含量估算")
w("> **本节为科研定量估算，不构成用药建议；存在较大不确定度。**")
w("")
w("---")
w("")
w("## 直接回答：会“好些”吗？")
w("")
w("> **会明显好一些——但“更低风险”不等于“安全/可自行联用”。**")
w(">")
w(f"> 把方剂降到 **1.5 g 蜜丸（≈{CRUDE:.2f} g 生药当量，约为常规日量的 1/4）**、利伐沙班用 **10 mg qd（低强度预防量）**，"
  "在定量模型中带来三点实质性改善：")
w("> 1. **药代(PK)相互作用信号大幅下降**：常规全量下 quercetin、licochalcone A 对 **BCRP** 的肠道 [I]/Ki 超阈；"
  "**减量到 1.5 g 后 quercetin 跌破阈值**（从 11.5 → 3.0），系统性 CYP3A4/P-gp/BCRP 抑制全部可忽略，"
  "仅剩 **licochalcone A** 一个（且其含量高度不确定）仍在阈值上。")
w("> 2. **药效(PD)抗栓负荷更小**：低剂量方剂递送的抗凝/抗血小板活性成分更少，与利伐沙班的叠加更弱。")
w(f"> 3. **甘草毒性风险很低**：该剂量仅含甘草 ~{gancao_g*1000:.0f} mg、甘草酸 ~{ga_mg:.1f} mg，"
  f"远低于诱发假性醛固酮增多症的 ~100 mg/日阈值。")
w(">")
w("> **但仍有不可忽视的残余风险**：利伐沙班本身抗凝，方中 **水蛭(水蛭素)** 是直接凝血酶抑制剂（本小分子模型未覆盖），"
  "**出血风险的叠加方向依旧存在**，且缺乏人体数据。故属 **“风险显著降低但非归零”**，仍需医师监护。")
w("")
w("---")
w("")
w("## 一、剂量换算与暴露估算")
w("")
w("| 项目 | 数值 | 说明 |")
w("|---|---|---|")
w(f"| 评估丸重 | 1.5 g | 同仁堂大蜜丸 3 g 之半粒 |")
w(f"| 生药当量 | ~{CRUDE:.2f} g | 大蜜丸 生药:蜂蜜≈1:1.1 稀释后 |")
w(f"| 相当于常规日量 | ~1/4 | 常规成人约 6 g 蜜丸/日(3 g bid) |")
w(f"| 利伐沙班 | 10 mg qd | 低强度/预防量；吸收~80–100%(不受food影响)，基线出血风险低于 15/20 mg 治疗量 |")
w("")
w("肠道浓度按 FDA 法 [I]gut = 递送摩尔量 / 250 mL；系统未结合 Cmax 用低生物利用度黄酮的量级假设估算。")
w("")
w("---")
w("")
w("## 二、药代动力学(PK)：对利伐沙班清除通路的抑制 [I]/Ki")
w("")
w("利伐沙班经 **CYP3A4** 代谢、是 **P-gp/BCRP** 底物。下表为 **1.5 g 剂量** 下各成分的估算：")
w("")
w("| 成分 | 递送量(mg) | 肠道[I](µM) | CYP3A4 Igut/Ki | P-gp Igut/Ki | **BCRP Igut/Ki** |")
w("|---|---|---|---|---|---|")
for c in order:
    r=half[c]
    def g(t):
        v=r.get(f"{t}_Igut/Ki"); return f"{v}" if isinstance(v,(int,float)) else "—"
    bcrp=r.get("ABCG2/BCRP_Igut/Ki")
    bmark=f"**{bcrp}** ⚠️" if isinstance(bcrp,(int,float)) and bcrp>=10 else g("ABCG2/BCRP")
    w(f"| {c} | {r['dose_mg']:.3f} | {r['Igut_uM']:.2f} | {g('CYP3A4')} | {g('ABCB1/P-gp')} | {bmark} |")
w("")
w(f"- **超过 FDA 筛查阈值(≥10)者**：1.5 g 剂量下仅 **{', '.join(f'{c}({t})' for c,t,_ in fh) if fh else '无'}**；"
  f"6 g 全量下为 **{', '.join(f'{c}({t})' for c,t,_ in ff)}**。")
w("- 系统性抑制（1+Cmax_u/Ki）在两种剂量下均 ≈1.0，可忽略——说明 **主要相互作用发生在肠道 BCRP**，且随剂量线性缩放。")
w("")
w("![剂量对 BCRP 相互作用的影响](../results/figures/fig7_dose_bcrp.png)")
w("")
w("> **PK 小结：** 减量到 1.5 g 使 quercetin 等跌破阈值，PK 相互作用从“多点超阈”降为“仅 licochalcone A 一点、且高度不确定”。"
  "考虑到 licochalcone A 主要存在于胀果甘草、常规乌拉尔甘草含量极低，**该方案下 PK 层面升高利伐沙班暴露的风险很可能较小**。")
w("")
w("---")
w("")
w("## 三、药效动力学(PD)：出血风险叠加——模型的关键盲点")
w("")
w("**必须强调：** 上述 PK 模型 **未覆盖** 出血风险叠加中最直接的一环——")
w("- **水蛭(制水蛭)中的水蛭素(hirudin)** 是直接凝血酶抑制剂（多肽，不在小分子库中）。1.5 g 蜜丸约含制水蛭 "
  f"~{CRUDE*(60/1320)*1000:.0f} mg；虽经炮制活性下降，但与利伐沙班(抗 FXa)构成 **双重抗凝**。")
w("- 方中 **quercetin/baicalein 等抗血小板(COX/LOX)** 作用亦随成分递送叠加。")
w("- 利伐沙班 10 mg qd 抗凝作用确定存在。")
w("")
w("→ 因此，**即便 PK 相互作用被剂量削弱，PD 层面的抗凝-抗血小板叠加仍使净出血风险高于单用利伐沙班**；"
  "只是相较“全量方剂 + 高剂量利伐沙班”，此低剂量组合的叠加幅度更小。")
w("")
w("---")
w("")
w("## 四、毒理学与其他安全性（该剂量下的量化）")
w("")
w(f"- **甘草/假性醛固酮增多症**：本剂量甘草酸 ~{ga_mg:.1f} mg（全量约 {ga_full:.1f} mg），"
  "远低于 ~100 mg/日风险阈值 → **该剂量下低风险**（长期每日服用仍需留意累积）。")
w("- **大黄蒽醌(emodin 等)**：单次递送 ~0.16 mg，急性毒性风险低；但 **长期使用** 蒽醌类与肝/肾负担、"
  "结肠黑变病、emodin 遗传毒性信号相关，不宜久服。")
w("- **动物药(水蛭/虻虫/蛴螬/土鳖虫)**：破血逐瘀，**孕妇禁用**；过敏体质注意。")
w("- 利伐沙班清除依赖肝肾 → **肝肾功能不全者暴露升高**，与本方联用时风险放大。")
w("")
w("---")
w("")
w("## 五、结论")
w("")
w("| 维度 | 全量方剂 + 高剂量利伐沙班 | **1.5 g 蜜丸 + 利伐沙班 10 mg qd** |")
w("|---|---|---|")
w("| PK(升高暴露) | quercetin+licochalcone A 超阈 | **仅 licochalcone A(不确定)**，显著减轻 |")
w("| PD(出血叠加) | 强 | **较弱但仍存在**(含水蛭素) |")
w("| 甘草毒性 | 需留意 | **低** |")
w("| 总体 | 风险较高 | **风险明显降低，但非归零** |")
w("")
w("> **一句话结论：** 用 **利伐沙班 10 mg qd + 半粒(1.5 g)大黄蛰虫丸** 相比高剂量组合 **确实“好一些”**——"
  "PK 相互作用大幅减弱、甘草毒性低；但因 **水蛭素等导致的出血叠加方向不变**、且 **无人体证据**，"
  "**仍须在专科医师/临床药师监护下、密切监测出血的前提下使用，不能作为自我药疗的安全依据**。")
w("")
w("**监测建议**：出血征象(瘀斑/黑便/血尿/Hb)、肝肾功能、（长期服用则查）血钾血压；避免叠加 NSAID/抗血小板药；"
  "高出血风险者与孕妇不联用。")
w("")
w("---")
w("")
w("## 六、方法学局限（务必知悉）")
w("")
w("- **成分含量为文献/药典量级估算**，游离苷元 vs 糖苷、licochalcone A 的甘草品种差异等带来 **数量级不确定度**。")
w("- **静态 [I]/Ki 筛查模型 ≠ PBPK 定量预测**；[I]gut 假设 250 mL 溶出、未计溶解度/首过/肠菌代谢(如槲皮素糖苷经肠菌释放苷元)。")
w("- **未定量水蛭素(PD 抗凝)**——出血叠加的主因之一，模型系统性低估 PD 风险。")
w("- 复方成分间相互作用、实际制剂含量未纳入。")
w("")
w("> **免责声明：** 本节为计算性科研分析，**不构成医疗建议**。抗凝药与中成药联用请务必咨询专业医师与临床药师并遵循循证指南。")
w("")
open("docs/DOSE_REGIMEN_DDI.md","w").write("\n".join(L))
print(f"wrote docs/DOSE_REGIMEN_DDI.md ({len(L)} lines)")
print(f"glycyrrhizic acid this dose ~{ga_mg:.1f} mg (full ~{ga_full:.1f} mg); leech ~{CRUDE*(60/1320)*1000:.0f} mg")
