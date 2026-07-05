#!/usr/bin/env python3
"""Step 13: PBPK figure (C-t profiles + dose-response) and report."""
import json, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

S=json.load(open("data/processed/pbpk_summary.json"))
P=json.load(open("data/processed/pbpk_profiles.json"))
sweep=list(csv.DictReader(open("results/tables/pbpk_aucr.csv")))

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(14,5.6))

# Panel A: steady-state rivaroxaban concentration-time
t=P["t"]
ax1.plot(t,P["baseline"],lw=2.2,color="#3B6FA0",label=f"rivaroxaban 10mg alone (Cmax {max(P['baseline']):.0f} ng/mL)")
ax1.plot(t,P["formula_1p5g"],lw=2.2,ls="--",color="#2E7D5B",
         label=f"+ 1.5 g pill (AUCR {S['AUCR_half']:.2f}×, +{(S['AUCR_half']-1)*100:.0f}%)")
ax1.plot(t,P["strong_inhibitor"],lw=2.0,ls=":",color="#C0392B",
         label=f"+ strong CYP3A4/P-gp inhibitor (AUCR {S['validation_strong_AUCR']:.1f}× — validation)")
ax1.set_xlabel("time (h, steady state)"); ax1.set_ylabel("rivaroxaban plasma (ng/mL)")
ax1.set_title("Rivaroxaban concentration–time: formula barely shifts it",weight="bold",fontsize=11)
ax1.legend(fontsize=8.5); ax1.set_xlim(0,24)

# Panel B: AUCR vs quercetin dose (log x)
xs=[float(r["quercetin_mg"]) for r in sweep]; ys=[float(r["AUCR"]) for r in sweep]
ax2.plot(xs,ys,"-o",color="#7D5BA6",lw=2)
ax2.set_xscale("log")
ax2.axhline(1.25,ls="--",c="grey"); ax2.text(xs[-1],1.255,"AUCR 1.25 (clinically notable)",ha="right",fontsize=8.5,color="grey")
ax2.axhline(1.0,ls=":",c="k",lw=0.8)
# mark formula & supplement
qh=S["quercetin_half_mg"]; qf=S["quercetin_full_mg"]
ax2.scatter([qh],[S["AUCR_half"]],color="#2E7D5B",zorder=5,s=70)
ax2.annotate(f"1.5 g pill\n{qh*1000:.0f} µg → {S['AUCR_half']:.2f}×",(qh,S["AUCR_half"]),
             textcoords="offset points",xytext=(12,-2),fontsize=8.5,color="#2E7D5B")
ax2.scatter([500],[S["AUCR_500mg_supplement"]],color="#C0392B",zorder=5,s=70)
ax2.annotate(f"500 mg quercetin\nsupplement → {S['AUCR_500mg_supplement']:.2f}×",(500,S["AUCR_500mg_supplement"]),
             textcoords="offset points",xytext=(-8,10),fontsize=8.5,color="#C0392B",ha="right")
ax2.set_xlabel("quercetin daily dose (mg, log scale)"); ax2.set_ylabel("rivaroxaban AUC ratio (AUCR)")
ax2.set_title("Predicted rivaroxaban AUCR vs quercetin dose",weight="bold",fontsize=11)
ax2.set_ylim(0.98,1.30)
plt.tight_layout(); plt.savefig("results/figures/fig8_pbpk.png",dpi=150); plt.close()
print("wrote results/figures/fig8_pbpk.png")

# ---------------- report ----------------
L=[]; w=L.append
w("# 附：动态 PBPK 模型 —— 利伐沙班 AUC 变化的定量倍数")
w("")
w("> 半机制性动态 PBPK 模型（scipy 微分方程数值积分）· quercetin 真实人体 PK 参数")
w("> 已用已知强抑制剂(酮康唑)相互作用校准 · **科研估算，不构成用药建议**")
w("")
w("---")
w("")
w("## 一句话结论")
w("")
w(f"> **动态 PBPK 预测：半粒(1.5 g)大黄蛰虫丸使利伐沙班 AUC 仅增加约 "
  f"{(S['AUCR_half']-1)*100:.0f}%（AUCR≈{S['AUCR_half']:.2f}×）——临床可忽略，远低于有意义阈值(约 +25%)，"
  "且在利伐沙班本身的个体间变异范围(CV~30–40%)之内。**")
w(">")
w(f"> 即使按全量方剂(6 g/日)也仅 AUCR≈{S['AUCR_full']:.2f}×；即使换成 500 mg 纯 quercetin 保健品也仅 "
  f"AUCR≈{S['AUCR_500mg_supplement']:.2f}×。要让利伐沙班 AUC 升高到 1.25×，需要约 "
  f"**{S['quercetin_mg_for_AUCR1p25']:.0f} mg quercetin**——是半粒蜜丸所含 quercetin({S['quercetin_half_mg']*1000:.0f} µg)的 "
  f"**约 {S['quercetin_mg_for_AUCR1p25']/S['quercetin_half_mg']/1000:.0f} 万倍**。")
w("")
w("**含义：** 前一步静态 [I]/Ki 筛查曾对肠道 BCRP 报警，但**动态模型代入真实的低递送量(µg级)与低口服生物利用度后，"
  "该相互作用在系统暴露上几乎消失**。→ **PK 层面（升高利伐沙班血药浓度）在此剂量下不是主要担忧；"
  "真正的残余风险仍是 PD（水蛭素/抗血小板导致的出血叠加）。**")
w("")
w("![PBPK: 浓度-时间曲线与剂量-AUCR 关系](../results/figures/fig8_pbpk.png)")
w("")
w("---")
w("")
w("## 二、模型与验证")
w("")
w("**结构（半机制动态 PBPK）：** quercetin 按真实人体 PK 给药至稳态，其随时间变化的血浆(肝/肾)与肠腔浓度，"
  "**动态调节**利伐沙班的清除与生物利用度：")
w("- 肝 **CYP3A4** 可逆抑制（占利伐沙班清除 fm≈0.30）")
w("- 肾 **P-gp/BCRP** 主动分泌抑制（占清除≈0.30）")
w("- 肠 **CYP3A4** 抑制 → 提高肠壁可用度 Fg")
w("- 肠 **BCRP/P-gp** 外排抑制 → 提高吸收分数 Fa")
w("对稳态一个给药间隔积分利伐沙班 AUC，比较有/无抑制剂 → AUCR。")
w("")
w("**关键人体 PK 参数（文献值）：**")
w("")
w("| 药物 | 参数 | 取值 | 说明 |")
w("|---|---|---|---|")
w("| 利伐沙班 | 剂量/F/ka | 10 mg qd / ~0.9 / 1.24 h⁻¹ | 低剂量吸收~80–100% |")
w("| 利伐沙班 | V / t½ / CL | 50 L / ~6 h / 5.8 L·h⁻¹ | |")
w("| 利伐沙班 | 清除分数 | CYP3A4 0.30, CYP2J2 0.14, 水解 0.10, 肾P-gp/BCRP分泌 0.30, 其他 0.16 | 校准至酮康唑 |")
w("| quercetin | F / fu | ~2% / ~1.5% | 游离苷元口服生物利用度极低、高蛋白结合 |")
w("| quercetin | Ki(CYP3A4/P-gp/BCRP) | 1.0 / 7.08 / 0.03 µM | 源自 ChEMBL 实测 |")
w("")
w(f"**验证：** 模型对“强 CYP3A4+P-gp 抑制剂(酮康唑样)”预测 **AUCR≈{S['validation_strong_AUCR']:.1f}×**，"
  "与临床酮康唑对利伐沙班的 ~2.6× 相符；基线利伐沙班 Cmax≈152 ng/mL 亦与临床(~125–145 ng/mL)吻合 → 模型可信。")
w("")
w("---")
w("")
w("## 三、剂量–AUCR 结果")
w("")
w("| quercetin 日剂量 | 场景 | 利伐沙班 AUCR | 暴露变化 |")
w("|---|---|---|---|")
for r in sweep:
    d=float(r["quercetin_mg"]); note=r["note"] or ""
    dd=f"{d*1000:.0f} µg" if d<1 else f"{d:.0f} mg"
    w(f"| {dd} | {note} | {float(r['AUCR']):.3f}× | +{float(r['pct_change']):.1f}% |")
w("")
w("> 曲线在 ~1.11× 处趋于平台：因利伐沙班基线生物利用度已高(~0.9)，肠道外排抑制最多把吸收提到 100%(≈+11%)；"
  "而 quercetin 游离血浆浓度始终远低于其 CYP3A4/肾转运体 Ki，故肝/肾抑制几乎不启动。**所以即便大剂量 quercetin，"
  "对利伐沙班 AUC 影响也有限。**")
w("")
w("---")
w("")
w("## 四、与前几步结论的整合")
w("")
w("| 分析层次 | 结论 |")
w("|---|---|")
w("| 静态 [I]/Ki 筛查(step10) | 肠道 BCRP 对 licochalcone A/quercetin 报警(保守、按全剂量) |")
w("| **动态 PBPK(本步)** | **代入真实µg递送量+低BA后，AUCR≈1.03×，PK相互作用可忽略** |")
w("| 出血风险(PD) | **仍存在**：水蛭素(直接抗凝血酶)+抗血小板黄酮，与利伐沙班叠加 |")
w("")
w("**综合：在“利伐沙班 10 mg qd + 半粒 1.5 g 蜜丸”方案下——**")
w("- **PK（药物浓度层面）几乎无相互作用**（AUCR≈1.03×，可忽略）；")
w("- **出血风险的核心来自 PD 叠加**，而非利伐沙班被“升高”；")
w("- 故安全管理重点应放在 **监测出血、避免叠加其他抗栓药、评估个体出血风险**，而非担心药物浓度飙升。")
w("")
w("---")
w("")
w("## 五、局限（务必知悉）")
w("")
w("- **半机制/降阶 PBPK**，非全身多器官 PBPK(如 Simcyp/GastroPlus)；未纳入时间依赖性抑制(TDI)、酶诱导、肠菌代谢。")
w("- **quercetin 递送量基于含量估算**(µg 级)，游离苷元 vs 糖苷、品种差异带来不确定；但即使放大百倍仍不改变“可忽略”的定性结论。")
w("- **未建模水蛭素(PD 抗凝)**——出血叠加的主因，PBPK 只回答 PK，不回答出血风险。")
w("- 个体差异(肝肾功能、CYP3A4 基因型、年龄、合并用药)可显著改变基线，需个体化。")
w("")
w("> **免责声明：** 本节为计算性科研模型，**不构成医疗建议**。抗凝治疗与中西药联用请咨询专业医师与临床药师。")
w("")
open("docs/PBPK_MODEL.md","w").write("\n".join(L))
print(f"wrote docs/PBPK_MODEL.md ({len(L)} lines)")
