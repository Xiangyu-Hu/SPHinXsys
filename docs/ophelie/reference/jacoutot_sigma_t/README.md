# Jacoutot σ(T) digitization prep（2026-08-11）

## 重要更正：CEP “Fig. 2” 不是 σ(T)

Jacoutot et al., *Chem. Eng. Process.* 47 (2008) / CES 搅拌论文里的 **Figure 2** 是 **EM–热流耦合算法流程图**，不是电导率曲线。

真正的 **σ(T)** 曲线在 **Jacoutot 博士论文**（HAL `tel-00121350`）：

| 图 | 公式 | 适用 |
|----|------|------|
| **III.12** | \(\log_{10}\sigma = 3.7921 - 3179.8/T\) (III.6) [C_2] | 静态 / 自然对流（French natural 主线） |
| **IV.10** | \(\log_{10}\sigma = 4.05726 - 3923.73/T\) [C_8] Uox2+RuO₂ | 机械搅拌工况 |

附录还写明 \(\sigma = 10^{3.7921-3179.8/T}\)，因此 **log = log10**。

## 为何不“扣像素点”

论文图是由上述拟合律画出的；直接用公式采样比 WebPlotDigitizer 像素扣点更可追溯、无读图误差。  
仓库内 CSV = 公式在离散温度上的采样点（便于审计 / 对照图）。

## 与 CEP Table 1 的差异（必须知情）

| 来源 | σ(1473 K) [S/m] |
|------|-----------------|
| CEP / CES **Table 1** 代表点 | **16**（范围 \(10^{-6}\)–30） |
| 论文公式 **III.6** | **≈43.0** |
| 论文公式 **IV.10** | **≈24.7** |
| III.6 给出 σ=16 的温度 | **≈1229 K** |

Table 1 的 16 是工况代表值；**不是** III.6 在 1473 K 的取值。  
Stage 3 以前用 `σ=16` 常量 / Arrhenius 锚定 16@1473，与论文 III.6 量级不同——换论文律后功率会变，属预期。

## 本目录文件

- `jacoutot_thesis_III12_sigma_T.csv` — 自然对流 / 静态主线
- `jacoutot_thesis_IV10_sigma_T.csv` — 搅拌 / Uox2
- 代码：`makeJacoutotThesisSigmaTemperatureLawIII12()` / `...IV10()`  
  `paper_digitized_sigma=1` 当选用论文律时

## 建议验收顺序

1. `--sigma-law=thesis-iii12` + 现有 σ(T) 外环（先 coil-only，再 +A_glass）  
2. 对照：`σ(1229K)≈16`，`σ(1473K)≈43`  
3. 观察 `P_joule` 相对旧 Arrhenius 的变化（应升高量级）  
4. 暂不把 Table1 的 16@1473 强行塞进 III.6

## HAL PDF

在线 thesis PDF 受 Anubis 防护，自动化下载常拿到 HTML challenge。  
人工打开：https://theses.hal.science/tel-00121350 （图 III.12 / IV.10 目视核对 CSV）。
