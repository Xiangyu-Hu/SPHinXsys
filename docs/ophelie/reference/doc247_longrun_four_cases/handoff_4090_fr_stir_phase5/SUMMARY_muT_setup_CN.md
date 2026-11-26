# Phase5 法国搅拌：变粘性（μ(T)）设置与结果摘要

**日期**：2026-11-26  
**口径**：ASSUMED 工程对照，不宣称文献 PASS  
**对应数据包**：同目录 `../phase5_fr_stir_muT10s_pack_20261126.tgz`  
**原始输出**：`../phase5_fr_stir_muT10s/{MU_const_rpm10,MU_varT_rpm10}/`

---

## 1. 算例目的

在**同一套搅拌基线**上，只改热流变开关，比较：

| 算例文件夹（包内） | 原始 case | 开关 | 粘性含义 |
|--------------------|-----------|------|----------|
| `01_MU_const_rpm10_N10rpm_mu4_t10s` | `MU_const_rpm10` | `--constant-thermal-rheology` | \(\mu\equiv 4\,\mathrm{Pa\cdot s}\)（钉在 \(T_0=1473\,\mathrm{K}\)） |
| `02_MU_varT_rpm10_N10rpm_muT_t10s` | `MU_varT_rpm10` | `--temperature-dependent-properties` | \(\mu(T),\,c_p(T),\,k(T)\) 按 CES 2008 Fig.2 表 |

其余物理设置两边相同，与前期 10 rpm / \(\Delta z=0\) 基线一致。

---

## 2. 共同启动与边界（两边相同）

| 项 | 设置 |
|----|------|
| 温度初场 | 数字化 **Fig.8a**：`SEED_STIR_Fig8a_axisym_ASSUMED_T_rz.csv`；\(T\in[972,1458]\,\mathrm{K}\)，体积均值约 **1319 K** |
| 焦耳热形状 | **Fig.10 Q**：`SEED_STIR_Fig10_Qshape_ASSUMED_rz.csv` + `--fig10-q-correct`，blend=0.8，every=1 |
| 功率 / EM | `--em-control=fixed-power`，`--target-absorbed-power=60000`；`--thermo-em-coupling=periodic`，`--em-update-interval=0.5` |
| 几何 / 桨 | `reload_HZ_z00`；`--rotation-rpm=10` |
| 热阻 | skull + free-crust（delta=0.005 m，free-crust-k=1） |
| 物理时间 | `--end-time=10`（约 1.67 转） |
| 弱可压 | \(c_0=8\,\mathrm{m/s}\)，默认 `acoustic_cfl=0.55` |

说明：**Fig.10 不是温度分布**，只是焦耳热 \(Q\) 的形状校正；温度 IC 是 Fig.8a。

---

## 3. 粘性如何随温度变化

### 3.1 不是“换一个常数系数”

开启 `--temperature-dependent-properties` 后，启用的是 Jacoutot et al., *Chem. Eng. Sci.* 63 (2008) **Fig.2 NUM** 的数字化表（代码中 OPEN_PARAMETER），**不是**作者给出的闭式公式。

- 对 \(\mu\)：节点间做 \(\log_{10}\mu\)–\(T\) 插值（`logarithmic_value_=true`）
- 锚点：Table-1，\(T=1473\,\mathrm{K}\)，\(\mu=4\,\mathrm{Pa\cdot s}\)
- 同时启用：\(c_p(T)\)、\(k(T)\)（同 Fig.2 数字化表）
- \(\sigma(T)\)：周期耦合里本来就会更新（常粘 / 变粘两边都有）

实现位置：`electromagnetic_ophelie_french_material_laws.h` 中  
`makeJacoutotCes2008Fig2ViscosityTemperatureLaw()` / `installJacoutotCes2008Fig2ThermalRheology()`。

### 3.2 数字化节点表 \(\mu(T)\) [Pa·s]

| \(T\) [K] | 300 | 550 | 800 | 1050 | 1200 | 1300 | 1400 | **1473** | 1550 | 1650 | 1800 |
|-----------|-----|-----|-----|------|------|------|------|----------|------|------|------|
| \(\mu\) | \(10^5\) | \(10^5\) | \(8\times10^4\) | 2000 | 80 | 20 | 8 | **4** | 2.5 | 1.5 | 1 |

### 3.3 与 Fig.8a 温场相关的插值示例

| \(T\) [K] | 972 | 1085 | 1200 | 1319 | 1400 | 1473 |
|-----------|-----|------|------|------|------|------|
| \(\mu\) [Pa·s] | \(\sim 6.3\times10^3\) | \(\sim 9.4\times10^2\) | 80 | \(\sim 17\) | 8 | 4 |

趋势：冷区极稠、热区变稀；相对常粘 \(\mu=4\)，偏冷区域有效粘性可高 1–3 个数量级。

---

## 4. \(\mu\) 如何进入动量方程（实现口径）

开启变物性后：

1. 对每个流体粒子按局部 \(T_i\) 查表得到 \(\mu_i(T_i)\)（粒子场会写）；
2. 水力粘性力使用的是 **体积平均**
   \[
   \mu_\mathrm{ref}=\langle\mu(T)\rangle_V
   \]
   再 `setReferenceViscosity`；
3. **不是**粒子对局部 \(\mu_{ij}\)。

日志对应：

```text
viscosity_mode=mean_mu_T_volume_average
property_modes: ... viscosity=mean_mu_T_volume_average
```

因此会出现：

- \(\mu(T_\mathrm{mean})\sim\mu(1320\,\mathrm{K})\sim 17\,\mathrm{Pa\cdot s}\)
- 但冷区 \(\mu\) 极大，把体积平均抬到本算例 **\(\mu_\mathrm{mean}\approx 224\,\mathrm{Pa\cdot s}\)**

即 \(\langle\mu(T)\rangle\gg\mu(\langle T\rangle)\)，来自 \(\mu(T)\) 强凸性，不是初场与基线不一致。

---

## 5. 本批数值结果（t ≈ 10 s）

| 量 | 常粘 `MU_const_rpm10` | 变粘 `MU_varT_rpm10` |
|----|----------------------|---------------------|
| \(\mu_\mathrm{min}/\mathrm{mean}/\mathrm{max}\) | 4 / 4 / 4 | \(\approx 4.1\) / **224** / \(\approx 7.8\times10^3\) |
| \(T_\mathrm{mean}\) | \(\sim 1320\,\mathrm{K}\) | \(\sim 1319\,\mathrm{K}\) |
| \(U_\mathrm{max}\)（终点附近） | \(\sim 0.34\,\mathrm{m/s}\) | 过程中长期 \(\sim 15\,\mathrm{m/s}\)（\(>c_0=8\)） |
| 墙时 | \(\sim 8\,\mathrm{min}\) | \(\sim 54\,\mathrm{min}\) |
| 退出 | `exit=0`，`diverged=0` | 同左；`validation_passed=0` |

变粘更慢的主因：有效 \(\mu\) 大 + `U_max` 高 → 声学时间步变小。  
`U_max\gg c_0` 表示弱可压假设局部不健康；速度场宜标 ASSUMED，不宜当论文级流场。

`validation_passed=0` 的 fail 原因包括：`q_transfer`、`energy_residual`、`variable_cp_not_validatable`。

---

## 6. 报告可用表述（可直接粘贴）

> 变粘算例在与基线相同的 Fig.8a 温场初值、Fig.10 焦耳热形状校正与 60 kW 周期电磁耦合下，启用 Jacoutot CES 2008 Fig.2 数字化 \(\mu(T)/c_p(T)/k(T)\) 表（\(\log\mu\)–\(T\) 插值，1473 K 锚于 \(\mu=4\,\mathrm{Pa\cdot s}\)）。水力粘性取体积平均 \(\langle\mu(T)\rangle\) 而非局部粒子对粘性；因冷区贡献，10 s 末 \(\mu_\mathrm{mean}\approx 224\,\mathrm{Pa\cdot s}\)，约为常粘对照（\(\mu=4\)）的 56 倍。该 run 完成指定物理时间且未触发发散退出，但全程 \(U_\mathrm{max}\) 显著高于人工声速 \(c_0=8\,\mathrm{m/s}\)，结果仅作工程演示对照，不宣称文献 PASS。

---

## 7. 画图数据位置（包内）

| 用途 | 路径 |
|------|------|
| 中面温度 | 各 case 下 `french_T_rz.csv`（列：`r_m,z_m,T_K`） |
| 中面速度 / 流线 | `derived_midplane_velocity_rz/*_U_rz_midplane.csv`（`r_m,z_m,Ur_mps,Uz_mps,Utheta_mps,count`） |
| 时间序列 | `french_stirring_em_monitor.csv`、`french_energy_budget.csv`（含 `mu_min/mean/max`） |
| Fig.8a 文献等值线 | `00_prototype_literature_seeds/STIR_Fig8a_contour_points.csv` |
| 给画图助手的短说明 | `GPT_PLOT_BRIEF.md` |

---

## 8. 关键命令备忘

变粘单 case（GPU0）：

```bash
bash /home/yongchuan/sphinxsys/docs/ophelie/reference/doc247_longrun_four_cases/overnight_6h/run_phase5_muT10s_gpu0_varT_only.sh
```

核心差异仅一行：

```bash
# 常粘
--constant-thermal-rheology

# 变粘 / 变物性
--temperature-dependent-properties
```
