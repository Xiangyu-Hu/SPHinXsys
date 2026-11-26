# Doc247 四套长时算例 — 就绪状态

**日期**：2026-11-25  
**分支意图**：生产长时挂机（非论文 PASS 宣称）

## 结论

| # | Case | 可执行体 | 长时挂机命令 | 生产身份成熟度 | 说明 |
|---|------|----------|--------------|----------------|------|
| 1 | 法国文献自然对流（最接近版） | `test_3d_ophelie_french_natural_convection_frozen_q` | **有**（600 s） | Doc246 证据到 ~60 s；~80 s 曾 U 门中止 | **C0_BASE**（JAC + coldwall + skull/crust + 50 kW）。无 doc227；`--doc224-u-blowup-gate=50`。**勿** F8SEED 长时；P65 仅能量对照。 |
| 2 | 法国文献搅拌 | `test_3d_ophelie_french_stirring_em` | **有**（1200 s / 10 h） | 长时脚手架已有；论文 Np–Re 口径未收口 | 60 kW + 10 rpm + Level-1 skull/crust。`Np_sph` 为工程候选，非论文冻结。 |
| 3 | 甲方自然对流 | frozen_q + `--client-nc-*` | **有**（600 s 工程挂机） | 仓库原仅 1 s 接线 PASS | ASSUMED 15 kW；电源 30 kW ≠ 已知效率。长跑矩阵此前未锁，本包给出可挂命令。 |
| 4 | 甲方搅拌 | `client_stirring_f1_heated_elastic`（弹性对照） / **`client_stirring_rigid_heated`（C2r 刚体，推荐长跑）** | **有**（C2r 默认 7200 s / 48 h，可 `--target-T-K=1400`） | C2r：规定刚体转动 + 规定热；软生产身份 | 弹性子步墙钟过长；冲 ~1400 K 用 C2r。热源仍非 EM。勿当论文 PASS。 |

## 法国搅拌是否存在？

**存在。** 路径：

- 源码：`tests/.../test_3d_ophelie_french_stirring_em/`
- 二进制：`build/.../test_3d_ophelie_french_stirring_em/bin/`
- Reload：`.../test_3d_ophelie_french_stirring_glass_relax/bin/reload/`（含 `GlassBody` + `Rotor`）
- README 已给生产长时默认（1200 s / 10 h 墙钟）

## 长时相对短时的软改动（必须）

法国 NC 在 Doc246 中约 **~80 s** 曾因 `U_max>0.05` 与/或 `doc227_viscous_force_age` 数值中止。长时版：

1. **不加** `--doc227-actual-timestep`
2. **加** `--doc224-u-blowup-gate=50`（仅防早期声速尖峰误杀；非有限值 U/T 仍会中止）
3. **不加** Fig.8 ASSUMED 冷壳全域初场（μ~1e4–1e5，无法长时）
4. 法国 NC **无** `--max-wall-hours`；靠 `--end-time=` + 人工停；搅拌两边有墙钟预算

## 身份边界（勿混）

- **法 NC ≠ 法搅拌**：NC 无桨；搅拌有机械功与桨几何。
- **法 ≠ 甲**：几何、功率、σ、桨机构不同。
- **甲搅拌 F1**：`heat-power` = 规定吸收 W，**不是** Biot–Savart EM。
- 本包 **不宣称** 文献 Fig.7/8 PASS；只保证可挂机、身份可复现。

命令：见同目录 `RUN_COMMANDS.md`。
