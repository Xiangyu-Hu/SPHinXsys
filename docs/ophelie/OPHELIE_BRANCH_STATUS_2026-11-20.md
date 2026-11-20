# OPHELIE 分支状态快照（2026-11-20）

> **分支：** `feature/EM-Biot-Savart`  
> **范围：** 法国自然对流 / 搅拌主线 + Doc201–227 工程修复 + 本轮仓库卫生  
> **不是文献定量复现声明。** 本地 round pack / `doc*_runs` 已不进 git；证据以本地包与本摘要为准。

---

## 1. 一句话结论

自然对流主线（`french_natural_*` + reduced-pressure）已形成可跑的 **EM → Q → 热 → 浮力 → WCSPH** 工程链，并完成 Doc221–227 一批关键稳定性修复（ForcePrior、真实时间步、局部 μ、周期 Couette、隐式黏性短程）；**Jacoutot 文献环流形态 / 绝对温度定量对比尚未通过。**  
搅拌主线（`french_stirring_*`）具备 **60 kW 标定 EM + 实验室系固定 Q 网格 + 运动学桨 + Boussinesq** 原型；**论文级耦合（σ(T)↔EM 迭代、冷坩埚/skull、Fig. 定量签字）仍未完成。**

---

## 2. 算例地图

| 算例 | 角色 | 成熟度（工程） |
|------|------|----------------|
| `test_3d_ophelie_french_natural_em` | Natural 几何 EM + 功率门 | 可跑；场精度仍开 |
| `test_3d_ophelie_french_natural_glass_relax` | Natural reload 网格 | 前置 |
| `test_3d_ophelie_french_natural_convection_frozen_q` | **自然对流主战场**：冻结/暖启 Q、热边、浮力、局部 μ、隐式 split | 工程链强；文献对比弱 |
| `test_3d_ophelie_french_natural_reduced_pressure` | reduced-pressure 生产路径；Doc180+ 全耦合入口 | 可跑；与 frozen_q 共用大量修复 |
| `test_3d_ophelie_french_complex_joule_to_heat_one_way` | σ(T)、单向热、Stage 3.x | 可跑 |
| `test_3d_ophelie_doc226_two_layer_couette` | 两层 / 常 μ Couette（SYCL+CK 周期 ghost） | Doc227 D 门控载体 |
| `test_3d_ophelie_french_stirring_glass_relax` | 搅拌三体 reload | 前置 |
| `test_3d_ophelie_french_stirring_em` | **搅拌主战场**：实验室系 Q 网格 + 桨 | 原型可跑；文献 Gate C 开 |
| `test_3d_ophelie_rh200_glass_em_stirring` | RH200 中期演示链 | 演示级 |
| `test_3d_ophelie_analytic_cylinder_skin` | EM 门控 / 解析柱 | Doc199–218 gates |

标准输入（进 git）：`reference/jacoutot_sigma_t/`、`reference/standard_data/`（Fig.7/8 PNG+CSV、Fig.8 ASSUMED T）。

---

## 3. 自然对流：已实现

### 3.1 物理 / 数值能力

- 复边通量 EM + φ（含 variational / P1 / D2 可选路径；文件已改名为 `phi_variational` / `phi_p1_projection` / `phi_d2_local_poly`）
- 焦耳 → 热：冻结粒子 Q、空间固定 Q、暖启 / Fig.8 ASSUMED T CSV
- 热扩散 + 法国式 Robin / 辐射边；顶部滑移 / paper-top 边界试验路径
- Boussinesq 浮力；**ForcePrior 生命周期**修复（不再错误清共享总力）
- 常 μ / **ParticleViscosity 局部 μ**；与 `Viscosity` 退化对照
- **ACTUAL_TIMESTEP** 约束调度（避免审计步长与真实推进脱节）
- reduced-pressure + **implicit_split** 隐式黏性（Doc227 F：可越过显式局部 μ 的早期 abort）
- Couette：CK 可读周期 ghost + 常 μ / 两层 μ 门控入口

### 3.2 Doc227 轮次门控（本地证据，包未入库）

| 门 | 结果（工程） |
|----|----------------|
| A 固定证据 | PASS |
| B ACTUAL_TIMESTEP | PASS |
| C 均匀 μ 退化 | PASS |
| D 真周期 Couette（邻域 + 常 μ + μ4/40） | PASS |
| E 法国局部 μ 短程 | 短程 PASS；**E4 ~5 s 在 t≈2.69 触 U_max 门 FAIL**；E4b 常 μ 3 s PASS |
| F 隐式局部 μ（reduced_pressure / implicit_split） | F0/F0b/F1/F1b PASS；F1b 3 s `U_max~0.0013`（对照 E4 ~0.052） |
| G 验收关闭 | 本轮工程验收关闭；**非文献复现** |

### 3.3 明确未完成（自然对流）

1. **文献环流 / Fig.7–8 定量对齐**（形态 + 速度尺度 + 温度绝对量）  
2. 去掉 `DIAGNOSTIC_CAP` / 固定温捷径后的长时稳态  
3. 全耦合：σ(T) → 周期 EM 更新 → Q 反馈（论文式迭代）  
4. 变 μ 显式路径在法国高黏度下的长时稳定性（隐式是缓解，不是文献签字）  
5. 方腔等独立热/浮力定量债（历史 VALIDATION_DEBT 仍在）

---

## 4. 搅拌：已实现

### 4.1 能力

- EM 预解 + 吸收功率标定（搅拌约定 ~60 kW）+ Picard 自感
- **实验室系欧拉 Q 网格**，每步插回运动粒子（相对 frozen_q 的关键差别）
- 运动学桨驱动（已自 Simbody 收敛到可记录 tip speed / 载荷的路径）
- Boussinesq + 法国热边约定；σ(T) 表可用 Jacoutot 数字化 CSV
- 周期热–EM 耦合基础设施（早期 Doc 层；搅拌论文 Fig.4 级迭代未签字）

### 4.2 明确未完成（搅拌）

1. Jacoutot CES 2008 **Gate C**：Fig. 温场 / 流场 / Np / Nu 定量  
2. 冷坩埚分段壁、skull、金属载体电流等论文几何  
3. 与 FLUENT↔OPHELIE 同构的 **σ↔EM 强迭代** 与长时（~1000 s）热平衡  
4. 搅拌工况下局部 μ / 隐式黏性与 Doc227 修复的系统并轨验证（主修复战场在 natural）

详见历史快照口径：[`09_FRENCH_LITERATURE_REPRODUCTION_STATUS_2026-09.md`](09_FRENCH_LITERATURE_REPRODUCTION_STATUS_2026-09.md)（日期较早，结论方向仍成立：原型强、定量弱）。

---

## 5. 共享基础设施（两 case 共用）

| 模块 | 状态 |
|------|------|
| Complex edge-flux + shared EMF / skin | Doc212–213 起生产路径可用；场 L2 仍 PARTIAL 级 |
| φ variational / P1 / D2 | 生产头已能力命名；CLI 仍保留历史 `doc21x_*` 开关名 |
| Analytic cylinder gates | 收入 `diagnostics/gates/` + umbrella |
| CK 时间步 / acoustic / solver 钩子 | Doc227 ACTUAL_TIMESTEP 相关 `src/shared_ck` 改动在本推送中 |
| 文档 | Docs 001–227 已并入 `archive/DOCS_*_MERGED.md`；单独 201–227 计划已删 |
| 运行产物 | `DOC*_ROUND_PACK*`、`doc*_runs/`、`packs/` **gitignore + 工作区已删** |

---

## 6. 建议下一步（按优先级）

1. **自然对流：** 在 F 隐式路径上放开诊断帽，做有文献可比窗口的长时跑；先形态门，再绝对 T / U。  
2. **Couette / 局部 μ：** 将 Doc227 D/F 结论固化为回归命令（CI 或脚本），防止回退。  
3. **搅拌：** 把 Doc227 黏性 / ForcePrior / 真实时间步修复并入 `french_stirring_em` 短回归，再谈论文 Fig.  
4. **文献 Gate：** 固定比较表（Fig.7/8 CSV 已在 `standard_data/`），禁止再改口径凑 PASS。

---

## 7. 本轮 git 推送内容（卫生 + 代码）

见同目录推送说明；要点：

- 合并归档 Docs 001–250；删除单独 201–213 追踪计划  
- `standard_data/`、welding 文本计划；忽略 packs/runs/worktrees  
- Doc214–227 相关代码、Couette 算例、诊断头重命名/分门  
- 本状态文档
