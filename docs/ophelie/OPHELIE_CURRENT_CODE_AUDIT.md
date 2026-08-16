# OPHELIE 当前代码审计与法国论文复现基线

审计日期：2026-08-04  
适用范围：`tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/` 及 French reduced、RH200、TEAM7 测试。

## 1. 当前可复用电磁链路

```text
多匝丝状线圈
  -> Biot-Savart A_coil / B_coil
  -> edge-flux 或 div-grad φ 求解
  -> E / J 重构
  -> Q = 0.5 * sigma * (|E_real|² + |E_imag|²)
  -> 玻璃体积积分功率
  -> 线圈电流标定到目标玻璃吸收功率
```

| 职责 | 主实现 | 状态 |
|---|---|---|
| 丝状多匝线圈和 `A_coil/B_coil` | `electromagnetic_ophelie_multiloop_source.h`、`electromagnetic_ophelie_biot_savart.*` | 可复用 |
| `phi` LHS/RHS 与 GMRES | `electromagnetic_ophelie_phi.*`、`electromagnetic_ophelie_phi_gmres.h` | 可复用；存在 div-grad 与 edge-flux 两条路径 |
| edge electromotive drop、守恒残差 | `electromagnetic_ophelie_edge_flux.h` | 主候选路径 |
| 边通量 LS 重构 `E/J` 与复数 `Q` | `electromagnetic_ophelie_edge_flux.h` | 主候选路径；需完成一致性验收 |
| 标量场后处理与场缩放 | `electromagnetic_ophelie_postprocess.*` | 可复用 |
| French reduced 调度 | `electromagnetic_ophelie_french_literature.h` | 可复用，但不是完整法国装置 |

## 2. 现有路线的精确定义

### 2.1 edge-flux

`ComputeOphelieEdgeFluxPhiRhsFromASrcCK`、`ComputeOphelieEdgeFluxResidualCK` 和
`ReconstructOphelieEdgeFluxElectricCurrentCK` 使用相同的 pair conductance 与 edge drop：

```math
e_{ij}=(phi_j-phi_i)+s_A omega ((A_i+A_j)/2) dot (x_j-x_i),
\qquad q_{ij}=-C_{ij}e_{ij}.
```

对复数场，实、虚两条 chain 分别求解，主 `JouleHeat` 复制自
`JouleHeatEdgeReconComplex`。图离散功率经 2026-08 Stage 1 Phase A 改为体积一致形式
`P_graph = Σ_i (Σ_j ¼ C e²) V_i`，并增加无向 `j>i` 审计 `P_undirected`；二者仍为诊断，
物理标定功率仍是 `P_recon`。详见 `docs/ophelie/08_EDGE_FLUX_UNDIRECTED_POWER_CLOSURE.md`。
法国主线的验收必须以 `P_recon` 为准，并记录 `P_graph/P_recon` 与 `P_undirected/P_graph`。

### 2.2 `A_ind`

`ComputeOphelieGlassSelfInducedBiotSavartCK` 已在 device kernel 中实现：

```math
A_ind(x_i)=mu_0/(4 pi) sum_{j != i} J_j V_j / sqrt(|x_i-x_j|²+epsilon²).
```

它是直接 `O(N²)` 求和；`stage2/electromagnetic_ophelie_self_induction.h` 提供
`A_coil + A_ind -> phi/J -> A_ind` Picard 调度。该路径在 French README 中仍标为
**experimental**：外层场比较和松弛含 host 同步，尚没有低电导率极限、标度律和圆柱交叉
参考的完整验收，不能作为已复现 OPHELIE 的依据。

### 2.3 French literature mode

`--literature-mode` 关闭后处理 field scaling，并先按单位电流求解，再按
`sqrt(P_target/P_raw)` 标定线圈电流。它是正确的固定玻璃吸收功率策略。

注意：未显式传入 current form 时，`applyFrenchLiteratureMode()` 仍选择 `DivGrad`。
因此 `literature_passed=1` 仅表示当前 selected operator 的 reduced internal gate 通过；
它不是“edge-flux 已冻结”，更不是“完整 OPHELIE 或法国论文复现”。

## 3. 与两篇 Jacoutot 论文的差距

| 论文要求 | 当前状态 | 必须动作 |
|---|---|---|
| 玻璃体积积分 EM | 已有 reduced cylinder | 冻结 edge-flux 基线并做独立验证 |
| `A_total=A_coil+A_glass+A_metal` | 仅 `A_coil` 主线；`A_glass` experimental；无 `A_metal` | 先验收 `A_glass` Picard，再实现 surface carriers |
| 冷坩埚分段/底板薄趋肤层 | 未实现 | 作为完整 OPHELIE 后期阶段 |
| 温度相关 `mu/cp/k/sigma` | 常数参数 | 单一、可追溯物性模块 |
| 50 kW 自然对流外层耦合 | 无 | `T -> sigma(T) -> EM -> Q -> thermo` 控制器 |
| 壁/底换热、液面辐射/对流、skull | RH200 壁面绝热；French reduced 无热边界 | 论文热边界与 `mu > 1e4 Pa s` skull 诊断 |
| 3D 温度方位平均到轴对称 EM | 未实现 | `(r,z)` 投影器 |
| 轴对称 `Q(r,z)` 回写 3D 搅拌 | 未实现 | `(r,z)` 插值器 |
| 10 rpm、60 kW、约 1000 s 的搅拌工况 | RH200 为约 100 rpm、50 kW 演示 | 建立独立 French stirred case |

## 4. RH200 与论文 case 的边界

`test_3d_ophelie_rh200_glass_em_stirring` 可复用：

- SYCL 3D 单相流动、刚体桨、温度变量、固定 Eulerian Q/J 网格采样；
- 体积功率和温度审计输出；
- VTP 动画与 reload 基础设施。

它不得作为论文复现主 case：其转速约 100 rpm，壁面是热绝缘，物性和 `sigma` 为常数，
EM 只在运行开头计算一次；这些都与 Jacoutot 3D 搅拌耦合模型不一致。

## 5. Stage 0 可重复基线

从 `build/` 运行。首次使用必须先以同一 `dp` 生成 reload：

```bash
cmake --build . --target test_3d_ophelie_edge_flux_power_uniform_field \
  test_3d_ophelie_french_reduced -j"$(nproc)"

./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_edge_flux_power_uniform_field/bin/\
test_3d_ophelie_edge_flux_power_uniform_field

FRENCH=./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/\
test_3d_ophelie_french_reduced

"$FRENCH" --relax=1 --dp=0.02 --glass-radius=0.325 --glass-height=0.50 \
  --state_recording=0

"$FRENCH" --reload=1 --dp=0.02 --glass-radius=0.325 --glass-height=0.50 \
  --literature-mode --ophelie-current-form=edge-flux --ophelie-edge-flux-complex=1 \
  --ophelie-edge-flux-normalization-mode=solver-local --state_recording=0
```

每次基线记录：

```text
particle_count, dp, frequency, sigma min/max/mean,
phi solver residual, phi equation residual,
edge residual reduction, edge antisymmetry residual,
P_graph, P_recon, P_target, Q min/max/mean,
A/J norms, output field finite counts, runtime.
```

禁止将 TEAM7 的高导体 plate 指标作为低导电玻璃 baseline 的 pass/fail gate；
但 TEAM7 必须保留为归一化、边界重构和高导体退化的回归测试。

## 6. Stage 0 退出条件

1. Uniform-field edge-flux power test 通过；
2. French reduced edge-flux 命令连续三次的关键标量在预定浮点容差内一致；
3. 保存上述每次运行的 console 与 CSV；
4. 默认/explicit operator、diagnostic-only 字段和 experimental `A_ind` 不再混淆；
5. 未达成前，不开始将 `sigma(T)`、流动或金属 carrier 混入现有主求解器。
