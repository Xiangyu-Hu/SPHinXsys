# Doc213 — LegacyFlux／edge EMF 趋肤场闭合（Cursor 计划）

日期：2026-11-07  
前置：Doc212（相位／边门 PASS；场 L2 ~5–6% PARTIAL；F5／F6 证明换 RHS 救场砸边）

## 1. 目标

在**生产路径不变**的前提下（shared complex `phi_rhs=legacy-flux` EMF；LHS=`legacy-pairwise`）：

1. 把圆柱弱趋肤中部 `L2_rel_J_mid` 从 ~5–6% **明显压下**（目标 ≤2%，理想 ≤1%），且  
2. **边门** `edge_rel_l2≤1%` 与 **相位 A1** 不回退。

禁止：把 `--doc212-complex-rhs=div-sigma-a|legacy-pairwise` 升为默认；禁止回退实虚分裂。

## 2. 已钉死事实（继承 Doc212）

| 事实 | 证据 |
|---|---|
| 场债在 LegacyFlux／EMF 离散，非 Picard／ghost／phi_rhs_ghost | `doc212_runs/E_power/EB_FIELD_ROOT_CAUSE.md` |
| 统一 DivΣA：场好、边 ~100% FAIL | F5 |
| 统一 pairwise：同交易 | F6 |
| 生产 shared EMF：边 ~0.23% PASS、场 ~5.5% | EA1 |

## 3. 工作顺序

| 步 | 内容 | 门 |
|---|---|---|
| **G0** | 基线三元组：no-self 生产／DivΣA；+ 生产 Picard+e1；粒子 CSV 含 `rhs_leg−rhs_div` | 复现 Doc212 数字 |
| **G1** | 径向／区域分解场债（rim vs interior；`rhs_leg_minus_div` vs \|J\| 误差相关） | 定位 EMF 哪一段 |
| **G2** | **只改**共享 EMF 装配（权重／A 插值／近壁通量），每刀同时报 edge+L2+A1 | 候选修复 |
| **G3** | 若 EMF 内无法兼得：设计**边一致**新 RHS（须过 e1），再考虑替换生产 | 仅当 G2 穷尽 |
| **P1 并行可选** | E-C 薄体绝缘短测（不阻塞 G0–G2） | EC 薄板 |
| **暂缓** | 法国长回归宣称、TEAM7、SPH-ε 再压 1% | — |

## 4. 运行约定

- GPU 命令由用户执行；**短期跑不进 tmux**。  
- 结果落 `docs/ophelie/reference/doc213_runs/`（gitignored `doc*_runs/`）。  
- 每批跑完后评审再改代码；不自动开法国／TEAM7。

## 5. 验收

- 生产 profile 不变旗下：`L2_rel_J_mid` 改善且 `edge_rel≤1%`、A1 PASS。  
- `ROUND213_REVIEW.md` + `ROOT_CAUSE_AND_FIX.md` 更新链。  
- 未达标不得改法国默认。
