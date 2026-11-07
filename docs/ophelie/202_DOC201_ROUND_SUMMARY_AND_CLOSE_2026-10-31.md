# Doc201 轮次总结与关闭

日期：2026-10-31。  
计划：[`201_OPHELIE_UNIFIED_COMPLEX_EM_AND_LIVE_BUOYANCY_CURSOR_PLAN_2026-10-30.md`](201_OPHELIE_UNIFIED_COMPLEX_EM_AND_LIVE_BUOYANCY_CURSOR_PLAN_2026-10-30.md)  
证据树：[`reference/doc201_runs/`](reference/doc201_runs/)（单层；无嵌套 PACKAGES）  
核对：[`reference/doc201_runs/DELIVERABLE_CHECKLIST.md`](reference/doc201_runs/DELIVERABLE_CHECKLIST.md)  
Brief：[`reference/doc201_runs/GPT_DISCUSSION_BRIEF_2026-10-30_DOC201.md`](reference/doc201_runs/GPT_DISCUSSION_BRIEF_2026-10-30_DOC201.md)

## 一句话

冻结命名统一电磁候选 `doc201_unified_legacy_flux_norm_off`（A1–U6 PASS）；A4 分列 PARTIAL；B1 分层 hold 改善并冻结 thaw→ForcePrior on；B2 真活声学顶 FAIL 且 dens 局部修停码。**本轮关闭；生产 DivΣA 默认未静默回切。**

## §13 六问（关闭答）

| 问题 | 答 |
| --- | --- |
| 相同生产入口能否处理 A、iA、A0+∇χ？ | **命名统一候选能**；生产混合入口 **不能** 宣称统一（A1 INVALID） |
| 是否仍按 case 切换 RHS / 拆梯度？ | 统一候选：**否**。生产默认仍 imag DivΣA ≠ real LegacyFlux |
| 场、局部连续、边界、能量是否共同相容？ | **否（统一候选）**：场阶梯 PASS；U1 ε_div≈14% FAIL；P1 field_rel FAIL |
| 分层残余最先从哪来？ | **声学推进首漂**；ForcePrior 双重计数放大 hold 残量 |
| 顶部是否真跑压力/密度/位置/声学壁？ | **是**（B2 live）；Doc200 冻结顶不算本项 PASS |
| 对法国弱对流 / 热耦合限制？ | EM 相位可用但 U1 连续与功率场分需解释；活顶切向未收口 → 长耦合仍缓 |

## 门禁摘要

| 线 | 状态 |
| --- | --- |
| A1–U6 统一候选 | **PASS** |
| A4 | **PARTIAL**（ε_div↔边权衡 + DivΣA/dual 已入账） |
| B1 hold | **改善**（5 s U_rms ~1.2e-5→9.8e-7） |
| B1 thaw 政策 | **默认 ForcePrior on**（thaw 关 owns） |
| B2 活声学顶 | **FAIL + dens 停码**（壁声学定位） |
| F1 法国固定 T 统一 Q | **未跑**（条件项） |

## 产物入口

- 进度 / 台账 / 合同：`reference/doc201_runs/{DOC201_PROGRESS,VALIDATION_LEDGER,UNIFIED_EM_CONTRACT}.*`
- 命令：`reference/doc201_runs/COMMANDS.md`、`scripts/`
- 下一轮选项：`reference/doc201_runs/NEXT_STAGE_DECISION.md`（**Doc201 已关闭**）

## 明确不做（关闭后仍成立）

- 不以 DivΣA 对照冒充统一候选 A4 PASS  
- 不以统一候选静默替换生产默认  
- 不把 0.2 s thaw U 分叉当自由对流校准  
- 不再叠 B2 dens 局部修  
