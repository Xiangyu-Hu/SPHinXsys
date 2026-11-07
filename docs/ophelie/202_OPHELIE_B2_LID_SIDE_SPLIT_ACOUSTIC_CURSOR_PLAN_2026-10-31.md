# Doc202：B2 盖/侧分体声学（非 dens）

日期：2026-10-31  
前置：Doc201 关闭（B2 dens 局部修停码；`ablate=wall` 定位壁声学）  
证据树：`docs/ophelie/reference/doc202_runs/`（单层；禁止嵌套 PACKAGES）

## 1. 目标

在**不再叠 dens 开关**的前提下，用可测假说分离「盖」与「侧/底」声学接触对活顶剪切的责任。

## 2. 已否决（不得重做）

| 项 | 证据 |
| --- | --- |
| 盖 free-slip dens 鬼元 + 一步力∥n | L2≈1.05 |
| H-dens-fs | L2≈1.05 |
| 全壁 free-slip dens/力∥n | L2≈0.99 |
| `ablate=acoustic` 作壁对照 | disp=0（无效） |

## 3. 本轮唯一主假说 H1（反臂）

**侧壁+底：整对跳过声学壁接触；盖：保留现有 free-slip 法向不可穿透（二步法向反射 + DissipativePJump·n；一步力∥n）。**

可测：

| 结果 | 判读 |
| --- | --- |
| core L2 → ~0.12（近 `ablate=wall`）且 uz/U0≤1%、disp>0 | 杀剪切的是**侧/底 Riemann**，不是盖 |
| L2 仍 ~1 | 盖接触仍是主因 → 下一刀 H1b（盖整对 continue，侧保留） |
| uz/U0 崩 | 盖法向不可穿透不够 → 记录失败，不叠 dens |

## 4. 明确不做

- dens / H-dens-fs / AcousticDensMode 再扫  
- 静默改生产 DivΣA / 法国长耦合  
- 冻结位置 / 每步 uz=0 投影  
- 无假说连改第二刀（两次无改善才停码）

## 5. 合同（沿用 Doc201 B2）

`french_natural_reduced_pressure --doc201-b2-live-top-shear`，`implicit_split`，`end=0.1 s`，`U0=0.1`，`H=0.1`，`R=0.2`，`dp=0.01`，`c0=20`，`μ=10`。  
门禁：core `L2_rel≤5%`，`amp_rel≤5%`，`uz_rms/U0≤1%`，`disp_rms>0`。  
新 CLI：`--doc202-b2-ablate-side-acoustic`。

## 6. 交付

1. `DOC202_PROGRESS.md`、`VALIDATION_LEDGER.csv`、`METRIC_CONTRACT.md`  
2. `H1_SIDE_ABLATE_REVIEW.md` + 运行目录  
3. `NEXT_STAGE_DECISION.md`  
4. 源码：壁接触 `AcousticWallContactOmit` + 法国标记/CLI
