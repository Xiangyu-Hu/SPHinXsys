# 第5阶段口径方案 — 法国搅拌矩阵（Fig8a 启动 / EM 纯算 / Q←Fig10 校准 / 60s·5s）

**日期**：2026-11-26  
**分支**：`feature/EM-Biot-Savart`  
**读者**：下一个 Cursor / 4090 操作者（无需重聊历史）  
**身份**：工程展示与综合对比脚手架；**不宣称论文 PASS / FRENCH_VALIDATED**

---

## 0. 给下一个 Cursor 的一句话

用 **法国搅拌** `test_3d_ophelie_french_stirring_em` 做可复现矩阵：每个甲方关心的变量类跑 **3 个水平**；一律 **STIR Fig.8a 数字化温度 IC**；**电磁场纯算不校准**；**熔体焦耳热 \(Q\) 在每次 EM→Q 更新后按 STIR Fig.10 形状 blend 校准**；每 case **end-time=60 s**，**state-record-interval=5 s**，中间帧也要能对比。

甲方合同表述（摘录）：

> 第5阶段（14个月内）：在充分模拟分析的基础上，完成不同**搅拌工况、电磁参数、结构参数**下的综合性能评估与设计对比，提交第三次计算分析报告。

本方案用**法国文献搅拌栈**先把“三变量×三水平”方法论跑通；甲方几何上的最终报告另用 Doc248 **C2r 刚体桨**承接（见文末）。二者**Reload / 图号 / 功率不可混**。

---

## 1. 已锁定的物理/展示口径（用户 2026-11-26 确认）

| 项 | 锁定 |
|----|------|
| 电磁场（电流标定、edge-flux/φ、σ(T) 周期刷新等） | **纯算，不按文献图校准** |
| 给流体用的焦耳热 \(Q\) | 每次 EM→Q 后与 **STIR Fig.10** 形状 **blend**（同 NC 的 Fig.6-correct 逻辑），总功率归一到该 case 的 \(P_\mathrm{abs}\) |
| 为何可以每步/每次刷新校准 | Fig.10（及 NC Fig.6）是文献**稳态/准稳态** \(Q\)；搅拌流动动态，但实验室系热源可近似准稳态。校准是为冲准稳态热闭合，**不是**声称 EM 预测了 Fig.10 |
| 给甲方看的焦耳热图 | 默认展示 **校准后 \(Q\)**；图注必须写 **ASSUMED / shape-corrected** |
| 给甲方看的电磁能力 | 另附 **未校准 EM** 功率/残差；禁止让人以为“电磁也按图调过” |
| 温度 IC | **仅** `SEED_STIR_Fig8a_axisym_ASSUMED_T_rz.csv`（搅拌 Fig.8a；禁止 NC Fig.8、禁止冷壳 750 K 映射） |
| 换热 BC | Level-1 **skull + free-crust**（与生产搅拌一致）；**不要** `--balance-heat-loss` 与 skull 同开 |
| 耦合 | `--thermo-em-coupling=periodic`（推荐 `--em-update-interval=0.5`） |
| Fig.11 等温旁路 | `--zero-q` / 关热 **不是**本矩阵；勿与本方案混跑同一 OUT |

NC 线（无桨）继续用：**Fig.8 T + EM 纯算 + \(Q\)←Fig.6 correct**（已有 `run_fr_nc_em_fig6_correct.sh`）。本文件主写**搅拌矩阵**。

---

## 2. 变量矩阵（每个变量 3 个 case）

基线（文献局部 Fig.8–10）：**10 rpm · 60 kW · skull-delta=0.005 m**。  
三类变量各扫 3 水平；**改一类时其余固定在基线** → **9 个正式 case**（若算力紧，可先跑每类的“低/基线/高”共 9 个，不要一上来全因子 27）。

### 2.1 搅拌工况 `SW`（rpm）

| case_id | rpm | \(P_\mathrm{abs}\) | skull-delta | 说明 |
|---------|-----|-------------------|-------------|------|
| `SW_rpm05` | **5** | 60 kW | 0.005 | 低搅拌 |
| `SW_rpm10` | **10** | 60 kW | 0.005 | 基线 |
| `SW_rpm20` | **20** | 60 kW | 0.005 | 高搅拌 |

### 2.2 电磁参数 `EM`（吸收功率）

文献附近水平：56 / 60 / 70 kW（见 `STIR_verified_parameters.csv`）。

| case_id | rpm | \(P_\mathrm{abs}\) | skull-delta |
|---------|-----|-------------------|-------------|
| `EM_P56` | 10 | **56 kW** | 0.005 |
| `EM_P60` | 10 | **60 kW** | 0.005 |
| `EM_P70` | 10 | **70 kW** | 0.005 |

CLI：搅拌侧用 `--target-absorbed-power=`（或文献包等价项）；**电流/场仍由求解器算**，只改功率目标。

### 2.3 结构参数 `ST`（本轮可跑的软结构代理）

真改桨浸深/桨径/坩埚比例需要 **重新 relax + 新 Reload**，本轮矩阵先用 **skull / free-crust 厚度** 作“结构隔热层”代理（仍写 ASSUMED）。真几何扫另开子阶段。

| case_id | rpm | \(P_\mathrm{abs}\) | skull-delta = free-crust-delta |
|---------|-----|-------------------|-------------------------------|
| `ST_d3` | 10 | 60 kW | **0.003 m** |
| `ST_d5` | 10 | 60 kW | **0.005 m** |
| `ST_d8` | 10 | 60 kW | **0.008 m** |

---

## 3. 单 case 运行约定（全部统一）

```text
二进制: test_3d_ophelie_french_stirring_em
Reload:  .../test_3d_ophelie_french_stirring_glass_relax/bin/reload   # 必须搅拌 Reload，禁止 NC Reload
IC_T:    docs/ophelie/reference/standard_data/doc247_seeds/SEED_STIR_Fig8a_axisym_ASSUMED_T_rz.csv
IC_Qshape: SEED_STIR_Fig10_Qshape_ASSUMED_rz.csv   # 若尚无 → 见 §5 缺口
end-time:              60
state-record-interval: 5          # → 约 t=0,5,10,...,60 共 ~13 帧，供中间态对比
em-update-interval:    0.5
Q-correct:             on, blend=0.8, every=1（每次 EM 更新后）
csv/screen:            建议 csv-every=10 screen-every=20
VTP:                   开（要中间态）；磁盘紧则只留 french_T_rz_t*.csv + monitor
max-wall-hours:        建议 ≥2（单 case）；整矩阵另排班
```

**输出目录建议**（gitignore 大文件；脚本与本 md 进 git）：

```text
docs/ophelie/reference/doc247_longrun_four_cases/overnight_6h/phase5_fr_stir_matrix_60s/
  SW_rpm05/  SW_rpm10/  SW_rpm20/
  EM_P56/    EM_P60/    EM_P70/
  ST_d3/     ST_d5/     ST_d8/
  MATRIX.md                 # 自动生成的退出码/墙钟/末态 T/Np 表
  figures_compare/          # 各 case 末态与中间帧 vs Fig8a / Fig10
```

每 case 至少保留：

- `french_stirring_em_monitor.csv`
- `french_q_rz.csv` / `french_T_rz.csv` + `trz_snaps/french_T_rz_t*.csv`（5 s 间隔）
- `french_thermo_em_coupling.csv`（EM 刷新记录）
- `exit_code.txt`、`run.log`（log 可本地、可 ignore）
- 可选 VTP：`out/GlassBody_ite_*.vtp`

---

## 4. 对比与验收（软门槛，非 PASS）

对每个 case、每个记录时刻 \(t=0,5,\ldots,60\)（至少末态）：

1. **T**：对折后的 STIR Fig.8a 数字化场 vs `french_T_rz_t*`（corr / RMSE；标 ASSUMED）  
2. **Q**：校准后 \(Q\) vs Fig.10 形状（若有 keeper）；同时存一份 raw EM \(P\) / 残差  
3. **搅拌**：\(N_p\)、\(U_\mathrm{rms}\)、转数；Fig.11 线仅趋势，正式 Fig.11 仍用等温 `--zero-q` 旁路  
4. **矩阵对比**：同 \(t\) 下三水平并排（SW / EM / ST 各一张图）

禁止：为贴线改 \(\mu\) 使 \(N_p\cdot Re=47\)；禁止把 NC Fig.6/8 种子接到搅拌 Reload。

---

## 5. 实现缺口（状态 2026-11-26）

| 项 | 状态 | 路径 / CLI |
|----|------|------------|
| Fig.10 Q shape seed | **已放 PLACEHOLDER**（`PLACEHOLDER_NOT_FIG10`；NC Fig.6→unit box） | `doc247_seeds/SEED_STIR_Fig10_Qshape_ASSUMED_rz.csv` + `.meta.json` |
| 搅拌 Fig.10-correct | **已实现**（EM 粒子 Q 上 blend+rescale，再 deposit 欧拉网格） | `--fig10-q-correct` `--fig10-q-correct-blend=` `--fig10-q-correct-every=` `--init-Q-rz-csv=`；共享头 `electromagnetic_ophelie_doc247_q_shape_correct.h` |
| 矩阵脚本 | **已有** | `02_fr_stir/run_fr_stir_phase5_matrix_60s.sh` |
| 后处理 | **已有** | `02_fr_stir/post_phase5_matrix_compare.py` → `MATRIX.md` + `figures_compare/` |

**仍欠（正式报告前）**：用文献扫描数字化 **真** STIR Fig.10 场，替换 PLACEHOLDER seed，并去掉 meta 里的 `PLACEHOLDER_NOT_FIG10`。在此之前 MATRIX / 第三次报告法国章必须标明 **PLACEHOLDER**。

---

## 6. 4090 建议跑法

```bash
cd /path/to/SPHinXsys
git checkout feature/EM-Biot-Savart && git pull
# 编译搅拌 +（若改了 correct）相关目标
cmake --build "$BUILD" --target test_3d_ophelie_french_stirring_em -j"$(nproc)"

# 缺口补齐后：
bash docs/ophelie/reference/doc247_longrun_four_cases/02_fr_stir/run_fr_stir_phase5_matrix_60s.sh
# 或 tmux 里分段：先 SW_*，再 EM_*，再 ST_*
```

单卡建议 **顺序跑 9 case**（避免抢 GPU）。墙钟粗估：此前 ~20 s 搅拌约数分钟级；60 s ×9 ≈ 可按“数小时级”排班（以本机 monitor 墙钟外推）。

并行策略（双卡/多机）：按 `SW_*` / `EM_*` / `ST_*` 三组拆三个会话，**同一 Reload 只读**。

---

## 7. 与现有产物的关系（防混）

| 路径/脚本 | 用途 | 是否本矩阵 |
|-----------|------|------------|
| `run_fr_stir.sh` | 单点生产搅拌（Fig8a + periodic + skull） | 基线模板；默认 end 已改短窗 |
| `run_fr_stir_fig11_*` / `zero-q` sweep | Fig.11 等温 / 软凑 Np | **否** |
| `run_fr_nc_em_fig6_correct.sh` | NC：Fig8 T + EM + Fig6 Q-correct | NC 线；图号不同 |
| `02_FR_STIR_Fig8a_EM` 10 s | 短窗展示 | 历史；矩阵用新 OUT 根目录 |
| Doc248 `04_run_C2r_rigid_stir.sh` | 甲方刚体桨冲 ~1400 K | 甲方线；**无 Fig8a** |

---

## 8. 甲方第5阶段如何映射（报告叙事）

1. **方法章（法国）**：本矩阵 9 case — 证明“搅拌/电磁/结构（代理）”三轴可扫，且 T 从 Fig.8a、\(Q\) 向 Fig.10 软闭合、电磁纯算。  
2. **应用章（甲方）**：同一套“三变量×三水平”换 Doc248 几何 + C2r；功率/转速用甲方数；**不要**写法文 Fig.8a 已在甲方坩埚上复现。  
3. **第三次计算分析报告**：法国矩阵图 + 甲方矩阵图分栏；每张图标注 ASSUMED / 非 PASS。

---

## 9. 下一个 Cursor / 4090 的任务清单（按序）

1. 读本文件 + `DECISIONS_LOCKED_20261125.md` + `french_stirring_em/README.md`。  
2. ~~补 PLACEHOLDER Fig.10 seed + Fig.10-correct + 矩阵/后处理脚本~~（已合入工作树；需 **compile** 搅拌目标）。  
3. 4090：**先 compile**，再冒烟 `CASES=SW_rpm10 END_PT=10`，再满矩阵 9×60 s。  
4. 后处理脚本在矩阵脚本末尾自动跑；检查 `MATRIX.md` 与 `PLACEHOLDER_NOT_FIG10` 标注。  
5. （可选高优先）数字化真 Fig.10 替换 PLACEHOLDER。  
6. 入库：脚本/seed/头文件/本 md；**不要** commit `overnight_6h` / `phase5_*` 大 CSV/VTP。

---

## 10. 关键路径速查

| 用途 | 路径 |
|------|------|
| 本方案 | `docs/ophelie/reference/doc247_longrun_four_cases/PLAN_PHASE5_FR_STIR_MATRIX_60s_20261126.md` |
| 4090 旧交接 | `.../WORK_SUMMARY_4090_HANDOFF_20261126.md` |
| 口径锁定 | `.../DECISIONS_LOCKED_20261125.md` |
| Fig8a seed | `docs/ophelie/reference/standard_data/doc247_seeds/SEED_STIR_Fig8a_axisym_ASSUMED_T_rz.csv` |
| 搅拌二进制 | `tests/.../test_3d_ophelie_french_stirring_em/` |
| 搅拌 Reload | `tests/.../test_3d_ophelie_french_stirring_glass_relax/bin/reload/` |
| NC Fig6-correct 参考实现 | `tests/.../test_3d_ophelie_french_natural_convection_frozen_q/` + `run_fr_nc_em_fig6_correct.sh` |
| 甲方刚体桨 | `docs/ophelie/reference/doc248_client_cases/04_run_C2r_rigid_stir.sh` |
