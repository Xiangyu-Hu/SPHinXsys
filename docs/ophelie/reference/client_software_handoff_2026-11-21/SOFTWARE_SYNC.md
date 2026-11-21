# 给软件机的最短同步清单（2026-11-21）

物理机已完成甲方 **自然对流 + 机械搅拌** 首轮短测门控。本文件是软件机接入时的唯一入口；细节见同目录 `HANDOFF.md`、`COMPILE_AND_RUN.sh`、`docs/ophelie/230_*.md`。

## 1. 同步什么

| 优先级 | 路径 | 用途 |
|--------|------|------|
| **必读** | `docs/ophelie/reference/client_software_handoff_2026-11-21/` | 本包：effective、summary、本清单 |
| **必读** | `effective_client_nc_sigma16.json` | NC 开发预设（冻结 μ + σ16 借用） |
| **必读** | `effective_client_stir.json` | 搅拌：`mechanics_mode=elastic_driven`，10 rpm，30 kW **电源输入**，电流 unknown |
| 证据 | `nc/*.json`、`stirring/*_summary.json` | 各门控 `passed=1` / `*_PASS=1` |
| 几何 | `tests/.../test_3d_ophelie_client_nc_glass_relax/data/` | NC CAD |
| 几何 | `tests/.../test_3d_ophelie_client_stirring_glass_relax/data/` | 搅拌 CAD + `client_stir_paddle_solver_m.stl` |
| 任务书 | `docs/ophelie/software/OPHELIE_DESKTOP_IMPLEMENTATION_CURSOR_PLAN_2026-11-19.md` | M1/M2 |
| 补充 | `docs/ophelie/228_*.md`、`229_*.md`、`230_*.md` | 弹性桨 / 材料 / 物理落地 |

**注意：** 软件机请拉取本分支最新提交（含甲方 case 源码、几何 STL、交接包）。不要停在旧基线 `2b35c44`。Reload/VTP/bin 不入库，需本机 relax 生成或使用已落盘 summary。

## 2. 门控状态（可接软件，勿夸大）

| 案 | 已过 | 软件可做 |
|----|------|----------|
| **NC** | 等温静水 + 受控 ΔT=20 冒烟 | M1：关搅拌项目、装配、粒子、静水/温差短跑、结果回放 |
| **搅拌** | E0/E1/E2短窗/F0/F1 | M2：`elastic_driven` 预设、桨轴/转速、等温短跑 + 规定加热回放 |

| 标签 | 值 | 含义 |
|------|-----|------|
| NC `passed` | 1 | 静水 / ΔT 功能短跑 |
| `PADDLE_HOLD_PASS` | 1 | E0 |
| `LOAD_RESPONSE_PASS` / `ELASTIC_PADDLE_ACTIVE` | 1 | E1，桨未锁刚体 |
| `DRIVEN_ELASTIC_ROTATION_SHORT_PASS` | 1 | E2 **短窗**（θ≈0.026 rad） |
| `ELASTIC_FSI_SHORT_PASS` | 1 | F0 等温弹性 FSI |
| `HEATED_ELASTIC_INTEGRATION_PASS` | 1 | F1 规定体积热（**非 EM**） |
| `full_revolution_pass` | **0** | 整圈未做 |
| 真实 EM Q | **未做** | 电流 unknown |

## 3. 软件机禁止事项

1. **禁止**把 E2 短窗标成「整圈 / 2π PASS」。
2. **禁止**把 F1 规定体积热标成「已计算 EM」或「30 kW 吸收功率已验证」。
3. **禁止**把测点碰到 ~1000℃ / ~1400℃ 标成合同定量 PASS（窗口未齐）。
4. **禁止**把刚体动画冒充 `elastic_driven`；界面应反映弹性桨已接通。
5. 30 kW 是**发电机输入**，不是法国 50/60 kW 吸收功率，勿照抄法国预设。

## 4. 求解器入口（CLI 复跑）

| 模式 | 目标 | Reload |
|------|------|--------|
| NC relax | `test_3d_ophelie_client_nc_glass_relax` | 自产 `bin/reload` |
| NC 推进 | `test_3d_ophelie_french_natural_convection_frozen_q` + `--client-nc-geometry --client-nc-material` | NC reload |
| 搅拌 relax | `test_3d_ophelie_client_stirring_glass_relax --bodies=all` | 自产 `bin/reload` |
| 搅拌 F0 | `test_3d_ophelie_client_stirring_f0_elastic_fsi` | 搅拌 reload |
| 搅拌 F1 | `test_3d_ophelie_client_stirring_f1_heated_elastic` | 同上；`--heat-power=` 为开发 Q |

一键命令：同目录 `COMPILE_AND_RUN.sh`。F0/F1 不锈钢亚步很慢，**软件演示优先用已落盘 summary/VTP**，不必每次重跑满窗。

## 5. 软件机开工顺序（建议）

1. 同步本包 + 甲方 case 源码/几何。  
2. 接 NC 示例项目（搅拌开关关）→ 读 `effective_client_nc_sigma16.json`。  
3. 接搅拌示例项目（开关开）→ 读 `effective_client_stir.json`，`mechanics_mode=elastic_driven`。  
4. 桌面：配置 → 预处理 → 启动求解器 → 回放已有短跑结果。  
5. 发行包验收时再按任务书 §2.1 / §7 复测；物理缺口（整圈、EM、定量 T）单列，不挡 M1/M2 首轮。

## 6. 仍留给物理机（不挡软件首轮）

- E2 整圈（∫ω dt = 2π）  
- 电流到位后真实 EM Q / σ 对照  
- 1000℃ 测点时间序列定量  
- 夹持 device 化 / 阻尼加速（可选）
