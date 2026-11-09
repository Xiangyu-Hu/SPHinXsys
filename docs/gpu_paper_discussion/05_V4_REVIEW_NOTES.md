# V4 计划阅读笔记与本地对齐（讨论稿）

日期：2026-11-03  
来源：`GPU_PAPER_DEVICE_TASKS_V4_2026-11-03.zip`  
状态：已通读；**尚未**把本机执行计划整份替换为 V4（待讨论确认）。

---

## 1. V4 相对 V2/本地粗采样的实质升级

V4 主要解决「有原则、无落地」：把频率写成**可执行命名规则**，并配 case 卡、配置矩阵、绑定脚本与合成样例。

| 能力 | V4 做法 |
|---|---|
| 曲线采样 | `uniform100` / `cycles100` / `event100` / `window51` / `industrial51` / `iterations100` |
| 场输出 | `transient5` / `cycle8` / `event8` / `evolution8` / `iterations8` / `industrial8`；**仅代表轨迹**，默认 5–8 组，≤10 |
| 绑定 | `tools/freeze_output_schedule.py`：先导冻结 T/S/P/[a,b]/W/I 后生成绝对时刻 JSON |
| 分层 | A：每配置 1 条验证；B/C：每配置×3 性能（关 VTP）+ 每例最多 1 条 showcase |
| 记录面 | 保留 V3 的约 282 字段 / 20 类指标（M01–M20）；新增 transfer/sync/profile/energy/missing 等 |
| 主端 | **4090 + 本机 64 核**（与本地共识一致）；外卡用已有 GPU 做 COMMON，不扩跨厂商/多 GPU |

合成校验（VALIDATION_REPORT）：均匀 100、四周期 100、事件 121、迭代 100、benchmark 51、工业流动+热等公式通过；**不证明** writer 已接到求解器。

---

## 2. 逐例：算什么 / 曲线 / VTP（摘自 11 + case_cards）

| ID | 曲线规则 | VTP（仅代表） | 备注 |
|---|---|---|---|
| A1 | uniform100 | transient5 @ GPU FP32 L2 | 稳态剖面另一次；增删与质量收支 |
| A2 | cycles100（约 4 周期×~25 点） | cycle8 | 先导密探针对照后可加密 |
| A3 | uniform100（或迭代轴） | transient5 | D4/D7 二选一 |
| A4 | event100（底 100 + 事件窗 +21） | event8 | **禁止**用稀疏点插值造峰值 |
| A5 | cycles100 | cycle8 | V4 仍在主矩阵；本地曾后置 |
| B1 | iterations100 | iterations8 @ L3 showcase | 固定迭代 vs 同质量耗时分开 |
| B2 | window51（benchmark） | evolution8 @ showcase | 性能关 VTP |
| B3 | window51 | cycle8 @ showcase | 同上 |
| C1–C3 | industrial51 / showcase 周期 100 | industrial8 | V4 仍主线；本地工业另分支 |

性能侧统一：计时**内存累计、结束写一次**；telemetry **1 s**；高成本邻居仅起/中/末。

---

## 3. 与本地已定约定的关系

| 本地约定（03/01） | V4 | 建议 |
|---|---|---|
| 只用 GPU 1 | 未写 | **继续强制** `CUDA_VISIBLE_DEVICES=1` |
| 助手不自跑，用户执行 + tmux | V4 假定 Cursor 可跑队列 | **继续本地执行权**；命令按 V4 口径写 |
| case 进 `extra_source_and_tests/.../gpu_paper_*` | 仍列原始 tests 路径作候选 | **复制后接 V4 writer/采样**，不改原始当跑点 |
| 结果进 `docs/gpu_paper_discussion/results/` | 通用 run 目录 | 本地归档路径保留，文件集按 V4 templates |
| C 组另分支 | V4 含完整 C | **本分支不跑 C**；契约/模板可预留 |
| A5 主线后置 | V4 含 A5 | **仍后置**；矩阵保留但不插队 |
| B1 几何后补 | 同 V4 | 不变 |
| 本地 `04_SAMPLING_POLICY` | 粗默认 | **频率以 V4 的 11/12/case 卡为准** |

冲突优先级（V4 自述）：**11 / 12 / case_cards > 旧 06/07 叙述**；本地运维以 **03_RUN_CONVENTIONS** 为准。

---

## 4. 写论文够不够 + 开跑策略（2026-11-03 更正）

V4 采样与字段面**够写论文**。开跑策略已与用户对齐（详见 `06_PRE_RUN_INSTRUMENTATION_GATE.md`）：

1. **正式跑第一个算例之前**，记录体系（含约 282 字段对应的表结构/writer/状态机）与物理 observer 尽量全部就绪，目标是**尽量不因漏记重跑**。  
2. 依赖工具/权限而不一定采得到的量：开跑前统一「尝试方式 + missing 口径」；采不到填 missing，**不**事后改规范。  
3. 采集频率与 VTP：**以 V4 为准**。  
4. **圆柱 + 柔性尾板 FSI2 后置**；工业三例另分支。  
5. ~~「先跑 case 再分期补字段」~~ **作废**——与不重跑目标冲突。

---

## 5. 已对齐 / 仍可确认

- [x] 频率与 VTP：V4  
- [x] 开跑前采集骨架 + 物理统计一次就绪；工具类统一 missing  
- [x] FSI2 后置；工业另分支  
- [ ] 外卡清单（3090 / 2080S / H20 等）现场确认  
- [ ] 下一步：先做公共 instrumentation + preflight，再接线第一个算例（二维 Poiseuille），命令仍由用户执行  
