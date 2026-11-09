# 开跑前门禁：采集体系一次就绪（避免重跑）

日期：2026-11-03  
权威采样/字段：V4（`11` / `12` / case 卡 / templates / parameter_dictionary）  
本地运维：`03_RUN_CONVENTIONS.md`（GPU1、用户执行、gpu_paper 副本、结果归档、tmux）

## 原则（已对齐）

1. **正式论文跑点开始前**，尽量把 V4 要求的记录面（约 282 字段所属的 CSV/JSON/manifest 体系）全部接好、口径统一。  
2. 第一个正式算例的意义是产出**可直接进论文链路的完整数据**，不是「先跑一版再补插桩」。  
3. 物理探针 / 收支 / 步数等 **observer 必须事先安插**；漏记会导致整条轨迹作废重跑——这是最高优先级。  
4. 依赖外部工具、权限或硬件而不一定采得到的量（如部分 Nsight 计数器、整机能耗）：**开跑前统一约定**  
   - 尝试路径（命令、工具、是否在代表短窗）  
   - 采不到时的 `measurement_status` / `missing_reason` 写法  
   - **不得**事后才发明缺失口径，也不得用 `0` 填未知。  
5. **圆柱 + 柔性尾板 FSI2（原 A5）后置**；工业三例另分支。其余算例共用同一套采集骨架。

「一次就绪」≠ 第一个算例上就要测出全部硬件计数器的真实数；  
「一次就绪」= **writer、表头、状态机、缺失规则、物理 observer、preflight 在开跑前全部就位**，正式跑只填数或填 missing，不再改插桩结构。

---

## 开跑前必须完成的工作包

### P0 — 公共记录骨架（所有算例共用）

- [ ] 按 V4 `templates/` 实现/对齐全部正式输出文件与列名（至少：manifest、run_summary、stage_timing、step_stats、body_stats、diagnostics、probes、telemetry、throughput、transfer_timing、sync_timing、memory_stats、missing_metrics、run_aggregate 等适用项）。  
- [ ] `measurement_status` 枚举统一：`measured` / `derived` / `estimated` / `unsupported` / `not_applicable` / `permission_denied` / `not_collected`。  
- [ ] `missing_metrics.csv`：字段、状态、原因、工具错误证据、补救动作；**禁止用 0 冒充未测**。  
- [ ] 计时账本：`wall_solver`、init/prewindow/output 边界、阶段互斥、传输/等待不重复相加（V4 `06`）。  
- [ ] preflight（V4 `13`）：缺必需字段或表头不通过 → **禁止进入正式长跑**。  
- [ ] runner：validate / benchmark / showcase / profile 模式分离；benchmark **关闭 VTP**；showcase 单独 run_id。  
- [ ] 本机：`CUDA_VISIBLE_DEVICES=1`、tmux `gp_*`、结果归档到 `docs/gpu_paper_discussion/results/`。

### P0 — 物理 observer（按算例，开跑该例前必须绿）

每个纳入本分支的算例，在**该例第一次正式跑之前**完成（可先做公共库，再逐例接线）：

| 算例（全称） | 必须事先接好的物理统计（摘要） | V4 曲线 / 代表场 |
|---|---|---|
| 二维 Poiseuille 通道流 | 中心速度、流量、压降、质量收支、N(t)、终态剖面 | uniform100 / L2 GPU：transient5 |
| 二维振荡梁 | 尾端位移速度、能量、步数/dt | cycles100 / cycle8 |
| 窗框传热（D4 或 D7） | 测点温度、内外热流、残差/热收支 | uniform100 / transient5 |
| 二维两相溃坝 | 压力探针、水位/前沿、峰值与到达时刻、相质量、机械能 | event100 / event8 |
| 三维复杂几何粒子生成/松弛 | 残差曲线、质量目标事件、起中末质量统计 | iterations100 / iterations8 |
| 三维溃坝性能 | 窗内最低诊断 51 点；showcase 完整曲线+场 | window51 / evolution8 |
| 三维扭转柱性能 | 端部运动、能量、关系重建次数等 | window51 / cycle8 |
| 圆柱+柔性尾板 FSI2 | （后置）尾端轨迹、载荷、子步保护等 | cycles100 / cycle8 |
| 工业齿轮箱 / 搅拌 / 储箱热晃荡 | （另分支）按 V4 industrial 规则 | industrial51 / industrial8 |

对每一例还要：真实测点坐标与单位、writer 字段名映射、`freeze_output_schedule` 绑定后的绝对时刻 JSON、`requested_time`/`actual_time`。

### P1 — 「可能采不到」但开跑前必须统一的工具类指标

在**任何正式算例之前**写好《工具采集尝试表》（可放 `results/_registry/tooling_attempt_policy.csv`）：

| 类别 | 计划尝试 | 失败时统一填写 |
|---|---|---|
| GPU kernel profile（occupancy、寄存器、cache、带宽等） | 仅 B2/B3 代表短窗；记录命令与原始报告路径 | `unsupported` / `permission_denied` + 错误日志摘要 |
| 传输 / 同步细分 | 能从现有 wait/队列拿则 measured；否则短窗 attempt | missing_reason，不阻塞主 wall_solver |
| 功率 / 温度 / 频率 / 利用率 | 外部 1 s 采样（本机 GPU1） | 覆盖不足则说明，不做伪峰值结论 |
| GPU/CPU 包能量积分 | 仅当功率采样覆盖有效时 derived | 否则 `not_collected` / `not_applicable` |
| Roofline / 消融 | 有开关与可信计数才做 | 无则 `not_applicable`，不新造求解器 |

原则：**结构先有、数值可空、理由必填**；正式跑中途不得改列名或新增必填物理列。

### P2 — 短测验收（仍不算正式论文点）

开跑前门禁用**短窗/粗档**验证（用户执行命令）：

- [ ] 每种模式真实写出完整表头 + ≥2 个物理采样点（空绑定不得放行）。  
- [ ] VTP 开关按模式正确；代表轨迹预计帧数与存储估计已打印。  
- [ ] 计时开关 on/off 配对短窗（B2/B3 代表）记录插桩开销比。  
- [ ] 工具类指标按上表尝试一次，生成首份 `missing_metrics` 样例（可来自短测）。  
- [ ] 通过后方可提交该算例的**第一次正式** validate/benchmark。

---

## 本分支算例顺序（全称）

1. 公共骨架 + 工具缺失政策（P0/P1）就绪  
2. 二维 Poiseuille → 二维振荡梁 → 窗框传热 → 二维两相溃坝  
3. 三维溃坝性能、三维扭转柱性能（含 group / CPU 线程）  
4. 三维复杂几何粒子生成/松弛（等你提供几何后）  
5. **后置**：圆柱 + 柔性尾板 FSI2  
6. **另分支**：工业齿轮箱、工业搅拌、两相储箱换热晃荡  

---

## 与此前「分期」说法的更正

此前「先接主性能再补 282 字段」**不符合**「尽量不重跑」目标，废止。  
正确表述：

- **开跑前**：记录体系、物理 observer、缺失口径、preflight **全部就绪**；  
- **开跑后**：只跑点、归档；缺工具的填 missing，**不**再改物理统计列定义；  
- 若开跑后发现物理量漏记 → 该配置必须重跑，并升 `experiment_revision`。
