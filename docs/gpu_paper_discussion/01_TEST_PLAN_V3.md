# GPU 论文测试计划 V3（本分支）

日期：2026-11-03  
基线：`GPU_PAPER_DEVICE_TASKS_V2`（2026-11-02）  
修订原因：主执行机改为本机 4090+64 核；工业 case 另分支；A5 后置；B1 几何后补。

比较用途标记沿用 V2：**V** 物理验证，**P** 精度，**R** 分辨率，**G** 工作组，**H** 硬件，**B** 阶段瓶颈，**C** 条件一致性。计时、manifest、吞吐率与加速比口径仍以 V2 的 `04_DATA_CONTRACT.md` 为准，除非本文件明确改写。

---

## 1. 范围与分工

### 1.1 本分支做 / 不做

| 范围 | 本分支 |
|---|---|
| A1–A4 精度与物理验证 | 做（本机 CPU FP64/FP32 + GPU FP32） |
| A5 FSI2 GPU 移植与验证 | **延后**：其余主线完成后再做 |
| B1–B3 五档规模、group、本机 CPU/GPU 性能 | 做 |
| C1–C3 工业齿轮箱 / 搅拌 / 热晃荡 | **不做**（另开分支；本计划只保留接口说明） |
| 跨其他 GPU（2080S / 3090 / H20）COMMON 配对 | 可选；本机先冻结 COMMON 与交接包 |
| 接触因子 / 覆盖率诊断 | 不作为本分支指标（工业分支再处理） |

### 1.2 设备角色（相对 V2 重映射）

| 设备 | V2 角色 | V3 角色 |
|---|---|---|
| **本机：4090 + 64 核 CPU** | （原分属 3090Ti 与 CPU64） | **主执行机**：公共工具、A 组、B 组、本机 CPU 基线与线程扫描、冻结包 |
| 3090 / 3090 Ti 服务器 | 主先导 | **可选跨卡**；不稳定则跳过或后补 |
| 2080S / H20（及其他卡） | 跨 GPU COMMON | 仍为跨硬件点；收到冻结包后只跑 COMMON×3 |
| 工业专用分支 / 机器 | C 组 | 工业 case 全部改到那里执行 |

本机双 4090：**固定使用 GPU 1**（`CUDA_VISIBLE_DEVICES=1`，UUID `GPU-6d9732cb-5aed-bfad-6d66-736689de9900`）；GPU 0 常被占用，禁止默认用 0。CPU 正式性能基线默认 **64 物理核**；线程扫描单独成表。

**执行权**：助手只改代码并给出完整命令；**由用户执行**。长跑必须用 tmux（会话名 `gp_*`）+ 文件日志，断线后可 `tmux attach` / `tail` / `scripts/check_runs.sh` 查看。详见 `03_RUN_CONVENTIONS.md`。

**源码与归档**：论文 case 复制到 `tests/extra_source_and_tests/{2d,3d}_examples/test_*_gpu_paper_*`，不在原始 `tests`/`tests_sycl` 算例上堆论文插桩；跑完结果归档到 `docs/gpu_paper_discussion/results/<CASE_ID>/`。

### 1.3 论文证据与案例对应（无偏向性，三块并重）

| 论文证据块 | 本分支案例 | 说明 |
|---|---|---|
| 正确性 / 精度 | A1–A4（A5 后补） | 解析解 / 文献 / COMSOL + FP64/FP32/GPU 互差 |
| 规模与实现性能 | B1–B3 | 五档、group、阶段、内存、本机 CPU↔GPU |
| 跨硬件可移植 | B 组 COMMON（+ 可选外卡） | 同输入同窗口同 group 政策 |
| 工业级端到端 | C1–C3（另分支） | 本分支只保证数据契约与交接接口兼容 |

---

## 2. 案例矩阵（本分支）

| ID | 案例 | 分辨率 | 本机任务 | 其他 GPU | 用途 | 优先级 |
|---|---|---|---|---|---|---|
| A1 | 2D mixed Poiseuille | L1/L2/L3 | 每档 CPU FP64、CPU FP32、GPU FP32 各 1 | 不重复 | V/P/R | P0 |
| A2 | 2D oscillating beam | L1/L2/L3 | 同 A1；L2 可选 GPU FP64×1 | 不重复 | V/P/R | P0 |
| A3 | 窗框传热 D4 优先（或成熟 D7） | L1/L2/L3 | 同 A2；D4/D7 只选一个 | 不重复 | V/P/R | P0 |
| A4 | 2D two-phase dambreak | L1/L2/L3 | 同 A1 | 不重复 | V/P/R | P0 |
| A5 | 2D cylinder + flexible tail FSI2 | L1/L2/L3 | GPU 移植 + 同 A1 | 不重复 | V/P/R | **P2 后置** |
| B1 | 3D 复杂几何单分辨率 level-set / 松弛 | L1–L5 | GPU 五档×3；L3 group 32/64/128/256×3；CPU FP32 五档×3 | COMMON×3/卡 | R/G/H/B | P1（几何后补） |
| B2 | 3D dambreak | L1–L5 | 同 B1 | COMMON×3/卡 | R/G/H/B | P1 |
| B3 | 3D twisting column | L1–L5 | 同 B1 | COMMON×3/卡 | R/G/H/B | P1 |
| C1–C3 | 工业三例 | COMMON | **本分支跳过** | 另分支规划 | H/B/C | 外置 |

说明：

- B1 指 **复杂几何的 level-set 粒子生成与单分辨率松弛**（V2 中的 B1 / 旧 P1）。几何体做到该 case 时再传入并冻结 hash；当前仓库可见 `barrel.stl` 等仅作占位/候选，非正式冻结几何。
- COMMON 优先 L3，须能在最小显存目标卡上跑通；本机 4090 显存充裕，不能单独把 COMMON 定得过大以致 2080S 无法配对。
- 工作组主基线暂定 **64**；扫描后若改统一设置，须升级 `experiment_revision` 并重跑受影响主性能点。

---

## 3. 本机工作量估算（上限）

不含短测、失败排查、可选 profile、工业 C 组。

| 类别 | 条数（约） | 备注 |
|---|---|---|
| A1–A4 精度 | 4×3×3 = **36** | CPU64/CPU32/GPU32 |
| A2/A3 可选 GPU FP64 | **0–2** | 仅中档 |
| B 五档 GPU 性能 | 3×5×3 = **45** | FP32，group=统一政策 |
| B L3 group 扫描 | 3×4×3 = **36** | 与 G=64 且条件全同时可复用对应 B 点，不重复计样本 |
| B 五档 CPU FP32 | 3×5×3 = **45** | 正式 64 核基线 |
| CPU 线程扫描 | B2/B3 × COMMON × {1,8,16,32,64} | 先导 1 次，正式点补至 3 次；SMT 可选 |
| A5 | **+约 9**（移植完成后） | 后置，不阻塞主线冻结 |
| **本机主线合计** | **约 160+**（A1–A4 + B） | 计划上限，按预算裁剪最细档 |

外卡每台（若启用）：B1–B3 COMMON×3 ≈ 9 次（本分支无 C 组时）。

---

## 4. 执行阶段（建议顺序）

### Phase 0 — 盘点与隔离（先做）

1. 记录本机双 4090 UUID、驱动、SYCL backend、CPU 拓扑 / NUMA、编译器与依赖。
2. 盘点源码路径与 git 状态；独立实验分支 / worktree，不 reset/clean 用户未提交工作。
3. 输出 `source_inventory.csv`（每例 repo/ref/SHA、CPU 路径、GPU 路径、参考数据、状态）。
4. 确认 A1–A4、B2、B3 可编译短测；A3 参考数据与 D4/D7 选择；B1 标记 `geometry_pending`。

### Phase 1 — 公共采集与构建（T1）

1. 隔离构建：CPU-CK FP32、CPU-CK FP64、GPU FP32；启动打印 `sizeof(Real)`、device、group、版本。
2. GPU FP64 仅作 A2/A3 可选；不通则标 `unsupported`，不阻塞主线。
3. runner / preflight / 窗口计时 / 日志解析 / 汇总脚本；修复已知 timer 闭合问题（如 Poiseuille 分段起点）。
4. dry-run：单 case 单配置可复现目录结构（见 V2 数据契约）。

**放行条件**：至少 1 个流体短测 + 1 个固体短测，CPU/GPU 均写出必要字段。

### Phase 2 — A 组精度（A1→A2→A3→A4；A5 跳过）

每例：粗档短测 → 三档正式轨迹。

| 顺序 | Case | 验证重点 | 阻塞处理 |
|---|---|---|---|
| 1 | A1 Poiseuille | 解析剖面、流量、质量收支 | 必修通，作为计时与精度模板 |
| 2 | A2 梁 | 频率/幅值/能量；可选 GPU64 | 参考不适用时换合适文献 |
| 3 | A3 窗框 | COMSOL；材料界面；稳态热流 | 无参考则标待补，不伪造；D4/D7 二选一 |
| 4 | A4 两相溃坝 | 界面、压力峰、守恒 | 对齐文献几何与无量纲 |
| 后 | A5 FSI2 | GPU 移植 + 标准尾板双向耦合 | **全部主线后再做** |

交付：`A_case_validation.csv`、误差图、参考来源、三档原始数据。

### Phase 3 — B 组规模 / group / 本机 CPU 配对

1. **B2、B3 先**：几何与源码已明确，先定 L3 短窗 → group 扫描 → 五档×3（GPU）→ 同窗 CPU FP32 五档×3。
2. **B1**：用户提供复杂几何后，冻结 STL/缩放/hash，确认体填充语义，再走与 B2/B3 相同流程。
3. 定 **COMMON**（优先 L3）：记录 N、邻居、显存；若外卡目标含 2080S，以其可运行为约束。
4. CPU 线程扫描：在 COMMON 上对 B2、B3 做 1/8/16/32/64；正式报告点补满 3 次。
5. 每例选一代表规模做 profile（必要时加最细档）；主加速比不用同步插桩时间。

交付：`resolution_table.csv`、`group_sensitivity.csv`、`B_scaling.csv`、`cpu_thread_scaling.csv`、阶段与内存曲线。

### Phase 4 — 冻结包（本分支非工业部分）

产出 `experiment_bundle`（A1–A4 + B1–B3）：

- 源码 ref / 补丁 / 依赖说明
- 输入、初态、resolution / COMMON / group / 窗口 / 计时政策
- run 清单、SHA256、少量参考结果
- `README_HARDWARE_HANDOFF`（仅已验证过的真实命令）

此后若改物理或计时边界：升 `experiment_revision`，受影响配对标待重跑。

### Phase 5 — 外卡 COMMON（可选）

2080S / 3090 / H20 等：每例 B 的 COMMON×3；禁止私改 dp/窗口/group。  
OOM 反馈协调端，或冻结额外小规模公共点并在本机补配对。

### Phase 6 — A5 FSI2（后置）

1. 以 CPU `fsi2` 为蓝本做 CK/SYCL 移植；禁止换成浮体/游鱼。
2. 固定流固 dp 比与真实尾板层数；统一初态。
3. 完成 L1–L3 × CPU64/CPU32/GPU32；纳入验证表与（如需要）附录性能点。

### Phase 7 — 工业 C 组（另分支，本文件仅接口）

另分支须遵守同一数据契约：COMMON、manifest、最低诊断、去掉专用 contact-factor 统计、跨设备同版本。  
本分支冻结包可预留 C 组目录结构，但不提交工业资产与运行。

---

## 5. 各 case 验证清单（摘要）

细节字段仍以 V2 `05_CASE_PARAMETER_TABLES.md` 为准。

### A1 Poiseuille
- 输入：通道、dp、rho/mu/c、Re、入口、探针
- 输出：u(y)、中心 u(t)、流量、质量收支、步数
- 比较：解析 V；三精度 P；三档 R

### A2 Oscillating beam
- 输入：几何、厚度层数、本构、激励、阻尼
- 输出：尾端位移/速度、频率、振幅、能量
- 比较：文献 V；P；R；可选 GPU64

### A3 Window frame
- 输入：材料区、导热、Robin、稳态阈值、COMSOL 配置
- 输出：测点温度、截面分布、内外热流、残差
- 比较：COMSOL V；P；界面细化 R

### A4 Two-phase dambreak
- 输入：几何、两相物性、重力、相间模型、探针
- 输出：界面、水位/前沿、压力峰与到达时间、质量/机械能
- 比较：实验/文献 V；P；R

### A5 FSI2（后置）
- 输入：圆柱+柔性尾板、流固 dp、耦合/子步、初态
- 输出：尾端轨迹、频率/相位、压力/黏性载荷
- 比较：标准 FSI 文献 V；P；R

### B1 Level-set + relaxation
- 固定：几何 hash、单分辨率、松弛规则
- 每档：dp、N、阶段时间、质量/迭代、内存
- 比较：R/G/H/B/C

### B2 3D dambreak
- 固定：物理窗口、物性、排序/关系政策
- 每档：N、邻居、solver 时间、步数、MPPS/MNIPS、内存
- 比较：R/G/H/B/C

### B3 Twisting column
- 固定：本构、初角速度、窗口
- 每档：固体 N、关系重建、结构步、能量/端部运动、时间与内存
- 比较：R/G/H/B/C

---

## 6. 加速比与发布规则（本分支强调）

1. 主加速比：`S = T_CPU_FP32 / T_GPU_FP32`，同 case / 输入 / 规模 / 窗口 / 算法 / 输出政策。
2. 本机 CPU 与本机 GPU 可直接配对；与外卡比较时表中必须列 host，并视为系统级对比。
3. 禁止：旧论文归一化时间与新结果混比；用线程扫描最快单次对 GPU 三次中位数不加说明比较；各卡私自最优 group 冒充统一基线。
4. 未跑完写 `pending`；OOM / 超预算写原因；不用 0 填未测。
5. A 组参考缺失时不把 CPU 结果冒充独立物理真值。

---

## 7. 近期下一步（进入执行时）

1. Phase 0：本机硬件与 `source_inventory` 盘点（命令仍交用户执行核对）。  
2. Phase 1：三套构建 + runner/preflight；先落地 A1 的 `test_2d_gpu_paper_a1_poiseuille` 副本。  
3. Phase 2：从 A1 粗档打通精度与计时模板（tmux + 归档）。  
4. 并行：确认 A3 的 D4/D7 与 COMSOL 数据是否齐全。  
5. B1 几何：做到该阶段再接收用户几何并冻结。  
6. A5、C 组按上文后置 / 外置，不插入当前关键路径。

---

## 8. 相对 V2 变更一览

| 项 | V2 | V3（本分支） |
|---|---|---|
| 主开发机 | 3090 Ti | 本机 4090（单 UUID）+ 本机 64 核 |
| CPU 正式基线 | 独立 64 核服务器 | 本机 64 核；外机仅可选复核 |
| 3090 服务器 | 先导 | 可选跨卡 |
| C1–C3 | 本包主线 | 另分支；本计划不跑 |
| A5 | A 组同步 | P2，主线完成后移植 |
| B1 几何 | 任务书内选定 | 用户后传，做到再冻 |
| 外卡任务 | B+C COMMON | 本分支仅 B COMMON（无 C） |
| 用哪张卡 | 单 UUID 即可 | **强制 GPU 1** |
| 谁跑仿真 | Cursor 本机可跑 | **用户执行**；助手只给命令 |
| Case 路径 | 可用原始 tests | 仅 `extra_source_and_tests` 下 `gpu_paper` 副本 |
| 结果位置 | 未统一 | `docs/gpu_paper_discussion/results/` |
| 断线 | 未强调 | tmux `gp_*` + `logs/` + `check_runs.sh` |
