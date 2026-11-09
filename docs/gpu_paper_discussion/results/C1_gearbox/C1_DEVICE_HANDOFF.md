# C1 工业齿轮箱 — 跨设备交接文档

> **给下一台机器上的 Cursor / 自己读**：打开本文件即可接续，不必重讲整段对话。  
> 最后更新：2026-11-09（本机 `sphinxsys-ml`，RTX 4090 GPU1，已完成 GPU FP32 ×3）

---

## 1. 一句话：我们在干什么

在 **SPHinXsys SYCL GPU 论文** 的工业算例组 **C1** 上，用 EinSIMO 真实减速器（齿轮箱搅油 + 传热）做 **可复现的性能基准**，并与 A 组（2D 验证）/ B 组（3D 规模）放在同一套归档体系里对比。

---

## 2. 为什么跑这个 case（目的）

| 问题 | 答案 |
|---|---|
| 论文里 C 组角色 | **工业规模 / 真实几何**，证明 SYCL-GPU 路径不只在理想基准上可用 |
| 为何选 EinSIMO gearbox | 已有真实 CAD、多刚体 + 油液 + 热模型；与学术 A/B 形成互补 |
| 对比什么 | **同一物理配置**下：GPU FP32 vs CPU FP32（及后续外卡 COMMON）；论文主指标是 **solver 墙时**（`wall_solver_s` / `timer_total`），**不是**输出/VTP 时间 |
| 不对比什么 | 不做 L1–L3 多分辨率矩阵（C 组只用 **一个 COMMON 分辨率**）；不把 relax 时间、VTP 写出时间算进性能 |
| 物理终点 | 冻结 **T = 1.05 s ≈ 100 转**（不是文章叙事里的 2.1 s） |
| Contact factor | **默认关闭**（`enable_contact_factor=0`）：不计算、不记录，避免干扰 flow+heat 主路径与计时 |

---

## 3. 仓库与分支（固定）

| 角色 | 路径（本机） | 远程 | 分支 | Git 推什么 |
|---|---|---|---|---|
| **C1 case + 交接文档（主推送）** | `/home/yongchuan/SPHinXsys-EinSIMO-gearbox-paper` | `Xiangyu-Hu/SPHinXsys-EinSIMO` | `yongchuan/reducer` | case 源码 + `docs/gpu_paper_c1/` |
| 论文 A/B case + 讨论文档 | `/home/yongchuan/SPHinXsys-gpu-paper` | `Xiangyu-Hu/SPHinXsys` | `paper/sycl-gpu-benchmarks` | case + MD/脚本；**不含** runs/reload |
| 计算结果 | 仅本机 `docs/gpu_paper_discussion/results/*/runs` 等 | — | — | **不推送**；SSH 打包下载 |

**原则**：C1 **不搬进** gpu-paper 里编译跑；在 EinSIMO 就地跑。数值结果只留本机，与 A/B 的 `results/` 目录并列存放，但 **不进 Git**。

换机时：

1. clone / pull 上述两个分支  
2. 设环境变量（见 §6）指向本机路径  
3. 编译 EinSIMO SYCL/CPU gearbox case  
4. 用已归档的 `reload/dp0p0012/` 或重新 `relax`

---

## 4. 冻结配置（务必一致）

见同目录 `C1_FROZEN_CONFIG.json`。摘要：

| 项 | 冻结值 |
|---|---|
| `case_id` | **8**（nominal / R-M01） |
| `resolution` / dp | **0.0012**（唯一 COMMON；tag=`dp0p0012`） |
| `end_time` | **1.05** s（≈100 shaft3 转） |
| `enable_contact_factor` | **0** |
| 热模型 | 开（flow + heat；viscous dissipation / oil–solid / solid conduction / wall convection） |
| 固体粒子 | **必须**该 dp 的 **relaxed reload**（Gear1/2/3 + Wall）；**禁止**用错分辨率的粗粒子 |
| 油粒子 | 运行时 **Lattice** 生成（不进 reload 存档） |
| relax 计时 | **单独进程**；**不进**论文 `timer_total` |
| 性能重复 | benchmark **×3**，取 **中位数** `wall_solver_s` |
| VTP | 末 1 转、每 0.25 转一帧（4 象限）；VTP 墙时从 `timer_total` **刨除**，仍记在 output 相关字段 |
| 本机 GPU | `CUDA_VISIBLE_DEVICES=1`（4090；勿用常被占的 GPU0） |

**粒子数（正式点，校验 reload 是否配对）**

| body | N |
|---|---:|
| oil | 755762 |
| gear1 | **78348** |
| gear2 | 211206 |
| gear3 | 293676 |
| wall | 491304 |
| **total** | **1830296** |

若 Gear1 ≈ 155 → **错了**，是旧粗网格；必须换 `reload/dp0p0012/`。

---

## 5. 本机已完成的计算（存档位置）

> **Git 政策**：远程只推 **case 源码 + 记录文档**；`runs/`、`reload/*.xml`、日志等计算结果 **不进 Git**，只留本机 / SSH 打包下载。



根目录：`docs/gpu_paper_discussion/results/C1_gearbox/`

```
C1_gearbox/
├── C1_DEVICE_HANDOFF.md          ← 本文件
├── C1_FROZEN_CONFIG.json
├── NOTES.md
├── SUMMARY_TABLE.csv             ← GPU×3 汇总
├── reload/dp0p0012/              ← 固体 reload（Gear*/Wall_rld.xml）
└── runs/
    ├── C1_case8_gpu_fp32_dp0p0012_benchmark_r01/
    ├── C1_case8_gpu_fp32_dp0p0012_benchmark_r02/
    └── C1_case8_gpu_fp32_dp0p0012_benchmark_r03/
```

日志：`docs/gpu_paper_discussion/logs/gp_C1_case8_gpu_fp32_dp0p0012_benchmark_r0{1,2,3}.log`  
矩阵总日志：`.../logs/gp_C1_gpu_matrix.log`

### GPU FP32 ×3（4090，dp=0.0012，T=1.05，CF=0）

| run | wall_solver_s | actual_t_end | EXIT |
|---|---:|---:|---|
| r01 | 26511.3 | 1.05 | 0 |
| r02 | 26477.2 | 1.05 | 0 |
| r03 | 26433.7 | 1.05 | 0 |
| **median** | **26477.2** | | |

每 run 目录应有：`run_manifest.json`、`run_summary.csv`、`runtime_breakdown.csv`、`stage_timing.csv`、`reducer_thermal_history.csv`、`exit_status.json`、observer/gear action `.dat`、`selected_snapshot_times.txt`。

### 已知缺口（换机时注意）

- `selected_snapshot_times.txt` 已写出末周期 4 个时刻，但 **归档目录里没有 `.vtp` 文件**（需后续查写出路径 / 再跑 showcase 或修归档 glob）。  
- **CPU FP32 ×3 尚未在本机跑完**（脚本已支持 `run_c1.sh cpu_fp32 ...`）。  
- 勿把 `results/B1_stanford_bunny/`（约 13G）当必拉内容；论文用精简的 `B1_particle_relaxation` 即可。

---

## 6. 换机后怎么直接跑同一配置

### 6.1 环境变量

```bash
export GPU_PAPER_ROOT=/path/to/SPHinXsys-gpu-paper          # 本仓库
export EINSIMO_ROOT=/path/to/SPHinXsys-EinSIMO-gearbox-paper
export CUDA_VISIBLE_DEVICES=1                              # 按该机空闲卡改
export RESOLUTION=0.0012
export CASE_ID=8
export ENABLE_CONTACT_FACTOR=0
```

脚本会读上述变量（见 `docs/gpu_paper_discussion/scripts/run_c1.sh`）。

### 6.2 依赖

- Intel oneAPI（`source /opt/intel/oneapi/setvars.sh`；注意它会把 `BIN_DIR` 设成 `bin64`，脚本已改用 `C1_BIN_DIR`）
- CUDA 设备 + 已编译的  
  `$EINSIMO_ROOT/build/tests/tests_extra/test_3d_gearbox_real_sycl/bin/test_3d_gearbox_real_sycl`  
  （CPU：`.../test_3d_gearbox_real_cpu/...`）
- EinSIMO 分支含 GPU-paper 补丁：`enable_contact_factor` 门控、manifest 字段、`GPU_PAPER_RUN_DIR` 归档、VTP 从 `timer_total` 刨除

### 6.3 命令

```bash
cd "$GPU_PAPER_ROOT"

# （可选）若该分辨率还没有固体 reload：
RESOLUTION=0.0012 bash docs/gpu_paper_discussion/scripts/run_c1.sh gpu_fp32 relax

# GPU 矩阵：smoke → benchmark×3（可 SKIP_SMOKE=1）
SKIP_SMOKE=1 RESOLUTION=0.0012 \
  bash docs/gpu_paper_discussion/scripts/run_c1_gpu_matrix.sh
```

单次：

```bash
RUN_SUFFIX=r01 bash docs/gpu_paper_discussion/scripts/run_c1.sh gpu_fp32 benchmark 1.05
RUN_SUFFIX=r01 bash docs/gpu_paper_discussion/scripts/run_c1.sh cpu_fp32 benchmark 1.05
```

### 6.4 从本 SSH 机下载存档（不经过 git 大文件时）

```bash
# 在服务器上打包（示例）
cd /home/yongchuan/SPHinXsys-gpu-paper
tar -czf /tmp/C1_gearbox_archive.tgz \
  docs/gpu_paper_discussion/results/C1_gearbox \
  docs/gpu_paper_discussion/logs/gp_C1_* \
  docs/gpu_paper_discussion/scripts/run_c1.sh \
  docs/gpu_paper_discussion/scripts/run_c1_gpu_matrix.sh

# 本机
scp user@host:/tmp/C1_gearbox_archive.tgz .
```

Git 侧：结果应已落在 `paper/sycl-gpu-benchmarks` 的 `docs/gpu_paper_discussion/results/C1_gearbox/`。

---

## 7. 必须记录的内容（论文 / 复现）

每正式 run 至少保留：

1. **身份**：`run_manifest.json`（case_id、dp、end_time、backend、device、particle_counts、CF、热开关、git/compiler）  
2. **性能**：`run_summary.csv` 的 `wall_solver_s`（主指标）；`runtime_breakdown.csv` / `stage_timing.csv`  
3. **物理/热**：`reducer_thermal_history.csv`；observer 温度与 gear action `.dat`（体状态）  
4. **退出**：`exit_status.json` + 日志末尾 `EXIT:0`  
5. **快照计划**：`selected_snapshot_times.txt`（+ 将来补齐的 VTP）  
6. **矩阵汇总**：×3 的 raw + **median**（见 `SUMMARY_TABLE.csv`）  
7. **reload 指纹**：`reload/dp0p0012/README.txt` + 粒子数与上表一致  

跨硬件对比时：**同一 dp / T / CF / reload**；只改设备与 backend；**不要**拿输出时间比快慢。

---

## 8. 给下一台 Cursor 的最短提示词（可复制）

```
请先读 SPHinXsys（branch paper/sycl-gpu-benchmarks）里
docs/gpu_paper_discussion/results/C1_gearbox/C1_DEVICE_HANDOFF.md
和 C1_FROZEN_CONFIG.json。

我们在做 GPU 论文工业算例 C1（EinSIMO gearbox，case_id=8，dp=0.0012，T=1.05s，
contact factor 关）。结果归档在 results/C1_gearbox，与 A/B 并列。
可执行代码在 EinSIMO 仓库 yongchuan/reducer，就地跑，用
docs/gpu_paper_discussion/scripts/run_c1.sh。

本机 GPU×3 已完成（median wall_solver≈26477s）。下一步按交接文档跑 CPU×3
或外卡 COMMON，并排查 benchmark 未归档到的 VTP。
```

---

## 9. 相关脚本与笔记

| 文件 | 用途 |
|---|---|
| `docs/gpu_paper_discussion/scripts/run_c1.sh` | relax / smoke / benchmark / showcase |
| `docs/gpu_paper_discussion/scripts/run_c1_gpu_matrix.sh` | GPU smoke + ×3 |
| `docs/gpu_paper_discussion/results/C1_gearbox/NOTES.md` | 短运行笔记 |
| `docs/gpu_paper_discussion/02_PROGRESS.md` | 全论文进度 |
| `docs/gpu_paper_discussion/03_RUN_CONVENTIONS.md` | 通用跑法约定 |
