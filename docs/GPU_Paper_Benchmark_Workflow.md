# SPHinXsys SYCL/GPU 论文 Benchmark 工作说明

分支：`paper/sycl-gpu-benchmarks`  
详细规划见：[SPHinXsys_SYCL_GPU_Paper_Benchmark_Plan.md](./SPHinXsys_SYCL_GPU_Paper_Benchmark_Plan.md)

本文档汇总：**论文要回答什么问题**、**Cursor/GPT 已完成的实现任务**、**当前仓库结构**、**如何编译与正式跑 GPU 数据**、**2080S / 3090 / 4090 采集流程**，以及**尚未完成项**。

---

## 1. 论文实验要回答的三个问题

| 主题 | 内容 | 本仓库负责范围 |
|------|------|----------------|
| **数值验证 / 后端一致性** | SYCL/GPU 是否与预期物理解一致；CPU 与 GPU 是否数值等价 | 六个 verification case + 可选 `OUTPUT=on` 回归轮 |
| **计算扩展性** | 粒子数增加时性能如何变化；2080S / 3090 / 4090 / CPU 如何对比 | 3D dam-break `s1`–`s6` scaling + 统一 `summary.csv` |
| **工业规模演示** | 复杂多物理真实算例 | **不在本分支**；后续在 SPHinXsys-EinSIMO 单独分支 |

当前主仓库工作聚焦前两项：**验证用例 + dam-break 扩展性**。

---

## 2. 六个核心 Verification Case

| ID | Case | 目录 | 分辨率（论文 campaign） | 默认 full 物理时间 |
|----|------|------|-------------------------|-------------------|
| V1 | Mixed Poiseuille 2D | `verification/poiseuille_2d` | coarse / standard / fine | 2.0 s |
| V2 | Dam-break 3D | `verification/dambreak_3d` | standard | 20.0 s |
| V3 | Fish FSI 2D | `verification/fish_fsi_2d` | standard | 1.7 s |
| V4 | Oscillating beam 2D | `verification/oscillating_beam_2d` | standard | 1.0 s |
| V5 | Twisting column 3D | `verification/twisting_column_3d` | standard | 0.5 s |
| V6 | Diffusion Neumann 2D | `verification/diffusion_neumann_2d` | coarse / standard / fine | 1.0 s |

Scaling（扩展性，与 V2 同一 executable）：

| 标签 | dp | 目标粒子量级（约） |
|------|-----|-------------------|
| s1 / standard | 0.050 | ~0.1–0.2 M |
| s2 | 0.040 | ~0.4–0.5 M |
| s3 | 0.032 | ~1 M |
| s4 | 0.025 | ~2 M |
| s5 | 0.020 | ~4 M |
| s6 | 0.016 | ~6–10 M（大卡可跑；2080S 可能 OOM，应如实记录） |

---

## 3. GPT / Cursor 任务清单与完成状态

对应规划文档 **Section 14 (Task A–G)** 及后续对齐工作。

| 任务 | 内容 | 状态 |
|------|------|------|
| **A** | 创建 `paper/sycl-gpu-benchmarks` 分支 | ✅ 已建分支 |
| **B** | 从 `tests/tests_sycl` 盘点并**独立复制**六个 case 到 `tests/gpu_paper_benchmarks` | ✅ 不改动原 tests_sycl |
| **C** | 共享 recorder：`summary.csv` + `environment.json` | ✅ `tests/paper_benchmarks_common/` |
| **D** | CLI：`--benchmark` `--dp` `--resolution` `--end-time` `--output` `--result-dir` `--run-id` 等 | ✅ |
| **E** | 六个 case 接入 benchmark；Poiseuille/Diffusion 三档分辨率 | ✅ |
| **F** | Dam-break `s1`–`s6` scaling | ✅ |
| **G** | 可复现脚本 + `collect_benchmark_results.py` | ✅ |
| **扩展** | 统一 schema：`mpps`、`particle_updates`、`solid_steps`、`peak_rss_kb`、组件计时占位列 | ✅ |
| **扩展** | `--benchmark` + `OUTPUT=off` 时仍做**轻量** observer/feature 采样 | ✅ |
| **扩展** | CPU 对称树 `tests/cpu_paper_benchmarks` + `SPHINXSYS_BUILD_CPU_PAPER_BENCHMARKS` | ✅ 结构就绪，**尚未跑 CPU campaign** |
| **扩展** | 正式 GPU campaign 脚本 `run_gpu_paper_campaign.sh`（48 runs，跑完自动退出） | ✅ |
| **待做** | CPU baseline（128 线程 host）正式数据 | ⏳ |
| **待做** | 数值核对专用轮（`OUTPUT=on`、非 benchmark、回归/DTW） | ⏳ 与性能轮分开 |
| **待做** | 峰值 GPU 显存、MNIPS/GPIPS（需库侧 API 或外部 profiling） | ⏳ 当前列为空/unavailable |
| **待做** | EinSIMO 工业算例 | ⏳ 另一仓库 |

---

## 4. 仓库目录结构

```text
tests/
├── paper_benchmarks_common/     # GPU/CPU 共用 recorder + CLI
│   ├── benchmark_config.h
│   └── benchmark_recorder.h
├── gpu_paper_benchmarks/        # GPU 论文 benchmark（SYCL build）
│   ├── verification/            # 六个 case
│   ├── scripts/                 # 启动与汇总脚本
│   └── curated/                 # 审阅后的汇总 CSV（可提交 git）
└── cpu_paper_benchmarks/        # CPU 对称树（host build，结构同 GPU）

build-sycl/benchmark_results/    # 原始逐 run 数据（gitignore，不提交）
docs/
├── SPHinXsys_SYCL_GPU_Paper_Benchmark_Plan.md   # 完整规划
└── GPU_Paper_Benchmark_Workflow.md              # 本文档
```

CMake 选项：

- `SPHINXSYS_BUILD_GPU_PAPER_BENCHMARKS=ON` → 目标 `gpu_paper_*`
- `SPHINXSYS_BUILD_CPU_PAPER_BENCHMARKS=ON` → 目标 `cpu_paper_*`

目标**不注册 CTest**；长跑由用户手动启动。

---

## 5. 编译（GPU）

```bash
cd /path/to/SPHinXsys   # 本 worktree，分支 paper/sycl-gpu-benchmarks

# 若未加载 oneAPI：
# source /opt/intel/oneapi/setvars.sh --force

cmake -S . -B build-sycl \
  -DSPHINXSYS_USE_SYCL=ON \
  -DSPHINXSYS_BUILD_GPU_PAPER_BENCHMARKS=ON

cmake --build build-sycl --target \
  gpu_paper_poiseuille_2d \
  gpu_paper_diffusion_neumann_2d \
  gpu_paper_dambreak_3d \
  gpu_paper_fish_fsi_2d \
  gpu_paper_oscillating_beam_2d \
  gpu_paper_twisting_column_3d -j
```

可执行文件位置：

```text
build-sycl/tests/gpu_paper_benchmarks/verification/<case>/bin/gpu_paper_<case>
```

---

## 6. 运行模式说明（重要）

### 6.1 性能 / 扩展性数据（写文章主表、图用）

- `--benchmark`（脚本默认通过 launcher 传入）
- `END_TIME_MODE=full`（各 case 使用上表物理结束时间）
- **`OUTPUT=off`**：关闭重 I/O（VTP、restart 等），避免 wall time / MPPS 被磁盘写放大；**不是**不做验证采样
- 轻量 observer/feature（鱼尾、位移、压力曲线等）在 benchmark 模式下**仍会写入** `output/` 并记入 `summary.csv`

### 6.2 数值核对 / 回归（另开一轮，可选）

- **不加** `--benchmark`（verification mode）
- **`OUTPUT=on`**：启用回归 `testResult()` / DTW、完整 observer 输出
- 建议 `REPETITIONS=1`，standard 分辨率即可
- 与性能 campaign **分开目录**，避免混淆

---

## 7. 正式 GPU 论文 Campaign（推荐）

**一条命令跑齐写文章用的 GPU 数据**（verification + dam-break scaling，默认各重复 3 次 = **48 runs**）。

### 7.1 RTX 2080 Super

```bash
cd /path/to/SPHinXsys

BUILD_DIR=/path/to/build-sycl \
ONEAPI_DEVICE_SELECTOR=cuda:gpu \
SPH_BENCH_DEVICE=rtx-2080s \
SPH_BENCH_HOST="$(hostname)" \
SPH_BENCH_OS="$(uname -sr)" \
REPETITIONS=3 \
END_TIME_MODE=full \
OUTPUT=off \
COLLECT=1 \
./tests/gpu_paper_benchmarks/scripts/run_gpu_paper_campaign.sh
```

### 7.2 RTX 3090 / RTX 4090

同一脚本，仅改设备标签（及该机器上的 `BUILD_DIR`）：

```bash
SPH_BENCH_DEVICE=rtx-3090 ./tests/gpu_paper_benchmarks/scripts/run_gpu_paper_campaign.sh
SPH_BENCH_DEVICE=rtx-4090 ./tests/gpu_paper_benchmarks/scripts/run_gpu_paper_campaign.sh
```

### 7.3 Campaign 覆盖范围

| 阶段 | 内容 | 次数（REPETITIONS=3） |
|------|------|----------------------|
| Phase 1 | Verification：上表 10 个 (case, resolution) 组合 | 10 × 3 = 30 |
| Phase 2 | Dam-break scaling：s1 … s6 | 6 × 3 = 18 |
| **合计** | | **48** |

脚本**前台运行**，结束后自动回到 bash（无需 Ctrl+C）。  
兼容旧名：`run_2080s_full_gpu_campaign.sh` → 调用同一 campaign。

### 7.4 输出位置

| 类型 | 路径 |
|------|------|
| 原始逐 run | `build-sycl/benchmark_results/<device>-paper-<UTC>/` |
| Campaign 日志 | 同上目录下 `campaign.log` |
| 汇总 CSV | `tests/gpu_paper_benchmarks/curated/<campaign-id>_all_runs.csv` |
| 重复统计 | `tests/gpu_paper_benchmarks/curated/<campaign-id>_repeat_stats.csv` |

每个 run 目录含：`summary.csv`、`environment.json`、`run.log`、`run.status`、`output/` 等。

### 7.5 查看进度（另开终端）

```bash
tail -n 50 build-sycl/benchmark_results/rtx-2080s-paper-*/campaign.log
find build-sycl/benchmark_results/rtx-2080s-paper-* -name summary.csv | wc -l
```

---

## 8. 分步运行（可选）

若不想一次跑 48 runs，可拆开：

```bash
# 仅 verification（30 runs @ REPETITIONS=3）
BUILD_DIR=build-sycl ONEAPI_DEVICE_SELECTOR=cuda:gpu SPH_BENCH_DEVICE=rtx-2080s \
  REPETITIONS=3 END_TIME_MODE=full OUTPUT=off \
  ./tests/gpu_paper_benchmarks/scripts/run_verification_gpu.sh

# 仅 dam-break scaling（18 runs @ REPETITIONS=3）
BUILD_DIR=build-sycl ONEAPI_DEVICE_SELECTOR=cuda:gpu SPH_BENCH_DEVICE=rtx-2080s \
  REPETITIONS=3 END_TIME_MODE=full OUTPUT=off \
  ./tests/gpu_paper_benchmarks/scripts/run_dambreak_scaling_gpu.sh
```

快速 schema 冒烟（短物理时间，10 runs × 1）：

```bash
./tests/gpu_paper_benchmarks/scripts/run_gpu_schema_smoke.sh
```

---

## 9. 单 case 手动运行示例

```bash
cd build-sycl/tests/gpu_paper_benchmarks/verification/dambreak_3d/bin

ONEAPI_DEVICE_SELECTOR=cuda:gpu \
SPH_BENCH_BACKEND=sycl \
SPH_BENCH_DEVICE=rtx-2080s \
./gpu_paper_dambreak_3d \
  --benchmark \
  --resolution=s3 \
  --end-time=20 \
  --output=off \
  --result-dir=/path/to/build-sycl/benchmark_results \
  --run-id=manual-s3-01
```

---

## 10. 结果汇总

Campaign 脚本在 `COLLECT=1` 时会自动调用；也可手动：

```bash
python3 tests/gpu_paper_benchmarks/scripts/collect_benchmark_results.py \
  build-sycl/benchmark_results/rtx-2080s-paper-<UTC> \
  tests/gpu_paper_benchmarks/curated/rtx-2080s-paper-<UTC>_all_runs.csv \
  --stats-output tests/gpu_paper_benchmarks/curated/rtx-2080s-paper-<UTC>_repeat_stats.csv
```

`summary.csv` 固定列含：`case`、`device`、`mpps`、`particle_updates`、`compute_seconds`、`peak_rss_kb` 等（共 49 个固定列 + case 特有 extra 列）。详见 `tests/paper_benchmarks_common/benchmark_recorder.h`。

---

## 11. 硬件采集矩阵（论文目标）

在同一 commit、同一脚本协议下，于各平台跑 **Section 7** 的正式 campaign：

| 平台 | `SPH_BENCH_DEVICE` 建议 | Build |
|------|-------------------------|-------|
| AMD 128 线程 CPU | `host_cpu` | `build-host`，`cpu_paper_*`，`SPHINXSYS_USE_SYCL=OFF` |
| RTX 2080 Super | `rtx-2080s` | `build-sycl` |
| RTX 3090 | `rtx-3090` | `build-sycl` |
| RTX 4090 | `rtx-4090` | `build-sycl` |

三台 NVIDIA GPU 用**相同** 48-run 协议；CPU 用对称 `tests/cpu_paper_benchmarks/scripts/`（待正式跑）。

大分辨率若 OOM：在 `summary.csv` 的 `status` / `failure_reason` 中记录，**不要**为 2080S 单独缩小 s6 定义。

---

## 12. Git 与数据提交策略

**建议提交：**

- `tests/gpu_paper_benchmarks/`、`tests/cpu_paper_benchmarks/`、`tests/paper_benchmarks_common/`
- `docs/SPHinXsys_SYCL_GPU_Paper_Benchmark_Plan.md`、`docs/GPU_Paper_Benchmark_Workflow.md`
- `tests/gpu_paper_benchmarks/curated/*_all_runs.csv`（审阅后的小体积汇总）

**不要提交：**

- `build-*/benchmark_results/` 下原始 run 树（已在 `.gitignore`）

推送到 GitHub 示例：

```bash
git add docs/GPU_Paper_Benchmark_Workflow.md \
        docs/SPHinXsys_SYCL_GPU_Paper_Benchmark_Plan.md \
        tests/gpu_paper_benchmarks tests/cpu_paper_benchmarks tests/paper_benchmarks_common \
        tests/CMakeLists.txt .gitignore

git status   # 确认无 build 产物、无 secrets
git commit -m "docs: add GPU paper benchmark workflow and shared recorder schema"
git push -u origin paper/sycl-gpu-benchmarks
```

---

## 13. 相关文档

| 文件 | 说明 |
|------|------|
| [SPHinXsys_SYCL_GPU_Paper_Benchmark_Plan.md](./SPHinXsys_SYCL_GPU_Paper_Benchmark_Plan.md) | 论文规划全文（图表、指标定义、Task A–G） |
| [tests/gpu_paper_benchmarks/README.md](../tests/gpu_paper_benchmarks/README.md) | Case 树与 CLI 说明 |
| [tests/gpu_paper_benchmarks/scripts/README.md](../tests/gpu_paper_benchmarks/scripts/README.md) | 环境变量与 launcher 细节 |

---

## 14. 下一步建议顺序

1. **2080S**：跑完 `run_gpu_paper_campaign.sh`（48 runs）→ 检查 `curated/*_all_runs.csv`
2. **3090 / 4090**：同协议各跑一遍
3. **CPU**：在 `build-host` 编译 `cpu_paper_*`，跑对称 campaign（待 baseline 数据）
4. **数值核对轮**（可选）：`OUTPUT=on`、非 benchmark，补 regression / L2 类指标
5. 根据 OOM 结果决定是否调整 s6 目标或论文中单独说明 2080S 上限
6. EinSIMO 工业算例：另开仓库分支，复用同一 summary schema
