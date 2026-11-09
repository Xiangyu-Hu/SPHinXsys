# GPU 论文实验讨论区

本目录存放相对 V2 设备任务书的**本地重规划、讨论结论与执行进度**。  
原始任务包仍保留在：`../GPU_PAPER_DEVICE_TASKS_V2_2026-11-02.zip`。

## 文件索引

| 文件 / 目录 | 用途 |
|---|---|
| [01_TEST_PLAN_V3.md](01_TEST_PLAN_V3.md) | 修订测试计划（测什么） |
| [02_PROGRESS.md](02_PROGRESS.md) | 执行进度、阻塞项 |
| [03_RUN_CONVENTIONS.md](03_RUN_CONVENTIONS.md) | **强制运行约定**（谁跑、GPU、路径、tmux） |
| [04_SAMPLING_POLICY.md](04_SAMPLING_POLICY.md) | 早期粗采样讨论（**频率以 V4 为准**） |
| [05_V4_REVIEW_NOTES.md](05_V4_REVIEW_NOTES.md) | V4 通读笔记与本地对齐讨论 |
| [06_PRE_RUN_INSTRUMENTATION_GATE.md](06_PRE_RUN_INSTRUMENTATION_GATE.md) | **开跑前门禁**：字段/observer/missing 一次就绪 |
| [07_USER_COMMANDS.md](07_USER_COMMANDS.md) | **请你执行的命令**（编译/短测/tmux） |
| [GPU_PAPER_CASE_DURATION_AND_PHYSICAL_OUTPUTS_V2_2026-11-03.md](GPU_PAPER_CASE_DURATION_AND_PHYSICAL_OUTPUTS_V2_2026-11-03.md) | **物理时长/监测/VTP 权威表** |
| [08_DURATION_OUTPUTS_V2_ALIGNMENT.md](08_DURATION_OUTPUTS_V2_ALIGNMENT.md) | 与本地执行的对齐说明 |
| `GPU_PAPER_DEVICE_TASKS_V4_2026-11-03.zip` | 权威任务包（采样/case 卡/模板） |
| `results/` | 各 case 归档数据（权威副本） |
| `logs/` / `sessions/` | 长跑日志与 tmux 元数据 |
| `scripts/check_runs.sh` | 重连后查看跑点状态 |

## 当前机器（本讨论默认主执行机）

- GPU：**只用 1 号卡** RTX 4090（UUID `GPU-6d9732cb-5aed-bfad-6d66-736689de9900`）；0 号常被占用
- CPU：AMD Ryzen Threadripper PRO 5995WX，64 物理核 / 128 逻辑线程
- 角色：主开发 + 本机 CPU/GPU 配对；助手给命令，**用户执行**

## 相对 V2 的关键变更（摘要）

1. 主开发与本机 CPU/GPU 配对改在本机完成；3090 服务器不稳定，降为可选跨卡点。
2. 工业 case（C1/C2/C3）不在本分支执行，另开分支处理。
3. A5（FSI2 GPU 移植）排到其余内容完成后再做。
4. B1 复杂几何待做到该 case 时由用户提供 STL/几何资产。
5. 论文 case 一律放在 `tests/extra_source_and_tests/.../test_*_gpu_paper_*`，不改原始 `tests` 主算例当跑点。
6. 跑完数据归档到本目录 `results/<CASE>/`；长跑用 tmux `gp_*` + 日志，断线可续看。

更新日期：2026-11-03。
