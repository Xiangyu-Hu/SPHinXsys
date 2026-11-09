# 运行与协作约定（强制）

日期：2026-11-03  
适用于本仓库 GPU 论文测试（本分支 A/B 组）。违反任一条前先改约定文档，不要 silently 绕过。

---

## 1. 谁执行仿真

- **助手只改代码、写配置、给出完整可复制命令，不自行启动正式/长跑 case。**
- 用户在终端执行；执行后把关键输出路径或退出码反馈即可。
- 短测命令同样发给用户执行（除非用户当次明确说“你直接跑这个短测”）。

---

## 2. GPU 绑定（本机双 4090）

| 项 | 值 |
|---|---|
| 使用设备 | **GPU 1**（不要用 0） |
| 型号 | NVIDIA GeForce RTX 4090 |
| UUID | `GPU-6d9732cb-5aed-bfad-6d66-736689de9900` |
| 原因 | GPU 0 常被占用 |

每条运行命令必须包含其一（优先环境变量，简单明确）：

```bash
export CUDA_VISIBLE_DEVICES=1
# 或单行：
CUDA_VISIBLE_DEVICES=1 <你的可执行文件与参数>
```

启动后应在日志中能核对实际 device（打印 GPU name / UUID / SYCL device）。若误落到 GPU 0，该次 run 作废重跑。

CPU-only 精度/性能跑点：**不要**设置占满两张卡的可见性；GPU 相关进程一律 `CUDA_VISIBLE_DEVICES=1`。

---

## 3. Case 源码位置（禁止改原始 tests 主算例当论文跑点）

| 用途 | 路径 |
|---|---|
| 论文专用 case 根目录 | `/home/yongchuan/SPHinXsys-gpu-paper/tests/extra_source_and_tests` |
| 2D | `.../2d_examples/test_2d_gpu_paper_<slug>/` |
| 3D | `.../3d_examples/test_3d_gpu_paper_<slug>/` |
| 共享实验工具（若需要） | `.../extra_src/` 下新增，名称带 `gpu_paper` |

约定：

1. 从 `tests/` 或 `tests/tests_sycl/` **复制**出独立目录，再加论文所需记录/计时/输出；**不在原始 case 上堆论文插桩。**
2. 目录名、目标名、源文件头注释标明：`GPU paper experiment (docs/gpu_paper_discussion)`。
3. CMake 目标名与文件夹一致，便于 `ninja`/`make` 单目标编译。
4. 已有 `extra_source_and_tests` 里非 `gpu_paper` 前缀的算例视为历史/其他用途，论文矩阵不直接当最终跑点（可参考，但要再复制一份 `gpu_paper` 版）。

建议命名对照：

| 计划 ID | 文件夹名（拟定） |
|---|---|
| A1 | `test_2d_gpu_paper_a1_poiseuille` |
| A2 | `test_2d_gpu_paper_a2_oscillating_beam` |
| A3 | `test_2d_gpu_paper_a3_window_frame` |
| A4 | `test_2d_gpu_paper_a4_two_phase_dambreak` |
| A5 | `test_2d_gpu_paper_a5_fsi2`（后置） |
| B1 | `test_3d_gpu_paper_b1_particle_relaxation` |
| B2 | `test_3d_gpu_paper_b2_dambreak` |
| B3 | `test_3d_gpu_paper_b3_twisting_column` |

---

## 4. 结果归档

每个 run 结束后，把该次运行目录（或约定导出物）归档到：

```text
/home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/results/<CASE_ID>/
```

建议结构：

```text
results/
  A1_poiseuille/
    runs/<run_id>/          # manifest、summary、logs、diagnostics…
    figures/                # 后处理图
    NOTES.md                # 该 case 备注
  A2_oscillating_beam/
  ...
  _registry/
    run_registry.csv        # 一行一次运行
sessions/
  <tmux_session>.meta.md    # 会话名、启动命令、日志路径
logs/
  <run_id>.stdout.log       # 可选：与 runs 内日志二选一或软链
```

规则：

- **不覆盖**：同 case 新跑用新 `run_id` 子目录。
- 失败/中断同样归档 `exit_status` 与日志。
- 原始 build 目录里的输出可以保留；**讨论区 `results/` 才是论文用权威副本。**
- 工业 C 组结果不进本分支该树（另分支）。

---

## 5. 断线可续看：tmux + 日志

本机已有 `tmux`。正式/长跑统一：

1. 固定会话名：`gp_<case>_<tag>`，例如 `gp_a1_l1_gpu32`。
2. 日志同时写入文件（便于不 attach 也能看）。
3. 在 `sessions/` 留下一行元数据（谁启动、命令、日志路径、GPU）。

### 5.1 推荐启动模板（发给用户的命令形态）

```bash
# ---- 用户执行：示例模板 ----
CASE=A1
TAG=l1_gpu32
SESSION=gp_${CASE}_${TAG}
ROOT=/home/yongchuan/SPHinXsys-gpu-paper
LOGDIR=$ROOT/docs/gpu_paper_discussion/logs
RUNDIR=$ROOT/docs/gpu_paper_discussion/results/${CASE}_*/runs   # 实际给出绝对路径
mkdir -p "$LOGDIR" "$RUNDIR_ACTUAL"

tmux new-session -d -s "$SESSION" "
  export CUDA_VISIBLE_DEVICES=1
  cd <build_or_case_workdir>
  <executable> <args> 2>&1 | tee $LOGDIR/${SESSION}.log
  echo EXIT:\$? | tee -a $LOGDIR/${SESSION}.log
"

# 元数据（可复制为 sessions/${SESSION}.meta.md）
```

### 5.2 重连后查看状态

```bash
# 列出论文相关会话
tmux ls | grep '^gp_'

# 进入查看（只读看可用 attach；不影响进程）
tmux attach -t gp_A1_l1_gpu32
# 分离：Ctrl-b d

# 不进入会话，只看尾部日志
tail -n 80 /home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/logs/gp_A1_l1_gpu32.log

# GPU 是否在 1 号卡上跑
nvidia-smi -i 1

# 辅助脚本
bash /home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/scripts/check_runs.sh
```

SSH 断开 **不会** 杀掉 tmux 内进程；重新登录后用上面命令即可。

若用户更习惯 `screen`，会话名同样用 `gp_*`，但默认文档与助手给的命令以 **tmux** 为准。

---

## 6. 助手交付命令时的检查清单

发给用户的每条运行说明应包含：

1. `CUDA_VISIBLE_DEVICES=1`
2. 工作目录（绝对路径）
3. 可执行文件与参数
4. tmux 会话名
5. 日志文件绝对路径
6. 预期结果归档目录
7. 如何确认仍在跑 / 已结束（`tmux ls`、`tail`、`EXIT:` 行）

---

## 7. 与计划文档的关系

- 测什么：`01_TEST_PLAN_V3.md`
- 做到哪：`02_PROGRESS.md`
- 怎么跑/放哪：本文  
原始 V2 数据契约（计时字段、manifest）仍然适用；路径与执行权以本文为准。
