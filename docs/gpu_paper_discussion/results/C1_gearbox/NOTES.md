# C1 工业齿轮箱（EinSIMO reducer）— 运行记录

最后更新：2026-11-09

> **跨设备交接请先读** [`C1_DEVICE_HANDOFF.md`](./C1_DEVICE_HANDOFF.md)  
> 冻结配置：[`C1_FROZEN_CONFIG.json`](./C1_FROZEN_CONFIG.json) · 汇总：[`SUMMARY_TABLE.csv`](./SUMMARY_TABLE.csv)

## 约定

| 项 | 值 |
|---|---|
| 源仓库 | EinSIMO `yongchuan/reducer`；就地跑 |
| 工况 | `case_id=8`，默认 **dp=0.0012**，T=**1.05 s≈100 转** |
| contact factor | **关** |
| **粒子** | **固体必须先 relax**（Gear/Wall）；按分辨率存 `reload/dp*/`；物理跑 **reload**。Oil 运行时 Lattice。 |
| relax 计时 | **单独进程**，不进论文 `timer_total` |
| VTP | benchmark 末周期 4 象限；VTP 墙时从 `timer_total` 刨除 |
| 性能 | ×3，取 **median** `wall_solver_s` |

## 重要：曾用错 reload

此前 `bin/reload` 是 **dp≈0.008 粗粒子**（Gear1≈155），与 case8/0.0012 正式点（Gear1=78348）不符。  
已换入匹配正式点的 **dp=0.0012** 固体 reload，并存档于  
`results/C1_gearbox/reload/dp0p0012/`。

## 流程

```bash
export GPU_PAPER_ROOT=/path/to/SPHinXsys-gpu-paper
export EINSIMO_ROOT=/path/to/SPHinXsys-EinSIMO-gearbox-paper
export CUDA_VISIBLE_DEVICES=1 RESOLUTION=0.0012

# 1) 某分辨率若还没有 relax 粒子：
bash $GPU_PAPER_ROOT/docs/gpu_paper_discussion/scripts/run_c1.sh gpu_fp32 relax

# 2) GPU 矩阵
SKIP_SMOKE=1 bash $GPU_PAPER_ROOT/docs/gpu_paper_discussion/scripts/run_c1_gpu_matrix.sh
```

## GPU 矩阵结果（dp=0.0012，本机 4090 GPU1）

| run_id | wall_solver_s | EXIT |
|---|---:|---|
| smoke | （通路径，不作性能） | — |
| benchmark r01 | 26511.3 | 0 |
| benchmark r02 | 26477.2 | 0 |
| benchmark r03 | 26433.7 | 0 |
| **median** | **26477.2** | |

## 状态

| 项 | 状态 |
|---|---|
| dp0p0012 reload 存档 | **就绪**（Gear1/2/3/Wall 与正式点一致） |
| GPU FP32 ×3 | **完成**（2026-11-09） |
| VTP 文件进归档 | **缺失**（有 `selected_snapshot_times.txt`，无 `.vtp`） |
| CPU FP32 ×3 | **未跑** |
