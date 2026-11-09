# B3 三维扭转柱 — 运行记录

最后更新：2026-11-07（L1–L5 × 三端 benchmark r01 矩阵已补全）

## 约定（Duration V2）

| 项 | 值 |
|---|---|
| 终点 | **End_Time=0.5 s** |
| 物性 | 源码：ρ=1100，E=17e6，**ν=0.45**，ω0=-300（与 ENOG ν=0.49 不同，须声明） |
| 分辨率 | `GPU_PAPER_DH_OVER_DP`=**PH/dp**：L1=5，L2=10（源码），L3=20，L4=40，L5=80 |
| benchmark | `B3_benchmark_T0p5_bound.json`（0–0.5 步长 0.01，51 点）；关 VTP；正式×3 另议 |
| validate/showcase | `B3_validate_T0p5_bound.json`（100 点）；VTP 仅 showcase GPU |
| 探针 | 端心 (PL,0,0) + 离轴 (PL,0.35,0)；截面转角拟合仍 missing |
| 二进制 | `test_3d_gpu_paper_b3_twisting_column` |

## 正式矩阵（benchmark r01）

|  | GPU FP32 | CPU FP32 | CPU FP64 |
|---|---|---|---|
| L1 | **ok** | **ok** | **ok** |
| L2 | **ok** | **ok** | **ok** |
| L3 | fail collapse | fail collapse | **ok** |
| L4 | fail collapse | fail collapse | fail collapse |
| L5 | fail collapse | fail collapse | fail collapse |

### 本轮补跑（L4/L5 CPU）

| run_id | 状态 | 墙时(s) | 塌缩约 t |
|---|---|---:|---|
| `B3_L4_cpu_fp32_benchmark_r01` | timestep_collapse | 46.5 | ~0.08–0.09（`acoustic_dt=nan`） |
| `B3_L4_cpu_fp64_benchmark_r01` | timestep_collapse | 85.9 | ~0.084（`-nan`） |
| `B3_L5_cpu_fp32_benchmark_r01` | timestep_collapse | 405 | ~0.048（`nan`） |
| `B3_L5_cpu_fp64_benchmark_r01` | timestep_collapse | 1270 | ~0.047（`acoustic_dt→0`） |

说明：L4/L5 即使用 **CPU64** 仍 collapse，故细档问题不只是 FP32 舍入；L3 仅 CPU64 能跑满 0.5 s。

## 状态

**15/15 矩阵格已有 r01 记录。**  
可用性能/正确性主表：**L1–L2 三端 + L3 CPU64**；L3 FP32 与 L4–L5 全端标 **blocked（timestep_collapse）**。×3 与阻尼/CFL 排查另议。
