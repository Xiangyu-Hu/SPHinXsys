# A3 窗框传热 D4 — 运行记录

最后更新：2026-11-04（L1–L3 正式 validate 矩阵完成）

## 约定（Duration V2）

| 项 | 值 |
|---|---|
| 变体 | **D4**（D7 备选，本分支先不做） |
| 终点 | **伪时间 End_Time=0.02**（非真实秒）；未收敛再考虑 0.04/0.08 |
| Schedule | `schedules/A3_validate_T0p02_bound.json`（0–0.02×100；VTP 0.004/0.008/0.012/0.016/0.02） |
| 分辨率 | `GPU_PAPER_DH_OVER_DP` = **scale**：L1=1→dp=1mm，L2=2→0.5mm，L3=4→0.25mm |
| VTP | 仅 GPU L2 |
| API | **经典 SPH**（非 CK）；仓库无现成 SYCL D4。GPU 构建能否通过以现场编译为准 |
| 二进制 | `test_2d_gpu_paper_a3_windows_frame_D4` |

## 正式 validate 矩阵

| run_id | 状态 | 墙时(s) | 备注 |
|---|---|---:|---|
| `A3_D4_L1_gpu_fp32_validate_r01` | **ok** EXIT:0 | 2.94 | scale=1；伪时间→0.02 |
| `A3_D4_L1_cpu_fp32_validate_r01` | **ok** EXIT:0 | 2.97 | L1 CPU 对照 |
| `A3_D4_L1_cpu_fp64_validate_r01` | **ok** EXIT:0 | 3.56 | L1 CPU 对照 |
| `A3_D4_L2_gpu_fp32_validate_r01` | **ok** EXIT:0 | 26.31 | scale=2；VTP=1；场归档 `runs/.../vtp/` |
| `A3_D4_L2_cpu_fp32_validate_r01` | **ok** EXIT:0 | 26.42 | L2 CPU 对照 |
| `A3_D4_L2_cpu_fp64_validate_r01` | **ok** EXIT:0 | 41.07 | L2 CPU 对照 |
| `A3_D4_L3_gpu_fp32_validate_r01` | **ok** EXIT:0 | 1191.4 | scale=4；伪时间略超 0.02（~0.022，块积分过冲） |
| `A3_D4_L3_cpu_fp32_validate_r01` | **ok** EXIT:0 | 1168.8 | L3 CPU 对照 |
| `A3_D4_L3_cpu_fp64_validate_r01` | **ok** EXIT:0 | 1984.5 | L3 CPU 对照；终时≈0.0200 |

## 状态

**A3 D4 正式 L1–L3 × GPU/CPU FP32/CPU FP64 validate 已完成。** 经典 SPH 路径上 GPU≈CPU32 墙时（未吃满 GPU）。COMSOL/L2D/Uf 仍 `not_collected`。下一案例：A4。

## SYCL (2026-11-05)
已改为 **DiffusionRelaxationCK + Robin CK**；GPU 构建通过。详见 `NOTES_SYCL.md`。

## CK 对齐重跑（2026-11-06）

| run_id | EXIT | wall_solver (s) | backend |
|---|---:|---:|---|
| A3_D4_L1_gpu_fp32_validate_ck01 | 0 | ≈2.88 | sycl_ck_gpu |
| A3_D4_L1_cpu_fp32_validate_ck01 | 0 | ≈5.41 | cpu_ck |
| A3_D4_L1_cpu_fp64_validate_ck01 | 0 | ≈5.93 | cpu_ck |

L1 三端 `status=ok`，api 均为 `DiffusionRelaxationCK_Robin_LocalIsotropicDiffusion`；GPU/CPU32 探针几乎一致，与 CPU64 差约 1e-3–2e-2 °C（伪时间 0.02 瞬态）。

| run_id | EXIT | wall_solver (s) |
|---|---:|---:|
| A3_D4_L2_gpu/c32/c64_validate_ck01 | 0 | ≈11.5 / 49.8 / 59.2 |
| A3_D4_L3_gpu/c32/c64_validate_ck01 | 0 | ≈83.9 / 719 / 813 |

L2–L3 亦全部 `ok`；GPU↔CPU32 探针一致（或机器精度内），相对 CPU64 仍为小差。**A3 CK 精度矩阵 L1–L3 × 三端完成。**
