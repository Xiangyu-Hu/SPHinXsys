# A2 二维梁振动 — 运行记录

最后更新：2026-11-04（L1–L3 正式 validate 矩阵完成）

## 约定（Duration V2）

| 项 | 值 |
|---|---|
| 物理终点 | **T = 2 s** |
| Schedule | `schedules/A2_validate_T2_bound.json`（0–2s×100 + 0.8–1.9s 窗×100 合并；VTP 8 帧） |
| 分辨率 | `GPU_PAPER_DH_OVER_DP` = **PH/dp**：L1=10，L2=20，L3=40（源码默认 PH/10） |
| VTP | 仅 GPU FP32 L2 |
| 对照 | ENOG 2023 §4.1；周期/幅值后处理，源算法≠论文算法须声明 |
| 二进制 | `test_2d_gpu_paper_a2_oscillating_beam` |

## 正式 validate 矩阵

| run_id | 状态 | 墙时(s) | 备注 |
|---|---|---:|---|
| `A2_L1_gpu_fp32_validate_r01` | **作废** | — | 2026-11-04：`acoustic_dt→0` 死循环（重算 LinearCorrectionMatrix） |
| `A2_L1_gpu_fp32_validate_r02` | **ok** EXIT:0 | 13.33 | T→2 s；PH/dp=10；VTP=0；tip_y≈±0.113 |
| `A2_L1_cpu_fp32_validate_r01` | **ok** EXIT:0 | 15.15 | 与 GPU L1 对照 |
| `A2_L1_cpu_fp64_validate_r01` | **ok** EXIT:0 | 15.92 | 与 GPU L1 对照 |
| `A2_L2_gpu_fp32_validate_r01` | **ok** EXIT:0 | 27.12 | PH/dp=20；VTP=1；代表场归档见 `runs/.../vtp/` |
| `A2_L2_cpu_fp32_validate_r01` | **ok** EXIT:0 | 45.71 | L2 CPU 对照 |
| `A2_L2_cpu_fp64_validate_r01` | **ok** EXIT:0 | 50.04 | L2 CPU 对照 |
| `A2_L3_gpu_fp32_validate_r01` | **ok** EXIT:0 | 57.01 | PH/dp=40；VTP=0 |
| `A2_L3_cpu_fp32_validate_r01` | **ok** EXIT:0 | 172.12 | L3 CPU 对照 |
| `A2_L3_cpu_fp64_validate_r01` | **ok** EXIT:0 | 210.03 | L3 CPU 对照 |

## 状态

**A2 正式 L1–L3 × GPU/CPU FP32/CPU FP64 validate 已完成。** tip_y 跨端一致（L1≈±0.113，L2/L3≈±0.110）。ENOG 周期/幅值后处理与应变能仍标 `not_collected`。下一案例：A3。

每次正式跑完成后在本表追加一行，结果进 `runs/<run_id>/`，勿覆盖。用新 `run_id`（如 `_r02`）重跑，勿覆盖作废目录除非你明确同意。
