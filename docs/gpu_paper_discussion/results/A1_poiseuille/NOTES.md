# A1 二维 mixed Poiseuille — 运行记录

最后更新：2026-11-04

## 约定

| 项 | 值 |
|---|---|
| 物理终点 | **T = 4 s**（Duration V2） |
| Schedule | `schedules/A1_validate_T4_bound.json`（uniform100 + VTP 0.2/0.5/1/2/3/4） |
| 分辨率 | `GPU_PAPER_DH_OVER_DP`：L1=20，L2=40，L3=80 |
| VTP | 仅 **GPU FP32 L2** 写代表场；其余 `GPU_PAPER_ENABLE_VTP=0` |
| 设备 | GPU：`CUDA_VISIBLE_DEVICES=1`；CPU：本机 Threadripper |
| 二进制 | `test_2d_gpu_paper_a1_poiseuille`（GPU=`build/`，CPU=`build_cpu_fp32/` / `build_cpu_fp64/`） |

## 正式 validate 矩阵（均 EXIT:0，速度剖面 15/15 通过）

| run_id | 墙时(s) | max 相对误差 | mean 相对误差 | 备注 |
|---|---:|---:|---:|---|
| `A1_L1_gpu_fp32_validate_r01` | 4.59 | 0.00338 | 0.00091 | |
| `A1_L1_cpu_fp32_validate_r01` | 17.87 | 0.00725 | 0.00508 | |
| `A1_L1_cpu_fp64_validate_r01` | 17.82 | 0.00731 | 0.00462 | |
| `A1_L2_gpu_fp32_validate_r01` | 9.80 | 0.00582 | 0.00444 | VTP 已拷到 `runs/.../vtp/`（12 文件） |
| `A1_L2_cpu_fp32_validate_r01` | 56.83 | 0.00784 | 0.00587 | |
| `A1_L2_cpu_fp64_validate_r01` | 57.70 | 0.00666 | 0.00492 | |
| `A1_L3_gpu_fp32_validate_r01` | 22.47 | 0.00670 | 0.00424 | |
| `A1_L3_cpu_fp32_validate_r01` | 210.94 | 0.00441 | 0.00293 | |
| `A1_L3_cpu_fp64_validate_r01` | 215.90 | 0.00399 | 0.00282 | |

日志：`docs/gpu_paper_discussion/logs/gp_A1_L{1,2,3}_*.log`（L1 GPU 历史名为 `gp_A1_L1_gpu32.log`）。

## 已知注意

- CPU FP32 全量 `extra_src` 遗留 UDF 在 icpx 下会挂；本分支默认只编 `shared/gpu_paper`。
- 运行脚本须用 `tee` 把 `N=`/`Time=` 打到终端；勿在 `set -u` 下 source oneAPI `setvars`。
- 进出口累计质量通量仍标 `not_collected`（见各 run 的 `missing_metrics.csv`）。

## 状态

**A1 正式 L1–L3 × GPU/CPU FP32/CPU FP64 validate 已完成。** 下一案例：A2。
