# A4 二维两相溃坝 — 运行记录

最后更新：2026-11-04（A4 L1–L3 正式 validate 矩阵完成）

## 约定（Duration V2）

| 项 | 值 |
|---|---|
| 终点 | **End_Time=10**（H=g=1 → tau=t=10） |
| Schedule | `schedules/A4_validate_T10_bound.json`（0–10×100 + 冲击窗 5.5–7.5×21；VTP 10 帧） |
| 分辨率 | `GPU_PAPER_DH_OVER_DP` = **LH/dp**：L1=20→dp=0.05，L2=40，L3=80 |
| VTP | 仅 GPU L2 |
| API | **CK / SYCL**（自 `tests_sycl/.../test_2d_two_phase_dambreak_ck`） |
| 几何 | **当前源码** DL=5.3H、DH=2H、探针 (DL,0.2)；**非**论文 5.366H/2.5H/0.19H |
| 二进制 | `test_2d_gpu_paper_a4_two_phase_dambreak` |

## 监测

- 压力：墙角 observer 核插值 `Pressure`
- 水机械能：`TotalMechanicalEnergyCK`
- 各相质量、水波前 `x_front`、条带水位 `h1`(x=LL/2)/`h2`(x=DL/2)、`vmax`、`rho_dev_max`
- 论文量规位置与几何匹配：`missing`（跨端一致性仍可跑）

## 正式 validate 矩阵

| run_id | 状态 | 墙时(s) | 备注 |
|---|---|---:|---|
| `A4_L1_gpu_fp32_validate_r01` | **ok** EXIT:0 | 6.55 | LH/dp=20；终时≈9.89→10 |
| `A4_L1_cpu_fp32_validate_r01` | **ok** EXIT:0 | 5.26 | L1 CPU 对照 |
| `A4_L1_cpu_fp64_validate_r01` | **ok** EXIT:0 | 5.27 | L1 CPU 对照 |
| `A4_L2_gpu_fp32_validate_r01` | **ok** EXIT:0 | 11.16 | LH/dp=40；VTP=1 |
| `A4_L2_cpu_fp32_validate_r01` | **ok** EXIT:0 | 15.53 | L2 CPU 对照 |
| `A4_L2_cpu_fp64_validate_r01` | **ok** EXIT:0 | 15.79 | L2 CPU 对照 |
| `A4_L3_gpu_fp32_validate_r01` | **ok** EXIT:0 | 27.94 | LH/dp=80 |
| `A4_L3_cpu_fp32_validate_r01` | **ok** EXIT:0 | 61.68 | L3 CPU 对照 |
| `A4_L3_cpu_fp64_validate_r01` | **fail** EXIT:139 | — | Time≈9.82 SIGSEGV；保留作失败样本 |
| `A4_L3_cpu_fp64_validate_r02` | **ok** EXIT:0 | 66.05 | 正式采用；终时≈9.998→10 |

## 状态

**A4 正式 L1–L3 × GPU/CPU FP32/CPU FP64 validate 已完成**（L3 CPU64 以 r02 为准；r01 segfault 不采信）。几何仍为 CK 源码非论文匹配；外部文献定量误差另冻配置后才报。
