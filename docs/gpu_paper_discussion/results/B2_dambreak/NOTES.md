# B2 三维溃坝 — 运行记录

最后更新：2026-11-08（CPU L3 ×3：r01/r02 ok；r03/r04=`non_finite_state`）

## 约定（Duration V2）

| 项 | 值 |
|---|---|
| 终点 | **End_Time=10**（H=g=1 → tau=t） |
| 几何 | CPC 附录匹配：DL=5.366、DH=2、DW=0.5；水柱 2×1×0.5；探针含 (DL,0.2,0.5 DW) |
| 分辨率 | `GPU_PAPER_DH_OVER_DP`=**LH/dp**：L1=10→0.1，L2=20→0.05，L3=40，L4=80，L5=160 |
| benchmark | `B2_benchmark_T10_bound.json`（0–10 步长 0.2，51 点）；**关 VTP**；正式×3 |
| validate/showcase | `B2_validate_T10_bound.json`（100+冲击加密）；VTP 仅 **showcase GPU** |
| 二进制 | `test_3d_gpu_paper_b2_dambreak` |

## 正式矩阵（待跑）

| run_id | 状态 | 墙时(s) | 备注 |
|---|---|---:|---|
| `B2_L2_gpu_fp32_smoke_r01` | **ok** EXIT:0 | 0.07 solver / 0.34 process | T=0.2；终时≈0.202；N_active=8000 |
| build gpu/cpu32/cpu64 | **ok** EXIT:0 | — | 三端二进制就绪 |
| `B2_L3_gpu_fp32_benchmark_r01` | **ok** EXIT:0 | **11.85** | LH/dp=40；终时≈9.62→10；先导可接受 |

## 状态

GPU L1–L5 / CPU L1–L3 首轮已 ok。L4 CPU 曾 `timestep_collapse`/假死（能量 NaN）。  
**2026-11-05**：L4 CPU `r02`（加水相 BoundingFromWallCK）在 t≈4.80 **EXIT:139 segfault**。  
GPU L4 `r01` 成功用的是**打补丁前**旧二进制，不能直接对比。  
已撤回单相 wall bounding（§2.7 针对轻相）；保留 dt/NaN/子步 abort。L4/L5 CPU 标为高分辨率不稳定，优先用 GPU 五档 + CPU≤L3。

**2026-11-08 CPU L3 ×3**：`r01/r02` ok（T→10）；`r03/r04` 均在 t≈8.6–9.2 以 `mech_e=nan` → `non_finite_state` 失败。属 FP32 溃坝偶发不稳，**不再重试**；该格性能以 r01/r02 中位数报，并注明仅×2。GPU L3 ×3 仍全 ok。

| `B2_L5_cpu_fp32_benchmark_r01` | **killed** EXIT:137 | — | 2026-11-05 06:30 UTC 用户中止；转 B3 |
