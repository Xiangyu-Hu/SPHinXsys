# B4 被动 NeoHookean 悬臂梁

最后更新：2026-11-07（CK/SYCL case 初建；物理阻尼暂缺）

## 来源与差异

| 项 | 值 |
|---|---|
| Classic | `tests/3d_examples/test_3d_passive_cantilever_neohookean` |
| GPU paper | `test_3d_gpu_paper_b4_passive_cantilever`（SPHSolver CK） |
| 材料 | NeoHookean：ρ=1265，E=5e4，ν=0.45 |
| 驱动 | 重力 g=9.8（`GravityForceCK`） |
| **阻尼** | Classic 有 `DampingProjectionInner(μ=50)`；**CK 无对应实现 → 当前省略**（meta `physical_damping=omitted_no_CK`） |
| 几何 | PL=0.1，PH=PW=0.04，SL=0.02；自由端探针 (PL,PH,PW) |

## 约定（Duration V2，仿 B3）

| 项 | 值 |
|---|---|
| 终点 | **End_Time=1.0 s**（smoke 默认 0.1） |
| 分辨率 | `PH/dp`：L1=**6**（源码），L2=12，L3=24，L4=48，L5=96 |
| 三端 | `gpu_fp32` / `cpu_fp32` / `cpu_fp64` |
| 二进制 | `test_3d_gpu_paper_b4_passive_cantilever` |
| GPU | `CUDA_VISIBLE_DEVICES=1` |
| 探针 | tip 位移/速度、`deflection_y`、动能 |

## 正式矩阵

|  | GPU FP32 | CPU FP32 | CPU FP64 |
|---|---|---|---|
| L1 smoke | **ok** | **ok** | **ok** |
| L1 bench | **ok** | **ok** | **ok** |
| L2 bench | **ok** | **ok** | **ok** |
| L3 bench | **ok** | **ok** | **ok** |
| L4 bench | **ok** | **ok** | **ok** |
| L5 bench | **ok** | **ok** | **ok** |

L1 smoke 三端 tip `deflection_y@0.1` 一致（≈-5.32e-2）。无物理阻尼时 T=1 位移会振荡，benchmark@1.0 的 defl≈-2.06e-2 属动态过程，勿当静挠度。

**L2**（N_lattice=5184，PH/dp=12，T=1）：三端 `status=ok`；tip `defl_y@1.0` ≈ −6.02e-2（CPU32/64/GPU 相对差 <0.3%）。`wall_solver`：CPU32=2.31 s，CPU64=2.84 s，GPU=1.33 s（GPU vs CPU32 ≈ **1.74×**；L1 时尚 GPU 慢于 CPU）。

**L3**（N_lattice=41472，PH/dp=24，T=1）：三端 `status=ok`；tip `defl_y@1.0` ≈ −5.30e-2（三端几乎重合）。`wall_solver`：CPU32=23.0 s，CPU64=31.2 s，GPU=3.77 s（GPU vs CPU32 ≈ **6.09×**）。

| 档 | N | CPU32 | CPU64 | GPU | GPU/CPU32 |
|---|---|---|---|---|---|
| L1 | 648 | 0.48 s | 0.52 s | 0.60 s | 0.81× |
| L2 | 5184 | 2.31 s | 2.84 s | 1.33 s | 1.74× |
| L3 | 41472 | 23.0 s | 31.2 s | 3.77 s | **6.09×** |
| L4 | 331776 | 395 s | 561 s | 62.0 s | **6.37×** |
| L5 | 2654208 | 7145 s | 9803 s | 1069 s | **6.69×** |

**L4**（N_lattice=331776，PH/dp=48，T=1）：三端 `status=ok`；tip `defl_y@1.0` ≈ −3.97e-2（三端相对差 <0.03%）。`wall_solver`：CPU32=395 s，CPU64=561 s，GPU=62.0 s。

**L5**（N_lattice=2654208，PH/dp=96，T=1）：三端 `status=ok`；tip `defl_y@1.0` ≈ −3.22e-2（CPU32/64 相对差 ≈0.02%）。`wall_solver`：CPU32=7145 s（≈2.0 h），CPU64=9803 s（≈2.7 h），GPU=1069 s（≈18 min；vs CPU32 ≈ **6.69×**）。

## 工作流

```bash
# GPU
cd /home/yongchuan/SPHinXsys-gpu-paper/build && cmake .. >/dev/null && ninja test_3d_gpu_paper_b4_passive_cantilever
# CPU
cd /home/yongchuan/SPHinXsys-gpu-paper/build_cpu_fp32 && cmake .. >/dev/null && ninja test_3d_gpu_paper_b4_passive_cantilever
cd /home/yongchuan/SPHinXsys-gpu-paper/build_cpu_fp64 && cmake .. >/dev/null && ninja test_3d_gpu_paper_b4_passive_cantilever

env RUN_SUFFIX=r01 bash /home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/scripts/run_b4.sh L1 cpu_fp32 smoke
```

## 状态

| 项 | 备注 |
|---|---|
| case 接线 | 已建 |
| 编译/短测 | **L1–L5 三端 bench 全齐** |
| 相对 classic | 缺物理阻尼；定量对比时须声明 |
