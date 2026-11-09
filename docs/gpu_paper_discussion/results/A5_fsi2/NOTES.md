# A5 圆柱+柔性尾板 FSI2 — 运行记录

最后更新：2026-11-07（**L1–L3 share 三端 validate 全齐**）

## 约定（Duration V2）

| 项 | 值 |
|---|---|
| 终点 | **End_Time=200**；smoke 用短 T |
| 分辨率 | D/dp：L1=10，L2=20，L3=40 |
| API | **SPHSolver + FSI CK** |
| Insert 粒子 | **relax → `./reload/Reload.xml` → Reload**（对齐经典 FSI2） |
| 边界 | **emitter + disposer**（替代 periodic+sponge） |
| 固体阻尼 | **不加** `StructureNumericalDamping` |
| VTP | validate + GPU L2 |
| GPU | `CUDA_VISIBLE_DEVICES=1` |

## 工作流

```bash
# 1) 松弛（每分辨率 / 每 build 目录各一次）
bash .../scripts/run_a5.sh L1 gpu_fp32 relax

# 2) 短测 / 正式（缺 Reload.xml 时脚本会自动先 relax）
bash .../scripts/run_a5.sh L1 gpu_fp32 smoke
```

环境变量：`GPU_PAPER_RELAX=1`、`GPU_PAPER_RELOAD=1`（默认 smoke/validate 为 1）；也可 `--relax` / `--reload`。

换 `D/dp` 或换 CPU/GPU 构建目录后需**重新 relax**（reload 与分辨率绑定）。

## 重要差异 vs 经典 `test_2d_fsi2`
- 通道周期 → emitter/disposer
- 密度：`CompressionSummation` + `DensityRegularization`
- 正式幅频须在 **175–200** 窗验证后再报

## 状态

| run_id | EXIT | 备注 |
|---|---:|---|
| A5_L1_gpu_fp32_relax_r01 | 0 | Reload.xml OK |
| A5_L1_gpu_fp32_smoke_r01（首跑） | 1 | Indicator 空指针；已修 |
| A5_L1_gpu_fp32_smoke_r01（修复后，≈03:47） | **0** | T→2.002；wall≈0.31 s；tip 有微动；N_water≈5.1k→5.4k；`insert_reload=1` |
| A5_L1_gpu_fp32_smoke_r02 | **0** | 同上可复现；wall≈0.31 s；`status=ok` |
| build cpu_fp32 / cpu_fp64 | **0** | 重编通过 |
| A5_L1_*_relax_r01（gpu/c32/c64） | **0** | 三端 Reload.xml 均有 |
| A5_L1_gpu_fp32_smoke_r03 | **0** | t=2.002；wall≈0.31 s |
| A5_L1_cpu_fp32_smoke_r01 | **0** | t=2.002；wall≈0.22 s |
| A5_L1_cpu_fp64_smoke_r01 | **0** | t=2.002；wall≈0.21 s |
| A5_L1_gpu_fp32_validate_ck01 | **0** | T→200.002；wall≈42.1 s；`status=ok`；有 tip=0 毛刺 |
| A5_L1_cpu_fp32_validate_ck01 | **0** | T→200.003；wall≈38.8 s；`status=ok`；有 tip=0 毛刺 |
| A5_L1_cpu_fp64_validate_ck01 | **0** | T→200.005；wall≈41.1 s；`status=ok`；有 tip=0 毛刺 |
| A5_L1_gpu_fp32_validate_ck02 | **0** | T→200.006；wall≈41.7 s；`tip_reject=0` |
| A5_L1_cpu_fp32_validate_ck02 | **0** | T→200.004；wall≈38.9 s；`tip_reject=0` |
| A5_L1_cpu_fp64_validate_ck02 | **0** | T→200.006；wall≈40.7 s；`tip_reject=0` |

## Tip 探针修复（2026-11-06）

根因：采样时重建 Lagrangian observer 接触 + 当前 Position CLL → 空邻域 → `(0,0)`。已改为不重建邻域 + `getObservedQuantity()` + 拒收近零。

**ck02 验证：三端 `tip_reject=0`，`tip_y==0` 次数=0，毛刺已消除。**

## L1 物理一致性（ck02，T=200）

**探针干净，但 GPU/CPU 轨迹仍非逐点一致。**

| 项 | 观察 |
|---|---|
| 通跑 | 三端 `EXIT:0`，`N_insert=596`，`N_water≈5.03k–5.06k` |
| 175–200 半峰幅值 | gpu≈0.214 / c32≈0.207 / c64≈0.119 |
| 175–200 相关 | gpu–c32 ≈ **−0.94**（近反相）；gpu–c64 ≈ −0.33 |
| 初值 | 三端各自 `Reload.xml`，L1 涡激振动易漂相位 |

下一步：共用同一 Reload 再对照，或升 L2；幅值量级 GPU≈CPU32 可作粗一致性参考，不能报逐点等价。

## 对齐流程（shared Reload）

三端各自 relax 的 `Reload.xml` sha256 不同 → 相位易漂。对齐时：

1. GPU L1 relax 并发布到 `shared_reload/L1/Reload.xml`
2. 三端 validate 均 `GPU_PAPER_SHARED_RELOAD=1`（脚本会拷入各自 `bin/reload/`）
3. 优先看 **gpu_fp32 ↔ cpu_fp32**（同精度）；cpu_fp64 作参考

```bash
env RUN_SUFFIX=share01 GPU_PAPER_PUBLISH_SHARED_RELOAD=1 \
  bash .../run_a5.sh L1 gpu_fp32 relax
env RUN_SUFFIX=share01 GPU_PAPER_SHARED_RELOAD=1 \
  bash .../run_a5.sh L1 gpu_fp32 validate
```

## L1 对齐结果（share01，同一 Reload sha256=`066e6d…`）

| run_id | EXIT | wall | tip_reject |
|---|---:|---:|---:|
| A5_L1_gpu_fp32_validate_share01 | 0 | 41.7 | 0 |
| A5_L1_cpu_fp32_validate_share01 | 0 | 38.7 | 0 |
| A5_L1_cpu_fp64_validate_share01 | 0 | 41.6 | 0 |

相对 ck02（各自 Reload）：

| 指标（gpu↔c32） | ck02 | share01 |
|---|---:|---:|
| t=0 tip_disp_y | 不一致 | **完全一致** `1.3804e-4` |
| 相关 [0,20] | 0.21 | **0.87** |
| 相关 [175,200] | **−0.94**（反相） | **+0.66** |
| late rms \|Δ\| | 0.257 | **0.093** |

**结论：** 共用 Reload 后 GPU↔CPU32 明显更齐（早段同步、晚窗同相）；仍非逐点等价（中段仍漂，L1+FP 混沌）。CPU64 幅值偏大，作参考即可。下一步可试 L2 shared，或短时窗 T≤20 作严格对齐验收。

## 对齐加强（align）

`share01` 显示 **t≈3.5 tip 仍齐，但 N_water 已差 ~10**（emitter/disposer 对 FP 敏感），随后 tip 分叉。

代码增加：
- `GPU_PAPER_DISABLE_PARTICLE_SORT=1`
- `GPU_PAPER_FIXED_ACOUSTIC_DT=0.003`（两端同一时间步）
- 脚本 `GPU_PAPER_ALIGN=1` → 自动 shared Reload + 上述两项；CPU 默认 `OMP_NUM_THREADS=1`

重编后跑 `RUN_SUFFIX=align01` + `GPU_PAPER_ALIGN=1` validate（沿用已有 `shared_reload/L1`，一般无需再 relax）。

### align01 结果（相对 share01）

三端 `EXIT:0`，`tip_reject=0`，同一 Reload；配置 `DISABLE_SORT=1`、`FIXED_ACOUSTIC_DT=0.003`。

| 指标（gpu↔c32） | share01 | align01 |
|---|---:|---:|
| 相关 [0,20] | **0.87** | −0.05 |
| 相关 [175,200] | **+0.66** | −0.51 |
| late rms \|Δ\| | **0.093** | 0.74 |
| late 半峰幅值 gpu/c32/c64 | 0.18/0.13/0.28 | **0.63/0.64/0.61**（幅值更齐，但振幅变大） |

**结论：** 固定 dt=0.003 + 关 sort **没有改善轨迹对齐**，反而早段相关崩溃；幅值三端更接近但振动更强，疑固定步长偏大。当前仍以 **share01（仅 shared Reload）为最佳对齐基线**；align 旋钮默认不要用于正式叙事。

### align02（dt=0.001）

三端 `EXIT:0`；墙时≈82/101/110 s。相对 share01：

| 指标（gpu↔c32） | share01 | align02 |
|---|---:|---:|
| 相关 [0,20] | **0.87** | −0.78 |
| 相关 [175,200] | **+0.66** | −0.47 |
| late 半峰 | 0.18/0.13/0.28 | 0.61/0.62/0.63（幅值齐但偏大） |
| 早段 ΔN_water | ~10–20 | **更小**（~1–7） |

N_water 更齐了，但 tip 轨迹仍不同步。**固定 dt 路线停用**；对齐基线保持 share01。

## 本机收束：L2 / L3（进行中）

正式叙事用 **shared Reload**（与 L1 `share01` 同策略）；跳过 align 固定 dt。顺序：GPU relax+publish → GPU smoke → 三端 validate(share)。

```bash
A5=/home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/scripts/run_a5.sh
# —— L2 ——
env RUN_SUFFIX=share01 GPU_PAPER_PUBLISH_SHARED_RELOAD=1 bash $A5 L2 gpu_fp32 relax
env RUN_SUFFIX=r01 bash $A5 L2 gpu_fp32 smoke
env RUN_SUFFIX=share01 GPU_PAPER_SHARED_RELOAD=1 bash $A5 L2 gpu_fp32 validate
env RUN_SUFFIX=share01 GPU_PAPER_SHARED_RELOAD=1 bash $A5 L2 cpu_fp32 validate
env RUN_SUFFIX=share01 GPU_PAPER_SHARED_RELOAD=1 bash $A5 L2 cpu_fp64 validate
# —— L3（L2 齐后再开）——
env RUN_SUFFIX=share01 GPU_PAPER_PUBLISH_SHARED_RELOAD=1 bash $A5 L3 gpu_fp32 relax
env RUN_SUFFIX=r01 bash $A5 L3 gpu_fp32 smoke
env RUN_SUFFIX=share01 GPU_PAPER_SHARED_RELOAD=1 bash $A5 L3 gpu_fp32 validate
env RUN_SUFFIX=share01 GPU_PAPER_SHARED_RELOAD=1 bash $A5 L3 cpu_fp32 validate
env RUN_SUFFIX=share01 GPU_PAPER_SHARED_RELOAD=1 bash $A5 L3 cpu_fp64 validate
```

| 档 | 状态 |
|---|---|
| L1 share01 | **ok** |
| L2 share01 | **ok** |
| L3 share01 | **ok**（见下） |

### L2 结果（share01，D/dp=20）

| run_id | EXIT | wall_solver | 备注 |
|---|---:|---:|---|
| A5_L2_gpu_fp32_relax_share01 | 0 | — | published `shared_reload/L2` sha256=`713131cb…` |
| A5_L2_gpu_fp32_smoke_r01 | 0 | 0.67 s | T→2.00 |
| A5_L2_gpu_fp32_validate_share01 | 0 | 85.5 s | T→200；tip_reject=0 |
| A5_L2_cpu_fp32_validate_share01 | 0 | 113 s | T→200；tip_reject=0 |
| A5_L2_cpu_fp64_validate_share01 | 0 | 127 s | T→200；tip_reject=0 |

gpu↔c32：`t=0` tip 接近（≈2.46e-4）；晚窗半峰幅值接近（≈0.81/0.80），但相关 [175,200]≈**−0.79**（相位仍漂，同 L1 量级问题）。正式通跑与幅值可比；不报逐点等价。

### L3 结果（share01，D/dp=40）

| run_id | EXIT | wall_solver | 备注 |
|---|---:|---:|---|
| A5_L3_gpu_fp32_relax_share01 | 0 | — | published `shared_reload/L3` sha256=`8f127b5c…` |
| A5_L3_gpu_fp32_smoke_r01 | 0 | 1.41 s | T→2.00 |
| A5_L3_gpu_fp32_validate_share01 | 0 | 191 s | T→200；tip_reject=0 |
| A5_L3_cpu_fp32_validate_share01 | 0 | 421 s | T→200；tip_reject=0 |
| A5_L3_cpu_fp64_validate_share01 | 0 | 520 s | T→200；tip_reject=0 |

gpu↔c32：`t=0` tip 一致（≈−9.3e-5）；晚窗半峰 ≈0.90/0.91，相关 [175,200]≈**+0.76**（明显好于 L2）。GPU vs CPU32 ≈ **2.2×**。

| 档 | GPU | CPU32 | CPU64 | GPU/CPU32 |
|---|---|---|---|---|
| L1 | 41.7 s | 38.7 s | 41.5 s | 0.93× |
| L2 | 85.5 s | 113 s | 127 s | 1.32× |
| L3 | 191 s | 421 s | 520 s | **2.20×** |
