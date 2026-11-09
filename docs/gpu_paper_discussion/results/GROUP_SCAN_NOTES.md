# B2 + B4 L3 group 扫描

最后更新：2026-11-08（**24/24 ok**）

## 目的
固定 L3，只改 SYCL `work_group_size` ∈ {32,64,128,256}×3，检验主表默认 G=64 是否合理。  
**B3 不做**（细档不稳定，文章可能用 B4 替代）。

## 结果摘要（median wall_solver）

| G | B2 L3 | B4 L3 | vs G=64 |
|---|---|---|---|
| 32 | 10.09 s | 4.09 s | 慢 ~5–7% |
| **64** | **9.58 s** | **3.81 s** | 基线 |
| 128 | 9.16 s | 3.79 s | B2 快 ~4%；B4 ≈持平 |
| 256 | 8.50 s | 3.66 s | 最快（B2 −11%；B4 −4%） |

重复性好（spread 多 ≤2%）。主表仍可用 **G=64**；文中可注 G=256 略优但非默认。  
（早期误跑的 `g1007` 仅作事故记录，不进主表。）

## 实现（不改 `src/`）
| 文件 | 作用 |
|---|---|
| `extra_src/shared/gpu_paper/gpu_paper_sycl_group.h` | 读 `GPU_PAPER_GROUP`，`setWorkGroupSize` |
| `run_b2.sh` / `run_b4.sh` | 可选 env；RUN_ID 含 `_g${GROUP}` |
| `run_group_scan_b2_b4.sh` | 批量扫描 |

## 跑法

```bash
# 先重编 GPU B2/B4（含 register 安全上限）
cd /home/yongchuan/SPHinXsys-gpu-paper/build && ninja test_3d_gpu_paper_b2_dambreak test_3d_gpu_paper_b4_passive_cantilever

# 清掉错误残留（若曾 export 过 GPU_PAPER_GROUP=1007）
unset GPU_PAPER_GROUP GROUPS

tmux new -s groupscan
bash /home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/scripts/run_group_scan_b2_b4.sh
```

单点示例：

```bash
env -u GPU_PAPER_GROUP CUDA_VISIBLE_DEVICES=1 GPU_PAPER_GROUP=128 RUN_SUFFIX=r01 \
  bash .../scripts/run_b2.sh L3 gpu_fp32 benchmark
```

G=64 的点可与已有 `B*_L3_gpu_fp32_benchmark_r0*`（未带 `_g64`）对照；扫描跑会写成 `*_benchmark_g64_r0*`。

## 事故记录（2026-11-08）
父 shell 里残留 **`GPU_PAPER_GROUP=1007`**（或旧脚本误用 `GROUPS`）→ B2 g1007 三次“碰巧”ok，B4 L3 因寄存器（80×1007 > 65536）`EXIT:134`。  
已改：扫描用 `SCAN_GROUPS`；`gpu_paper_sycl_group.h` 硬顶 **512**。
