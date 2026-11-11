# A/B/C VTP 存放说明（本机 / SSH，不进 Git）

计算结果（含 VTP）按政策 **不推送 GitHub**，只留在本机：

`docs/gpu_paper_discussion/results/<CASE>/runs/<run_id>/`

## 当前本机盘点（2026-11-10）

| 系列 | VTP 情况 | 典型路径 |
|---|---|---|
| A1 / A2 / A3(D4) | **有少量** validate/showcase 归档 | `.../runs/*/vtp/*.vtp` |
| A4 / A5 | **results 里基本无 VTP**；部分残留在 `build_*/.../bin/output*` | 未统一归档进 results |
| B1 | **有**（bunny particles，体积很大 ~13G runs） | `results/B1_stanford_bunny/runs/*/particles/*.vtp` |
| B2 / B3 / B4 | **正式 benchmark 默认关 VTP**（`run_b*.sh` 里 `VTP=0`；GPU showcase 才开） | `runs/*/snapshots/` 多为空 |
| C1 | **无 VTP 文件**（仅有 `selected_snapshot_times.txt`） | GPU×3 未写出/未归档到 runs |

打包下载：`archives/GPU_PAPER_ABC_{SLIM,FULL}_2026-11-10.tgz` 含上述 `results/` 里已有的 VTP（FULL 含 B1 bunny）。

外卡复现 B 系列性能：**不依赖**这些 VTP；用 case 源码 + `run_b*.sh` 即可。
