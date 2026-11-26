# Handoff: French stirring Phase5（本机 4090 / 隔离树）

**分支**：`feature/ophelie-fr-stir-phase5-4090`（本机搅拌后续只提交到此分支）  
**共享分支**：`feature/EM-Biot-Savart` — 请勿再往其上叠加本机搅拌专属改动

| 文件 | 说明 |
|------|------|
| [MANUAL_CN.md](MANUAL_CN.md) | **说明书**（改动范围、跑法、μ(T)、冲突规避） |
| [PULL_ON_OTHER_PC_CN.md](PULL_ON_OTHER_PC_CN.md) | **另一台如何 fetch，且不动 EM-Biot-Savart** |
| [APPLY_ON_OTHER_MACHINE_CN.md](APPLY_ON_OTHER_MACHINE_CN.md) | **到另一台如何审阅 / 合入** |
| [SUMMARY_muT_setup_CN.md](SUMMARY_muT_setup_CN.md) | 变粘设置与结果摘要（报告用） |
| [GPT_PLOT_BRIEF.md](GPT_PLOT_BRIEF.md) | 画图数据包说明 |
| [PLAN_PHASE5_snapshot.md](PLAN_PHASE5_snapshot.md) | Phase5 计划快照 |
| `seeds/` | Fig.8a T + **真 Fig.10 Q** 种子（本机数字化） |
| `scripts/` | 数字化 / 高度扫 / μ(T) 跑批脚本副本 |
| `figures/` | Fig.10 数字化 QA 图 |

大结果包仍在本机 gitignore 路径（不入库）：

`docs/ophelie/reference/doc247_longrun_four_cases/overnight_6h/`
