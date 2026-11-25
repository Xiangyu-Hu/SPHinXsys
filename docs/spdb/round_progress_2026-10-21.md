# SPDB 本轮进度（2026-10-21）

工作区：`/home/ubuntu/software/SPHinXsys-SPDB`，分支 `feature/spdb-thermo-viscoplastic`。
计划：[`SPDB_Next_Round_Cursor_Plan_2026-10-17.md`](SPDB_Next_Round_Cursor_Plan_2026-10-17.md)。
新代码只在 `tests/extra_source_and_tests/`，未改 `src/`。默认未 git commit / push。

## 已完成

| 项 | 结果摘要 | 记录 |
|----|----------|------|
| Ren 2014 材料点 | 核对点相对误差 < 1e-9；10 MPa、7200 s 常速参考 3.09×10⁻⁵ | `material_point_results.md` |
| 无孔 CK 60 s / 600 s | 600 s 芯部非弹性 2.203×10⁻⁶；与标量回映对齐 | `ck_hold_results.md` |
| 物理步长减半 | 20×30 s 相对 10×60 s 差 **2.36%**（略超建议 ≤2%） | `ck_hold_results.md` |
| E = 100 GPa 齐平接触 1/2/4 | 穿透 0.707 / 0.363 / 0.186 µm | `elastic_contact_E100GPa_results.md` |
| 有槽弹性接触 | 穿透 0.906 µm；弹性 Δ 净高 −2.87 µm | `slot_contact_results.md` |
| 有槽保压 60 s | 穿透漂移 7.5 pm；保压 Δ 净高 36 pm | `slot_hold_results.md` |
| 有槽 10×60 s | 退出码 0，墙钟 3173 s。保压 Δ 净高 +1.19 nm；穿透漂移 +94 pm。下/上平均非弹性 3.03×10⁻⁶ / 2.30×10⁻⁶ | `slot_hold_results.md`、`results/slot_hold_dt60_n10/` |
| 有槽 120×60 s | **完成**。墙钟 24835 s（约 6.9 h）。`hold_physical_s=7200`，`aborted=0`。保压 Δ 净高 +0.066 nm；穿透漂移 −0.18 nm。下/上平均非弹性 1.69×10⁻⁵ / 1.42×10⁻⁵ | `slot_hold_results.md`、`results/slot_hold_dt60_n120/` |

材料标注：E = 100 GPa 等为 **demo_assumption**，不是 Ren 2014 模量，也不是 950℃ 工艺预测。

## 关键实现选择

- 无孔保压：CK `StructureIntegration1stHalf`。
- 齐平/有槽接触：旧 API `ContactForce`。库内 `RepulsionForceCK` 在齐平面上不可用。
- 有槽保压：`ContactForce` + `ContactScaledRen2014` + `commit_ren2014_physical_increment`；松弛中非弹性冻结。
- 空腔指标：粒子中心净高/跨度/顶壁挠度；**不宣称体积**。

## 成本（实测）

- 有槽 60 s：约 9.5 min
- 有槽 600 s（10×60 s）：约 53 min
- 有槽 7200 s（120×60 s）：**约 6.9 h**（早先按 600 s 外推约 10.5 h，偏保守）

## 未完成 / 未验证

- 重启快照（需保存物理保压时钟，不只 VTK 粒子）
- 论文图 6 读点对照
- SYCL / GPU 对照（CPU 完成 ≠ GPU 完成）
- 有槽 900℃ / 1000℃ / 20 MPa 对照
- git commit（需用户明确要求）

## 可复现入口

```bash
# 有槽 60 s smoke
docs/spdb/scripts/run_slot_hold_smoke.sh

# 有槽 600 s 预演
HOLD_STEPS=10 PHYSICAL_DT=60 docs/spdb/scripts/run_slot_hold_smoke.sh

# 有槽 7200 s（后台 + 状态文件）
docs/spdb/scripts/run_slot_hold_background.sh
docs/spdb/scripts/check_slot_hold_status.sh
```

构建说明见 [`baseline_env.md`](baseline_env.md)。交接入口：[`README.md`](README.md)。
