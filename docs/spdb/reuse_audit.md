# SPDB 程序温度 Demo 复用审计

日期：2026-11-16
工作区：`/home/ubuntu/software/SPHinXsys-SPDB`
分支：`feature/spdb-thermo-viscoplastic`（相对 `origin/master` behind，本地基线 SHA 与先前一致）
构建：`build-cpu-icpx`（icpx，`SPHINXSYS_USE_SYCL=OFF`）可用。

## 已有且本轮复用

| 项 | 文件 / 类 | 状态 |
|----|-----------|------|
| Ren2014 材料点 | `extra_src/shared/spdb_arrhenius.*` | 已验证；包内 `material_reference.py` 通过 |
| 三维本构 + 物理提交 | `Ren2014X20Cr13Solid`，`commit_ren2014_physical_increment` | 有槽 7200 s 已跑通；本轮加 `rate_scale` 门控 |
| 接触保压 | `ContactForce` + `ContactScaledRen2014`，`spdb_slot_hold.cpp` | 齐平接触可用；不用 `RepulsionForceCK` |
| BC 模式 | 下端轴向支撑、上端恒力、径向自由、无摩擦 | 复用；载荷改为 5 MPa |
| 解析双槽诊断 | `SlotMode::Both` | 仅对照；**主工况用甲方 STL** |
| STL 粒子范式 | 库内 `TriangleMeshShapeSTL` + LevelSet + Lattice | 库有，旧 SPDB case 未接 |

## 本轮必须新建

| 项 | 说明 |
|----|------|
| `test_3d_spdb_programmed_temperature` | 独立 case，不覆盖旧 geometry/slot_hold |
| STL 双槽粒子 | `lower_m.stl` / `upper_m.stl`，scale=1，单位 m |
| 程序温度 T(t) | CSV 分段线性；全粒子均匀赋温 |
| 热膨胀 | demo α=12e-6；参考 T0=293.15 K |
| 850–900℃ 门函数 | smoothstep 乘到非弹性速率 |
| 双槽贯通空腔指标 | 名义净高 20 mm（非旧盲槽 11 mm） |
| F0–F3 对照驱动 | 宏步 10–30 s，落点 300/600/…/1800 |

## 明确不做 / 缺项

- 不改 `src/`；不求解传热、电场、石墨、焊合、相变
- 本轮以 CPU 参考为准；SYCL/GPU 标缺项
- `E/ν/ρ/α` 为 `demo_assumption`，不是 Ren2014 或实测 20Cr13 全温物性
- 旧有槽 7200 s 结果保留可回归，不改写时间标签

## 输入包

`docs/spdb/program_temperature_demo_2026-11-16/`（自 zip 解压）。算例 `input/` 另有运行时拷贝。
