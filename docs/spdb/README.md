# SPDB 工作交接入口（给新 Chat 先读）

**最后更新：** 2026-10-27
**工作区：** `/home/ubuntu/software/SPHinXsys-SPDB`
**分支：** `feature/spdb-thermo-viscoplastic`
**基线 SHA：** `4e686ee357de1cb798ac47eb760dd79e2fb407ff`（与上游 `origin/master` 一致）
**状态：** 有槽保压 **60 s / 600 s / 7200 s** 均已通过（CPU）。7200 s 墙钟约 6.9 h；保压 Δ 净高约 0.07 nm，远小于弹性闭合约 3 µm。见 [`slot_hold_results.md`](slot_hold_results.md)。**尚未 git commit / push**

> 新对话请先读本文件，再按下方“必读清单”深入。不要从旧 OPHELIE / Maxwell / level-set 分支切出；不要把结果写成 950℃ 20Cr13 真实预测。

---

## 1. 项目是什么

脉冲电流辅助扩散焊（SPDB）数值验证：20Cr13 上下圆柱 + 下侧盲方槽高温保压变形，后续扩展电热耦合。
开发总计划见：[`../SPDB_Cursor_Development_Plan.md`](../SPDB_Cursor_Development_Plan.md)

| 阶段 | 内容 | 本轮 |
|------|------|------|
| P0 | 库能力审计 | ✅ |
| P1 | 参数化几何 + VTK/体积校核 | ✅ |
| + | 无孔弹性压缩 smoke | ✅（短时，非验收级） |
| + | 无槽单轴名义 10 MPa，1 mm | ✅ 力/芯部应力；端面应变 4.15%。见 [`elastic_axial_results.md`](elastic_axial_results.md) |
| + | 同上，0.5 mm | ✅ 端面应变 1.89%，内部标距 1.04%。见 [`elastic_axial_dp0.5mm_results.md`](elastic_axial_dp0.5mm_results.md) |
| + | 上下接触，1 mm，E = 2 GPa | ✅ 力平衡约 1e-8；最大穿透 33.23 µm。见 [`elastic_contact_results.md`](elastic_contact_results.md) |
| + | 同上，demo E = 100 GPa，刚度 1/2/4 | ✅ 最大穿透 0.707 / 0.363 / 0.186 µm；材料缩短几乎不变。见 [`elastic_contact_E100GPa_results.md`](elastic_contact_E100GPa_results.md) |
| 资料 | Ren 2014 Arrhenius 材料点 | ✅ 核对点相对误差 < 1e-9。单轴加载、保压、卸载已通过。见 [`material_point_results.md`](material_point_results.md) |
| + | 接触刚度 1/2/4 | ✅ 穿透下降，轴向应变变化 0.014%。见 [`elastic_contact_stiffness_results.md`](elastic_contact_stiffness_results.md) |
| + | 无孔 CK 保压，dp = 1 mm，60 s 与 600 s | ✅ 600 s（10×60 s）芯部 2.203×10⁻⁶。Δt = 30 s 同累计时间得 2.255×10⁻⁶，相对差 2.36%（略超 ≤2%）。见 [`ck_hold_results.md`](ck_hold_results.md) |
| + | 有槽弹性接触，E = 100 GPa | ✅ 穿透 0.906 µm；弹性 Δ 净高 −2.87 µm。见 [`slot_contact_results.md`](slot_contact_results.md) |
| + | 有槽保压，dp = 1 mm，60 / 600 / 7200 s | ✅ 7200 s：Δ 净高 +0.066 nm；下/上平均非弹性 1.69×10⁻⁵ / 1.42×10⁻⁵。见 [`slot_hold_results.md`](slot_hold_results.md) |
| P2+ | 温度/载荷对照、电热/LLM、SYCL | ❌ 未做 |

---

## 2. 仓库与目录关系

| 路径 | 角色 |
|------|------|
| `/home/ubuntu/software/SPHinXsys` | 原主库 worktree，`master`，**勿覆盖用户未提交改动** |
| `/home/ubuntu/software/SPHinXsys-SPDB` | 本项目 worktree（`git worktree add`），专做 SPDB |
| 共用 `.git` | 同一远程：`git@github.com:Xiangyu-Hu/SPHinXsys.git` |

创建方式回顾：

```bash
cd /home/ubuntu/software/SPHinXsys
git fetch origin
git worktree add -b feature/spdb-thermo-viscoplastic \
  /home/ubuntu/software/SPHinXsys-SPDB origin/master
```

---

## 3. 本轮已改动的文件（未提交）

```
docs/SPDB_Cursor_Development_Plan.md          # 开发总计划副本
docs/spdb/                                    # 本交接目录（审计/环境/结果）
tests/extra_source_and_tests/extra_src/shared/spdb_arrhenius.h
tests/extra_source_and_tests/extra_src/shared/spdb_arrhenius_ck.hpp
tests/extra_source_and_tests/extra_src/shared/spdb_arrhenius.cpp
tests/extra_source_and_tests/extra_src/for_3D_build/spdb_geometry.h
tests/extra_source_and_tests/extra_src/for_3D_build/spdb_ren2014_solid.h
tests/extra_source_and_tests/extra_src/for_3D_build/spdb_ren2014_solid.hpp
tests/extra_source_and_tests/extra_src/for_3D_build/spdb_ren2014_solid.cpp
tests/extra_source_and_tests/extra_src/for_3D_build/spdb_ck_axial.h
tests/extra_source_and_tests/3d_examples/test_3d_spdb_geometry/
  CMakeLists.txt
  spdb_geometry.cpp                           # 算例入口
```

后续新代码继续放在 `extra_src`（可复用源码）和 `extra_source_and_tests/3d_examples`（三维算例）。不要再往 `tests/3d_examples` 加 SPDB 文件。

构建产物（本地，通常不入库）：`build-cpu-icpx/`、`build-gnu/`（GNU 因 `-Werror`+Eigen 失败）

---

## 4. 案例怎么跑

构建目录：`build-cpu-icpx`（`icpx`，**SYCL=OFF**，CPU 参考）

```bash
cd /home/ubuntu/software/SPHinXsys-SPDB
source /opt/intel/oneapi/setvars.sh
cmake -G Ninja -B build-cpu-icpx \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=icx \
  -DCMAKE_CXX_COMPILER=icpx \
  -DSPHINXSYS_USE_SYCL=OFF \
  -DSPHINXSYS_BUILD_2D_EXAMPLES=OFF \
  -DSPHINXSYS_BUILD_OPTIMIZATION_EXAMPLES=OFF \
  -DSPHINXSYS_BUILD_PYTHON_INTERFACE=OFF \
  -DSPHINXSYS_BUILD_UNIT_TESTS=OFF \
  -DSPHINXSYS_BUILD_EXTRA_SOURCE_AND_TESTS=ON \
  -DSPHINXSYS_BUILD_MODULES=OFF \
  -DSPHINXSYS_3D=ON \
  -DSPHINXSYS_BUILD_3D_EXAMPLES=ON \
  -DSPHINXSYS_BUILD_TESTS=ON
cmake --build build-cpu-icpx --target test_3d_spdb_geometry -j$(nproc)

cd build-cpu-icpx/tests/extra_source_and_tests/3d_examples/test_3d_spdb_geometry/bin
./test_3d_spdb_geometry --geometry_only              # 下侧盲槽
./test_3d_spdb_geometry --geometry_only --no_slot    # 无槽
./test_3d_spdb_geometry --elastic_no_slot            # 合成弹性短时 smoke
./test_3d_spdb_geometry --elastic_axial              # 无槽名义 10 MPa，dp=1 mm，合成 E=2 GPa
./test_3d_spdb_geometry --elastic_axial --dp 0.0005 --end_time 0.004
./test_3d_spdb_geometry --elastic_contact --dp 0.001 --end_time 0.004
./test_3d_spdb_geometry --elastic_contact --dp 0.001 --end_time 0.004 --contact_stiffness_scale 4
./test_3d_spdb_geometry --material_point
./test_3d_spdb_geometry --ck_hold --dp 0.001 --physical_dt 60 --damping_factor 0.125
./test_3d_spdb_geometry --ck_hold --dp 0.001 --physical_dt 60 --hold_steps 10 --damping_factor 0.125
./test_3d_spdb_geometry --slot_contact --dp 0.001 --contact_stiffness_scale 1 --end_time 0.0015
./test_3d_spdb_geometry --slot_hold --dp 0.001 --physical_dt 60 --hold_steps 1 --contact_stiffness_scale 1
```

几何（SI，z 为加载轴）：

- 下圆柱：`x²+y²≤0.02²`，`z∈[-0.02,0]`；可选减去盲槽 `|x|,|y|<0.005`，`z∈[-0.01,0]`
- 上圆柱：`x²+y²≤0.02²`，`z∈[0,0.02]`
- 粗粒距 `dp=1 mm`；材料 `LinearElasticSolid`，标注 **synthetic**，不是 950℃ 20Cr13

实现要点：`GeometricShapeCylinder`（局部轴 x，用 `FromTwoVectors` 转到 z）+ `ComplexShape::subtract<GeometricShapeBox>`；骨架参考 `test_3d_muscle_soft_body_contact`。

---

## 5. 关键验证结果（摘要）

| 检查项 | 结果 |
|--------|------|
| 有槽下体 / 上体粒子数 | 24280 / 25280（差 1000 ≈ 槽体积/`dp³`） |
| 体积相对误差 | ≈0.6%（建议门槛 ≤2%，通过） |
| 槽 AABB 内下体粒子 | 0（槽向上开口） |
| 无槽上下体粒子 | 均为 25280 |
| 弹性 smoke | 跑通，约 8.2 s；**非** 10 MPa / 准静态验收 |
| 单轴 10 MPa，1 mm | 反力相对误差 3×10⁻⁸；芯部应力 0.12%；端面轴向应变 4.15% |
| 单轴 10 MPa，0.5 mm | 反力相对误差 1.8×10⁻⁷；芯部应力 0.47%；端面轴向应变 1.89%；内部标距 1.04% |
| 上下接触，1 mm，E = 2 GPa | 接触合力互消 1.5×10⁻⁹；下端面反力相对误差 2.8×10⁻⁸；最大穿透 0.033 mm（3.3% 粒距）；两体轴向应变误差约 3.66% |
| 上下接触，1 mm，E = 100 GPa | 罚函数接触，刚度 1/2/4。最大穿透 0.707 / 0.363 / 0.186 µm；材料缩短约 1.837 µm，1 倍到 4 倍只变 0.005%。传力不平衡优于 10⁻⁹ |
| 有槽接触，E = 100 GPa | 最大穿透 0.906 µm；弹性 Δ 净高 −2.87 µm；顶壁挠度 −0.60 µm。体积未宣称 |
| 有槽保压 60 s | 保压段穿透漂移 7.5×10⁻¹² m；Δ 净高 3.6×10⁻¹¹ m。下/上平均非弹性 3.46×10⁻⁷ / 2.54×10⁻⁷。反力相对误差 0.0025% |
| 有槽保压 10×60 s | 退出码 0，墙钟 3173 s。终了下/上平均非弹性 3.03×10⁻⁶ / 2.30×10⁻⁶；保压 Δ 净高 1.19 nm；穿透漂移 94 pm。这是 600 s，不是 7200 s |
| 有槽保压 120×60 s | **完成**。墙钟 24835 s。`hold_is_7200s=1`，`aborted=0`。保压 Δ 净高 +0.066 nm；穿透漂移 −0.18 nm。下/上平均非弹性 1.69×10⁻⁵ / 1.42×10⁻⁵ |
| 接触刚度 1/2/4 | 最大穿透 33.2 / 17.6 / 9.2 µm。两体材料缩短保持 183 µm，1 倍到 4 倍只变 0.014%。力平衡仍优于 10⁻⁷ |
| 无孔 CK，E = 100 GPa，60 s | 阻尼 0.125 与 0.25：端面 ε_z 相对胡克解 3.77% / 3.71%，芯部 Cauchy_zz 误差 0.21% / 0.26%。芯部非弹性应变 2.255×10⁻⁷ / 2.265×10⁻⁷，第二次松弛漂移为 0。保压时钟是 60 s，不是伪时间 |
| 无孔 CK，10×60 s | 退出码 0，墙钟 474 s。最终芯部非弹性 2.203×10⁻⁶，与该步标量回映差 0.004%。端面应变 −9.845×10⁻⁵，反力相对误差 0.32%。这是 600 s，不是 7200 s |
| 无孔 CK，20×30 s | 同一 600 s。芯部非弹性 2.255×10⁻⁶。相对 10×60 s 变化 2.36%，绝对值差 5.2×10⁻⁸。略超建议 ≤2% 门槛 |
| SYCL/GPU | 本机有 RTX 3090 Ti，**本轮未验证 SYCL** |
| GNU g++ | 配置可过，编译因 Eigen `maybe-uninitialized` + `-Werror` 失败 |

详细数字与日志：[`first_round_results.md`](first_round_results.md)、[`elastic_axial_results.md`](elastic_axial_results.md)、[`elastic_axial_dp0.5mm_results.md`](elastic_axial_dp0.5mm_results.md)、[`elastic_contact_results.md`](elastic_contact_results.md)、[`elastic_contact_stiffness_results.md`](elastic_contact_stiffness_results.md)、[`elastic_contact_E100GPa_results.md`](elastic_contact_E100GPa_results.md)、[`slot_contact_results.md`](slot_contact_results.md)、[`slot_hold_results.md`](slot_hold_results.md)、[`material_point_results.md`](material_point_results.md)、[`ck_hold_results.md`](ck_hold_results.md)、[`results/`](results/)
截面图：[`results/section_xz_with_slot.svg`](results/section_xz_with_slot.svg)、[`results/section_yz_with_slot.svg`](results/section_yz_with_slot.svg)

---

## 6. 库能力结论（P0 摘要）

详见 [`capability_audit.md`](capability_audit.md)。

| 需求 | 状态 |
|------|------|
| 圆柱/方盒/布尔减 | 可用 |
| Lattice 固体粒子 | 可用 |
| 弹性应力更新 + 动态接触 + 阻尼 | 可用（准静态接触**未验证**） |
| 高温蠕变本构 | **库中没有**（`ViscousPlasticSolid` ≠ creep） |
| E(T) / 热力耦合 | **没有** |
| 脉冲电热 | **没有**（后续阶段） |
| Ren 2014 Arrhenius | 材料点已实现。10 MPa 保压是低于 0.01 s⁻¹ 的外推，单点应变约 3×10⁻⁵ |
| CK 塑性 | 仅 `HardeningPlasticSolid` 有较完整证据 |

---

## 7. 硬约束（新 Chat 必须遵守）

1. 不擅自填入 20Cr13 高温蠕变常数；synthetic/proxy 必须标注，不能称为实际工况预测。
2. 第一阶段不做冶金结合、等离子体、扩散焊强度。
3. 导电焦耳热 ≠ 感应加热 A–φ；不要误搬电磁焊接感应模型。
4. 默认不 push；未确认写权限不要推远程。
5. 不清理/覆盖 `/home/ubuntu/software/SPHinXsys` 上用户未提交工作。
6. 每阶段验证通过再进入下一阶段；失败保留诊断。

---

## 8. 建议下一轮任务（按优先级）

1. 有槽 7200 s 已通过（墙钟约 6.9 h）。下一步优先：重启快照（含物理时钟）或 900–1000℃ / 20 MPa 对照。
2. 孔高必须把弹性闭合（约 −3.16 µm）与保压增量（7200 s 约 +0.07 nm）分开。`RepulsionForceCK` 不能用于齐平保压。
3. 无孔步长敏感性：30 s 相对 60 s 差 2.36%，略超 ≤2%。
4. 论文图 6 读点、SYCL 对照仍缺。监控：`docs/spdb/scripts/check_slot_hold_status.sh`。
5. E(950℃) 仍是 demo_assumption，不是 Ren 2014。本轮进度：[`round_progress_2026-10-21.md`](round_progress_2026-10-21.md)。

---

## 9. 必读清单（新 Chat 阅读顺序）

1. **本文件** `docs/spdb/README.md`（你在这里）
2. [`../SPDB_Cursor_Development_Plan.md`](../SPDB_Cursor_Development_Plan.md) — 总计划与验收边界
2b. [`SPDB_Next_Round_Cursor_Plan_2026-10-17.md`](SPDB_Next_Round_Cursor_Plan_2026-10-17.md) — 文献本构下一轮计划
3. [`capability_audit.md`](capability_audit.md) — P0 真实类名/缺口
4. [`baseline_env.md`](baseline_env.md) — 构建环境
5. [`first_round_results.md`](first_round_results.md) — 首轮数字结果
6. [`material_20cr13_arrhenius.md`](material_20cr13_arrhenius.md) — Ren 2014 完整系数
6b. [`material_point_results.md`](material_point_results.md) — 材料点核对结果
7. [`elastic_axial_results.md`](elastic_axial_results.md) — 无槽 10 MPa，1 mm
8. [`elastic_axial_dp0.5mm_results.md`](elastic_axial_dp0.5mm_results.md) — 同一工况，0.5 mm
9. [`elastic_contact_results.md`](elastic_contact_results.md) — 1 mm 上下接触
10. [`elastic_contact_stiffness_results.md`](elastic_contact_stiffness_results.md) — 接触刚度 1/2/4
11. 源码：`tests/extra_source_and_tests/extra_src/shared/spdb_arrhenius.{h,cpp}`、`extra_src/for_3D_build/spdb_geometry.h`、算例 `3d_examples/test_3d_spdb_geometry/spdb_geometry.cpp`

---

## 10. 给新 Chat 的启动提示（可复制）

```
请阅读 docs/spdb/README.md 及其中链接的计划/审计/结果文档。
工作区已是 SPHinXsys-SPDB，分支 feature/spdb-thermo-viscoplastic。
P0/P1、弹性 smoke、无槽 10 MPa 单轴（1 mm 与 0.5 mm）和 1 mm 上下接触已完成且未提交。
0.5 mm 端面应变误差 1.89%，内部标距 1.04%。接触最大穿透是粒距的 3.3%，力平衡到约 1e-8。
Ren 2014 材料点已通过。无孔 CK 保压已做到 60 s 和 600 s（10×60 s，芯部非弹性 2.203×10⁻⁶）。E = 100 GPa 齐平接触要用罚函数实测，不要用 RepulsionForceCK。不要把 600 s 或 10 MPa、7200 s 的外推写成孔塌陷预测。
```
