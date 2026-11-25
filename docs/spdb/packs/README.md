# SPDB 第一轮工作包

日期：2026-10-17。工作区 `/home/ubuntu/software/SPHinXsys-SPDB`，分支 `feature/spdb-thermo-viscoplastic`，基线 `4e686ee35`。打包时这些文件尚未 git commit。

本包是该轮全部记录、案例源码和已保存的计算结果。材料均为合成弹性，不是 950℃ 20Cr13 预测。甲方 Arrhenius 热压缩资料只归档，未进入求解器。

## 目录

| 路径 | 内容 |
|------|------|
| `docs/SPDB_Cursor_Development_Plan.md` | 开发总计划 |
| `docs/spdb/` | 交接说明、能力审计、环境、材料记录、各次结果说明、日志与 VTK |
| `tests/3d_examples/test_3d_spdb_geometry/` | 打包当时的案例源码。现行代码已迁到 `tests/extra_source_and_tests/` |

## 已完成的计算

| 工况 | 要点 |
|------|------|
| 几何，1 mm，有槽 / 无槽 | 体积相对误差约 0.6%；有槽下体比上体少 1000 粒 |
| 无槽弹性冒烟 | 短时动态松弛，不是 10 MPa 验收 |
| 无槽单轴 10 MPa，1 mm | 反力误差 3×10⁻⁸；芯部应力 0.12%；端面轴向应变误差 4.15% |
| 无槽单轴 10 MPa，0.5 mm | 端面轴向应变误差 1.89%；内部标距 1.04%；反力误差 1.8×10⁻⁷ |
| 上下接触，1 mm，同一名义力 | 接触合力互消 1.5×10⁻⁹；最大穿透 0.033 mm（粒距的 3.3%） |

名义力 12566.37 N。单轴与接触使用合成模量 2 GPa，以便在 1 mm 粒距下能看到缩短量。详细说明见 `docs/spdb/README.md` 及其中链接的结果文档。

## 复现入口

在案例构建目录的 `bin` 下：

```text
./test_3d_spdb_geometry --geometry_only
./test_3d_spdb_geometry --geometry_only --no_slot
./test_3d_spdb_geometry --elastic_no_slot
./test_3d_spdb_geometry --elastic_axial
./test_3d_spdb_geometry --elastic_axial --dp 0.0005 --end_time 0.004
./test_3d_spdb_geometry --elastic_contact --dp 0.001 --end_time 0.004
```

构建命令见 `docs/spdb/baseline_env.md`。

压缩包：`spdb_round_2026-10-17.tar.gz`（29 MB，86 个条目）。

```text
sha256  3afd14b4ba42fc34cb02b8f45da599dfd15e145a79c92a3b05611222faeddd09
```
