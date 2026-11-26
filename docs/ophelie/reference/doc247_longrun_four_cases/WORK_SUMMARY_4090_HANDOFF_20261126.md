# Doc247/248 工作总结 — 4090 交接（2026-11-26）

分支：`feature/EM-Biot-Savart`  
身份：展示/工程短窗与长跑脚手架；**不宣称论文 PASS**。

## 1. 本轮在做什么

| 线 | 目标 | 现状 |
|----|------|------|
| 法国 NC | Fig.8 T + Fig.6 Q 展示；EM→Fig6 blend | prescribed 与 α=0.8/60s 短窗可用；长时 450s 类挂机已放弃优先 |
| 法国搅拌 | Fig.8a 主；Fig.10 核；**Fig.11 等温 Np–Re** | 10 s 耦合软点有；Fig.11 需 **3 个常物性点**（5/10/20 rpm） |
| 甲方 NC | Doc248 C0 | 脚本包已有，他机继续 |
| 甲方搅拌 | 原 C2 **弹性**太慢；新增 **C2r 刚体桨** 冲 ~1400 K | 源码+脚本已加（见下） |

## 2. Git 追踪口径（瘦身）

**跟进仓库**

- 源码 / CMake / Doc247–248 **脚本与锁定文档**
- `standard_data/doc247_seeds/`（Fig8 T、Fig6 Q shape、STIR Fig8a）
- `doc248_client_cases/configs/*.json` + run/COMPILE 脚本
- 本交接 `WORK_SUMMARY_*.md`

**不跟进（已在 .gitignore）**

- `doc247_longrun_four_cases/overnight_6h/`（大体积 CSV/VTP/log）
- `doc247_french_digitization_kit/`（完整数字化栅格；keeper 在 `standard_data`）
- `doc248_client_cases/runs/`
- `**/out/`、`**/trz_snaps/`、`**/run.log`

本地对比图包可留在机器上：`overnight_6h/acceptance_pack_short/`（ignored）。

## 3. 4090 上建议并行的算例

```bash
cd /path/to/SPHinXsys
git checkout feature/EM-Biot-Savart && git pull

# 编译（按机改 BUILD）
bash docs/ophelie/reference/doc248_client_cases/COMPILE.sh
cmake --build "$BUILD" --target \
  test_3d_ophelie_french_stirring_em \
  test_3d_ophelie_client_stirring_rigid_heated \
  -j"$(nproc)"
```

**会话 A — 法国 Fig.11（3 case 一条脚本）**

```bash
tmux new -s fig11
cd /path/to/SPHinXsys
bash docs/ophelie/reference/doc247_longrun_four_cases/02_fr_stir/run_fr_stir_fig11_isothermal_points.sh
# Ctrl-b d
```

**会话 B — 甲方刚体搅拌冲 1400 K**

```bash
tmux new -s c2r
cd /path/to/SPHinXsys
bash docs/ophelie/reference/doc248_client_cases/01_generate_particles.sh   # 若尚无 Reload
END_TIME=7200 MAX_WALL_HOURS=48 TARGET_T_K=1400 HEAT_POWER=15000 \
  bash docs/ophelie/reference/doc248_client_cases/04_run_C2r_rigid_stir.sh
```

可选再开：Doc248 C0 NC（`02_run_C0_nc.sh`）；法国 NC EM-blend 仅作短窗展示，不必再挂 450 s 搅拌长跑。

## 4. 新 case：甲方刚体桨（C2r）

| | |
|--|--|
| 源码 | `tests/.../test_3d_ophelie_client_stirring_rigid_heated/` |
| 二进制 | `test_3d_ophelie_client_stirring_rigid_heated` |
| 跑法 | `docs/ophelie/reference/doc248_client_cases/04_run_C2r_rigid_stir.sh` |
| 配置 | `configs/C2r_rigid_effective.json` |
| 力学 | **整桨规定刚体转动**（非 Simbody 法国桨；非 NeoHookean 弹性） |
| 热 | 规定均匀体积热（默认 15 kW ASSUMED）；**非 EM** |
| 停机 | `--target-T-K=1400`（体积平均 T 达标可早停）；或 `--end-time` / `--max-wall-hours` |

与 C2 弹性对比：弹性固体子步墙钟不可接受时，用 C2r 做升温长跑。

## 5. 验收图优先级（法国）

1. NC Fig.6 Q / Fig.8 T（已有短窗展示）
2. STIR Fig.8a（10 s 对折对比，短窗）
3. STIR Fig.11 — **等温三点**（本交接主缺口）
4. Fig.10 仅 Qmax 核对；Fig.8c/d 辅助；Fig.16 不做

## 6. 注意

- 法 ≠ 甲；NC ≠ 搅拌；勿混 Reload。
- Fig.11 必须常物性等温；热耦合 10 s 软点不能当验收。
- GPU：先 `nvidia-smi`；本开发机曾出现驱动不可用。
- 默认 `run_fr_stir.sh` 已改回短窗 `END=10`（不再默启 450 s）。
