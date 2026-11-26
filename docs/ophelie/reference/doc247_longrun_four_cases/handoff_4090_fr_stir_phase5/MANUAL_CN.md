# 本机法国搅拌 Phase5 说明书（隔离 handoff）

**机器**：4090 工作站（sphinxsys-ml）  
**日期**：2026-11-26  
**口径**：ASSUMED 工程演示；不宣称文献 PASS  
**专用分支**：`feature/ophelie-fr-stir-phase5-4090`  
**共享分支**：`feature/EM-Biot-Savart`（另一台可继续在此分支改自己的内容）

---

## 0. 为什么要隔离

本机在法国搅拌上做了较多工作（Fig.10 真数字化、转速/桨高 10 s 扫、μ(T) 对照）。  
共享分支 `feature/EM-Biot-Savart` 上另一台电脑也有修改。若继续把本机后续改动直接提交到共享分支，容易在同一路径（尤其种子 CSV、`.cpp`）产生冲突。

**做法：**

1. 已推上共享分支的 **Fig.10 Q-correct 代码能力保留**（默认关闭，见 §2），**不 revert**。
2. 本机后续文档、真种子、跑批脚本、说明书只放在本目录，并只提交到 **`feature/ophelie-fr-stir-phase5-4090`**。
3. 工作区对共享路径上的未提交脏改动，在拷入本 handoff 后 **restore 还原**，避免本机继续脏改共享种子。

另一台若需要本机搅拌成果：merge 或 cherry-pick `feature/ophelie-fr-stir-phase5-4090` 即可。

---

## 1. 本机工作范围一览

| 类别 | 内容 | 落点 |
|------|------|------|
| 代码（已在共享分支 tip `a502e16d8`） | Fig.10 Q 形状校正头文件 + stirring_em CLI | 见 §2 |
| 真 Fig.10 数字化 | 替换早期 PLACEHOLDER（NC Fig.6 重映射） | 本目录 `seeds/`、`scripts/digitize_stir_fig10.py` |
| 转速扫 10 s | 5 / 10 / 20 rpm | 本机 `overnight_6h/phase5_fr_stir_rpm10s/`（gitignore） |
| 桨高扫 10 s | Δz = −20 / 0 / +20 mm | 同左 `phase5_fr_stir_height10s/` |
| 六算例打包 | T / 流线 / midplane U CSV | `phase5_fr_stir_6cases_pack_20261126.tgz` |
| μ(T) 对照 | 常粘 vs 变粘，10 rpm × 10 s | `phase5_fr_stir_muT10s/` + 本目录 `SUMMARY_muT_setup_CN.md` |

---

## 2. 已合入共享分支的代码能力（默认关闭）

### 2.1 新头文件（相对独立）

`tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_doc247_q_shape_correct.h`

在 EM→Q 之后，按种子 Q 形状做 blend，再 rescale 到 `target_absorbed_power`。

### 2.2 接线

`tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_em/test_3d_ophelie_french_stirring_em.cpp`

| CLI | 默认 | 含义 |
|-----|------|------|
| `--fig10-q-correct` | 关 | 打开 Fig.10 Q 校正 |
| `--fig10-q-correct-blend=` | 0.8 | 形状权重 ∈[0,1] |
| `--fig10-q-correct-every=` | 1 | 每 N 次 EM 更新校正一次 |
| `--init-Q-rz-csv=` | — | Q 形状种子路径 |
| `--init-T-rz-csv=` | — | 温度初场（本机用 Fig.8a） |
| `--temperature-dependent-properties` | 关 | 打开 μ/cp/k(T) |
| `--constant-thermal-rheology` | 开 | 钉在 \(T_0\) 常值 |

**对其他算例影响**：不传 `--fig10-q-correct` 时行为与原先一致。  
**冲突风险**：另一台若也改了同一个 `test_3d_ophelie_french_stirring_em.cpp`，合并时可能只在该文件冲突；新头文件一般可干净合入。

---

## 3. 本机推荐物理口径（10 s 工程对照）

两边对照（常粘 / 变粘）除流变开关外应一致：

| 项 | 设置 |
|----|------|
| 温度 IC | **Fig.8a** → `seeds/SEED_STIR_Fig8a_axisym_ASSUMED_T_rz.csv` |
| 焦耳热形状 | **Fig.10 Q** → `seeds/SEED_STIR_Fig10_Qshape_ASSUMED_rz.csv` + `--fig10-q-correct` blend=0.8 |
| 功率 | fixed-power **60 kW** |
| 耦合 | periodic，`em-update-interval=0.5` |
| 几何 | `reload_HZ_z00`（Δz=0 基线桨高） |
| 热阻 | skull + free-crust（delta=0.005） |
| 转速 | 10 rpm（矩阵另含 5/20；桨高扫另含 ±20 mm） |
| 结束时间 | 10 s |

注意：**Fig.10 不是温度分布**，只是 \(Q\) 形状；温度来自 Fig.8a。

---

## 4. 本机脚本（本目录 `scripts/`）

| 脚本 | 用途 |
|------|------|
| `digitize_stir_fig10.py` | 数字化文献 Fig.10 → Q 种子 |
| `run_fr_stir_phase5_matrix_60s.sh` | Phase5 矩阵脚手架（含 60 s 口径；本机后续多用 10 s） |
| `run_fr_stir_height_relax_and_10s.sh` | 桨高 relax + 10 s |
| `run_phase5_muT10s_gpu0.sh` | 常粘 + 变粘 + 变粘@20rpm（GPU0） |
| `run_phase5_muT10s_gpu0_varT_only.sh` | **仅变粘** `MU_varT_rpm10`（推荐） |

示例（tmux + GPU0，仅变粘）：

```bash
tmux new -s phase5_muT_gpu0   # 或 attach
export ONEAPI_DEVICE_SELECTOR=cuda:0
bash /home/yongchuan/sphinxsys/docs/ophelie/reference/doc247_longrun_four_cases/handoff_4090_fr_stir_phase5/scripts/run_phase5_muT10s_gpu0_varT_only.sh
```

脚本内种子路径仍指向仓库标准 `doc247_seeds/`（共享路径上的已提交 PLACEHOLDER 版）。  
若要用本 handoff 里的**真 Fig.10**，把脚本里的 `SEED_Q=` 改为本目录：

```text
.../handoff_4090_fr_stir_phase5/seeds/SEED_STIR_Fig10_Qshape_ASSUMED_rz.csv
```

（本机跑 μ(T) 时工作区种子已是真 Fig.10；restore 共享路径后请改用上述 handoff 路径。）

---

## 5. 变粘性 μ(T)（摘要）

详情见同目录 [SUMMARY_muT_setup_CN.md](SUMMARY_muT_setup_CN.md)。

要点：

- 不是改一个常数 μ，而是 Jacoutot CES 2008 **Fig.2** 数字化表 + \(\log\mu\)–\(T\) 插值；同时开 \(c_p(T),k(T)\)。
- 水力用 **体积平均** \(\langle\mu(T)\rangle\)，不是粒子对局部 μ。
- Fig.8a 场偏冷 → 本算例终点 \(\mu_\mathrm{mean}\approx 224\,\mathrm{Pa\cdot s}\)（常粘对照为 4）。
- `MU_varT_rpm10`：`exit=0`，`diverged=0`，但过程中 `U_max≃15 > c_0=8`，流场标 ASSUMED。

---

## 6. 本机大数据包位置（不入库）

根目录（gitignore）：

```text
/home/yongchuan/sphinxsys/docs/ophelie/reference/doc247_longrun_four_cases/overnight_6h/
```

| 包 / 目录 | 内容 |
|-----------|------|
| `phase5_fr_stir_6cases_pack_20261126.tgz` | 转速+桨高 六算例 + midplane T/流线 |
| `phase5_fr_stir_muT10s_pack_20261126.tgz` | 常粘 vs 变粘画图包 |
| `phase5_fr_stir_rpm10s/`、`height10s/`、`muT10s/` | 原始 case 输出 |

曾上传 Google Drive（若仍保留）：`google-drive:file transfer/sphinxsys_ophelie/` 一类路径（以当时 rclone 日志为准）。

---

## 7. 与另一台协作：冲突规避清单

| 建议 | 说明 |
|------|------|
| 另一台继续用 `feature/EM-Biot-Savart` | 改自己的非搅拌或其它文件 |
| 本机搅拌后续只提交 `feature/ophelie-fr-stir-phase5-4090` | 本 handoff 树 |
| 需要搅拌能力时再 merge handoff 分支 | 预期冲突热点：`test_3d_ophelie_french_stirring_em.cpp` |
| 不要在两边同时改 `SEED_STIR_Fig10_*` 共享路径 | 真种子以本 `seeds/` 为准 |
| 不要 revert 共享分支上的 Q-correct | 默认关；其它 case 不受影响 |

---

## 8. 关键代码 / 文档路径速查

```text
# Q-correct（共享分支已有）
tests/.../electromagnetic_ophelie_doc247_q_shape_correct.h
tests/.../test_3d_ophelie_french_stirring_em/test_3d_ophelie_french_stirring_em.cpp

# 本 handoff
docs/ophelie/reference/doc247_longrun_four_cases/handoff_4090_fr_stir_phase5/

# 本机结果（gitignore）
docs/ophelie/reference/doc247_longrun_four_cases/overnight_6h/
```

---

## 9. 报告可用一句话

> 本机在与 10 rpm 基线相同的 Fig.8a 温场初值、Fig.10 Q 形状校正与 60 kW 周期电磁耦合下，完成了转速/桨高短时矩阵，并对照开启 CES 2008 Fig.2 数字化 μ(T)/cp(T)/k(T)（水力为体积平均 μ）；相关说明书、真 Fig.10 种子与跑批脚本隔离在 `handoff_4090_fr_stir_phase5` 与分支 `feature/ophelie-fr-stir-phase5-4090`，以避免与共享分支上其它机器的并行修改冲突。
