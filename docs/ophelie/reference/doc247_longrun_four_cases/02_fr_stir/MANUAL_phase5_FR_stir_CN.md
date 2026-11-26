# 法国搅拌 Phase5 说明书（已合入本机）

**来源**：`origin/feature/ophelie-fr-stir-phase5-4090` @ `df3185e38`（文件合入，未 merge 该分支）  
**口径**：ASSUMED 工程演示；不宣称文献 PASS  
**桨机构**：Simbody **刚体规定转动**（不是甲方 C2 NeoHookean 弹性桨）

---

## 1. 合入内容

| 类别 | 落点 |
|------|------|
| 真 Fig.10 Q 种子 | `standard_data/doc247_seeds/SEED_STIR_Fig10_Qshape_ASSUMED_*`（旧 PLACEHOLDER 备份为 `*_PLACEHOLDER_from_Fig6_*`） |
| Fig.8a T 种子 | 同目录（与 handoff 一致） |
| 数字化脚本/图 | `standard_data/doc247_french_digitization/` |
| Q-correct C++ | 共享分支已有（`--fig10-q-correct`，默认关） |
| 跑批脚本 | 本目录 `run_fr_stir*.sh` |

---

## 2. 推荐物理口径（10 s）

| 项 | 设置 |
|----|------|
| 温度 IC | Fig.8a → `SEED_STIR_Fig8a_axisym_ASSUMED_T_rz.csv` |
| 焦耳热形状 | 真 Fig.10 Q + `--fig10-q-correct` blend=0.8 |
| 功率 | fixed-power **60 kW** |
| 耦合 | periodic，`em-update-interval=0.5` |
| 热阻 | skull + free-crust（delta=0.005） |
| 转速 | 10 rpm（矩阵另含 5/20） |

注意：**Fig.10 不是温度分布**，只是 \(Q\) 形状。

---

## 3. 两案对照（常粘 / 变粘）

| 脚本 | 开关 | 含义 |
|------|------|------|
| `run_fr_stir_muT_const10s.sh` | `--constant-thermal-rheology` | \(\mu\equiv 4\,\mathrm{Pa\cdot s}\) |
| `run_fr_stir_muT_varT10s.sh` | `--temperature-dependent-properties` | \(\mu(T),c_p(T),k(T)\) CES Fig.2 |

详情见 [SUMMARY_muT_setup_CN.md](SUMMARY_muT_setup_CN.md)。

```bash
cd /home/ubuntu/SPHinXsys
export ONEAPI_DEVICE_SELECTOR=cuda:0
bash docs/ophelie/reference/doc247_longrun_four_cases/02_fr_stir/run_fr_stir_muT_const10s.sh
bash docs/ophelie/reference/doc247_longrun_four_cases/02_fr_stir/run_fr_stir_muT_varT10s.sh
```

基线（常粘 + Fig10 Q-correct）：

```bash
bash docs/ophelie/reference/doc247_longrun_four_cases/02_fr_stir/run_fr_stir.sh
```

桨高扫：`run_fr_stir_height_relax_and_10s.sh`  
矩阵：`run_fr_stir_phase5_matrix_60s.sh`

---

## 4. 变粘摘要

- Jacoutot CES 2008 Fig.2 数字化表；\(\log\mu\)–\(T\) 插值；1473 K 锚 \(\mu=4\)。
- 水力用体积平均 \(\langle\mu(T)\rangle\)，不是粒子对局部 μ。
- Fig.8a 偏冷 → 终点 \(\mu_\mathrm{mean}\approx 224\,\mathrm{Pa\cdot s}\)；过程中可有 \(U_\mathrm{max}>c_0\)，流场标 ASSUMED。

---

## 5. 报告可用一句话

> 在 Fig.8a 温场初值、真 Fig.10 Q 形状校正与 60 kW 周期电磁耦合下，用刚体规定转动桨完成转速/桨高短时矩阵，并对照常粘与 CES 2008 Fig.2 变粘 \(\mu(T)\)（体积平均）；结果为工程演示，不宣称文献 PASS。
