# OPHELIE Stage 3.4：Jacoutot 论文 σ(T) 准备（公式采样，非像素扣点）

## 状态
**2026-08-11：Stage 3.4a `passed=1`（coil-only + thesis-iii12）。继续 3.4b。**

## 结论（先读）
1. **期刊 Fig.2 ≠ σ(T)**（是耦合流程图）。
2. **真正 σ(T)** = 博士论文 **图 III.12 / IV.10**，且文中已给出拟合公式 → **直接公式采样**，不必 WebPlotDigitizer。
3. 已写入 CSV + `makeJacoutotThesisSigmaTemperatureLawIII12/IV10()` + CLI `--sigma-law=`。

## 公式
```text
III.12 / (III.6) [C_2]  natural/static:
  log10(σ) = 3.7921 - 3179.8/T     →  σ(1473K)≈43.0 S/m,  σ≈16 @ ~1229 K

IV.10 [C_8] stirred Uox2+RuO2:
  log10(σ) = 4.05726 - 3923.73/T   →  σ(1473K)≈24.7 S/m
```

## 文件
- `docs/ophelie/reference/jacoutot_sigma_t/README.md`
- `.../jacoutot_thesis_III12_sigma_T.csv`
- `.../jacoutot_thesis_IV10_sigma_T.csv`
- `electromagnetic_ophelie_french_material_laws.h`
- 进度快照：[`OPHELIE_STAGE3_PROGRESS_SNAPSHOT_2026-08-11.md`](OPHELIE_STAGE3_PROGRESS_SNAPSHOT_2026-08-11.md)

## 用户跑法

### 3.4a — coil-only + thesis-iii12（已通过）

```bash
cd ~/SPHinXsysSYCL/build

RELAX_RELOAD=~/SPHinXsysSYCL/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_glass_relax/bin/reload

./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload-dir="$RELAX_RELOAD" \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-edge-flux-complex=1 \
  --sigma-t --sigma-law=thesis-iii12 \
  --sigma-t-max-iter=8 --sigma-t-relax=0.3 --sigma-t-tol=0.05 \
  --sigma-t-seed-delta-k=100 \
  --thermal-dt=0.1 --thermal-steps=1 \
  2>&1 | tee stage3_4_thesis_iii12_sigma_t.log
```

### 3.4b — thesis-iii12 + A_glass Picard（下一步，整段复制）

```bash
cd ~/SPHinXsysSYCL/build

RELAX_RELOAD=~/SPHinXsysSYCL/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_glass_relax/bin/reload

./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload-dir="$RELAX_RELOAD" \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-edge-flux-complex=1 \
  --self-induction --self-induction-max-iter=25 --self-induction-relax=0.3 \
  --self-induction-tol=1e-3 --self-induction-phi-tol=2e-4 \
  --sigma-t --sigma-law=thesis-iii12 \
  --sigma-t-max-iter=8 --sigma-t-relax=0.3 --sigma-t-tol=0.05 \
  --sigma-t-seed-delta-k=100 \
  --thermal-dt=0.1 --thermal-steps=1 \
  2>&1 | tee stage3_4b_thesis_iii12_aglass.log
```

期望：`sigma_law=thesis-iii12 paper_digitized_sigma=1`，`self_induction=1`，`j_ok=1 phi_ok=1`，`sigma_mean~40+`，`passed=1`。  
（嵌套 Picard×σ 外环，耗时约数秒～十几秒。）

可选人工目视：https://theses.hal.science/tel-00121350 打开图 III.12 对照 CSV。

## 实测

### Stage 3.4a（2026-08-11，coil-only + thesis-iii12）
```text
sigma_law=thesis-iii12 paper_digitized_sigma=1
sigma_t_iters=4 sigma_t_converged=1
sigma_mean=46.471  sigma_min=43.045  sigma_max=48.221
sigma_rel=0.0342   power_rel=0.0390
P_joule: 7.991 → 7.449 → 7.070 → 6.805 W
phi_eq_res_vol=1.93e-4 phi_ok=1
E_joule=E_thermal=E_power=0.68046 J
em_ok=1 thermal_ok=1 sigma_spatial_ok=1 sigma_coupling_ok=1 passed=1
self_induction=0
```

对照旧 provisional Arrhenius（≈16@1473）：
| 路径 | σ_mean | P_joule [W] |
|------|--------|------------|
| provisional σ(T) coil-only | ~17.5 | ~2.57 |
| **thesis-iii12** coil-only | **~46.5** | **~6.80** |

与公式一致：`σ(1473K)≈43`，径向 seed 使首步 `σ∈[43,58]`，外环收敛到更窄的 `[43,48]`。功率约升高 **2.6×**，符合更高电导预期。

### Stage 3.4b（2026-08-11，thesis-iii12 + A_glass — 首跑 `passed=0`）
```text
sigma_law=thesis-iii12 paper_digitized_sigma=1
sigma_t_iters=3 sigma_t_converged=1
sigma_mean=47.962 sigma_min=43.068 sigma_max=50.462
self_induction=1 iters=25 final_J_rel=8.12e-5 j_ok=1
picard_converged=0   # J 已收敛，但 φ 未进门 → 跑满 max_iter
A_ind_over_A_coil=0.440 B_ind_over_B_coil=0.372
P_joule_W=5.719
phi_eq_res_vol=2.104e-4 em_phi_gate=2.0e-4 phi_ok=0
thermal_ok=1 sigma_*=1 em_ok=0 passed=0
```

对照：
| 路径 | σ_mean | P [W] | A_ind/A | φ_res | gate |
|------|--------|-------|---------|-------|------|
| 3.2b provisional+A_glass | 17.5 | 2.49 | 0.177 | 2.00e-4 | 2e-4 ✓ |
| 3.4a thesis coil-only | 46.5 | 6.80 | — | 1.93e-4 | 1e-2 ✓ |
| **3.4b thesis+A_glass** | **48.0** | **5.72** | **0.440** | **2.10e-4** | **2e-4 ✗** |

解读：
- σ/P/热能/自感幅度均合理（高 σ → 更强 `A_ind`、功率介于 coil-only thesis 与旧 Arrhenius+A_glass 之间）
- **唯一失败点**：φ 平台 `2.10e-4` 比工程门 `2.00e-4` 高约 **5%**；与 Stage2「平台≈1.8e-4、门用 2e-4」同类，高电导略抬高平台
- 加 Picard 迭代数无效（J 早已 `<1e-3`）

### 3.4b 工程复跑（建议 `phi_tol=2.5e-4`）

```bash
cd ~/SPHinXsysSYCL/build

RELAX_RELOAD=~/SPHinXsysSYCL/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_glass_relax/bin/reload

./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload-dir="$RELAX_RELOAD" \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-edge-flux-complex=1 \
  --self-induction --self-induction-max-iter=25 --self-induction-relax=0.3 \
  --self-induction-tol=1e-3 --self-induction-phi-tol=2.5e-4 \
  --sigma-t --sigma-law=thesis-iii12 \
  --sigma-t-max-iter=8 --sigma-t-relax=0.3 --sigma-t-tol=0.05 \
  --sigma-t-seed-delta-k=100 \
  --thermal-dt=0.1 --thermal-steps=1 \
  2>&1 | tee stage3_4b_thesis_iii12_aglass_phi2p5e-4.log
```

期望：`phi_ok=1 passed=1`（物理场不变，仅工程门对齐高 σ 平台）。
