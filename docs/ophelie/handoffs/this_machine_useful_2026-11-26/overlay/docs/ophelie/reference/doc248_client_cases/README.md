# Doc248 甲方 NC（C0）+ 搅拌（C2 / C3 / C4 / C5）交接包

与 Doc247 法国线隔离。粒子**一律本机重新 relax 生成**，不附带现成 Reload。

## 搅拌路径（勿混比材料/机构/热源结论）

| Case | 角色 | 桨机构 / 热源 | 输出 |
|------|------|---------------|------|
| **C2** | 工程挂机 / 数值稳跑 | 软化弹性（默认 E=5e8）+ 体热源 | `runs/C2_stir/` |
| **C3** | 甲方不锈钢桨材料声明 | 弹性 E=200 GPa + 体热源 | `runs/C3_stir_stainless/` |
| **C4** | 长时流/热观效 | **刚体规定转动** + 体热源 | `runs/C4_stir_rigid/` |
| **C5** | 30 kW EM 标定 + 法式散热/环流 | `french_stirring_em` host-spin + **EM fixed-power 30 kW** + skull/crust（默认 `SKULL_DELTA=0.030`，净加热目标 ~18–20 kW） | `runs/C5_stir_em30/` |

共用搅拌 Reload；脚本与 `OUT` 独立。C2/C3 不可写成对方材料结论；C4 **不是** `elastic_driven`，仅用于长时看搅拌流场与规定热升温趋势。C5 **不是** C0 的 15 kW ASSUMED 吸收，也**不是** C4 假热源：电源 30 kW 仅作工程 EM 标定目标（`power_identity=generator_input_used_as_EM_scaling_target`），不宣称耦合效率 100%。

## 测点

- 主测点：`client_midradius_probe` — `r=0.5 R`，液面下 **0.12 m**，方位 **+X（ASSUMED）**
- 旧 CAD 热电偶 tip：单独命名，勿与 midradius 混比
- 映射：`configs/geometry_and_probe_mapping.json`

## 在本机

```bash
bash docs/ophelie/reference/doc248_client_cases/COMPILE.sh
```

## 拷到 GPU 机后

```bash
# 1) 编（若未带二进制）
bash docs/ophelie/reference/doc248_client_cases/COMPILE.sh

# 2) 从零生成粒子（NC 现为 H=0.230 m + 线圈随液位 remap；会删旧 reload）
bash docs/ophelie/reference/doc248_client_cases/01_generate_particles.sh

# 3) C0 自然对流（H=230 mm）
# 推荐：细热舌+轴向分层 IC + Doc246 上偏 Q（现有二进制即可；线圈压缩需先 COMPILE）
# OUT=docs/ophelie/reference/doc248_client_cases/runs/C0_nc_hybrid_Qup_10 \
# END_TIME=10 TARGET_ABSORBED=18000 H_SIDE=75 H_FREE=10 EMISSIVITY=0.4 \
# INIT_T_RZ=docs/ophelie/reference/doc248_client_cases/configs/client_nc_hybrid_tongue_H230_ASSUMED_T_rz.csv \
# Q_UPPER_KAPPA=2.5 \
# bash docs/ophelie/reference/doc248_client_cases/02_run_C0_nc.sh
# 重建后可再加：COIL_Z_FRAC_LO=0.55 COIL_Z_FRAC_HI=0.95
# Fig.6 映射 Q 形状参考（尚不能直接加载）：configs/client_nc_fig6_mapped_H230_ASSUMED_Q_rz.csv
bash docs/ophelie/reference/doc248_client_cases/02_run_C0_nc.sh

# 4a) C2 软化弹性挂机
bash docs/ophelie/reference/doc248_client_cases/03_run_C2_stir.sh

# 4b) C3 不锈钢弹性短窗探伤
bash docs/ophelie/reference/doc248_client_cases/04_run_C3_stir_stainless.sh

# 4c) C4 刚体桨（先短窗，再长跑看效果）
bash docs/ophelie/reference/doc248_client_cases/05_run_C4_stir_rigid.sh
END_TIME=60 MAX_WALL_HOURS=24 bash docs/ophelie/reference/doc248_client_cases/05_run_C4_stir_rigid.sh
# END_TIME=600 MAX_WALL_HOURS=48 bash ...

# 4d) C5：30 kW EM 标定 + skull/crust + 甲方搅拌几何/物性（默认 60 s 看净 ΔT）
END_TIME=60 bash docs/ophelie/reference/doc248_client_cases/06_run_C5_stir_em30.sh
```

输出在 `runs/`（gitignore）。C2/C3 弹性固体子步多、墙钟很长；C4 跳过弹性子步，长跑墙钟通常短得多。C5 走周期 EM，墙钟取决于 EM 刷新与 SPH 热耦合。

## 目标二进制

| Case | Relax | Run |
|------|-------|-----|
| C0 NC | `test_3d_ophelie_client_nc_glass_relax` | `frozen_q --client-nc-*` |
| C2 软化弹性 | `client_stirring_glass_relax --bodies=all` | F1 + `--youngs-modulus=` |
| C3 不锈钢弹性 | 同上 Reload | F1（不传 youngs → 200 GPa） |
| C4 刚体桨 | 同上 Reload | F1 + `--rigid-rotor` |
| C5 EM 30 kW | 同上 Reload | `french_stirring_em --client-stir-*` + fixed-power 30000 |
