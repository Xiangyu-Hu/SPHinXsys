# 230｜甲方算例物理落地（材料 / 几何 / 预设）

日期：2026-11-20  
依据：Doc228、Doc229、软件任务书、甲方几何（已放入对应 relax `data/`）与回复。  
本轮**不包含** Linux GUI、Windows、弹性桨 E0–F1 完成声明。

## 1. 已落地代码

| 文件 | 作用 |
|------|------|
| `electromagnetic_ophelie_client_material_laws.h` | 甲方 `mu=exp(-4+8500/T)`、`sigma=exp(4.5-8500/T)`；冻结 μ；法国 σ=16 借用路径；公式复算 CSV |
| `electromagnetic_ophelie_client_geometry.h` | 两套 STL 角色、AABB、线圈四环、名义热电偶、CAD→求解器变换说明 |
| `electromagnetic_ophelie_client_case_params.h` | NC/搅拌 development 预设；30 kW **电源输入**账本；effective JSON |
| `test_3d_ophelie_client_case_audit` | 主机审计入口：公式门控 + STL 清单 + geometry/effective JSON |
| `OphelieClientNcGlassCylinderShape` / `OphelieClientNcCrucibleWallShape` | 甲方 NC **参数化**液腔与数值壁（不用壳 STL 填流体） |
| `test_3d_ophelie_client_nc_glass_relax` | 甲方 NC 液腔粒子生成 + SYCL-CK relax + Reload；写出测点/材料 meta |
| `applyClientNcMeltToFrenchReduced` / `applyClientNcDevelopmentToLiterature` | 把甲方熔体/材料映射进现有 French-reduced knobs |
| `frozen_q` CLI `--client-nc-geometry` / `--client-nc-material` / `--zero-q` skip-EM | 甲方 Reload 等温静水入口 |

## 2. 材料口径（已编码）

- **作废**：`pow(10, 4.5-8500/T)` 旧电导率读法及 1473 K≈0.0536 / 1800 K≈0.599 旧表。
- **现行甲方公式**：自然指数；1473 K → μ≈5.874、σ≈0.280681（相对误差门 1e-4）。
- **NC 开发默认**：`french_table1_sigma16_borrowed` + **冻结** μ(1000℃)≈14.532 Pa·s。
- **对照预设**：同几何 / 同频率的 `client_ln_sigma_provisional`（供后续同电流未归一功率对照；电流仍未知，本轮不跑假 kA）。
- ρ=2700 暂定；cp/k/β 借用法国表1并标记 borrowed。
- **不改**法国回归材料 ID / 法国算例默认值。

## 3. 几何口径（已编码）

- STL 单位 mm；CAD 竖直 +Y；推荐变换：坩埚外底轴心→原点后 `(x,y,z)→(x,-z,y)/1000`。
- 自然对流无桨；搅拌桨偏心历史值约 −67 mm（须复核后才硬编码驱动轴）。
- 热电偶尖端 = AABB 下端名义点，`status=unconfirmed_sensor_location`。
- 水冷管 / 线圈展示；零面积三角形保留原件并记入说明。

## 4. 功率口径（已编码）

- `generator_input_power_W=30000`（confirmed）。
- `coil_current_A` / peak-vs-RMS = **unknown**（不填 1000 A）。
- 禁止把法国 50/60 kW 吸收功率直接复制为甲方预设。

## 5. 如何复跑审计

已跑通一次，产物在：

`docs/ophelie/reference/client_case_audit_2026-11-20/`

```bash
# 配置并编译 test_3d_ophelie_client_case_audit 后：
./test_3d_ophelie_client_case_audit \
  --repo-root /path/to/SPHinXsys \
  --output-dir docs/ophelie/reference/client_case_audit_2026-11-20
```

期望：`client_formula_check.csv`、两套 `*_geometry.json` / `*_stl_manifest.csv`、三份 `effective_*.json`、`AUDIT_SUMMARY.md`；进程退出码 0。

## 6. 编译与运行（请本机执行）

### 6.1 材料/几何审计（无需 GPU 长跑）

```bash
cd /home/ubuntu/SPHinXsys/build
cmake --build . --target test_3d_ophelie_client_case_audit -j$(nproc)

./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_case_audit/bin/test_3d_ophelie_client_case_audit \
  --repo-root /home/ubuntu/SPHinXsys \
  --output-dir /home/ubuntu/SPHinXsys/docs/ophelie/reference/client_case_audit_2026-11-20
```

### 6.2 甲方 NC 液腔 relax（需 SYCL/GPU 构建）

```bash
cd /home/ubuntu/SPHinXsys/build
cmake --build . --target test_3d_ophelie_client_nc_glass_relax -j$(nproc)

cd tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_nc_glass_relax/bin
./test_3d_ophelie_client_nc_glass_relax \
  --dp=0.010 \
  --melt-height=0.180 \
  --relax-steps=1000 \
  --relax-vtp-every=100 \
  --output-meta-dir=./client_nc_relax_meta
```

可选：`--client-ln-sigma` 切换甲方 ln-σ 对照材料；默认法国 σ=16 借用 + 冻结 μ(1000℃)。

检查：`client_nc_relax_meta/{melt_geometry.json,effective-case.json,probe_nominal.csv,README.md}`，以及 Reload / VTP。

液深 `0.180 m` 为开发暂定，甲方确认后改 `--melt-height=`。

## 7. 等温静水短跑（Reload → frozen_q）

`test_3d_ophelie_french_natural_convection_frozen_q` 已接：

- `--client-nc-geometry`：甲方 NC 参数化液腔（R≈0.1415、H 默认 0.180、底 z≈0.001、4 匝线圈、100 kHz）
- `--client-nc-material`：冻结 μ(1000℃) + σ16 借用；关闭浮力；**不要**用法 50 kW 标定
- `--zero-q`：**跳过 EM Phase A**，`frozen_q=0`（等温静水）
- `--melt-height=`：与 relax 同口径（仅配合 `--client-nc-geometry`）

不要单独用 `--glass-height=`（会把 center 重置为 H/2，与甲方底面 z=0.001 冲突）；需要改高度用 `--melt-height=` 或同时设 `--glass-center=`。

```bash
cd /home/ubuntu/SPHinXsys/build
# 若刚改过源码，先 reconfigure 再编（未知 target 时）：cmake .
cmake --build . --target test_3d_ophelie_french_natural_convection_frozen_q -j$(nproc)

FROZEN_BIN=tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/bin
RELOAD=../test_3d_ophelie_client_nc_glass_relax/bin/reload   # 相对 FROZEN_BIN

cd "$FROZEN_BIN"
./test_3d_ophelie_french_natural_convection_frozen_q \
  --client-nc-geometry \
  --client-nc-material \
  --zero-q \
  --zero-buoyancy \
  --dp=0.010 \
  --melt-height=0.180 \
  --reload-dir="$RELOAD" \
  --end-time=0.50 \
  --hydro-settle-time=0.15 \
  --state-record-interval=0.05 \
  --summary-json=./client_nc_isothermal_hydro_summary.json
```

绝对路径也可：`--reload-dir=/home/ubuntu/SPHinXsys/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_nc_glass_relax/bin/reload`（目录内须有 `Reload.xml`）。

期望日志含：`[ophelie][client-nc] geometry`、`[ophelie][client-nc] material`、`[ophelie][stage4.1] zero_q: skip EM`；静水后 `escaped≈0`、`U_max ≪ c0`。

短冒烟可用 `--end-time=0.05`；正式静水验收用上面的 `0.50`。未显式 `--reload-dir=` 且带 `--client-nc-geometry` 时，会优先找 `../test_3d_ophelie_client_nc_glass_relax/bin/reload`。

**已知坑（已修）**：
1. `--zero-q` 跳过 Doc219 spatial-q coverage abort（否则 `n_uncovered=N`）。
2. `--zero-q` 在 settle 后补清 `DensityChangeRate/Force/Displacement`；Doc224 blowup 门改为 `0.05·c0`（约 0.66 m/s），避免静水残速 ~0.05 误杀。
3. `material_check passed=0` 在甲方冻结 μ 路径上是预期日志，不中止。

### 7.1 本机结果（2026-11-21，exit 0）

| 量 | 值 |
|----|----|
| `passed` / `smoke_passed` | 1 / 1 |
| `escaped` | 0/11232 |
| end `U_max` / `U_rms` | 0.0135 / 0.0038（settle 后衰减） |
| `P_joule` / `em_updates` | 0 / 0 |
| `boussinesq` | 0（`--zero-buoyancy`） |
| 备注 | Robin 冷却使 `dT` 缓增到 ~3.7 K（无扩散、无 Q）；等温静水门看 escaped + U 衰减，不钉绝对 T |

产物：`bin/client_nc_isothermal_hydro_summary.json`、`output/GlassBody_ite_*.vtp`。

## 8. 受控温差短冒烟（已跑）

冻结 μ≈14.5 下 Ra 被埋；本步只验：**温差场可施加 + 有/无浮力对照不炸**，**不**宣称可观测自然对流。

### 8.1 本机结果（2026-11-21，两跑均 exit 0）

| 跑 | boussinesq | dT | U_max | U_rms | escaped | passed |
|----|------------|-----|-------|-------|---------|--------|
| A `client_nc_dT20_boussinesq_summary.json` | 1（β=1e-5） | 19.82 | 0.01350 | 0.00381 | 0 | 1 |
| B `client_nc_dT20_zero_buoyancy_summary.json` | 0 | 19.82 | 0.01350 | 0.00379 | 0 | 1 |

结论：冻结径向 ΔT=20 K 稳定；A/B 速度几乎相同 → **在冻结 μ 开发预设下浮力信号不可分辨**（与 `buried_defer_finer_mesh` 一致）。等温静水 + 受控温差冒烟均 PASS。

## 9. Doc228 S0/T0：用仓库现成算例写报告（不自造通用力学 case）

用户口径（2026-11-21）：**纯模块独立验证**不必新写悬臂/导热 case；用 SPHinXsys 自带固体/流体/传热算例即可。自写的 `test_3d_ophelie_doc228_s0_cantilever` **作废**（且 `DistributingPointForces` 要壳 `Thickness`，体积梁会 `Thickness is nullptr`）。

S0/T0 报告证据（本机为 **SYCL 构建**：上游 `tests/3d_examples` 在 SYCL 下不编入；用 `extra_source` 里已有算例）：

```bash
cd /home/ubuntu/SPHinXsys/build

# S0 固体：现成 LG 被动悬臂（extra，可 SYCL 编）
cmake --build . --target test_3d_passive_cantilever_LG -j$(nproc)
cd tests/extra_source_and_tests/3d_examples/test_3d_passive_cantilever_LG/bin
./test_3d_passive_cantilever_LG

# T0 传热：OPHELIE 已有 1D 导热解析门（不是新造框架）
cd /home/ubuntu/SPHinXsys/build
cmake --build . --target test_3d_ophelie_thermal_1d_conduction -j$(nproc)
cd tests/extra_source_and_tests/3d_examples/test_3d_ophelie_thermal_1d_conduction/bin
./test_3d_ophelie_thermal_1d_conduction --summary-json=./t0_thermal_summary.json
```

非 SYCL 构建还可补跑上游 `test_3d_mr_cantilever_beam` / `test_3d_passive_cantilever`。  
不锈钢 304 预设头 `electromagnetic_ophelie_client_stainless_steel.h` 留给 **甲方桨 E0+**。

## 10. 甲方搅拌路线（对齐法国搅拌 case 结构）

用户口径（2026-11-21）：**按法国搅拌 case 装配甲方电磁+搅拌**；桨改为弹性力学；材料用开发说明书（Doc228/229）。不是另起一套通用力学玩具。

| 层 | 法国参考 | 甲方目标 |
|----|----------|----------|
| Relax | `french_stirring_glass_relax` → Glass/Rotor/Wall Reload | 同结构；熔体参数化液腔−桨；桨=甲方 STL |
| EM+流动 | `french_stirring_em`（host-spin 刚桨） | 同 Phase A/B；桨=`elastic_driven` + 304 |
| 材料 | 法国表1 / 60 kW | Doc229 冻结μ/σ预设；30 kW **电源输入**；电流 unknown |

头文件：`electromagnetic_ophelie_client_stirring_geometry.h`（脚手架）。  
E0（桨单独 hold）是接入全案前的 Doc228 固体隔离步，**不是**最终交付形态。

### 10.1 E0 桨 hold（源码树路径；先修 BoundingBox 后再跑）

```bash
SRC=/home/ubuntu/SPHinXsys
cd "$SRC/build"
cmake .
cmake --build . --target test_3d_ophelie_doc228_e0_paddle_hold -j$(nproc)

ST_RELAX=$SRC/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_glass_relax
# CAD + 预处理桨在 stirring relax/data/（CMake 会拷到 E0 bin/data）
python3 "$ST_RELAX/preprocess_client_paddle_stl.py"
BINE0=$SRC/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_doc228_e0_paddle_hold/bin

OUT=$SRC/docs/ophelie/reference/elastic_paddle_runs_2026-11-21/E0_hold
mkdir -p "$OUT"
cd "$BINE0"
./test_3d_ophelie_doc228_e0_paddle_hold \
  --stl=./data/client_stir_paddle_solver_m.stl \
  --dp=0.004 --solid-steps=100 \
  --out-dir="$OUT" --summary-json=summary.json
```

甲方 CAD：`client_nc_glass_relax/data/`、`client_stirring_glass_relax/data/`（已从 `docs/ophelie/software/*stl` 迁入）。  
其后：client stirring relax → F0 等温弹性 FSI。

### 10.2 E0 本机结果（2026-11-21，exit 0）

| 量 | 值 |
|----|-----|
| `PADDLE_HOLD_PASS` | 1 |
| n / n_clamp | 10119 / 649（夹持约 6.4%，非整桨） |
| E | 2e11 Pa（304，未降 E） |
| max_disp / vonMises | 0 / 0（无载等温 hold） |
| WARN | `dp=0.004` 叶片可能欠分辨；E1 应力前加密 |

产物：`docs/ophelie/reference/elastic_paddle_runs_2026-11-21/E0_hold/summary.json`

### 10.3 F0 等温弹性 FSI 短跑（编译 + 运行；本机 GPU 跑）

前置：`test_3d_ophelie_client_stirring_glass_relax --bodies=all` 已写出 Reload（`n_rotor≈4332`，`dp≈0.005`）。

```bash
SRC=/home/ubuntu/SPHinXsys
cd "$SRC/build"
cmake .
cmake --build . --target test_3d_ophelie_client_stirring_f0_elastic_fsi -j$(nproc)

RELAX_RELOAD=$SRC/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_glass_relax/bin/reload
BINF0=$SRC/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_f0_elastic_fsi/bin
# CMake 已从 stirring relax/data 拷贝 client_stir_paddle_solver_m.stl → bin/data
mkdir -p "$BINF0/output"

OUT=$SRC/docs/ophelie/reference/elastic_paddle_runs_2026-11-21/F0_isothermal_fsi
mkdir -p "$OUT"
cd "$BINF0"
./test_3d_ophelie_client_stirring_f0_elastic_fsi \
  --reload-dir="$RELAX_RELOAD" \
  --rotor-stl=./data/client_stir_paddle_solver_m.stl \
  --end-time=0.05 --omega-ramp=0.05 --c0=8 \
  --clamp-top-band=0.05 \
  --max-wall-hours=2 \
  --out-dir="$OUT" --summary-json=summary.json
```

要点：桨=`elastic_driven` + 304（E=2e11，不降 E）；关 EM / Q / Boussinesq；重力保留；夹持=桨顶 50 mm 暂定带；监控 CSV=`output/client_f0_elastic_fsi_monitor.csv`。不锈钢固体声速高 → 固体亚步多，墙钟可能较长。

### 10.4 F0 本机结果（2026-11-21，exit 0）

| 量 | 值 |
|----|-----|
| `ELASTIC_FSI_SHORT_PASS` | 1 |
| t / wall | 0.050 s / ≈41 s |
| n_glass / n_rotor / n_clamp | 43917 / 4332 / 278（夹持≠整桨） |
| E | 2e11 Pa（未降） |
| f_pres / f_visc / Tq | ≈1.79 N / 0.44 N / −0.012 Nm（双向载荷非零） |
| tip / e_off | ≈0.071 m/s / ≈0.22 mm |
| solid_steps | 70759 |

产物：`docs/ophelie/reference/elastic_paddle_runs_2026-11-21/F0_isothermal_fsi/summary.json`

### 10.5 F1 加热弹性 FSI（规定体积热源；非 EM）

Doc228：F0 后恢复稳定热源路径；本轮选 **规定体积热**（电流未知，不得标成“已算 EM”）。常黏度、无浮力、无桨热应变；热电偶 Observer 采样 T。

```bash
SRC=/home/ubuntu/SPHinXsys
cd "$SRC/build"
cmake .
cmake --build . --target test_3d_ophelie_client_stirring_f1_heated_elastic -j$(nproc)

RELAX_RELOAD=$SRC/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_glass_relax/bin/reload
BINF1=$SRC/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_f1_heated_elastic/bin
OUT=$SRC/docs/ophelie/reference/elastic_paddle_runs_2026-11-21/F1_heated_integration
mkdir -p "$OUT" "$BINF1/output"
cd "$BINF1"
./test_3d_ophelie_client_stirring_f1_heated_elastic \
  --reload-dir="$RELAX_RELOAD" \
  --rotor-stl=./data/client_stir_paddle_solver_m.stl \
  --end-time=0.2 --omega-ramp=0.05 --c0=8 \
  --heat-power=1e5 \
  --clamp-top-band=0.05 \
  --max-wall-hours=2 \
  --out-dir="$OUT" --summary-json=summary.json
```

默认 `P=1e5 W` 为短窗可见 ΔT 的开发热功率（非 30 kW 电源输入）。看 `HEATED_ELASTIC_INTEGRATION_PASS`、`energy_rel_err`、`probe_T_K`。

### 10.6 E1 非零载荷响应（无流体；不锈钢 E）

Doc228：夹持保留，对自由区施加两级已知总力（质量加权写入 `ForcePrior`），检查可分辨弹性偏移、应力更新、夹持静止及近似比例关系。不降低 E。

```bash
SRC=/home/ubuntu/SPHinXsys
cd "$SRC/build"
cmake .
cmake --build . --target test_3d_ophelie_doc228_e1_paddle_load -j$(nproc)

BINE1=$SRC/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_doc228_e1_paddle_load/bin
OUT=$SRC/docs/ophelie/reference/elastic_paddle_runs_2026-11-21/E1_load_response
mkdir -p "$OUT" "$BINE1/output"
cd "$BINE1"
./test_3d_ophelie_doc228_e1_paddle_load \
  --dp=0.004 --solid-steps=5000 --clamp-top-band=0.05 \
  --force1=1.5e4 --force2=3e4 \
  --out-dir="$OUT" --summary-json=summary.json
```

看 `LOAD_RESPONSE_PASS`、`ELASTIC_PADDLE_ACTIVE`、`e_off_ratio` / `ratio_rel_err`。粗 `dp` 下不宣称应力精度。

### 10.7 E1 本机结果（2026-11-21，exit 0）

| 字段 | 值 |
|------|-----|
| `LOAD_RESPONSE_PASS` | 1 |
| `ELASTIC_PADDLE_ACTIVE` | 1 |
| F1 / F2 | 1.5e4 / 3.0e4 N（+Y，自由区质量加权） |
| e_off_free 1/2 | ≈17.3 / 34.5 mm |
| e_off_ratio / force_ratio | ≈1.991 / 2（rel_err≈0.45%） |
| sigma/E 1/2 | ≈0.0043 / 0.0086 |
| clamp_e_off | 0 |
| wall_s | ≈48 |

产物：`docs/ophelie/reference/elastic_paddle_runs_2026-11-21/E1_load_response/summary.json`

### 10.8 E2 规定旋转短窗（无流体；10 rpm）

Doc228：夹持区平滑转速斜坡，先 0.05–0.1 s 冒烟。检查轴、符号、约束与固体亚步；弹性偏移相对刚体参考轨迹。整圈 2π 不在本短窗宣称。

```bash
SRC=/home/ubuntu/SPHinXsys
cd "$SRC/build"
cmake .
cmake --build . --target test_3d_ophelie_doc228_e2_paddle_spin -j$(nproc)

BINE2=$SRC/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_doc228_e2_paddle_spin/bin
OUT=$SRC/docs/ophelie/reference/elastic_paddle_runs_2026-11-21/E2_driven_rotation
mkdir -p "$OUT" "$BINE2/output"
cd "$BINE2"
./test_3d_ophelie_doc228_e2_paddle_spin \
  --dp=0.004 --end-time=0.05 --omega-ramp=0.05 --clamp-top-band=0.05 \
  --out-dir="$OUT" --summary-json=summary.json
```

看 `DRIVEN_ELASTIC_ROTATION_SHORT_PASS`、`theta_rad`、`sign_ok`、`max_clamp_err_m`。

### 10.9 E2 本机结果（2026-11-21，exit 0）

| 字段 | 值 |
|------|-----|
| `DRIVEN_ELASTIC_ROTATION_SHORT_PASS` | 1 |
| `full_revolution_pass` | 0（短窗不宣称） |
| t / θ / ω | 0.05 s / ≈0.026 rad / ≈1.047 rad/s（10 rpm） |
| max_e_off / clamp_err | ≈0.33 mm / 0 |
| tip_u / sign_ok | ≈0.068 m/s / 1 |
| solid_steps / wall_s | 73645 / ≈346 |

产物：`docs/ophelie/reference/elastic_paddle_runs_2026-11-21/E2_driven_rotation/summary.json`

## 11. 仍未做

1. 电流到位后：σ16 vs ln-σ 未归一 Q 对照；再换真实 EM Q。  
2. Doc228 E2 整圈（∫ω dt = 2π；本短窗不宣称）。  
3. 测点时间序列 + 1000℃ 窗口（不得钉测点温度冒充 PASS）。  
4. ~~两案软件交接包~~ → 见 `docs/ophelie/reference/client_software_handoff_2026-11-21/`（HANDOFF.md + COMPILE_AND_RUN.sh）。整圈 / EM Q / 1000℃ 定量仍开放。

## 12. 给软件机的契约提示

读取 `effective_*.json` 的 `material_id`、`power.generator_input_power_W`、`mechanics_mode`。  
自然对流短跑与搅拌 Doc228 门控（E0/E1/E2 短窗、F0、F1）本机已过；软件机可接两案预设，但不得把短窗标成整圈 PASS，也不得把规定体积热标成已算 EM。
