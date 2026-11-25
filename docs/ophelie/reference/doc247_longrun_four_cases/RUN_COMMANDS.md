# Doc247 四套长时算例 — 编译与运行命令

物理机 `/home/ubuntu/SPHinXsys`。先一次性编译四个目标，再按需挂机。

## 0. 一次性编译

```bash
SRC=/home/ubuntu/SPHinXsys
BUILD=$SRC/build
cmake --build "$BUILD" --target \
  test_3d_ophelie_french_natural_convection_frozen_q \
  test_3d_ophelie_french_stirring_em \
  test_3d_ophelie_client_stirring_f1_heated_elastic \
  -j"$(nproc)"
```

公共路径：

```bash
SRC=/home/ubuntu/SPHinXsys
BUILD=$SRC/build
OUT_ROOT=$SRC/docs/ophelie/reference/doc247_longrun_four_cases
FROZEN=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/bin
FR_STIR=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_em/bin
CL_STIR=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_f1_heated_elastic/bin
RELOAD_FR_NC=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_glass_relax/bin/reload
RELOAD_FR_STIR=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_glass_relax/bin/reload
RELOAD_CL_NC=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_nc_glass_relax/bin/reload
RELOAD_CL_STIR=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_glass_relax/bin/reload
SEED=$SRC/docs/ophelie/reference/doc243_runs/03_thermal_initial_guess/SEED_coldwall_skin_ASSUMED_T_rz.csv
TABLES=$SRC/docs/ophelie/reference/doc239_runs/solver_tables
export ONEAPI_DEVICE_SELECTOR=cuda:gpu
export OMP_NUM_THREADS=8
export TBB_NUM_THREADS=8
```

---

## 1. 法国文献自然对流（最接近文献 · 长时）

**身份**：Doc246 `C0_BASE` = JAC_THESIS_RAW_BORROWED + coldwall 皮初场 + R10c skull/crust + 固定玻璃吸收 **50 kW** + paper flat slip top。  
**非** F8SEED 全域冷壳；**非** Q1_UPPER；**非** P65（P65 更贴能量残差，但整体环流/体平均仍以 C0 为生产锚）。

**长时软改动**：无 `--doc227-actual-timestep`；`--doc224-u-blowup-gate=50`；`end-time=600`（可改 1200）。

```bash
OUT=$OUT_ROOT/01_fr_nc/FR_NC_C0_LONGRUN_600s
mkdir -p "$OUT/out" "$OUT/figures"
cd "$FROZEN"
rm -rf output output_backup && mkdir -p output

./test_3d_ophelie_french_natural_convection_frozen_q \
  --reload-dir="$RELOAD_FR_NC" \
  --dp=0.015 \
  --doc223-paper-flat-slip-top \
  --em-control=fixed-power --target-absorbed-power=50000 \
  --thermo-em-coupling=periodic --em-update-interval=0.5 --aind=off \
  --enable-thermal-diffusion \
  --temperature-dependent-properties \
  --doc226-local-particle-viscosity \
  --viscosity-mode=implicit_split \
  --material-dataset=JAC_THESIS_RAW_BORROWED \
  --material-tables-dir="$TABLES" \
  --material-evidence-dir="$OUT" \
  --init-T-rz-csv="$SEED" \
  --skull-resistance \
  --skull-delta-side=0.005 --skull-delta-bottom=0.005 \
  --skull-k=1 \
  --skull-h-cool-side=300 --skull-h-cool-bottom=35 \
  --free-crust-resistance --free-crust-delta=0.005 --free-crust-k=1 \
  --doc224-u-blowup-gate=50 \
  --end-time=600.0 --hydro-settle-time=0.05 \
  --state-record-interval=30.0 \
  --meridional-average-start=60.0 \
  --boundary-face-t-mode=local_dt_dz \
  --summary-json="$OUT/summary.json" \
  2>&1 | tee "$OUT/run.log"

cp -f output/french_energy_budget.csv \
      output/french_thermo_em_coupling.csv \
      output/french_meridional_timeseries.csv \
      output/french_meridional_ur_uz_mean.csv \
      output/french_T_rz.csv \
      output/french_spatial_stats.csv \
      "$OUT/" 2>/dev/null || true
cp -a output/*.vtp "$OUT/out/" 2>/dev/null || true
```

启动应见：`dataset=JAC_THESIS_RAW_BORROWED`；`h_eff=120` / `29.787`；`with_skull_and_free_crust`。  
中止若仍出现：把 gate 再放宽或把 `end-time` 分段重开（无 checkpoint 时只能改初场/身份）。

---

## 2. 法国文献搅拌（长时 · 已有生产体）

**身份**：文献搅拌几何 + 10 rpm + 文献玻璃吸收功率（case 内置，约 60 kW）+ Round 3D Level-1 skull + free-crust。  
**勿**与 `--balance-heat-loss` 同开 skull。默认已按长时设计：`end-time=1200`、`max-wall-hours=10`。

```bash
OUT=$OUT_ROOT/02_fr_stir/FR_STIR_LONGRUN_1200s
mkdir -p "$OUT/out"
cd "$FR_STIR"
rm -rf output && mkdir -p output

./test_3d_ophelie_french_stirring_em \
  --reload-dir="$RELOAD_FR_STIR" \
  --rotation-rpm=10 \
  --thermo-em-coupling=periodic --em-update-interval=0.5 \
  --em-control=fixed-current --aind=off \
  --skull-resistance --skull-delta=0.005 \
  --free-crust-resistance --free-crust-delta=0.005 --free-crust-k=1 \
  --end-time=1200 --max-wall-hours=10 \
  --state-record-interval=30 \
  --summary-json="$OUT/summary.json" \
  2>&1 | tee "$OUT/run.log"

cp -f output/*.csv "$OUT/" 2>/dev/null || true
cp -a output/*.vtp "$OUT/out/" 2>/dev/null || true
```

备选（旧热闭合，**不要**与 skull 同开）：

```bash
./test_3d_ophelie_french_stirring_em \
  --reload-dir="$RELOAD_FR_STIR" \
  --balance-heat-loss \
  --end-time=1200 --max-wall-hours=10
```

---

## 3. 甲方自然对流（长时）

**身份**：`--client-nc-geometry --client-nc-material`（σ=16 借用）；平坦滑移顶；**ASSUMED** 玻璃吸收 **15 kW**（电源确认 30 kW，效率未确认 → 勿抄法国 50/60 kW）。  
长时：无 doc227；`--doc224-u-blowup-gate=50`；`end-time=600`。

```bash
OUT=$OUT_ROOT/03_client_nc/CLIENT_NC_ASSUMED15kW_LONGRUN_600s
mkdir -p "$OUT/out"
cd "$FROZEN"
rm -rf output output_backup && mkdir -p output

./test_3d_ophelie_french_natural_convection_frozen_q \
  --client-nc-geometry --client-nc-material \
  --dp=0.010 --melt-height=0.180 \
  --reload-dir="$RELOAD_CL_NC" \
  --doc223-paper-flat-slip-top \
  --em-control=fixed-power --target-absorbed-power=15000 \
  --thermo-em-coupling=periodic --em-update-interval=0.5 --aind=off \
  --enable-thermal-diffusion \
  --doc224-u-blowup-gate=50 \
  --end-time=600.0 --hydro-settle-time=1.0 \
  --state-record-interval=30.0 \
  --meridional-average-start=60.0 \
  --boundary-face-t-mode=local_dt_dz \
  --summary-json="$OUT/summary.json" \
  2>&1 | tee "$OUT/run.log"

cp -f output/french_energy_budget.csv \
      output/french_thermo_em_coupling.csv \
      output/french_meridional_timeseries.csv \
      output/french_T_rz.csv \
      "$OUT/" 2>/dev/null || true
cp -a output/*.vtp "$OUT/out/" 2>/dev/null || true
```

---

## 4. 甲方搅拌（工程挂机 · 弹性桨 F1）

**身份**：`client_stirring_f1_heated_elastic`；Reload 含 `GlassBody`+`Rotor`；桨 STL 在 bin/data；转速由 case（确认 **10 rpm**）；热源为 **规定体积热** `--heat-power=`（此处 ASSUMED **15000 W**，与甲 NC 对齐，**非** EM）。  
**成熟度**：仓库尚无与 `french_stirring_em` 对等的甲搅拌生产长跑身份；F0/F1 历史多为短窗诊断且 ρ/U 质量常 FAIL。本命令仅工程挂机，**不授验收**。  
墙钟：`--max-wall-hours=10`；物理：`--end-time=600`。

```bash
OUT=$OUT_ROOT/04_client_stir/CLIENT_STIR_F1_ASSUMED15kW_LONGRUN_600s
mkdir -p "$OUT"
cd "$CL_STIR"
rm -rf output && mkdir -p output

./test_3d_ophelie_client_stirring_f1_heated_elastic \
  --reload-dir="$RELOAD_CL_STIR" \
  --rotor-stl=./data/client_stir_paddle_solver_m.stl \
  --end-time=600 --omega-ramp=0.05 --c0=8 \
  --heat-power=15000 --clamp-top-band=0.05 \
  --max-wall-hours=10 \
  --out-dir="$OUT" --summary-json=summary.json \
  2>&1 | tee "$OUT/run.log"
```

若只要机械响应、热功率暂用开发默认：把 `--heat-power=15000` 改成 `--heat-power=1e5`（演示值，勿当甲方实测）。

---

## 建议挂机顺序

1. 法 NC 600 s（最贵、最易数值中止）  
2. 甲 NC 600 s  
3. 法搅拌 1200 s（自带墙钟退出）  
4. 甲搅拌 600 s  

四套可并行若有多卡；单卡请串行。跑完把对应 `run.log` / `summary.json` 回传即可核对。
