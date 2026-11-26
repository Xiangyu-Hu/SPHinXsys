# Doc248 甲方 NC（C0）+ 搅拌（C2 / C2r）交接包

与 Doc247 法国线隔离。粒子**一律本机重新 relax 生成**，不附带现成 Reload。

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

# 3) C0 自然对流（STANDARD：H=230、GPT IC、18 kW、κ=3.5、coil 68–97%H、
#    H_SIDE=70 H_FREE=10 EMISSIVITY=0.4；默认 end=60 s）
# 短烟：END_TIME=10 bash .../02_run_C0_nc.sh
bash docs/ophelie/reference/doc248_client_cases/02_run_C0_nc.sh

# 4a) C2 弹性搅拌（默认 END=5s、ASSUMED soft E=2e8；不锈钢 E=2e11 会在 ~0.3s tip 跑飞）
# Rotor VTP 含 VonMisesStress / VonMisesStrain / GreenLagrangeStrain / DeformationGradient
END_TIME=5 MAX_WALL_HOURS=24 bash docs/ophelie/reference/doc248_client_cases/03_run_C2_stir.sh
# 更长：END_TIME=60 …；若要硬钢短烟：YOUNGS_E=2e11 OMEGA_RAMP=0.5 END_TIME=0.5 …

# 4b) C2r 刚体搅拌（冲 ~1400 K 的生产长跑；推荐）
END_TIME=7200 MAX_WALL_HOURS=48 TARGET_T_K=1400 \
  bash docs/ophelie/reference/doc248_client_cases/04_run_C2r_rigid_stir.sh
```

输出在 `runs/`（gitignore）。C2 弹性固体子步多，墙钟会很长；**升温长跑请用 C2r 刚体桨**。

C2 / C2r 默认挂 **ASSUMED** 初场（法国 Fig.8a T + Fig.10 Q 映射到甲方 R/H）：

- `docs/ophelie/reference/standard_data/doc248_seeds/SEED_CLIENT_STIR_Fig8a_T_rz.csv`
- `docs/ophelie/reference/standard_data/doc248_seeds/SEED_CLIENT_STIR_Fig10_Q_rz.csv`

可用 `INIT_T_CSV=` / `INIT_Q_CSV=` 覆盖。重编后才生效：

```bash
cmake --build "$BUILD" --target test_3d_ophelie_client_stirring_f1_heated_elastic \
  test_3d_ophelie_client_stirring_rigid_heated -j
```

## 目标二进制

| Case | Relax | Run |
|------|-------|-----|
| C0 NC | `test_3d_ophelie_client_nc_glass_relax` | `frozen_q --client-nc-*` |
| C2 搅拌（弹性） | `test_3d_ophelie_client_stirring_glass_relax --bodies=all` | `client_stirring_f1_heated_elastic` |
| **C2r 搅拌（刚体）** | 同上 | **`client_stirring_rigid_heated`** |
