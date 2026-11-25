# Doc248 甲方 NC（C0）+ 搅拌（C2）交接包

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

# 2) 从零生成粒子（NC + 搅拌，会删旧 reload）
bash docs/ophelie/reference/doc248_client_cases/01_generate_particles.sh

# 3) C0 自然对流（默认 60 s；可用 END_TIME= 改）
bash docs/ophelie/reference/doc248_client_cases/02_run_C0_nc.sh

# 4) C2 弹性搅拌（慢；默认 end=60 s / max-wall=24 h）
END_TIME=600 MAX_WALL_HOURS=48 bash docs/ophelie/reference/doc248_client_cases/03_run_C2_stir.sh
```

输出在 `runs/`（gitignore）。C2 弹性固体子步多，墙钟会很长，属预期。

## 目标二进制

| Case | Relax | Run |
|------|-------|-----|
| C0 NC | `test_3d_ophelie_client_nc_glass_relax` | `frozen_q --client-nc-*` |
| C2 搅拌 | `test_3d_ophelie_client_stirring_glass_relax --bodies=all` | `client_stirring_f1_heated_elastic` |
