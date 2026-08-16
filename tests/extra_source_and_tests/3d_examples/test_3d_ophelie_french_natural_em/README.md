# test_3d_ophelie_french_natural_em

法国 natural-convection 玻璃圆柱的 **第一阶段电磁**：complex edge-flux，无 `A_glass` 自感，无热/流动。

默认几何 / 激励（期刊 + EREBUS 补充）：

| 量 | 默认 | 备注 |
|---|---|---|
| `R` | 0.250 m | 期刊 / EREBUS fill |
| `H` | 0.185 m | 同上 |
| `f` | 282 kHz | EREBUS 体系补充 |
| `N_turn` | 7 | 博士论文/装置补充，非期刊正文硬数值 |
| `coil_R` | 0.285 m | Stage1.5 coil geometry (EREBUS inductor) |
| `sigma` | 16 S/m | Table 1 @ 1473 K |
| `target_P` | 50 kW | target-power 模式 |

## 推荐流程

先跑 relax case 得到 `Reload.xml`，再跑本 case：

```bash
cd ~/SPHinXsysSYCL/build
cmake --build . --target test_3d_ophelie_french_natural_glass_relax \
  test_3d_ophelie_french_natural_em -j$(nproc)

# 1) relax
cd tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_glass_relax/bin
./test_3d_ophelie_french_natural_glass_relax --dp=0.015 --relax-steps=1000

# 2) EM（fixed current，验证用）
RELAX_RELOAD=$PWD/reload
cd ../../test_3d_ophelie_french_natural_em/bin
./test_3d_ophelie_french_natural_em \
  --reload=1 --reload-dir="$RELAX_RELOAD" \
  --excitation-mode=current --coil-current=1.0

# 3) EM（标定到 50 kW）
./test_3d_ophelie_french_natural_em \
  --reload=1 --reload-dir="$RELAX_RELOAD" \
  --excitation-mode=target-power --target-power=50000
```

也可把 `Reload.xml` 拷到本 case 的 `./reload/`，然后省略 `--reload-dir`。

## CLI

本 case：

| Flag | 含义 |
|------|------|
| `--excitation-mode=current` | 固定 `I_per_loop`（默认） |
| `--excitation-mode=target-power` | 标定线圈电流使 `P≈target` |
| `--reload=1` / `--reload-dir=` | 读取 relax 粒子 |
| `--lattice` / `--reload=0` | 不读 reload（仅调试） |

复用 French / OPHELIE 通用选项：

| Flag | 含义 |
|------|------|
| `--dp=` `--glass-radius=` `--glass-height=` | 几何 |
| `--coil-radius=` `--coil-num-loops=` `--coil-current=` | 线圈 |
| `--coil-segments-per-loop=` `--ampere-turns=` | 线圈细分 / 安匝 |
| `--frequency=` `--sigma=` `--target-power=` | EM 物性 / 目标功率 |

## 输出

- console：`P_raw`、`P_recon_edge`、`I_peak` / `I_rms`、φ residual、edge residual reduction
- `./output/GlassBody_ite_0000000000.vtp`：`A`、`E/J` edge recon、`Q`、`phi`

## 本阶段明确不做

- `A_glass` Picard 自感
- `sigma(T)` / 热边界 / 自然对流
- 冷坩埚 / 底板金属表面电流
