# test_3d_ophelie_french_natural_glass_relax

法国 natural-convection 玻璃圆柱：几何 → SYCL level-set → lattice 粒子 → SYCL-CK relax。

几何：

- `TriangleMeshShapeCylinder`
- `R = 0.250 m`
- `H = 0.185 m`
- 轴线 `+z`，`z ∈ [0, H]`
- `dp = 0.015 m`

实现对齐：

- `tests/tests_sycl/3d_examples/test_3d_particle_relaxation_single_resolution_sycl`
- `tests/tests_sycl/3d_examples/test_3d_taylor_bar_sycl`

关键步骤：

```text
defineBodyLevelSetShape(par_ck).correctLevelSetSign().writeLevelSet()
→ Lattice
→ RandomizeParticlePositionCK
→ SYCL-CK loop:
    cell linked list
    update inner relation
    KernelGradientIntegral + LevelsetKernelGradientIntegral
    RelaxationScalingCK
    PositionRelaxationCK
    LevelsetBounding
→ ReloadParticleIOCK
```

## 编译与运行

需要 `SPHINXSYS_USE_SYCL=ON`。

```bash
cd ~/SPHinXsysSYCL/build
cmake --build . --target test_3d_ophelie_french_natural_glass_relax -j$(nproc)

cd tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_glass_relax/bin
./test_3d_ophelie_french_natural_glass_relax

# 冒烟
./test_3d_ophelie_french_natural_glass_relax --dp=0.02 --relax-steps=200 --relax-vtp-every=50

# 正式
./test_3d_ophelie_french_natural_glass_relax --dp=0.015 --relax-steps=1000
```

## CLI

| Flag | 默认 | 含义 |
|------|------|------|
| `--dp=` | `0.015` | 粒子间距 |
| `--glass-radius=` | `0.25` | 玻璃半径 |
| `--glass-height=` | `0.185` | 玻璃高度 |
| `--mesh-resolution=` | `20` | `TriangleMeshShapeCylinder` 网格分辨率 |
| `--relax-steps=` | `1000` | relax 步数 |
| `--relax-vtp-every=` | `100` | VTP 间隔；`0` 关闭 |
| `--state_recording=` | SPHSystem | 透传给 `handleCommandlineOptions` |

## 输出

- `./output/GlassBody_ite_*.vtp`
- `./reload/Reload.xml`
