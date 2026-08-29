# test_3d_ophelie_french_stirring_glass_relax

法国搅拌工况的粒子分布生成：玻璃熔体、搅拌桨、坩埚壁三个体的 relax。

输出的 `reload/Reload.xml` 是 `test_3d_ophelie_french_stirring_em` 的前置输入。

## 几何

CAD 交付时熔体轴沿 **-y**，但所有共享的法国工具（线圈叠放、圆柱 level set、
Robin/辐射面标记、phi 边界法向）都假设 **+z** 圆柱。`data/rotate_cad_to_z_up.py`
把网格一次性转到 z-up 并把熔体底面移到 z = 0：

| 文件 | 说明 |
|---|---|
| `glass-y.stl` / `stirring_paddle_2.stl` | CAD 原件（y-up，mm），只作为转换脚本的输入 |
| `glass-z.stl` | **运行用。** 熔体，z 圆柱 R = 250 mm，z ∈ [0, 210] mm |
| `stirring_paddle_2_z.stl` | **运行用。** 非对称搅拌桨，轴在 z 上，最薄叶片约 8 mm |

`--geometry-scale=0.001`（mm → m），无平移，STL 按导出坐标直接用。

坩埚壁不是 STL，是由 `OphelieFrenchNaturalCrucibleWallShape` 按熔体包围盒解析生成的
壳体（侧环 + 底板，顶部敞开），厚度 `wall_thickness_factor × dp`。

## 分辨率与 level set

搅拌桨叶片比坩埚薄得多，而且两个体都要分辨它：转子**就是**叶片，玻璃则要挖出一个
叶片形状的空腔。由此定下两个设置。

- `dp = 0.005` — `dp = 0.010` 时一个粒子就跨过整片叶片，转子的 lattice 会塌成约 96 个粒子。
- 两个体的 level-set refinement 都取 `5` — 参照 sphinxsys-einsimo `test_3d_stirring`，
  `ratio = dp / (0.1 × 叶片厚度) = 0.005 / 0.0008 ≈ 6`。取 1 时 level-set 网格等于 `dp`，
  `cleanLevelSet()` 会把空腔抹掉，玻璃出来是个实心圆柱。

`dp = 0.005` 下的实测粒子数：`Rotor` 2192、`GlassBody` 322396、`WallBoundary` 91976。

## 运行

```bash
cd build
ninja test_3d_ophelie_french_stirring_glass_relax

cd tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_glass_relax/bin
./test_3d_ophelie_french_stirring_glass_relax --bodies=all
```

约需 8 分钟，其中大部分是玻璃体 level set 的 `UpdateKernelIntegrals`。

输出：

- `./reload/Reload.xml` — `GlassBody` + `Rotor` + `WallBoundary`（约 46 MB）
- `./output/*_ite_*.vtp` — 在 ParaView 里检查粒子分布

`Reload.xml` 是所有体共用的一个文件，每次运行都从头重写，所以
`--bodies=rotor` 之类的会产出流动算例用不了的 reload。这些模式只用于单独检查某个体，
正式跑之前必须 `--bodies=all`。

## CLI

- `--bodies=all|glass|rotor|wall`（逗号分隔的子集也可以）
- `--glass-stl=./input/glass-z.stl`
- `--rotor-stl=./input/stirring_paddle_2_z.stl`
- `--geometry-scale=0.001`
- `--dp=0.005`
- `--rotor-level-set-refinement=5`
- `--glass-level-set-refinement=5`
- `--wall-thickness-factor=4`
- `--wall-mesh-resolution=40`
- `--relax-steps=1000`
- `--relax-vtp-every=100`
