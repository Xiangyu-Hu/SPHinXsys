# OPHELIE French stirring — CAD meshes

CMake 会把本目录整个拷到 `bin/input/`，所以下面的路径在运行时都是 `./input/<name>`。

| 文件 | 作用 |
|---|---|
| `glass-z.stl` | **运行用。** 熔体，z 圆柱 R = 250 mm，z ∈ [0, 210] mm |
| `stirring_paddle_2_z.stl` | **运行用。** 搅拌桨实体，轴在 z 上 |
| `glass-y.stl` | CAD 原件（y-up），只作为转换脚本输入 |
| `stirring_paddle_2.stl` | CAD 原件（y-up），同上 |
| `rotate_cad_to_z_up.py` | y-up → z-up 转换脚本 |

单位 **mm**，运行时用 `--geometry-scale=0.001` 换算成米（默认值）。不做自动平移，
STL 按导出坐标直接使用。

## 为什么要转坐标

CAD 交付时熔体轴沿 -y，而共享的法国 EM/热工具（线圈叠放、圆柱 level set、
Robin/辐射面标记、phi 边界法向）全部按 +z 圆柱写死。脚本做一次 +90° 绕 x 轴旋转，
并把熔体底面平移到 z = 0。几何改了就重跑：

```bash
python3 rotate_cad_to_z_up.py
```

生成的 z-up 文件要同步拷贝到 `test_3d_ophelie_french_stirring_em/data/`。

## relax 之后

`reload/Reload.xml` 包含三个体：

- `GlassBody` — `glass-z` 减去 `stirring_paddle_2_z`
- `Rotor` — 搅拌桨实体
- `WallBoundary` — 解析生成的坩埚壳（不来自 STL）
