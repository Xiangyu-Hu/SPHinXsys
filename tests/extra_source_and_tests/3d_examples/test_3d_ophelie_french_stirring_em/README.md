# test_3d_ophelie_french_stirring_em

法国冷坩埚玻璃熔炉：感应加热 + 机械搅拌 + 自然对流的耦合算例。

对应文献：Jacoutot et al., *Chem. Eng. Process.* 47 (2008) 449-455（搅拌工况），
几何与感应器参数取 EREBUS 装置，电导率取 Jacoutot 博士论文（HAL tel-00121350）图 IV.10。

## 计算流程

**Phase A — 电磁预解（一次）**

在 relax 出来的玻璃粒子上求解 EM（complex edge-flux + self-induction Picard），
把线圈电流标定到 60 kW 吸收功率，然后把 Joule 功率密度 `Q [W/m^3]`
用 CIC 沉积到一张**固定的欧拉网格**上（默认 `h = dp = 5 mm`，105×105×47）。

**Phase B — 流动 + 传热（长时间推进）**

WCSPH 熔体 + Simbody 驱动的搅拌桨 + Boussinesq 浮力 + 法国 Robin/辐射边界条件。
每个 advection step 从欧拉网格上**重新插值** Q 到当前粒子位置。

用欧拉网格而不是像 `test_3d_ophelie_french_natural_convection_frozen_q`
那样把 Q 冻结在粒子上，是搅拌工况的关键：感应器固定在实验室坐标系里，
而熔体被搅拌桨输运着穿过这个空间固定的热源场。冻结在粒子上会让热源跟着流体一起转。

## 前置条件

先跑 relax 生成三个体的 `Reload.xml`：

```bash
cd build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_glass_relax/bin
./test_3d_ophelie_french_stirring_glass_relax --bodies=all
```

得到 `GlassBody`（玻璃减去搅拌桨，322396 个粒子）、`Rotor`（搅拌桨，2192 个）、
`WallBoundary`（坩埚壳，91976 个）。

## 参数

| 量 | 值 | 来源 |
|---|---|---|
| 熔体 | z 向圆柱 R = 0.25 m, H = 0.21 m | `data/glass-z.stl` |
| ρ, μ, cp, k, β | 2800 kg/m³, 4 Pa·s, 1150 J/(kg·K), 4 W/(m·K), 3e-5 1/K | CEP 2008 Table 1 + 搅拌工况 |
| T₀ | 1473 K | 同上 |
| σ(1473 K) | 24.75 S/m，`log10σ = 4.05726 − 3923.73/T` | 论文图 IV.10 |
| 感应器 | 7 匝, R = 0.285 m, 叠高 0.230 m（与熔体同心）, 300 kHz | EREBUS |
| 吸收功率 | 60 kW | 搅拌工况 |
| 搅拌 | 绕 z 轴 10 rpm（1.047 rad/s） | 搅拌工况 |
| h_side / h_bottom / h_free | 300 / 35 / 20 W/(m²·K) | 法国生产工况约定 |
| ε, T_cool / T_amb | 0.8, 300 K | 同上 |

`dp = 5 mm`，`c0 = 8 m/s`。

c0 的上限是稳定性给的，不是精度给的：玻璃、搅拌桨、坩埚壁是分别 relax 出来的，
接触面上有轻微重叠，由此产生的压力尖峰按 c0² 增长。实测 c0 = 10 会在 50 个声学步内
发散到 U_max ~ 3e7，c0 = 8 可以稳定长跑。代价是 0.21 m 的熔体柱有约 3% 的静水压缩，
这部分由自由面跟踪吸收（见下）。

## 已知的建模简化

**热平衡不闭合。** 按文献系数原样计算，1473 K 时边界总散热约 149 kW，
而感应输入只有 60 kW，其中自由面辐射就占 43 kW。真实冷坩埚在水冷壁上会长出
一层凝固的玻璃壳（skull），自由面上会结壳，两者都起隔热作用，本模型没有。
不加处理的话熔体会以约 0.6 K/s 的速率净冷却，1200 s 之后凉几百 K，没有意义。

- `--balance-heat-loss`：把 h_side / h_bottom / h_free / ε 统一乘一个系数
  （实测 0.4028），使 T₀ 时刻总散热正好等于 60 kW。相对分配比例仍按文献，
  缩放后 h_side = 121 W/(m²·K)，仍落在 docs/ophelie 记录的 [100, 400] 敏感性区间内。
  长时间生产算例推荐用这个。实测 8 s 物理时间内总散热 59.5 kW，T_mean 只漂 0.16 K。
- `--t-min=`（默认 300 K）：温度下限。边界壳层在没有 skull 的情况下会一直被抽热，
  不设下限会跑到 0 K 以下。

**自由面标记跟踪液面。** `setupOphelieThermalFrenchNaturalBoundaryFaces` 是按到名义
圆柱顶面的距离来标记辐射面的。WCSPH 是可压缩的，熔体在重力下会沉几个 dp，搅拌还会晃，
固定平面在 1 秒内就会丢掉大部分辐射粒子——而辐射占了热收支的三分之二。
本算例每次重标记前先用粒子 z 坐标的 99.5 分位数量出当前液面，再传给标记函数。

**没有 SPH 热传导。** 与 `frozen_q` 生产算例一致，只有 Q 源项 + 边界条件。
熔体热扩散系数 α = 1.24e-6 m²/s，在几百秒的物理时间内扩散长度只有几个 dp。

**搅拌桨绝热。** 文献没有给桨的热物性，桨不参与传热。

**边界面标记会随流动失效。** Robin/辐射边界是按粒子到名义圆柱面的距离标记的，
搅拌会把粒子搬走。默认每 50 个 advection step 重新标记一次（`--bc-retag-every=`）。

## 运行

```bash
cd build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_em/bin
./test_3d_ophelie_french_stirring_em \
    --reload-dir=../../test_3d_ophelie_french_stirring_glass_relax/bin/reload \
    --balance-heat-loss
```

Phase A 约需 3 分钟（EM self-induction Picard + 坩埚壁 level set）。Phase B 按
`--max-wall-hours=10` 的预算跑，到点干净退出并写最后一帧。

实测推进速率约 0.039 s 物理时间 / s 墙钟（322396 个流体粒子），
10 小时约合 1400 s 物理时间 ≈ 230 转，所以默认 `--end-time=1200`（200 转）
会在预算内自然跑完。

常用选项：

| 选项 | 默认 | 说明 |
|---|---|---|
| `--end-time=` | 1200 s | 物理时间上限（200 转） |
| `--max-wall-hours=` | 10 | 墙钟预算，到点停 |
| `--c0=` | 8 | 人工声速，超过 9 会发散 |
| `--balance-heat-loss` | 关 | 缩放边界系数使 T₀ 热平衡 |
| `--h-side=` `--h-bottom=` `--h-free=` | 文献值 | 单独覆盖对流系数 |
| `--t-min=` `--t-max=` | 300 / 2500 K | 温度钳制 |
| `--q-grid-spacing=` | dp | 欧拉 Q 网格间距 |
| `--bc-retag-every=` | 50 | 边界面重标记间隔（advection step） |
| `--state-record-interval=` | 10 s | VTP 输出间隔 |
| `--no-state-recording` | — | 关闭 VTP |
| `--rotation-rpm=` | 10 | 搅拌转速 |
| `--no-boussinesq` | — | 关闭浮力（纯机械搅拌） |
| `--no-self-induction` | — | EM 不做 Picard 自感迭代 |

## 输出

- `output/StirringJouleHeatGrid.vti` — 欧拉 Q 场，可在 ParaView 里直接看感应加热分布
- `output/GlassBody_*.vtp` — 熔体粒子（Pressure / Temperature / JouleHeat）
- `output/RotorProxy_*.vtp` — 搅拌桨表面网格
- `output/french_stirring_em_monitor.csv` — 时间、转数、U_max、T_mean、T_max、墙钟

`--state-record-interval=10` 时约输出 120 帧，每帧 48 MB，合计约 5.8 GB。
10 rpm 下一转是 6 s，想看清混合过程可以调到 `--state-record-interval=3`（约 23 GB，先看磁盘）。

结束时打印 `passed=1` 需要同时满足：EM 功率标定误差 < 1%、
Q 网格重采样相对 L2 < 0.2 且没有粒子掉出网格、没有发散、速度有限且非零、边界散热为正。
主循环里有发散保护：声学 dt 变成非有限值或塌缩到 1e-9 以下就带诊断信息退出，
不会像之前那样在 NaN 物理时间上空转。
