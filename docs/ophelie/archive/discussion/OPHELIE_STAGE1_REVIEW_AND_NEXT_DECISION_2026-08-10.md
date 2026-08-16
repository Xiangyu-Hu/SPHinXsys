# OPHELIE French Natural Stage 1 代码审查与下一步决策
**日期：2026-08-10**

## 0. 最终决策

本轮工作整体方向正确，**旧的 edge-flux 功率灾难性错误 B1 可以认为已经解决**：

- 旧 `P_graph/P_recon ~ 1e5` 已被定位为 pointwise SPH Laplace pair coefficient 被错误当成 Watt conductance，缺少粒子体积语义；
- Phase A 后 French case 变为 `P_graph/P_recon ≈ 0.355`；
- `P_undirected/P_graph = 1`；
- uniform-field MMS 中 `P_recon/P_exact = 1`，`P_graph/P_exact ≈ 0.783`；
- complex real/imag 两条链均通过。

**但我不建议现在直接宣布 Stage 1 完全封版并立即进入大规模 A_glass。**
原因不是 `0.355` 本身，而是审查中发现了两个比 `0.355` 更重要的配置/一致性问题：

1. **当前 French run 的粒子数与声明的 `dp=0.01` 明显不一致，疑似 reload 实际由约 `dp=0.015` 生成，却用 `dp=0.01` 的 SPHSystem / smoothing length 读取。**
2. **线圈几何仍是 provisional：代码 `coil_R=0.300 m`，而 EREBUS 补充参数为直径 `570 mm`，即 `R=0.285 m`；当前 auto-z/inset 方案也没有体现 `230 mm` 的线圈总高度。**

因此建议增加一个非常短的 **Stage 1.5 Closeout**。把上述配置问题和弱形式 edge-current 诊断修正后，**不再追求强行让 `P_graph/P_recon → 1`，直接进入 A_glass Picard。**

---

# 1. 本轮代码中符合预期的内容

## 1.1 French cylinder case 已建立

`test_3d_ophelie_french_natural_glass_relax`：

- `R=0.250 m`
- `H=0.185 m`
- 使用 `TriangleMeshShapeCylinder`
- `LevelSetShape`
- lattice particles
- `RandomizeParticlePositionCK`
- SYCL/CK relaxation
- `ReloadParticleIOCK`

这个方向正确，符合我们此前“不需要 STL，直接 SPHinXsys 内建几何 + relax + reload”的决策。

## 1.2 French EM case 已强制 complex edge-flux

`test_3d_ophelie_french_natural_em.cpp` 明确：

```cpp
params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
params.edge_flux_complex_ = true;
params.enable_self_induction_ = false;
```

因此这个新 French Natural EM case 不会因为 generic `literature-mode` 默认 DivGrad 而走错 operator。

## 1.3 fixed-current 与 target-power 两种激励模式已经存在

这是正确设计。

- fixed-current：用于电磁真实性、标度律和 A_glass 验证；
- target-power：用于 French 50 kW / 60 kW 论文复现。

不要删掉 fixed-current 模式。

## 1.4 Phase A 的 volume semantics 修正是正确方向

当前：

\[
C_{ij}=\sigma_{ij}\alpha_{ij}V_j
\]

本质仍是 pointwise SPH Laplacian coefficient，而不是 Siemens。

因此：

\[
\sum_j \frac14 C_{ij}|e_{ij}|^2
\]

应首先解释为局部离散能量密度量，再乘 \(V_i\) 才得到 particle Watt contribution。

当前实现：

\[
P_{\rm graph}
=
\sum_i
V_i
\sum_j
\frac14 C_{ij}|e_{ij}|^2
\]

以及无向审计：

\[
P_{\rm undirected}
=
\sum_{j>i}
\frac12
\alpha_{ij}V_iV_j|e_{ij}|^2
\]

在数学上是相容的；因此得到：

\[
P_{\rm undirected}/P_{\rm graph}=1
\]

是合理结果。

旧的 `~1e5` 已不再是 blocker。

## 1.5 LS reconstructed E/J/Q 仍作为 production 主链是合理的

当前 production：

\[
E_{\rm recon}
\rightarrow
J_{\rm recon}=\sigma E_{\rm recon}
\rightarrow
Q_{\rm recon}
=
\frac12\sigma
\left(
|E_{\rm real}|^2+|E_{\rm imag}|^2
\right)
\]

target-power 使用：

\[
P_{\rm recon}
=
\sum_iQ_{{\rm recon},i}V_i
\]

进行标定。

目前不建议切换成 `P_graph` 标定。

---

# 2. 发现的首要问题：dp / reload 很可能不一致

当前 snapshot 写：

```text
R = 0.25 m
H = 0.185 m
dp = 0.01 m
n = 10474
```

解析圆柱体积：

\[
V
=
\pi(0.25)^2(0.185)
=
0.0363247\ {\rm m^3}
\]

若真的是 3D single-resolution lattice 且：

\[
dp=0.01\ {\rm m}
\]

名义粒子数数量级应约：

\[
N
\approx
\frac{0.0363247}{(0.01)^3}
\approx
36325
\]

而不是 10474。

反过来，由 10474 个粒子推算的 cubic-equivalent spacing：

\[
dp_{\rm eff}
=
\left(
\frac{V}{N}
\right)^{1/3}
\approx
0.01514\ {\rm m}
\]

而：

\[
\frac{V}{(0.015)^3}
\approx10763
\]

与 10474 非常接近。

**因此目前高度怀疑：当前 Reload.xml 实际来自约 `dp=0.015` 的 relax run，而 EM case 仍以 `french.dp=0.010` 创建 SPHSystem。**

这会导致：

- particle coordinate spacing
- `VolumetricMeasure`
- reference smoothing length
- cell linked-list/support radius
- pairwise kernel derivative
- boundary width

之间不一致。

这甚至可能是 French case：

\[
P_{\rm graph}/P_{\rm recon}=0.355
\]

比 uniform MMS 的 0.783 更低的重要原因之一。

## 必须立即增加 reload consistency audit

EM case 读取 reload 后，在真正做 EM 之前打印：

```text
declared_dp
n_particles
sum_particle_volume
analytic_volume
volume_rel_error
mean_particle_volume
dp_from_mean_volume = cbrt(mean_particle_volume)
nearest_neighbor_mean
nearest_neighbor_p50
nearest_neighbor_min/max
reference_smoothing_length
mean_neighbor_count
min/max_neighbor_count
```

至少 gate：

\[
\left|
\frac{\sum_iV_i-V_{\rm exact}}
{V_{\rm exact}}
\right| < 3\%
\]

以及：

\[
|dp_{\rm eff}-dp_{\rm declared}|/dp_{\rm declared}<5\%
\]

如果当前 reload 确实是 `dp=0.015`，有两种选择：

### 快速开发方案（推荐）

先正式把开发 baseline 设置为：

```text
dp = 0.015
```

重新用同样 `dp=0.015` 建立 SPHSystem 并读对应 reload。

优点：

- n ~ 1e4；
- A_glass 的 O(N²) Picard 成本更可控；
- 最快推进全耦合。

之后再做：

```text
dp = 0.0125
dp = 0.010
```

收敛检查。

### 或正式重生成 dp=0.010 reload

这样粒子数应接近 3.5–3.6 万，A_glass O(N²) 成本会显著增加。

当前“尽快给甲方出全耦合结果”的目标下，我建议先用 **dp=0.015 作为开发级 baseline**，不要假装它是 0.010。

---

# 3. 第二个必须修正的问题：线圈几何

当前代码：

```cpp
french.coil.num_loops = 7;
french.coil.loop_radius = 0.300; // provisional
french.coil.segments_per_loop = 64;
french.coil.use_cell_centered_loops = true;
french.coil_z_inset = 0.02;
french.auto_coil_z = true;
```

但是 EREBUS 补充参数：

```text
Inductor diameter = 570 mm
number of turns = 7
total height = 230 mm
frequency = 282 kHz
```

因此 filament 等效圆环半径应首先改为：

\[
R_{\rm coil}=0.285\ {\rm m}
\]

而不是 0.300 m。

更重要的是：**线圈 stack 总高度应该体现 230 mm，而不是把 7 匝自动塞进 185 mm 玻璃高度内部。**

Natural glass：

\[
H_{\rm glass}=185\ {\rm mm}
\]

而线圈总高度：

\[
H_{\rm coil}=230\ {\rm mm}
\]

所以线圈必然在轴向上超出液位范围。

## 推荐实现

增加明确参数：

```text
coil_radius = 0.285
coil_num_loops = 7
coil_stack_height = 0.230
coil_center_z = configurable
```

若当前资料还找不到精确每匝 z 坐标，则第一版：

- 7 匝等间距；
- 将 230 mm stack 的中心与 glass bath 中心对齐；
- 文档明确标记：`equal-spacing / center-aligned supplemental assumption`。

必须打印全部：

```text
loop 0: R, z
loop 1: R, z
...
loop 6: R, z
```

不要只打印 `N_turn/coil_R`。

---

# 4. 对 Cursor 提出的 “Phase B 对称 G” 的决定

**现在不要直接把 `C_ij` 强制改成 `C_ij=C_ji`。**

原因：

当前 pointwise coefficient：

\[
C_{ij}
=
\alpha_{ij}V_j
\]

其中：

\[
\alpha_{ij}
=
\sigma_{ij}\,w_{ij}
\]

在 kernel / harmonic-mean 对称的情况下：

\[
\alpha_{ij}=\alpha_{ji}
\]

所以本来就有：

\[
C_{ij}\ne C_{ji}
\]

当 \(V_i\ne V_j\) 时这是正常的，因为它是 **pointwise Laplace row coefficient**。

真正有 conservation/current units 意义的 weak edge current 应定义为：

\[
\boxed{
I_{ij}^{w}
=
V_iq_{ij}
=
-V_iC_{ij}e_{ij}
=
-\alpha_{ij}V_iV_je_{ij}
}
\]

反向：

\[
I_{ji}^{w}
=
-\alpha_{ji}V_jV_ie_{ji}
=
-I_{ij}^{w}
\]

因此正确的 antisymmetry 应检查：

\[
\boxed{
I_{ij}^{w}+I_{ji}^{w}\approx0
}
\]

而不是当前的：

\[
q_{ij}+q_{ji}\approx0
\]

## 当前 q_antisym diagnostic 应修改

当前代码比较：

```cpp
q_ij = -c_ij * edge_drop_ij;
q_ji = -c_ji * edge_drop_ji;
antisym = q_ij + q_ji;
```

对于不同粒子体积，这不是严格正确的守恒量。

改成：

```cpp
I_ij = Vol_i * q_ij;
I_ji = Vol_j * q_ji;
antisym = I_ij + I_ji;
```

或直接构造：

\[
G_{ij}^{w}
=
\alpha_{ij}V_iV_j
\]

然后：

\[
I_{ij}^{w}=-G_{ij}^{w}e_{ij}
\]

## 是否以后可以用 symmetric weak operator？

可以。

因为：

\[
V_iR_i=0
\]

与：

\[
R_i=0
\]

对 \(V_i>0\) 是等价的。

所以将来若想得到真正对称的：

\[
B^TGB
\]

并尝试 PCG，可以用 weak row-scaled system：

\[
G_{ij}^{w}=\alpha_{ij}V_iV_j
\]

但这是后续 solver optimization，不应阻塞 A_glass。

---

# 5. `P_graph/P_recon = 0.355` 是否阻塞 Stage 2？

**我的决定：不把它作为硬 blocker。**

原因：

`P_graph` 与 `P_recon` 目前不是同一个离散 functional。

### Graph

\[
P_{\rm graph}
=
\sum_e
\frac12G_e^w|e_e|^2
\]

是 weak SPH discrete energy。

### Reconstruction

LS reconstruction 先从 edge directional drops 重构：

\[
E_i
\]

再算：

\[
P_{\rm recon}
=
\sum_i
\frac12\sigma_i|E_i|^2V_i
\]

它是 continuum-field estimator。

只要 reconstruction 不是严格由同一个 variational functional 得到，就不能要求两者在有限粒子分辨率下严格相等。

uniform-field MMS 已经说明：

```text
P_recon/P_exact = 1
P_graph/P_exact = 0.78285
```

因此强行要求：

```text
P_graph/P_recon < 5% difference
```

目前并没有充分数学依据。

## 正确的下一步 gate

不要问：

> `P_graph` 是否等于 `P_recon` 到 1%？

应该问：

> 在 geometry/dp/coil 都一致以后，`P_recon`、J reconstruction 和 normalized Q distribution 是否随分辨率收敛？

所以：

- production calibration 继续使用 `P_recon`；
- thermal source 继续使用 `Q_recon`；
- A_glass 使用 `J_recon`；
- `P_graph` 保留为 weak-energy audit；
- `P_undirected/P_graph≈1` 继续作为离散 identity gate。

若 dp refinement 后 `P_graph/P_recon` 收敛到一个稳定的 0.4–0.8，而 `P_recon/Q` 自身收敛，不要再浪费主线时间强行把它压到 1。

---

# 6. 当前 `passed=1` 不能等价于 “Stage 1 fully accepted”

French case 的最终 gate 现在基本只有：

```cpp
fields_ok = finite
power_ok = target within 5%
passed = fields_ok && power_ok
```

这太宽松。

建议增加：

```text
--strict-stage1
```

至少检查：

- reload/dp consistency
- total-volume error
- phi solver residual
- edge residual reduction
- weak-current antisymmetry
- no nonfinite Q
- no negative reconstructed Q
- P_undirected/P_graph identity
- current scaling law
- coil segment convergence status

target-power 的 5% gate 也太宽，因为这里是解析 \(I^2\) calibration，建议：

\[
|P-P_{\rm target}|/P_{\rm target}<10^{-3}
\]

作为 numerical calibration gate。

---

# 7. uniform-field MMS gate 也建议收紧

当前：

```text
0.5 < P_recon/P_exact < 2
0.25 < P_graph/P_exact < 4
```

对 regression 太宽。

既然现在实际得到：

\[
P_{\rm recon}/P_{\rm exact}=1
\]

建议：

```text
abs(P_recon/P_exact - 1) < 1e-3
abs(P_undirected/P_graph - 1) < 1e-5
```

`P_graph/P_exact` 暂时只设稳定范围，不要求逼近 1。

---

# 8. Stage 1.5：下一轮 Cursor 应先完成的最小修补

## Task 1 — Reload / spacing audit 【必须】

1. 在 relax 输出中写入/打印：
   - dp
   - particle count
   - sum Vol
   - analytic volume
   - volume relative error
2. EM 读取 reload 后重新打印同样信息。
3. 计算 `dp_eff = cbrt(mean Vol)`。
4. 打印 neighbor-count。
5. 若 reload 与 requested dp 不一致，直接 fail，禁止继续 EM。

建议先确认当前 10474 粒子到底是不是 dp=0.015。

## Task 2 — 修正 EREBUS filament coil 【必须】

改：

```text
coil_R: 0.300 -> 0.285 m
```

增加：

```text
coil_stack_height = 0.230 m
```

记录 7 个 z。

若 exact axial offset 不可得，先 center-align，并标记补充假设。

## Task 3 — 改 weak-current continuity audit 【必须】

定义：

\[
I_{ij}^w=V_iq_{ij}
\]

新增：

```text
weak_edge_antisym_rel_l2
global_weak_kcl
```

不要再把 directed `q_ij+q_ji` 当成最终 current-conservation gate。

## Task 4 — Coil quadrature convergence 【很快】

fixed current, fixed sigma：

```text
segments = 32 / 64 / 128
```

比较：

- `P_recon`
- `Q_max`
- normalized `Q(r,z)`
- A norm

64 vs 128：

```text
P_recon relative difference < 1%
normalized-Q L2 difference < 1–2%
```

通过后固定 64 或 128。

## Task 5 — 最小 dp convergence 【必须，但不要过度】

为了快速推进：

```text
dp = 0.015
dp = 0.0125
```

先两档。

如成本允许再：

```text
dp = 0.010
```

全部使用各自匹配的 relax reload。

fixed-current `I=1 A/loop`, `sigma=16`, corrected coil geometry。

比较：

- `P_recon`
- `P_graph/P_recon`
- `Q_outer/Q_center`
- `Q_max`
- normalized Q map
- phi residual
- runtime

若 finest two：

```text
P_recon change < 5%
Q_outer/Q_center change < 5%
normalized Q shape stable
```

即可进入 A_glass。

**不要为了 `P_graph/P_recon=1` 卡住。**

---

# 9. Stage 2：随后立即进入 A_glass Picard

Stage 1.5 通过后，下一主任务就是玻璃 self-induction。

先固定：

```text
T = 1473 K
sigma = 16 S/m
excitation = fixed-current
base_em = off
crucible_em = off
coil_surface = off
```

求：

\[
A_{\rm total}=A_{\rm coil}+A_{\rm glass}
\]

\[
A_{\rm glass}(x_i)
=
\frac{\mu_0}{4\pi}
\sum_{j\ne i}
\frac{J_jV_j}
{\sqrt{|x_i-x_j|^2+\epsilon^2}}
\]

## Picard

\[
J^{n}
\rightarrow
A_{\rm glass}^{n}
\rightarrow
\phi^{n+1}
\rightarrow
J^{n+1}
\]

建议先对 `A_glass` 欠松弛：

\[
A_{\rm glass}^{n+1}
=
(1-\alpha)A_{\rm glass}^{n}
+
\alpha A_{\rm new}
\]

初始：

```text
alpha = 0.3–0.5
max_iter = 50
```

inner phi solve 在 A_glass 阶段建议收紧到至少约：

```text
phi solver relative residual <= 1e-6
```

否则 `1e-4` inner error 与 outer Picard target 同量级，不利于判断真正收敛。

## Picard convergence

同时记录：

\[
r_A=
\frac{\|A_g^{n+1}-A_g^n\|}
{\|A_{\rm coil}+A_g^{n+1}\|}
\]

\[
r_J=
\frac{\|J^{n+1}-J^n\|}
{\|J^{n+1}\|}
\]

\[
r_P=
\frac{|P^{n+1}-P^n|}
{P^{n+1}}
\]

第一版 gate：

```text
r_A < 1e-4
r_J < 1e-4
r_P < 1e-4
```

---

# 10. A_glass 必做验证

## 10.1 Low-conductivity limit

例如：

```text
sigma = 0.16 / 1.6 S/m
```

应看到：

\[
A_{\rm glass}/A_{\rm coil}\rightarrow0
\]

并且 self-induction 与 coil-only 结果接近。

## 10.2 正式 sigma=16 S/m

报告：

```text
||A_glass||/||A_coil||
max(|A_glass|/|A_coil|)
P_glass self-ind / coil-only
Q_outer/center change
Q_peak location change
```

## 10.3 Current scaling

固定 sigma：

\[
I_2=2I_1
\]

即使包含线性 self-induction，也应该近似：

\[
A_2=2A_1
\]

\[
J_2=2J_1
\]

\[
P_2=4P_1
\]

## 10.4 Relaxation-parameter independence

测试：

```text
alpha = 0.2 / 0.4 / 0.6
```

只允许 convergence speed 改变，最终解应相同。

## 10.5 性能

当前 direct Biot-Savart 是：

\[
O(N^2)
\]

所以必须记录：

```text
N
A_glass kernel ms
phi solve ms
one Picard cycle ms
total Picard time
GPU
```

这将决定后面 thermal transient 的 `EM update interval`。

---

# 11. Stage 3：A_glass 通过后立刻做 sigma(T)

当前 `electromagnetic_ophelie_french_material_laws.h` 的数据层设计是正确的，但目前真正安装的只是 1473 K 的单点 constant law：

```text
rho = 2750
beta = 1e-5
mu = 4
cp = 1150
k = 4
sigma = 16
```

因此：

> **sigma(T) infrastructure 已有，但 French sigma(T) physics 还没有完成。**

论文明确要求温度相关的 viscosity / cp / thermal conductivity / electrical conductivity，并且 EM–thermal coupling 是通过：

\[
\sigma(T)
\leftrightarrow
Q_{\rm Joule}
\]

进行。

## 快速实施顺序

先只上：

\[
\sigma(T)
\]

不要同时上全部 `mu/k/cp`。

链：

\[
T
\rightarrow
\sigma(T)
\rightarrow
A_{\rm glass}+\phi+J
\rightarrow
Q
\rightarrow
T
\]

并对 sigma 使用 under-relaxation。

等这个稳定，再加入：

\[
\mu(T),\quad k(T),\quad c_p(T)
\]

---

# 12. Stage 4：French natural full coupling

之后加入：

- Joule heating
- thermal conduction
- side wall effective cooling
- bottom cooling
- free-surface convection
- free-surface radiation
- Boussinesq buoyancy
- natural convection

形成：

\[
T
\rightarrow
\sigma
\rightarrow
EM
\rightarrow
Q
\rightarrow
T
\rightarrow
buoyancy
\rightarrow
u
\rightarrow
T
\]

Natural 论文复现使用 target-power：

```text
P_glass = 50 kW
```

每次 EM 更新：

1. 当前 T -> sigma(T)
2. self-induced EM Picard 收敛
3. 获得 raw Q distribution
4. 按 target 重新标定 coil current
5. 将最终 Q 放入 thermal solver

传感器此时仍只做 virtual probes。

---

# 13. Stage 5：Stirred case

Natural 完整跑通后再加入：

```text
H = 0.215 m
P = 60 kW
N = 10 rpm
```

复用同一 EM/thermal solver。

再加入：

- stirrer geometry
- rigid rotation
- stirrer thermal boundary
- 3D fluid motion
- Boussinesq buoyancy

不要此时再重写 EM。

---

# 14. 底板与线圈 surface current 仍后置

当前甲方允许忽略 segmented crucible EM。

所以顺序保持：

```text
Stage 1.5 closeout
-> A_glass
-> sigma(T)
-> full material laws
-> natural full coupling
-> stirred full coupling
-> base surface current
-> coil surface current
```

底板和线圈面积分都不阻塞当前交付。

---

# 15. 给 Cursor 的明确回答

### Q1：`P_graph/P_recon=0.355` 能否进入 A_glass？

**可以，但先完成 Stage 1.5 的 dp/reload consistency、coil geometry 和 weak-current diagnostic。**
不需要先把 0.355 压到 1 或 `<5%`。

### Q2：是否做 Phase B，强制 `C_ij=C_ji`？

**现在不做。**
先用：

\[
I_{ij}^{w}=V_iq_{ij}
\]

建立真正的 weak-current antisymmetry/KCL diagnostic。
将来若为 PCG/对称系统需要，可把整个 equation row 乘 \(V_i\) 得到 symmetric weak form，而不是直接篡改 pointwise `C_ij`。

### Q3：production power 用谁？

继续：

```text
P_recon / JouleHeatEdgeReconComplex
```

作为 calibration 和 thermal source。

```text
P_graph / P_undirected
```

作为 weak discrete-energy audit。

### Q4：下一主功能是什么？

**A_glass Picard。**

但先修：
1. reload/dp mismatch；
2. coil R/height；
3. weak current audit；
4. 32/64/128 coil convergence；
5. 最小 dp convergence。

这些应是一轮很短的 closeout，不要继续无限打磨 Stage 1。
