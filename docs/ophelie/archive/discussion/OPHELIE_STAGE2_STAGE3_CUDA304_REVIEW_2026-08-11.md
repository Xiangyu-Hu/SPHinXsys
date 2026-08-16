# OPHELIE Stage 2 / Stage 3 审查与 CUDA 304 故障定位计划
**日期：2026-08-11**

## 0. 决策摘要

本轮总体进展是正向的。

### Stage 1.5
此前要求的两个关键修正已经落地：

- reload 当前确实是 `dp=0.015 m`：
  - `n = 10474`
  - 每粒子 `VolumetricMeasure = 3.375e-6 = 0.015^3 m^3`
  - `sum(V) = 0.03534975 m^3`
  - 解析圆柱体积 `0.036324665 m^3`
  - 体积误差约 `-2.684%`
- 线圈几何已改为：
  - `R_coil = 0.285 m`
  - `N_turn = 7`
  - `z_min = -0.0225 m`
  - `z_max = 0.2075 m`
  - stack span = `0.230 m`
  - 与 185 mm 玻璃浴中心对齐
- weak-current audit 已改为使用：
  \[
  I_{ij}^{w}=V_iq_{ij}
  \]
  而不是简单比较 `q_ij + q_ji`。

因此此前 Stage 1.5 的核心配置问题已基本关闭。

### Stage 2
`A_glass` complex edge-flux Picard 已经进入正确工作区间：

```text
final_J_rel = 1.44154e-05
P_joule_W = 2.2428
A_ind_over_A_coil = 0.160174
B_ind_over_B_coil = 0.135285
phi_eq_res_vol = 1.78261e-4
```

`J_rel` 收敛很好；`phi_eq_res_vol` 没过 `1e-4` gate，但当前 inner GMRES 默认 tolerance 本身就是 `1e-4`，因此这个 fail 更像 gate/solver-tolerance 设计问题，不是 Picard 物理解发散。

Stage 2 可以作为 engineering pass 继续推进；正式 closeout 时只需增加一个更紧的 inner GMRES tolerance（建议 `1e-6`）重跑一次。

### Stage 3 CUDA 错误
当前：

```text
CUDA_ERROR_OPERATING_SYSTEM = 304
OS call failed
```

**优先判断为软件 / SYCL-UR / host-device 生命周期或上下文管理问题，而不是硬件故障。**

当前最强软件嫌疑是：

> Stage 3 在 EM GPU pipeline 已运行之后，才动态 `registerStateVariable` 新的 Temperature / DeltaT / thermal auxiliary fields，再立即创建/同步 device delegated data。

第二嫌疑是：

```cpp
DelegatedData(ExecutionPolicy{})
```

通过一个临时新构造的 execution policy 去确保 device allocation，可能引入不一致的 queue/context/device 生命周期。

不要先升级 CUDA/driver。先用最小二分测试和 API trace 把第一次失败的 native call 定位出来。

---

# 1. Stage 2 代码审查

## 1.1 Picard 主线符合预期

complex edge-flux self induction：

\[
A_{\rm total}=A_{\rm coil}+A_{\rm glass}
\]

并通过：

\[
J^n
\rightarrow
A_{\rm glass}^n
\rightarrow
\phi^{n+1}
\rightarrow
J^{n+1}
\]

迭代。

当前代码中：

```cpp
ComputeOphelieGlassSelfInducedBiotSavartCK
CombineOphelieCoilAndInducedVectorPotentialCK
execOphelieComplexEdgeFluxSolveReconAndPower
```

组合顺序合理。

欠松弛当前作用在 reconstructed complex J 上：

```cpp
hostRelaxVecdFieldTowardPrevious(...)
```

当前 `relax=0.3` 下 J 已稳定到 `1.44e-5`。

## 1.2 当前 Stage 2 fail 不应被解释为 self-induction 失败

运行设置：

```text
self-induction-tol = 1e-3
self-induction-phi-tol = 1e-4
```

但全局：

```cpp
phi_gmres_tolerance_ = 1.0e-4;
```

最终：

```text
phi_imag_rel_res = 1.72e-4
phi_real_rel_res = 1.39e-4
phi_eq_res_vol = 1.78e-4
```

因此 outer gate 与 inner solve accuracy 同量级，甚至要求 outer equation audit 比 inner solver 实际完成度更严格。

### 建议

新增 CLI：

```text
--phi-gmres-tol=
```

Stage 2 formal closeout：

```text
phi_gmres_tol = 1e-6
self_induction_phi_tol = 2e-4 或 1e-4
```

先看 `phi_eq_res_vol` 是否随 inner tolerance 继续下降。

不要为了这个 gate 停止 thermal unit-test 开发。

---

# 2. Stage 2 结果的物理含义

当前：

\[
\frac{\|A_{\rm ind}\|}{\|A_{\rm coil}\|}
\approx0.160
\]

\[
\frac{\|B_{\rm ind}\|}{\|B_{\rm coil}\|}
\approx0.135
\]

说明在 `sigma=16 S/m`、282 kHz 下，玻璃 self-induction 不是零阶小量。

当前 `P_joule=2.2428 W` 是 reference current 条件下的 raw load response，不应与后续 50 kW target power 混淆。

之后 target-power 模式只是把收敛后的 spatial Q distribution 统一缩放到 French 论文对应总玻璃焦耳功率。

仍需补做但不阻塞 Stage 3 runtime 修复的 Stage 2 regression：

```text
low-sigma limit
I -> 2I scaling
relaxation factor independence
```

---

# 3. Stage 3 当前实际上是什么

`test_3d_ophelie_french_complex_joule_to_heat_one_way` 当前是：

\[
\text{coil-only complex edge-flux EM}
\rightarrow
Q
\rightarrow
\text{one-way thermal}
\]

因为 `OphelieParameters::enable_self_induction_` 默认是 false，而这个 case 没有打开它，也没有调用 Stage 2 Picard pipeline。

因此它目前是一个 **thermal handoff regression test**，不是最终的：

\[
A_{\rm coil}+A_{\rm glass}
\rightarrow
Q
\rightarrow
T
\]

这并非错误，但文档必须说清楚。

此外本次失败命令没有：

```text
--use-literature-thermal
```

所以 thermal material 仍是 prototype：

```text
rho = 2500
cp = 1200
k = 1
T0 = 300 K
```

这对 unit test 可以接受，但不是 French literature run。

---

# 4. CUDA 304 为什么不像硬件错误

Stage 2 在同一 CUDA/SYCL 环境下可以运行重得多的 self-induction / EM pipeline。

Stage 3：

- `thermal_steps=3` 会失败；
- 缩成 `thermal_steps=1` 仍相同失败；
- 没有 NaN / divergence 日志；
- 错误码不是 `CUDA_ERROR_ILLEGAL_ADDRESS (700)`；
- 错误是 `CUDA_ERROR_OPERATING_SYSTEM (304)`。

因此第一优先应检查：

- CUDA/UR API 调用；
- queue/context/device 一致性；
- device buffer 创建；
- host-device sync；
- field registration / allocation 时机；
- runtime loader。

而不是先认定 GPU 显存、硬件损坏或数值发散。

---

# 5. 当前最强嫌疑：late registration of thermal fields

Stage 3 main 先执行：

```cpp
RegisterOphelieGlassFields register_glass(...);
AssignOphelieGlassSigmaCK assign_sigma(...);
assign_sigma.exec();              // 已经进入 device path
```

随后：

```cpp
runFrenchReducedEmPipeline(...)   // 大量 device work
```

然后在 `runFrenchReducedEmThenJouleHeatOneWay()` 中才：

```cpp
registerOphelieJouleHeatTemperatureField(...)
```

该函数内部：

```cpp
particles.registerStateVariable<Real>("OphelieThermalDeltaT", 0);
particles.registerStateVariable<Real>("Temperature", t_initial);
...
syncVariableToDevice(...)
```

diffusion 路径还会在 EM 之后注册：

```text
ThermalLaplaceT
ThermalConductivity
ThermalBoundaryMask
```

这和现有 EM 的安全模式不一致：EM field 是在第一批 device dynamics 之前集中注册的。

### 决策

把 thermal fields 的 registration 全部移到 **第一次 device execution 之前**。

建议新建：

```cpp
RegisterOphelieThermalFields
```

在：

```cpp
RegisterOphelieGlassFields
```

之后、`assign_sigma.exec()` 之前执行。

即使 thermal diffusion 关闭，也可以先注册最基本：

```text
Temperature
OphelieThermalDeltaT
```

如果 diffusion 启用，则其 aux fields 也必须在首个 offload 前注册。

---

# 6. 第二嫌疑：临时 ExecutionPolicy 的 DelegatedData

当前辅助函数：

```cpp
template <class ExecutionPolicy, typename DataType>
inline void ensureOphelieVariableDelegatedOnDevice(...)
{
    ...
    variable->DelegatedData(ExecutionPolicy{});
}
```

这里人为创建一个新的：

```cpp
ExecutionPolicy{}
```

来触发 delegated allocation。

这值得怀疑。

正式 dynamics 本身会通过自己的 execution policy / delegated data path 管理 device data；额外用临时 policy 构造 device allocation，可能导致 queue/context/device 生命周期与真实 kernel path 不一致。

### 决策

优先删除这个手工 `ensure...`。

让：

```cpp
StateDynamics<MainExecutionPolicy,...>
```

沿框架正常取得 delegated data。

如果确实必须预创建 device storage，则必须复用与实际 dynamics 相同的 execution context，而不是新构造临时 policy。

---

# 7. 第三个嫌疑：UR loader/runtime

Stage2/Stage3 CMake 都手工搜索并 link：

```text
/opt/intel/oneapi/compiler/*/lib/libur_loader.so*
```

并优先找：

```text
libur_loader.so.0.11
```

由于两个 case CMake 相同而 Stage2 可跑，它不太像第一根因。

但 Stage3 可能第一次触发 Stage2 没用到的 memory/sync API，所以仍需检查 Stage2 与 Stage3：

```bash
ldd <binary> | grep -E 'sycl|ur|cuda'
```

必须完全一致。

长期建议不要由每个 case 手工挑 UR loader；最好统一交给整个 SPHinXsys/SYCL toolchain 配置管理。

---

# 8. 下一轮必须按二分法定位，不要盲改环境

## A. 加 stage markers

所有 marker 用：

```cpp
std::cerr << "... " << std::endl;
```

确保 flush。

至少：

```text
[S3] after system construction
[S3] after reload
[S3] after system configuration
[S3] after EM field registration
[S3] after thermal field registration
[S3] before assign_sigma
[S3] after assign_sigma
[S3] before EM
[S3] after EM
[S3] before Q handoff
[S3] after Q handoff
[S3] before thermal sync
[S3] after thermal sync
[S3] before thermal kernel
[S3] after thermal kernel
[S3] before reductions
[S3] after reductions
```

在关键 GPU phase 后加 queue/runtime 可用的 `wait_and_throw()`，使 asynchronous error 在真正发生的 phase 抛出。

## B. 建四个最小路径

### T0 — thermal-only

同一 reload：

- 不跑 EM；
- early-register thermal fields；
- synthetic constant Q；
- 跑一个 `T += Q dt/(rho cp)`。

若 T0 失败：thermal field/device path 本身有问题。

若 T0 通过：thermal kernel 基本无问题。

### E0 — EM-only Stage3 skeleton

保持 Stage3 初始化和 EM，但完全不 register thermal fields、不做 thermal。

若 E0 通过：EM 本身没问题。

### E1 — early-register thermal fields + EM only

不跑 thermal kernel。

若 E1 通过：early registration 不破坏 EM。

### E2 — early-register + EM + one thermal step

这是最小完整 handoff。

四个测试可以非常快地把问题从“环境大海捞针”缩到某一小段。

---

# 9. 运行时追踪命令

先看实际设备 selector：

```bash
sycl-ls --ignore-device-selectors
sycl-ls --verbose
```

然后按 `sycl-ls` 输出的 CUDA selector 固定到单卡，例如实际输出若允许：

```bash
ONEAPI_DEVICE_SELECTOR=cuda:0 ...
```

或对应它显示的精确 selector。

建议：

```bash
unset ONEAPI_DEVICE_SELECTOR
SYCL_PI_TRACE=2 ./test_3d_ophelie_french_complex_joule_to_heat_one_way ... \
  2>&1 | tee stage3_sycl_trace.log
```

若需要更完整：

```bash
SYCL_PI_TRACE=-1 ...
```

重点看最后一个成功的 UR call 和紧接着返回 304 的 API。

同一 selector 下：

```text
Stage2 -> Stage3 -> Stage2
```

连续跑。

如果前后 Stage2 都正常而 Stage3 稳定报 304，更强地指向 Stage3 代码路径，而非 GPU 硬件随机故障。

---

# 10. 基础环境检查

只做诊断，不先升级：

```bash
nvidia-smi
nvidia-smi -q -d ECC,ERROR,COMPUTE
ulimit -l
ls -l /dev/nvidia*
```

比较：

```bash
ldd test_3d_ophelie_french_self_induction_picard | grep -E 'sycl|ur|cuda'
ldd test_3d_ophelie_french_complex_joule_to_heat_one_way | grep -E 'sycl|ur|cuda'
```

如可用：

```bash
compute-sanitizer --tool memcheck ./test_3d_ophelie_french_complex_joule_to_heat_one_way ...
```

如果在任何 CUDA kernel memory violation 之前就直接返回 304，则继续聚焦 runtime / buffer / context，而不是数值 kernel illegal access。

---

# 11. exception 后的 segfault 不要当成第一根因

日志顺序：

```text
sycl::exception
CUDA 304
terminate
Segmentation fault
```

当前应首先定位抛出 `sycl::exception` 的 native API。

后面的 segfault 很可能是异常退出/析构阶段的二次错误。

建议在 Stage3 顶层临时加入：

```cpp
try
{
    ...
}
catch (const sycl::exception &e)
{
    std::cerr << "[S3][SYCL] " << e.what() << std::endl;
    return 2;
}
catch (const std::exception &e)
{
    std::cerr << "[S3][STD] " << e.what() << std::endl;
    return 3;
}
```

以避免原始错误被 secondary crash 淹没。

---

# 12. 修复后 Stage 3 的正确推进顺序

先证明：

```text
coil-only complex edge-flux
-> Q
-> one thermal step
```

energy closure 正确。

再把 Stage 2 self induction 真正接进 Stage3：

\[
A_{\rm coil}+A_{\rm glass}
\rightarrow
Q
\rightarrow
T
\]

然后才开始：

\[
\sigma(T)
\]

的真正反馈：

\[
T
\rightarrow
\sigma(T)
\rightarrow
EM
\rightarrow
Q
\rightarrow
T
\]

不要在 CUDA 304 未定位时同时引入 `sigma(T)`，否则故障维度会再次扩大。

---

# 13. 代码/文档清理

### French geometry defaults

`electromagnetic_ophelie_french_reduced_geometry.h` 仍保留旧 legacy 默认/注释：

```text
650 mm
300 kHz
Rglass=0.325
coil_R=0.40
8 loops
```

尽管 Stage2/3 CLI 已正确 override，仍建议建立明确 preset：

```text
LegacyReducedDemo
FrenchNaturalEREBUS
FrenchStirredEREBUS
```

避免后续漏参数时误跑 legacy 几何。

### Stage3 README

当前 README 仍有旧粒子数量/几何描述，需要同步到：

```text
R=0.25
H=0.185
dp=0.015
n=10474
```

并明确：

> 当前 Stage3 是 coil-only EM -> thermal handoff regression，不是 self-induced full-coupling case。

---

# 14. 给 Cursor 的明确决定

1. **不要先更新 CUDA/driver。**
2. Stage2 允许 engineering 放行；正式 closeout 时增加 `phi-gmres-tol=1e-6` 重跑一次。
3. Stage3 CUDA 304 先按软件路径查。
4. 第一修改：所有 thermal state variables 在第一次 device kernel 前注册。
5. 第二修改：去掉/重构 `DelegatedData(ExecutionPolicy{})` 临时 policy 路径。
6. 建 T0/E0/E1/E2 四个最小回归测试。
7. 用 `SYCL_PI_TRACE=2` 和单 GPU selector 定位第一条失败 API。
8. Stage3 one-way 通过后，再接入 Stage2 `A_glass`。
9. `A_glass -> Q -> T` 通过后，再进入 `sigma(T)` 双向耦合。
10. 暂时不要把 runtime 错误和 French thermal physics 同时改，保持变量隔离。

## 最重要判断

当前没有证据说明 RTX GPU 硬件本身有问题。

更符合现有证据的解释是：

\[
\boxed{
\text{Stage3 新增 thermal fields 的 device registration/sync/context path}
}
\]

触发了 Stage2 从未走过的 CUDA/UR 运行时操作。

先把第一次返回 304 的 native API 准确抓出来，再决定是否需要环境层修改。
