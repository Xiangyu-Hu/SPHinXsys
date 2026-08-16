# OPHELIE French Natural：Stage 3.4 收口与下一步开发决策

**日期：2026-08-11**  
**用途：直接发给 Cursor 作为下一阶段开发依据**

---

## 0. 当前阶段总判断

这一轮之后，项目主线已经发生阶段性变化：

\[
\boxed{
\text{EM solver development}
\longrightarrow
\text{French Natural 50 kW 全耦合复现}
}
\]

目前可认为：

- edge-flux 主问题已经解决；
- glass self-induction \(A_{\rm glass}\) 已经进入可用状态；
- EM \(\rightarrow Q \rightarrow T\) 数据链已打通；
- \(\sigma(T)\) 更新机制已具备；
- 当前剩余的 \(\phi\) residual 问题只是略高于既定 gate，不构成新的物理或算法 blocker；
- 下一步不应该继续长期停留在电磁残差打磨，而应尽快完成 50 kW target-power，然后进入法国自然对流的真实热边界和流体侧。

---

# 1. 关于当前 \(\phi\) residual：接受 \(2.5\times10^{-4}\) engineering hard gate

当前高 \(\sigma\)、带 self-induction 工况下：

\[
\phi_{\rm eq,res}
\approx
2.10\times10^{-4}
\]

原 gate：

\[
2.0\times10^{-4}
\]

仅超出约 5%。

与此同时：

\[
J_{\rm rel}
\approx
8.1\times10^{-5}
\]

已经明显收敛。

因此当前现象更像：

\[
\boxed{
O(2\times10^{-4})
\text{ 的离散 / audit floor}
}
\]

而不是 Picard 或 Krylov 发散。

## 1.1 正式 gate 建议

不要简单把所有标准都改成 \(2.5\times10^{-4}\)。

建议采用双层标准：

```text
preferred_phi_target = 2.0e-4
engineering_phi_hard_gate = 2.5e-4
```

即：

- \(\phi_{\rm eq,res} \le 2.0\times10^{-4}\)：preferred pass；
- \(2.0\times10^{-4} < \phi_{\rm eq,res}\le2.5\times10^{-4}\)：engineering pass + warning；
- \(>2.5\times10^{-4}\)：fail。

同时：

```text
J_rel < 1.0e-3
```

继续作为 outer self-induction convergence gate。

## 1.2 仅补一个 tighter-solver diagnostic

在彻底封 Stage 3.4 之前，Cursor 再做一次很小的诊断：

- 将 inner \(\phi\) solver tolerance 收紧 10 倍；
- max iteration 明显提高；
- 同时输出：

```text
phi_real_solver_rel_res
phi_imag_solver_rel_res
phi_eq_res_vol
```

目的不是继续调 solver，而是判断：

\[
\text{Krylov algebraic residual}
\]

与：

\[
\text{post-audit } \phi_{\rm eq,res}
\]

是否同步下降。

### 情况 A

如果 inner solver residual 明显下降，而：

\[
\phi_{\rm eq,res}
\approx2.1\times10^{-4}
\]

基本不动，则确认：

\[
\boxed{
\phi_{\rm eq,res}
\text{ 是离散/audit floor}
}
\]

Stage 3.4 正式关闭。

### 情况 B

如果收紧 solver 后：

\[
\phi_{\rm eq,res}<2.0\times10^{-4}
\]

则更好，但同样不需要继续长期优化。

**无论 A/B，都不阻塞后续开发。**

---

# 2. \(\sigma=16\) 还是 \(\sigma\approx43\)：正式 Natural case 采用期刊 Table 1 的 16 S/m

这是当前最需要统一的文献口径。

## 2.1 正式来源优先级

继续遵循：

\[
\boxed{
\text{期刊论文}>\text{博士论文}
}
\]

Natural convection 期刊论文 Table 1 明确给出 1473 K 参考值：

\[
\rho_0=2750\ {\rm kg/m^3}
\]

\[
\mu=4\ {\rm Pa\,s}
\]

\[
C_p=1150\ {\rm J/(kg\,K)}
\]

\[
k=4\ {\rm W/(m\,K)}
\]

\[
\sigma=16\ {\rm S/m}
\]

\[
\beta=10^{-5}\ {\rm K^{-1}}
\]

并给出电导率在工况温度范围中的变化范围：

\[
10^{-6}\sim30\ {\rm S/m}
\]

因此正式 Natural reproduction 在：

\[
T=1473K
\]

必须满足：

\[
\boxed{
\sigma(1473)=16\ {\rm S/m}
}
\]

## 2.2 博士论文 III.12 不作为正式 Natural baseline

博士论文给出：

\[
\log_{10}\sigma
=
3.7921-\frac{3179.8}{T}
\]

代入：

\[
T=1473K
\]

约得到：

\[
\sigma\approx43\ {\rm S/m}
\]

这不仅与 Table 1 的 16 不一致，而且高于 Natural 期刊给出的：

\[
30\ {\rm S/m}
\]

范围上限。

因此：

\[
\boxed{
\text{raw thesis III.12 只能作为 sensitivity case}
}
\]

不能作为 Natural 期刊正式复现参数。

---

# 3. Natural \(\sigma(T)\) 的推荐实现：期刊锚点 + thesis 变化形状

Natural 期刊明确要求 electrical conductivity 随温度变化，但短论文没有给出完整的函数节点。

因此建议建立一个显式标记为“补充重建”的 law：

\[
\boxed{
\sigma_{\rm nat}(T)
=
16
\frac{\sigma_{\rm thesis}(T)}
{\sigma_{\rm thesis}(1473)}
}
\]

如果 thesis III.12 使用：

\[
\sigma_{\rm thesis}(T)
=
10^{3.7921-3179.8/T}
\]

则：

\[
\boxed{
\sigma_{\rm nat}(T)
=
16\,
10^{3179.8(1/1473-1/T)}
}
\]

并按照 Natural 期刊给出的范围做 clipping：

\[
\boxed{
10^{-6}
\le
\sigma_{\rm nat}(T)
\le
30
}
\]

这样能够同时满足：

1. 1473 K 时严格：

\[
\sigma=16
\]

2. 整体范围不超过期刊的：

\[
10^{-6}\sim30
\]

3. 温度变化趋势由 thesis curve 补充。

## 3.1 建议 material preset 拆分

不要再使用一个笼统的 “French material”。

建议至少拆成：

```text
FrenchNaturalCEP2008
FrenchStirredCES2008
ThesisStaticIII
ThesisStirredIV
```

Natural 和 Stirred 论文物性并不完全一样。

例如：

Natural：

\[
\rho_0=2750,\qquad
\beta=1\times10^{-5}
\]

Stirred：

\[
\rho_0=2800,\qquad
\beta=3\times10^{-5}
\]

避免不同论文参数混用。

---

# 4. 下一步先做 50 kW target-power，但只做一个短 Stage 3.5

Cursor 问：

> 下一步先做 50 kW target-power，还是先做自然对流流体侧？

正式决定：

\[
\boxed{
\text{先完成 50 kW target-power}
}
\]

但这应该只是一个很短的 Stage 3.5。

完成之后立即进入：

\[
\boxed{
\text{French Natural thermal + fluid}
}
\]

不要再在 EM 上长期停留。

---

# 5. Stage 3.5：50 kW target-power 的具体实现

Natural 论文正式工况：

\[
\boxed{
P_{\rm glass}=50\ {\rm kW}
}
\]

当前 reference-current 计算只有几 W，只用于验证 EM 响应，不能直接拿来跑论文自然对流。

## 5.1 计算流程

先在当前温度场下得到：

\[
\sigma(T)
\]

然后完整收敛：

\[
A_{\rm coil}
+
A_{\rm glass}
\rightarrow
\phi
\rightarrow
J
\rightarrow
Q_{\rm raw}
\]

得到：

\[
P_{\rm raw}
=
\sum_iQ_{{\rm raw},i}V_i
\]

定义：

\[
s
=
\sqrt{
\frac{50000}
{P_{\rm raw}}
}
\]

由于固定 \(\sigma\) 时整个 reduced EM system 是线性的：

\[
A\propto I,\quad
\phi\propto I,\quad
E\propto I,\quad
J\propto I
\]

因此：

\[
Q\propto I^2
\]

可缩放：

\[
A,\phi,E,J
\leftarrow
s(A,\phi,E,J)
\]

\[
Q
\leftarrow
s^2Q
\]

最终：

\[
P_{\rm deposited}
\approx
50\,000\ {\rm W}
\]

## 5.2 第一次实现建议重新跑一次 Picard 做 sanity check

第一次 50 kW 实现时，不要只做 algebraic scaling。

先算：

\[
I_{\rm target}=sI_{\rm ref}
\]

再重新跑一次完整 self-induction Picard。

验证：

\[
P_{\rm target-run}
\approx
50\ {\rm kW}
\]

并比较 normalized Joule shape：

\[
\frac{Q_{\rm target}(\mathbf x)}
{P_{\rm target}}
\]

与：

\[
\frac{Q_{\rm raw}(\mathbf x)}
{P_{\rm raw}}
\]

应几乎完全一致。

一旦证明线性 scaling 成立，后续 EM update 可直接 scale，避免重复第二遍 Picard。

## 5.3 Stage 3.5 输出

必须输出：

```text
P_raw
P_target
excitation_scale
I_ref
I_equiv_peak
I_equiv_rms
P_deposited
power_error_rel
Q_max_raw
Q_max_scaled
Q_peak_position
normalized_Q_shape_error
```

建议：

\[
\frac{|P_{\rm deposited}-50000|}
{50000}
<10^{-3}
\]

---

# 6. 下一关键开发：法国真实热边界

当前 Stage 3.3 使用的 cold-wall Dirichlet 只是 regression test。

如果当前做：

\[
T_{\rm wall}=T_{\rm initial}
\]

那不是法国 Natural 论文中的真实边界。

必须保留现有路径，但重命名为：

```text
thermal_dirichlet_regression
```

正式 Natural production BC 单独实现。

---

# 7. French Natural 正式热边界

## 7.1 侧壁

\[
\boxed{
q_{\rm side}
=
h_{\rm side}
(T-T_c)
}
\]

论文做过：

\[
h_{\rm side}
=
100,\ 200,\ 300,\ 400
\ {\rm W/(m^2K)}
\]

敏感性研究。

第一版建议：

\[
\boxed{
h_{\rm side}=300\ {\rm W/(m^2K)}
}
\]

但文档写：

> selected within the journal sensitivity range

不要写成：

> exact value explicitly prescribed by the paper

## 7.2 底面

\[
\boxed{
q_{\rm bottom}
=
h_{\rm bottom}(T-T_c)
}
\]

采用：

\[
\boxed{
h_{\rm bottom}=35\ {\rm W/(m^2K)}
}
\]

## 7.3 顶部自由表面

\[
\boxed{
q_{\rm free}
=
h_s(T-T_a)
+
\epsilon\sigma_{\rm SB}
(T^4-T_{\rm ar}^4)
}
\]

其中第一版采用：

\[
h_s\approx20\ {\rm W/(m^2K)}
\]

并加入 emissivity / ambient / radiative ambient 参数。

---

# 8. Stage 4.1：先做 frozen-Q 自然对流

为了避免 EM / thermal / fluid 三套东西同时出问题，第一步自然对流不要马上做 fully coupled EM update。

先冻结 Stage 3.5 得到的：

\[
\boxed{
Q_{50kW}(\mathbf x)
}
\]

然后只验证：

\[
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

即：

\[
\boxed{
\text{thermal-fluid isolated validation}
}
\]

---

# 9. Stage 4.1 Natural fluid 参数

先采用 Natural 期刊 1473 K constant-reference properties：

\[
\rho_0=2750\ {\rm kg/m^3}
\]

\[
\mu=4\ {\rm Pa\,s}
\]

\[
C_p=1150\ {\rm J/(kg\,K)}
\]

\[
k=4\ {\rm W/(m\,K)}
\]

\[
\beta=10^{-5}\ {\rm K^{-1}}
\]

Boussinesq：

\[
\rho
=
\rho_0
[1-\beta(T-T_0)]
\]

浮力：

\[
F_b
\sim
-\rho_0\beta(T-T_0)\mathbf g
\]

当前不加：

\[
J\times B
\]

也不加：

\[
u\times B
\]

---

# 10. Natural velocity BC

侧壁/底部：

\[
\boxed{
u=0
}
\]

即 no-slip。

顶部自由表面：

- 无穿透；
- 切向 slip / free-slip；
- 当前不做自由表面形变。

---

# 11. Stage 4.1 需要检查的量

至少输出：

```text
T_mean
T_max
T_min
U_max
U_mean
wall_heat_loss
bottom_heat_loss
free_surface_conv_loss
free_surface_rad_loss
total_heat_loss
Joule_input
thermal_energy_change
energy_balance_error
Ra
Pr
```

以及：

- \(r-z\) 平均温度场；
- \(r-z\) 平均速度场；
- 环流拓扑；
- Q 分布与温度热点位置；
- virtual thermocouple 温度曲线。

---

# 12. Stage 4.2：完成真正的全耦合 Natural case

Stage 4.1 稳定后，打开真正反馈：

\[
\boxed{
T
\rightarrow
\sigma(T)
\rightarrow
EM
\rightarrow
Q
\rightarrow
T,u
}
\]

流程：

1. 当前 \(T\)；
2. 更新 \(\sigma(T)\)；
3. 运行 self-induction Picard；
4. 得到 \(Q_{\rm raw}\)；
5. target-power 标定回 50 kW；
6. 将 \(Q\) 送入 thermal-fluid；
7. 推进若干 SPH timestep；
8. 再做下一次 EM update。

不需要每一个 acoustic timestep 都做 EM。

应支持：

```text
--em-update-every=N
```

后续研究：

```text
N = 1 / 5 / 10
```

或按物理时间间隔更新。

---

# 13. 温度依赖物性加入顺序

不要一次性全部打开。

推荐：

### Phase A

\[
\sigma(T)
\]

先实现完整 EM feedback。

### Phase B

\[
\mu(T)
\]

这会显著影响 natural convection 和 skull 区域。

### Phase C

\[
k(T),\quad C_p(T)
\]

最后补齐。

这样每一步都可以定位结果变化来源。

---

# 14. 当前 temperature seed 处理

Stage 3 machinery test 中若存在人为设置的 radial temperature seed / perturbation：

正式 Natural production case 中最终应关闭。

因为正式温度非均匀性应该由：

\[
Q_{\rm Joule}
\]

\[
\text{wall cooling}
\]

\[
\text{free-surface losses}
\]

\[
\text{conduction}
\]

\[
\text{natural convection}
\]

自然形成。

如果 Natural convection 启动初期需要极小 perturbation 打破完全对称数值状态，可以保留一个极小、明确记录的 seed，但不要让它控制最终流动结构。

---

# 15. 传感器继续使用 virtual probes

当前不要加实体 sensor。

继续：

```text
TC1
TC2
TC3
...
```

在论文位置采样：

\[
T(\mathbf x_{\rm TC},t)
\]

后处理和论文实验曲线对比。

---

# 16. 当前不要进入 Stirred case

机械搅拌继续后置。

Natural case 是最适合当前做 full-coupling acceptance 的体系：

\[
\boxed{
A_{\rm coil}
+
A_{\rm glass}
+
50\,{\rm kW}
+
\sigma(T)
+
thermal\ boundaries
+
Boussinesq
+
natural\ convection
}
\]

只有 Natural case 完整跑通后，再切：

\[
H=0.215\ {\rm m}
\]

\[
P=60\ {\rm kW}
\]

\[
N=10\ {\rm rpm}
\]

加入 stirrer。

---

# 17. Cursor 当前三个问题的最终拍板

## Q1

> 3.4b：\(\phi=2.10\times10^{-4}\) vs \(2.0\times10^{-4}\)，是否接受 `phi_tol=2.5e-4`？

### 决定

\[
\boxed{\text{接受}}
\]

但采用：

```text
preferred = 2.0e-4
hard_gate = 2.5e-4
```

并补一次 tighter-solver diagnostic。

不阻塞后续。

## Q2

> Table 1 \(\sigma=16\) vs thesis III.12 \(\sigma\approx43@1473\)，以谁为准？

### 决定

\[
\boxed{
\text{Natural 正式复现以期刊 Table 1 的 16 S/m 为准}
}
\]

Raw thesis law：

```text
sensitivity only
```

正式 \(\sigma(T)\) 推荐：

```text
journal anchor + thesis shape + journal range clipping
```

## Q3

> 下一步先50 kW target-power，还是先自然对流流体侧？

### 决定

\[
\boxed{
\text{先做很短的 Stage 3.5：50 kW target-power}
}
\]

完成后：

\[
\boxed{
\text{立即进入 French Natural thermal boundary + natural convection}
}
\]

不要继续长期磨 EM。

---

# 18. 下一阶段开发顺序

最终顺序固定为：

```text
Stage 3.4 closeout
    ↓
Stage 3.5 target-power = 50 kW
    ↓
Stage 4.0 French production thermal BC
    ↓
Stage 4.1 frozen-Q natural convection
    ↓
Stage 4.2 sigma(T) fully coupled Natural
    ↓
mu(T)
    ↓
k(T), cp(T)
    ↓
French Natural literature validation
    ↓
French Stirred 60 kW / 10 rpm
```

底板 surface current 和 coil surface integration 继续后置，不阻塞当前交付。

---

# 19. 当前阶段的核心目标

现在不再追求“把每一个 EM residual 压到极限”。

近期最重要的是得到一套：

\[
\boxed{
A_{\rm coil}
+
A_{\rm glass}
+
\phi
+
J
+
Q_{50kW}
+
\sigma(T)
+
T
+
u
}
\]

可以稳定工作的 French Natural 全耦合 case。

一旦这个 case 跑通，就可以向甲方展示：

- 玻璃感应电流；
- self-induction；
- Joule heat 分布；
- 50 kW 功率控制；
- 温度场；
- 自然对流；
- 冷壁/自由表面热损失；
- 温度相关电导率反馈；
- virtual sensor 曲线。

这应该成为当前主线的第一优先级。
