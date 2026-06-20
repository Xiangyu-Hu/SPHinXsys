# OPHELIE-like French Reduced: Phi Discretization Consistency, Boundary Treatment, and Next-Stage Development Recommendations

> This document consolidates analysis of `phi_eq_res≈0.49` saturation, Neumann boundary corrections, corrected gradient, SPHinXsys/SYCL code structure, and next-stage Cursor development tasks. Goal: clear roadmap for follow-on coding, avoiding scattered CLI toggle trial-and-error.

---

## 0. Overall assessment of current issues

Current French reduced OPHELIE-like pipeline:

```text
multiloop Biot coil
  -> A_src / B_src
  -> φ solve: div(σ grad φ) = -ω div(σ A)
  -> E / J / JouleHeat
  -> literature-mode: P_raw ≈ 50 kW
```

Continuous equation:

\[
\nabla\cdot(\sigma\nabla\phi) = -\omega\nabla\cdot(\sigma \mathbf A)
\]

Postprocess:

\[
\mathbf E_{imag} = -\nabla\phi - \omega\mathbf A
\]

\[
\mathbf J_{imag} = \sigma\mathbf E_{imag}
\]

\[
Q_{Joule}=\frac{|\mathbf J|^2}{2\sigma}
\]

Observed behavior:

```text
phi_eq_res_vol ≈ 0.48–0.49 saturated
RHS zero-mean projection almost ineffective
RHS one-sided Neumann flux almost ineffective
Grad postprocess Neumann can drop Jn to ~1e-8 but destroys divJ / literature acceptance
Standalone corrected gradient worsens phi_eq_res from ~0.486 to ~0.685
```

Core judgment:

**This is not simply insufficient GMRES iterations or a single missing boundary flux term. The main issue is SPH discrete operator inconsistency: LHS, RHS, grad postprocess, and divJ diagnostic do not form a compatible `D/G` discrete pair.**

That is, the current matrix-free equation is:

\[
L_h\phi=b_h
\]

Where:

\[
L_h = D_h(\sigma G_h)
\]

\[
b_h = -\omega D_h(\sigma A)
\]

If `D_h`, `G_h`, `DivSigmaA`, `DivSigmaGrad`, `E/J` postprocess, and `DivJ` diagnostic are not the same compatible discrete system:

- GMRES internal residual can decrease;
- physical residual or `phi_eq_res_vol` may stall at a high floor;
- fixing boundary or gradient alone may improve a local metric but break overall `divJ`.

---

## 1. Fix engineering issue first: CLI prefix offset bug

Before mathematical/physical judgment, fix the CLI prefix length bug in `electromagnetic_ophelie_cli.h`.

Current code has patterns like:

```cpp
std::strncmp(arg, "--phi-gmres-max-outer-iter=", 25)
params.phi_gmres_max_outer_iterations_ = atoi(arg + 25);
```

But the string:

```text
--phi-gmres-max-outer-iter=
```

actual length is **27**, not 25.

High-risk items already found:

```text
--phi-gmres-max-outer-iter=   should be 27, may use 25
--phi-gmres-restart=          should be 20, may use 18
--phi-gmres-eq-res-tol=       should be 23, may use 22
--phi-eq-res-gate=            should be 18, may use 17
```

This causes:

```text
atoi / atof reads from wrong position, e.g. from "er=...", "t=...", or "=..."
final parameters may become 0 or wrong values
GMRES restart, outer iter, eq_res_tol, gate judgment polluted
```

### 1.1 Recommended fix

Do not hand-write magic numbers. Add helper:

```cpp
inline bool startsWith(const char* arg, const char* key)
{
    return std::strncmp(arg, key, std::strlen(key)) == 0;
}

inline const char* valueAfter(const char* arg, const char* key)
{
    return arg + std::strlen(key);
}
```

Then all parsing:

```cpp
const char* key = "--phi-gmres-restart=";
if (startsWith(arg, key))
{
    params.phi_gmres_restart_dimension_ =
        static_cast<UnsignedInt>(std::atoi(valueAfter(arg, key)));
}
```

### 1.2 Must print final params at startup

Every French reduced test run should print effective parameters:

```text
[ophelie] final_params:
  literature_mode = ...
  target_power = ...
  phi_gauge_penalty = ...
  phi_solver = ...
  phi_gmres_restart_dimension = ...
  phi_gmres_max_outer_iterations = ...
  phi_gmres_eq_res_tolerance = ...
  phi_eq_res_gate = ...
  phi_lhs_operator = ...
  phi_rhs_operator = ...
  phi_boundary_mode = ...
  phi_rhs_project_zero_mean = ...
  phi_boundary_grad_neumann = ...
  phi_boundary_lhs_grad_neumann = ...
  phi_gradient_correction = ...
  state_recording = ...
```

Otherwise actual case vs command line cannot be confirmed.

---

## 2. Second audit item: does RHS finalize actually enter GMRES solve?

Current pipeline flow:

```cpp
setupOpheliePhiImagRhsFromASrc(...);
finalizeOpheliePhiImagRhsHost(...);
solvePhiImag(...);
```

Must confirm:

**Does `solvePhiImag()` internally call `setupOpheliePhiImagRhsProblem()` or `setupOpheliePhiImagRhsFromASrc()` again?**

If solver regenerates RHS without calling again:

```cpp
finalizeOpheliePhiImagRhsHost(...)
```

then:

```text
--phi-rhs-project-zero-mean
--phi-boundary-mode=one-sided-neumann
```

may enter pre-solve diagnostic only, not the actual GMRES solve.

### 2.1 Recommended interface refactor

Clarify RHS setup/finalize and solve lifecycle:

```cpp
setupPhiRhs(...);
finalizePhiRhs(...);
solvePhiImagWithCurrentRhs(...);
```

Or unify finalize inside solver:

```cpp
solvePhiImag(...)
{
    setupOpheliePhiImagRhsFromASrc(...);
    finalizeOpheliePhiImagRhsHost(...);
    gmresSolveWithFinalizedRhs(...);
}
```

Do not allow external pipeline and internal solver each to generate RHS once.

### 2.2 Required diagnostics

Before GMRES in `solvePhiImag()`, output or record:

```text
rhs_l2_before_gmres
rhs_sum_before_gmres
rhs_min_before_gmres
rhs_max_before_gmres
rhs_hash_or_checksum_before_gmres
```

Record after finalize and before GMRES to ensure:

```text
RHS projection / Neumann correction actually changed the rhs used for GMRES
```

---

## 3. Q1 judgment: phi_eq_res≈0.49 is mainly not a GMRES issue

Current evidence:

1. RHS zero-mean projection almost unchanged `phi_eq_res`;
2. RHS Neumann flux almost unchanged `phi_eq_res`;
3. French real Biot boundary `|n·σA|≈1e-10` — boundary source term itself tiny;
4. Grad postprocess fixes `Jn` but destroys `divJ`;
5. Standalone corrected gradient worsens residual.

Priority judgment:

```text
phi_eq_res≈0.49 more likely:
  - LHS/RHS discrete range mismatch
  - grad/div incompatible
  - operator inconsistency from boundary kernel truncation
not GMRES/preconditioner insufficient.
```

Before concluding, complete two audits:

```text
1. CLI parameters actually take effect
2. finalized RHS actually enters GMRES solve
```

---

## 4. Discrimination tests to add

### 4.1 Operator determinism test

For same field `phi`, apply LHS twice:

```text
L1 = L(phi)
L2 = L(phi)
```

Check:

\[
\frac{\|L_1-L_2\|}{\|L_1\|} < 10^{-10}
\]

If not satisfied, possible:

```text
device/host sync issue
temporary variable pollution
relation update inconsistency
matvec nondeterminism
```

### 4.2 Operator linearity test

Random fields `x, y`, check:

\[
L(x+y)-L(x)-L(y)
\]

and:

\[
L(\alpha x)-\alpha L(x)
\]

Acceptance recommendation:

```text
linearity_add_rel   < 1e-8
linearity_scale_rel < 1e-8
```

If fail, GMRES mathematical premise invalid.

### 4.3 Discrete self-consistency MMS

Key test for next stage.

Do not start with continuous analytic `A`. Build strictly compatible RHS from discrete operators.

Flow:

1. Given manufactured `phi_exact`.
2. With current `ComputeOphelieScalarPhiGradientCK`:

   \[
   G_h\phi_{exact}
   \]

3. Construct:

   \[
   A_h = -\frac{1}{\omega}G_h\phi_{exact}
   \]

4. With current `DivSigmaA` generate RHS:

   \[
   b_h=-\omega D_h(\sigma A_h)
   \]

5. With current LHS:

   \[
   L_h\phi_{exact}=D_h(\sigma G_h\phi_{exact})
   \]

6. Check:

   \[
   \|L_h\phi_{exact}-b_h\|
   \]

If this test cannot approach 0, current code still has:

```text
LHS/RHS sign, scaling, fields, D/G chain inconsistent
```

If near 0 but real Biot A residual still ~0.49, then real Biot A RHS has range projection error vs current discrete LHS.

### 4.4 Small dense matrix range test

Small particle count, e.g.:

```text
N < 1000
```

Explicitly assemble:

\[
L_h
\]

Dense least-squares or SVD:

\[
\min_\phi \|L_h\phi-b_h\|
\]

Judge:

```text
If dense least-squares minimum residual also ~0.49: range issue.
If dense solve drops low but matrix-free GMRES does not: solver/preconditioner/matvec implementation issue.
```

---

## 5. Q2 judgment: why grad postprocess destroys divJ

Current grad postprocess:

\[
\nabla\phi_i \leftarrow \nabla\phi_i - n_i(n_i\cdot\nabla\phi_i-g_n)
\]

forces normal component of `grad_phi` after solve.

But GMRES solves:

\[
D_h(\sigma G_h\phi)=b_h
\]

Postprocess computes:

\[
J=\sigma[-P(G_h\phi)-\omega A]
\]

where `P` is manual normal projection.

Then original divergence:

\[
D_hJ
\]

These three are not one compatible discrete system:

```text
GMRES sees G_h
E/J postprocess uses P(G_h)
DivJ diagnostic uses D_h
```

So MMS can make `Jn` excellent (MMS checks normal boundary only) but French real Biot destroys bulk `divJ`.

### 5.1 Conclusion

`--phi-boundary-grad-neumann=1` must not be production route.

At most:

```text
diagnostic / MMS-only tool
```

Must not enter:

```text
literature-mode default
French reduced acceptance
subsequent mainline solver
```

---

## 6. How ghost should enter operators, not postprocess

Boundary Neumann:

\[
n\cdot\nabla\phi = g_n
\]

In OPHELIE φ equation:

\[
g_n=-\omega n\cdot A
\]

### 6.1 Weak-form boundary flux route

For diffusion / Poisson:

\[
\nabla\cdot(\sigma\nabla\phi)=b
\]

integral form:

\[
\int_\Omega \nabla\cdot(\sigma\nabla\phi)dV
=\int_{\partial\Omega}\sigma\nabla\phi\cdot n\,dS
\]

Boundary contribution approximately:

\[
S_i\sigma_i g_{n,i}
\]

Current use:

```cpp
rhs_i += (V_i / d_b) * sigma_i * g_n;
```

equivalent to:

\[
S_i\approx \frac{V_i}{d_b}
\]

Theoretically valid, but on French case:

\[
|n\cdot\sigma A|\approx10^{-10}
\]

so `g_n` is tiny; RHS Neumann will not significantly improve residual.

### 6.2 Ghost-in-gradient route

If using ghost, do not fix `grad_phi` after solve; ghost contribution must enter `ComputeOphelieScalarPhiGradientCK` and `DivSigmaGrad` matvec.

For boundary particle `i`, outward normal `n_i`, boundary distance `d_i`, mirror ghost:

\[
x_g=x_i+2d_i n_i
\]

Neumann:

\[
\frac{\phi_g-\phi_i}{2d_i}=g_n
\]

so:

\[
\phi_g=\phi_i+2d_i g_n
\]

Add ghost to gradient/divergence pairwise sums:

\[
G_h\phi_i
=
\sum_{j\in\Omega}(\phi_i-\phi_j)w_{ij}
+ (\phi_i-\phi_g)w_{ig}
\]

Key point:

```text
ghost contribution must enter GMRES every operator apply
cannot fix only in E/J postprocess
```

But on French case boundary `g_n≈0`, ghost mainly completes kernel support, not nonzero flux. Without reasonable ghost volume/kernel weight design, may not immediately improve `phi_eq_res`.

---

## 7. Q3 judgment: why standalone corrected gradient worsens

Current P2 corrected gradient fixes only half:

```text
LHS: Div_uncorrected(σ Grad_corrected(φ))
RHS: Div_uncorrected(σ A)
E/J:  uses Grad_corrected(φ)
DivJ diagnostic: still mainly Div_uncorrected(J)
```

Not compatible discretization.

### 7.1 Corrected operators must be paired

B-matrix corrected gradient:

\[
G_c\phi_i
=
\sum_j(\phi_i-\phi_j)C_i\nabla W_{ij}V_j
\]

Compatible divergence:

\[
D_cF_i
=
\sum_j(F_i-F_j)\cdot C_i\nabla W_{ij}V_j
\]

Unified use:

\[
L_c\phi=D_c(\sigma G_c\phi)
\]

\[
b_c=-\omega D_c(\sigma A)
\]

\[
J=\sigma(-G_c\phi-\omega A)
\]

\[
\nabla\cdot J = D_cJ
\]

### 7.2 Minimum consistent set

If enabling corrected operator, at minimum fix together:

```text
1. grad_phi
2. DivSigmaGrad LHS gradient and divergence
3. DivSigmaA RHS divergence
4. E/J postprocess grad_phi
5. DivJ diagnostic divergence
```

Fixing only one worsens residual or divJ.

### 7.3 CLI recommendation

Do not use `--phi-gradient-correction=1` as mainline production switch.

Add paired switch:

```text
--phi-compatible-correction=0|1
```

Meaning:

```text
G_c / D_c / L_c / b_c / J / divJ switched as a set
```

Old:

```text
--phi-gradient-correction=1
```

diagnostic only.

---

## 8. Corrected gradient matrix inversion risk

Current `InvertOpheliePhiGradLinearCorrectionMatrixCK` logic is risky:

```cpp
const Real determinant = linear_correction_[index_i].determinant();
const Real det_sqr = SMAX(tikhonov_alpha_ - determinant, Real(0));
const Matd inverse = inverseTikhonov(linear_correction_[index_i], SqrtEps);
const Real weight = determinant / (determinant + det_sqr + TinyReal);
linear_correction_[index_i] = weight * inverse + (1 - weight) * Identity;
```

If determinant negative and `tikhonov_alpha_=0`:

```text
det_sqr = -determinant
denominator = determinant + (-determinant) + TinyReal ≈ TinyReal
weight ≈ determinant / TinyReal
```

Huge negative weight may amplify local gradient, especially dangerous near boundary particles.

### 8.1 Recommended conservative version first

```cpp
const Real det = B.determinant();
if (!std::isfinite(det) || std::abs(det) < det_min || det <= 0.0)
{
    C = Identity;
}
else
{
    C = inverseTikhonov(B, eps);
}
```

Also output diagnostics:

```text
B_det_min
B_det_max
B_det_negative_count
B_det_small_count
B_condition_proxy_max
fallback_identity_count
```

If `det_negative_count` is large, correction matrix definition, sign convention, or local particle configuration has issues. Do not run French acceptance with it.

---

## 9. Q4 judgment: French boundary is not current main bottleneck

From current data:

```text
French real Biot boundary |n·σA| ≈ 1e-10
RHS Neumann correction magnitude ≈ 1e-8
RHS Neumann almost no effect on phi_eq_res / divJ
Jn_post_rel baseline ≈ 0.0038, already not large
```

Main bottleneck more likely:

```text
kernel truncation + grad/div incompatible + RHS/LHS range mismatch
```

not missing physical Neumann boundary source.

### 9.1 Next priority

```text
1. Fix CLI and RHS lifecycle
2. Operator determinism / linearity
3. Discrete self-consistency MMS
4. Paired corrected D/G operator
5. Then ghost-in-gradient or virtual shell diagnostic
```

Do not keep optimizing `Jn_post` as primary target.

---

## 10. Q5 acceptance strategy

Current stage primary acceptance:

```text
P_raw ≈ 50 kW
fields finite
divJ_L2_red ≥ 1.25
phi_eq_res does not worsen
Qmax/Qmean reasonable
```

Where:

```text
P_raw       -> literature-mode power scale closure
divJ_L2_red -> whether phi correction truly improves current continuity
phi_eq_res  -> diagnostic and gate; do not optimize alone
Jn_post     -> boundary diagnostic, not main target
```

Do not sacrifice:

```text
divJ_L2_red / literature_passed
```

to achieve:

```text
Jn_post: 1e-3 -> 1e-8
```

---

## 11. SPHinXsys/SYCL code structure assessment

Overall:

```text
Large framework matches SPHinXsys/SYCL CK style;
but new boundary / gradient correction remains experimental, not suitable for default merge into mainline.
```

### 11.1 What works well

Main dynamics use:

```text
LocalDynamics
Interaction<Inner<>>
UpdateKernel / InteractKernel
StateDynamics<ExecutionPolicy, ...>
InteractionDynamicsCK<ExecutionPolicy, ...>
DiscreteVariable
DelegatedData(ex_policy)
```

Matches SPHinXsys CK style.

### 11.2 Structural issues needing remediation

#### 11.2.1 Host/device path naming must be clear

Some functions named like dynamics but actually:

```text
device -> host
host projection
host -> device
```

Recommended split naming:

```cpp
applyOpheliePhiBoundaryGradNeumannProjectionHost(...)
applyOpheliePhiBoundaryGradNeumannProjectionCK(...)
```

French default should not use grad projection; host version diagnostic only.

#### 11.2.2 Correction matrix should not recompute every matvec

B-matrix correction depends only on geometry and neighbors. For static glass particles, precompute once before GMRES:

```text
precomputePhiGradientCorrectionMatrixOnce()
applyCorrectedGradientUsingPrecomputedMatrix()
```

Do not repeat in every `execOphelieScalarPhiGradient()` or GMRES matvec:

```text
compute matrix -> invert matrix -> compute grad
```

High cost and sync/state pollution risk.

#### 11.2.3 UR loader duplication in CMake

Multiple tests manually find:

```cmake
libur_loader.so
```

Short term OK; later extract unified CMake helper.

#### 11.2.4 CLI must refactor

Most urgent engineering item. CLI parse bug may pollute all subsequent comparisons.

---

## 12. Next sprint development plan

Next sprint: do not open A_ind, do not add thermal coupling, do not try new boundary toggles. Three task categories only.

---

### Sprint 1: parameter and RHS lifecycle audit

#### Task 1.1 Fix CLI prefix bug

Files:

```text
electromagnetic_ophelie_cli.h
```

Tasks:

```text
all strncmp + arg offset -> strlen helper
remove magic numbers
add final params print
```

Acceptance:

```text
after passing gmres_restart / max_outer / eq_res_tol / gate on command line, final params display correctly
```

#### Task 1.2 Audit solvePhiImag RHS lifecycle

Files:

```text
electromagnetic_ophelie_phi.h / .hpp
electromagnetic_ophelie_phi_gmres.h
french_literature pipeline solvePhiImag call sites
```

Tasks:

```text
confirm finalized RHS not overwritten inside solve
change interface to solvePhiImagWithCurrentRhs or unify finalize inside solve
```

Acceptance:

```text
RHS checksum: consistent after finalize and before GMRES
with rhs_project_zero_mean, rhs_sum before GMRES actually changes
with one-sided-neumann, rhs_l2 before GMRES changes or confirm ~0 physical reason
```

#### Task 1.3 Add operator determinism / linearity test

Add:

```text
test_3d_ophelie_phi_operator_linear_consistency
```

Acceptance:

```text
repeat_apply_rel < 1e-10
linearity_add_rel < 1e-8
linearity_scale_rel < 1e-8
```

---

### Sprint 2: discrete self-consistency MMS

Add:

```text
test_3d_ophelie_phi_discrete_self_consistency
```

Core flow:

```text
assign phi_exact
compute G_h(phi_exact)
set A = -G_h(phi_exact) / omega
compute RHS = -omega * D_h(sigma A)
compute LHS = D_h(sigma G_h phi_exact)
measure ||LHS - RHS||
```

Acceptance:

```text
discrete_self_consistency_rel < 1e-8 ~ 1e-6
```

If fail, fix LHS/RHS sign, scaling, field sync before Sprint 3.

---

### Sprint 3: paired corrected D/G operator

Prerequisite: Sprint 1/2 pass.

Tasks:

1. Fix corrected matrix inversion fallback;
2. Precompute correction matrix;
3. Add compatible divergence `D_c`;
4. Switch as set:

```text
G_c(phi)
D_c(sigma G_c phi)
D_c(sigma A)
E/J using G_c(phi)
DivJ diagnostic using D_c(J)
```

Add CLI:

```text
--phi-compatible-correction=0|1
```

Do not use old:

```text
--phi-gradient-correction=1
```

as mainline acceptance.

Acceptance:

```text
operator linearity still passes
discrete self-consistency still passes
French case phi_eq_res does not worsen
divJ_L2_red not below baseline
P_raw ≈ 50 kW
```

---

## 13. Routes to abandon or freeze temporarily

Do not invest further for now:

```text
1. grad postprocess Neumann default on
   reason: destroys divJ; diagnostic/MMS only.

2. standalone LHS grad Neumann projection
   reason: MMS worsened; current implementation ineffective.

3. standalone --phi-gradient-correction=1
   reason: gradient only, not divergence/RHS/diagnostic; proven to worsen.

4. literature default Neumann on
   reason: French boundary n·σA≈0; RHS Neumann no benefit.

5. virtual shell
   reason: complex geometry unsuitable as main route; main cause not boundary source.

6. A_ind / thermal coupling
   reason: phi discrete consistency unsettled; A_ind/thermal would stack issues.
```

---

## 14. Direct execution manifest for Cursor

Execute in order:

```text
Step 1:
Fix electromagnetic_ophelie_cli.h.
All prefixes use strlen helper.
Add final params print.

Step 2:
Open solvePhiImag / GMRES implementation.
Confirm whether RHS is re-setup.
If yes, change to solveWithCurrentRhs or unify finalize inside solve.
Add RHS checksum / l2 / sum diagnostics.

Step 3:
Add test_3d_ophelie_phi_operator_linear_consistency.
Validate repeat apply, additivity, scaling.

Step 4:
Add test_3d_ophelie_phi_discrete_self_consistency.
Build discrete-compatible RHS with A = -G_h(phi_exact)/omega.
Require L(phi_exact) + omega D(sigma A) near 0.

Step 5:
Fix corrected gradient matrix inversion.
det<=0 / det too small / nonfinite -> fallback identity.
Output det diagnostics.

Step 6:
Implement paired corrected divergence:
D_c(sigma G_c phi), D_c(sigma A), D_c(J) switched as set.
Add --phi-compatible-correction=0|1.

Step 7:
Rerun French baseline and compatible-correction case:
compare phi_eq_res, divJ_L2_red, Jn_post, P_raw, Qmax/Qmean, runtime.
```

---

## 15. Recommended regression table

After each change run at least:

```text
A. French baseline
--reload=1 --literature-mode --state_recording=0

B. RHS projection
--reload=1 --literature-mode --phi-rhs-project-zero-mean=1 --state_recording=0

C. Neumann RHS only
--reload=1 --literature-mode --phi-boundary-mode=one-sided-neumann --state_recording=0

D. Compatible correction
--reload=1 --literature-mode --phi-compatible-correction=1 --state_recording=0
```

Record:

```text
case_label
n_glass
P_raw
Q_min / Q_mean / Q_max / Q_rms
phi_eq_res_vol
GMRES iter / final residual
divJ_L2_L0
divJ_L2_phi
divJ_L2_red
Jn_pre_rel
Jn_post_rel
rhs_sum
rhs_l2
operator_linearity_passed
discrete_self_consistency_rel
runtime
literature_passed
```

---

## 16. Stage conclusion

Current blocker is not that French OPHELIE-like reduced direction is wrong, but:

```text
phi correction under SPH particle discretization needs a compatible D/G operator system.
```

Current baseline although:

```text
phi_eq_res≈0.49
```

still achieves:

```text
divJ_L2_red≈2.06
P_raw≈50 kW
literature_passed=1
```

So phi correction already gives physical improvement; do not reject the whole chain because of `phi_eq_res` floor alone.

Next stage: avoid scattered toggle trial-and-error; proceed along:

```text
CLI/RHS lifecycle audit
-> operator linearity
-> discrete self-consistency MMS
-> paired corrected D/G operator
```
