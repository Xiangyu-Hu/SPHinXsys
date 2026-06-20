# OPHELIE Phi Discretization Operator Follow-On Remediation: Physical Credibility Validation of D and Dc, Biot A Divergence Tests, and Edge-Flux Projection Development Plan

## 0. Purpose of This Document

The current French reduced OPHELIE-like EM chain is already connected:

```text
multiloop Biot → A_src/B_src → φ correction → E/J/JouleHeat → literature 50 kW
```

But a key theoretical and numerical issue remains open:

```text
Uncorrected D and corrected Dc give completely different numerical divergence interpretation for real Biot A_src.
Volume relative difference of D(σA) and Dc(σA) is about 97%; cosine about 0.38.
There is not yet evidence that Dc(A) is closer to continuous physical divergence than D(A).
```

This is not fixable by tuning CLI switches. Next steps must be “operator validation” to answer:

```text
For a real vector field A, which of D(A) and Dc(A) is closer to continuous ∇·A?
For Biot–Savart A_src, is the current residual floor a projection issue of the physical field itself, or numerical noise from SPH discrete divergence?
```

Goals for Cursor are clear:

1. No longer default to `--phi-compatible-correction=1`.
2. No longer substitute discrete self-consistent MMS for continuous vector-field divergence validation.
3. Add continuous vector-field divergence MMS validating D and Dc separately.
4. Add discrete divergence diagnostics for Biot A, including dp convergence, interior/boundary decomposition, axisymmetry decomposition.
5. Design and implement edge-conservative φ projection operator as next-generation mainline candidate.

---

## 1. Why the Previous Issue Was Not Resolved

Previous remediation already proposed vector-field divergence MMS and edge-flux projection, but did not mark them as blocking items; Cursor prioritized clearly encodable work: CLI fixes, RHS lifecycle, operator linearity, discrete self-consistency MMS, paired Gc/Dc.

Cursor completed these usefully:

```text
CLI prefix bug fixed.
RHS overwritten by GMRES entry fixed.
RHS fingerprint proves finalize / gmres_entry / post_solve consistency.
Uncorrected and compatible operators pass linearity tests.
Discrete self-consistent MMS rel ≈ 3e-7.
```

But “discrete self-consistent MMS” only proves code self-consistency, not that D or Dc is more accurate for real Biot A continuous divergence.

Current MMS construction:

```text
Given φ_exact
Use current discrete gradient G_h or G_c to generate A_h = -G_h(φ_exact)/ω
Use same D_h or D_c to generate RHS
Check L_h(φ_exact) and RHS consistency
```

This test passes easily because RHS is manufactured by the same discrete operator. It catches sign errors, field overwrite, host/device sync errors, but not:

```text
For an externally given real vector potential A_src, e.g. Biot–Savart A, is D(A_src) or Dc(A_src) close to continuous ∇·A_src?
```

Therefore the remaining issue is not failure to execute simple tasks, but insufficiently specific, high-priority “continuous vector-field divergence MMS” in prior docs. This document fills that gap.

---

## 2. Current Key Facts and Judgment

### 2.1 Current Facts

French reduced baseline:

```text
phi_eq_res_vol ≈ 0.486
divJ_L2_red ≈ 2.06
P_raw ≈ 50000 W
literature_passed = 1
```

Compatible Gc/Dc:

```text
phi_eq_res_vol ≈ 0.745
divJ_L2_red ≈ 1.34
P_raw ≈ 50000 W
literature_passed = 0
```

D and Dc on real Biot A:

```text
rhs_unc_vs_compatible_vol ≈ 0.97
rhs_cosine_unc_compatible ≈ 0.38
```

B matrix looks healthy:

```text
det_min ≈ 0.075
det_negative = 0
fallback_identity = 0
condition_proxy_max ≈ 4.65
```

### 2.2 Current Judgment

1. `Dc` is not simply a “more accurate D”; it is another divergence discretization.
2. `D(A)` and `Dc(A)` differ greatly on real Biot A, indicating different numerical divergence interpretation for non-potential vector potentials.
3. Without continuous divergence MMS or analytic benchmark, cannot say `Dc(A)` is more physical than `D(A)`.
4. Only production route remains baseline uncorrected D/G.
5. Next step is not more compatible correction tuning but validating D/Dc divergence accuracy on continuous vector fields.

---

## 3. Core Mathematical Issue

French OPHELIE-like φ correction continuous form:

```math
\nabla\cdot[\sigma(-\nabla\phi - \omega A)] = 0
```

Equivalent to:

```math
\nabla\cdot(\sigma\nabla\phi) = -\omega\nabla\cdot(\sigma A)
```

Discrete:

```math
L_h\phi = b_h
```

Where:

```math
L_h = D_h(\sigma G_h), \qquad b_h = -\omega D_h(\sigma A)
```

With corrected operator:

```math
L_c = D_c(\sigma G_c), \qquad b_c = -\omega D_c(\sigma A)
```

Issue:

```text
D_c can be self-consistent for potential fields from G_c(φ).
But Biot A is not from G_c(φ); it is a real vector potential with curl.
Must separately validate whether D_c(A) gives accurate physical divergence for general vector fields.
```

Especially note:

```math
\nabla\times A = B
```

Biot A is typically a solenoidal-curl field. Testing only `A = -G_c(φ)/ω` is insufficient—that A is pure gradient or nearly so.

---

## 4. First Priority: Add Continuous Vector-Field Divergence MMS

### 4.1 New Test Name

Recommend adding:

```text
tests/extra_source_and_tests/3d_examples/test_3d_ophelie_vector_divergence_mms/
```

Or:

```text
test_3d_ophelie_phi_vector_divergence_mms
```

### 4.2 Test Goal

For given analytic vector field `A_exact(x,y,z)`, compute:

```text
D_unc(A_exact)
D_c(A_exact)
```

Compare with continuous analytic divergence:

```text
divA_exact = ∇·A_exact
```

### 4.3 Required Field 1: Constant Field

```math
A(x,y,z) = (1, 0, 0)
```

Continuous divergence:

```math
\nabla\cdot A = 0
```

Purpose: validate no spurious divergence on constant vector field.

Acceptance:

```text
interior L2 error near machine/discrete noise.
boundary error may be larger but should decrease with dp refinement.
Both D and Dc should pass.
```

### 4.4 Required Field 2: Linear Divergent Field

```math
A(x,y,z) = (x, y, z)
```

Continuous divergence:

```math
\nabla\cdot A = 3
```

Purpose: validate first-order linear consistency.

Acceptance:

```text
D_c should theoretically beat D, especially in interior.
If D_c does not beat D on this case, B-matrix correction implementation or sign has issues.
```

### 4.5 Required Field 3: Pure Rotational Divergence-Free Field

```math
A(x,y,z) = (-y, x, 0)
```

Continuous divergence:

```math
\nabla\cdot A = 0
```

Purpose: key case—Biot A has curl; potential-field tests alone are insufficient.

Acceptance:

```text
Both D(A) and Dc(A) should give near-zero divergence in interior.
If Dc gives obvious false divergence while D does not, Dc cannot handle Biot A.
If both give large divergence, current SPH divergence has significant discrete error on rotational fields; phi RHS carries numerical noise.
```

### 4.6 Required Field 4: Quadratic Divergent Field

```math
A(x,y,z) = (x^2, y^2, z^2)
```

Continuous divergence:

```math
\nabla\cdot A = 2x + 2y + 2z
```

Purpose: check higher-order error and dp convergence.

### 4.7 Optional Field 5: Biot-Like Rotational Field with Analytic Zero Divergence

Smooth vector potential, e.g.:

```math
A(x,y,z) = (-y f(r,z), x f(r,z), 0)
```

Where:

```math
r^2=x^2+y^2
```

If `f=f(r,z)`, axisymmetric toroidal vector potential with continuous divergence:

```math
\nabla\cdot A = 0
```

Closer to annular coil A structure in cylindrical region.

Recommend:

```math
f(r,z)=\exp[-\alpha(r^2+z^2)]
```

or:

```math
f(r,z)=1+0.1r^2+0.2z^2
```

Acceptance:

```text
Key case for whether D/Dc suits Biot-like rotational A.
```

---

## 5. Output Metrics Must Split Interior and Boundary

Each vector MMS case must output:

```text
case_name
n_total
n_interior
n_boundary
D_unc_l2_all
D_unc_l2_interior
D_unc_l2_boundary
D_unc_linf_all
D_c_l2_all
D_c_l2_interior
D_c_l2_boundary
D_c_linf_all
D_c_vs_D_l2
D_c_vs_D_cosine
mean_exact_div
mean_D_unc
mean_D_c
global_integral_exact = sum(div_exact * V)
global_integral_unc = sum(D_unc * V)
global_integral_c = sum(D_c * V)
```

Recommended boundary definition:

```text
boundary_shell = distance_to_boundary < 2 dp
interior = distance_to_boundary >= 2 dp
```

Or consistent with current `interior_vol_frac≈0.75` definition.

---

## 6. dp Convergence Required

Not only `dp=0.02`. At minimum:

```text
dp = 0.04, 0.03, 0.02, 0.015
```

If cost allows:

```text
dp = 0.01
```

Each dp outputs CSV:

```text
output/ophelie_vector_divergence_mms.csv
```

CSV columns:

```text
dp, case_name, operator_kind, n_total, n_interior,
l2_all, l2_interior, l2_boundary, linf_all,
global_integral_error, mean_error, cosine_with_exact
```

Acceptance logic:

```text
If D_c clearly beats D on linear field but loses on rotational divergence-free field, D_c suits scalar-gradient consistency only, not Biot A RHS.
If both D and D_c are poor on rotational divergence-free field, main phi RHS noise is vector divergence discretization; need conservative edge-flux projection.
If D is more stable than D_c on Biot-like rotational case, baseline uncorrected remains production.
```

---

## 7. Second Priority: Biot A Discrete Divergence Diagnostics

Add diagnostic function or test:

```text
test_3d_ophelie_biot_a_divergence_diagnostics
```

Or in French reduced add:

```text
--phi-biot-divergence-diagnostics=1
--phi-biot-divergence-csv=PATH
```

### 7.1 Quantities to Output

For real Biot A_src:

```text
D_unc(σA)
D_c(σA)
```

Output:

```text
l2_all
l2_interior
l2_boundary
linf_all
volume_integral
mean
rms
D_unc_vs_D_c_l2
D_unc_vs_D_c_cosine
```

Normalized metric:

```math
\frac{\|D(\sigma A)\|}{\|\sigma A\|/L}
```

Where `L` is glass cylinder radius or height scale.

### 7.2 Axisymmetric Decomposition

If particles are cylindrical glass, decompose A in cylindrical components:

```text
A_r, A_theta, A_z
```

For ideal annular coil, axisymmetric setup, main vector potential is usually azimuthal `A_theta`. Continuous cylindrical divergence:

```math
\nabla\cdot A = \frac{1}{r}\frac{\partial(rA_r)}{\partial r}
              + \frac{1}{r}\frac{\partial A_\theta}{\partial \theta}
              + \frac{\partial A_z}{\partial z}
```

If field is approximately axisymmetric with mainly `A_theta(r,z)`:

```math
\frac{1}{r}\frac{\partial A_\theta}{\partial \theta} \approx 0
```

Continuous divergence should be near 0.

If SPH `D(A)` or `D_c(A)` gives large divergence on approximate pure azimuthal axisymmetric A, RHS contains large numerical divergence, not physical.

### 7.3 Additional Recommended Output

```text
A_r_rms
A_theta_rms
A_z_rms
A_theta_fraction = A_theta_rms / A_rms
axisymmetry_error = rms(A_theta variation over theta bins)
```

If `A_theta_fraction` is high and `axisymmetry_error` is low but `D(A)` is large, the divergence operator creates nonphysical divergence on rotational fields / boundary truncation.

---

## 8. Third Priority: Edge-Conservative φ Projection Operator

If vector divergence MMS shows D/Dc unreliable on Biot-like rotational fields, turn to edge-conservative projection.

Core idea: do not use two-step form:

```text
G(phi) → D(sigma G(phi))
D(sigma A)
```

Define current flux directly on each particle pair so LHS and RHS use the same pair conductance.

### 8.1 Continuous Target

Current continuity:

```math
\nabla\cdot J=0
```

Where:

```math
J = \sigma(-\nabla\phi - \omega A)
```

Edge direction approximation:

```math
\nabla\phi\cdot e_{ij} \approx \frac{\phi_j-\phi_i}{r_{ij}}
```

```math
A\cdot e_{ij} \approx \bar A_{ij}\cdot e_{ij}
```

Edge current flux:

```math
F_{ij} = C_{ij}\left[(\phi_i-\phi_j) - \omega\bar A_{ij}\cdot(x_i-x_j)\right]
```

Or equivalent sign form; confirm via MMS.

Per-particle equation:

```math
\sum_j F_{ij}=0
```

Expand:

```math
\sum_j C_{ij}(\phi_i-\phi_j)
= \omega\sum_j C_{ij}\bar A_{ij}\cdot(x_i-x_j)
```

Exact sign per code current convention; validate with manufactured test.

### 8.2 Pair Conductance Cij

`C_ij` must satisfy:

```text
C_ij = C_ji
C_ij >= 0
Same C_ij for LHS and RHS
```

Recommend reusing existing pairwise Laplace / diffusion weights, e.g.:

```math
C_{ij} = V_j \sigma_{ij} \frac{\partial W_{ij}}{\partial r}\frac{1}{r_{ij}}
```

Or SPHinXsys current diffusion operator form; do not invent unrelated weights.

Conductivity average: harmonic mean:

```math
\sigma_{ij}=\frac{2\sigma_i\sigma_j}{\sigma_i+\sigma_j+\epsilon}
```

Constant σ should reduce to constant σ.

### 8.3 Why Edge-Flux May Resolve phi_eq_res Floor

Two-step method issue:

```text
LHS uses D(G(phi))
RHS uses D(A)
For real Biot A, numerical divergence of D(A) may mainly come from kernel truncation and vector-field discretization error.
```

Edge-flux benefits:

```text
LHS and RHS defined on same edge.
Pair fluxes appear in pairs; stronger conservation.
RHS more likely in range of LHS.
Closer to stable SPH diffusion / Laplace structure.
```

### 8.4 New CLI

Recommend:

```text
--phi-projection-operator=div-grad|compatible-div-grad|edge-flux
```

Where:

```text
div-grad             current baseline D/G
compatible-div-grad  current Gc/Dc, diagnostic only
edge-flux            new pairwise conservative projection
```

Do not multiply confusing combination switches.

### 8.5 Recommended New Files

```text
electromagnetic_ophelie_phi_edge_flux.h
```

Main classes:

```text
ComputeOpheliePhiEdgeFluxLhsCK
ComputeOpheliePhiEdgeFluxRhsFromA CK
ApplyOpheliePhiEdgeFluxOperatorCK
```

Or matrix-free apply first:

```text
apply_edge_flux(phi) -> lhs
setup_edge_flux_rhs(A) -> rhs
```

### 8.6 Edge-Flux Acceptance

MMS first:

1. Discrete self-consistent MMS:

```text
Given phi_exact and A = -edge_gradient(phi_exact)/omega
Require L(phi_exact)-rhs≈0
```

2. Continuous vector-field MMS:

```text
For A=(1,0,0), A=(x,y,z), A=(-y,x,0), Biot-like A_theta field evaluate RHS numerical divergence.
```

3. French reduced comparison:

```text
baseline vs edge-flux
```

Output:

```text
phi_eq_res_vol
divJ_L2_red
Jn_post_rel
P_raw
Qmax
Qmean
Qmax/Qmean
runtime
```

Recommended pass conditions:

```text
P_raw ≈ 50000 W
divJ_L2_red >= baseline or at least >= 1.5
phi_eq_res_vol <= baseline or not clearly worse
Q distribution still strong outside, weak inside, no numerical hotspots
```

---

## 9. On Hodge/Helmholtz Projection

Not recommended as production route for Biot A Hodge/Helmholtz projection now.

Reasons:

```text
φ equation itself performs current-continuity projection.
Hodge projection on A before φ may double-project.
Changes real induction field structure, especially near nonuniform σ or boundaries.
```

Use as diagnostic only:

```text
--phi-a-hodge-diagnostic=1
```

Output:

```text
A_before divergence
A_after divergence
||A_after - A_before|| / ||A_before||
P_total before/after
Q distribution before/after
```

Judge how much of `D(A)` is numerical divergence; do not enter default literature mode.

---

## 10. Current Acceptance Strategy Adjustment

Do not optimize `phi_eq_res` alone. Primary acceptance at current stage:

```text
1. P_raw ≈ target power, default 50 kW
2. divJ_L2_red >= 1.25, preferably >= 1.5 or >= 2.0
3. fields finite, no NaN/Inf
4. JouleHeat spatial distribution reasonable: strong outside, weak inside, strong mid-height
5. current scaling correct: I×2 → A/B/E/J×2, Q×4
```

Keep `phi_eq_res` as diagnostic gate:

```text
phi_eq_res_vol < 0.65 for baseline literature
```

But do not sacrifice:

```text
divJ_L2_red
P_raw
Q distribution
```

to chase:

```text
phi_eq_res < 0.1
```

---

## 11. Cursor Next-Round Execution Manifest

### Step 1: Freeze Baseline Regression

Command:

```bash
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced \
  --reload=1 --literature-mode --state_recording=0
```

Also run:

```bash
--phi-rhs-project-zero-mean=1
--phi-boundary-mode=one-sided-neumann
--target-power=25000
--target-power=50000
--target-power=100000
```

Output CSV:

```text
output/ophelie_french_baseline_regression.csv
```

Columns:

```text
case, dp, n, target_power, P_raw, phi_eq_res, divJ_L2_red, Jn_post_rel,
Qmax, Qmean, Qmax_over_Qmean, maxJ, runtime, literature_passed
```

### Step 2: Clean CLI Recognition Functions

In `isOphelieTestCommandLineOption()`, replace magic-length `strncmp` with helper:

```cpp
ophelieCliStartsWith(arg, "--option=")
```

### Step 3: Add Vector-Field Divergence MMS

Implement `test_3d_ophelie_vector_divergence_mms`.

Must include:

```text
A=(1,0,0)
A=(x,y,z)
A=(-y,x,0)
A=(x^2,y^2,z^2)
A=(-y f(r,z), x f(r,z), 0) optional but strongly recommended
```

Test:

```text
D_unc(A)
D_c(A)
```

Output L2/Linf errors for interior/boundary/all and dp convergence.

### Step 4: Add Biot A Divergence Diagnostic

In French reduced add:

```text
--phi-biot-divergence-diagnostics=1
--phi-biot-divergence-csv=PATH
```

Output:

```text
D_unc(σA), D_c(σA), D_unc_vs_D_c, A_r/A_theta/A_z, axisymmetry_error
```

### Step 5: Based on Steps 3/4, Decide Edge-Flux Projection

If D/Dc unreliable on rotational divergence-free or Biot-like A_theta fields, implement:

```text
--phi-projection-operator=edge-flux
```

Experimental first; do not replace baseline.

---

## 12. Explicit Stop Items

Next round do not:

```text
1. Continue tuning --phi-compatible-correction=1 as production.
2. Put Hodge projection in default route.
3. Default-enable grad Neumann postprocess.
4. Connect A_ind.
5. Connect thermal coupling.
6. Add more random CLI combinations.
```

---

## 13. One-Sentence Conclusion

The unresolved core is not “Cursor failed simple instructions” but insufficient distinction between two MMS types:

```text
Discrete self-consistent MMS: validates code self-consistency.
Continuous vector-field divergence MMS: validates whether D(A) or Dc(A) is close to physical ∇·A.
```

Cursor completed the first. Next must complete the second. Only after proving D or Dc divergence reliability on solenoidal-curl fields, especially Biot-like A_theta, can we choose baseline D/G, switch to Dc/Gc, or develop edge-conservative φ projection.
