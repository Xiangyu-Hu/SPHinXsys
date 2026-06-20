# OPHELIE-like French Reduced Branch: Current Status, GPU Path, A_ind Self-Induction, and Next-Step Remediation Plan

This document guides Cursor to continue `electromagnetic_ophelie` / `test_3d_ophelie_french_reduced`.  
This round update focuses on:

```text
1. Accurate positioning of current French reduced OPHELIE-like branch;
2. How coil excitation is applied today;
3. Whether host-side Krylov + device apply_operator must be GPU-ized now;
4. How to implement device-resident GMRES/Krylov;
5. What remains for full simplified OPHELIE without skin depth, sigma(T), water cooling, two-way thermo-EM;
6. How to add A_ind = K[J_glass] self-induction feedback;
7. Next-step code remediation, tests, docs, and customer demo route.
```

---

## 0. Current Route Positioning

Main line is not full PDE A-phi or TEAM7 benchmark main line, but:

```text
French reduced OPHELIE-like particle induction module
```

Current reduced chain:

```text
code-generated multi-loop coil line source
    -> Biot-Savart A_src / B_src on glass particles
        -> PhiImag correction in glass
            -> E / J / JouleHeat
                -> optional power scaling to 50 kW
                    -> later one-way thermal source
```

Explicitly out of scope:

```text
1. full PDE A-phi solver;
2. TEAM7 as main validation;
3. STL geometry;
4. segmented cold-crucible EM;
5. metal coil/crucible/base skin-depth surface model;
6. sigma(T);
7. water cooling;
8. thermal-EM two-way iteration;
9. real helical coil / lead wires / impedance;
10. full FLUENT-OPHELIE thermo-EM iteration from French literature.
```

To add progressively:

```text
1. device-resident Krylov / GMRES;
2. A_ind = K[J_glass] self-induction feedback;
3. JouleHeat -> thermal one-way coupling.
```

---

## 1. What Is Already Done

### 1.1 Geometry and source

French reduced case uses:

```text
GlassBody:
    analytic cylinder
    relaxed / reload capable
    EM active

CoilSource:
    code-generated multi circular loop line source
    no coil particles required for EM

CoilVisualBody:
    optional visual body only

CrucibleWallVisualBody:
    optional visual / future thermal wall only
    no segmentation
    no EM
```

### 1.2 How coil excitation is applied

No solid coil particles in EM.  
Excitation current is on **code-generated multiloop centerline segments**.

Each loop discretized to line quadrature points:

```cpp
source_pos[k]       // source position x_k
source_moment[k]    // I_loop * dl_k, unit A*m
```

Biot-Savart:

```math
A(x_i) =
\frac{\mu_0}{4\pi}
\sum_k
\frac{I_k \Delta l_k}{|x_i-x_k|}
```

```math
B(x_i) =
\frac{\mu_0}{4\pi}
\sum_k
\frac{I_k \Delta l_k \times (x_i-x_k)}
{|x_i-x_k|^3}
```

That is:

```text
Current source = filament current moment I dl
```

Not:

```text
J dV on coil volume particles
K dS on coil surface particles
unknown coil conductor current
```

For ParaView coil view, add `CoilVisualBody`; it does not participate in Biot or EM solve.

---

## 2. Current Phi / Literature-Mode Status

This round key closure:

```text
1. GradPhi sign validated by test_3d_ophelie_phi_gradient_linear;
2. Do not flip ComputeOphelieScalarPhiGradientCK (phi_i - phi_j);
3. LegacyFlux RHS confirmed ≈ negative of DivSigmaA RHS;
4. DivSigmaGrad LHS + DivSigmaA RHS is current correct main line;
5. reload French reduced literature-mode:
       divJ_red ≈ 1.74
       literature_passed = 1
6. eq_res ≈ 0.57 is discrete solvability / projection floor for real Biot A_src; do not chase <0.1 with more GMRES iterations alone.
```

### 2.1 Defaults to solidify

```text
literature-mode default:
    LHS = DivSigmaGrad
    RHS = DivSigmaA
```

### 2.2 Paths to deprecate

```text
LegacyFlux:
    diagnostic / regression only
    not for real run
```

Must state in README and stdout:

```text
LegacyFlux ≈ -DivSigmaA
Deprecated for production/literature-mode
```

---

## 3. Customer Demo Scope

### 3.1 Can present

Recommended demo command:

```bash
test_3d_ophelie_french_reduced \
  --reload=1 \
  --no-phi \
  --target-power=50000 \
  --state_recording=1
```

Show:

```text
1. analytic cylinder glass geometry;
2. relaxed/reloaded particles;
3. external multiloop coil source;
4. Biot-Savart A/B;
5. induced E/J;
6. JouleHeat field;
7. total JouleHeat scaled to 50 kW;
8. ParaView visualization.
```

### 3.2 Internal presentation

```bash
test_3d_ophelie_french_reduced \
  --reload=1 \
  --literature-mode \
  --state_recording=1
```

Scope:

```text
Reduced literature-mode internal acceptance passed.
```

### 3.3 Cannot claim

Do not say:

```text
1. full OPHELIE reproduced;
2. Jacoutot / French paper fully reproduced;
3. self-induced A_ind has been solved;
4. metal skin-depth / surface conductor model included;
5. sigma(T) and thermal-electromagnetic iteration included;
6. water cooling included.
```

---

## 4. Must Host-Side Krylov + Device apply_operator Change Immediately?

### 4.1 Current status

GMRES/Krylov roughly:

```text
large operator apply:
    device / SYCL CK

Krylov vector operations:
    host-side or partially host-side

small Hessenberg / Givens:
    host
```

French reduced particle count ~:

```text
N_glass ≈ 16k
```

Phi solved once or few times; therefore:

```text
host-side Krylov + device apply_operator does not block customer demo or reduced literature acceptance.
```

### 4.2 Can Krylov move to device now?

Yes. No fundamental blocker.  
Engineering architecture optimization, not physics.

Recommended: not “all control flow on GPU,” but:

```text
Hybrid device-resident GMRES:
    Krylov vectors on device
    dot/norm/axpy/scale/copy on device
    apply_operator on device
    preconditioner on device
    Hessenberg / Givens / stop logic on host
    only scalar reductions copied back to host
```

---

## 5. Device-Resident GMRES/Krylov Implementation Plan

### 5.1 Recommended new files

```text
electromagnetic_ophelie_device_vector_ops.h
electromagnetic_ophelie_device_vector_ops.hpp

electromagnetic_ophelie_device_reductions.h
electromagnetic_ophelie_device_reductions.hpp

electromagnetic_ophelie_device_gmres_workspace.h
electromagnetic_ophelie_device_gmres_workspace.hpp

electromagnetic_ophelie_device_gmres_solver.h
electromagnetic_ophelie_device_gmres_solver.hpp
```

### 5.2 Device vector kernels

Implement:

```cpp
x = 0
y = x
y = alpha * x
y += alpha * x
y = alpha * x + beta * y
w = x - y
x *= alpha
```

Recommended interface:

```cpp
class OphelieDeviceVectorOps
{
public:
    void zero(VariableData<Real> &x);
    void copy(VariableData<Real> &y, VariableData<Real> &x);
    void scale(VariableData<Real> &x, Real alpha);
    void axpy(VariableData<Real> &y, Real alpha, VariableData<Real> &x);
    void xpay(VariableData<Real> &y, VariableData<Real> &x, Real beta);
    void combine(VariableData<Real> &y, Real alpha, VariableData<Real> &x, Real beta);
};
```

Adjust to SPHinXsys variable API.

### 5.3 Device reductions

GMRES needs:

```text
weighted dot:
    dot(x,y) = sum_i x_i*y_i*Vol_i

weighted norm:
    norm(x) = sqrt(dot(x,x))

max abs:
    max_i |x_i|
```

Recommended interface:

```cpp
Real weightedDot(DeviceVector &x, DeviceVector &y, DeviceVector &vol);
Real weightedNorm(DeviceVector &x, DeviceVector &vol);
Real maxAbs(DeviceVector &x);
```

Points:

```text
1. reduction on device;
2. only scalar result to host each time;
3. no full Krylov vector D2H per step.
```

### 5.4 Workspace design

GMRES restart m needs:

```text
V[0..m]      Krylov basis
Z[0..m]      preconditioned basis if right-preconditioned
W            apply_operator temporary
R            residual
TMP          general temporary
DIAG_INV     preconditioner
```

For scalar PhiImag, each vector is one `Real` particle field.

Memory estimate:

```text
m = 30
N = 16066
(31 vectors) * N * 8 bytes ≈ 4 MB
```

Even N=1e6:

```text
31 * 1e6 * 8 bytes ≈ 248 MB
```

Still acceptable.

### 5.5 GMRES main flow

Keep host small matrices:

```text
H[m+1][m]
cs[m]
sn[m]
g[m+1]
```

Large vectors on device.

Pseudo flow:

```text
r = b - A x                         device
beta = ||r||                        scalar to host
v0 = r / beta                       device

for k = 0..m-1:
    z_k = M^{-1} v_k                device
    w = A z_k                       device

    for j = 0..k:
        h[j,k] = dot(w, v_j)        device reduction -> host scalar
        w -= h[j,k] * v_j           device axpy

    h[k+1,k] = norm(w)              device reduction -> host scalar
    v[k+1] = w / h[k+1,k]           device scale

    apply Givens rotations          host small ops
    update residual estimate        host scalar

solve small least squares            host
x += sum_j y_j z_j                  device axpy loop
```

### 5.6 Acceptance

Device GMRES must match host-GMRES:

```text
1. phi_mms:
       same eq_res order
       same divJ_red trend

2. phi_biot_rhs_solvability:
       reload eq_res≈0.574
       divJ_red≈1.74

3. french_reduced literature:
       literature_passed=1
       same P_raw / JouleHeat / divJ_red within tolerance

4. no full-vector copy per Arnoldi step:
       only scalar reductions copied to host
```

### 5.7 Priority

```text
P0: solidify current reduced literature-mode and customer demo
P1: device vector ops + reductions
P2: device-resident GMRES
P3: use device GMRES for self-induction / large 3D
```

Do not let device GMRES block customer demo.

---

## 6. Gap to Simplified OPHELIE

With explicit exclusions:

```text
skin depth
sigma(T)
water cooling
segmented crucible EM
thermal-electromagnetic two-way coupling
```

Core remaining gap:

```text
A_ind = K[J_glass] self-induced vector potential
```

Current:

```math
A \approx A_{src}
```

More complete simplified OPHELIE-like:

```math
A = A_{src} + A_{ind}
```

Where:

```math
A_{ind}(x_i)
=
\frac{\mu_0}{4\pi}
\sum_{j \in glass}
\frac{J_j V_j}{|x_i - x_j|}
```

Then:

```math
J = -\sigma(\nabla\phi + i\omega A)
```

With:

```math
\nabla \cdot J = 0
```

---

## 7. A_ind = K[J_glass] Implementation Plan

### 7.1 Fix/confirm complex channel

Frequency-domain variables:

```text
JReal, JImag
AIndReal, AIndImag
BIndReal, BIndImag
```

Biot-Savart is linear real operator:

```text
JReal -> AIndReal, BIndReal
JImag -> AIndImag, BIndImag
```

Do not integrate `JImag` into `AIndReal`.

### 7.2 Source current moment from glass

Per glass particle:

```text
moment_real_j = JReal_j * Vol_j
moment_imag_j = JImag_j * Vol_j
```

Units:

```text
J [A/m^2] * Vol [m^3] = A*m
```

### 7.3 Glass-to-glass Biot

Recommended files:

```text
electromagnetic_ophelie_glass_self_induction.h
electromagnetic_ophelie_glass_self_induction.hpp
```

Core operator:

```cpp
ComputeOphelieGlassSelfInducedBiotSavartCK
```

Or host reference first, then CK. Formula:

```cpp
for each target glass i:
    AIndReal = 0
    AIndImag = 0
    BIndReal = 0
    BIndImag = 0

    for each source glass j, j != i:
        r = x_i - x_j
        inv_r = 1 / sqrt(r2 + eps2)
        inv_r3 = inv_r / r2

        mr = JReal[j] * Vol[j]
        mi = JImag[j] * Vol[j]

        AIndReal += coeff * mr * inv_r
        AIndImag += coeff * mi * inv_r

        BIndReal += coeff * cross(mr, r) * inv_r3
        BIndImag += coeff * cross(mi, r) * inv_r3
```

Note:

```text
B direction is still current_moment cross r
```

### 7.4 Combine fields

Add:

```cpp
ASrcReal = ACoilReal + AIndReal
ASrcImag = ACoilImag + AIndImag
BSrcReal = BCoilReal + BIndReal
BSrcImag = BCoilImag + BIndImag
```

Or rename:

```text
ATotalReal / ATotalImag
BTotalReal / BTotalImag
```

If code uses `ASrc` heavily, keep `ASrc = total A` and output:

```text
ACoilReal
AIndReal
ASrcReal
```

### 7.5 Picard self-consistent iteration

Minimal iteration:

```text
Precompute ACoil, BCoil.

Initialize:
    AInd = 0
    ATotal = ACoil
    solve Phi/J using ATotal

for k = 0..max_iter:
    compute AInd_new = K[J^k]
    ATotal_new = ACoil + AInd_new

    solve Phi^{k+1} using ATotal_new
    compute E/J/Joule^{k+1}

    under-relax:
        J^{k+1} = (1-alpha) J^k + alpha J_new
        AInd^{k+1} = (1-alpha) AInd^k + alpha AInd_new

    check:
        rel_J = ||J^{k+1} - J^k|| / ||J^k||
        rel_Q = |P^{k+1} - P^k| / P^k
```

Suggested defaults:

```text
--self-induction=1
--self-induction-max-iter=10
--self-induction-relax=0.3
--self-induction-tol=1e-3
```

### 7.6 One-way induced field first

Before full Picard, one-way diagnostic:

```text
1. Use ACoil only to solve Phi/J.
2. Compute AInd/BInd from J.
3. Output AInd/BInd magnitude relative to ACoil/BCoil.
4. Do not feed AInd back yet.
```

Answers:

```text
Is the glass-induced field large enough to matter?
```

Metrics:

```text
||AInd|| / ||ACoil||
||BInd|| / ||BCoil||
P_Joule with ACoil only
P_Joule with ACoil + AInd one-way
```

### 7.7 Complexity warning

Glass-to-glass all-to-all:

```text
O(N_glass^2)
```

For N≈16k:

```text
~2.6e8 interactions per self-induction Biot apply
```

Feasible on GPU for diagnostics; too slow on CPU/host.  
Therefore:

```text
A_ind must be GPU/SYCL CK for realistic runs.
```

Bigger performance issue than Krylov dot/axpy.

---

## 8. A_ind and Device Krylov Relationship

```text
device Krylov:
    improves Phi solve performance
    mostly O(N * m^2) vector ops
    engineering architecture improvement

A_ind:
    adds missing simplified OPHELIE physics
    O(N^2) long-range interaction
    much more expensive physically
```

Recommended order:

```text
1. solidify current reduced chain;
2. implement device vector ops / reductions;
3. implement device GMRES;
4. implement one-way A_ind GPU Biot;
5. implement Picard self-induction.
```

If time tight:

```text
Priority one-way A_ind GPU Biot
```

More directly completes simplified OPHELIE physics chain.

---

## 9. JouleHeat -> Thermal One-Way Coupling

User excludes two-way thermo-EM; next thermal is one-way:

```text
EM solved once
JouleHeat fixed or periodically updated
Thermal solver uses Q as source
No sigma(T) feedback
```

Minimal validation:

```text
test_3d_ophelie_joule_to_heat_one_way
```

Equation:

```math
\frac{dT}{dt} = \frac{Q}{\rho c_p}
```

For uniform Q:

```math
\Delta T = \frac{Q \Delta t}{\rho c_p}
```

Outputs:

```text
HeatingRate = JouleHeat / (rho * cp)
Temperature
Total internal energy increase
Compare with integral(Q dt)
```

After demo fields stable; before full self-induction if customer wants thermal visualization.

---

## 10. Cursor Next-Step Task Manifest

### P0: Docs and default route solidification

```text
1. Update OPHELIE_PHI_TEST_UPDATE_AFTER_GRADIENT_RESULT.md;
2. Update OPHELIE_PHI_OPERATOR_AUDIT.md;
3. Update FRENCH_LITERATURE_MODE.md;
4. State clearly:
   DivSigmaGrad + DivSigmaA = literature default;
   LegacyFlux = deprecated diagnostic only;
   eq_res_gate = 0.65;
   divJ_red_gate = 1.25;
   literature_passed=1 means reduced internal acceptance, not full OPHELIE.
```

### P1: GPU Krylov infrastructure — **closed (2026-06-02)**

```text
Done (scaffold):
  electromagnetic_ophelie_phi_device_vector_ops.h
  test_3d_ophelie_phi_device_vector_ops --reload=1 --gmres-parity + --krylov-arnoldi
  biot --gmres-device-ops-parity

Done (SphinxSys SYCL native, 2026-06-02):
  electromagnetic_ophelie_phi_krylov_ck.h
    ophelieVolWeightedDotParticleFields → particle_reduce + LoopRangeCK<ParallelDevicePolicy, SPHBody>
    OphelieKrylovAxpyFieldsCK / OphelieKrylovSubtractScaledFieldsCK / OphelieKrylovCopyFieldCK + StateDynamics
    PhiKrylovWorkspace/Basis/ScratchA/B on glass particles (register_fields)
  Krylov USM shared storage; dot/norm stage via hostAssign + syncVariableToDevice + particle_reduce
  subtract/scale/copy remain USM parallel_for (shared memory; validated)
  reload parity: device_ops ~10 s, device_krylov ~28 s @40 outer; rel_err ~7e-6

Optional later:
  subtract/normalize fully via CK without USM parallel_for
  trim hostAssign staging per dot (persistent particle-field basis columns)
  french_reduced end-to-end with --phi-gmres-device-ops=1
```

French reduced `--literature-mode --self-induction --reload=1` (2026-06-02):
  P calibrate → 50 kW OK; self_ind_J_rel≈1.9e-5; A_ind/A_coil≈0.371;
  literature_passed=0 (by design); passed=1; divJ_gate=0 (Picard path skips L0 baseline).

### P2: A_ind one-way diagnostic — **implemented (2026-06-02)**

```text
1. Fixed ComputeOphelieGlassSelfInducedBiotSavartCK:
   JReal -> AIndReal/BIndReal, JImag -> AIndImag/BIndImag.
2. New test: test_3d_ophelie_french_aind_diagnostic
   helper: electromagnetic_ophelie_aind_diagnostic.h
3. Reload n=16066 sample (coil-only phi, no feedback):
   A_ind/A_coil ≈ 0.37, B_ind/B_coil ≈ 0.28
   max_J_imag ≈ 82 A/m^2, max_J_real = 0 (harmonic formulation)
```

→ **Picard self-induction is justified** (A_ind not negligible); next P3.

### P3: A_ind Picard self-consistency — **implemented & reload-validated (2026-06-02)**

```text
1. runFrenchReducedSelfInductionPicard + test_3d_ophelie_french_self_induction_picard
2. Lattice n=262: J_rel≈1.5e-3 in 2 outer iters; A_ind/A_coil≈0.43
3. Reload n=16066 (--reload=1, max_iter=5, relax=0.3):
   self_induction_iters=2, final_J_rel≈4.2e-6, converged=1, passed=1
   phi_eq_res_vol≈0.574 (GMRES floor, unchanged)
   P_joule_W≈4.723, A_ind/A_coil≈0.3708, B_ind/B_coil≈0.279
4. Reload Picard metrics ≈ one-way diagnostic (same J fixed point):
   feeding K[J] into A_src barely moves J at this φ accuracy → not a bug;
   A_ind vol-norm ratio still ~37% of coil (physically non-negligible).
5. --self-induction wired on test_3d_ophelie_french_reduced (experimental; not literature_passed).
```

### P4: JouleHeat -> thermal one-way

```text
1. Add HeatingRate = Q/(rho*cp).
2. Add one-step thermal verification.
3. Add optional fixed-Q heat diffusion.
4. Prepare customer visual outputs.
```

---

## 11. Recommended Implementation Order

```text
Step 1:
    Freeze current French reduced literature-mode and docs.

Step 2:
    Prepare customer demo output:
        no-phi/scaled and literature-mode/internal.

Step 3:
    Implement device vector ops and reductions.

Step 4:
    Implement device-resident GMRES and verify against host-GMRES.

Step 5:
    Implement one-way A_ind GPU Biot.

Step 6:
    Implement Picard A_ind self-consistency.

Step 7:
    Implement one-way thermal source.
```

If customer demo urgent:

```text
Step 1 -> Step 2 -> Step 7
```

If OPHELIE physics completeness urgent:

```text
Step 1 -> Step 5 -> Step 6
```

If SYCL architecture completeness urgent:

```text
Step 1 -> Step 3 -> Step 4
```

---

## 12. One-Line Instruction for Cursor

```text
French reduced OPHELIE-like branch already has source-field + phi + JouleHeat reduced literature internal acceptance. Next: do not revisit GradPhi/LegacyFlux or TEAM7/STL. Solidify docs and tests first, then parallel two enhancement lines: (1) device-resident GMRES/Krylov—GPU Krylov vector ops and reductions; (2) A_ind=K[J_glass] one-way diagnostic and Picard self-consistency—the key gap to simplified OPHELIE EM chain without skin-depth, sigma(T), water cooling, and two-way thermo-EM.
```

---

## 13. Final Goals

Short term:

```text
Customer demo:
    multiloop coil source
    glass A/B/E/J/JouleHeat
    50 kW scaling
    ParaView output
```

Midterm reduced OPHELIE:

```text
A = A_src + A_ind
Phi/J self-consistent
JouleHeat output
```

Engineering optimization:

```text
device-resident GMRES
GPU A_ind long-range Biot
one-way thermal coupling
```

Long term:

```text
passive crucible EM
surface/skin-depth
sigma(T)
thermal-electromagnetic coupling
real French geometry
```
