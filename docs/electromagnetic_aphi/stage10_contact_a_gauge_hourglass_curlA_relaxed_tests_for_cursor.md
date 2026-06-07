# Stage 10 Contact / A-Gauge / Hourglass / Irregular Particle Test Execution Plan

This document is intended for Cursor execution. It consolidates the latest assessment of the current Stage 10 Contact multi-body A--phi solver, and supplements four new key questions:

1. Whether the current operators may introduce zero-energy modes / hourglass modes;
2. How `B = curl A` should be computed, and whether SPHinXsys kernel/particle form should be used;
3. Whether testing should be extended from regular box particle distributions to sphere/cylinder/relaxed particle distributions;
4. How Stage 10.7 A-divergence penalty should be advanced at this stage.

---

## 0. Overall Conclusion

The current Contact operator and coupled multi-body GMRES main line are largely in place, but that does not mean:

```text
All operators are free of zero-energy modes / hourglass modes;
Coulomb gauge div A = 0 is already satisfied;
B-field / curlA post-processing is complete;
Validation already covers real particle discretization and relaxed particle distributions.
```

Therefore the next phase should advance four diagnostic lines in parallel:

```text
Line A: A-divergence penalty sign/source/PC diagnostics;
Line B: post-solve diagnostics for divA and B=curlA;
Line C: hourglass / near-null-mode / checkerboard perturbation tests;
Line D: relaxed sphere/cylinder particle distribution tests.
```

It is not recommended to proceed directly to:

```text
Contact A-penalty integration;
three-body TEAM7 λ_A sweep;
cold-crucible four-body;
10D thermal coupling;
large helper merge.
```

---

# Part I. Direct Answers to Current Key Questions

## 1. Do all current operators avoid introducing zero-energy modes?

That cannot be claimed.

We have already validated many things:

```text
Contact vs Inner / monolithic equivalence;
fused apply vs debug assemble;
RHS consistency;
coupled GMRES convergence;
Joule/E/J observable alignment;
phi gauge penalty effective against raw field drift.
```

These show:

```text
Operator plumbing is basically correct;
Contact splitting does not obviously break the current scaffold;
the coupled GMRES main line can run end-to-end.
```

But these do not strictly prove:

```text
No hourglass mode;
no unphysical near-null longitudinal mode;
no lattice checkerboard mode;
no discrete zero-energy mode under relaxed particle disorder.
```

In particular, the current issues already expose:

```text
With phi gauge penalty off, Krylov residual is good but raw field drifts severely;
divA_relative is about O(1) in two-body MMS, about O(0.2) globally in three-body;
after adding A-divergence penalty, divA drops but Krylov does not converge.
```

All of this indicates:

```text
The current system still has open issues at the gauge / nullspace / conditioning level.
```

So the conclusion is:

> We cannot claim that all current operators are naturally free of hourglass/zero-energy modes. Dedicated hourglass / near-null-mode diagnostics are needed.

---

## 2. Which hourglass / zero-mode issues should we test?

Three categories are recommended.

### 2.1 Gauge null mode

Including:

```text
constant phi mode;
gradient-like A mode;
divergence-dominated A mode;
bodywise phi offset mode;
multi-body interface longitudinal mode.
```

Already available:

```text
phi_gauge_penalty sweep
divA diagnostic
```

Still needed:

```text
A-side gauge / divA Rayleigh quotient;
gradient-mode response test;
bodywise gauge offset test.
```

### 2.2 Lattice checkerboard / high-frequency mode

On regular box particle distributions, high-frequency alternating modes may be weakly responded to by some SPH pairwise second-derivative operators, forming hourglass-like modes.

Recommended constructions:

```text
checkerboard scalar phi;
checkerboard A_x/A_y/A_z;
random high-frequency perturbation;
smooth low-frequency reference perturbation.
```

Compare:

```text
||K(checkerboard)|| / ||checkerboard||
||K(random_high_freq)|| / ||random_high_freq||
||K(smooth)|| / ||smooth||
Rayleigh quotient <u, K u> / <u, u>
```

If the operator response to checkerboard mode is abnormally small, hourglass risk exists.

### 2.3 Particle-disorder mode

On relaxed / jittered / spherical-geometry particle distributions, irregular kernel support may expose issues invisible on regular lattices.

Recommended tests:

```text
regular box;
jittered box;
relaxed sphere/cylinder;
relaxed two-body sphere/cylinder Contact;
```

Compare:

```text
MMS error;
Krylov convergence;
divA_relative;
checkerboard/random perturbation response;
E/J/Joule observable.
```

---

## 3. The statement "curl B equals A" needs correction

The strict electromagnetic relation is:

```text
B = curl A
```

Not:

```text
curl B = A
```

If discussing Ampère's law, in magnetostatics/frequency domain one may have:

```text
curl H = J
B = μH
```

But in the current A--phi solver, the most direct magnetic-field post-processing should be:

```text
B_real = curl(A_real)
B_imag = curl(A_imag)
```

Then:

```text
H_real = ν B_real
H_imag = ν B_imag
```

To verify the current relation, one can also check:

```text
curl H ≈ J_total
```

But the first step should be `B = curl A`.

---

## 4. How should `curl A` be computed?

Grid finite differences are not recommended, and it should not be written as post-processing unrelated to particle/kernel form.

SPHinXsys kernel/particle gradient form should be used.

### 4.1 Recommended implementation

Reuse the current gradient diagnostic pipeline:

```text
A field
  -> SPH/kernel gradient tensor grad A
  -> curl A from skew part of grad A
```

Let:

```text
gradA[m,n] = ∂A_m / ∂x_n
```

Then:

```text
B_x = ∂A_z/∂y - ∂A_y/∂z = gradA[z,y] - gradA[y,z]
B_y = ∂A_x/∂z - ∂A_z/∂x = gradA[x,z] - gradA[z,x]
B_z = ∂A_y/∂x - ∂A_x/∂y = gradA[y,x] - gradA[x,y]
```

### 4.2 Should existing SPHinXsys curl/vorticity operators be used?

If SPHinXsys already has mature vorticity/curl CK operators, their style and kernel gradient conventions can be referenced, but do not directly black-box the fluid vorticity operator onto the A-field unless confirmed:

```text
It supports arbitrary vector variables;
it supports real/imag variable pairs;
it supports Inner and Contact;
it uses the same correction / kernel gradient convention;
it can be registered and output under the current electromagnetic variable system.
```

A safer approach:

```text
Add AphiComputeCurlACK or AphiMagneticFluxDensityCK;
internally use SPHinXsys kernel gradient/relation data structures;
compute curl separately for A_real/A_imag;
output B_real/B_imag;
optionally output H=νB.
```

### 4.3 Is Inner + Contact needed?

Short term:

```text
Inner diagnostic first;
then Contact diagnostic.
```

Final multi-body magnetic-field post-processing should use:

```text
Inner + Contact gradient of A
```

Otherwise B near interfaces will only see in-body neighbors and miss cross-body support domains.

---

## 5. Should sphere/cylinder/relaxed particle distribution tests be used?

Yes.

A major reason for the concern that the method is "not based on particle discretization, not based on kernel form" may be that most current tests use regular box lattices, which can make it look like a grid/finite-difference scaffold.

Therefore at least one relaxed particle test must be added:

```text
sphere/cylinder geometry;
level-set based particle generation;
particle relaxation;
reload relaxed particles;
kernel-based Inner/Contact relation;
MMS or analytic field;
output convergence / residual / divA / curlA / hourglass metrics.
```

This is important to demonstrate that the method is particle/kernel-based.

Recommended priority:

```text
1. relaxed sphere single-body Inner MMS;
2. relaxed two-body sphere/cylinder Contact interface MMS;
3. relaxed three-body TEAM7-like simplified geometry;
```

Do not start with complex TEAM7 geometry; first build a relaxed sphere/cylinder with controllable exact field.

---

# Part II. Current Assessment of Stage 10.7 A-Divergence Penalty

## 6. Key results of current Stage 10.7

Already implemented:

```text
LhsA += λ_A · grad(div A)
```

In Inner prototype:

```text
λ_A = off:
    converged = true
    outer ≈ 2
    true_rel ≈ 4e-6
    div_A_relative ≈ 0.578

λ_A = 50--100:
    converged = false
    outer = 100
    true_rel ≈ 0.5--0.8
    div_A_relative ≈ 0.03--0.05

λ_A = 300--1000:
    converged = false
    divA rebounds
```

Key assessment:

```text
divA can drop, but GMRES does not converge;
divA on non-converged solutions cannot serve as success evidence.
```

Therefore:

```text
Current A-divergence penalty prototype cannot yet be integrated into Contact;
no three-body λ_A sweep;
not for cold-crucible use.
```

---

## 7. Current `+grad(divA)` may not be a penalty

If the penalty energy is:

```text
E = 1/2 ∫ (div A)^2 dV
```

Its variation corresponds to:

```text
- grad(div A)
```

Because:

```text
<A, -grad(div A)> = ||div A||²
```

is positive under appropriate boundary conditions.

Whereas:

```text
+ grad(div A)
```

is closer to:

```text
curl curl A = grad(div A) - ∇²A
```

It is not a positive divergence penalty and may introduce or amplify nullspace, making GMRES harder to converge.

Therefore the next step is not to keep tuning λ_A, but to run sign/energy diagnostics first.

---

## 8. Cursor must-do: sign / energy diagnostic

Add test:

```text
test_3d_aphi_ck_inner_a_divergence_penalty_sign_energy_diagnostic
```

For a random or manufactured A-field, compute:

```text
D = grad(div A)
```

Output:

```text
inner_A_plus_graddiv_A  = <A, +D>
inner_A_minus_graddiv_A = <A, -D>
div_A_norm2             = ||div A||²
ratio_plus              = inner_A_plus_graddiv_A / div_A_norm2
ratio_minus             = inner_A_minus_graddiv_A / div_A_norm2
```

Expected:

```text
ratio_minus ≈ +1
ratio_plus  ≈ -1
```

If this holds, the true penalty should be changed to:

```text
LhsA -= λ_A · grad(div A)
```

And naming should be clarified:

```text
use_coulomb_divergence_penalty
```

Do not continue calling `+grad(divA)` a penalty.

---

## 9. Current impressed-current source also needs checking

Current Stage 10.7 uses an impressed-current box. If the source current is not divergence-free, it may be incompatible with the Coulomb gauge constraint.

Real coil current should satisfy:

```text
div J_s = 0
```

A local box source with abrupt truncation usually has large:

```text
div J_s
```

Therefore add test:

```text
test_3d_aphi_ck_impressed_current_divergence_diagnostic
```

Output:

```text
source_J_L2
source_div_J_L2
source_div_J_relative
source_div_J_Linf
source_region_count
```

If `source_div_J_relative` is large, the current impressed-current box is not suitable as the core acceptance test for A-divergence penalty.

A solenoidal source test should be added, for example:

```text
J_s = curl C
```

or a closed loop / ring current source.

---

## 10. PC is currently insufficient

Cursor has tried:

```text
d_a += λ_A · laplace_a_diag
```

But `grad(div A)` is a directional longitudinal operator; in Fourier space it is similar to:

```text
k k^T
```

It is not a scalar Laplace diagonal.

Therefore scalar diagonal PC cannot adequately represent:

```text
component coupling;
longitudinal-only stiffness;
cross-component contribution.
```

Next steps:

```text
finite-difference local diagonal diagnostic;
3x3 A-block diagonal approximation;
penalty_to_laplace_diag_ratio output.
```

Do not continue relying only on:

```text
λ_A * laplace_diag
```

to conclude the penalty is infeasible.

---

# Part III. B=curlA and divA Diagnostic Execution Plan

## 11. Add AphiCurlA / B diagnostic

### 11.1 Suggested new classes

```cpp
AphiComputeCurlACK
AphiComputeMagneticFluxDensityCK
AphiComputeMagneticFieldHCK
```

First version can be diagnostic helper only:

```cpp
AphiCurlADiagnosticHelper
```

### 11.2 Input variables

```text
AReal = (A_x_real, A_y_real, A_z_real)
AImag = (A_x_imag, A_y_imag, A_z_imag)
```

### 11.3 Output variables

```text
BReal = curl(AReal)
BImag = curl(AImag)

optional:
HReal = nu * BReal
HImag = nu * BImag
```

### 11.4 Output metrics

```text
B_real_L2
B_imag_L2
B_L2
B_Linf
bodywise_B_L2
regionwise_B_L2

optional:
curl_H_L2
curl_H_minus_J_L2
```

### 11.5 Verification case

Construct manufactured A:

```text
A = (-y/2, x/2, 0)
```

Then:

```text
B = (0, 0, 1)
div A = 0
```

Or:

```text
A = curl Psi
```

Naturally satisfies:

```text
div A = 0
```

Used to verify:

```text
curlA is correct;
divA is small;
kernel gradient convention is correct.
```

---

## 12. Keep divA diagnostic

Already available:

```text
div_A_relative = ||div A|| / ||grad A||
```

Recommended continued output:

```text
global_div_A_relative
bodywise_div_A_relative
regionwise_div_A_relative
plate_div_A_relative
interface_band_div_A_relative
core_div_A_relative
```

Note:

```text
global divA should not be the sole criterion;
three-body TEAM7 should focus on plate/conductor regions.
```

Reason:

```text
air region may occupy large volume;
smooth or weak A field in air affects global ratio;
physically, conductor/plate/melt/crucible regions matter most.
```

---

# Part IV. Hourglass / Near-Null-Mode Test Plan

## 13. Add test: operator mode response diagnostic

Add:

```text
test_3d_aphi_ck_operator_mode_response_diagnostic
```

### 13.1 Mode types

Construct the following perturbations:

```text
smooth low-frequency mode:
    sin(pi x) sin(pi y) sin(pi z)

checkerboard mode:
    (-1)^(i+j+k)

random high-frequency mode:
    random per particle with zero mean

gradient mode:
    A = grad psi

solenoidal mode:
    A = curl Psi
```

### 13.2 For each mode compute

```text
K_mode_norm = ||K(u_mode)||
mode_norm = ||u_mode||
response_ratio = K_mode_norm / mode_norm

Rayleigh:
    rq = <u_mode, K(u_mode)> / <u_mode, u_mode>
```

### 13.3 Assessment

If:

```text
checkerboard response_ratio << smooth response_ratio
```

or:

```text
checkerboard Rayleigh quotient is near 0
```

hourglass / near-null mode may be present.

If:

```text
gradient mode response is abnormally small
```

gauge/nullspace needs further treatment.

---

## 14. Add relaxed particle test

## 14.1 Single-body relaxed sphere MMS

Add:

```text
test_3d_aphi_ck_inner_relaxed_sphere_mms
```

Workflow:

```text
1. Create sphere level set;
2. Generate particles;
3. Use SPHinXsys particle relaxation;
4. Reload relaxed particles;
5. Run A--phi MMS / operator apply diagnostic;
6. Output residual/error/divA/curlA/hourglass mode response.
```

### 14.2 Two-body relaxed sphere/cylinder Contact MMS

Add:

```text
test_3d_aphi_ck_contact_relaxed_two_body_interface_mms
```

Optional geometries:

```text
sphere split by plane;
cylinder split by plane;
inner conductor sphere + outer air shell;
```

Goal:

```text
Prove Contact operators are not only valid on regular box lattices;
prove kernel support / relaxed disorder still works.
```

### 14.3 Acceptance metrics for relaxed geometry

```text
coupled_converged
global_true_rel
max_bodywise_true_rel
interface_band_rel
div_A_relative
curlA_error
left/right field diagnostic
Joule/E/J observable if applicable
```

Do not pursue very strict accuracy at first; focus on:

```text
no blow-up;
reasonable trends;
error converges with dp/refinement;
no obvious hourglass mode.
```

---

# Part V. Answers to Four Questions for Cursor

## Q1: Continue A-penalty or switch to Coulomb projection?

Answer:

```text
Continue A-penalty, but fix sign/source/PC first.
Projection deferred.
```

Reason:

```text
projection requires multi-body Poisson, chi boundary conditions, and phi gauge-invariant update;
too early at this stage.
```

Next steps:

```text
sign/energy diagnostic;
source divergence diagnostic;
corrected-sign penalty;
PC improvement.
```

---

## Q2: How to calibrate λ_A?

Answer:

```text
Do not calibrate from divA on non-converged solutions;
do not directly use absolute λ_A=50~100.
```

Use relative scale:

```text
eta_A = ||λ_A GDiv|| / ||LaplaceA||
```

or:

```text
eta_A = median(|λ_A D_graddiv|) / median(|D_laplace|)
```

And output:

```text
penalty_to_laplace_diag_ratio
```

Only on solutions where:

```text
GMRES converged
true_rel acceptable
```

evaluate divA reduction.

---

## Q3: How to integrate Contact?

Answer:

```text
Do not integrate Contact for now.
Close Inner corrected penalty first.
```

Future Contact should implement:

```text
A field
  -> Inner + Contact gradient
  -> divA owner value
  -> Inner + Contact grad(divA)
  -> A-equation penalty contribution
```

Not only:

```text
per-body independent Inner penalty
```

otherwise interface longitudinal modes are missed.

---

## Q4: How to tier acceptance?

Answer:

```text
Tier-1:
    Krylov convergence + E/J/Joule observable

Tier-2:
    phi raw field + divA + curlA + hourglass diagnostic
```

A-penalty success criteria:

```text
1. GMRES converged;
2. true_rel meets target;
3. divA drops on converged solution;
4. divA drops in conductor/plate regions;
5. E/J/Joule not degraded;
6. λ_A sweep has stable interval;
7. hourglass mode response does not worsen.
```

---

# Part VI. Concrete Execution Order

## 15. Do immediately

```text
1. sign/energy diagnostic for grad(divA)
2. impressed-current source divergence diagnostic
3. B=curlA diagnostic helper
4. hourglass / mode response diagnostic
5. relaxed sphere Inner MMS
```

## 16. Then do

```text
6. corrected-sign A-divergence penalty
7. better PC / graddiv diagonal diagnostic
8. solenoidal source test
9. A-penalty λ_A relative-scale sweep
```

## 17. After that

```text
10. Contact A-penalty pipeline
11. three-body TEAM7 λ_A sweep
12. relaxed two-body Contact geometry
13. cold-crucible four-body EM scaffold
```

## 18. Continue to defer

```text
Coulomb projection
10D thermal coupling
helper large merge
real TEAM7 benchmark finalization
```

---

# Part VII. External Messaging

## 19. Can say

```text
Contact operator and coupled multi-body GMRES have reached engineering-scaffold maturity.
Phi gauge drift has been identified and controlled by phi_gauge_penalty.
divA diagnostic reveals that A-side Coulomb gauge is still open.
The next stage focuses on A-gauge diagnostics, curlA/B post-processing, and relaxed-particle validation.
```

## 20. Cannot say

```text
All operators are free of hourglass modes.
The solver enforces div A = 0.
The current A-divergence penalty is successful.
The current B-field post-processing is complete.
The current tests already prove particle/kernel discretization on relaxed geometries.
```

---

## 21. Final Decision

Current optimal route:

```text
Do not directly integrate Contact A-penalty;
do not switch directly to projection;
do not keep tuning λ_A only;
do not enter cold-crucible four-body yet.

Do first:
    sign/energy;
    source divergence;
    curlA/B diagnostic;
    hourglass mode response;
    relaxed sphere/cylinder particle tests.
```

Core assessment:

> The current divA issue is real, but `LhsA += λ_A grad(divA)` cannot yet be called a successful Coulomb-gauge penalty. Sign, source compatibility, PC, and particle-disorder robustness must be confirmed first. Meanwhile, kernel-based `B=curlA` post-processing and relaxed-geometry tests should be added to address the concern about whether the method is truly based on particle/kernel discretization.
