# Stage 10.11 Execution Plan: Boundary divA, B=curlA, E/J/H Validation, Projection, and eta_A Explanation

This document guides Cursor to continue Stage 10 A--phi solver mainline. It answers key questions:

1. Large boundary divA residual — can SPH boundary treatment help?
2. Are A/phi equations accurate but div A recovery fails at boundary?
3. Does divA issue mean A field itself is wrong?
4. Must we validate `B = curl A`?
5. Must we validate `H` or `curl H`?
6. Must we validate electric field `E`?
7. Do these quantities have theoretical solutions?
8. Do projection now?
9. What is `eta_A`, what does it do?
10. How should Cursor execute next?

---

## 0. Overall Conclusion

Stage 10.10 results show:

```text
1. A--phi solver mainline and pairwise A-divergence penalty route are promising;
2. core-region pairwise divA already very good for Az2D / CrossSine3D divergence-free fields;
3. current large global divA mainly from boundary layer and unsuitable global normalization;
4. this cannot simply equate to "entire A field wrong";
5. but cannot ignore boundary divA; future Contact interface / cold-crucible walls have similar incomplete support;
6. next must validate B=curlA, E/J/Joule, H=nu B / curlH, and build exact/reference;
7. projection not next mainline for now; do boundary-aware diagnostic and SPH boundary treatment first;
8. eta_A is A-divergence penalty relative strength, not physical parameter, not equal to lambda_phi.
```

Current mainline recommendation:

```text
Keep Contact baseline:
    lambda_phi = 100
    use_a_divergence_penalty = false

Continue Inner A-penalty research:
    eta_A = 0.1 primary research candidate
    eta_A = 0.2 optional
    eta_A = 0.3 upper research bound
```

Do not now directly:

```text
1. Enable Contact A-penalty;
2. Do cold-crucible four-body;
3. Do 10D thermal coupling;
4. Do projection;
5. Claim Coulomb gauge globally satisfied.
```

---

# Part I. Can Large Boundary divA Borrow SPH Boundary Treatment?

## 1. Yes, and Should Borrow SPH Boundary Treatment Ideas

Current large boundary divA likely from:

```text
incomplete kernel support;
missing particles outside boundary;
gradient / divergence operator consistency reduced at boundary;
unsuitable normalization denominator;
A field at box boundary not matching divA operator boundary assumptions.
```

Similar to common SPH boundary issues: near free surface, walls, interfaces, kernel support truncated; interior formula directly applied gives large gradient/divergence error.

Can borrow following SPH treatments.

---

## 2. Recommended Boundary Diagnostics and Corrections

### 2.1 Do Not Fix Directly First; Quantify Boundary Layer

Add/keep output:

```text
core_mask
boundary_mask
distance_to_boundary
neighbor_count
kernel_sum
linear_correction_condition
support_deficiency
DivAAbs
DivAReal
DivAImag
```

Partition metrics:

```text
core_divA_rel
boundary_divA_rel
boundary_divA_energy_fraction
core_particle_count
boundary_particle_count
```

Goal:

```text
Confirm divA concentrates in boundary layer only;
Confirm core already near theoretical 0;
Confirm boundary divA varies with boundary layer thickness.
```

### 2.2 Boundary Exclusion / Core-Only Gate

For manufactured validation, use:

```text
core-only hard gate
boundary diagnostic warning
```

Reason:

```text
Many MMS exact fields lack matching SPH boundary conditions;
Incomplete support at boundary causes natural error;
Using global divA as hard gate wrongly rejects correct interior operator.
```

But boundary still matters separately as open issue.

### 2.3 Ghost / Dummy / Buffer Particles

Can borrow SPHinXsys wall/dummy particle idea; supplement outside boundary:

```text
ghost particles
dummy particles
buffer particles
mirror particles
```

To restore kernel support.

For divA/curlA diagnostic, try:

```text
1. Add ghost/buffer layer outside box;
2. Assign ghost A from analytic exact field;
3. Recompute divA/curlA;
4. Compare whether boundary divA decreases.
```

If boundary divA drops significantly after ghost, issue is mainly support deficiency, not A field error.

### 2.4 Boundary Correction / Renormalization

If linear correction matrix exists, test:

```text
uncorrected pairwise divA
linear-corrected gradient trace divA
renormalized kernel gradient divA
B-corrected divA
```

But note:

```text
diagnostic can use corrected;
penalty apply and PC must be same family.
```

Cannot have:

```text
penalty apply uses corrected;
PC uses pairwise;
```

Or Stage 10.7 non-convergence returns.

### 2.5 Manufactured Boundary Condition Matching

Many exact fields do not satisfy natural boundary assumptions at box boundary. Construct MMS better suited to boundaries:

```text
A field decays to 0 near boundary;
or use periodic-like field;
or validate only in core region.
```

E.g. envelope:

```text
w(x,y,z) = x(1-x)y(1-y)z(1-z)
```

Construct:

```text
A = curl(Psi), Psi = w * smooth function
```

So A smoother/nearer zero at boundary; helps distinguish boundary support artifact from operator error.

---

## 3. Is It "A/phi Solve Accurate, but divA Recovery Fails at Boundary"?

More accurate statement:

```text
A/phi equation residual and MMS/observable validated in several cases;
post-solve divA diagnostic good in core;
global divA large mainly from boundary layer and normalization;
therefore cannot simply say A/phi solution wrong, nor divA diagnostic entirely wrong.
```

Current issue more like:

```text
boundary divA diagnostic / gauge consistency not fully closed.
```

Not:

```text
entire A/phi solver wrong.
```

But cannot ignore boundary because real problems have important boundaries/interfaces.

---

## 4. Does divA Issue Prove A Field Itself Wrong?

Not necessarily; depends on case.

### 4.1 If divA Mainly at Boundary, Core Very Small

Usually means:

```text
A field in interior may be correct;
boundary diagnostic affected by support deficiency / boundary condition;
B/E/J/Joule may still be reasonable.
```

Cannot directly judge entire A field wrong.

### 4.2 If divA Large in Core Too

More serious; may mean:

```text
A contains non-physical longitudinal component;
Coulomb gauge not satisfied;
Laplace-form A solution may not match curl-curl physics;
E/J/Joule may be polluted.
```

### 4.3 Cannot Judge A Physical Accuracy from divA Alone

Must also check:

```text
B = curl A
E = -j omega A - grad(phi)
J = sigma E
Joule = 0.5 sigma |E|^2
H = nu B
curl H ≈ J_total
```

If:

```text
divA boundary large, but B/E/J/Joule vs exact/reference good;
```

Mainly boundary gauge diagnostic artifact.

If:

```text
divA large, and B/E/J/Joule also wrong;
```

A field physics problematic; need further correction.

---

# Part II. Must We Validate B, H, E, J, Joule?

## 5. Must Validate B = curl A

Yes.

Reason:

```text
A itself has gauge DOF;
B = curl A is physical magnetic flux density, insensitive to gradient gauge mode;
If A longitudinal component problematic, B may still be correct or expose transverse error.
```

Therefore `B=curlA` is important post-processing to judge A field physical meaning.

### 5.1 How to Compute B

Use SPHinXsys kernel-based gradient, not grid finite differences.

First compute:

```text
gradA[m,n] = ∂A_m / ∂x_n
```

Then:

```text
B_x = ∂A_z/∂y - ∂A_y/∂z
B_y = ∂A_x/∂z - ∂A_z/∂x
B_z = ∂A_y/∂x - ∂A_x/∂y
```

Do separately for real / imag:

```text
B_real = curl(A_real)
B_imag = curl(A_imag)
```

### 5.2 Recommended Implementation

Add diagnostic helper:

```text
AphiCurlADiagnosticHelper
AphiComputeCurlACK
AphiMagneticFluxDensityDiagnostic
```

Output:

```text
BReal
BImag
BAbs
B_L2_error_vs_exact
B_core_error
B_boundary_error
```

---

## 6. Must We Validate H?

Yes, but after B.

In current eddy-current A--phi setting:

```text
B = curl A
H = nu B
```

Here `nu` is reluctivity, `nu = 1 / mu`.

If `nu` constant, H is B scaling; after B validated, H usually follows.

If `nu` spatially varying or multi-material:

```text
H = nu(x) B
```

Must validate:

```text
H correct near material interface;
H_tangential / B_normal interface physics reasonable;
```

### 6.1 Do Not Say "curlB Gives H"

Strict relation is not:

```text
curl B = H
```

But:

```text
H = nu B
curl H = J_total
```

Where `J_total` includes:

```text
source current J_s
induced current J_ind = sigma E
```

Exact form depends on current frequency-domain equations and whether displacement current ignored.

### 6.2 Optional curlH Validation Later

With `H`, can validate:

```text
curl H ≈ J_s + J_ind
```

Or in some regions:

```text
curl H ≈ J
```

Closer to Maxwell equations than H alone.

---

## 7. Must We Validate E?

Yes.

Electric field is direct source of Joule heating:

```text
E = -j omega A - grad(phi)
```

In real/imag split, roughly:

```text
E_real = omega * A_imag - grad(phi_real)
E_imag = -omega * A_real - grad(phi_imag)
```

Current Joule uses:

```text
q = 0.5 sigma (|E_real|^2 + |E_imag|^2)
```

If E wrong, J and Joule wrong.

### 7.1 E Validation Metrics

Do not look at single real/imag component only. Use uniformly:

```text
E_combined_L2 = sqrt( ∫ (|E_real|² + |E_imag|²) dV )
```

Error:

```text
E_combined_error_vs_exact
```

Similarly:

```text
J_combined_error_vs_exact
Joule_error_vs_exact
```

### 7.2 Why E_combined

In phasor problems, some exact fields may have only real or imag component. Single near-zero component gives false relative error explosion.

Therefore use combined norm uniformly.

---

## 8. Do These Quantities Have Theoretical Solutions?

Yes; manufactured solutions can construct them.

### 8.1 Az2D Exact Field

Set:

```text
A = (0, 0, sin(pi x) sin(pi y))
phi = 0
```

Then:

```text
div A = 0
```

Magnetic field:

```text
B_x = ∂A_z/∂y = pi sin(pi x) cos(pi y)
B_y = -∂A_z/∂x = -pi cos(pi x) sin(pi y)
B_z = 0
```

Electric field:

```text
E = -j omega A
```

If A real phasor:

```text
E_real = 0
E_imag = -omega A_real
```

Current:

```text
J = sigma E
```

Joule:

```text
q = 0.5 sigma omega^2 |A|^2
```

H:

```text
H = nu B
```

If nu constant, direct scaling.

### 8.2 CrossSine3D Exact Field

Set:

```text
A_x = sin(pi y) sin(pi z)
A_y = sin(pi z) sin(pi x)
A_z = sin(pi x) sin(pi y)
phi = 0
```

Then:

```text
divA = 0
```

B analytically:

```text
B_x = ∂A_z/∂y - ∂A_y/∂z
    = pi sin(pi x) cos(pi y) - pi cos(pi z) sin(pi x)

B_y = ∂A_x/∂z - ∂A_z/∂x
    = pi sin(pi y) cos(pi z) - pi cos(pi x) sin(pi y)

B_z = ∂A_y/∂x - ∂A_x/∂y
    = pi sin(pi z) cos(pi x) - pi cos(pi y) sin(pi z)
```

E/J/Joule likewise analytic.

### 8.3 Use These Exact Fields For

Validate:

```text
A error
divA error
B=curlA error
E error
J error
Joule error
H error
curlH vs J consistency
```

These are true absolute validation.

---

# Part III. Try Projection Now?

## 9. What Is Projection

Coulomb projection similar to incompressible-flow pressure projection:

```text
solve Laplace(chi) = div A
A_new = A_old - grad chi
```

Theoretically:

```text
div A_new = div A_old - Laplace chi = 0
```

In frequency-domain A--phi, for gauge consistency may need:

```text
phi_new = phi_old + j omega chi
```

Exact sign must align with current phasor convention.

---

## 10. Not Recommended to Do Projection Immediately

Reasons:

```text
1. Projection needs extra Poisson solve;
2. Multi-body Contact projection needs Contact Poisson;
3. chi boundary conditions unclear;
4. phi update must preserve E;
5. current boundary divA mainly support artifact; projection may over-correct;
6. pairwise penalty not finished boundary/exact observable validation;
7. projection now adds more unknowns.
```

Current recommendation:

```text
Projection: design document only, not code mainline.
```

Only enter projection prototype when:

```text
1. Inner pairwise/B-corrected divA diagnostic both show unacceptable core divA;
2. A-divergence penalty cannot lower divA without breaking E/J/Joule;
3. divA still high after boundary treatment;
4. Contact or real TEAM7 must strictly satisfy divA=0.
```

---

# Part IV. What Is eta_A?

## 11. Definition of eta_A

`eta_A` is relative strength parameter for A-divergence penalty.

Not a physical parameter added directly to equations; used to auto-compute actual `lambda_A` from current discrete operator scale.

Current idea:

```text
eta_A = median(|lambda_A * D_graddiv|) / median(|D_laplace|)
```

i.e.:

```text
typical diagonal strength of A-divergence penalty / typical diagonal strength of A-Laplace
```

Therefore:

```text
eta_A = 0.1
```

Roughly means:

```text
A-divergence penalty typical strength about 10% of A-Laplace.
```

### 11.1 Relation Between eta_A and lambda_A

```text
lambda_A = eta_A * median(|D_laplace|) / median(|D_graddiv|)
```

Code may invert via helper.

### 11.2 Why Not Absolute lambda_A

Because absolute `lambda_A` depends on:

```text
dp
kernel
particle volume
operator scaling
material parameters
dimension
```

Not directly comparable across cases.

So say:

```text
eta_A = 0.1
```

Not:

```text
lambda_A = 100
```

### 11.3 eta_A Unrelated to lambda_phi

`lambda_phi` is:

```text
zero-order gauge penalty in phi equation:
LhsPhi += lambda_phi * phi
```

Locks phi constant gauge.

`eta_A/lambda_A` is:

```text
second-derivative divergence penalty in A equation:
A equation += lambda_A * [-grad(divA)]
```

Controls A divergence.

No direct relation between the two.

---

# Part V. Cursor Next Execution Plan

## 12. Stage 10.11-A: Boundary divA Specialization

Add/extend test:

```text
test_3d_aphi_ck_pairwise_diva_boundary_diagnostic
```

### 12.1 Output

```text
boundary_width
core_shell
neighbor_count
kernel_sum
support_deficiency
linear_correction_condition
pairwise_divA_abs
Bcorrected_divA_abs
core_divA_rel_gradDen
boundary_divA_rel_gradDen
boundary_energy_fraction
```

### 12.2 VTP Output

```text
DivAReal
DivAImag
DivAAbs
CoreMask
BoundaryMask
NeighborCount
KernelSum
SupportDeficiency
```

### 12.3 Test Content

For exact fields:

```text
Az2D
CrossSine3D
Sinusoidal3DCurlPsi
```

Do:

```text
core_shell sweep
boundary_width sweep
dp sweep
```

### 12.4 Judgment

If boundary divA drops significantly with ghost/buffer/core exclusion:

```text
boundary support artifact confirmed
```

If not:

```text
check div operator or boundary condition.
```

---

## 13. Stage 10.11-B: B/E/J/Joule/H Exact Observable Validation

Extend:

```text
test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms
```

### 13.1 Must Output

```text
A_error
divA_core_error
divA_boundary_error
B_error
E_combined_error
J_combined_error
Joule_error
H_error
```

### 13.2 Optional Output

```text
curlH_error_vs_J
```

If implementing curlH, note:

```text
curlH = J_total
```

Clarify what `J_total` includes per current case.

### 13.3 Exact Fields

At minimum support:

```text
Az2D
CrossSine3D
```

---

## 14. Stage 10.11-C: SPH Boundary Treatment Attempts

Try in priority order:

### 14.1 Core-Only Validation

Hard gate limited to core first.

### 14.2 Ghost/Buffer Analytic Particles

In MMS case, add analytic ghost particles outside boundary with exact A/phi for divA/curlA diagnostic.

### 14.3 Mirror Particles

For simple box, try mirror particles and mirror A values.

### 14.4 Boundary Correction Factor

Renormalization from kernel sum / support deficiency.

### 14.5 Do Not Change Main Solver First

These as diagnostic / post-processing tests first; do not directly change production A--phi solve.

---

## 15. Stage 10.11-D: Projection Design Document

Add document:

```text
CURSOR_APHI_STAGE10_PROJECTION_DESIGN_NOTES.md
```

Design only, no implementation.

Must answer:

```text
1. chi equation: which operator for Laplace chi = divA?
2. What are chi boundary conditions?
3. A_new = A - grad chi uses which gradient?
4. Update phi as phi_new = phi + j omega chi?
5. Does E remain unchanged?
6. How for Contact multi-body?
7. Projection post-solve correction or embedded in iteration?
```

---

# Part VI. Acceptance and External Wording

## 16. Can Say Now

```text
The A-divergence issue is now localized: the pairwise divA is accurate in the core region for non-degenerate divergence-free MMS fields, while the global divA metric is dominated by boundary-layer contributions and normalization effects.
```

Summary:

```text
A divergence issue localized: for non-degenerate divergence-free MMS, pairwise divA correct in core; global divA mainly from boundary layer and normalization.
```

## 17. Cannot Say Now

```text
The Coulomb gauge is globally enforced.
Boundary divA is solved.
A-penalty is production ready.
Contact A-penalty can be enabled.
Projection is required and implemented.
```

## 18. Most Accurate Wording Now

```text
Stage 10.10/10.11 establishes the core-region validity of the pairwise divA diagnostic and identifies boundary support deficiency as the dominant remaining issue. The next step is boundary-aware divA/curlA/E/J/H validation before enabling Contact A-penalty.
```

---

# Part VII. Final Route

## 19. Next Priority

```text
1. boundary divA diagnostic;
2. B=curlA exact validation;
3. E/J/Joule exact validation;
4. H=nuB and optional curlH validation;
5. SPH boundary treatment experiments;
6. projection design only;
7. Contact A-penalty later.
```

## 20. Contact Unblocking Conditions

Only consider Contact A-penalty when:

```text
1. core divA validated;
2. boundary divA explained or has treatment strategy;
3. B/E/J/Joule exact validation passed;
4. eta_A=0.1 does not cause unacceptable exact observable deviation;
5. projection not needed yet or design exists;
6. production baseline still stable.
```

---

# 21. Final Decision

Core judgment:

```text
Large boundary divA can borrow SPH boundary treatment;
currently more like boundary support / diagnostic issue, not entire A field wrong;
but must confirm physical solution accuracy via B/E/J/Joule/H exact validation;
projection: design only for now, no implementation;
eta_A is A-divergence penalty relative strength; eta_A=0.1 means penalty about 10% of A-Laplace typical strength;
next focus: boundary-aware diagnostic and exact observable validation.
```
