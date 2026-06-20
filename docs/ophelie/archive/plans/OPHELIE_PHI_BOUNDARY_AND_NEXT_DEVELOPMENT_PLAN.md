# OPHELIE-like French Reduced: Phi Boundary Treatment, Accuracy Bottlenecks, and Next-Stage Development Plan

> Purpose: for Cursor follow-on development. This document explains why the current French reduced OPHELIE-like EM chain still needs glass-domain boundary treatment even though the air domain is not computed; why an outer dummy shell is not recommended as priority; and how to reference SPHinXsys fluid/thermal boundary ideas to progressively improve the `phi` equation, `divJ`, and JouleHeat accuracy.

---

## 0. Current positioning

The current route is neither full OPHELIE nor full Maxwell PDE, but a French reduced OPHELIE-like particle-integral EM chain on SPHinXsys/SYCL:

```text
multiloop filament coil (I dl, Biot-Savart)
  -> A_src, B_src on glass particles
  -> solve phi correction in glass domain
  -> E, J, JouleHeat
  -> target total Joule power, e.g. 50 kW
```

Current outputs on cylindrical glass geometry:

```text
A_src / B_src
phi / grad_phi
E / J / divJ
JouleHeat
P_total = sum(Q_i V_i)
```

The largest numerical issue:

```text
phi_eq_res_vol ≈ 0.45--0.55 saturated
GMRES outer iteration increase barely lowers it further
reconstructed divJ improves somewhat but remains too large
```

This looks more like discrete compatibility, boundary flux, and operator consistency issues than insufficient GMRES iterations alone.

---

## 1. Why boundary issues remain when the air domain is not computed

In OPHELIE-like reduced EM, the air domain need not be a PDE computation domain because:

1. Coil source fields `A_src/B_src` come directly from Biot-Savart integration;
2. Glass-induced current fields `A_ind/B_ind` can also come from Biot-Savart integration over glass particle current;
3. Therefore magnetic vector potential and field need not be solved on an air volume grid.

But the `phi` equation still solves an elliptic correction problem inside the glass domain. Continuous form:

```math
\nabla \cdot \left[\sigma\left(\nabla \phi + \omega A\right)\right] = 0
```

Equivalent to:

```math
\nabla \cdot (\sigma \nabla \phi)
= -\omega \nabla \cdot (\sigma A)
```

Where the current reduced route usually first takes:

```math
A \approx A_{src}
```

Later self-induction route takes:

```math
A = A_{src} + A_{ind}, \quad A_{ind}=K[J_{glass}]
```

This `phi` equation is defined only on glass domain `Omega_glass`, so boundary conditions are required. Not computing air does not mean glass boundary has no conditions.

Physical boundary condition: no current leaves the outer glass boundary:

```math
n \cdot J = 0
```

Since:

```math
J = \sigma E,
\quad
E_{imag} = -\nabla \phi_{imag} - \omega A_{real}
```

the boundary condition is inhomogeneous Neumann:

```math
n \cdot \sigma(\nabla \phi + \omega A) = 0
```

i.e.:

```math
n \cdot \sigma \nabla \phi = -\omega \; n \cdot \sigma A
```

If `sigma` is scalar and nonzero, locally:

```math
\frac{\partial \phi}{\partial n} = -\omega\; n\cdot A
```

So the boundary issue is not "air domain must not be computed", but:

```text
Glass-domain phi equation needs no-normal-current boundary conditions.
SPH kernel truncation at boundary: without boundary flux or ghost state, DivSigmaGrad and DivSigmaA are hard to discretize compatibly.
```

---

## 2. What are dummy shell / ghost particles? Are they suitable for complex geometry?

### 2.1 What dummy shell is

Dummy shell or ghost shell usually means placing one or more layers of non-physical particles outside the real geometry boundary to complete SPH kernel support. They do not represent real air or real conductor; they help inner glass particles' gradient, divergence, or Laplace operators near the boundary approach full support.

For example on a cylinder, virtual points can be placed outward along the normal with Neumann-satisfying `phi_ghost`.

### 2.2 Why dummy shell is not recommended as first choice

For cylinders or regular geometry, dummy shell is easy to generate; for real cold crucible, complex STL, segmented geometry, or geometry with coil/stirrer, generating high-quality outer layers is not robust. Main issues:

```text
1. Outer normals discontinuous on complex geometry; shell particles overlap or miss coverage;
2. Concave corners, narrow gaps, thin walls: shell particle quality hard to guarantee;
3. Every CAD change requires shell generation and relax;
4. Increases contact/relation and output complexity;
5. Dummy particles easily mistaken for physical particles, polluting JouleHeat or P_total.
```

Therefore dummy shell can serve as diagnostic or fallback, but should not be the current main route.

---

## 3. Can we reference SPHinXsys fluid one-sided Riemann boundary ideas?

Reference the **idea**, not the form directly.

SPHinXsys fluid one-sided Riemann boundary essence:

```text
Need not generate full fluid domain outside wall;
Pass wall normal, wall state, and ghost/extrapolated state to construct boundary flux;
Impose no-penetration, slip/no-slip, etc. through discrete flux.
```

The `phi` equation is not hyperbolic fluid; no Riemann solver needed; but similar boundary-state construction applies:

```text
Do not generate physical air particles outside glass;
For boundary-near glass particles, use level set normal or geometry normal for one-sided Neumann boundary state;
Add n·sigma(grad_phi + omega A)=0 as boundary flux to phi operator or RHS.
```

This route suits complex geometry better than dummy shell because it uses existing geometry/level-set boundary normals rather than extra volumetric shell particles.

---

## 4. Can we reference SPHinXsys thermal diffusion boundary form?

Yes, and this is the more recommended main route.

Mathematically `phi` is close to diffusion/Poisson:

```math
\nabla\cdot(\sigma\nabla\phi)=b
```

Thermal conduction:

```math
\nabla\cdot(k\nabla T)=source
```

So `phi` LHS can reference SPHinXsys thermal diffusion operators; boundary conditions analogize to thermal Dirichlet or Neumann.

Physical notes:

### 4.1 `phi` should not default to Dirichlet boundary

Unless glass boundary has electrodes or specified potential, do not simply set:

```math
\phi = 0 \quad \text{on boundary}
```

That forces electric field incorrectly; usually wrong for cold-crucible glass melt boundary.

Dirichlet is better for gauge lock, e.g.:

```text
mean(phi)=0
or fix one reference particle phi=0
```

Used to remove constant null space, not as physical boundary condition.

### 4.2 `phi` physical boundary is closer to Neumann

In cold-crucible reduced glass-only EM, physical boundary is no normal current:

```math
n \cdot J = 0
```

Corresponding inhomogeneous Neumann:

```math
n \cdot \sigma \nabla\phi = -\omega n\cdot\sigma A
```

Reference thermal specified heat-flux boundary:

```math
n\cdot k\nabla T = q_n
```

Here:

```math
n\cdot\sigma\nabla\phi = q_{phi,n}
= -\omega n\cdot\sigma A
```

---

## 5. Recommended implementation route: one-sided Neumann + compatibility diagnostics before complex dummy shell

### Overall recommendation

Do not prioritize generic complex-geometry dummy shell now. Proceed in this order:

```text
P0: strengthen phi boundary diagnostics and RHS compatibility diagnostics
P1: one-sided Neumann boundary correction without outer particles
P2: corrected gradient / linear consistency option
P3: virtual shell as cylinder/simple-geometry diagnostic only
P4: then A_ind self-induction and thermal one-way/weak coupling
```

---

## 6. P0: diagnostics to add first

### 6.1 RHS compatibility / zero-mean projection

Neumann-type Poisson problems need RHS compatibility. Discretely check at least:

```math
\sum_i b_i V_i \approx 0
```

Where:

```math
b_i = -\omega\;\nabla_h\cdot(\sigma A)_i
```

Recommended diagnostics:

```text
phi_rhs_volume_integral = sum(rhs_i * V_i)
phi_rhs_volume_mean     = sum(rhs_i * V_i) / sum(V_i)
phi_rhs_l2_before_proj
phi_rhs_l2_after_proj
phi_eq_res_before_rhs_projection
phi_eq_res_after_rhs_projection
```

Simplest projection:

```math
b_i^{proj} = b_i - \frac{\sum_j b_jV_j}{\sum_j V_j}
```

Development requirement:

```text
Add CLI:
  --phi-rhs-project-zero-mean=0|1

Default recommendation:
  diagnostic output first without changing results by default;
  enable in literature-mode after comparison confirms benefit.
```

Acceptance:

```text
If eq_res drops clearly, compatibility is a main factor;
If almost unchanged, main issue is boundary flux or operator range mismatch.
```

---

### 6.2 Boundary particle identification and normal diagnostics

No outer particles needed first, but must know which glass particles are boundary-near and what the outward normal is.

Recommended fields:

```text
PhiBoundaryFlag       scalar, 0/1
PhiBoundaryNormal     vector
PhiBoundaryDistance   scalar
PhiKernelDeficiency   scalar
JnBoundary            scalar
```

Boundary identification — pick one of two methods initially:

#### Method A: level set / shape normal

If `TriangleMeshShapeCylinder + LevelSet` provides signed distance and normal, prioritize it. Cursor should check SPHinXsys LevelSetShape / level-set data interfaces.

Logic:

```text
if abs(signed_distance_i) < c * dp or kernel_deficiency_i > threshold:
    boundary_flag_i = 1
    normal_i = outward_normal_from_level_set
```

#### Method B: kernel deficiency

If level-set normal is inconvenient, use insufficient neighbor support:

```text
kernel_support_sum_i = sum_j W_ij V_j
kernel_deficiency_i = reference_support - kernel_support_sum_i
```

Large kernel deficiency likely means boundary-near particle; but normal is less reliable than level set — use as auxiliary only.

### 6.3 Output boundary normal current

Need not only bulk `divJ` but boundary normal current leakage:

```math
J_n = n\cdot J = n\cdot\sigma(-\nabla\phi-\omega A)
```

Output:

```text
JnBoundary_min/max/mean/l2
JnBoundary_rel = ||Jn||_boundary / ||J||_boundary
```

Acceptance target:

```text
After phi correction, JnBoundary_rel should drop clearly vs no-phi;
If divJ improves but Jn does not, boundary condition did not enter phi solve.
```

---

## 7. P1: one-sided Neumann boundary correction without outer geometry particles

### 7.1 Physical boundary condition

For boundary particles, impose:

```math
n\cdot\sigma\nabla\phi = -\omega n\cdot\sigma A
```

Let:

```math
g_n = \frac{\partial\phi}{\partial n} = -\omega n\cdot A
```

If `sigma` is scalar, it cancels on both sides; if later `sigma(T)` varies spatially, keep flux form:

```math
q_{phi,n}= -\omega n\cdot\sigma A
```

### 7.2 Ghost state form without explicit particles

For boundary particle `i`, construct virtual ghost point at distance `d_b` along outward normal locally in the operator, not as SPH particle, not in neighbor relation.

If ghost at:

```math
x_g = x_i + d_b n_i
```

Neumann:

```math
\frac{\phi_g - \phi_i}{d_b} = g_n
```

then:

```math
\phi_g = \phi_i + d_b g_n
```

Similar to fluid one-sided ghost state, but scalar Neumann for `phi`.

### 7.3 How to enter LHS/RHS

Do not output ghost as real particle or include in JouleHeat integral. In `DivSigmaGrad` operator, add local boundary term for boundary particles:

```text
L_phi_i += boundary_coeff_i * (phi_g - phi_i)
```

Since:

```math
phi_g - phi_i = d_b g_n
```

this mainly enters RHS; full Robin/Dirichlet ghost would also change LHS.

First version as RHS correction:

```text
rhs_i += boundary_flux_correction_i
```

Points:

```text
1. Implement as diagnostic first; do not default affect literature-mode;
2. Apply only to boundary_flag_i=1 particles;
3. Sign must be validated by MMS / 1D slab test;
4. Correction must not change P_total definition; P_total still integrates glass particles only.
```

### 7.4 CLI design

Recommended:

```text
--phi-boundary-mode=none|one-sided-neumann|virtual-shell-diagnostic
--phi-boundary-normal-source=level-set|kernel-deficiency|analytic-cylinder
--phi-boundary-distance-factor=0.5
--phi-boundary-diagnostics=1
```

First-stage default:

```text
--phi-boundary-mode=none
```

Explicitly enable in tests:

```text
--phi-boundary-mode=one-sided-neumann
```

Consider literature-mode default only after validation stabilizes.

---

## 8. P2: corrected gradient / linear consistency

Current `grad_phi`, `div_sigma_A`, `DivSigmaGrad` boundary and non-uniform particle errors may be large. Add corrected gradient option:

```math
\nabla \phi_i
= B_i^{-1}\sum_j (\phi_j-\phi_i)\nabla W_{ij} V_j
```

Where:

```math
B_i = \sum_j (x_j-x_i)\otimes \nabla W_{ij} V_j
```

Development recommendation:

```text
Add CLI:
  --phi-gradient-correction=0|1

Add tests:
  test_3d_ophelie_phi_gradient_linear
  test_3d_ophelie_phi_grad_div_laplace_consistency
  test_3d_ophelie_phi_solve_manufactured_A
```

Acceptance:

```text
linear phi gradient error drops significantly;
MMS grad/div/laplace consistency improves;
French real Biot case phi_eq_res_vol and divJ_L2_red do not worsen, preferably improve.
```

Note: corrected gradient may break strict symmetry; GMRES is acceptable but record conditioning and runtime.

---

## 9. P3: virtual shell diagnostic only, not main route

On cylinder case, keep virtual shell diagnostic to answer:

```text
If kernel support at boundary is completed, does phi_eq_res drop clearly?
```

Not the main route for complex geometry.

Recommended mode:

```text
--phi-boundary-mode=virtual-shell-diagnostic
```

Limits:

```text
1. analytic cylinder or simple box/slab only;
2. shell particles do not output JouleHeat, not in P_total;
3. compare one-sided Neumann vs support completion only;
4. do not require arbitrary STL support.
```

---

## 10. P4: subsequent A_ind and thermal coupling

After boundary and `phi` stabilize, proceed:

```math
A = A_{src} + A_{ind}, \quad A_{ind}=K[J_{glass}]
```

Current A_ind/A_coil magnitude is already non-negligible; treat as separate self-induction acceptance.

Recommended:

```text
--self-induction=0|1
--self-induction-max-iter=N
--self-induction-relax=alpha
```

Acceptance output:

```text
A_ind/A_src
B_ind/B_src
J_picard_rel_change
P_total_picard_rel_change
divJ_L2_red
phi_eq_res_vol
```

Thermal coupling one-way first:

```math
\rho C_p \frac{dT}{dt}=Q_{Joule}
```

No natural convection or complex cooling initially; energy closure only:

```math
\Delta E_{thermal} \approx \int P_{Joule}(t) dt
```

Then add cooling-wall Neumann/Robin thermal boundary step by step.

---

## 11. Tests to add or update

### 11.1 Basic boundary condition test

Add:

```text
test_3d_ophelie_phi_neumann_slab
```

Goal: regular slab with analytic A and phi, known boundary:

```math
\partial_n \phi = -\omega n\cdot A
```

Check:

```text
grad_phi_error
boundary_Jn_error
phi_eq_res_vol
```

### 11.2 Cylinder analytic normal test

Add:

```text
test_3d_ophelie_phi_neumann_cylinder
```

Analytic cylinder normals:

```text
side: n = (x, y, 0) / sqrt(x^2+y^2)
top:  n = (0,0,1)
bottom: n = (0,0,-1)
```

Validates one-sided Neumann sign and boundary Jn.

### 11.3 French reduced regression

Update:

```text
test_3d_ophelie_french_reduced
```

Add output:

```text
phi_rhs_volume_mean
phi_rhs_projection_ratio
JnBoundary_L2
JnBoundary_rel
phi_boundary_mode
phi_gradient_correction
```

Compare four groups:

```text
A. baseline: no rhs projection, no boundary correction
B. rhs projection only
C. one-sided Neumann only
D. rhs projection + one-sided Neumann + optional corrected gradient
```

Acceptance focus:

```text
P_total still near target power;
divJ_L2_red improves;
JnBoundary_rel improves;
phi_eq_res_vol does not worsen, preferably improves;
JouleHeat spatial distribution has no unphysical spikes.
```

---

## 12. Current Dirichlet / Neumann conclusions

### Dirichlet acceptable when

```text
1. gauge lock: fix mean(phi)=0 or one reference particle phi=0;
2. manufactured solution test: known phi boundary for operator validation;
3. real physics with electrode/specified potential boundary.
```

### Dirichlet not recommended when

```text
French reduced cold crucible glass-only EM main chain; do not default phi=0 on glass boundary.
```

Reason: pins potential boundary incorrectly; not equivalent to no current through boundary.

### Neumann priority when

```text
No current leaves outer glass boundary, i.e. n·J=0.
```

i.e.:

```math
n\cdot\sigma\nabla\phi = -\omega n\cdot\sigma A
```

This is the more natural boundary for current OPHELIE-like reduced glass-domain `phi` equation.

---

## 13. Cursor execution manifest

### Step 1: diagnostics first

1. Add RHS compatibility output;
2. Add `--phi-rhs-project-zero-mean=0|1`;
3. Add boundary flag / normal / JnBoundary output;
4. Update French reduced CSV/console summary.

### Step 2: boundary no-normal-current

1. Add `--phi-boundary-mode=none|one-sided-neumann`;
2. Prefer level-set normal; if interface unclear, French cylinder uses analytic normal first;
3. Add Neumann boundary flux correction in phi operator/RHS;
4. Validate sign with slab/cylinder MMS.

### Step 3: consistency enhancement

1. Add optional corrected gradient;
2. Run gradient linear test;
3. Run manufactured A test;
4. Run French reduced baseline comparison.

### Step 4: diagnostic virtual shell only

1. Cylinder support only;
2. Exclude from P_total;
3. Not main route for complex geometry;
4. Judge theoretical eq_res improvement from boundary support completion.

### Step 5: A_ind and thermal

1. Self-induction mode independent gate;
2. One-way JouleHeat to heat;
3. Energy conservation test;
4. Later sigma(T) and thermal-electromagnetic Picard.

---

## 14. Stage assessment

Do not simply explain `phi_eq_res≈0.5` as insufficient GMRES. More reasonable:

```text
current SPH DivSigmaGrad / DivSigmaA on real Biot A and bounded glass domain
still lacks sufficiently strict Neumann boundary flux closure and discrete compatibility control.
```

Do not invest immediately in complex outer dummy shell — not robust for real complex geometry. Prefer SPHinXsys boundary ideas:

```text
Use level-set/geometry normal for one-sided Neumann boundary state,
reference scalar diffusion/heat-transfer boundary flux implementation,
close n·J=0 inside glass domain directly,
add RHS projection and boundary Jn diagnostics.
```

This matches OPHELIE-like glass-only physics and future complex STL geometry on SPHinXsys/SYCL.
