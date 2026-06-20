# SPHinXsys/SYCL French OPHELIE-like Particle-Integral Electromagnetic Heating Branch Development Plan

## 0. Current Decision Conclusion

The goal of this branch is not to deliver the full PDE-type A–phi solver currently in progress, but to implement a **French-literature OPHELIE-style particle-integral induction heating module** within the SPHinXsys/SYCL framework.

Recommendation: create a new branch from the **SPHinXsys mainline / clean SYCL branch**, rather than continuing development from the current A–phi branch:

```bash
cd ~/sphinxsys   # or the current clean, compilable SYCL SPHinXsys worktree
git status
git checkout -b yongchuan/french-ophelie-particle
```

This branch only implements the electromagnetic heat-source chain required for client project closure:

```text
External equivalent AC coil source
    -> Biot-Savart particle integral yields A_src / B_src inside conductive glass
        -> Internal electric scalar potential V/phi correction in glass
            -> E, J, JouleHeat
                -> JouleHeat as input to the heat equation
```

This branch temporarily does not include:

```text
Full A-phi PDE solver
Air-domain far-field boundary
Coulomb gauge / divA projection
Contact A-penalty
Real metal-surface skin-depth treatment
Cold-crucible water-cooling heat exchange
Temperature-dependent conductivity sigma(T)
Thermo-electromagnetic bidirectional iteration
Coil internal potential / impedance / skin effect
Segmented cold-crucible metal induction current
```

If the client accepts “a simplified particle/SYCL implementation based on the French-literature OPHELIE integral approach,” this route can serve as closure code; the full A–phi solver can be preserved as the follow-on JCP paper mainline.

---

## 1. French-Literature Electromagnetic Model and Correspondence with This Branch

The core electromagnetic relations in the French literature can be written in frequency domain as:

```math
\hat{\mathbf{j}}
=
-\sigma\left(\nabla \hat{V} + i\omega \hat{\mathbf{A}}\right)
```

```math
\hat{\mathbf{A}}(\mathbf{x})
=
\frac{\mu_0}{4\pi}
\int
\frac{\hat{\mathbf{j}}(\mathbf{x}')}{|\mathbf{x}-\mathbf{x}'|}
dV'
```

```math
Q_{\mathrm{Joule}}
=
\frac{|\hat{\mathbf{j}}|^2}{2\sigma}
```

Where:

- `V` is the electric scalar potential, playing the same role as `phi` in our discussion;
- `A` is the magnetic vector potential;
- `j` is the induced current density;
- `Q_Joule` is the Joule heat source entering the heat equation;
- The frequency-domain convention uniformly uses `exp(i omega t)`, therefore `dA/dt -> i omega A`.

This branch first implements a simplified version:

```math
\hat{\mathbf{A}}
\approx
\hat{\mathbf{A}}_{\mathrm{src}}
```

Where `A_src` is given only by the external coil prescribed current source passed through Biot-Savart particle integration. If closer to OPHELIE is needed later, add the glass induced-current back-reaction:

```math
\hat{\mathbf{A}}
=
\hat{\mathbf{A}}_{\mathrm{src}}
+
K[\hat{\mathbf{j}}_{\mathrm{glass}}]
```

---

## 2. Geometry and Physical Simplifications

### 2.1 Recommended Bodies

The first version needs only three bodies:

```text
GlassBody
    Conductive glass/melt, discretized as SPH particle volume, participates in EM induction.

CoilSourceBody
    TEAM7-like thick annular cylinder, as prescribed AC current source.
    Source and visualization only; no solve for coil internal potential/current distribution.

CrucibleWallBody / OuterWallBody
    Optional. Liquid wall surface or visualization only; does not participate in EM computation.
```

The first version does not build a real helical coil, thin shell, or surface EM treatment. Coil geometry can be represented as a thick annular cylinder:

```text
Rin_coil
Rout_coil
H_coil
z_center
```

Prescribe azimuthal current density on coil particles:

```math
\mathbf{J}_{\mathrm{src}} = J_0 \mathbf{e}_{\theta}
```

```math
\mathbf{e}_{\theta}
=
\left(
-\frac{y}{r},
\frac{x}{r},
0
\right),
\quad
r=\sqrt{x^2+y^2}
```

The equivalent current density can initially be taken as:

```math
J_0 =
\frac{N_{\mathrm{turn}} I_0}
{(R_{\mathrm{out}}-R_{\mathrm{in}})H_{\mathrm{coil}}}
```

Where `N_turn * I0` is the equivalent ampere-turns. For the first version, power normalization is more important; `I0` can initially be set to 1 A.

### 2.2 Recommended Default Parameters

```text
glass_radius = 0.325 m     # French-literature cold-crucible diameter 650 mm
glass_height = 0.50 m      # First version can be simplified initially
coil_inner_radius = 0.38 m
coil_outer_radius = 0.45 m
coil_height = 0.55 m
frequency = 300 kHz
sigma_glass = 16 S/m       # French natural-convection literature order of magnitude near 1473 K
target_joule_power = 50 kW or 60 kW
```

---

## 3. Staged Implementation Roadmap

## Stage 0: Geometry and VTP Output

Goal: confirm the client can see French cold-crucible-style geometry.

Implementation:

```text
GlassBody: cylinder
CoilSourceBody: thick annular cylinder
Optional WallBody: cylindrical wall/bottom
```

Acceptance:

```text
In ParaView, external coil surrounds internal glass; geometry is clear.
```

---

## Stage 1: Initialize Equivalent Coil Current Source

Register and initialize on `CoilSourceBody`:

```text
JSrcReal : Vecd
JSrcImag : Vecd
```

First version phasor current is real:

```math
\hat{I}=I_0+0i
```

Therefore:

```text
JSrcReal = J0 * e_theta
JSrcImag = 0
```

Acceptance:

```text
In VTP, JSrcReal is tangential about the z axis.
```

---

## Stage 2: Biot-Savart Long-Range Particle Integral for A_src and B_src

Register on `GlassBody`:

```text
ASrcReal : Vecd
ASrcImag : Vecd
BSrcReal : Vecd
BSrcImag : Vecd
```

First version computes only the real source field:

```text
ASrcImag = 0
BSrcImag = 0
```

For each glass particle `i`, accumulate over all coil source particles `j`:

```math
\mathbf{A}_{\mathrm{src},i}
=
\frac{\mu_0}{4\pi}
\sum_{j\in \mathrm{coil}}
\frac{\mathbf{J}_{\mathrm{src},j} V_j}
{\sqrt{|\mathbf{x}_i-\mathbf{x}_j|^2+\epsilon^2}}
```

```math
\mathbf{B}_{\mathrm{src},i}
=
\frac{\mu_0}{4\pi}
\sum_{j\in \mathrm{coil}}
\mathbf{J}_{\mathrm{src},j}V_j
\times
\frac{\mathbf{x}_i-\mathbf{x}_j}
{(|\mathbf{x}_i-\mathbf{x}_j|^2+\epsilon^2)^{3/2}}
```

The kernel here is the Biot-Savart Green function:

```math
G_A(r)=1/r,\quad G_B(r)=r/r^3
```

Not the SPH smoothing kernel `W(r,h)`.

### Why a Standard SPH Neighbor List Cannot Be Used

A standard SPH neighbor list is a compact-support local neighborhood, suitable only for:

```text
grad
div
Laplacian
heat diffusion
phi correction
```

Biot-Savart is a long-range integral; in principle every coil particle contributes to every glass particle. If a standard cutoff is used:

```math
\frac{1}{r}
\rightarrow
\frac{1}{r}\mathbf{1}_{r<r_c}
```

the magnetic field is hard-truncated at the cutoff. Central glass particles beyond the cutoff from the coil would incorrectly yield:

```text
ASrc ≈ 0
BSrc ≈ 0
E ≈ 0
J ≈ 0
Q ≈ 0
```

Therefore the first version should explicitly be written as:

```text
for each glass particle:
    loop over all coil source particles
```

That is, a `glass × coil` global source summation. It is still particle discretization, just not SPH compact-kernel discretization.

### Whether a Cell Linked List Is Needed

Not for the first version.

Later optimization can use near/far splitting or cell aggregation, but do not simply enlarge the standard neighbor-list cutoff.

---

## Stage 3: Level 0 — No Phi Correction; Compute E/J/Q Directly

Register on `GlassBody`:

```text
Sigma : Real
EReal : Vecd
EImag : Vecd
JReal : Vecd
JImag : Vecd
JouleHeat : Real
```

Using the `exp(i omega t)` convention:

```math
\hat{\mathbf{E}}
=
-i\omega \hat{\mathbf{A}}_{\mathrm{src}}
```

Split into real/imag:

```math
\mathbf{E}_{r}
=
\omega \mathbf{A}_{i}
```

```math
\mathbf{E}_{i}
=
-\omega \mathbf{A}_{r}
```

First version `A_i=0`, therefore:

```text
EReal = 0
EImag = -omega * ASrcReal
JReal = 0
JImag = sigma * EImag
JouleHeat = 0.5 * sigma * |EImag|^2
```

More general code form:

```cpp
EReal =  omega * ASrcImag;
EImag = -omega * ASrcReal;

JReal = sigma * EReal;
JImag = sigma * EImag;

JouleHeat = 0.5 * sigma * (dot(EReal, EReal) + dot(EImag, EImag));
```

Then normalize to target total power:

```math
P_{\mathrm{raw}}
=
\sum_i Q_i V_i
```

```math
Q_i
\leftarrow
Q_i
\frac{P_{\mathrm{target}}}{P_{\mathrm{raw}}}
```

Acceptance:

```text
A/B/E/J nonzero
Q >= 0
P can be normalized to 50 kW or 60 kW
When I0 doubles, unnormalized total Q power scales by roughly 4x
```

---

## Stage 4: Level 1 — Add Internal Glass Phi/V Correction

This step is closer to the French formula:

```math
\hat{\mathbf{j}}
=
-\sigma
\left(
\nabla\hat{V}
+
i\omega\hat{\mathbf{A}}_{\mathrm{src}}
\right)
```

Or using the `phi=V` notation:

```math
\hat{\mathbf{E}}
=
-\nabla\hat{\phi}
-
i\omega\hat{\mathbf{A}}_{\mathrm{src}}
```

Current continuity requires:

```math
\nabla\cdot \hat{\mathbf{j}} = 0
```

Yielding:

```math
\nabla\cdot
\left[
\sigma
\left(
\nabla\hat{\phi}
+
i\omega\hat{\mathbf{A}}_{\mathrm{src}}
\right)
\right]
=0
```

i.e.:

```math
\nabla\cdot(\sigma\nabla\hat{\phi})
=
-\nabla\cdot(i\omega\sigma\hat{\mathbf{A}}_{\mathrm{src}})
```

Split into real/imag:

```math
\nabla\cdot(\sigma\nabla\phi_r)
=
\nabla\cdot(\omega\sigma A_i)
```

```math
\nabla\cdot(\sigma\nabla\phi_i)
=
-\nabla\cdot(\omega\sigma A_r)
```

First version `A_i=0`, so solve only `PhiImag` initially:

```math
\nabla\cdot(\sigma\nabla\phi_i)
=
-\nabla\cdot(\omega\sigma A_r)
```

Then:

```cpp
EReal = -GradPhiReal + omega * ASrcImag;
EImag = -GradPhiImag - omega * ASrcReal;

JReal = sigma * EReal;
JImag = sigma * EImag;

JouleHeat = 0.5 * sigma * (dot(EReal, EReal) + dot(EImag, EImag));
```

If only `PhiImag` is enabled:

```cpp
EReal = 0.0;
EImag = -GradPhiImag - omega * ASrcReal;
JImag = sigma * EImag;
JouleHeat = 0.5 * sigma * dot(EImag, EImag);
```

### Numerical Difficulties of Phi Correction

This part can borrow from completed content on the current A–phi branch.

Difficulties include:

```text
1. Scalar Poisson/Laplace pairwise operator
2. Complex real/imag split
3. Matrix-free iterator
4. Constant nullspace handling
5. Harmonic average for variable sigma
```

For project-closure first version, constant sigma and solving only `PhiImag` significantly reduces difficulty.

### Constant Nullspace of Phi

Because only `grad(phi)` matters:

```math
\phi + C
```

is equivalent to `phi`. The constant nullspace must be handled. Optional approaches:

```text
Option A: pin one reference particle phi_ref = 0
Option B: add small penalty lambda_phi * phi
Option C: subtract volume-weighted mean after each iteration
```

First version recommendation: use small `lambda_phi` penalty or reuse the current A–phi branch approach.

---

## Stage 5: Optional OPHELIE-like Self-Induced Iteration

To get closer to French OPHELIE, add the back-reaction magnetic vector potential from glass induced current:

```math
A = A_{\mathrm{src}} + A_{\mathrm{ind}}
```

```math
A_{\mathrm{ind},i}
=
\frac{\mu_0}{4\pi}
\sum_{j\in \mathrm{glass}}
\frac{J_j V_j}{|x_i-x_j|}
```

Iteration flow:

```text
Precompute A_src from coil source.

Initialize J^0 using A_src only.

For k = 0, 1, 2, ...
    Compute A_ind^k = K[J^k] on EM particles
    A^k = A_src + A_ind^k
    Solve phi^{k+1}:
        div(sigma grad phi) = -div(i omega sigma A^k)
    Update:
        E^{k+1} = -grad(phi^{k+1}) - i omega A^k
        J^{k+1} = sigma E^{k+1}
    Check:
        ||J^{k+1} - J^k|| / ||J^k||
End
Compute JouleHeat.
```

This step is computationally heavy; recommend as an enhanced version only, not blocking Stages 0–4.

---

## 4. Acceleration Strategies

### 4.1 Most Recommended: Separate EM Particles from Fluid Particles

French literature already separates EM grid and thermal-flow grid, then interpolates. We should do the same:

```text
GlassFluidBody
    Real fluid/thermal particles; count can be large.

GlassEMBody
    Coarse EM particles for Biot-Savart / phi / J / Q only.

CoilSourceBody
    Coil source particles or source integration points.
```

Flow:

```text
CoilSourceBody + GlassEMBody
    -> solve A, phi, J, Q on coarse EM particles
        -> interpolate Q to GlassFluidBody
            -> heat/flow solver uses Q
```

Recommendation:

```text
N_fluid: 100k~500k
N_EM: 5k~20k
N_coil: 1k~5k
```

Even with self-induced `A_ind`, this avoids direct full-fluid-particle `N^2`.

### 4.2 Source Splitting

Decompose:

```math
A = A_{src} + A_{ind}
```

Where:

```text
A_src = K[J_coil] computed once.
A_ind = K[J_glass] updated each iteration.
```

Do not recompute coil -> glass every OPHELIE iteration.

### 4.3 Near/Far Splitting

For glass -> glass long-range integration:

```text
near cells:
    direct particle summation

far cells:
    cell aggregated as equivalent source points
```

First-order far-field approximation:

```math
M_c = \sum_{j\in c} J_j V_j
```

```math
A_i^{far}
\approx
\sum_{c\in far}
\frac{\mu_0}{4\pi}
\frac{M_c}{|x_i-x_c|}
```

Simpler than FMM; suitable for later optimization.

### 4.4 Axisymmetric EM

If the client case temporarily does not require stirrer effects on EM, the closest French acceleration is:

```text
2D/axisymmetric EM source
    -> Q(r,z)
        -> rotate/interpolate to 3D SPH
```

In the stirring paper, French authors also use 2D axisymmetric EM + 3D thermohydrodynamics mixed coupling.

---

## 5. Recommended Code File Structure

Recommend creating:

```text
src/shared/electromagnetic_dynamics/
    electromagnetic_ophelie_particle_parameters.h
    electromagnetic_ophelie_particle_variables.h
    electromagnetic_ophelie_particle_source.h
    electromagnetic_ophelie_particle_source.hpp
    electromagnetic_ophelie_particle_biot_savart.h
    electromagnetic_ophelie_particle_biot_savart.hpp
    electromagnetic_ophelie_particle_phi.h
    electromagnetic_ophelie_particle_phi.hpp
    electromagnetic_ophelie_particle_postprocess.h
    electromagnetic_ophelie_particle_postprocess.hpp

tests/extra_source_and_tests/test_3d_french_ophelie_particle/
    test_3d_french_ophelie_particle.cpp
    CMakeLists.txt
```

Can merge into fewer files first and split after the pipeline runs.

---

## 6. Variables That Must Be Registered

### CoilSourceBody

```text
JSrcReal : Vecd
JSrcImag : Vecd
```

### GlassEMBody / GlassBody

```text
Sigma : Real

ASrcReal : Vecd
ASrcImag : Vecd

BSrcReal : Vecd
BSrcImag : Vecd

PhiReal : Real
PhiImag : Real

GradPhiReal : Vecd
GradPhiImag : Vecd

EReal : Vecd
EImag : Vecd

JReal : Vecd
JImag : Vecd

JouleHeat : Real

DivJReal : Real    # diagnostic
DivJImag : Real    # diagnostic
```

First-round minimum output:

```text
JSrcReal
ASrcReal
BSrcReal
EImag
JImag
JouleHeat
Sigma
```

---

## 7. Content to Borrow/Copy from the Current A–phi Branch

Do not copy the full A–phi solver to the client branch. Borrow only the following minimal components.

### 7.1 Scalar Pairwise Laplace / Poisson Operator

Purpose:

```text
Solve phi correction:
div(sigma grad phi) = RHS
```

Find files/classes on the current A–phi branch related to:

```text
electromagnetic_aphi_matrix_free_operators.h/.hpp
electromagnetic_aphi_matrix_free_residuals.h/.hpp
electromagnetic_aphi_matrix_free_solver.h/.hpp
electromagnetic_aphi_matrix_free_pairwise_graph.h/.hpp
```

Extract:

```text
scalar variable-coefficient Laplace apply
pairwise graph / neighbor loop
harmonic sigma average
matrix-free residual computation
```

Do not extract:

```text
vector A operator
divA regularization
curl-curl / B curl diagnostics
A contact penalty
full coupled A-phi block residual
```

### 7.2 Complex Real/Imag Data Organization

Purpose:

```text
PhiReal/PhiImag
EReal/EImag
JReal/JImag
```

Borrow:

```text
real/imag split approach from complex scalar Helmholtz / scalar Poisson tests
```

Requirement:

```text
Do not use std::complex inside device kernels; continue storing Real arrays / Vecd arrays separately.
```

### 7.3 GMRES / Jacobi / Block-Jacobi Iterators

Purpose:

```text
Scalar linear system for phi correction
```

Can borrow:

```text
restarted GMRES
Jacobi or block-Jacobi preconditioner
residual norm reporting
max iteration / tolerance handling
```

First version can also use simplified Jacobi / Richardson / CG-like variants as long as the pipeline runs.

### 7.4 Phi Constant Nullspace Handling

Purpose:

```text
phi + C nullspace
```

Borrow from the current A–phi branch:

```text
lambda_phi gauge penalty
pin reference particle
mean removal
```

First version recommendation:

```text
small lambda_phi penalty
```

or:

```text
pin one particle phi = 0
```

### 7.5 VTP Output and Test Framework

Borrow:

```text
existing TEAM7-like scaffold in tests/extra_source_and_tests
field registration and output style from the current A–phi branch
```

Output fields must be easy to inspect in ParaView:

```text
ASrcReal
BSrcReal
EImag
JImag
JouleHeat
```

### 7.6 Diagnostic Tools

Borrow:

```text
L2 norm
relative residual
sum(Q*Vol)
scaling check
divJ diagnostic
```

At minimum implement:

```text
P_raw = sum(Q_i V_i)
P_scaled = sum(Q_scaled_i V_i)
min/max/avg JouleHeat
||J||_2
||divJ||_2 / ||J||_2
```

---

## 8. Content Not to Copy from the A–phi Branch

These belong to the paper mainline; not recommended for the client branch:

```text
full vector A unknown solver
Coulomb gauge / divA projection
eta_A regularization
lambda_A contact A-penalty
B=curlA corrected curl dual-track
air padding sensitivity
multi-body contact A-phi coupled GMRES
TEAM7 absolute 1e-5 validation scaffold
```

The client branch should stay lightweight, engineering-oriented, and easy to explain:

```text
OPHELIE-like particle integral induction module
```

---

## 9. Update Recommendations for Previous Starter Code

The previous `french_coil_source_starter.zip` needs these updates:

```text
1. Rename from FrenchCoilSource to OphelieParticle / OphelieLike.
2. Clearly create branch from SPHinXsys mainline, not dependent on A–phi branch.
3. Biot-Savart section must state global source summation, not SPH neighbor relation.
4. Level 1 phi correction: solve PhiImag only first, because ASrc is real-only in first version.
5. Add optional GlassEMBody / GlassFluidBody separation architecture description.
6. Add optional Stage 5 OPHELIE-like self-induced iteration.
7. Delete or weaken full A–phi solver wording to avoid Cursor misusing vector A/gauge/divA.
```

---

## 10. Cursor First-Round Task Manifest

Use the content below directly as first-round development tasks.

### Task 1: Create Test Case

```text
tests/extra_source_and_tests/test_3d_french_ophelie_particle/
```

Create:

```text
GlassBody: cylinder
CoilSourceBody: thick annular cylinder
```

No fluid time integration yet; run one EM source-field computation and output VTP.

### Task 2: Register Variables

Register:

```text
JSrcReal / JSrcImag on CoilSourceBody

Sigma
ASrcReal / ASrcImag
BSrcReal / BSrcImag
EReal / EImag
JReal / JImag
JouleHeat
on GlassBody
```

### Task 3: Initialize Coil Current Source

Implement CK kernel:

```cpp
JSrcReal[i] = J0 * e_theta(pos[i]);
JSrcImag[i] = Vecd::Zero();
```

### Task 4: Implement Coil -> Glass Biot-Savart Global Kernel

Each glass particle loops over all coil particles:

```cpp
for glass_i:
    A = 0
    B = 0
    for coil_j:
        r = x_glass[i] - x_coil[j]
        A += coeff * Jsrc[j] * Vol[j] / |r|
        B += coeff * cross(Jsrc[j] * Vol[j], r) / |r|^3
```

Do not use compact-support neighbor list.

### Task 5: Compute Level 0 E/J/Q

```cpp
EReal = omega * ASrcImag;
EImag = -omega * ASrcReal;

JReal = sigma * EReal;
JImag = sigma * EImag;

JouleHeat = 0.5 * sigma * (|EReal|^2 + |EImag|^2);
```

### Task 6: Power Normalization

```cpp
P_raw = sum(JouleHeat[i] * Vol[i]);
JouleHeat[i] *= P_target / P_raw;
```

### Task 7: VTP Output and Log

Output:

```text
ASrcReal
BSrcReal
EImag
JImag
JouleHeat
Sigma
```

Log output:

```text
N_glass
N_coil
frequency
sigma
P_raw
P_scaled
min/max JouleHeat
```

---

## 11. Second-Round Task: Phi Correction

After first-round figures, have Cursor borrow the scalar solver from the A–phi branch:

```math
\nabla\cdot(\sigma\nabla\phi_i)
=
-\nabla\cdot(\omega\sigma A_r)
```

Solve only `PhiImag`:

```text
PhiImag
GradPhiImag
```

Then update:

```cpp
EImag = -GradPhiImag - omega * ASrcReal;
JImag = sigma * EImag;
Q = 0.5 * sigma * |EImag|^2;
```

Added diagnostics:

```text
DivJImag
||DivJImag|| / ||JImag||
```

---

## 12. Final Delivery Narrative Recommendation

Describe to the client as:

```text
This module is based on the OPHELIE integral-method idea from French cold-crucible literature,
implemented in the SPHinXsys/SYCL particle framework for external AC coil source,
Biot-Savart magnetic vector potential / magnetic field integration,
internal electric scalar potential correction in conductive melt,
induced current and Joule heat source computation.

This stage uses an equivalent thick annular coil and constant-conductivity conductive melt,
without metal skin-depth surface model, water-cooled wall, temperature-dependent conductivity,
or thermo-electromagnetic feedback.
Joule heat source can serve as input for subsequent heat conduction / flow computation.
```

This aligns with French literature without delivering a full A–phi PDE solver.
