# OPHELIE Edge-Flux Solver Development Plan and Legacy Code Organization Guide

> Cursor development document  
> Project: SPHinXsysSYCL / French reduced OPHELIE-like cold-crucible EM solver  
> Current decision: shift from `div-grad baseline` residual-floor tuning to **SPH pairwise edge-flux current formulation**.  
> Important constraint: “edge” in edge-flux means **SPH neighbor particle pair i-j**, not external finite-volume grid, FEM edge, or host CSR graph. All production paths must stay SPHinXsys / SYCL CK style.

---

## 0. Purpose of This Document

Guide next-stage OPHELIE-like solver development. We can already run the French reduced demo:

```text
multiloop Biot A_src/B_src
→ phi correction
→ E/J/JouleHeat
→ 50 kW literature calibration
```

But the `div-grad baseline` route long has:

```text
phi_eq_res_vol ≈ 0.48–0.49
```

Even after CLI fixes, RHS lifecycle fixes, Neumann boundary, B-matrix corrected gradient, compatible D/G, vector divergence MMS, edge-flux diagnostics, this residual floor is not truly resolved.

Therefore this stage is not tuning the demo further but developing an **edge-flux current formulation** closer to French OPHELIE:

```text
phi edge equation
→ pairwise edge current q_ij
→ edge-consistent JouleHeat
→ particle J reconstruction
→ later A_ind = K[J]
→ later Picard self-induction OPHELIE-like solver
```

---

## 1. Why Shift from div-grad Baseline to Edge-Flux

### 1.1 Current div-grad Baseline Status

Current production route:

```text
--phi-projection-operator=div-grad
```

Rough results:

```text
phi_eq_res_vol ≈ 0.486
divJ_L2_red   ≈ 2.06
P_raw         ≈ 50000 W
literature    = PASS
```

Useful for current particle-level postprocess `E/J/divJ/JouleHeat`; keep as demo / fallback / reference.

Residual stuck near 0.49. We already ruled out or weakened:

```text
1. CLI parameters not taking effect;
2. RHS overwritten after finalize by GMRES entry;
3. insufficient GMRES outer iterations alone;
4. RHS volume mean compatibility;
5. one-sided Neumann boundary as dominant source;
6. simple corrected gradient fix;
7. basic D_unc divergence error on linear/rotational fields.
```

Likely cause: elimination form

\[
L_h \phi = D_h(\sigma G_h\phi),
\]

\[
b_h = -\omega D_h(\sigma A_{src})
\]

has discrete range / compatibility floor on real Biot vector potential \(A_{src}\). Real \(A_{src}\) is solenoidal-curl, not from same \(G_h\phi\). Small patches on this route unlikely to push `phi_eq_res` very low.

### 1.2 French OPHELIE Core Is Not “Post-Process J After Phi”

French OPHELIE core:

\[
\nabla\cdot \mathbf j = 0,
\]

\[
\mathbf j = -\sigma\left(\nabla V + \frac{\partial \mathbf A}{\partial t}\right),
\]

\[
\mathbf A(M)=\frac{\mu_0}{4\pi}\int_V \frac{\mathbf j}{r}\,dV.
\]

Current density \(\mathbf j\) is a core unknown, not last-step particle gradient postprocess from an incompatible phi solve.

Current `div-grad baseline` is equivalent in continuous math to eliminating \(\mathbf j\), but not necessarily in SPH discretization:

```text
phi equation uses one LHS/RHS set;
E/J postprocess uses another gradient;
divJ diagnostic may use a third divergence.
```

Source of current contradictions.

### 1.3 Edge-Flux Is SPH Version of “Current as Core Variable”

Edge-flux: do not treat current as postprocess only; define discrete current flux \(q_{ij}\) on each SPH neighbor pair with pairwise current conservation as the equation.

Edge means:

```text
SPH neighbor pair i-j
```

Not:

```text
FVM mesh edge
FEM edge
CSR graph edge
host adjacency list
external grid
```

Still SPH: particle i loops inner-relation neighbors j, using SPH kernel / kernel gradient / particle volume / smoothing length / material coefficients for pairwise flux.

---

## 2. Hard Implementation Constraints

### 2.1 Forbidden Implementation Approaches

Cursor must not use:

```text
1. external finite-volume grid;
2. FEM/mesh edge;
3. CSR sparse matrix graph;
4. host-side std::vector<std::vector<int>> adjacency;
5. pure CPU pair loop as production;
6. bypass SPHinXsys InnerRelation;
7. custom graph parallel to SPH neighbor relation;
8. bypass SYCL CK backend for edge-flux.
```

### 2.2 Required Implementation Approach

Production code must:

```text
1. use SPH particles;
2. use SPHinXsys cell linked list / inner relation;
3. use LocalDynamics / Interaction<Inner<>>;
4. use InteractionDynamicsCK or StateDynamics CK style;
5. use UpdateKernel / InteractKernel;
6. use DelegatedData / DiscreteVariable;
7. run on SYCL GPU backend;
8. compute q_ij on-the-fly in neighbor loop;
9. pair weight C_ij from SPH kernel or existing pairwise Laplace weight.
```

---

## 3. Positioning Old vs New Routes

### 3.1 Preserve but No Longer Primary

Keep:

```text
--phi-projection-operator=div-grad
```

For:

```text
1. demo fallback;
2. regression baseline;
3. comparison with edge-flux;
4. current French reduced visualization reference.
```

No longer invest heavily in:

```text
phi_eq_res 0.49 → 0.1
```

Unless new theoretical evidence appears.

### 3.2 Freeze as Diagnostic

Preserve source and tests but not production:

```text
1. compatible-div-grad / G_c-D_c;
2. phi-gradient-correction alone;
3. grad Neumann postprocess;
4. one-sided Neumann as default route;
5. legacy edge-flux phi-only production;
6. singular_toroidal_legacy field as method conclusion basis.
```

Do not delete; archive for papers or debug explaining route change.

---

## 4. Legacy Code, Test, and Document Archive Requirements

Goal: do not delete old work but avoid interfering with new mainline.

### 4.1 Recommended Archive Directory

At repo root or under `docs/ophelie/`:

```text
docs/ophelie/archive_phi_divgrad_residual_floor/
```

Recommend placing or copying:

```text
OPHELIE_PHI_DEVELOPMENT_LOG.md
OPHELIE_PHI_DISCRETIZATION_NEXT_SPRINT_PLAN.md
OPHELIE_DIVERGENCE_OPERATOR_VALIDATION_AND_EDGE_FLUX_PLAN.md
OPHELIE_GUIDE_MISUNDERSTANDING_AUDIT_AND_CURSOR_IMPLEMENTATION_PLAN.md
OPHELIE_PHI_BOUNDARY_AND_NEXT_DEVELOPMENT_PLAN.md
```

If not moving originals, add README in archive:

```text
README_PHI_DIVGRAD_ARCHIVE.md
```

Describing diagnostic history of old route.

### 4.2 Recommended Archive README Content

README must state:

```text
1. div-grad baseline kept as fallback/reference;
2. compatible-div-grad, boundary Neumann, grad projection are diagnostic;
3. current residual floor ~ 0.49;
4. no longer primary goal to lower div-grad phi_eq_res;
5. new mainline is SPH edge-flux current formulation.
```

### 4.3 Do Not Mass-Move Code Initially

Do not move many `.h/.hpp/.cpp` at once—breaks CMake. Steps:

```text
Step 1: document archive;
Step 2: CLI mark diagnostic-only;
Step 3: new edge-flux module under existing OPHELIE shared directory;
Step 4: preserve old tests, not main CI;
Step 5: reorganize CMake targets after edge-flux stable.
```

### 4.4 Old Test Naming and Comments

Preserve old tests; mark in README or test output:

```text
[diagnostic-only]
```

Including:

```text
test_3d_ophelie_phi_neumann_slab
test_3d_ophelie_phi_neumann_cylinder
test_3d_ophelie_phi_compatible_operator_mms
test_3d_ophelie_vector_divergence_mms
```

Vector divergence MMS still valuable for operator diagnostic, not edge-flux production acceptance itself.

---

## 5. Edge-Flux Solver Mathematical Form

### 5.1 Basic Variables

Particle variables:

```text
phi_i       : scalar electric potential correction, imaginary component
A_i         : real vector potential, initially A_src_real
sigma_i     : electrical conductivity
V_i         : particle volume
x_i         : particle position
```

Neighbor-pair variables:

```text
r_ij        = x_i - x_j
A_ij_avg    = 0.5 * (A_i + A_j)
sigma_ij    = harmonic average of sigma_i and sigma_j
C_ij        = SPH pair conductance, must reuse existing pairwise Laplace/legacy weight
edge_drop   = potential + inductive electric drop along pair
q_ij        = pairwise current flux
```

### 5.2 Edge Electric Drop

Match existing legacy-pairwise sign; confirm with legacy operator. Recommend:

\[
e_{ij} = (\phi_i-\phi_j) + \omega \bar{\mathbf A}_{ij}\cdot(\mathbf x_i-\mathbf x_j).
\]

Or overall sign per code convention, but must ensure:

```text
LHS pairwise term and RHS A-flux term use exactly the same C_ij and sign convention.
```

### 5.3 Pair Current Flux

\[
q_{ij} = - C_{ij} e_{ij}.
\]

Must satisfy:

\[
q_{ij}=-q_{ji}.
\]

If computed on-the-fly from i’s neighbor loop only without storing q_ji, validate pair antisymmetry via diagnostic.

### 5.4 Conservation Equation

Per particle:

\[
R_i = \sum_{j\in\mathcal N(i)} q_{ij}=0.
\]

New φ equation. Still scalar φ solve, but built on pairwise current flux unlike old div-grad.

### 5.5 Joule Heat Edge Energy Form

Time-averaged Joule power per pair:

\[
P_{ij}=\frac{1}{2}C_{ij}e_{ij}^2.
\]

Half pair power to each particle; volumetric heat source on particle i:

\[
Q_i^{edge}=\frac{1}{V_i}\sum_j \frac{1}{2}P_{ij}
=\frac{1}{V_i}\sum_j \frac{1}{4}C_{ij}e_{ij}^2.
\]

Total power:

\[
P_{edge}=\sum_i Q_i^{edge}V_i.
\]

Acceptance must include:

```text
I × 2  =>  P_edge × 4
P_edge finite
Q_edge strong outside, weak inside, no obvious checkerboard noise
```

Note: graph edge energy may not equal physical Joule power—see power-fix plan; use reconstructed power for calibration in production.

### 5.6 Particle E/J Reconstruction

Primary current is \(q_{ij}\); for visualization and later \(A_{ind}=K[J]\), reconstruct particle \(E_i,J_i\).

Along pair direction:

\[
E_i\cdot\hat{r}_{ij}\approx \frac{e_{ij}}{|r_{ij}|}.
\]

Least squares:

\[
M_i E_i = b_i,
\]

\[
M_i=\sum_j w_{ij}\hat r_{ij}\otimes\hat r_{ij},
\]

\[
b_i=\sum_j w_{ij}\left(\frac{e_{ij}}{|r_{ij}|}\right)\hat r_{ij}.
\]

Then:

\[
J_i^{edge}=\sigma_iE_i.
\]

If \(M_i\) ill-conditioned, fallback and output:

```text
EdgeReconCondition
EdgeReconFallback
```

---

## 6. New Modules and File Recommendations

### 6.1 New Files

Under:

```text
tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/
```

Add:

```text
electromagnetic_ophelie_edge_current.h
electromagnetic_ophelie_edge_current.hpp
```

Optional split:

```text
electromagnetic_ophelie_edge_current_fields.h
electromagnetic_ophelie_edge_current_diagnostics.h
electromagnetic_ophelie_edge_joule_heat.h
```

First stage: limit file count; one main header and implementation file.

### 6.2 Register Fields

In field names and register fields add:

```text
EdgeCurrentResidual
EdgeCurrentResidualAbs
EdgeCurrentResidualLevel0
JouleHeatEdge
PowerEdgeParticle
EEdgeRecon
JEdgeRecon
EdgeReconCondition
EdgeReconFallback
EdgeDropRms
EdgeCurrentRms
```

Or with prefix per naming style:

```text
OphelieEdgeCurrentResidual
OphelieJouleHeatEdge
OphelieEEdgeRecon
OphelieJEdgeRecon
```

Keep VTP readable.

### 6.3 New CLI

Add:

```text
--ophelie-current-form=particle-gradient|edge-current
```

Meaning:

```text
particle-gradient : old div-grad baseline particle E/J/JouleHeat;
edge-current      : new edge-flux current / JouleHeatEdge / JEdgeRecon.
```

Preserve old:

```text
--phi-projection-operator=div-grad|compatible-div-grad|edge-flux
```

New mainline should avoid ambiguous switch combinations. Logic:

```text
if --ophelie-current-form=edge-current:
    phi projection operator must be edge-flux
    JouleHeat output uses JouleHeatEdge
    production acceptance uses edge-current metrics
else:
    use old particle-gradient baseline
```

If conflicting parameters, print warning or reject.

---

## 7. CK Class Design Recommendations

### 7.1 Pair Conductance Helper

Extract unified helper so φ LHS, RHS, q_ij, JouleHeatEdge share same \(C_{ij}\).

Add:

```cpp
inline Real computeOphelieEdgePairConductance(
    const Real sigma_i,
    const Real sigma_j,
    const Vecd& r_ij,
    const Real dW_ij,
    const Real Vol_j,
    const Real h,
    const OphelieParameters& params);
```

Or reuse existing legacy pairwise Laplace function. Key:

```text
Do not copy a similar-looking but not identical C_ij.
```

### 7.2 ComputeOphelieEdgeCurrentResidualCK

Compute level0 and post-phi edge current residual.

Pseudocode:

```cpp
class ComputeOphelieEdgeCurrentResidualCK
    : public LocalDynamics,
      public Interaction<Inner<>>
{
public:
    class InteractKernel
    {
    public:
        void interact(size_t index_i, Real dt = 0.0)
        {
            Real residual_i = 0.0;
            Real abs_sum_i = 0.0;

            const Vecd xi = pos_[index_i];
            const Real phi_i = phi_[index_i];
            const Vecd A_i = a_real_[index_i];
            const Real sigma_i = sigma_[index_i];

            for each neighbor j:
                Vecd rij = xi - xj;
                Real Cij = computeOphelieEdgePairConductance(...);
                Vecd Aavg = 0.5 * (A_i + A_j);
                Real edge_drop = (phi_i - phi_j) + omega * Aavg.dot(rij);
                Real qij = -Cij * edge_drop;
                residual_i += qij;
                abs_sum_i += abs(qij);

            edge_residual_[index_i] = residual_i;
            edge_residual_abs_[index_i] = abs_sum_i;
        }
    };
};
```

### 7.3 ComputeOphelieEdgeJouleHeatCK

Edge energy Joule heat.

Pseudocode:

```cpp
class ComputeOphelieEdgeJouleHeatCK
    : public LocalDynamics,
      public Interaction<Inner<>>
{
public:
    class InteractKernel
    {
    public:
        void interact(size_t index_i, Real dt = 0.0)
        {
            Real Qi_power = 0.0;

            for each neighbor j:
                Real Cij = ...;
                Real edge_drop = ...;
                Real Pij = 0.5 * Cij * edge_drop * edge_drop;

                // half of pair power assigned to particle i
                Qi_power += 0.5 * Pij;

            joule_heat_edge_[index_i] = Qi_power / Vol_i;
            power_edge_particle_[index_i] = Qi_power;
        }
    };
};
```

Note double-counting: i-loop and j-loop both visit same pair; half to i and half to j is reasonable. Total:

```text
P_total_edge = sum_i power_edge_particle_i
```

Should match pair total energy.

### 7.4 ReconstructOphelieEdgeElectricCurrentCK

Least-squares reconstruct particle \(E_i,J_i\) from edge drops.

Pseudocode:

```cpp
class ReconstructOphelieEdgeElectricCurrentCK
    : public LocalDynamics,
      public Interaction<Inner<>>
{
public:
    class InteractKernel
    {
    public:
        void interact(size_t index_i, Real dt = 0.0)
        {
            Matd M = Matd::Zero();
            Vecd b = Vecd::Zero();

            for each neighbor j:
                Vecd rij = xi - xj;
                Real r = norm(rij);
                Vecd e = rij / r;
                Real edge_drop = ...;
                Real directional_E = edge_drop / r;
                Real wij = choose_reconstruction_weight(Cij or kernel weight);

                M += wij * outer(e, e);
                b += wij * directional_E * e;

            if invertible(M):
                E_i = inverse(M) * b;
                fallback = 0;
            else:
                E_i = Vecd::Zero();
                fallback = 1;

            J_i = sigma_i * E_i;
        }
    };
};
```

---

## 8. New Acceptance / Diagnostics

### 8.1 Edge-Current Acceptance

Add:

```text
edge_res_l2_level0
edge_res_l2_post_phi
edge_res_red = edge_res_l2_level0 / edge_res_l2_post_phi
q_antisymmetry_error
P_total_edge
P_total_recon
P_edge_over_recon
Q_edge_max
Q_edge_mean
Q_edge_max_over_mean
J_edge_max
E_edge_max
recon_fallback_fraction
```

### 8.2 Initial Pass Criteria

First stage not too strict:

```text
edge_res_red > 100
q_antisymmetry_error < 1e-8 relative
P_total_edge finite
I×2 → P_total_edge×4
Q_edge field finite
Q_edge strong outside, weak inside
recon_fallback_fraction < 5%
P_edge / P_recon within 0.5–2 as initial target; 0.25–4 OK for diagnostic
```

Production calibration should use `P_total_recon`; see power-fix plan.

### 8.3 Do Not Use Old Particle divJ as Primary Edge-Current Fail Metric

For edge-current formulation, old particle `divJ_L2_red` is reference only. Primary:

```text
edge_res_red
P_edge / P_recon
Q_edge
J_edge_reconstruction quality
```

---

## 9. Stage 1 Development Steps

### Step 1: Organize Old Route, Create Archive README

Output:

```text
docs/ophelie/archive_phi_divgrad_residual_floor/README.md
```

Commit once done.

### Step 2: Add Edge-Current Fields and CLI

Implement:

```text
--ophelie-current-form=particle-gradient|edge-current
```

Default:

```text
particle-gradient
```

### Step 3: Implement Pair Conductance Helper

Must reuse legacy pairwise actual weights. Cursor must find current `legacy-pairwise` and `legacy-flux` pair weight formulas—do not rewrite by guess.

### Step 4: Implement EdgeCurrentResidualCK

Residual only first; no JouleHeat.

Acceptance:

```text
after edge-flux phi solve, edge_res_red > 100
```

### Step 5: Implement EdgeJouleHeatCK

Output `JouleHeatEdge` and `PowerEdgeParticle`.

Acceptance:

```text
P_total_edge finite
current scaling: I×2 => P_edge×4
```

### Step 6: Implement Edge E/J Reconstruction CK

Output:

```text
EEdgeRecon
JEdgeRecon
EdgeReconFallback
EdgeReconCondition
```

### Step 7: French Edge-Current Case

Add or extend:

```text
test_3d_ophelie_french_reduced --ophelie-current-form=edge-current
```

Three-way comparison:

```text
particle-gradient baseline
edge-current new route
compatible diagnostic
```

### Step 8: VTP Visualization

VTP at minimum:

```text
JouleHeatEdge
JEdgeRecon
EEdgeRecon
EdgeCurrentResidual
```

---

## 10. Stage 2: A_ind Self-Induction Route — Not Immediate

After Stage 1 success:

```text
A_ind = K[JEdgeRecon]
```

One-way diagnostic first, then Picard:

```text
A^k = A_src + A_ind^k
solve edge-current φ/q
reconstruct JEdge
A_ind^{k+1}=K[JEdge]
under-relax
repeat
```

Do not do A_ind in Stage 1.

---

## 11. Stage 3: Full J-Phi Coupled System — Not Immediate

If edge-current + Picard insufficient, consider full coupled system:

\[
\begin{bmatrix}
I+\sigma\omega K & \sigma G \\
D & 0
\end{bmatrix}
\begin{bmatrix}
J\\
\phi
\end{bmatrix}
=
\begin{bmatrix}
-\sigma\omega A_{src}\\
0
\end{bmatrix}.
\]

Highest cost: J three components, phi scalar, global Biot K[J], complex nonsymmetric matrix-free GMRES, GPU O(N^2) or far-field acceleration.

Not next sprint target.

---

## 12. Stop / Go Criteria

### Go: Continue Edge-Current Production

Satisfy:

```text
1. edge_res_red > 100;
2. P_recon current scaling correct (see power-fix plan);
3. JouleHeatEdgeRecon distribution physically reasonable;
4. P_recon used for calibration; P_graph diagnostic;
5. JEdgeRecon not noisy;
6. stable on SYCL GPU backend;
7. more self-consistent than div-grad baseline without breaking French 50 kW calibration.
```

### Stop: Freeze Edge-Current

Any serious failure:

```text
1. P_edge and P_recon differ by order of magnitude and unexplained;
2. JEdgeRecon widespread fallback or strong noise;
3. JouleHeatEdge nonphysical hotspots/checkerboard;
4. must introduce external mesh/CSR/CPU loop to continue;
5. SYCL CK implementation cost out of control;
6. cannot maintain current scaling.
```

If stop, return to:

```text
particle-gradient baseline + A_ind diagnostic
```

But no longer blindly tune div-grad residual.

---

## 13. Short Execution Summary for Cursor

```text
We now officially start edge-flux solver development.
Edge in edge-flux is SPH neighbor pair i-j, not external grid.
All implementation must stay SPHinXsys SYCL CK style.

Old div-grad baseline preserved as fallback/reference; no longer primary phi_eq_res goal.
compatible-div-grad, Neumann, grad correction, old edge-flux phi-only all archived as diagnostic.

New mainline:
1. Define edge_drop and q_ij with SPH pair conductance C_ij;
2. Solve φ with sum_j q_ij=0;
3. Get JouleHeatEdge from q_ij / edge energy (graph diagnostic);
4. Least-squares reconstruct EEdge/JEdge; use P_recon for calibration;
5. Later JEdge for A_ind=K[J].

Stage 1 forbidden: A_ind, thermal coupling, external grid, CSR, host-only loop.

Stage 1 acceptance:
edge_res_red>100, P_recon finite and calibrated, I×2→P_recon×4, Q distribution reasonable, JEdgeRecon not noisy.
```

---

## 14. Final Roadmap

```text
Stage 0: Archive old phi-divgrad diagnostics
Stage 1: SPH edge-current φ/q/JouleHeatEdge/JEdgeRecon
Stage 2: A_ind = K[JEdge] one-way diagnostic
Stage 3: Picard self-induction OPHELIE-like solver
Stage 4: TEAM7 / reference benchmark (see docs/TEAM7-reference/)
Stage 5: JouleHeat → thermal one-way
Stage 6: σ(T) and thermal-electromagnetic coupling
```

---

## 15. Key Conclusion

Current issue is not “demo insufficient” but `div-grad baseline` no longer suitable as true OPHELIE-like solver mainline. Edge-flux worth trying but must be strictly:

```text
SPH pairwise edge-current formulation
```

Not external grid methods.

If Stage 1 succeeds, closer to French OPHELIE core:

```text
Current is core variable; current continuity is the equation; Joule heat from same current discretization.
```

Mainline for developing a usable OPHELIE solver.
