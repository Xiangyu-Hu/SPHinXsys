# OPHELIE Complex Edge-Flux Solver Plan

**Purpose.** This document records the decision to move from the real-only edge-flux prototype to a full complex edge-flux formulation, explains why a short A_ind phase audit is still useful, and gives Cursor a concrete implementation plan for code changes, tests, and acceptance criteria.

---

## 0. Executive decision

We should build the **complex edge-flux OPHELIE-like solver**.

The current real-only edge-flux route has already proven the key SPH/SYCL idea:

```text
SPH neighbor-pair edge drop
→ edge current residual
→ edge-reconstructed E/J
→ JouleHeatEdgeRecon
→ 50 kW calibration
```

This is now much stronger than the old `div-grad` baseline. The remaining issue is that the real-only formulation only covers the chain:

```text
A_real → phi_imag → E_imag/J_imag → Joule heat
```

For true OPHELIE-like self-induction, the induced current must generate a complex vector potential:

```text
J_real  → A_ind_real
J_imag  → A_ind_imag
```

and both components must feed back into the electric field. Therefore, if accuracy is the goal, the correct next solver is:

```text
complex edge-flux:
  A_real, A_imag
  phi_real, phi_imag
  E_real, E_imag
  J_real, J_imag
  JouleHeat = 0.5 * sigma * (|E_real|^2 + |E_imag|^2)
```

The A_ind phase audit is **not** a substitute for complex edge-flux. It is a short, mandatory sanity test to verify the sign convention, phase convention, and Biot-Savart source phase before implementing Picard self-induction. It should not delay the complex solver; it should be folded into the first complex-edge-flux sprint.

---

## 1. Why not continue with real-only edge-flux?

The real-only edge-flux formulation has worked well for externally prescribed real coil vector potential:

```text
A_src_real → phi_imag → J_imag → JouleHeatEdgeRecon
```

This is sufficient for a reduced demonstration where the source field is prescribed and self-induction is ignored. However, once we compute:

```text
A_ind = K[J]
```

we cannot assume that the feedback remains in `A_real`. In phasor notation, the Biot-Savart operator itself does not multiply by `i`; it preserves the phase of the current phasor:

```text
A = K[J]
A_ind_real = K[J_real]
A_ind_imag = K[J_imag]
```

If the coil source is stored as `A_src_real`, the externally driven eddy current is mainly `J_imag`. Therefore the first induced vector potential will generally be mainly:

```text
A_ind_imag = K[J_imag]
```

A solver that only reads `A_real` in the edge drop cannot fully feed this back. It may incorrectly conclude that A_ind has no effect, simply because the feedback is in the imaginary component which the current edge-drop ignores.

Therefore, for a real OPHELIE-like solver, **complex edge-flux is required**.

---

## 2. Why still do an A_ind phase audit?

The phase audit is not about avoiding complexity. It is about preventing us from implementing a complex solver with the wrong sign or wrong phase convention.

The audit answers these questions before Picard iteration:

1. Does `K[JEdgeRecon]` write the induced vector potential into `A_ind_real`, `A_ind_imag`, or both?
2. Is the dominant induced component consistent with the phasor convention used in the edge-drop formula?
3. Are the stored field names and VTP outputs aligned with the mathematical convention?
4. Will a real-only active-A test be meaningful, or does the feedback immediately require `A_imag`?
5. Does `JEdgeRecon` produce an A_ind magnitude comparable to previous diagnostic estimates, e.g. `A_ind/A_src ~ O(0.3–0.4)`?

The phase audit should be a short diagnostic step, not a separate research branch.

**Decision:** perform the phase audit as Stage 3.0, then proceed directly to complex edge-flux implementation.

---

## 3. Harmonic convention and complex equations

Use the same physical starting point as the French OPHELIE formulation:

\[
\mathbf j = -\sigma\left(\nabla V + \frac{\partial \mathbf A}{\partial t}\right),
\qquad
\frac{\partial \mathbf A}{\partial t}= i\omega \mathbf A.
\]

Let:

\[
\mathbf A = \mathbf A_r + i\mathbf A_i,
\qquad
\phi = \phi_r + i\phi_i.
\]

Then:

\[
\mathbf E = -\nabla\phi - i\omega \mathbf A.
\]

Separating real and imaginary components:

\[
\mathbf E_r = -\nabla\phi_r + \omega \mathbf A_i,
\]

\[
\mathbf E_i = -\nabla\phi_i - \omega \mathbf A_r.
\]

and:

\[
\mathbf J_r = \sigma \mathbf E_r,
\qquad
\mathbf J_i = \sigma \mathbf E_i.
\]

The time-averaged Joule heat for peak-value phasors is:

\[
Q = \frac{1}{2}\sigma\left(|\mathbf E_r|^2 + |\mathbf E_i|^2\right)
  = \frac{1}{2\sigma}\left(|\mathbf J_r|^2 + |\mathbf J_i|^2\right).
\]

This generalizes the current real-only implementation, which effectively computes only:

\[
Q_i = \frac{1}{2}\sigma |\mathbf E_i|^2.
\]

---

## 4. Complex edge-drop definitions

The current edge-flux code uses the imaginary-chain edge drop:

\[
e_{ij}^{(i)} = (\phi_{i,j} - \phi_{i,i}) + \omega \bar{\mathbf A}_{r,ij}\cdot(\mathbf x_j-\mathbf x_i).
\]

This is correct for:

\[
\mathbf E_i = -\nabla\phi_i - \omega \mathbf A_r.
\]

For the real component:

\[
\mathbf E_r = -\nabla\phi_r + \omega \mathbf A_i.
\]

The corresponding edge drop is:

\[
e_{ij}^{(r)} = (\phi_{r,j} - \phi_{r,i}) - \omega \bar{\mathbf A}_{i,ij}\cdot(\mathbf x_j-\mathbf x_i).
\]

Therefore the complex edge-flux solver needs two scalar edge-drop solves:

```text
Imaginary equation:
  edge_drop_imag = (phi_imag_j - phi_imag_i)
                   + omega * A_real_avg · (x_j - x_i)

Real equation:
  edge_drop_real = (phi_real_j - phi_real_i)
                   - omega * A_imag_avg · (x_j - x_i)
```

Both use the same SPH pair conductance `C_ij`.

---

## 5. Sign tests that Cursor must implement

Do not infer the complex edge-drop sign from the old vector-divergence MMS `sign_alpha`. The edge-flux sign must be verified by direct edge-drop-zero tests.

### 5.1 Imaginary-chain sign test

Set:

\[
\mathbf A_r = \mathbf A_0,
\qquad
\phi_i = -\omega \mathbf A_0\cdot\mathbf x.
\]

Then:

\[
\phi_{i,j}-\phi_{i,i} = -\omega \mathbf A_0\cdot(\mathbf x_j-\mathbf x_i),
\]

so:

\[
e_{ij}^{(i)}=0.
\]

Test output:

```text
edge_drop_imag_zero_l2
edge_drop_imag_zero_linf
edge_drop_imag_zero_mean_abs
```

### 5.2 Real-chain sign test

Set:

\[
\mathbf A_i = \mathbf A_0,
\qquad
\phi_r = +\omega \mathbf A_0\cdot\mathbf x.
\]

Then:

\[
\phi_{r,j}-\phi_{r,i} = +\omega \mathbf A_0\cdot(\mathbf x_j-\mathbf x_i),
\]

so:

\[
e_{ij}^{(r)}=0.
\]

Test output:

```text
edge_drop_real_zero_l2
edge_drop_real_zero_linf
edge_drop_real_zero_mean_abs
```

Initial tolerance:

```text
l2   < 1e-6 to 1e-8, depending on Real precision and GPU backend
linf < 1e-5 initially acceptable for float/SYCL tests
```

If either sign test fails, do not proceed to A_ind or Picard.

---

## 6. Immediate P0 fixes before complex implementation

### 6.1 Fix `test_3d_ophelie_french_aind_diagnostic` CLI

Current issue: `test_3d_ophelie_french_aind_diagnostic` does not parse general OPHELIE CLI options. Therefore:

```text
--ophelie-current-form=edge-flux
--ophelie-use-a-total-for-edge-flux
--coil-current-scale
--target-power
```

may be ignored.

Cursor must connect the same OPHELIE CLI path used in `test_3d_ophelie_french_reduced.cpp`:

```cpp
OphelieParameters params;
OphelieFrenchReducedCaseParams french;
applyFrenchReducedDefaults(params, french);

OphelieTestCliOptions cli_options;
StdVec<std::string> filtered_arguments =
    filterOphelieTestCommandLine(ac, av, params, cli_options);

if (cli_options.literature_mode)
{
    OphelieFrenchLiteratureProfile literature_profile;
    applyFrenchLiteratureMode(params, cli_options, literature_profile);
}

finalizeOphelieCurrentFormConfiguration(params);
syncFrenchReducedToParameters(french, params);
applyOphelieCoilCurrentScale(french, params);
logOphelieFinalParams(params, cli_options);
```

The exact helper names may differ; Cursor should follow the current French reduced test implementation.

### 6.2 Fix edge-flux RHS setup in A_ind diagnostic

Current issue: in `runFrenchReducedAIndOneWayDiagnostic()`, the edge-flux branch calls `solvePhiImagWithCurrentRhs()` without clearly setting up and finalizing the edge-flux RHS first.

Cursor must ensure the sequence is:

```cpp
zero_phi.exec();

setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, glass_inner, names, params);
finalizeOpheliePhiImagRhsHost(particles, names, params, nullptr, french.dp, nullptr);

(void)solvePhiImagWithCurrentRhs<ExecutionPolicy>(glass_body, glass_inner, names, params);

applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, glass_inner, names, params);
execOphelieEdgeFluxPostPhiPipeline<ExecutionPolicy>(glass_body, glass_inner, names, params);
syncOphelieEdgeReconToPrimaryEJQ<ExecutionPolicy>(glass_body, names);
```

The key rule is: **never call `solvePhiImagWithCurrentRhs()` before the current RHS field is initialized and finalized.**

---

## 7. Stage 3 implementation plan: complex edge-flux

### Stage 3.0 — Phase audit and sign validation

Deliverables:

1. Fixed A_ind diagnostic CLI.
2. Correct edge-flux RHS setup in A_ind diagnostic.
3. A_ind phase output:

```text
A_ind_real_norm
A_ind_imag_norm
A_ind_real/A_src_real
A_ind_imag/A_src_real
B_ind_real_norm
B_ind_imag_norm
B_ind_real/B_src_real
B_ind_imag/B_src_real
max_A_ind_real
max_A_ind_imag
max_B_ind_real
max_B_ind_imag
```

4. Direct edge-drop sign tests for real and imaginary chains.

Acceptance:

```text
current_form=edge-flux is actually used in A_ind diagnostic
J source = JEdgeRecon
edge_drop_zero tests pass
A_ind real/imag phase is clearly reported
```

This is not a stopping point; it is the setup for complex edge-flux.

---

### Stage 3.1 — Register complex edge-flux fields

Required fields, or their current equivalents if already present:

```text
PhiReal
PhiImag
ATotalReal
ATotalImag
ASrcReal
ASrcImag
AIndReal
AIndImag
EEdgeReconReal
EEdgeReconImag
JEdgeReconReal
JEdgeReconImag
JouleHeatEdgeReconComplex
EdgeFluxResidualReal
EdgeFluxResidualImag
EdgeDropRealStats
EdgeDropImagStats
```

For compatibility with downstream modules:

```text
EImag / JImag / JouleHeat
```

may remain aliases of the dominant production fields, but documentation must state clearly what they mean in edge-flux mode.

Recommended production outputs in VTP:

```text
EEdgeReconReal
EEdgeReconImag
JEdgeReconReal
JEdgeReconImag
JouleHeatEdgeReconComplex
PhiReal
PhiImag
ATotalReal
ATotalImag
AIndReal
AIndImag
```

---

### Stage 3.2 — Refactor edge-flux kernels to accept component configuration

Current edge-flux kernels are mostly hardwired to the imaginary chain:

```text
A_real → phi_imag → E_imag/J_imag
```

Cursor should implement a generic component interface instead of duplicating large code blocks.

Suggested component descriptor:

```cpp
struct OphelieEdgeFluxComponent
{
    std::string phi_field;
    std::string rhs_field;
    std::string lhs_field;
    std::string active_a_field;
    std::string edge_residual_field;
    std::string e_recon_field;
    std::string j_recon_field;

    // +1 for imag equation using A_real:
    //   edge_drop = phi_j - phi_i + omega A_real · dx
    // -1 for real equation using A_imag:
    //   edge_drop = phi_j - phi_i - omega A_imag · dx
    Real a_sign;
};
```

Then use:

```text
Imag component:
  phi_field     = PhiImag
  active_a      = ATotalReal
  a_sign        = +1
  e_recon       = EEdgeReconImag
  j_recon       = JEdgeReconImag

Real component:
  phi_field     = PhiReal
  active_a      = ATotalImag
  a_sign        = -1
  e_recon       = EEdgeReconReal
  j_recon       = JEdgeReconReal
```

Generic edge-drop:

```cpp
edge_drop = (phi_j - phi_i)
          + component.a_sign * omega * a_avg.dot(xj_minus_xi);
```

---

### Stage 3.3 — Solve both scalar systems

The complex edge-flux solver should run two scalar solves:

```text
1. Imaginary solve:
   ATotalReal → PhiImag → EEdgeReconImag/JEdgeReconImag

2. Real solve:
   ATotalImag → PhiReal → EEdgeReconReal/JEdgeReconReal
```

Important implementation point:

Current GMRES and field names may be hardwired to `PhiImag`, `PhiRhsImag`, `PhiLhsImag`. Cursor must either:

1. Refactor the scalar edge-flux solve to accept a component descriptor and field names; or
2. Duplicate the minimal field set for `PhiReal` and add `solvePhiRealWithCurrentRhs` as a temporary bridge.

Preferred solution: field-name parameterization.

Acceptance:

```text
Imag solve alone reproduces current H edge-flux v4 behavior.
Real solve with ATotalImag=0 gives near-zero PhiReal/JReal.
Direct sign tests pass for both components.
```

---

### Stage 3.4 — Complex Joule heat

Compute:

\[
Q = \frac{1}{2}\sigma(|E_r|^2 + |E_i|^2).
\]

In code:

```cpp
joule_heat_complex[i] = 0.5 * sigma[i] *
    (E_edge_real[i].squaredNorm() + E_edge_imag[i].squaredNorm());
```

Total power:

```cpp
P_complex = sum_i joule_heat_complex[i] * Vol_i;
```

The current `P_recon_edge` should become:

```text
P_recon_complex
```

For real-only source with `ATotalImag=0`, this should match the current H edge-flux v4 result.

Acceptance:

```text
P_complex ≈ old P_recon_edge when AImag=0
I×2 → P_complex×4
Q_complex nonfinite_count=0
Q_complex negative_count=0
Q_outer/center > 1 for French reduced case
```

---

### Stage 3.5 — A_ind one-way with complex fields

After complex edge-flux works without self-induction:

1. Solve with:

```text
ATotalReal = ASrcReal
ATotalImag = 0
```

2. Reconstruct:

```text
JEdgeReconReal
JEdgeReconImag
```

3. Compute:

```text
AIndReal = K[JEdgeReconReal]
AIndImag = K[JEdgeReconImag]
```

4. Build:

```text
ATotalReal = ASrcReal + AIndReal
ATotalImag = ASrcImag + AIndImag
```

5. Re-solve complex edge-flux once.

Output:

```text
P_complex_coil_only
P_complex_total_A
AIndReal/ASrcReal
AIndImag/ASrcReal
BIndReal/BSrcReal
BIndImag/BSrcReal
Qmax change
Qmean change
edge_res_red_real
edge_res_red_imag
```

Do not do Picard yet.

---

### Stage 3.6 — Picard self-induction prototype

Only after one-way complex A_ind is stable:

```text
A_total^k = A_src + A_ind^k
solve complex edge-flux
reconstruct J_real/J_imag
A_ind^{k+1} = K[J]
under-relax A_ind
repeat
```

Suggested initial controls:

```text
--ophelie-aind-picard=1
--ophelie-aind-picard-max-iter=5
--ophelie-aind-relax=0.3
--ophelie-aind-rel-tol=1e-3
```

Acceptance:

```text
A_ind update norm decreases
P_complex remains finite
edge residuals remain good
Q distribution remains physical
no NaN/Inf
```

---

## 8. Tests Cursor must add or update

### 8.1 Complex edge-drop sign test

Extend existing `test_3d_ophelie_edge_flux_sign`:

```text
imag chain: AReal=A0, PhiImag=-omega*A0·x → edge_drop_imag≈0
real chain: AImag=A0, PhiReal=+omega*A0·x → edge_drop_real≈0
```

### 8.2 Complex uniform-field power test

Extend `test_3d_ophelie_edge_flux_power_uniform_field`:

Case A:

```text
AReal=0, AImag=0, PhiReal=-E0·x, PhiImag=0
```

Expected:

```text
E_real=E0
P_exact=0.5*sigma*|E0|^2*V
```

Case B:

```text
Phi=0, AReal=A0, AImag=0
```

Expected:

```text
E_imag=-omega*A0
P_exact=0.5*sigma*omega^2*|A0|^2*V
```

Case C:

```text
Phi=0, AReal=0, AImag=A0
```

Expected:

```text
E_real=+omega*A0
P_exact=0.5*sigma*omega^2*|A0|^2*V
```

### 8.3 French complex regression

Run:

```text
--ophelie-current-form=edge-flux
--ophelie-edge-flux-complex=1
--target-power=50000
```

Expected with no A_ind:

```text
P_complex≈50000
J_real≈0 or very small
J_imag comparable to old H v4 JEdgeRecon
edge_res_imag good
edge_res_real near trivial if AImag=0
```

### 8.4 A_ind phase audit test

Run true edge-flux current form and output A_ind real/imag ratios. Do not accept results unless CLI log confirms:

```text
current_form=edge-flux
J source=JEdgeRecon complex fields
```

---

## 9. What Cursor must not do now

Do not do these until the complex edge-flux and one-way A_ind tests are complete:

```text
1. No thermal coupling.
2. No A_ind Picard before one-way complex A_ind.
3. No further div-grad residual repair.
4. No compatible-div-grad or Neumann production revival.
5. No P_graph kappa calibration.
6. No monolithic J-phi solve yet.
7. No external mesh, CSR graph, host adjacency, or CPU-only edge loops.
```

The solver must remain SPH neighbor-pair based and SPHinXsys SYCL compatible.

---

## 10. Required log language to avoid confusion

In edge-flux complex mode, logs should avoid old div-grad gate language.

Use:

```text
edge_flux_complex_continuity_passed
edge_res_red_real
edge_res_red_imag
q_antisym_rel_l2_real
q_antisym_rel_l2_imag
P_complex_recon
Q_complex_spatial_gate
```

Avoid or clearly mark:

```text
divJ_L2_red                     diagnostic only
particle_divJ_continuity_passed diagnostic only
P_graph_edge                    diagnostic only, not physical Joule power
```

---

## 11. Direct instruction summary for Cursor

```text
We will proceed to complex edge-flux.

Do not treat the A_ind phase audit as an alternative to complex edge-flux.
It is Stage 3.0 sanity validation for signs and field phases.

Immediate tasks:
1. Fix test_3d_ophelie_french_aind_diagnostic CLI so --ophelie-current-form=edge-flux is actually parsed.
2. Fix runFrenchReducedAIndOneWayDiagnostic edge branch so RHS is setup/finalized before solvePhiImagWithCurrentRhs.
3. Add A_ind phase audit outputs: AIndReal/Imag and BIndReal/Imag norms and ratios.
4. Add real-chain edge-drop sign test.
5. Implement complex edge-flux fields and generic component kernels:
   imag: PhiImag driven by ATotalReal with +omega sign;
   real: PhiReal driven by ATotalImag with -omega sign.
6. Compute complex Joule heat:
   Q = 0.5*sigma*(|EReal|^2 + |EImag|^2).
7. Verify complex edge-flux reduces to current H v4 when AImag=0.
8. Only after one-way complex A_ind is stable, consider Picard.
```

---

## 12. Final position

The user is right to question why we would remain in real-only edge-flux if the physical system is complex. The correct answer is:

```text
We should not remain real-only.
We should move to complex edge-flux now.
```

The A_ind phase audit is still useful, but only as a quick validation step for sign and phase conventions. It should be part of the complex-edge-flux transition, not a reason to postpone it.

The new main line is:

```text
SPH/SYCL complex edge-flux
→ complex E/J reconstruction
→ complex Joule heat
→ one-way A_ind = K[J_complex]
→ Picard self-induction
```

This is the route most consistent with the French OPHELIE spirit, where current and scalar potential are central unknowns and the vector potential is obtained by Biot-Savart coupling.
