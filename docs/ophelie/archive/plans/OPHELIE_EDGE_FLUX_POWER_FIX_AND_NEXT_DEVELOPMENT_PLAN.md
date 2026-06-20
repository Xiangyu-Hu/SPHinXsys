# OPHELIE Edge-Flux Solver: Current Status, Power Definition Issue, and Next Development Plan

> Purpose: This document consolidates this round’s analysis conclusions, current blockers, code structure review, power/JouleHeat definition fix plan, and next-step development instructions for Cursor.  
> Current stage focus is not continuing to tune old `div-grad baseline` residuals, but making the new SPH/SYCL `edge-flux` solver a closed **edge current → reconstructed E/J → Joule heat → power calibration** chain.

---

## 1. Overall Assessment

This round’s edge-flux direction is correct. Code largely keeps SPHinXsys/SYCL particle neighbor-pair form—not external finite-volume grid, CSR graph, host-only adjacency, or discretization outside SPH kernel-sum.

But the blocker is clear:

**Edge-flux φ equation / edge residual already solves well; what failed is treating graph edge energy directly as physical Joule power.**

Current H production key results:

```text
P_total_edge  = 50000.2 W
P_total_recon = 0.429548 W
P_edge_over_recon = 116402
```

This indicates `P_total_edge` from

\[
P_{ij}=\frac12 C_{ij}e_{ij}^2
\]

is not trustworthy physical Joule power. It is a graph/Laplace discrete energy quantity and cannot serve as the 50 kW calibration target.

Physical power for Stage 1 should come from edge-drop least-squares reconstructed E/J:

\[
Q_i^{\mathrm{edge-recon}}
=
\frac12\sigma_i |E_i^{\mathrm{edge}}|^2,
\]

\[
P_{\mathrm{recon}}
=
\sum_i Q_i^{\mathrm{edge-recon}}V_i.
\]

Therefore next steps do not reject edge-flux but fix **power ownership / calibration target / JouleHeat output**.

---

## 2. What This Round’s Edge-Flux Already Completed

### 2.1 Edge-Flux Is Already SPH/SYCL Form

New code mainly in:

```text
electromagnetic_ophelie_edge_flux.h
electromagnetic_ophelie_edge_flux.hpp
```

Core classes:

```cpp
ComputeOphelieEdgeFluxPhiRhsFromASrcCK
ComputeOphelieEdgeFluxResidualCK
ComputeOphelieEdgeFluxJouleHeatCK
ReconstructOphelieEdgeFluxElectricCurrentCK
```

Based on:

```cpp
Interaction<Inner<>>
InteractionDynamicsCK
DelegatedData(ex_policy)
SPH neighbor loop
```

Requirements met:

```text
Must use SPH particles;
Must use SPHinXsys cell linked list / inner relation;
Must use Inner neighbor pairs;
Must run on SYCL CK backend;
No external mesh;
No CSR graph matrix;
No host-only adjacency;
No CPU-only pair loop.
```

Direction is correct; continue preserving it.

---

### 2.2 New Edge-Drop Form Already Wired

Current edge-drop:

\[
e_{ij}
=
(\phi_j-\phi_i)
+
\omega \bar{A}_{ij}\cdot(x_j-x_i),
\]

Where:

\[
\bar{A}_{ij}=\frac{A_i+A_j}{2}.
\]

Code:

```cpp
return (phi_j - phi_i) + omega * a_avg.dot(xj_minus_xi);
```

Not old legacy-flux:

\[
g_{ij}\cdot(A_i-A_j),
\]

but pairwise electromotive drop. Correct direction.

Emphasize:

```text
legacy-flux RHS is divergence-form φ equation RHS;
new edge-flux is pairwise current flux / electromotive drop formulation;
do not conflate them.
```

---

### 2.3 Edge Residual Solves Well

H production log:

```text
phi_eq_res_vol = 8.0e-5
edge_res_red_l2 ≈ 10920
```

Edge-flux equation is easy to solve with strong residual reduction. Versus old div-grad baseline:

```text
div-grad baseline phi_eq_res ≈ 0.486
```

Edge-flux is clearly better at φ equation level.

Therefore edge-flux did not fail overall. More precisely:

```text
edge-flux φ equation / edge residual: PASS
edge-flux physical Joule power closure: FAIL
edge-flux production: FAIL for now
```

---

## 3. Where It Actually Blocks

### 3.1 Key Failed Metrics

H production:

```text
P_total_edge  = 50000.2 W
P_total_recon = 0.429548 W
P_edge_over_recon = 116402
```

Graph edge energy calibrated to 50 kW; reconstructed physical power only 0.43 W—about 116000× gap.

Not small error or acceptable reconstruction error:

**`P_total_edge` cannot serve as physical Joule power.**

---

### 3.2 Why `P_total_edge` Is So Wrong

Graph energy uses:

\[
P_{ij}
=
\frac12 C_{ij} e_{ij}^2,
\]

Each particle gets half pair energy:

\[
Q_i^{\mathrm{edge}}
=
\frac{1}{V_i}
\sum_j \frac12 P_{ij}
=
\frac{1}{V_i}
\sum_j \frac14 C_{ij}e_{ij}^2.
\]

Code similar to:

```cpp
const Real p_ij = Real(0.5) * c_ij * edge_drop * edge_drop;
power_i += Real(0.5) * p_ij;
joule_heat_edge_[index_i] = power_i / vol_i;
```

Issue:

\[
C_{ij}
\]

From:

```cpp
pairwiseNegativeLaplaceWeight(...)
```

Old pairwise Laplace discrete coefficient. Usable for φ residual, not automatically physical conductance for Joule heat.

More specifically:

```text
C_ij can make φ equation self-consistent;
C_ij used in C_ij e_ij² as physical Joule heat needs separate validation;
current P_edge/P_recon≈1e5 shows that validation fails.
```

Hence H production shrinks coil current heavily:

```text
ampere_turns_eff = 1.93578
max_BSrc = 3.7e-6 T
max_EImag = 0.9745 V/m
P_particle ≈ 0.43 W
```

Old baseline needed:

```text
ampere_turns_eff ≈ 660
max_BSrc ≈ 1.27e-3 T
P_raw ≈ 50000 W
```

Graph energy formula overestimates power; 50 kW calibration shrinks coil current ~300×.

---

## 4. Why H Production `production_literature_passed=0` Is Correct

H log:

```text
edge_flux_literature_acceptance:
edge_res_gate=1
P_edge_finite=1
P_edge_over_recon_soft=0
phi_res=1
P_raw=1
fields=1
literature_passed=1
production_literature_passed=0
```

Reasonable judgment.

Although:

```text
edge residual is good;
P_total_edge calibrated to 50 kW;
```

But:

```text
P_edge_over_recon_soft = 0
```

Graph edge power and reconstructed physical power are inconsistent—not production.

Recommend clearer log naming:

```text
literature_passed=1
```

Easily misread as edge-flux production passed. Add or replace with:

```text
edge_stage1_residual_passed=1
edge_power_closure_passed=0
production_literature_passed=0
```

Correct H case status:

```text
edge-flux residual solver: PASS
edge-flux physical power closure: FAIL
edge-flux production: FAIL
```

---

## 5. Fix Direction: Change Stage 1 Primary Power Definition

### 5.1 No Longer Use `P_total_edge` for 50 kW Primary Calibration

`P_total_edge` is graph energy diagnostic, not physical Joule power. Do not use it for French 50 kW calibration short term.

Stage 1 primary calibration target:

\[
P_{\mathrm{recon}}
=
\sum_i \frac12 \sigma_i |E_i^{\mathrm{edge}}|^2V_i.
\]

Where:

```text
E_i^edge from edge-drop least-squares reconstruction;
J_i^edge = sigma_i E_i^edge;
Q_i^edge-recon = 0.5 * sigma_i * |E_i^edge|².
```

H case:

```text
P_total_recon = 0.429548 W
P_particle    = 0.429642 W
```

Edge reconstruction matches old particle-field energy scale—reconstruction units are more credible.

---

### 5.2 New Stage 1 Calibration Strategy

When:

```text
--ophelie-current-form=edge-flux
```

Use:

```text
P_total_recon
```

As `joule_power_raw` / calibration power, not `P_total_edge`.

If code has:

```cpp
if (ophelieUseEdgeFluxElectromotiveRhs(params) && result.joule_power_edge > 0.0)
{
    result.joule_power_raw = result.joule_power_edge;
}
```

Change to:

```cpp
if (ophelieUseEdgeFluxElectromotiveRhs(params) && result.edge_power.p_total_recon > 0.0)
{
    result.joule_power_raw = result.edge_power.p_total_recon;
}
```

Also keep clear fields:

```cpp
result.joule_power_graph_edge = result.edge_power.p_total_edge;
result.joule_power_recon_edge = result.edge_power.p_total_recon;
```

Avoid ambiguous `joule_power_edge`.

---

### 5.3 Expected Results After Fix

Rerun H edge-flux production:

```text
P_recon ≈ 50000 W
ampere_turns_eff ≈ hundreds of A
max_BSrc ≈ 1e-3 T
P_graph_edge / P_recon may still be large, warning only
edge_res_red > 100
production_literature_passed may become 1
```

If after fix:

```text
P_recon still small;
ampere_turns_eff still only a few A;
```

Calibration still uses graph power—not fully fixed.

---

## 6. Keep `P_total_edge` but Clarify Role

Preserve, but mark diagnostic only.

Rename from:

```text
P_total_edge
```

To:

```text
P_graph_edge
```

Or log as:

```text
P_total_edge_graph_diagnostic
```

Avoid mistaking for physical power.

Current:

```text
P_edge/P_recon≈1.16e5
```

Record in dev log:

```text
pairwise Laplace coefficient C_ij cannot be used directly as physical conductance for Joule heat.
```

Future graph energy normalization:

\[
P_{\mathrm{graph,norm}}
=\kappa P_{\mathrm{graph}},
\]

\(\kappa\) from uniform-field MMS:

\[
P_{\mathrm{exact}}
=\frac12\sigma |E_0|^2V.
\]

Later research, not current production mainline.

---

## 7. Tests to Add or Fix

### 7.1 Direct Edge-Drop Sign Test

Current `test_3d_ophelie_edge_flux_sign` focuses on residual—insufficient. Residual is neighbor sum; pair-level sign issues may cancel.

Add direct edge-drop statistics:

```text
edge_drop_l2
edge_drop_linf
edge_drop_mean_abs
```

Setup:

\[
A=A_0,
\]

\[
\phi=-\omega A_0\cdot x.
\]

Each pair:

\[
e_{ij}
=(\phi_j-\phi_i)+\omega A_0\cdot(x_j-x_i)
\approx 0.
\]

Acceptance:

```text
edge_drop_linf < 1e-6  or reasonable small amount under float/SYCL
edge_drop_l2   < 1e-8~1e-6
```

Note: do not rely on residual alone.

---

### 7.2 Uniform E Power Test

Add:

```text
test_3d_ophelie_edge_flux_power_uniform_field
```

#### Case 1: Pure Potential Electric Field

Set:

\[
A=0,
\]

\[
\phi=-E_0\cdot x.
\]

Exact:

\[
E=E_0,
\]

\[
P_{\mathrm{exact}}
=\frac12\sigma |E_0|^2 V.
\]

Check:

```text
P_recon / P_exact
P_graph / P_exact
```

Expected:

```text
P_recon/P_exact ≈ 1;
P_graph/P_exact may be very large, diagnostic only.
```

#### Case 2: Pure Vector-Potential Induction Field

Set:

\[
\phi=0,
\]

\[
A=A_0.
\]

Then:

\[
E=-\omega A_0,
\]

\[
P_{\mathrm{exact}}
=\frac12\sigma\omega^2 |A_0|^2 V.
\]

Compare:

```text
P_recon / P_exact
P_graph / P_exact
```

Cleaner than French case for edge reconstruction physical power.

---

### 7.3 I×2 Scaling Test: Use `P_recon`

Old acceptance:

```text
I×2 → P_edge×4
```

Now:

```text
I×2 → P_recon×4
```

`P_graph` can still be logged, diagnostic only.

---

## 8. Concrete Code Change Plan

### Step 1: Preserve Edge Residual; Do Not Change φ Solve

These classes are correct—do not change for power issue:

```cpp
ComputeOphelieEdgeFluxPhiRhsFromASrcCK
ComputeOphelieEdgeFluxResidualCK
```

Edge residual already drops significantly.

---

### Step 2: Fix Power Ownership

In `runFrenchReducedEmPipeline()` or French pipeline, find raw power logic when current form is edge-flux.

Old logic likely:

```cpp
if (ophelieUseEdgeFluxElectromotiveRhs(params) && result.joule_power_edge > 0.0)
{
    result.joule_power_raw = result.joule_power_edge;
}
```

Change to:

```cpp
if (ophelieUseEdgeFluxElectromotiveRhs(params) && result.edge_power.p_total_recon > 0.0)
{
    result.joule_power_raw = result.edge_power.p_total_recon;
}
```

Add clearer fields:

```cpp
result.joule_power_graph_edge = result.edge_power.p_total_edge;
result.joule_power_recon_edge = result.edge_power.p_total_recon;
```

Deprecate or rename `joule_power_edge` if present.

---

### Step 3: Fix Acceptance

Current acceptance may resemble:

```cpp
acceptance.p_total_edge = joule_power_edge;
acceptance.power_edge_finite = std::isfinite(joule_power_edge) && joule_power_edge > 0.0;
acceptance.power_raw_ok = abs(joule_power_edge - target) / target < tol;
acceptance.production_literature_passed = acceptance.passed && acceptance.power_edge_over_recon_soft_ok;
```

Recommend:

```cpp
acceptance.p_total_recon = power_metrics.p_total_recon;
acceptance.p_graph_edge = power_metrics.p_total_edge;

acceptance.power_recon_finite =
    std::isfinite(power_metrics.p_total_recon) && power_metrics.p_total_recon > 0.0;

acceptance.power_raw_ok =
    abs(power_metrics.p_total_recon - target) / target < tol;

acceptance.graph_power_diagnostic_ok =
    std::isfinite(power_metrics.p_total_edge) && power_metrics.p_total_edge > 0.0;
```

Do not use:

```text
P_edge/P_recon ∈ [0.5,2]
```

As production gate.

Use warning:

```cpp
acceptance.graph_power_matches_recon_warning =
    ratio >= 0.5 && ratio <= 2.0;
```

Log:

```text
P_graph_over_recon_warning=0
```

But do not fail `production_literature_passed` for that.

New production gate:

```text
edge_res_red > 100
P_recon ≈ 50 kW
P_recon finite
phi residual ok
fields finite
I scaling ok
```

---

### Step 4: Split `JouleHeatEdge` Field

Current:

```text
JouleHeatEdge
```

Misleading—it stores graph energy heat, not edge-reconstructed physical heat.

Split:

```text
JouleHeatEdgeGraph
JouleHeatEdgeRecon
```

Minimum short term:

```text
1. Rename or mark current JouleHeatEdge as graph diagnostic in log/VTP;
2. Add JouleHeatEdgeRecon;
3. VTP priority output JouleHeatEdgeRecon.
```

`JouleHeatEdgeRecon`:

\[
Q_i^{\mathrm{edge-recon}}
=\frac12\sigma_i |E_i^{\mathrm{edge}}|^2.
\]

If reconstruction already has:

```text
EEdgeRecon
JEdgeRecon
```

In reconstruction CK:

```cpp
joule_heat_edge_recon_[i] = 0.5 * sigma_i * e_edge_i.squaredNorm();
power_edge_recon_particle_[i] = joule_heat_edge_recon_[i] * Vol_i;
```

Or separate CK:

```cpp
ComputeOphelieEdgeFluxReconstructedJouleHeatCK
```

---

### Step 5: Rerun H Production

Command similar:

```bash
test_3d_ophelie_french_reduced \
  --reload=1 \
  --literature-mode \
  --state_recording=0 \
  --ophelie-current-form=edge-flux \
  --phi-edge-flux-diagnostics=1
```

After fix expected:

```text
P_recon ≈ 50000 W
ampere_turns_eff ≈ hundreds of A
max_BSrc ≈ 1e-3 T
P_graph_over_recon still large, warning only
edge_res_red > 100
production_literature_passed may become 1
```

If still:

```text
P_recon ≪ 50 kW
ampere_turns_eff ≈ 1~2 A
```

Calibration still uses graph power.

---

## 9. Code Structure Review and Recommendations

### 9.1 What Works Well

```text
1. electromagnetic_ophelie_edge_flux.h/.hpp centralizes edge-flux mainline;
2. core operators use CK / Inner interaction, not external graph;
3. --ophelie-current-form=edge-flux is main switch;
4. --phi-projection-operator=edge-flux vs current form already distinguished in logs;
5. archive direction correct; archive_phi_divgrad_residual_floor/README.md useful.
```

Preserve these.

---

### 9.2 Code Structure Issues to Fix

#### 9.2.1 Diagnostics Should Not Include CLI

`electromagnetic_ophelie_edge_flux_diagnostics.h` may still have:

```cpp
#include "electromagnetic_ophelie_cli.h"
```

Unclean. Diagnostics should not depend on CLI.

Move from CLI:

```cpp
phiProjectionOperatorKindName(...)
inferOpheliePhiProjectionOperatorKind(...)
ophelieCurrentFormKindName(...)
```

To:

```text
electromagnetic_ophelie_parameters.h
```

Or new lightweight header:

```text
electromagnetic_ophelie_projection_operator.h
```

CLI and diagnostics both include base header.

---

#### 9.2.2 Legacy Diagnostic Should Not Reuse `DivJImag`

Legacy edge flux residual may reuse:

```text
names.div_j_imag
```

As temporary residual container.

OK as historical diagnostic temporarily; structurally unsafe. Add:

```text
LegacyEdgeFluxResidual
EdgeFluxResidualImag
```

Avoid overwriting real `DivJImag`.

---

#### 9.2.3 `P_total_edge` Naming Must Change

Misleading name. Change to:

```text
P_graph_edge
P_recon_edge
```

Or log explicitly:

```text
P_total_edge_graph_diagnostic
P_total_edge_reconstructed_physical
```

---

#### 9.2.4 Sign Test Should Be Direct Edge-Drop Test

If sign test only checks residual, add pair-level edge drop stats:

```text
edge_drop_l2
edge_drop_linf
edge_drop_mean_abs
```

Avoid false pass from neighbor-sum cancellation.

---

## 10. Current Edge-Flux Status Summary

### 10.1 Completed

```text
- SPH/SYCL edge-drop RHS;
- edge residual CK;
- edge residual reduction > 10000;
- edge E/J least-squares reconstruction;
- edge graph JouleHeat diagnostic;
- --ophelie-current-form=edge-flux main switch;
- French pipeline wired;
- directory archive and README started.
```

### 10.2 Not Passed

```text
- graph edge power and reconstructed physical power differ by 1e5;
- 50 kW calibration currently uses wrong power source;
- production_literature_passed=0.
```

### 10.3 Current Bottleneck

```text
Not φ solver;
Not edge residual;
Not SPH/SYCL structure;
But JouleHeat / power physical definition and calibration target.
```

---

## 11. Direct Development Instructions for Cursor

Copy-paste:

```text
Continue edge-flux direction but fix Stage 1 primary power definition.

Conclusion:
- edge residual / phi solve succeeds;
- P_total_edge = Σ 0.25*C_ij*edge_drop² cannot be physical Joule power;
- P_edge/P_recon≈1.16e5 shows pairwise Laplace weight C_ij is not direct physical conductance for Joule heat;
- 50 kW calibration must not use P_total_edge;
- Stage 1 primary power switches to P_recon from edge-reconstructed E/J.

Modify in order:

1. Preserve ComputeOphelieEdgeFluxPhiRhsFromASrcCK and ComputeOphelieEdgeFluxResidualCK; do not change φ solve.

2. For current_form=edge-flux, change calibration power from p_total_edge to p_total_recon.

3. Rename or mark P_total_edge as P_graph_edge_diagnostic; no longer production gate.

4. Add JouleHeatEdgeRecon:
   Q_i = 0.5 * sigma_i * |EEdgeRecon_i|²
   P_recon = Σ Q_i V_i
   VTP output JouleHeatEdgeRecon.

5. Acceptance:
   edge_res_red > 100
   P_recon≈50kW
   P_recon finite
   phi residual ok
   fields finite
   P_graph/P_recon warning only, not production fail.

6. Add direct edge-drop sign test:
   A=A0, phi=-omega*A0·x
   each edge_drop ≈ 0, not residual only.

7. Add uniform field power test:
   A=0, phi=-E0·x
   exact P=0.5*sigma*|E0|²*V
   validate P_recon/P_exact≈1;
   P_graph/P_exact diagnostic only.

8. Rerun H edge-flux production.
   Expect ampere_turns_eff back to hundreds of A;
   P_recon≈50kW;
   P_graph/P_recon may still be large but non-blocking.
```

---

## 12. Next-Step Order

Recommended next round:

```text
Step 1: Fix power ownership
  - calibration target → P_recon;
  - P_graph_edge diagnostic only;
  - adjust acceptance gate.

Step 2: Add JouleHeatEdgeRecon field
  - VTP output reconstructed Joule heat;
  - separate graph heat and recon heat.

Step 3: Add direct edge-drop sign test
  - A=A0, phi=-omega*A0·x;
  - pair-level edge_drop ≈ 0.

Step 4: Add uniform field power test
  - A=0, phi=-E0·x;
  - phi=0, A=A0;
  - validate P_recon/P_exact.

Step 5: Rerun H/A/G comparison
  - H: full edge-flux production;
  - A: div-grad fallback;
  - G: phi-only diagnostic.

Step 6: If H passes, consider Stage 2 A_ind diagnostic
  - do not connect thermal coupling before Stage 1 fix;
  - do not treat P_graph_edge as physical quantity.
```

---

## 13. Final Conclusion

Edge-flux solver main direction is correct:

```text
SPH/SYCL structure correct;
edge-drop RHS wired correctly;
edge residual excellent;
φ equation residual already ~1e-4;
```

Failure is not the φ solver but:

```text
treating graph edge energy as physical Joule power.
```

Do not reject edge-flux because `P_edge/P_recon≈1e5`. Immediately switch Stage 1 primary power to edge-reconstructed E/J Joule heat, i.e. `P_recon`, and demote `P_graph_edge` to diagnostic.

If after fix:

```text
P_recon≈50kW;
edge_res_red>100;
JouleHeatEdgeRecon distribution reasonable;
I×2→P_recon×4;
```

edge-flux can continue as OPHELIE-like solver mainline. Next stage:

```text
A_ind = K[J_edge_reconstructed]
Picard self-induction
TEAM7 benchmark
JouleHeat → thermal one-way
```
