# Stage 10 Contact Latest Progress, Phi Gauge Penalty, and A Divergence Constraint Diagnostic Execution Plan

This document is for Cursor execution. It summarizes the latest Stage 10 Contact multi-body A--phi development results and adds a new key question:

> Has the current Laplace-form A--phi solver implemented Coulomb gauge correction for `∇·A = 0`?

Conclusions:

1. Contact operators and coupled multi-body GMRES mainline are largely connected;
2. The main root cause of two-body MMS left-body field error has shifted from "Contact/GMRES error" to "insufficient phi gauge penalty calibration";
3. `phi_gauge_penalty` fixes the gauge of `phi`, not the divergence constraint on `A`;
4. Current code has not implemented projection, penalty, or explicit constraint for `∇·A = 0`;
5. Next phase should do `Stage 10.6 divA diagnostic` first, then decide whether to add A-divergence penalty or projection;
6. Defer cold-crucible four-body and 10D thermal coupling; close gauge/penalty/three-body TEAM7 parameter sensitivity first.

---

## 0. Current One-Sentence Status

Contact is no longer at the early inner-comparison stage. It has entered:

```text
multi-body Contact engineering scaffold validation
```

But still cannot say:

```text
real TEAM7 benchmark is complete;
Coulomb gauge ∇·A=0 is satisfied;
can directly enter cold-crucible thermal coupling.
```

More accurate status:

```text
Contact operator layer: largely complete;
coupled multi-body GMRES: usable as current mainline;
three-body TEAM7-like scaffold: engineering passed;
phi gauge penalty: root cause largely confirmed, needs calibration;
A divergence-free gauge: not implemented, needs immediate diagnostic;
real TEAM7 benchmark: later stage.
```

---

# Part I. Latest Contact / Penalty Progress Summary

## 1. Completed and Largely Trustworthy Content

### 1.1 Contact Operator Layer

The following Contact operator layer can be considered largely connected:

```cpp
AphiPairwiseLaplaceCK<Contact<>>
AphiGradPhiCouplingCK<Contact<>>
AphiDivSigmaACouplingCK<Contact<>>
AphiComputeScalarPhiGradientCK<Contact<>>
AphiApplyCK<Contact<>>
AphiComputeBlockJacobiDiagonalCK<Contact<>>
```

Passed validation types include:

```text
pairwise Laplace equivalence
apply vs monolithic
fused apply equivalence
block Jacobi diagonal equivalence
interface RHS consistency
three-body apply / observable vs monolithic
```

Conclusion:

```text
The main issue is no longer whether Contact pairwise operators work,
but gauge constraints, penalty calibration, bodywise residual, and real benchmark formalization.
```

---

### 1.2 Workspace Contamination Has Been Fixed

Previously found:

```text
If single-body Contact GMRES Arnoldi uses full Contact apply,
it reads neighbor body Krylov workspace, e.g. GMRESZk,
causing current body Krylov operator contamination.
```

Mathematically, single-body block-GS GMRES should use:

```text
A_ii z_i
```

Not:

```text
A_ii z_i + A_ij z_j
```

Therefore the current fix direction is correct:

```text
single-body Contact GMRES:
    Arnoldi apply uses block-diagonal Contact apply

coupled multi-body GMRES:
    Arnoldi apply uses full Contact apply
```

Current mainline should be explicit:

```text
Production / benchmark mainline:
    AphiMultiBodyContactGMRESSolverCK
    coupled multi-body GMRES
    full Inner + full Contact apply
    global Arnoldi / global dot / global norm / global update

Diagnostic / fallback:
    AphiGMRESContactSolverCK
    single-body Contact GMRES
    block-diagonal Arnoldi
    optional polish / block-GS diagnostic
```

Do not revert to legacy block-GS as default production path.

---

### 1.3 Three-Body TEAM7-Like Scaffold Passed Engineering Prototype

Current three-body setup:

```text
air body
coil body
plate/conductor body
```

Current results:

```text
mono_plate_joule_power      = 0.136195
coupled_plate_joule_power   = 0.136228
plate_joule_power_gap       ≈ 2.44e-4
plate_j_L2_gap              ≈ 9.58e-5
plate_joule_max_gap         ≈ 3.87e-4

coupled_outer               = 9
coupled_converged           = 1
coupled_global_true_rel     ≈ 3.06e-4
coupled_max_bodywise_rel    ≈ 1.51e-4
air_true_rel                ≈ 1.51e-4
coil_true_rel               ≈ 1.06e-4
plate_true_rel              ≈ 3.97e-5
```

Judgment:

```text
Three-body Contact scaffold E/J/Joule observables align well with monolithic scaffold.
This shows Contact splitting and coupled GMRES mainline did not break plate Joule / current observables.
```

But note:

```text
This is still TEAM7-like impressed-current scaffold,
not formal TEAM7 benchmark.
```

---

## 2. Latest Conclusions on Phi Gauge Penalty

The most critical new evidence this round comes from `phi_gauge_penalty_sweep`.

### 2.1 Penalty Sweep Result Summary

Representative results:

```text
lambda_phi = 0:
    converged        = 1
    outer            = 3
    global_true_rel  ≈ 1.79e-6
    left_true_rel    ≈ 1.57e-6
    left_continuous  ≈ 88.36
    phi_real_mean    ≈ 0.829676

lambda_phi = 1:
    converged        = 1
    outer            = 5
    global_true_rel  ≈ 3.90e-6
    left_true_rel    ≈ 6.15e-6
    left_continuous  ≈ 0.0523

lambda_phi = 10:
    converged        = 1
    outer            = 2
    global_true_rel  ≈ 4.99e-6
    left_true_rel    ≈ 4.71e-5
    left_continuous  ≈ 0.0205

lambda_phi = 100:
    converged        = 1
    outer            = 2
    global_true_rel  ≈ 2.09e-7
    left_true_rel    ≈ 6.09e-6
    left_continuous  ≈ 0.00265
```

Conclusion:

```text
At lambda_phi = 0, Krylov residual is good but phi/raw field drifts severely;
At lambda_phi = 100, left continuous error drops from O(1%) to about 0.26%;
This suggests the previous left-field issue was likely insufficient phi gauge regularization.
```

---

### 2.2 Meaning of Phi-Only Test

In phi-only case:

```text
lambda_phi = 10:
    left_continuous ≈ 4.6e-5

lambda_phi = 0:
    left_continuous ≈ 1.12
```

This shows:

```text
phi equation itself has obvious constant/null gauge drift;
phi_gauge_penalty is effective for locking phi constant DOF.
```

Also:

```text
global Krylov convergence does not automatically guarantee comparable raw phi field.
```

---

## 3. Can Stage 10.5 Be Closed?

Suggested status change:

```text
Stage 10.5:
    Root cause largely identified as phi gauge regularization.
    Status: partially resolved.
    Remaining: penalty calibration and observable sensitivity.
```

Suspicions that can be closed:

```text
Contact pairwise operator fundamentally wrong;
coupled multi-body GMRES mainline wrong;
interface band residual out of control;
RHS assembly fundamentally wrong;
polish sweep insufficient;
monolithic-vs-split is the only cause.
```

Still cannot close:

```text
whether lambda_phi = 100 should be default;
whether lambda_phi should scale with dp or Laplace diagonal;
whether lambda_phi affects three-body TEAM7 E/J/Joule;
whether A-side gauge / divA has hidden issues;
whether field error should use gauge-aligned metric.
```

---

## 4. Should Default `lambda_phi` Be Changed to 100?

Recommendation:

```text
Short term: can change MMS/scaffold default test parameters to lambda_phi = 100.
Long term: do not treat fixed 100 as final physical default for all problems.
```

Reason:

```text
At current dp ≈ 0.1, 1 / dp^2 = 100,
so lambda_phi = 100 likely matches discrete Laplace scale.
```

More stable long-term parameterization:

```text
lambda_phi = C_phi / dp^2
```

Or:

```text
lambda_phi = C_phi * median(|LaplacePhiDiagonal|)
```

Recommended Cursor action for now:

```text
A/B:
    lambda_phi = 10 vs 100

sweep:
    lambda_phi = 10, 30, 100, 300
```

Then decide default from three-body E/J/Joule sensitivity.

---

## 5. Acceptance Logic: Tier-1 / Tier-2 Must Be Separated

### 5.1 Tier-1: Solver/Operator/Observable Acceptance

These can enter main `passed`:

```text
coupled_converged = 1
global_true_rel < tolerance or engineering tolerance
max_bodywise_true_rel < bodywise tolerance
interface_band_rel reasonable
discrete_mms_defect absolute value reasonable
E/J/Joule observable aligned with monolithic/scaffold reference
```

### 5.2 Tier-2: Field/Gauge Diagnostic

These should not directly gate mainline `passed`:

```text
left_continuous raw max-norm vs exact
raw phi mean error
raw A/phi block Linf
contact_vs_exact raw field error
contact_vs_mono raw field error
```

Reason:

```text
raw phi contains gauge DOF;
lambda_phi = 0 already proves: Krylov residual can be good while raw field drifts severely.
```

Recommended output:

```text
tier1_solver_passed
tier1_observable_passed
tier2_field_diagnostic_passed
field_diagnostic_warn
passed
```

Where:

```text
passed = tier1_solver_passed && tier1_observable_passed
```

Tier-2 temporarily as warning / diagnostic.

---

## 6. Positioning of Polish

Current polish sweep:

```text
sweep = 0:
    global_true_rel ≈ 4.99e-6
    left_true_rel   ≈ 4.71e-5
    left_continuous ≈ 2.05%

sweep = 1:
    global_true_rel ≈ 1.81e-6
    left_true_rel   ≈ 2.47e-5
    left_continuous ≈ 0.95%

sweep = 2:
    global_true_rel ≈ 1.04e-6
    left_true_rel   ≈ 1.41e-5
    left_continuous ≈ 0.96%

sweep = 3:
    global_true_rel ≈ 6.75e-7
    left_true_rel   ≈ 9.20e-6
    left_continuous ≈ 0.97%
```

Conclusion:

```text
polish improves global/bodywise residual,
but left_continuous plateaus around ~1%;
therefore polish is not the fundamental fix for left field issue.
```

Recommendation:

```text
production default:
    coupled GMRES only

diagnostic optional:
    --enable_contact_polish
    --polish_sweeps=N

Do not let polish mask phi gauge penalty calibration issues.
```

---

# Part II. Key New Question: Is `∇·A = 0` Correction Implemented?

## 7. Background: Why Care About `∇·A = 0`

Current approach uses Laplace-form A--phi equations, not direct curl-curl form.

Classic relation:

```text
∇ × ∇ × A = ∇(∇·A) - ∇²A
```

With Coulomb gauge:

```text
∇·A = 0
```

Then:

```text
∇ × ∇ × A = -∇²A
```

Therefore Laplace-form A equation implicitly requires handling gauge, especially:

```text
∇·A = 0
```

Otherwise vector Laplace form is not strictly equivalent to curl-curl form.

Early discussion mentioned incompressible-flow-like pressure Poisson correction:

```text
A* -> solve chi Poisson -> A = A* - grad chi
```

So that:

```text
∇·A = 0
```

This is Coulomb gauge projection.

---

## 8. What Already Exists in Current Code

### 8.1 `div_sigma_a_coupling`

Current code has something like:

```cpp
AphiDivSigmaACouplingCK
```

Adding to phi equation:

```text
j omega div(sigma A)
```

Derived from current continuity:

```text
div J = 0
J = sigma E
E = -j omega A - grad phi
```

So it controls:

```text
div[ sigma(-j omega A - grad phi) ] = 0
```

Or in discrete form:

```text
div(sigma A) + div(sigma grad phi) terms
```

But it is not:

```text
div A = 0
```

Especially in:

```text
spatially varying sigma
air / coil / plate / melt / crucible multi-body
conductive / non-conductive interface
```

`div(sigma A)` and `div A` are completely different.

---

### 8.2 `phi_gauge_penalty`

Current code has:

```cpp
AphiPhiGaugePenaltyCK
```

Form:

```text
L_phi += lambda_phi * phi
```

This term fixes:

```text
phi constant gauge / null mode
```

Not:

```text
A divergence-free gauge
```

Results show it is very important, but it cannot replace `∇·A = 0`.

---

### 8.3 `AphiVectorGradientDivergenceCK`

Current code has something like:

```cpp
AphiVectorGradientDivergenceCK
```

Function:

```text
divergence = trace(grad A)
```

i.e. compute:

```text
div A
```

But from current design, it mainly appears in diagnostic / debug assembly paths, e.g.:

```text
div_a_real_diagnostic
div_a_imag_diagnostic
AphiDivSigmaAConstSigmaDiagnosticCouplingCK
```

It is not currently a production solver constraint term.

---

## 9. What Is Not Implemented

Currently no mechanism seen for:

```text
A <- A - grad chi projection
solve Laplace chi = div A
phi <- phi + j omega chi
```

Also no:

```text
A equation LHS += lambda_A * grad(div A)
```

Also no:

```text
A divergence penalty
Coulomb gauge projection
A gauge correction step
A-side nullspace removal
```

Clear conclusion:

```text
Current code implements phi gauge penalty,
but not A divergence-free gauge constraint.
```

---

## 10. Will This Affect Current Results?

### 10.1 For Contact Scaffold

Current Contact vs monolithic three-body TEAM7 observables are already good, showing:

```text
Contact splitting did not break current Laplace-form scaffold E/J/Joule observables.
```

So this does not overturn current Contact work.

### 10.2 For Real TEAM7 / Cold Crucible

If seriously benchmarking real TEAM7 or cold-crucible physics later, must answer:

```text
Does current solution satisfy or approximately satisfy div A = 0?
If not, are B/E/J/Joule still reliable?
Is A-gauge penalty or projection needed?
```

Otherwise risks:

```text
A contains non-physical longitudinal component;
A point values do not match FEM reference;
B = curl A may be OK, but A/phi raw field interpretation unreliable;
E/J/Joule may be polluted in some high-contrast/complex geometries;
cold-crucible melt/crucible Joule may mix gauge mode effects.
```

---

# Part III. New Stage 10.6: A Gauge / divA Diagnostic

## 11. Stage 10.6 Overall Goal

Without changing current main equations, first diagnose whether current solution is close to:

```text
div A = 0
```

And judge whether to add later:

```text
A-divergence penalty
```

Or:

```text
Coulomb gauge projection
```

---

## 12. Task 10.6-A: Organize divA Diagnostic Helper

### 12.1 Goal

Organize existing or similar:

```cpp
AphiVectorGradientDivergenceCK
```

Into reusable diagnostic helper.

Suggested naming:

```cpp
AphiComputeDivACK
```

Or test helper:

```cpp
AphiDivADiagnosticHelper
```

### 12.2 Input

Compute divergence for:

```text
A_real = (A_x_real, A_y_real, A_z_real)
A_imag = (A_x_imag, A_y_imag, A_z_imag)
```

### 12.3 Output Variables

Suggested registration:

```text
DivAReal
DivAImag
DivAAbs
DivARelative
```

If not adding particle variables, at least output reduction metrics.

### 12.4 Reduction Metrics

Global:

```text
div_A_real_L2
div_A_imag_L2
div_A_L2
div_A_Linf
grad_A_L2
div_A_relative = div_A_L2 / grad_A_L2
```

Bodywise:

```text
body_div_A_L2
body_div_A_Linf
body_grad_A_L2
body_div_A_relative
```

Regionwise:

```text
left_div_A_relative
right_div_A_relative
air_div_A_relative
coil_div_A_relative
plate_div_A_relative
interface_band_div_A_relative
core_div_A_relative
```

---

## 13. Task 10.6-B: Connect divA to Existing Tests

At minimum connect to:

```text
test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured
test_3d_aphi_ck_contact_phi_gauge_penalty_sweep_diagnostic
test_3d_aphi_ck_gmres_contact_team7_three_body_scaffold
```

Also suggested:

```text
test_3d_aphi_ck_contact_left_field_error_diagnostic
```

---

## 14. Task 10.6-C: Relationship Between divA and Penalty

In `phi_gauge_penalty_sweep` add output:

```text
lambda_phi
global_true_rel
left_continuous
phi_mean_error
div_A_relative
left_div_A_relative
right_div_A_relative
```

Focus on:

```text
lambda_phi = 0, 1, 10, 100
```

Judgment:

```text
If lambda_phi improves left_continuous but divA unchanged:
    phi gauge is independent; A gauge still open.

If lambda_phi also lowers divA:
    A/phi coupling may indirectly stabilize A divergence, but still record independently.

If divA always small:
    temporarily no A penalty needed; record evidence only.

If divA large:
    enter Stage 10.7 A-divergence penalty.
```

---

## 15. Task 10.6-D: Relationship Between divA and E/J/Joule

In three-body TEAM7 output:

```text
plate_div_A_relative
plate_E_L2
plate_J_L2
plate_Joule_power
plate_Joule_L2
```

Relationship judgment:

```text
If divA large but E/J/Joule stable vs monolithic:
    current observable insensitive to gauge pollution, but real benchmark still needs statement.

If divA large and E/J/Joule also change significantly with penalty/source/sigma:
    must consider A-divergence penalty or projection.

If divA small and E/J/Joule stable:
    current Laplace-form scaffold can continue.
```

---

## 16. Preliminary divA Threshold Suggestions

Do not set hard pass/fail yet. Suggested diagnostic levels:

```text
div_A_relative < 1e-2:
    good; no A correction needed for now.

1e-2 <= div_A_relative < 1e-1:
    needs attention; keep warning; continue parameter sweep.

div_A_relative >= 1e-1:
    high risk; enter A-divergence penalty or projection design.
```

These thresholds are engineering initial values, not final theoretical standards.

---

# Part IV. If divA Has Issues, How to Fix

## 17. Option A: Diagnostic Only

Applicable when:

```text
div_A_relative is small;
B/E/J/Joule observables stable;
real TEAM7 reference comparison shows no gauge issue.
```

Advantages:

```text
no change to main equations;
does not break current passed tests;
lowest implementation cost.
```

---

## 18. Option B: A-Divergence Penalty

If divA is not small, prioritize:

```text
A equation LHS += lambda_A * grad(div A)
```

i.e. penalty on longitudinal mode.

### 18.1 New Parameters

```cpp
bool use_a_divergence_penalty;
Real a_divergence_penalty;
```

Suggested placement:

```cpp
AphiLhsTermFlags
AphiCouplingParameters
```

Or same parameter structure managing `phi_gauge_penalty`.

### 18.2 New Operators

```cpp
AphiComputeDivACK<Inner<>>
AphiComputeDivACK<Contact<>>

AphiGradDivAPenaltyCK<Inner<>>
AphiGradDivAPenaltyCK<Contact<>>
```

Or fused into:

```cpp
AphiApplyCK
```

### 18.3 Notes

This term is not a physical material term; label as:

```text
Numerical Coulomb-gauge regularization
```

Must do:

```text
lambda_A sweep
divA reduction
E/J/Joule sensitivity
TEAM7 observable sensitivity
```

---

## 19. Option C: Coulomb Gauge Projection

More complete but complex:

```text
solve Laplace chi = div A
A <- A - grad chi
phi <- phi + j omega chi
```

Advantages:

```text
closest to incompressible-flow projection;
clear physical interpretation;
directly suppresses divA.
```

Disadvantages:

```text
extra Poisson solve needed;
multi-body Contact Poisson needed;
chi boundary conditions complex;
phi update must preserve E gauge invariance;
high implementation cost.
```

Recommendation:

```text
Do not do projection directly for now.
Do divA diagnostic first;
if divA issue severe, try A-divergence penalty first.
```

---

# Part V. Next Execution Plan

## 20. Immediate: Penalty A/B and Stage 10.6 divA

### Step 1: MMS Default Penalty A/B

Run:

```text
lambda_phi = 10
lambda_phi = 100
```

Tests:

```text
test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured
test_3d_aphi_ck_contact_left_field_error_diagnostic
test_3d_aphi_ck_contact_two_body_polish_sweep_diagnostic
```

Output:

```text
global_true_rel
max_bodywise_true_rel
interface_band_rel
discrete_defect
left_continuous
left_block_L2_after_mean_sub
left_E/J/Joule error if available
```

---

### Step 2: Three-Body TEAM7 Penalty Sweep

Run:

```text
lambda_phi = 10, 30, 100, 300
```

Test:

```text
test_3d_aphi_ck_gmres_contact_team7_three_body_scaffold
```

Output:

```text
coupled_global_true_rel
max_bodywise_true_rel
air_true_rel
coil_true_rel
plate_true_rel

plate_joule_power
plate_E_real_max
plate_J_L2
plate_Joule_L2
plate_joule_max
plate_solution_block_max
```

Judgment:

```text
If plate_Joule_power / J_L2 change < 0.5%--1% within 30--300:
    lambda_phi = 100 can be current scaffold default.

If change significant:
    lambda_phi needs scaling; cannot fix at 100.
```

---

### Step 3: Implement divA Diagnostic

Add:

```text
Stage 10.6 A-gauge diagnostic
```

Minimum integration:

```text
two-body MMS
phi penalty sweep
three-body TEAM7 scaffold
```

Output:

```text
div_A_real_L2
div_A_imag_L2
div_A_L2
div_A_Linf
grad_A_L2
div_A_relative
bodywise_div_A_relative
regionwise_div_A_relative
```

---

### Step 4: Update Documentation

Update:

```text
CURSOR_APHI_STAGE10_CONTACT_PROGRESS.md
CURSOR_APHI_STAGE10_5_PHI_GAUGE_DIAGNOSTIC_RECORD.md
```

Add:

```text
CURSOR_APHI_STAGE10_6_A_GAUGE_DIVERGENCE_DIAGNOSTIC_RECORD.md
```

Record must include:

```text
current has phi gauge penalty;
current has no A divergence-free correction;
divA diagnostic results;
whether to enter A-divergence penalty.
```

---

## 21. Deferred Content

Before completing penalty A/B + divA diagnostic, not recommended:

```text
cold-crucible four-body EM scaffold
10D Joule -> heat coupling
large helper merge
common GMRES driver large refactor
real TEAM7 geometry development
```

Reason:

```text
If A gauge / divA issue not diagnosed clearly,
subsequent four-body cold-crucible melt/crucible Joule interpretation will be unstable.
```

---

# Part VI. Cursor Execution Checklist

## 22. Code Task List

### Task A: Keep and Improve Phi Gauge Penalty

1. Confirm all paths include `lambda_phi * phi`:
   ```text
   fused apply
   debug apply / assembly
   block Jacobi diagonal
   coupled GMRES path
   contact path
   inner path
   ```
2. Add comment:
   ```text
   phi_gauge_penalty is numerical gauge regularization, not a physical term.
   ```
3. Support CLI or test option:
   ```text
   --phi_gauge_penalty
   --use_phi_gauge_penalty
   ```

### Task B: Penalty Sensitivity Tests

1. Two-body MMS:
   ```text
   lambda_phi = 10, 100
   ```
2. Three-body TEAM7:
   ```text
   lambda_phi = 10, 30, 100, 300
   ```
3. Output solver + field + E/J/Joule metrics.

### Task C: Implement divA Diagnostic

1. Organize or reuse:
   ```cpp
   AphiVectorGradientDivergenceCK
   ```
2. Add reduction helper:
   ```text
   div_A_L2
   div_A_Linf
   div_A_relative
   bodywise_div_A_relative
   regionwise_div_A_relative
   ```
3. Connect to:
   ```text
   two-body MMS
   penalty sweep
   three-body TEAM7
   left field diagnostic
   ```

### Task D: Add Stage 10.6 Documentation

Add:

```text
CURSOR_APHI_STAGE10_6_A_GAUGE_DIVERGENCE_DIAGNOSTIC_RECORD.md
```

Record:

```text
1. No A projection currently;
2. No A-divergence penalty currently;
3. Only phi gauge penalty currently;
4. divA diagnostic results;
5. whether A penalty needed.
```

### Task E: Decide Stage 10.7 Based on divA Results

If:

```text
div_A_relative < 1e-2
```

Then:

```text
record evidence; no A correction for now.
```

If:

```text
div_A_relative >= 1e-2
```

Then:

```text
keep warning; continue sweep.
```

If:

```text
div_A_relative >= 1e-1
```

Then enter:

```text
Stage 10.7 A-divergence penalty prototype
```

---

# Part VII. Current External Wording Recommendations

## 23. Can Say

```text
The Contact operator layer has been largely validated against monolithic inner scaffolds.
The coupled multi-body Contact GMRES path is operational.
Three-body TEAM7-like impressed-current scaffold agrees with the monolithic scaffold in plate E/J/Joule observables.
The previous left-field diagnostic issue is now largely attributed to insufficient phi-gauge regularization.
```

## 24. Cannot Say

```text
The real TEAM7 benchmark is fully validated.
The Coulomb gauge condition div A = 0 has been enforced.
The current solver is ready for thermal coupling.
The solution_block_max is magnetic field.
The phi gauge penalty also enforces div A = 0.
```

## 25. Most Accurate Wording

```text
Stage 10 Contact has reached the multi-body engineering-scaffold stage.
The remaining immediate issues are phi-gauge penalty calibration and A-gauge/divA diagnostics.
After these diagnostics, the solver can proceed to three-body parameter sweeps and then cold-crucible four-body EM scaffolds.
```

---

# 26. Final Decision

Current recommended route:

```text
Now:
    do not revert block-GS;
    do not enter 10D;
    do not merge helpers;
    do not rush cold-crucible four-body.

Do immediately:
    lambda_phi = 100 A/B;
    three-body lambda_phi sweep;
    Stage 10.6 divA diagnostic.

If divA small:
    fix current Laplace-form scaffold route;
    do three-body dp/sigma/source sweep;
    then cold-crucible four-body EM.

If divA large:
    do A-divergence penalty first;
    projection later.
```

Core conclusion:

> This penalty sweep proves the left-body raw field issue is mainly insufficient phi gauge regularization; but phi gauge penalty cannot replace A Coulomb gauge. Current code has no `∇·A=0` correction yet, so next step must add Stage 10.6 divA diagnostic, then decide whether to add A-divergence penalty or projection.
