# Stage 10 Contact Multi-Body A--phi: Current Status Review and Next-Step Execution Plan

This document guides Cursor in continuing development of the SPHinXsys SYCL/CK multi-body Contact A--phi solver. It focuses on answering four questions:

1. How far the current Contact work has progressed;
2. How much remains before Contact is "done";
3. Where it is blocked, whether it can be resolved, and how;
4. Whether current completion is relative to the inner scaffold, or already close to a real TEAM7 physics solver.

---

## 0. One-Sentence Conclusion

Current Contact development has passed the most critical operator-layer threshold: `Contact<>` Laplace / grad-phi / div-sigma-A / Joule / fused apply are largely usable, single-body GMRES workspace contamination has been diagnosed and fixed, and coupled multi-body GMRES can run two-body MMS and three-body TEAM7-like scaffold.

However:

> Current completion is mainly "multi-body Contact equivalence and engineering prototype relative to inner / monolithic scaffold," not "real TEAM7 physics benchmark fully complete."

More specifically:

- **Contact operator layer**: largely complete;
- **Contact coupled GMRES solver prototype**: usable, but still needs bodywise residual, field diagnostic, and more parameter sweeps;
- **Three-body TEAM7-like scaffold**: engineering prototype passes; plate Joule aligns with monolithic;
- **Real TEAM7 validation**: not complete; still missing real TEAM7 geometry / reference FEM/MFEM / literature benchmark / physical quantity comparison;
- **Cold-crucible Contact**: four-body real structure not started; only prior box/scaffold experience;
- **10D thermal coupling**: still not recommended immediately; should wait until Contact three-body/four-body EM is stable.

---

## 1. Tiered Assessment of Current Completion

### 1.1 Completion relative to inner / monolithic scaffold

If the goal is:

> Migrate the originally inner-only / monolithic region-tag A--phi operators to true multi-body Contact, and verify consistency with the monolithic scaffold.

Then current completion is already high.

Can be considered done:

1. `Contact<>` pairwise Laplace aligned with inner;
2. `Contact<>` grad-phi / div-sigma-A / Joule apply aligned with monolithic;
3. fused apply and block Jacobi diagonal pass equivalence tests;
4. workspace contamination confirmed;
5. single-body Contact GMRES Arnoldi changed to block-diagonal;
6. coupled multi-body GMRES can do full Contact apply;
7. three-body TEAM7-like scaffold plate Joule/E/J aligned with monolithic;
8. metric naming changed from `*_field_max` to `*_solution_block_max`, avoiding mislabeling A/phi solution block as magnetic field.

At the "compare inner / monolithic scaffold" level, current status:

```text
Contact operator layer: 80%--90% complete
Contact coupled-GMRES prototype: 60%--75% complete
three-body scaffold vs monolithic: preliminary pass
```

These percentages are not mathematically precise; they are engineering maturity judgments.

### 1.2 Completion relative to real TEAM7 problem

If the goal is:

> Be able to seriously claim the current A--phi Contact solver can quantitatively compute the real TEAM7 benchmark.

Then current status is insufficient.

Reasons:

1. Current three-body TEAM7 is still TEAM7-like scaffold, not real TEAM7 geometry;
2. Main comparison target is monolithic inner / region-tag reference, not external FEM/MFEM or literature benchmark;
3. Still need real geometry, coil/plate/air dimensions, material parameters, frequency, current, boundary conditions, probe-line data;
4. Still need to output and compare physical quantities:
   - magnetic vector potential A;
   - electric field E;
   - current density J;
   - Joule power density q;
   - total Joule power;
   - B = curl A, once implemented;
   - possibly H = nu B;
5. Still need dp sweep, sigma contrast sweep, frequency/source sweep;
6. Root cause of ~1% left-body field-level max-norm error in current two-body MMS not fully explained.

At the "real TEAM7 benchmark" level, current status:

```text
Real TEAM7 solver readiness: 35%--50%
```

This does not mean code cannot advance; it means current results still belong to:

> solver infrastructure and scaffold validation, not final physical benchmark validation.

---

## 2. Closed Items

### 2.1 Contact operator layer is largely closed

The following Contact operators can be considered Stage 10 Contact baseline capability:

```cpp
AphiPairwiseLaplaceCK<Contact<>>
AphiGradPhiCouplingCK<Contact<>>
AphiDivSigmaACouplingCK<Contact<>>
AphiComputeScalarPhiGradientCK<Contact<>>
AphiApplyCK<Contact<>>
AphiComputeBlockJacobiDiagonalCK<Contact<>>
```

These operators have passed:

- pairwise Laplace equivalence;
- apply vs monolithic;
- fused apply equivalence;
- block Jacobi diagonal equivalence;
- interface RHS consistency diagnostic;
- three-body apply vs monolithic.

Conclusion:

> The main problem is no longer the Contact pairwise operator itself, but multi-body solver semantics, convergence criteria, field-level MMS interpretation, and real benchmark formalization.

---

### 2.2 workspace contamination has been correctly identified

Previous issue:

```text
When single-body Contact GMRES Arnoldi uses full Contact apply,
it reads neighbor body GMRESZk / Krylov workspace.
```

Mathematically, single-body block-GS Krylov apply should only use:

```text
A_ii z_i
```

Not:

```text
A_ii z_i + A_ij z_j
```

Therefore changing single-body Contact GMRES Arnoldi to block-diagonal apply is correct.

Current route should be explicit:

```text
full Contact apply:
    for coupled GMRES;
    for residual;
    for RHS;
    for K(solution).

block-diagonal Contact apply:
    only for single-body GMRES / block-GS Arnoldi apply.
```

---

### 2.3 coupled multi-body GMRES should be the main line

Currently available:

```cpp
aphi_multibody_contact_gmres_ck.h/.hpp
```

For coupled multi-body GMRES.

This is the correct route. Future production / benchmark / customer-like simulation should use coupled multi-body GMRES as the main line, not legacy block-GS.

Recommended formal positioning:

```text
Production main path:
    AphiMultiBodyContactGMRESSolverCK
    coupled multi-body GMRES
    full Inner + Contact apply over all bodies
    global Arnoldi / global dot / global norm / global update

Diagnostic / fallback path:
    AphiGMRESContactSolverCK
    single-body Contact GMRES
    block-diagonal Arnoldi
    optional block-GS polish
```

---

### 2.4 Three-body TEAM7-like scaffold already has engineering significance

Current three-body setup:

```text
air body
coil body
plate/conductor body
```

Can run coupled GMRES, and plate Joule / E / J align with monolithic scaffold; plate Joule gap is about 0.02%.

This result is very important. It shows:

1. Multi-body Contact apply does not obviously fail in three-body structure;
2. coupled GMRES works across bodies;
3. plate observable is not broken by Contact splitting;
4. current code has basis to enter more realistic TEAM7-like / cold-crucible scaffold.

But it is still scaffold, not formal TEAM7 validation.

---

## 3. Open Issues

### 3.1 Two-body MMS left-body max-norm error ~1%

Current key data:

```text
coupled only:
    converged = 1
    outer = 2
    global_true_rel ≈ 5e-6
    left_continuous ≈ 2.0%
    right_continuous ≈ 1e-6

coupled + 1 polish:
    converged = 1
    outer = 2
    global_true_rel ≈ 1.8e-6
    left_continuous ≈ 0.96%
    right_continuous ≈ 8e-7

polish sweep 0→3:
    left_continuous ≈ 2% → 1% → 0.97% → 0.98%
```

Assessment:

> This does not look like simple insufficient iteration, because it plateaus after polish; nor like Contact operator completely wrong, because interface band / RHS consistency / three-body Joule all behave normally.

More likely causes:

1. gauge / nullspace / constant offset;
2. max-norm oversensitive to local boundary/interface particles;
3. manufactured exact not fully consistent with split-body Contact discrete operator;
4. left body is conductive/stiff region, more sensitive to field variables;
5. bodywise residual masked by global average;
6. body order / contact direction / material asymmetry causing field-level asymmetry.

Conclusion:

```text
Left max-norm 1% should not block three-body TEAM7 main line.
But it must be listed as Stage 10.5 open issue for separate resolution.
```

---

### 3.2 global residual insufficient to represent all bodies converged

Current coupled GMRES has:

```text
global_true_rel ≈ 1e-6
```

But left body:

```text
left_true_rel ≈ 2.5e-5
```

indicating global norm may mask bodywise issues.

Therefore must add or strengthen:

```text
bodywise_true_rel
bodywise_rhs_norm
bodywise_solution_norm
bodywise_residual_norm
max_bodywise_true_rel
```

Recommended output for all future multi-body solver tests:

```text
global_true_rel
max_bodywise_true_rel
bodywise_true_rel per body
bodywise_rhs_norm per body
bodywise_solution_norm per body
```

---

### 3.3 No real B/H output yet

E/J/Joule metric helpers exist, but without:

```text
B = curl(A)
H = nu B
```

do not claim magnetic field is already output.

Current `solution_block_max` is only:

```text
A_real
A_imag
phi_real
phi_imag
```

block norm/max. It is not B, not E, not J.

For magnetic field plots later, must add:

```cpp
AphiComputeMagneticFluxDensityCurlACK
AphiComputeMagneticFieldHCK
```

Output:

```text
B_real
B_imag
B_L2
B_max
H_real
H_imag
H_L2
H_max
```

---

### 3.4 Real TEAM7 still lacks external reference

Current three-body TEAM7-like scaffold aligned with monolithic only shows Contact splitting does not break the current discrete model.

It cannot prove:

```text
Current A--phi discretization is correct for real TEAM7 benchmark.
```

Real TEAM7 validation still needs:

1. real TEAM7 geometry;
2. mesh/particle resolution study;
3. external FEM/MFEM reference;
4. literature benchmark;
5. plate Joule / induced current / B-field profile comparison;
6. consistent frequency/source/current parameters;
7. consistent boundary conditions.

---

## 4. How Far Contact Is from "Done"

Contact "done" should be split into three levels.

### 4.1 Level 1: Contact infrastructure done

Definition:

```text
Contact operators, Contact apply, Contact coupled GMRES, three-body scaffold all run stably.
```

Current status:

```text
Near complete, about 75%--85%.
```

Remaining:

1. complete bodywise residual output;
2. mark left MMS 1% as diagnostic and do preliminary localization;
3. confirm three-body TEAM7 stable under different dp / sigma / source;
4. update progress documentation;
5. keep helpers unmerged but clearly categorized.

This level is the short-term target.

### 4.2 Level 2: Contact engineering solver ready

Definition:

```text
Usable as engineering solver for subsequent cold-crucible four-body scaffold and 10D heat-source preprocessing.
```

Current status:

```text
About 55%--70%.
```

Remaining:

1. three-body TEAM7 dp sweep;
2. sigma contrast sweep;
3. source/frequency sweep;
4. global + bodywise residual convergence criteria;
5. E/J/Joule observable convergence;
6. optional B=curl(A);
7. four-body cold-crucible contact scaffold:
   - air;
   - coil;
   - crucible;
   - melt;
8. melt / crucible Joule power diagnostic;
9. bodywise residual for all four bodies.

Only after this level should fixed sigma/nu 10D Joule → heat equation one-way thermal coupling be considered.

### 4.3 Level 3: TEAM7 / paper-grade validation ready

Definition:

```text
Suitable for serious physical validation in papers or benchmark chapters.
```

Current status:

```text
About 30%--45%.
```

Remaining:

1. real TEAM7 geometry;
2. external FEM/MFEM reference;
3. mesh/particle resolution convergence;
4. full B/E/J/Joule observable set;
5. line profile / cross-section profile;
6. plate total Joule power;
7. error tables;
8. physical parameter tables;
9. consistency statement with literature benchmark;
10. SPHinXsys Contact multi-body computation flowchart/pseudocode.

This level cannot be completed in the current sprint; it needs a separate phase.

---

## 5. Can Current Blockers Be Resolved?

Yes, and they do not appear to be "wrong direction" blockers.

### 5.1 left MMS 1% can be resolved or explained

This issue can likely be explained or converged via:

1. gauge/offset decomposition;
2. L2 vs Linf decomposition;
3. core vs interface band decomposition;
4. E/J/Joule observable comparison;
5. material swap / body order swap;
6. split-vs-split comparison;
7. bodywise residual-based stopping criteria;
8. possibly gauge pinning / bodywise mean constraint.

Even if left solution block max-norm cannot reach `1e-5`, as long as E/J/Joule and interface flux are correct, engineering use need not be blocked. Cold-crucible heating ultimately cares about E/J/Joule, not some gauge-dependent max-norm of A/phi.

### 5.2 true TEAM7 validation takes time, but route is clear

Real TEAM7 is not done with one code change. Main work:

- build geometry;
- build reference;
- unify parameters;
- compare physical quantities;
- run convergence study.

This is engineering effort, not evidence that Contact architecture is infeasible.

### 5.3 many helpers is not a root problem

Many helper files exist because the test matrix is complex. Merging is not recommended now. Reorganize after Contact milestone is fully complete.

---

## 6. Stage 10.5: left MMS Field Error Diagnostic Tasks

Cursor next step: create or extend diagnostics specifically for the left-body 1% issue.

### 6.1 Task 10.5-A: gauge / mean-offset diagnostic

Add output:

```text
left_A_mean_error
left_phi_mean_error
left_A_L2_error
left_phi_L2_error
left_A_Linf_error
left_phi_Linf_error

left_A_L2_after_mean_subtraction
left_phi_L2_after_mean_subtraction
left_A_Linf_after_mean_subtraction
left_phi_Linf_after_mean_subtraction
```

Assessment:

- If error drops sharply after mean-subtraction, mainly gauge/constant mode;
- If not, field shape error.

### 6.2 Task 10.5-B: E/J/Joule field error diagnostic

Add:

```text
left_E_L2_error
left_E_Linf_error
left_J_L2_error
left_J_Linf_error
left_Joule_power_error
left_Joule_L2_error
left_Joule_Linf_error
```

Assessment:

- If A/phi error ~1% but E/J/Joule error small, issue is gauge-biased, does not affect main heat source;
- If E/J/Joule also ~1%, further investigate Contact discretization.

### 6.3 Task 10.5-C: core/interface/boundary partitioned error

Split left body into:

```text
left_core_region
left_interface_band
left_outer_boundary_band
```

Output separately:

```text
L2_error
Linf_error
particle_count
max_error_location
```

Assessment:

- If error only in interface/boundary Linf, should not be main pass criterion;
- If core L2 also high, issue is more serious.

### 6.4 Task 10.5-D: material/body/order swap

Run three control cases:

```text
Case 1:
    original left/right materials and body order

Case 2:
    swap material labels

Case 3:
    swap body construction order / contact direction order
```

Assessment:

- Error follows left body: possible body/order/contact direction issue;
- Error follows conductor material: possible material stiffness / gauge issue;
- Error follows interface normal direction: check sign convention.

### 6.5 Task 10.5-E: split-vs-monolithic restricted comparison

Do not compare only analytic exact. Add:

```text
monolithic_inner_solution
contact_split_solution
restricted_monolithic_solution_on_left_particles
```

Compare:

```text
contact_left - analytic_exact
contact_left - restricted_monolithic_left
monolithic_left - analytic_exact
```

Assessment:

- If monolithic also deviates from analytic exact but contact is close to monolithic, analytic MMS is not split discrete exact;
- If contact deviates from monolithic, continue investigating Contact operator;
- If contact E/J/Joule is close to monolithic, engineering observable can be released first.

---

## 7. Three-Body TEAM7 Next-Step Reinforcement Tasks

Three-body TEAM7 has preliminary pass, but needs reinforcement.

### 7.1 Task T7-A: bodywise residual

In three-body TEAM7 scaffold output:

```text
global_true_rel
air_true_rel
coil_true_rel
plate_true_rel
max_bodywise_true_rel
```

Suggested pass criteria:

```text
global_true_rel < 1e-5 ~ 5e-5
max_bodywise_true_rel < 1e-4
```

Engineering phase can relax slightly, but cannot rely on global only.

### 7.2 Task T7-B: observable convergence

Output:

```text
plate_solution_block_max
plate_E_L2
plate_J_L2
plate_Joule_power
plate_Joule_max
coil_rhs_norm
coil_solution_block_max
air_solution_block_max
```

Compare:

```text
monolithic inner scaffold
split Contact scaffold
```

Suggested acceptance:

```text
plate_Joule_power gap < 1% initially
plate_E_L2 gap < 5%
plate_J_L2 gap < 5%
plate_Joule_max gap < 10%--20% diagnostic only
```

If current plate Joule gap ~0.02%, that item is already very good, but more parameter points need verification.

### 7.3 Task T7-C: dp sweep

At least:

```text
dp = coarse
dp = medium
dp = fine
```

Record:

```text
particle_count per body
GMRES outer / restart
global_true_rel
bodywise_true_rel
plate_Joule_power
plate_E_L2
plate_J_L2
runtime
```

Purpose:

- assess whether observables are stable with dp under Contact multi-body;
- assess whether three-body vs monolithic error converges with dp.

### 7.4 Task T7-D: sigma contrast sweep

At least:

```text
sigma_plate / sigma_air = 10^0
10^1
10^2
10^3
10^4
```

or reasonable range from current nondimensional parameters.

Output:

```text
condition/stiffness indicator
GMRES outer
global_true_rel
plate_Joule_power
plate_J_L2
bodywise residual
```

Purpose:

- verify Contact solver stability under high contrast;
- prepare early for cold-crucible crucible/melt high contrast.

### 7.5 Task T7-E: source/frequency sweep

Small range:

```text
omega = low / medium / high
J_source = weak / baseline / strong
```

Check:

```text
linear scaling of fields
quadratic scaling of Joule with source amplitude
solver robustness
```

Joule should scale quadratically with source amplitude; good sanity check.

---

## 8. When to Start Cold-Crucible Four-Body Contact

10D thermal coupling is not recommended immediately, but four-body EM scaffold can start after three-body TEAM7 reinforcement.

Suggested conditions:

```text
1. Three-body TEAM7 bodywise residual output complete;
2. Three-body plate Joule/E/J stable vs monolithic at at least 2--3 parameter points;
3. left MMS 1% at least preliminarily classified as gauge/max-norm/discrete mismatch;
4. helpers not merged yet, but documentation clear.
```

Four-body setup:

```text
air body
coil body
crucible body
melt body
```

First version only:

```text
EM solve
E/J/Joule output
bodywise residual
melt_Joule_power
crucible_Joule_power
coil_solution_block_max
air_solution_block_max
```

Do not do yet:

```text
sigma(T)
nu(T)
thermal feedback
Lorentz force
melting
fluid coupling
```

---

## 9. Why Are There Many Helper Files, and Is That Normal?

### 9.1 helpers are not SPHinXsys core physics types

These helpers are not core physics/dynamics types in the SPHinXsys framework like `LocalDynamics`, `InteractionDynamicsCK`, `StateDynamics`.

Helpers are essentially test support code for reusing:

- geometry construction;
- body/relation setup;
- MMS exact field;
- RHS assembly;
- GMRES options;
- diagnostic metrics;
- E/J/Joule observable;
- CSV / log output;
- monolithic vs split comparison;
- workspace contamination;
- readback sync;
- polish sweep.

You may not have noticed before because ordinary SPHinXsys examples often fit in one `.cpp` without this much cross-case reuse.

### 9.2 many helpers currently is reasonable

Current A--phi test matrix is complex:

```text
inner vs contact
monolithic vs split
two-body MMS
three-body TEAM7
single-body GMRES
coupled multi-body GMRES
workspace contamination
readback sync
polish sweep
E/J/Joule observables
benchmark naming
```

Without helpers, each test `.cpp` would repeat hundreds of lines, harder to maintain later.

Therefore:

```text
Many helpers is not an error.
Merging not recommended now.
But helpers must not enter production umbrella include.
```

### 9.3 boundary between production code and helpers

Production solver/operator should go in:

```text
src/shared/particle_dynamics/electromagnetic_dynamics/
```

Including:

```text
aphi_multibody_contact_gmres_ck.*
aphi_gmres_contact_solver_ck.*
aphi_matrix_free_contact_solve_ck.*
aphi_matrix_free_operator_ck.*
aphi_laplace_ck.*
aphi_grad_phi_coupling_ck.*
aphi_div_sigma_a_coupling_ck.*
aphi_joule_heating_ck.*
aphi_block_jacobi_preconditioner_ck.*
```

Test helpers should go in test directories:

```text
tests/extra_source_and_tests/3d_examples/.../helpers/
```

Or common test helper directory:

```text
tests/extra_source_and_tests/3d_examples/common/aphi/
```

### 9.4 not merging files now is correct

left MMS 1% not fully explained; three-body TEAM7 still needs reinforcement. Merging helpers or abstracting GMRES driver now risks new bugs.

Suggested order:

```text
Complete Contact solver semantics first;
then three-body TEAM7 reinforcement;
then Stage 10.5 left MMS field error attribution;
then start cold-crucible four-body EM scaffold;
finally reorganize helpers and abstract common driver.
```

---

## 10. Near-Term Execution Order

Cursor can execute in this order.

### Step 1: Update progress documentation

Update:

```text
CURSOR_APHI_STAGE10_CONTACT_PROGRESS.md
```

Must include:

```text
1. workspace contamination confirmed;
2. single-body Contact GMRES Arnoldi changed to block-diagonal;
3. coupled multi-body GMRES becomes main line;
4. three-body TEAM7 plate Joule gap ~0.02%;
5. left MMS max-norm ~1% is open diagnostic;
6. polish sweep 0→3 plateaus;
7. helpers not merged is deliberate decision;
8. current completion is scaffold/contact infrastructure, not real TEAM7 validation.
```

### Step 2: Add bodywise residual output to coupled GMRES

All multi-body Contact solver tests output:

```text
global_true_rel
max_bodywise_true_rel
bodywise_true_rel per body
bodywise_rhs_norm per body
bodywise_solution_norm per body
```

Including:

```text
two-body MMS
three-body TEAM7 scaffold
future cold-crucible four-body scaffold
```

### Step 3: Run Stage 10.5 left MMS diagnostic

Add or extend:

```text
test_3d_aphi_ck_contact_left_field_error_diagnostic
```

At minimum include:

```text
mean offset removal
L2 vs Linf
core/interface/boundary split
E/J/Joule error
material swap
body order swap
split-vs-monolithic restricted comparison
```

Goal is not immediately reaching left max-norm `1e-5`, but classifying it as:

```text
gauge/offset
local Linf
Contact split discrete mismatch
body/order bug
material stiffness
```

### Step 4: Reinforce three-body TEAM7 scaffold

On existing three-body TEAM7 coupled GMRES, add:

```text
dp sweep
sigma contrast sweep
source amplitude sweep
frequency sweep
bodywise residual
plate observable convergence
```

Pass criteria:

```text
plate Joule/E/J stable and consistent with monolithic scaffold;
bodywise residual shows no obvious anomaly;
source amplitude vs Joule quadratic scaling correct.
```

### Step 5: Start cold-crucible four-body EM scaffold

When conditions met, construct:

```text
air
coil
crucible
melt
```

Phase one EM only:

```text
A/phi solve
E/J/Joule output
bodywise residual
melt_Joule_power
crucible_Joule_power
```

Do not enter thermal coupling.

### Step 6: File reorganization and helper merge

Do last:

```text
1. move helpers to common test helper directory;
2. merge diagnostic helpers;
3. mark deprecated legacy block-GS;
4. abstract common GMRES driver;
5. sync documentation.
```

Do not do large-scale refactor before left MMS and three-body TEAM7 are stable.

---

## 11. Can We Say Contact Is Done Now?

Suggested wording:

### Can say

```text
Contact operator layer has been largely validated.
The coupled multi-body Contact GMRES prototype is operational.
Three-body TEAM7-like scaffold agrees with the monolithic scaffold in plate Joule/E/J observables.
The previous single-body GMRES workspace contamination issue has been identified and corrected by using a block-diagonal Arnoldi apply.
```

### Cannot say

```text
The full Contact A--phi solver is completely validated for real TEAM7.
The current TEAM7 scaffold is already a publication-grade benchmark.
The left-body MMS field error is fully understood.
The solver is ready for thermal coupling without further EM diagnostics.
The current solution_block_max is magnetic field.
```

### Most accurate project status

```text
Stage 10 Contact has completed operator and coupled solver prototype validation,
and is moving from scaffold equivalence into engineering validation.
Next steps: reinforce bodywise residual, three-body parameter sweeps, and left MMS field diagnostic.
Real TEAM7 benchmark and cold-crucible production case remain later phases.
```

---

## 12. Final Route Decision

Recommended route:

```text
Current:
    do not revert to legacy block-GS;
    do not rush helper merge;
    do not enter 10D thermal coupling directly.

Short term:
    bodywise residual + left MMS diagnostic + three-body TEAM7 sweep.

Mid term:
    cold-crucible four-body EM scaffold.

Later:
    fixed sigma/nu Joule -> thermal one-way coupling.

Final:
    external FEM/MFEM TEAM7 benchmark + B/E/J/Joule physical validation.
```

Core assessment:

> Current Contact work is not as early as "only inner comparison," but has not reached "real TEAM7 physics problem fully done." It is in the middle layer: multi-body Contact infrastructure and engineering scaffold are largely in place; next step is to push from scaffold toward benchmark-grade validation.
