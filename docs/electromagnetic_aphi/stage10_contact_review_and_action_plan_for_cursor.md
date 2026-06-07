# Stage 10 Contact Multi-Body A--phi Development Review and Follow-Up Execution Plan

This document summarizes the development status of Stage 10 Contact in the `electromagnetic_dynamics-contact` branch: what has passed, current blockers, code change suggestions, and follow-up Sprint execution order. Goal: let Cursor execute directly from this document with minimal extra route decisions.

---

## 0. Overall Conclusion

Current Contact progress is effective; direction is largely correct.

Can already consider:

- `Contact<>` pairwise Laplace / grad-phi / div-sigma-A / Joule / fused apply operator layer largely connected;
- Sprint 1--3 can serve as Contact operator foundation acceptance;
- Sprint 4 two-body interface MMS proves Contact operator and RHS assembly have no fundamental error;
- but current `block Gauss--Seidel + single-body Contact GMRES` is still only a transitional diagnostic path, not the final production solver.

What is truly not closed:

1. Multi-body coupled solver semantics;
2. Contact neighbor body variable state sync between device/host;
3. Single-body GMRES neighbor Krylov workspace contamination;
4. block-GS global residual and full sweep convergence control;
5. Three-body TEAM7-like scaffold plate Joule / E / J metric acceptance;
6. File organization and test helper placement.

Therefore do not enter formal 10D thermal coupling next; first add **Sprint 4.5** to fix block-GS Contact GMRES key issues, then enter Sprint 5 three-body TEAM7-like scaffold.

---

## 1. Current Development Progress Assessment

| Phase | Current Status | Conclusion |
|---|---:|---|
| Contact Laplace vs Inner | Complete | Pass |
| Contact Grad/Div/Joule apply | Complete | Pass |
| Contact fused apply + block Jacobi | Complete | Largely pass |
| Two-body Contact interface MMS | Passed but `converged=0` | Proof-of-concept, not production solver |
| Contact zero/pin diagnostic | Found clear issue | Valuable, not fully resolved |
| Three-body TEAM7-like scaffold | Not started | Next, but add Sprint 4.5 first |
| Coupled multi-body GMRES | Not done | Needed for formal route |
| E/J/Joule/B metric naming | Not unified | Handle before Sprint 5 |

Current stage can be stated as:

> Contact operator layer is mostly validated, but the multi-body A--phi solver layer is not yet finalized.

---

## 2. Is Completed Work Acceptable?

### 2.1 Sprint 1--3: Acceptable

Current Contact operator style overall matches SPHinXsys SYCL/CK conventions:

- Uses `LocalDynamics` / `Interaction<Inner>` / `Interaction<Contact>`;
- Kernels get variables via `DelegatedData(ex_policy)`;
- Contact contributions use owner-body gather, writing only current body `i` particles;
- No atomic writes to neighbor bodies;
- Uses `InteractionDynamicsCK` / `StateDynamics` / `ReduceDynamicsCK` for parallel organization;
- Inner and Contact specializations share same physical operator class system.

Following Contact operators can be considered largely acceptable at operator layer:

```cpp
AphiPairwiseLaplaceCK<Contact<>>
AphiGradPhiCouplingCK<Contact<>>
AphiDivSigmaACouplingCK<Contact<>>
AphiComputeScalarPhiGradientCK<Contact<>>
AphiApplyCK<Contact<>>
AphiComputeBlockJacobiDiagonalCK<Contact<>>
```

### 2.2 Sprint 4: Diagnostic Pass, Solver Not Formally Closed

Current Sprint 4 results:

```text
mono:
    outer = 2
    discrete_defect ≈ 5e-6
    continuous ≈ 1.8e-6

contact:
    outer = 30
    discrete_defect ≈ 4.6e-8
    continuous ≈ 4.9e-6
    true_rel ≈ 1.4e-5
    converged = 0
```

This shows:

- Contact operator itself likely has no fundamental error;
- RHS assembly is not the main issue;
- Two-body interface MMS solution can be recovered;
- but current `block-GS + single-body Contact GMRES` convergence control is not fully closed.

Therefore do not write Sprint 4 as:

> Contact GMRES production solver converged.

Write instead:

> Contact MMS passes the diagnostic accuracy criterion, but the single-body Contact GMRES path is not yet a fully converged production solver.

---

## 3. Current Main Blockers

### 3.1 After Zeroing Contact Neighbor Body, Pin Still Reads Wrong

Current diagnostic:

```text
After assembly u = exact, not zeroed:
    |Au - b| / |b| ≈ 1e-11    ✅

Zero both sides → pin r_hat:
    ≈ 0.64                    ❌

Zero left only, right keeps exact:
    ≈ 0.0016                  ⚠️
```

This shows the issue is not interface discretization formula, but:

> Contact neighbor body variable state is inconsistent on host/device/contact apply side.

Short-term strategy:

- In block-GS, **do not zero Contact neighbor body solution block**;
- During current body solve, neighbor body solution should be current block-GS iterate / boundary state;
- Only allow zeroing current solving body's residual/search/workspace.

Long-term strategy:

- Must find why pin back to host/device after zero still makes Contact apply read wrong;
- This should not be treated as normal framework behavior;
- At least need minimal Contact readback test to clarify source.

### 3.2 Single-Body Contact GMRES May Be Contaminated by Neighbor Krylov Workspace

This is the highest-priority fix needed now.

Current `AphiGMRESContactSolverCK` in Arnoldi uses:

```cpp
AphiApplyContactDynamicsBundle<ExecutionPolicy> apply_z_to_w(
    body_, inner_, contact_, z_j, names_.search, ...
);
```

This full Contact apply reads:

- Current body `z_j`;
- Contact neighbor body `z_j`.

But for single-body block-GS, current body Krylov vector should only represent current block perturbation. Mathematically, current body GMRES should solve:

```text
A_ii z_i
```

Not:

```text
A_ii z_i + A_ij z_j
```

If neighbor body `GMRESZk` is not explicitly zeroed, or was just used by another body's GMRES, current body Arnoldi operator is contaminated.

Typical risk:

1. Left-body GMRES runs first; left `GMRESZ0, GMRESZ1, ...` written;
2. Right-body GMRES runs next;
3. Right Contact apply reads left `GMRESZk`;
4. If left `GMRESZk` nonzero, right Krylov operator is not single-block operator;
5. Causes high outer count, `converged=0`, or recursive vs true residual mismatch.

Therefore single-body Contact GMRES should not continue using full Contact apply directly as Arnoldi apply.

---

## 4. Required Fix: Sprint 4.5

Before Sprint 5 three-body TEAM7-like scaffold, add Sprint 4.5.

---

### 4.5-A: Add Workspace Contamination Test

Purpose:

> Verify whether single-body Contact GMRES Arnoldi apply is affected by neighbor body Krylov workspace.

Suggested new test:

```cpp
test_3d_aphi_ck_contact_gmres_workspace_contamination.cpp
```

Test logic:

1. Construct two-body Contact interface case;
2. Set current body `GMRESZ0` to known analytic vector;
3. Set neighbor body `GMRESZ0` to zero; run `apply_z_to_w`, record `w_zero_neighbor`;
4. Set neighbor body `GMRESZ0` to random or constant nonzero; rerun `apply_z_to_w`, record `w_dirty_neighbor`;
5. Compare:

```text
||w_dirty_neighbor - w_zero_neighbor|| / ||w_zero_neighbor||
```

Expected:

- If still using full Contact apply, test likely fails;
- Failure shows single-body GMRES contaminated by neighbor Krylov workspace;
- After fix, difference should be near machine precision or at least below `1e-12 ~ 1e-10`.

Acceptance standard suggestion:

```text
workspace_contamination_rel_diff < 1e-10
```

In double precision with reduction order effects, can relax to:

```text
workspace_contamination_rel_diff < 1e-8
```

---

### 4.5-B: Add Block-Diagonal Contact Apply

Purpose:

> Single-body GMRES Arnoldi apply should only apply current body diagonal block `A_ii`, not read neighbor body Krylov input.

Suggested new class:

```cpp
template <class DataDelegationType>
class AphiApplyContactBlockDiagonalCK;
```

Suggested bundle:

```cpp
template <class ExecutionPolicy>
class AphiApplyContactBlockDiagonalDynamicsBundle;
```

For `AphiGMRESContactSolverCK` internal `apply_z_to_w`.

#### 4.5-B-1. Difference Between Full and Block-Diagonal Contact Apply

Full Contact apply for residual / RHS / `K(solution)`:

```text
K_i(x_i, x_j) = A_ii x_i + A_ij x_j
```

Block-diagonal Contact apply for single-body GMRES Arnoldi:

```text
K_i^{diag}(z_i) = A_ii z_i
```

Therefore:

- Residual still uses full Contact apply;
- RHS / MMS apply(exact) still uses full Contact apply;
- Only `apply_z_to_w` uses block-diagonal Contact apply.

#### 4.5-B-2. Contact Laplace Block-Diagonal Form

Full pairwise Laplace equivalent to:

```text
w_ij (u_i - u_j)
```

Block-diagonal apply keeps only current body unknown:

```text
w_ij u_i
```

Implementation:

- Still loop contact neighbors;
- Still compute pairwise weight;
- Do not read neighbor input block `u_j`;
- Accumulate only diagonal contribution to owner body output.

#### 4.5-B-3. Grad-Phi Coupling Block-Diagonal Form

Full form:

```text
sigma_i * g_ij * (phi_i - phi_j)
```

Block-diagonal form:

```text
sigma_i * g_ij * phi_i
```

Note:

- Do not read contact neighbor `phi_j`;
- Keep same sign convention as current full grad-phi operator;
- Can verify sign via finite difference vs full operator diagonal part.

#### 4.5-B-4. Div-Sigma-A Coupling Block-Diagonal Form

Full form:

```text
omega * sigma_ij * g_ij · (A_i - A_j)
```

Block-diagonal form:

```text
omega * sigma_ij * g_ij · A_i
```

Note:

- `sigma_ij` still per current harmonic/face average rule;
- Do not read neighbor body `A_j`;
- Keep same sign convention as full div-sigma-A operator.

#### 4.5-B-5. A/phi Mass-Like Term

For local term:

```text
j omega sigma A
```

Already body-local; no Contact neighbor. Block-diagonal apply same as full apply.

#### 4.5-B-6. Replacement Location in GMRES

In `AphiGMRESContactSolverCK`:

- Residual / true residual / RHS consistency:
  - Keep full `AphiApplyContactDynamicsBundle`
- Arnoldi `apply_z_to_w`:
  - Change to `AphiApplyContactBlockDiagonalDynamicsBundle`

i.e.:

```cpp
// residual: full coupled operator for current body
AphiApplyContactDynamicsBundle<ExecutionPolicy> apply_solution_to_lhs(...);

// Arnoldi: diagonal block operator only
AphiApplyContactBlockDiagonalDynamicsBundle<ExecutionPolicy> apply_z_to_w(...);
```

Do not use full Contact apply directly on `GMRESZk` for Arnoldi unless forcing all neighbor body `GMRESZk` to zero. Block-diagonal apply is clearer algorithmically.

---

### 4.5-C: Add Minimal Contact Readback / Zero-Pin Sync Test

Purpose:

> Confirm Contact apply reads neighbor variables consistent with host/device current state.

Suggested new test:

```cpp
test_3d_contact_variable_readback_sync.cpp
```

No A--phi or GMRES needed. Minimal Contact variable read/write only.

Steps:

1. Construct left body / right body;
2. On right body register scalar variable, e.g.:
   ```cpp
   NeighborProbePhi
   ```
3. Use device dynamics to write right body `NeighborProbePhi` to known value;
4. Left body reads right body `NeighborProbePhi` via Contact, writes to left body `ReadNeighborPhi`;
5. Host checks left body `ReadNeighborPhi` equals known value;
6. Test separately:

```text
Case A: device write → contact read, no finalizeLoadIn
Case B: device write → finalizeLoadIn → contact read
Case C: host write → finalizeLoadIn → contact read
Case D: device zero → device copy/pin → contact read
Case E: host zero → host pin → finalizeLoadIn → contact read
```

Must clarify:

- Whether `finalizeLoadIn(execution::par_device)` overwrites device with stale host data;
- After device-side copy, should not call host-to-device sync;
- Whether Contact kernel uses stale delegated data pointer;
- Whether Contact relation update is unrelated to variable buffer sync.

Expected:

- After host modify, correct host-to-device sync must make contact read visible;
- After device modify, do not call host-to-device load that overwrites device;
- If Contact still reads old value after device write, `DelegatedData` or contact variable pointer has deeper issue.

---

## 5. block-GS Solver Contract

Write block-GS multi-body solve rules as documentation or code comments to avoid repeated mistakes.

Suggested contract:

### 5.1 Solution Block Rules

```text
1. During single-body solve, only modify current solving body solution block.
2. Do not zero Contact neighbor body solution block.
3. Neighbor body solution represents current block-GS iterate / boundary state.
4. Initially can set unified initial value for all bodies, e.g. zero or exact/MMS initial guess.
5. After each block sweep, all body solutions together form current global iterate.
```

### 5.2 RHS Block Rules

```text
1. Each body stores its own RHS independently.
2. In MMS, RHS = K(exact) must be assembled with all bodies in exact solution state.
3. For coupled RHS, must apply exact state on all bodies, not single-side assembly.
4. RHS must not be overwritten by zero operations in block-GS.
```

### 5.3 Residual / LHS / Search Block Rules

```text
1. Current body residual/search/lhs/workspace can be freely zeroed.
2. Contact neighbor body solution can be read.
3. Contact neighbor body Krylov workspace must not participate in single-body Arnoldi.
4. If using full contact apply as Arnoldi apply, must ensure all neighbor workspace zero.
5. Prefer block-diagonal contact apply to avoid reading neighbor Krylov workspace entirely.
```

### 5.4 Contact Relation Rules

```text
1. Must update relation after body position / cell linked list changes.
2. Pure variable changes should not rely on update relation to fix.
3. Variable host/device sync direction must be explicit:
   - host modify → host-to-device load;
   - device modify → do not overwrite device with stale host;
   - device-to-host copy only before output.
4. Contact apply should read current device variable state.
```

### 5.5 Convergence Criteria Rules

```text
1. block-GS cannot look only at per-body residual.
2. After each block sweep must record global true residual:
   global_norm = sqrt(sum over bodies ||r_body||^2)
3. Also record bodywise true residual to locate stuck body.
4. Physical quantities like Joule power / E / J should converge with block sweep.
```

---

## 6. How to Handle `converged=0 / outer=30 / true_rel≈1.4e-5`

Current result acceptable as Sprint 4 diagnostic, not as production solver acceptance.

Suggested test distinguishes two flags:

```text
contact_converged = GMRES solver actually reached requested tolerance
accepted_by_diagnostic = final true residual and MMS error satisfy relaxed diagnostic criteria
```

Current result can be:

```text
contact_converged = 0
accepted_by_diagnostic = 1
```

Suggested temporary diagnostic thresholds:

```text
final_true_relative_residual < 10 * tolerance
discrete_mms_defect < 10 * tolerance
continuous_error < 5e-5
```

But final production thresholds should require:

```text
contact_converged = 1
final_true_relative_residual < tolerance or specified engineering tolerance
outer < max_outer
```

After Sprint 4.5 fix, rerun with:

```text
max_outer = 30, 50, 80
restart_dimension = 30, 50, 80
```

Observe whether true residual naturally drops below `1e-5`.

---

## 7. Sprint 5: Three-Body TEAM7-Like Scaffold

Sprint 5 should continue, but after Sprint 4.5.

### 7.1 Three-Body Setup

Suggested construction:

```text
air body
coil body
plate/conductor body
```

Minimum functionality:

```text
coil body:
    impressed current / source region

plate body:
    solve induced A/phi response
    compute E, J, Joule

air body:
    participates in EM field propagation
```

### 7.2 Acceptance Order

Do not run full GMRES first. Accept in this order.

#### Step 5.1: Apply-Level Monolithic vs Split

Setup:

```text
monolithic inner air/coil/plate region-tag case
split Contact air + coil + plate case
```

For same analytic/separable solution field, compare:

```text
K(X)
Laplace contribution
grad-phi contribution
div-sigma-A contribution
Joule source
```

Suggested metrics:

```text
core-away relative L2 diff < 1e-3 ~ 1e-2
interface-band relative diff < 1e-2
regional integrated Joule gap < 1% ~ 3%
```

#### Step 5.2: Joule / Current Diagnostics

Output and compare in plate region:

```text
plate_Joule_power
plate_Joule_max
plate_Joule_L2
plate_J_L2
plate_E_L2
plate_solution_block_max
```

Suggested thresholds:

```text
plate_Joule_power gap < 1% ~ 3%
plate_J_L2 gap < 5%
plate_Joule_max gap < 10% ~ 20%
```

Note:

- `max` metrics easily affected by local particle error; not sole pass criterion;
- `integrated Joule power` and `L2` better as primary metrics.

#### Step 5.3: Solver-Level block-GS

If using block-GS:

```text
for sweep = 1..max_sweeps:
    solve air body with neighbors fixed
    solve coil body with neighbors fixed
    solve plate body with neighbors fixed
    compute global true residual over all bodies
    compute plate Joule power
```

Output:

```text
global_true_rel
bodywise_true_rel_air
bodywise_true_rel_coil
bodywise_true_rel_plate
plate_Joule_power
plate_Joule_power_change_per_sweep
```

Suggested acceptance:

```text
global_true_rel < 5e-4 or 1e-4
plate_Joule_power relative change over last 2 sweeps < 1e-3 ~ 1e-2
bodywise residual no body stuck at O(1)
```

---

## 8. Metric Naming Unification

Must unify before Sprint 5 to avoid scaffold rework.

### 8.1 Must Rename

Current names like:

```text
conductor_field_max
melt_field_max
crucible_field_max
coil_field_max
```

If counting `A_real, A_imag, phi_real, phi_imag` solution block max, must change to:

```text
conductor_solution_block_max
melt_solution_block_max
crucible_solution_block_max
coil_solution_block_max
```

Reason:

> These are not magnetic field B, nor electric field E, nor current density J.

### 8.2 Physical Metrics to Add

Suggested in Joule post-processing:

```text
E_real_max
E_imag_max
E_L2

J_real_max
J_imag_max
J_L2

Joule_power
Joule_max
Joule_L2
```

Output by region/body:

```text
plate_E_L2
plate_J_L2
plate_Joule_power
plate_Joule_max

coil_rhs_norm
coil_solution_block_max

air_solution_block_max
global_true_rel
bodywise_true_rel
```

### 8.3 Future B/H Metrics

If `B = curl(A)` not implemented, do not output or name as magnetic field.

Can add separately later:

```text
B_real_max
B_imag_max
B_L2

H_real_max
H_imag_max
H_L2
```

Only discuss magnetic field spatial distribution after implementation.

---

## 9. Code Structure and File Organization Suggestions

Many source files currently, but no large refactor immediately. Small-step organization preferred; do not break tests.

### 9.1 Keep Current Physical Operator File Structure

Following structure is reasonable; do not force merge:

```text
aphi_laplace_ck.*
aphi_grad_phi_coupling_ck.*
aphi_div_sigma_a_coupling_ck.*
aphi_joule_heating_ck.*
aphi_matrix_free_operator_ck.*
aphi_block_jacobi_preconditioner_ck.*
```

Principle:

> Same physical operator file provides Inner and Contact specializations.

Better than splitting inner/contact into parallel file sets.

### 9.2 Suggested Test Helper Moves

If following files are under production electromagnetic dynamics directory, move to test helper directory:

```text
aphi_contact_test_helpers.h
aphi_contact_gmres_test_helpers.h
aphi_contact_interface_diagnostic_helpers.h
aphi_test_device_sync.h
aphi_lhs_test_helpers.h
aphi_gmres_test_helpers.h
aphi_gmres_benchmark_helpers.h
aphi_gmres_robustness_sweep_helpers.h
aphi_scalar_phi_diagnostic_helpers.h
```

Suggested target:

```text
tests_extra/electromagnetic_dynamics/test_helpers/
```

Or:

```text
tests_extra/test_3d_aphi_ck/helpers/
```

Principle:

> Production dynamics only holds reusable solver/operator; test case geometry, diagnostics, CSV output, MMS helpers should not pile in production directory long term.

### 9.3 `aphi_contact_jacobi_debug_ck` Should Be Renamed or Merged

If used by `AphiGMRESContactSolverCK`, it is not a debug file.

Suggested rename:

```text
aphi_contact_block_jacobi_preconditioner_ck.*
```

Or merge into:

```text
aphi_block_jacobi_preconditioner_ck.*
```

Because Contact specialization and Contact diagonal contribution belong to block Jacobi preconditioner.

### 9.4 Inner/Contact GMRES Merge Not Urgent

Current:

```text
aphi_gmres_solver_ck.*
aphi_gmres_contact_solver_ck.*
```

May have high duplication, but do not merge now. Reasons:

- Contact solver block-diagonal apply not fixed yet;
- Coupled multi-body GMRES not designed yet;
- Merging now may solidify wrong abstractions.

After Sprint 4.5 and Sprint 5, consider abstracting to:

```cpp
template <class ExecutionPolicy, class ApplyBundleType, class JacobiBundleType>
class AphiGMRESDriverCK;
```

Inner and contact only provide different:

```text
apply_solution_to_lhs
apply_krylov_vector_to_search
compute_jacobi_diagonal
```

### 9.5 Future Unified Wrapper

If currently have:

```text
aphi_matrix_free_solve_ck.*
aphi_matrix_free_contact_solve_ck.*
```

Can later unify common solve driver; keep two wrappers for now.

---

## 10. Follow-Up Route Recommendations

Suggested route:

```text
Sprint 4.5:
    workspace contamination test
    block-diagonal contact apply
    minimal Contact readback / zero-pin sync test
    rerun two-body interface MMS

Sprint 5:
    three-body TEAM7-like contact scaffold
    monolithic vs split apply-level comparison
    plate E/J/Joule diagnostics
    block-GS solver-level global residual

Sprint 6:
    coupled multi-body GMRES prototype
    global Arnoldi over all bodies
    global dot/norm/update
    full Inner + Contact apply for all bodies

Sprint 7:
    physical cold-crucible contact scaffold
    air / coil / crucible / melt
    segmented crucible or slit wall if possible
    fixed sigma/nu Joule distribution

Stage 10D:
    one-way thermal coupling
    fixed sigma/nu
    Joule -> heat equation RHS
    energy conservation check
```

Do not enter now:

```text
sigma(T)
nu(T)
Lorentz force
melting/phase change
fully coupled EM-thermal-fluid
```

Leave for 10E or later.

---

## 11. Cursor Direct Execution Task List

Execute in order; do not jump to formal 10D.

### Task 1: Add Workspace Contamination Test

Suggested file:

```text
test_3d_aphi_ck_contact_gmres_workspace_contamination.cpp
```

Requirements:

- Two-body Contact case;
- Current body `GMRESZ0` fixed;
- Neighbor body `GMRESZ0` zero and nonzero respectively;
- Compare `apply_z_to_w`;
- If fails now, record as expected failure;
- Must pass after block-diagonal apply fix.

### Task 2: Implement Block-Diagonal Contact Apply

Add:

```cpp
AphiApplyContactBlockDiagonalCK
AphiApplyContactBlockDiagonalDynamicsBundle
```

Requirements:

- For single-body GMRES Arnoldi;
- Do not read Contact neighbor body Krylov input;
- Contact Laplace accumulates diagonal term only;
- Contact grad-phi accumulates owner phi term only;
- Contact div-sigma-A accumulates owner A term only;
- Local mass-like term unchanged.

### Task 3: Modify Contact GMRES

In `AphiGMRESContactSolverCK`:

- Residual / true residual continue full Contact apply;
- RHS consistency continue full Contact apply;
- `apply_z_to_w` use block-diagonal Contact apply;
- Zero current body GMRES workspace before each body solve;
- Do not rely on neighbor workspace zero for correctness.

### Task 4: Add Minimal Contact Sync Test

Suggested file:

```text
test_3d_contact_variable_readback_sync.cpp
```

Requirements:

- No A--phi dependency;
- Only test Contact reading neighbor variables;
- Cover device write, host write, zero, pin, finalizeLoadIn combinations;
- Clarify whether zero/pin issue is helper misuse or underlying Contact/device cache issue.

### Task 5: Rerun Two-Body Interface MMS

Output required:

```text
contact_converged
accepted_by_diagnostic
outer
true_rel
discrete_defect
continuous_error
defect_ratio
continuous_ratio
```

Target:

```text
contact_converged = 1
true_rel < 1e-5
outer < max_outer
```

If still not reached, record reason and output block-GS residual history.

### Task 6: Enter Sprint 5 Three-Body Scaffold

Construct:

```text
air body
coil body
plate body
```

Apply-level first:

```text
monolithic inner vs split Contact
```

Then Joule:

```text
plate_Joule_power
plate_J_L2
plate_E_L2
plate_Joule_max
```

Finally solver:

```text
global_true_rel
bodywise_true_rel
plate_Joule_power convergence over block sweeps
```

### Task 7: Unify Metric Naming

Change immediately:

```text
*_field_max
```

To:

```text
*_solution_block_max
```

And add:

```text
*_E_L2
*_J_L2
*_Joule_power
*_Joule_max
```

Do not call solution block max magnetic field.

### Task 8: Organize File Structure

Light organization first:

- Move test helpers to test directory;
- Rename or merge `aphi_contact_jacobi_debug_ck`;
- Do not merge inner/contact GMRES yet;
- Abstract common GMRES driver after Sprint 5 stable.

---

## 12. Final Assessment

Contact development has not gone off track. Correct assessment:

```text
Sprint 1--3:
    Contact operator layer acceptable.

Sprint 4:
    Contact interface MMS diagnostic passed, but solver not production.

Current blockers:
    multi-body solver semantics,
    zero/pin device state,
    neighbor Krylov workspace contamination,
    missing global residual.

Next:
    Sprint 4.5 fix solver semantics first,
    then Sprint 5 three-body TEAM7-like,
    then coupled multi-body GMRES and 10D thermal coupling.
```

Most important execution principle:

> Contact operators are largely connected, but multi-body A--phi has not formed a true coupled solver yet. Current block-GS can serve as transitional route, but must add block-diagonal Contact apply, workspace hygiene, global residual, and minimal Contact sync test before entering three-body TEAM7-like and cold-crucible scaffold.
