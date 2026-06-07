# Stage 10 A--phi / Contact Pre-Discussion Organization and Demo Execution Plan

This document is for Cursor execution. Goal: before the afternoon discussion with the advisor, organize current A--phi / Contact development from "many source files, many test cases, many diagnostics" into a clear presentation package.
This demo does not aim to cover all latest blockers, and should not use `divA` / A-divergence penalty open issues as the mainline. Focus on:

1. Operators are kernel/particle form, not grid finite differences;
2. Matrix-free operator can stably apply;
3. GMRES solver converges based on matrix-free apply;
4. Contact multi-body version can reproduce monolithic / inner scaffold;
5. Can already output E/J/Joule visualization results;
6. Will continue relaxed particles / sphere/cylinder / divA / curlA enhanced validation.

---

## 0. Principles for This Organization

### 0.1 No Large Refactor Before Afternoon

Current source has many files, but no large moves/merges before afternoon because:

```text
1. Contact / GMRES / diagnostic just stabilized;
2. Large helper moves easily break include / CMake;
3. Advisor cares more about "method and results credibility" than final directory structure;
4. Code organization can first be "indexed, grouped, presentation-ready"; not "destructive merge".
```

Therefore two levels:

```text
Level A: Lightweight organization required before afternoon
    - Build file map;
    - Keep few demo cases;
    - Prepare unified output;
    - Write clear purpose and expected results per case;
    - Do not or minimally touch production source.

Level B: Real merge after afternoon
    - Helper classification moves;
    - Diagnostic helper merges;
    - Common GMRES driver abstraction;
    - Mark obsolete cases deprecated.
```

### 0.2 Demo Should Not Cover All Latest Blockers

Afternoon mainline suggested:

```text
Operator reliability → matrix-free solver → Contact multi-body → TEAM7-like observable → follow-up relaxed-particle validation
```

Do not put mainline on:

```text
A-divergence penalty currently not converging;
whether to do projection;
PC not matched yet;
too many helper files.
```

These can be "follow-up work" or when advisor asks.

---

# Part I. Source File Organization Plan

## 1. Why Current Source Files Look Messy

Many files not because code is bad, but because Stage 10 simultaneously includes:

```text
1. production operator;
2. matrix-free GMRES solver;
3. Contact multi-body solver;
4. Joule/E/J observable;
5. phi gauge penalty;
6. divA diagnostic;
7. A-divergence penalty prototype;
8. workspace/readback/interface diagnostics;
9. two-body MMS;
10. three-body TEAM7 scaffold;
11. benchmark helper;
12. historical debug cases.
```

Mixed together, directory looks scattered. Need "explainable grouping" first.

---

## 2. Suggested New File Map Before Afternoon

New document:

```text
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/STAGE10_APHI_FILE_MAP_FOR_DISCUSSION.md
```

Content structure below.

### 2.1 Production / Reusable Operator and Solver

Files to present as "core code":

```text
extra_src/shared/electromagnetic_dynamics/
    aphi_coupling_modes_ck.h
        Parameter switches: phi_gauge_penalty, a_divergence_penalty, etc.

    aphi_matrix_free_operator_ck.h/.hpp
        Matrix-free fused apply main entry.

    aphi_laplace_ck.h/.hpp
        A / phi Laplace pairwise operator.

    aphi_grad_phi_coupling_ck.h/.hpp
        sigma grad(phi) coupling in A equation.

    aphi_div_sigma_a_coupling_ck.h/.hpp
        div(sigma A) coupling in phi equation for current continuity.

    aphi_phi_gauge_penalty_ck.h/.hpp
        phi gauge regularization: LhsPhi += lambda_phi * phi.

    aphi_joule_heating_ck.h/.hpp
        E/J/Joule post-processing.

    aphi_block_jacobi_preconditioner_ck.h/.hpp
        block Jacobi / approximate preconditioner.

    aphi_gmres_solver_ck.h/.hpp
        inner / single-body GMRES base.

    aphi_multibody_contact_gmres_ck.h/.hpp
        coupled multi-body Contact GMRES mainline.

    aphi_gmres_contact_solver_ck.h/.hpp
        single-body Contact GMRES / diagnostic / polish path.

    aphi_matrix_free_contact_solve_ck.h/.hpp
        Contact matrix-free solve wrapper.
```

### 2.2 Diagnostics / Not Yet Production

Do not present as afternoon mainline; only as "ongoing research issues":

```text
    aphi_gradient_divergence_debug_ck.*
        Pipeline computing gradA / divA.

    aphi_a_divergence_penalty_pipeline.h
        Stage 10.7 A-divergence penalty prototype.
        Still in sign/source/PC diagnostic stage; not production solver.

    aphi_div_a_diagnostic_helpers.h
        divA metric helper; diagnostic.
```

### 2.3 Test Helpers

No need to explain each to advisor; just say:

```text
These avoid repeating geometry, MMS, output, and metrics across test cases.
They are not SPHinXsys core framework types.
```

Suggested functional labels:

```text
Geometry / setup:
    aphi_contact_test_helpers.h

Two-body MMS:
    aphi_contact_gmres_test_helpers.h

Three-body TEAM7-like:
    aphi_team7_contact_test_helpers.h

Observable:
    aphi_em_observable_helpers.h

Diagnostic:
    aphi_contact_workspace_contamination_helpers.h
    aphi_contact_readback_sync_helpers.h
    aphi_contact_interface_diagnostic_helpers.h
    aphi_div_a_diagnostic_helpers.h

Generic:
    aphi_gmres_test_helpers.h
    aphi_lhs_test_helpers.h
    aphi_test_device_sync.h
```

---

## 3. Safe Lightweight Merge Before Afternoon

### 3.1 Can Merge Docs, Not Production Code

Before afternoon can merge/add:

```text
STAGE10_APHI_FILE_MAP_FOR_DISCUSSION.md
STAGE10_APHI_TEACHER_DEMO_CASES.md
STAGE10_APHI_TEACHER_DEMO_SUMMARY.md
```

Not recommended before afternoon:

```text
aphi_matrix_free_operator_ck.*
aphi_multibody_contact_gmres_ck.*
aphi_gmres_contact_solver_ck.*
aphi_block_jacobi_preconditioner_ck.*
```

Moving these easily breaks build.

### 3.2 Optional Helper Merge, Not Required

If Cursor has time, can merge diagnostic helpers only, not main solver:

```text
Before merge:
    aphi_contact_workspace_contamination_helpers.h
    aphi_contact_readback_sync_helpers.h
    aphi_contact_interface_diagnostic_helpers.h

After merge:
    aphi_contact_diagnostic_helpers.h
```

If time tight, skip. Afternoon demo does not depend on this merge.

---

# Part II. Recommended Demo Case Set for Afternoon

## 4. Demo Case Overall Principles

Do not show all cases. At most 4 main cases + 1 optional backup.

Recommended mainline:

```text
Case 1: Operator equivalence / Contact apply vs monolithic
Case 2: Two-body interface MMS + coupled GMRES
Case 3: phi gauge penalty sweep
Case 4: Three-body TEAM7-like scaffold, show E/J/Joule VTP
Optional: relaxed sphere/cylinder geometry test (if time before afternoon)
```

Not for main demo:

```text
A-divergence penalty λ_A sweep;
source divergence diagnostic;
projection;
helper merge;
cold-crucible four-body;
10D thermal coupling.
```

---

## 5. Case 1: Contact Apply vs Monolithic

### 5.1 Recommended Test

```text
test_3d_aphi_ck_contact_apply_vs_monolithic
```

If VTP output inconvenient, use:

```text
test_3d_aphi_ck_contact_fused_apply_equivalence
```

As stdout-only operator equivalence test.

### 5.2 Demo Purpose

Tell advisor:

```text
Same A--phi operator can get consistent matrix-free apply on monolithic inner body and split Contact multi-body.
```

This shows:

```text
Contact pairwise accumulation correct;
Inner/Contact kernel form consistent;
Not assembled from explicit large matrix.
```

### 5.3 Metrics to Watch

stdout:

```text
lhs_max_diff
lhs_L2_diff
relative_diff
passed
```

If component metrics exist, prefer:

```text
LaplaceA diff
LaplacePhi diff
gradPhi coupling diff
divSigmaA coupling diff
Joule diff
```

### 5.4 Output Needed

Cursor should save:

```text
stage10_teacher_demo/case01_contact_apply_vs_monolithic/stdout.txt
stage10_teacher_demo/case01_contact_apply_vs_monolithic/summary.md
```

VTP optional. Case mainly shows numerical consistency; visualization not required.

---

## 6. Case 2: Two-Body Flux-Matched Interface MMS + Coupled GMRES

### 6.1 Recommended Test

```text
test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured
```

### 6.2 Demo Purpose

Tell advisor:

```text
We not only verify operator apply, but put matrix-free operator in coupled GMRES;
Under two-material interface, RHS / solution / interface band all close.
```

This case shows:

```text
1. Variable-material interface handled;
2. Contact multi-body GMRES converges;
3. interface-band residual very small;
4. bodywise residual output available;
5. phi gauge penalty controls raw phi drift.
```

### 6.3 Metrics to Watch

stdout focus:

```text
contact_converged = 1
contact_outer
contact_global_true_rel
contact_max_bodywise_true_rel
contact_interface_band_rel
contact_discrete_defect
contact_left_continuous
contact_right_continuous
passed = 1
```

Explain to advisor:

```text
global_true_rel:
    GMRES equation residual.

max_bodywise_true_rel:
    Worst body residual in multi-body; prevents global norm masking local issues.

interface_band_rel:
    Residual near interface; shows Contact interface has no obvious error.

left/right continuous:
    Manufactured exact field error; raw field affected by gauge; diagnostic.
```

### 6.4 VTP Output Suggestion

If VTP exists, ensure at least:

```text
A_real / A_imag or solution block
Phi_real / Phi_imag
RHS / LHS / residual
MaterialId or BodyRegion
ExactSolution
ErrorNorm or ErrorMagnitude
InterfaceBand marker
```

If missing, Cursor can temporarily add minimum fields:

```text
SolutionBlockNorm
ResidualNorm
MaterialId
InterfaceBand
```

Output directory:

```text
stage10_teacher_demo/case02_two_body_interface_mms/vtp/
stage10_teacher_demo/case02_two_body_interface_mms/stdout.txt
stage10_teacher_demo/case02_two_body_interface_mms/summary.md
```

---

## 7. Case 3: Phi Gauge Penalty Sweep

### 7.1 Recommended Test

```text
test_3d_aphi_ck_contact_phi_gauge_penalty_sweep_diagnostic
```

### 7.2 Demo Purpose

Present not as "problem" but as:

```text
A--phi equations have gauge DOF;
We identified and controlled phi constant gauge;
After lambda_phi=100 raw phi/field drift significantly reduced.
```

Shows team checks gauge stability, not only GMRES residual.

### 7.3 Metrics to Watch

Focus:

```text
lambda_phi = 0:
    global_true_rel very small, but left_continuous huge

lambda_phi = 100:
    global_true_rel very small, left_continuous significantly reduced
```

Explanation:

```text
Good Krylov residual does not mean raw field necessarily comparable;
phi gauge penalty is necessary numerical gauge regularization;
phi penalty is not divA=0; can mention in follow-up work.
```

### 7.4 Suggested Output Table

Cursor generates:

```text
stage10_teacher_demo/case03_phi_gauge_penalty/phi_penalty_summary.csv
stage10_teacher_demo/case03_phi_gauge_penalty/summary.md
```

CSV columns:

```text
lambda_phi
converged
outer
global_true_rel
left_true_rel
left_continuous
phi_real_mean
phi_L2_mean_sub
div_A_relative
```

If `div_A_relative` available, add to show:

```text
phi penalty controls phi gauge but barely lowers divA.
```

Do not expand too much; avoid pulling discussion into current blockers.

---

## 8. Case 4: Three-Body TEAM7-Like Scaffold

### 8.1 Recommended Test

```text
test_3d_aphi_ck_gmres_contact_team7_three_body_scaffold
```

### 8.2 Demo Purpose

Best afternoon "engineering result":

```text
air / coil / plate three-body Contact;
coupled multi-body GMRES;
impressed current source;
plate outputs E/J/Joule;
Contact plate Joule/J observables align with monolithic scaffold.
```

Shows:

```text
1. Moved from single inner to multi-body Contact;
2. Not only MMS; can run TEAM7-like physical scaffold;
3. Has visualization results;
4. Outputs are E/J/Joule for subsequent electromagnetic heating.
```

### 8.3 Metrics to Watch

stdout focus:

```text
mono_plate_joule_power
coupled_plate_joule_power
plate_joule_power_gap

mono_plate_J_L2
coupled_plate_J_L2
plate_j_L2_gap

mono_plate_joule_max
coupled_plate_joule_max
plate_joule_max_gap

coupled_converged
coupled_outer
coupled_global_true_rel
coupled_max_bodywise_true_rel
air_true_rel
coil_true_rel
plate_true_rel
```

Talking points:

```text
plate_joule_power_gap very small:
    Contact three-body aligns with monolithic scaffold.

plate_J_L2_gap very small:
    Current density distribution also aligns.

bodywise residual:
    No single body masked by global norm.
```

### 8.4 VTP Output Suggestion

Best case for ParaView.

Ensure VTP per body:

```text
air.vtp
coil.vtp
plate.vtp
```

Or equivalent multi-body VTP.

Suggested fields:

```text
SolutionBlockNorm
AReal / AImag
PhiReal / PhiImag
EReal
EImag
JReal
JImag
JouleHeating or JoulePowerDensity
MaterialId / BodyId
SourceCurrent or ImpressedCurrent
```

If field names inconsistent, Cursor lists actual names in `summary.md`.

ParaView display order:

```text
1. Show three-body geometry: air / coil / plate;
2. Hide air or set transparent;
3. Show coil source current;
4. Show JouleHeating on plate;
5. Then show |J| or J_L2;
6. Optional SolutionBlockNorm.
```

Output directory:

```text
stage10_teacher_demo/case04_team7_three_body_scaffold/vtp/
stage10_teacher_demo/case04_team7_three_body_scaffold/stdout.txt
stage10_teacher_demo/case04_team7_three_body_scaffold/summary.md
```

---

## 9. Optional Case 5: Relaxed Sphere / Cylinder Particle Test

### 9.1 Required Before Afternoon?

Not required. If time insufficient, do not force.

But advisor previously questioned:

```text
Is code not based on particle discretization?
Is it only regular-lattice grid-like tests?
```

Therefore prepare relaxed geometry plan; even if afternoon only shows geometry and particle distribution.

### 9.2 Recommended Minimal Case

Add or prepare:

```text
test_3d_aphi_ck_inner_relaxed_sphere_mms
```

If sphere inconvenient, use cylinder:

```text
test_3d_aphi_ck_inner_relaxed_cylinder_mms
```

### 9.3 Purpose

Explain:

```text
Particles from level set geometry;
After particle relaxation;
Computation still uses SPH kernel neighbors;
Not dependent on regular Cartesian lattice.
```

### 9.4 Output

VTP at minimum:

```text
ParticlePosition
RelaxedParticleDistribution
BodyRegion
SolutionBlockNorm or analytic scalar field
KernelNeighborCount or particle_volume
```

If MMS not ready, at least output relaxed particles:

```text
stage10_teacher_demo/case05_relaxed_sphere_particles/vtp/
```

### 9.5 Follow-Up Formal Metrics

Add later:

```text
MMS error
GMRES residual
divA_relative
curlA/B diagnostic
hourglass mode response
```

Afternoon can say:

```text
We are preparing to extend current regular-box scaffold to relaxed sphere/cylinder particles to address particle discretization and kernel-based validation.
```

---

# Part III. Afternoon Demo Material Directory

## 10. Suggested Unified Directory for Cursor

Add:

```text
stage10_teacher_demo/
```

Structure:

```text
stage10_teacher_demo/
    README_FOR_TEACHER_DISCUSSION.md
    RUN_COMMANDS.md

    case01_contact_apply_vs_monolithic/
        stdout.txt
        summary.md

    case02_two_body_interface_mms/
        stdout.txt
        summary.md
        vtp/

    case03_phi_gauge_penalty/
        stdout.txt
        phi_penalty_summary.csv
        summary.md

    case04_team7_three_body_scaffold/
        stdout.txt
        summary.md
        vtp/

    case05_relaxed_sphere_particles_optional/
        stdout.txt
        summary.md
        vtp/
```

---

## 11. What `README_FOR_TEACHER_DISCUSSION.md` Should Include

Suggested one page.

### 11.1 Project Goal

```text
Build SPHinXsys particle/kernel-based matrix-free A--phi electromagnetic solver
for subsequent cold-crucible electromagnetic heating and multi-body Contact scenarios.
```

### 11.2 Currently Completed

```text
1. Inner / Contact pairwise operator;
2. Matrix-free fused apply;
3. GMRES solver;
4. phi gauge regularization;
5. Joule / E / J post-processing;
6. two-body interface MMS;
7. three-body TEAM7-like scaffold.
```

### 11.3 Current Demo Cases

List 4 cases:

```text
Case 1: operator equivalence
Case 2: two-body interface MMS
Case 3: phi gauge penalty
Case 4: three-body TEAM7-like Joule/J observable
```

### 11.4 Follow-Up Work

```text
1. divA / Coulomb gauge diagnostic;
2. curlA/B post-processing;
3. relaxed sphere/cylinder particle tests;
4. true TEAM7 reference;
5. cold-crucible four-body scaffold;
6. one-way thermal coupling.
```

---

## 12. Suggested `RUN_COMMANDS.md` Content

Cursor writes actual commands. Template:

```bash
cd build

ninja \
  test_3d_aphi_ck_contact_apply_vs_monolithic \
  test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured \
  test_3d_aphi_ck_contact_phi_gauge_penalty_sweep_diagnostic \
  test_3d_aphi_ck_gmres_contact_team7_three_body_scaffold

CUDA_VISIBLE_DEVICES=1 ./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_contact_apply_vs_monolithic/bin/test_3d_aphi_ck_contact_apply_vs_monolithic \
  > ../stage10_teacher_demo/case01_contact_apply_vs_monolithic/stdout.txt

CUDA_VISIBLE_DEVICES=1 ./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured/bin/test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured \
  > ../stage10_teacher_demo/case02_two_body_interface_mms/stdout.txt

CUDA_VISIBLE_DEVICES=1 ./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_contact_phi_gauge_penalty_sweep_diagnostic/bin/test_3d_aphi_ck_contact_phi_gauge_penalty_sweep_diagnostic \
  > ../stage10_teacher_demo/case03_phi_gauge_penalty/stdout.txt

CUDA_VISIBLE_DEVICES=1 ./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_gmres_contact_team7_three_body_scaffold/bin/test_3d_aphi_ck_gmres_contact_team7_three_body_scaffold \
  > ../stage10_teacher_demo/case04_team7_three_body_scaffold/stdout.txt
```

Note: if tests default output to `./output`, after each case copy:

```bash
cp -r output ../stage10_teacher_demo/caseXX_*/vtp
```

Or modify test to support `--output_dir`.

---

# Part IV. How to Present Each Case

## 13. Case 1: Operators Are Not Ad Hoc

Suggested talk:

```text
We first do not solve real EM problem; we verify matrix-free operator itself.
Same discrete operator apply on monolithic inner body and split Contact multi-body
should give consistent LHS.
This test checks pairwise kernel operator consistency.
```

Advisor may ask:

```text
Is this a matrix?
```

Answer:

```text
Not explicit matrix. We do not assemble sparse matrix.
GMRES only calls matrix-free apply each iteration, computing K(x) via particle neighbors and kernel gradients.
```

---

## 14. Case 2: Interface and GMRES Already Run

```text
This case is two-material interface MMS.
RHS from manufactured solution, solved with coupled GMRES.
We look at global residual, bodywise residual, and interface band residual.
```

Focus:

```text
Not only overall residual;
Also each body;
Also interface band;
Avoid hidden Contact interface errors.
```

---

## 15. Case 3: Gauge Issue Identified and Controlled

```text
A--phi equations have gauge DOF; phi can have constant drift.
We fix this DOF via phi_gauge_penalty.
At lambda_phi=0 residual can be good but raw phi/field drifts;
At lambda_phi=100 raw field stabilizes significantly.
```

Note:

```text
Not a physical material parameter; numerical gauge regularization.
```

Do not proactively expand:

```text
divA still open
A-divergence penalty not converging yet
```

Unless advisor asks.

---

## 16. Case 4: Can Already Output Physical Quantities

```text
This case is TEAM7-like three-body scaffold: air / coil / plate.
Impressed current in coil; induced E/J/Joule in plate.
We compare Contact multi-body vs monolithic scaffold plate Joule/J.
```

Focus:

```text
Engineering observables:
E field, current density, Joule heating.
```

Advisor may ask:

```text
Is this already real TEAM7?
```

Answer:

```text
Not full real TEAM7 benchmark yet.
TEAM7-like scaffold validates multi-body Contact solver and E/J/Joule chain.
Real TEAM7 still needs real geometry, external FEM/MFEM reference, resolution convergence.
```

---

# Part V. Hourglass / Zero-Energy Mode Answer Prep

## 17. If Advisor Asks About Hourglass / Zero-Energy Mode

Suggested answer:

```text
We have validated operator equivalence, MMS residual, Contact interface residual, and physical observables.
These show current operators did not expose obvious spurious modes in these tests.
But strictly, cannot claim no hourglass modes from these alone.
Next step: mode-response diagnostic:
checkerboard mode, random high-frequency mode, gradient mode, solenoidal mode.
```

Do not say:

```text
All operators guarantee no hourglass.
```

Say:

```text
We are extending from residual validation to near-null-mode / particle-disorder validation.
```

---

# Part VI. `B = curl A` Answer Prep

## 18. If Advisor Asks About Magnetic Field B

Answer:

```text
Current main outputs are A/phi, E, J, Joule.
B = curl A is next post-processing item.
We will use SPHinXsys kernel-based gradient of A for curl A,
not grid finite differences.
```

Formula:

```text
gradA[m,n] = ∂A_m / ∂x_n

B_x = gradA[z,y] - gradA[y,z]
B_y = gradA[x,z] - gradA[z,x]
B_z = gradA[y,x] - gradA[x,y]
```

Note:

```text
If SPHinXsys has vorticity/curl operator, can reference its kernel gradient convention;
but for A_real/A_imag and Contact multi-body, we will implement AphiComputeCurlACK.
```

---

# Part VII. Relaxed Particles Answer Prep

## 19. If Advisor Questions Regular Box Particles

Answer:

```text
Current box case is scaffold for manufactured solution and operator equivalence.
Not final validation.
We will add sphere/cylinder relaxed particle case:
particles from level set, particle relaxation, same A--phi operator / solver validation with kernel neighbors.
```

Suggested figure:

```text
relaxed sphere/cylinder particle distribution
```

If full MMS not ready before afternoon, at least prepare particle distribution VTP.

---

# Part VIII. Cursor Concrete Execution Checklist

## 20. Required Tasks

### Task 1: Create Demo Directory

```bash
mkdir -p stage10_teacher_demo
```

Create:

```text
README_FOR_TEACHER_DISCUSSION.md
RUN_COMMANDS.md
```

### Task 2: Generate File Map

Create:

```text
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/STAGE10_APHI_FILE_MAP_FOR_DISCUSSION.md
```

Per Part I of this document.

### Task 3: Prepare Four Demo Cases

Build and run:

```text
test_3d_aphi_ck_contact_apply_vs_monolithic
test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured
test_3d_aphi_ck_contact_phi_gauge_penalty_sweep_diagnostic
test_3d_aphi_ck_gmres_contact_team7_three_body_scaffold
```

Each case save:

```text
stdout.txt
summary.md
VTP output if available
```

### Task 4: TEAM7-Like Case Must Have Viewable VTP

At least output on plate:

```text
JouleHeating / JoulePowerDensity
JReal / JImag or CurrentDensity
EReal / EImag
SolutionBlockNorm
```

### Task 5: Optional Relaxed Particles

If time permits, add or run:

```text
test_3d_aphi_ck_inner_relaxed_sphere_particles
```

At least output:

```text
relaxed particle distribution VTP
```

MMS not required before afternoon.

---

## 21. Do Not Do For Now

```text
Do not large-scale merge production source;
Do not connect Contact A-divergence penalty;
Do not do cold-crucible four-body;
Do not do 10D thermal coupling;
Do not use A-divergence penalty current results as demo mainline;
Do not treat helper merge as pre-afternoon goal.
```

---

# Part IX. Afternoon Presentation Suggestions

## 22. Presentation Order

Suggested 15--20 minute flow:

```text
1. Background: why A--phi matrix-free SPH;
2. Operators: pairwise kernel operator + matrix-free apply;
3. Solver: GMRES without matrix assembly;
4. Interface: two-body MMS validation;
5. Contact: three-body TEAM7-like E/J/Joule output;
6. Gauge: phi gauge controlled;
7. Follow-up: divA/curlA/relaxed particles/real TEAM7.
```

## 23. One-Sentence Positioning

```text
We are presenting core-chain validation of A--phi SPH matrix-free solver,
not final real TEAM7 benchmark.
Current results show operator, GMRES, Contact, and E/J/Joule output chain largely usable;
next: Coulomb-gauge divA, B=curlA, and relaxed particle validation.
```

---

## 24. Final Decision

Most important before afternoon:

```text
Show less, show stable results;
Merge less, build file map first;
Discuss less blockers, discuss validation chain;
Must have VTP, especially TEAM7-like plate Joule/J;
Prepare to explain box is scaffold, relaxed particles is next step.
```
