# CURSOR — Stage 10.13 A-phi Electromagnetic Solver Next-Step Work Instructions (Detailed)

> For Cursor Auto mode. Execute strictly per this document; do not expand scope on your own.  
> Repository root assumed to be `/home/yongchuan/sphinxsys`, branch assumed to be `feature/electromagnetic`.  
> All new/modified content continues in the `tests/extra_source_and_tests` system; do not change the old aphi route; do not break SPHinXsys SYCL CK style.

---

## 0. Current stage conclusion — read first

Stage 10.11 / 10.12 conclusions are largely frozen:

1. **A-phi main solve chain is connected**: matrix-free GMRES, A/phi, E, J, Joule manufactured-solution validation is basically usable.
2. **Inner region core divA is already good**: global divA anomaly mainly from boundary/support deficiency, not core pairwise operator error.
3. **Contact two-body A-divergence penalty research gate is closed**: effective approach is:
   - main operator `K(A, phi)` still uses `Inner + Contact`;
   - but A-divergence penalty stencil must use `InnerOnly`;
   - penalty `divA/grad(divA)` must not use Contact neighbors across split-body interface.
4. **production default still must not enable lambda_A**:
   - `production lambda_A = off`;
   - `eta_A = 0.1` only as research primary;
   - when `eta_A >= 0.2`, must disable `polish_sweeps`, otherwise block-GS polish destroys the solution.
5. **projection still deferred**: do not implement post-solve projection. No clear definition yet of `A_new = A - grad chi` and `phi_new = phi + i omega chi` synchronized boundary conditions under Contact/multibody.
6. **B/H are now computable and displayable, but not as high-precision production hard gate**:
   - B=curlA postprocessing must use B-corrected / LinearCorrectionMatrix + LinearGradient route;
   - uncorrected pairwise curl only for comparison diagnostic, not pass/fail criterion;
   - Az2D/CrossSine B error ~3.5% @ dp=0.1, refining to dp=0.075 ~2%, indicating discrete curl accuracy issue.
7. **cold crucible showcase case may be built**, but must be explicitly labeled `demonstration / pipeline smoke test / not quantitative validation`; do not write that physical validation is complete.

---

## 1. Absolute prohibitions

Cursor must not:

1. Change `lambda_A` to production default on.
2. Implement post-solve projection.
3. Change Contact A-penalty default stencil back to `InnerContact`.
4. Use `curlH ~= J` as current div-free MMS hard gate.
5. Use `uncorrected pairwise curl` as B=curlA hard gate.
6. Write cold crucible demo as quantitative validation.
7. Refactor GMRES main solver.
8. Change SPHinXsys core framework; only advance in `tests/extra_source_and_tests` and corresponding electromagnetic helpers.
9. Bypass CK/SYCL style with CPU-only special cases, unless host-side diagnostic statistics only.
10. Introduce new third-party libraries.

---

## 2. Concept clarification: why solve uses pairwise Laplace but B=curlA postprocessing uses B-corrected curl

### 2.1 Do not confuse two goals

Pairwise Laplace in the solver is the PDE operator for constructing the matrix-free equation:

```text
K(A, phi) = RHS
```

It emphasizes:

- pairwise symmetry;
- conservation / consistency with SPH Laplacian formulation;
- consistency with matrix-free GMRES/PC discrete system;
- suitable for `Laplace A`, `grad phi`, `div(sigma A)`, etc.

B=curlA is a postprocessing observable:

```text
B = curl A
```

It emphasizes first-derivative reconstruction accuracy, especially first-order consistency for exact linear vector fields. Direct uncorrected pairwise gradient for curl lacks local kernel moment linear correction; in current tests, Linear2D fields that should be easy still show ~O(1) relative error. Therefore it cannot be B hard gate.

### 2.2 "B" in B-corrected is not magnetic field B

`B-corrected` here is SPHinXsys kernel correction / linear correction matrix semantics, not magnetic field `B`.

Current helper postprocessing route should remain:

```cpp
InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>>> linear_correction_matrix(inner);
InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd>>> field_gradient(inner, field_name);
StateDynamics<MainExecutionPolicy, AphiVectorGradientCurlCK> curl_kernel(body, gradient_name, curl_name);
```

That is:

1. compute `LinearCorrectionMatrix` first;
2. then `LinearGradient<Inner<Vecd>>` for corrected vector field gradient;
3. finally `AphiVectorGradientCurlCK` takes curl from gradient tensor.

This is essentially SPHinXsys corrected gradient/linear gradient approach; A-phi reuses it for B=curlA diagnostic/observable postprocessing.

### 2.3 Standard wording for this stage

Use externally or in documentation:

```text
The matrix-free A-phi solve still uses the pairwise Laplace-form operator. The magnetic field B=curlA is an observable reconstructed from the solved vector potential. For this first-derivative reconstruction, the corrected SPHinXsys LinearGradient route is used as the diagnostic/reference path. The uncorrected pairwise curl is retained only as a diagnostic comparison and must not be used as a hard accuracy gate.
```

---

## 3. Concept clarification: current Contact test is two bodies, not one body with two regions

Current Contact test should be understood as: **two independent SPHBody / SolidBody, finding each other's particles via Contact relation**.

Not two regions defined in the same body.

Typical structure:

```cpp
SolidBody left_body(sph_system, makeShared<AphiHalfSpaceBoxShape>("LeftBody", left_center, halfsize));
SolidBody right_body(sph_system, makeShared<AphiHalfSpaceBoxShape>("RightBody", right_center, halfsize));

Inner<> left_inner(left_body);
Inner<> right_inner(right_body);

Contact<> left_to_right(left_body, {&right_body});
Contact<> right_to_left(right_body, {&left_body});
```

Then in apply/GMRES, execute separately for left/right body:

```cpp
AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_left(
    left_body, left_inner, left_to_right, ...);

AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_right(
    right_body, right_inner, right_to_left, ...);
```

So it truly tests SPHinXsys Contact relation:

- left body finds right body particles via `left_to_right`;
- right body finds left body particles via `right_to_left`;
- each body has its own particles, variables, materials, inner relation and contact relation.

This is the meaning of the Contact route. If same body with two regions, only material field discontinuity can be tested, not multibody neighbor handling under SPHinXsys Contact relation.

---

## 4. Concept clarification: what full Maxwell MMS is — not solving full Maxwell equations

### 4.1 Meaning of MMS

MMS = Method of Manufactured Solutions.

It is not a new physical model, nor must it solve full Maxwell equations. The approach is to artificially construct an analytic solution, e.g.:

```text
A_exact(x,y,z), phi_exact(x,y,z)
```

substitute into current discrete/continuous equation, derive RHS. Then let the numerical solver solve this RHS, check recovery of `A_exact, phi_exact`, and further check derived quantities E/J/Joule/B/H, etc.

### 4.2 Why full Maxwell MMS is needed

Current div-free MMS mainly guarantees:

```text
div A = 0
A, phi -> E
E -> J
J -> Joule
B = curl A
H = nu B
```

But it does not necessarily guarantee:

```text
curl H = J
```

Because `curlH = J` is the Maxwell-Ampere relation; it must be satisfied together by constructed A/phi/material/RHS. Ordinary div-free A field does not necessarily satisfy this constraint.

So full Maxwell MMS does not mean writing a complete Maxwell solver, but constructing manufactured fields that simultaneously satisfy or explicitly test:

```text
B_exact = curl A_exact
H_exact = nu B_exact
J_exact = curl H_exact
E_exact = J_exact / sigma
E_exact = -i omega A_exact - grad phi_exact
```

If such a self-consistent analytic field set can be constructed, `curlH≈J` can become hard gate. Otherwise current P5 `curlH≈J` can only be diagnostic/informational.

### 4.3 Is full Maxwell MMS required now?

Not highest priority currently.

Suggested order:

1. cold crucible demo/pipeline smoke test first;
2. then B=curlA dp refinement;
3. then full Maxwell MMS, making `curlH≈J` stricter validation.

---

## 5. Current work overview for Cursor

Execute in order; do not skip steps.

| Priority | Work package | Goal | Required |
|---|---|---|---|
| P0 | Documentation freeze and test gate cleanup | Avoid Auto Cursor misjudgment later | Required |
| P1 | cold crucible demonstration case | Show staged pipeline to client | Required |
| P2 | demo output fields and VTP/VTU visualization | Displayable A/phi/B/E/J/Joule | Required |
| P3 | B=curlA dp refinement enhancement | Check if B decreases with resolution | Recommended |
| P4 | Contact InnerOnly PC consistency supplemental test | Prevent PC being changed wrong later | Recommended |
| P5 | full Maxwell MMS design document | Not implemented yet, design only | Optional |
| P6 | three-body Contact A-penalty | Deferred, do not do | Forbidden now |
| P7 | projection | Deferred, do not do | Forbidden now |

---

# P0 — Documentation freeze and test gate cleanup

## P0.1 Add master handoff document

New file:

```text
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_13_CHATGPT_HANDOFF.md
```

Content must include:

```text
production lambda_A = off
research eta_A = 0.1
Contact A-penalty stencil = InnerOnly
main Contact operator = Inner + Contact
eta_A >= 0.2 => polish_sweeps = 0
projection = deferred
boundary divA = open diagnostic / support deficiency
B=curlA = use corrected LinearGradient route, uncorrected pairwise curl only diagnostic
E/J/Joule = validated in MMS, suitable for demo pipeline
cold crucible demo = pipeline smoke test, not quantitative validation
```

## P0.2 Update misleading status in old records

Check and update:

```text
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_11_BOUNDARY_VALIDATION_RECORD.md
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_12_B_H_OBSERVABLE_RECORD.md
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_12_CONTACT_A_PENALTY_RECORD.md
```

Ensure no misleading wording such as:

```text
Contact eta>0 open
```

Replace with more accurate:

```text
Contact two-body A-penalty research gate closed for eta=0 and eta=0.1 with InnerOnly stencil; eta=0.2 is usable with polish_sweeps=0; three-body A-penalty remains deferred.
```

## P0.3 Strengthen C4/C5 test gates

Check tests:

```text
tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_contact_a_divergence_penalty_inner_only_stencil_diagnostic/
tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_contact_a_divergence_penalty_upper_eta_solver_diagnostic/
```

Requirements:

1. C4 must hard-check:
   - `eta0_regression_ok == true`
   - `research_improved == true`
   - `InnerOnly eta=0.1 mms_global_rel_l2 <= 1e-6` or use current test threshold;
   - `InnerContact eta=0.1` should be clearly worse than `InnerOnly eta=0.1`, ratio at least `1e3`.
2. C5 must hard-check:
   - `eta=0.2, polish_sweeps=0` `global_true_rel` finite;
   - `global_true_rel <= 1e-4`, if current record stable at `~7.8e-6`, may use stricter `1e-4`;
   - `eta=0.2, polish_sweeps=1` diagnostic only, do not require fail, but if fail print explanation.

Reproduction commands:

```bash
cd /home/yongchuan/sphinxsys/build
ninja test_3d_aphi_ck_contact_a_divergence_penalty_inner_only_stencil_diagnostic \
      test_3d_aphi_ck_contact_a_divergence_penalty_upper_eta_solver_diagnostic

./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_contact_a_divergence_penalty_inner_only_stencil_diagnostic/bin/test_3d_aphi_ck_contact_a_divergence_penalty_inner_only_stencil_diagnostic

./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_contact_a_divergence_penalty_upper_eta_solver_diagnostic/bin/test_3d_aphi_ck_contact_a_divergence_penalty_upper_eta_solver_diagnostic
```

Acceptance: both tests output `passed=1`.

---

# P1 — Add cold crucible demonstration case

## P1.1 Goal

Add a "cold crucible shape showcase case". Goal is not accurate electromagnetic field validation, but to show the client the current A-phi solve chain can run:

```text
geometry/material/frequency/source setup
-> A-phi solve or impressed-field assignment
-> compute E
-> compute J
-> compute Joule heat
-> compute B=curlA diagnostic
-> compute H=nuB diagnostic
-> write VTP/VTU output for ParaView
```

Must clearly write in README and terminal output:

```text
This is a demonstration / smoke-test case, not a quantitative validation case.
```

## P1.2 New directory

New directory:

```text
tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_cold_crucible_demo/
```

At minimum include:

```text
CMakeLists.txt
test_3d_aphi_ck_cold_crucible_demo.cpp
README.md
```

Optional new helper:

```text
tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_cold_crucible_demo_helpers.h
```

If helper content is large, must be in helper; do not put all logic in cpp.

## P1.3 Geometry design

Start with simple 3D axis-aligned / cylindrical approximation; do not pursue real cold crucible slit geometry.

Suggested geometry:

1. `AirBody`: outer air domain, a box or cylinder-like body.
2. `CrucibleWallBody`: vessel/crucible wall, conductor or low-conductor material acceptable.
3. `MeltBody`: inner melt region, conductive fluid/conductor material.
4. `CoilBody` or impressed source region: outer coil region. May skip real coil particles initially, use analytic impressed A or impressed J source.

Minimum viable version may start with only:

```text
AirBody + MeltBody + analytic impressed coil source
```

Then second step add:

```text
CrucibleWallBody
```

Do not start with complex multibody contact; easy to get stuck on geometry and boundary.

## P1.4 Recommended MVP route

For fastest client showcase, recommend `impressed analytic A` version first, not immediate real outer coil solve.

### MVP-A: no full coil solve, directly give analytic impressed A

Construct vector potential around z-axis, e.g.:

```text
A = A_theta(r,z) e_theta
```

In Cartesian coordinates:

```text
A_x = -A_theta(r,z) * y / (r + eps)
A_y =  A_theta(r,z) * x / (r + eps)
A_z = 0
```

where:

```text
r = sqrt(x^2 + y^2)
A_theta(r,z) = A0 * exp(-((r - r_coil)^2 / w_r^2 + z^2 / w_z^2))
```

Given frequency `omega`, if `phi=0`:

```text
E = -i omega A
J = sigma E
Joule = 0.5 * sigma * |E|^2
B = curl A
H = nu B
```

This is enough to show:

- A/B in air;
- E/J/Joule in conductive melt;
- clear J/Joule difference under different material sigma;
- ParaView can display magnetic field, current and Joule heat.

Note: this route demonstrates postprocessing pipeline, not complete coil-fluid coupling solve.

### MVP-B: then connect A-phi solve

If existing A-phi solver RHS setup runs stably, add solve mode after MVP-A:

```text
--mode solve
```

Flow:

```text
construct RHS from manufactured/impressed source
run GMRES A-phi solve
compute E/J/Joule/B/H
write output
```

If solve mode unstable, do not block MVP-A. MVP-A must complete first.

## P1.5 Material parameter suggestions

Use numerically stable nondimensional parameters first; do not use real physical properties directly.

Suggested:

```text
omega = 1.0e2 or 1.0e3
sigma_air = 0.0 or very small, e.g. 1e-8
sigma_melt = 1.0
sigma_wall = 0.1 or 1.0
nu_air = 1.0
nu_melt = 1.0
nu_wall = 1.0
A0 = 1.0e-3 or adjust by output magnitude
```

Output must print:

```text
max|A|, max|B|, max|E|, max|J|, max Joule
integral Joule over MeltBody
integral Joule over AirBody
```

Acceptance criteria:

1. All max values finite;
2. MeltBody `max|J| > 0`;
3. MeltBody `integral Joule > 0`;
4. AirBody if `sigma_air=0`, then `J/Joule` should be near 0;
5. Program outputs `passed=1`.

## P1.6 Contact relation design

If demo uses multiple bodies, use real Contact relation; do not fake Contact with two regions in same body.

Minimum Contact structure:

```cpp
Inner<> air_inner(air_body);
Inner<> melt_inner(melt_body);

Contact<> air_contact(air_body, {&melt_body});
Contact<> melt_contact(melt_body, {&air_body});
```

If adding wall:

```cpp
Contact<> air_contact(air_body, {&wall_body, &melt_body});
Contact<> wall_contact(wall_body, {&air_body, &melt_body});
Contact<> melt_contact(melt_body, {&wall_body, &air_body});
```

Note: demo MVP-A if only analytic impressed field + local postprocess, may not depend on Contact apply temporarily; but solve mode must use Contact apply.

## P1.7 Output fields

At minimum output these fields, names consistent with existing Aphi variable names:

```text
A_real, A_imag
Phi_real, Phi_imag
B_real, B_imag
H_real, H_imag
E_real, E_imag
J_real, J_imag
JouleHeatSource
Sigma
Nu
BodyId or MaterialId
```

If existing writer only outputs SPHinXsys state variables, ensure these variables are registered as state variables and synced to device/host correctly.

## P1.8 CMake integration

After adding test directory, update parent CMake so:

```bash
cd /home/yongchuan/sphinxsys/build
ninja test_3d_aphi_ck_cold_crucible_demo
./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_cold_crucible_demo/bin/test_3d_aphi_ck_cold_crucible_demo
```

Acceptance:

```text
passed=1
```

And ParaView-readable files found in output directory.

## P1.9 README must be clear

`test_3d_aphi_ck_cold_crucible_demo/README.md` must include:

```text
Purpose
-------
This case demonstrates the current A-phi electromagnetic post-processing pipeline for a cold-crucible-like geometry. It is not a quantitative validation case.

What is demonstrated
--------------------
- analytic/impressed vector potential A around a crucible-like domain
- B=curlA diagnostic using corrected LinearGradient route
- frequency-domain electric field E = -i omega A - grad phi
- induced current J = sigma E in conductive melt/wall regions
- Joule heat source q = 0.5 sigma |E|^2
- ParaView output for A/B/E/J/Joule

What is not demonstrated yet
----------------------------
- fully validated coil-current electromagnetic solve
- real cold crucible slit geometry
- wall dummy/support-complete boundary treatment
- quantitative TEAM7/cold-crucible benchmark accuracy
- fluid-thermal two-way coupling
```

---

# P2 — Demo visualization output enhancement

## P2.1 Output ParaView-friendly scalar/vector

Besides complex fields, suggest additional magnitude scalars:

```text
A_magnitude
B_magnitude
E_magnitude
J_magnitude
JouleHeatSource
```

So client showcase can view color maps directly without ParaView Calculator.

## P2.2 Recommended ParaView visualization

Add to README:

```text
Suggested ParaView visualization:
1. Load all body VTP/VTU outputs.
2. Color MeltBody by JouleHeatSource.
3. Add Glyph for J_real or J_magnitude if vector output is available.
4. Color AirBody by B_magnitude and set opacity to 0.2-0.4.
5. Show CrucibleWallBody as semi-transparent gray if present.
```

## P2.3 Acceptance

After demo run must print similar:

```text
test_3d_aphi_ck_cold_crucible_demo max_A=... max_B=... max_E=... max_J=... max_Joule=... melt_Joule_integral=... passed=1
```

---

# P3 — B=curlA dp refinement enhancement

## P3.1 Goal

Currently known:

```text
Az2D/CrossSine: B_err ~3.5% at dp=0.1
Az2D: B_err ~2% at dp=0.075
```

Next add:

```text
dp = 0.05
```

At least for Az2D. If too slow, Az2D only, skip CrossSine3D.

## P3.2 Files

Modify:

```text
tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_curl_b_dp_refinement_diagnostic/
tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_curl_b_dp_refinement_diagnostic_helpers.h
```

## P3.3 Acceptance criteria

Do not require B_err to 1e-6. Current goal is only to confirm convergence with dp.

Suggested hard gate:

```text
B_err(dp=0.075) < B_err(dp=0.10)
B_err(dp=0.05) <= B_err(dp=0.075) * 1.2
B_err(dp=0.05) < 0.02  // if slightly higher in practice, may change to informational, but must record reason
```

If dp=0.05 too slow or insufficient memory, may be optional diagnostic only, but must note in record.

---

# P4 — Contact InnerOnly graddiv PC consistency supplemental test

## P4.1 Goal

Current C2 tested Contact graddiv PC consistency, but focused on InnerContact pipeline. Now research default is InnerOnly penalty stencil, so suggest supplemental test confirming:

```text
Contact main operator: Inner + Contact
A-penalty stencil: InnerOnly
PC graddiv block: only Inner contribution
```

## P4.2 New test directory

Add:

```text
tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_contact_inner_only_graddiv_pc_consistency_diagnostic/
```

May reuse helper:

```text
tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_contact_graddiv_pc_diagnostic_helpers.h
```

If new helper needed, name:

```text
tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_contact_inner_only_graddiv_pc_consistency_helpers.h
```

## P4.3 Check content

1. Construct left/right two-body Contact case.
2. Set options:

```cpp
options.use_a_divergence_penalty = true;
options.a_divergence_penalty = lambda_a;
options.contact_a_divergence_penalty_stencil = AphiContactADivergencePenaltyStencilMode::InnerOnly;
```

3. Finite-difference or direct apply comparison on probe vector.
4. Verify PC diagonal does not include Contact graddiv contribution.
5. Core region error should be near original Inner PC magnitude.

## P4.4 Acceptance

Output:

```text
inner_only_pc_core_max_diff=...
contact_graddiv_skipped=1
passed=1
```

---

# P5 — full Maxwell MMS design document, not implemented yet

## P5.1 New design document

Add:

```text
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_13_FULL_MAXWELL_MMS_DESIGN.md
```

Design only, no implementation.

## P5.2 Must explain

Document must clearly state:

```text
MMS = Method of Manufactured Solutions.
It is not a new full Maxwell solver.
The purpose is to construct exact A, phi, B, H, E, J fields that are mutually consistent so that curlH≈J can be tested as a hard gate.
```

## P5.3 Suggested formulas

List relationships to satisfy:

```text
B_exact = curl A_exact
H_exact = nu B_exact
J_exact = curl H_exact
E_exact = J_exact / sigma
E_exact = -i omega A_exact - grad phi_exact
Joule_exact = 0.5 sigma |E_exact|^2
```

If all relationships cannot be satisfied simultaneously, must clarify source term or treat some as diagnostic only.

---

# P6 — three-body Contact A-penalty: deferred

Do not do three-body A-penalty currently.

Enter three-body only when all conditions met:

1. P0 complete;
2. P1 cold crucible demo at least MVP-A complete;
3. P4 InnerOnly PC consistency complete;
4. ChatGPT/user explicitly requests three-body gate.

Fixed configuration when starting three-body:

```text
main K = Inner + Contact
A-penalty stencil = InnerOnly
eta_A = 0.1 first
eta_A = 0.2 optional with polish_sweeps=0
production lambda_A = off
```

---

# P7 — projection: continue deferred

Do not implement currently:

```text
A_new = A - grad chi
phi_new = phi + i omega chi
```

Reasons:

1. chi interface condition under Contact/multibody undefined;
2. projection may destroy verified E/J/Joule;
3. boundary divA mainly support deficiency; projection may inject boundary artifact into full field;
4. no wall dummy/real TEAM7 boundary to support projection validation currently.

Only retain or extend design documents; no production code.

---

## 8. Final regression test checklist

After P0/P1 complete, run at least:

```bash
cd /home/yongchuan/sphinxsys/build

ninja test_3d_aphi_ck_contact_a_divergence_penalty_pairwise_operator_diagnostic \
      test_3d_aphi_ck_contact_graddiv_pc_consistency_diagnostic \
      test_3d_aphi_ck_contact_a_divergence_penalty_two_body_mms \
      test_3d_aphi_ck_contact_a_divergence_penalty_inner_only_stencil_diagnostic \
      test_3d_aphi_ck_contact_a_divergence_penalty_upper_eta_solver_diagnostic \
      test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms \
      test_3d_aphi_ck_curl_b_dp_refinement_diagnostic \
      test_3d_aphi_ck_h_nu_b_curl_h_j_observable_diagnostic \
      test_3d_aphi_ck_cold_crucible_demo
```

Run individually:

```bash
./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_contact_a_divergence_penalty_pairwise_operator_diagnostic/bin/test_3d_aphi_ck_contact_a_divergence_penalty_pairwise_operator_diagnostic
./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_contact_graddiv_pc_consistency_diagnostic/bin/test_3d_aphi_ck_contact_graddiv_pc_consistency_diagnostic
./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_contact_a_divergence_penalty_two_body_mms/bin/test_3d_aphi_ck_contact_a_divergence_penalty_two_body_mms
./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_contact_a_divergence_penalty_inner_only_stencil_diagnostic/bin/test_3d_aphi_ck_contact_a_divergence_penalty_inner_only_stencil_diagnostic
./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_contact_a_divergence_penalty_upper_eta_solver_diagnostic/bin/test_3d_aphi_ck_contact_a_divergence_penalty_upper_eta_solver_diagnostic
./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms/bin/test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms
./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_curl_b_dp_refinement_diagnostic/bin/test_3d_aphi_ck_curl_b_dp_refinement_diagnostic
./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_h_nu_b_curl_h_j_observable_diagnostic/bin/test_3d_aphi_ck_h_nu_b_curl_h_j_observable_diagnostic
./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_cold_crucible_demo/bin/test_3d_aphi_ck_cold_crucible_demo
```

Acceptance:

```text
all existing gates still passed=1
cold crucible demo outputs passed=1
cold crucible demo outputs ParaView-readable files
production lambda_A not enabled
projection not implemented
```

---

## 9. Records to write after completion

Add or update:

```text
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_13_COLD_CRUCIBLE_DEMO_RECORD.md
```

Must include:

1. Modified file list;
2. New test list;
3. Demo geometry description;
4. Material parameters;
5. source/impressed A formulas;
6. Output fields;
7. Terminal passed=1 summary;
8. ParaView output paths;
9. Explicit statement: not quantitative validation;
10. Next-step suggestions.

---

## 10. Final status summary template for user/ChatGPT

After completion, provide in record:

```text
Stage 10.13 completed the cold-crucible-like A-phi electromagnetic demonstration case. The case uses analytic/impressed vector potential to represent outer coil excitation, outputs A/phi/B/H/E/J/Joule, and can display air-domain magnetic field, conductive melt induced current and Joule heat distribution in ParaView. The case is a pipeline smoke test, not quantitative validation. Existing MMS gates not broken; Contact A-penalty remains research-only, production lambda_A off; projection continues deferred.
```
