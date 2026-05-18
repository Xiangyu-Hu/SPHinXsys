# EM A-phi Electromagnetic Heating Archive

## Purpose

This file is the single archival note for the electromagnetic-heating work under
`tests/extra_source_and_tests`. It replaces the temporary Chinese notes,
presentation drafts, scan scripts, and ad hoc experiment records that were used
during prototyping.

The goal of the retained tree is simple:

- keep the code and tests that are needed to build and continue the work;
- keep the minimum geometry assets required by the TEAM7-related examples;
- keep one English summary that explains what was implemented, what was
  validated, and what still needs follow-up.

## What Is Intentionally Kept

### Shared implementation

The following shared files are part of the retained A-phi electromagnetic
heating implementation:

- `extra_src/shared/legacy_aphi_archive/baselines/electromagnetic_aphi_global_implicit_solver.*`
- `extra_src/shared/legacy_aphi_archive/baselines/electromagnetic_aphi_laplace_eigen.*`
- `extra_src/shared/legacy_aphi_archive/baselines/electromagnetic_aphi_operator_comparison_eigen.*`
- `extra_src/shared/legacy_aphi_archive/ck_operators/electromagnetic_component_hessian_ck.*`
- `extra_src/shared/aphi_case_support/electromagnetic_multiturn_coil_drive.*`
- `extra_src/shared/aphi_case_support/electromagnetic_team7_aphi_dynamics.*`
- `extra_src/shared/aphi_case_support/electromagnetic_team7_aphi_frequency_dynamics.*`
- `extra_src/shared/aphi_case_support/em_aphi_residual_diagnostics_snippet.hpp`
- `extra_src/shared/em_matched_stencil_body_relations.h`

These files provide:

- the original TEAM7-oriented A-phi dynamics path;
- the frequency-domain extensions;
- the Laplace-structured global complex prototype;
- operator-comparison helpers;
- source-current and diagnostics utilities used by the retained examples.

### Retained examples and tests

The retained EM-heating examples fall into three groups.

#### 1. Operator and field-chain validation

- `test_3d_em_dual_body_linear_a`
- `test_3d_em_frequency_dual_body_linear_a`
- `test_3d_em_phi_gradient_linear`
- `test_3d_em_frequency_phi_only_linear`
- `test_3d_em_frequency_joule_chain`
- `test_3d_em_frequency_a_solver_manufactured`
- `test_3d_em_frequency_a_solver_nonzero_curlcurl_manufactured`
- `test_3d_em_frequency_nonzero_curlcurl_manufactured`

These are the compact validation cases used to isolate:

- `A -> curl(A) -> curl(nu B)`;
- `phi -> grad(phi)`;
- `A, phi -> E -> J -> Joule heat`;
- manufactured frequency-domain solver behavior.

#### 2. Laplace-structured global prototype cases

- `test_3d_em_aphi_laplace_eigen_manufactured`
- `test_3d_em_aphi_operator_mode_comparison`
- `test_3d_em_aphi_laplace_eigen_physical_source`
- `test_3d_em_aphi_laplace_eigen_team7_like`
- `test_3d_em_aphi_laplace_eigen_team7_stl_em_only`

These cases were created as an independent comparison track instead of another
local patch to the strong-form curl-curl route. They are the main retained
prototype set for the global complex Laplace-structured A-phi solver.

#### 3. Main TEAM7 scaffold

- `particle_generation_em/particle_generation_em.cpp`
- `particle_generation_em/CMakeLists.txt`
- `particle_generation_em/data/coil.stl`
- `particle_generation_em/data/plate.stl`

This remains the main geometry-rich TEAM7 scaffold for the original
particle-based route.

## Validation Ladder Summary

The retained work established a layered validation ladder instead of jumping
directly from code changes to TEAM7 curve fitting.

### Step 1: Magnetic operator sanity

Use:

- `test_3d_em_dual_body_linear_a`
- `test_3d_em_frequency_dual_body_linear_a`

Purpose:

- verify discrete `curl(A)`;
- verify discrete `curl(nu curl(A))`;
- compare merged-body and split-body behavior across an internal interface.

### Step 2: Scalar-potential path sanity

Use:

- `test_3d_em_phi_gradient_linear`
- `test_3d_em_frequency_phi_only_linear`

Purpose:

- verify `grad(phi)` on its own;
- verify the frequency-domain scalar-potential source and relaxation path before
  mixing it back into the full coupled solve.

### Step 3: Frequency-domain post-processing sanity

Use:

- `test_3d_em_frequency_joule_chain`

Purpose:

- verify sign conventions and reconstruction for `E`, `J`, and Joule heat when
  `A` and `phi` are prescribed.

### Step 4: Global-complex prototype sanity

Use:

- `test_3d_em_aphi_laplace_eigen_manufactured`
- `test_3d_em_aphi_operator_mode_comparison`

Purpose:

- validate assembly of the global complex sparse system;
- compare Laplace-structured and weak-inspired curl-curl operator responses;
- study gauge, scaling, equilibration, and mode sensitivity in a controlled
  setting.

### Step 5: Geometry-carrying cases

Use:

- `test_3d_em_aphi_laplace_eigen_physical_source`
- `test_3d_em_aphi_laplace_eigen_team7_like`
- `test_3d_em_aphi_laplace_eigen_team7_stl_em_only`
- `particle_generation_em`

Purpose:

- move from manufactured/operator tests toward realistic geometry, current
  driving, material jumps, and TEAM7-like setup decisions.

## Main Technical Conclusions So Far

### Original TEAM7-oriented frequency A-phi route

The retained investigation supports the following summary:

- the post-processing chain `A, phi -> E -> J -> Joule heat` is not the main
  unresolved issue;
- the harder problem is still the quality of the solved electromagnetic field,
  especially the magnetic operator path and coil-air transmission in the
  geometry-rich setting;
- previous debugging identified important implementation-level fixes, including
  operator execution order and block-sweep timestep normalization;
- after those fixes, stability improved substantially, but reference-quality
  accuracy for TEAM7 remains unfinished.

### Laplace-structured global prototype

The retained investigation supports the following summary:

- the global complex Laplace-structured A-phi prototype is viable as an
  independent comparison path;
- at validation-friendly scales, the manufactured case recovers the solution
  with small error and acceptable linear residuals;
- at larger `sigma * frequency`, raw performance degrades mainly because of
  matrix scaling and solver-conditioning issues rather than immediate evidence
  that the formulation itself is invalid;
- simple diagonal equilibration materially improves difficult manufactured
  stress cases;
- material jumps, especially in `nu`, remain a priority for continued
  validation;
- the STL TEAM7 EM-only case is already useful for structure and scaling checks,
  but not yet a trusted physical benchmark.

## Practical Notes For Future Work

### Minimal retained entry points

The most useful entry points after cleanup are:

- `particle_generation_em` for the original TEAM7 scaffold;
- `test_3d_em_aphi_laplace_eigen_manufactured` for controlled global-prototype
  verification;
- `test_3d_em_aphi_operator_mode_comparison` for operator-energy comparisons;
- `test_3d_em_aphi_laplace_eigen_team7_stl_em_only` for the small STL-based
  three-body EM-only prototype;
- the linear/made-up validation cases for isolating operator behavior.

### What was removed on purpose

The cleanup removes items that were useful during iteration but should not be
treated as durable source assets:

- temporary Chinese progress notes and status reports;
- slide drafts and presentation artifacts;
- scan scripts used for one-off sweeps;
- temporary logs;
- regenerated or heavyweight input references not required to build the retained
  cases.

If a deleted sweep becomes important again, it should be recreated as a focused,
documented workflow rather than restored as a pile of exploratory scripts.

## Current Recommendation

Treat the retained tree as a compact GitHub checkpoint:

- code and cases are preserved;
- geometry assets required by the STL/TEAM7 examples are preserved;
- exploratory notes are condensed into this file;
- the next contributor should rebuild any new sweep script only when the exact
  validation question is known.
