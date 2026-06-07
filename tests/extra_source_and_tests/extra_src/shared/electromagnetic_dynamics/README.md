# A-phi Matrix-Free (CK / SYCL) Source File Guide

This directory contains the Stage 10 **matrix-free Laplace-type A-phi** solver and related test helper code developed on the SPHinXsys particle/kernel framework.

## Directory Structure (2026-05 Cleanup)

```
electromagnetic_dynamics/
├── README.md                          ← this file
├── all_electromagnetic_dynamics_ck.h  ← production operator/solver unified entry
│
├── (top level) production operators and solvers (Inner + Contact paired classes merged into same files)
│   aphi_coupling_modes_ck.h
│   aphi_matrix_free_operator_ck.*
│   aphi_laplace_ck.*                    ← Inner<> / Contact<> template parameter
│   aphi_grad_phi_coupling_ck.*
│   aphi_div_sigma_a_coupling_ck.*
│   aphi_phi_gauge_penalty_ck.*
│   aphi_a_divergence_penalty_ck.*
│   aphi_reaction_ck.*
│   aphi_joule_heating_ck.*
│   aphi_block_jacobi_preconditioner_ck.* ← includes AphiComputeBlockJacobiContactDynamicsBundle
│   aphi_gmres_solver_ck.*               ← includes AphiGMRESSolverCK + AphiGMRESContactSolverCK
│   aphi_gmres_workspace_ck.*
│   aphi_multibody_contact_gmres_ck.*
│   aphi_matrix_free_solve_ck.*          ← includes AphiMatrixFreeSolveCK + AphiMatrixFreeContactSolveCK
│   ...
│
├── alternate_krylov/                  ← Stage 8 off-main-line BiCGStab / PCG (not legacy)
├── diagnostics/                       ← debug pipelines and diagnostic helpers
├── test_helpers/                      ← test geometry/MMS/metrics (deletable when done)
└── benchmark/                         ← MMS / TEAM7 / cold crucible case definitions

2026-06 pre-commit: deleted 12 top-level one-line `#include` aliases; includes unified to canonical files (`aphi_gmres_solver_ck.*`, `aphi_matrix_free_solve_ck.*`, `diagnostics/aphi_assemble_lhs_debug_ck.*`, `aphi_block_jacobi_preconditioner_ck.*`; BiCGStab/PCG only in `alternate_krylov/`).
```

## Three-Layer Responsibilities

| Layer | Location | Description |
|---|---|---|
| **Production** | top-level `.h/.hpp` pairs | particle/kernel operators, matrix-free apply, Krylov/Contact solvers |
| **Diagnostics** | `diagnostics/` | grad/div debug, assemble debug, divA/A-gauge metrics, sweep helpers |
| **Test helpers** | `test_helpers/` | two-body/three-body geometry, MMS drivers, EM observables, device sync |
| **Benchmark** | `benchmark/` | manufactured solutions, TEAM7 canonical, cold crucible scaffold data |

## Related Documentation

- Doc index (3d_examples): `tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/README.md`
- Doc index (repo root): `docs/electromagnetic_aphi/README.md`
- File map: `tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/STAGE10_APHI_FILE_MAP_FOR_DISCUSSION.md`
- Test index: `tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/STAGE10_APHI_TEST_CASE_INDEX.md`
- TEAM7 native record: `tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_17_TEAM7_NATIVE_RECORD.md`
- MR plan: `tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_TEAM7_SMALL_AIR_MULTIRES_PLAN.md`
- Pre-commit cleanup record: `tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_17_19_PRECOMMIT_CLEANUP_RECORD.md`

## Conventions

- **OPERATOR** files: participate in fused matrix-free apply or Contact GMRES main path.
- **DIAGNOSTIC** files: debug assemble or post-solve metrics only; do not change production default behavior.
- Before modifying production operators, run `test_3d_aphi_ck_fused_apply_vs_debug_lhs` and related Contact equivalence tests.
