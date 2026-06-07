# Stage 10.17–10.19 Pre-Commit Cleanup Record

> Date: 2026-06-05  
> Branch: `feature/electromagnetic`  
> Plan: `STAGE10_17_PRECOMMIT_CLEANUP_AND_CLI_FIX_PLAN.md` (see docs/electromagnetic_aphi/)

## Objectives

Before the milestone commit: fix known CLI bugs, unify reload/metadata, clean up legacy-route files, **without changing** the main numerical path of matrix-free Laplace A-φ + SYCL CK + Contact GMRES + TEAM7 native.

## A. CLI / Metadata Fixes

### Bug: `--team7-case=` off-by-one

**File:** `extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_team7_native_geometry_config.h`

- Original `arg.substr(12)` for `--team7-case=` (length 13) passed `=uniform-legacy-dp6`, causing preset to fail silently.
- All `--team7-*` prefix parsing changed to `team7NativeCliSuffix(arg, prefix)`.

### New / Improved

| Item | Description |
|----|------|
| `multires-small-dp6-l1` preset | Verified MR levels=1 case |
| `uniform-small-dp6` preset | small air + uniform dp6 |
| `team7_case_preset` field | Recorded after CLI successfully applies preset |
| `team7NativeAuditReloadMetadata` | `exit(1)` if `--team7-reload-case` and `reload_case_id` in `team7_native_geometry.txt` disagree |
| `team7NativeReloadCaseLookupOrder` | Alias staging such as `uniform_legacy_dp6` ↔ `team7_native_uniform_legacy_dp6_l0` |
| `printTeam7NativeGeometryConfig` | Outputs `team7_case=` |

### Scripts

- `run_team7_uniform_legacy_bz_vtp.sh`: probe adds `--team7-case=uniform-legacy-dp6` and echoes case id.

## B. Alias Header Consolidation (12 files deleted)

After updating includes to point to canonical paths, deleted top-level aliases:

- `aphi_gmres_contact_solver_ck.*` → `aphi_gmres_solver_ck.*`
- `aphi_matrix_free_contact_solve_ck.*` → `aphi_matrix_free_solve_ck.*`
- `aphi_contact_assemble_lhs_debug_ck.*` → `diagnostics/aphi_assemble_lhs_debug_ck.*`
- `aphi_contact_jacobi_debug_ck.*` → `aphi_block_jacobi_preconditioner_ck.*`
- `aphi_bicgstab_solver_ck.*` / `aphi_pcg_solver_ck.*` → `alternate_krylov/` only (top-level aliases removed)

## C. Legacy Route Directory Deletion

After `rg` confirmed no `test_3d_aphi_ck_*` dependencies, deleted:

```text
extra_src/shared/aphi_sphinxsys/
extra_src/shared/legacy_aphi_archive/
extra_src/shared/aphi_case_support/
3d_examples/legacy_aphi_archive/
```

(Legacy tests such as `test_3d_em_*` were already deleted in the working tree.)

## D. Retained

- All `electromagnetic_dynamics/` production + `test_3d_aphi_ck_*`
- TEAM7 native P0–P4 tests and `particle_generation_em`
- `alternate_krylov/` (off main path, not deleted for now)
- `diagnostics/` research-oriented helpers

## E. Minimum Gate Tests (to run before committing this record)

```text
test_3d_aphi_ck_team7_native_geometry_audit
test_3d_aphi_ck_team7_native_source_rhs_audit
test_3d_aphi_ck_team7_native_vacuum_source_b_sanity
test_3d_aphi_ck_team7_native_geometry_bz_reference_probe
test_3d_aphi_ck_fused_apply_vs_debug_lhs
test_3d_aphi_ck_contact_apply_vs_monolithic
test_3d_aphi_ck_contact_fused_apply_equivalence
test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured
test_3d_aphi_ck_joule_uniform_field_analytic_verification
```

Results in **Gate Run Record** below (filled in by agent).

## F. Explicitly Out of Scope for This Stage

- big-air MR, GMRES tuning, Joule thermal coupling, profile <10%, A-penalty research

## Gate Run Record

> 2026-06-05, sequential execution under `build/`; P3 uses `--team7-case=uniform-legacy-dp6 --team7-reload-case=uniform_legacy_dp6`.

| Test | Result | Notes |
|------|------|------|
| `test_3d_aphi_ck_team7_native_geometry_audit` | PASS | `passed=1` |
| `test_3d_aphi_ck_team7_native_source_rhs_audit` | PASS | `passed=1` |
| `test_3d_aphi_ck_team7_native_vacuum_source_b_sanity` | PASS | `final_true_rel≈0.006` |
| `test_3d_aphi_ck_team7_native_geometry_bz_reference_probe` | PASS | `team7_case=uniform-legacy-dp6`, `Bz≈7.10 mT`, `profile L2 real≈13%` |
| `test_3d_aphi_ck_fused_apply_vs_debug_lhs` | PASS | |
| `test_3d_aphi_ck_contact_apply_vs_monolithic` | PASS | |
| `test_3d_aphi_ck_contact_fused_apply_equivalence` | PASS | |
| `test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured` | PASS | |
| `test_3d_aphi_ck_joule_uniform_field_analytic_verification` | PASS | |

**Build fix:** `team7NativeAuditReloadMetadata` must be placed after `loadTeam7NativeGeometryMetadata` (otherwise undeclared identifier).
