# Stage 10 A-phi Source File Map (For Discussion)

> Path root: `tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/`
>
> As of 2026-06 organized into **production / diagnostics / test_helpers / benchmark** subdirectories; top-level **one-line alias** headers see `electromagnetic_dynamics/README.md` (consolidation pending GPT discussion).

---

## 1. Production — Operators and Solvers (Afternoon Demo Main Line)

These files represent the reusable SPHinXsys kernel-form A-phi solver stack.

| File | Responsibility |
|---|---|
| `aphi_coupling_modes_ck.h` | Coupling modes and switches: `phi_gauge_penalty`, `a_divergence_penalty`, etc. |
| `aphi_matrix_free_operator_ck.*` | **Main entry**: fused matrix-free `K(input)` apply |
| `aphi_laplace_ck.*` | A / phi pairwise Laplace |
| `aphi_grad_phi_coupling_ck.*` | `σ grad(φ)` coupling in A equation |
| `aphi_div_sigma_a_coupling_ck.*` | `div(σA)` coupling in φ equation (current continuity) |
| `aphi_phi_gauge_penalty_ck.*` | φ gauge regularization: `LhsPhi += λ_φ φ` |
| `aphi_reaction_ck.*` | Local reaction: `LhsARe += -ωσAIm`, etc. |
| `aphi_joule_heating_ck.*` | Post-processing: E / J / Joule |
| `aphi_block_jacobi_preconditioner_ck.*` | Block Jacobi / approximate preconditioner; **includes** `AphiComputeBlockJacobiContactDynamicsBundle` |
| `aphi_gmres_solver_ck.*` | inner / single-body GMRES; **includes** `AphiGMRESContactSolverCK` |
| `aphi_gmres_workspace_ck.*` | Krylov basis variable registration |
| `aphi_multibody_contact_gmres_ck.*` | coupled multi-body Contact GMRES main line |
| `aphi_matrix_free_solve_ck.*` | matrix-free solve wrapper; **includes** `AphiMatrixFreeContactSolveCK` |
| `diagnostics/aphi_assemble_lhs_debug_ck.*` | debug assemble; **includes** `AphiAssembleLhsDebugContactDynamicsBundle` |
| `aphi_block_vector_ops_ck.*` | block vector dot/norm, etc. |
| `aphi_block_zero_ck.*` | Zero / Copy / Residual |
| `aphi_variables_ck.*` | Variable registration |
| `aphi_field_names_ck.h` | Naming conventions |
| `aphi_gmres_solver_ck.*` | **production GMRES** (Inner + Contact share `detail::aphiGMRESSolve`) |
| `alternate_krylov/aphi_pcg_solver_ck.*` | PCG (Stage 8 scalar φ sanity, off main line) |
| `alternate_krylov/aphi_bicgstab_solver_ck.*` | BiCGStab (Stage 8 comparison research, off main line) |
| `aphi_krylov_diagnostics_ck.h` | Krylov shared finite/breakdown helpers |
| `all_electromagnetic_dynamics_ck.h` | production unified include |

### Stage 10.7 Prototype Operators (Still in research, **not afternoon main line**)

| File | Status |
|---|---|
| `aphi_a_divergence_penalty_ck.*` | Operator implemented: `LhsA -= λ_A grad(div A)` |
| `diagnostics/aphi_a_divergence_penalty_pipeline.h` | debug assemble orchestration |
| `diagnostics/aphi_inner_a_divergence_penalty_sweep_helpers.h` | λ_A sweep + GMRES parameter scan |

**Current conclusion**: fused vs debug match; divA improves after 3×3 block PC, but GMRES still stalls ~4–8% rel under strict tol. **Not closed; continue research.**

---

## 2. Diagnostics — Debug and Metrics (Not Shown as Production)

Directory: `diagnostics/`

| File | Purpose |
|---|---|
| `diagnostics/aphi_assemble_lhs_debug_ck.*` | Sub-term debug assemble (Inner + Contact same file) |
| `aphi_contact_assemble_lhs_debug_ck.*` (top-level alias) | Forward → `diagnostics/aphi_assemble_lhs_debug_ck.*` |
| `aphi_gmres_contact_solver_ck.*` (top-level alias) | Forward → `aphi_gmres_solver_ck.*` |
| `aphi_matrix_free_contact_solve_ck.*` (top-level alias) | Forward → `aphi_matrix_free_solve_ck.*` |
| `diagnostics/aphi_gradient_divergence_debug_ck.*` | gradA / divA debug pipeline |
| `aphi_curl_a_diagnostic_ck.*` | curl(A) manufactured-solution diagnostic |
| `aphi_a_divergence_penalty_pipeline.h` | A-divergence penalty debug orchestration |
| `aphi_div_a_diagnostic_helpers.h` | divA L2/L∞ and other metrics |
| `aphi_a_gauge_diagnostic_helpers.h` | divA energy, sign energy, etc. |
| `aphi_graddiv_block_diagnostic_helpers.h` | grad-div block diagonal FD comparison |
| `aphi_inner_a_divergence_penalty_sweep_helpers.h` | inner λ_A sweep |
| `aphi_contact_phi_gauge_sweep_helpers.h` | Contact φ gauge λ_φ sweep |
| `aphi_contact_left_field_error_helpers.h` | left-body field error decomposition |
| `aphi_contact_bodywise_residual_helpers.h` | per-body residual |
| `aphi_contact_workspace_contamination_helpers.h` | workspace contamination detection |
| `aphi_contact_readback_sync_helpers.h` | device readback sync |
| `aphi_contact_interface_diagnostic_helpers.h` | interface band metrics |
| `aphi_scalar_phi_diagnostic_helpers.h` | scalar φ diagnostics |
| `aphi_gmres_robustness_sweep_helpers.h` | GMRES robustness sweep |

### TEAM7 Native (Stage 10.17+, current benchmark main line)

| File | Purpose |
|---|---|
| `diagnostics/aphi_team7_native_geometry_config.h` | air preset, MR levels, reload_cases CLI |
| `diagnostics/aphi_team7_native_team7_shapes.h` | Adaptive STL shapes |
| `diagnostics/aphi_team7_native_reload_geometry_helpers.h` | three-body Contact solve, VTP, Joule pipeline |
| `diagnostics/aphi_team7_native_coil_source_helpers.h` | centerline tangent source J₀=NI/A |
| `diagnostics/aphi_team7_native_bz_reference_probe_helpers.h` | Bz probe vs muFEM CSV |
| `diagnostics/aphi_team7_native_mr_probe_helpers.h` | MR adaptive kernel probe sampling |

Record: `CURSOR_APHI_STAGE10_17_TEAM7_NATIVE_RECORD.md`.

---

## 3. Test Helpers — Test Use Only

Directory: `test_helpers/`

| File | Purpose |
|---|---|
| `aphi_contact_test_helpers.h` | two-body geometry / Contact setup |
| `aphi_contact_gmres_test_helpers.h` | two-body MMS + coupled GMRES |
| `aphi_team7_contact_test_helpers.h` | three-body TEAM7 scaffold |
| `aphi_em_observable_helpers.h` | E/J/Joule output and comparison |
| `aphi_gmres_test_helpers.h` | generic GMRES test driver |
| `aphi_gmres_benchmark_helpers.h` | benchmark run wrapper |
| `aphi_lhs_test_helpers.h` | LHS comparison helpers |
| `aphi_test_device_sync.h` | host/device sync |

**Note**: These are not SPHinXsys core framework types; they exist only to avoid duplicated code across 70+ tests.

---

## 4. Benchmark — Standard Case Definitions

Directory: `benchmark/`

| File | Purpose |
|---|---|
| `aphi_benchmark_case_ck.*` | generic MMS benchmark (includes solenoidal J assignment) |
| `aphi_team7_canonical_case_ck.h` | TEAM7-like canonical geometry/materials |
| `aphi_cold_crucible_case_ck.*` | cold crucible four-body scaffold (follow-up) |

---

## 5. Recommended File Subset for Afternoon Discussion

**Must cover (4–6 files)**:

1. `aphi_matrix_free_operator_ck.*` — kernel-form operator
2. `aphi_multibody_contact_gmres_ck.*` — Contact GMRES
3. `aphi_phi_gauge_penalty_ck.*` — φ gauge control
4. `aphi_joule_heating_ck.*` — physical observables

**Optional deep dive**:

5. `aphi_block_jacobi_preconditioner_ck.*` — preconditioner
6. `test_helpers/aphi_contact_gmres_test_helpers.h` — how two-body MMS drives tests

**Not recommended for main-line coverage**:

- `diagnostics/aphi_inner_a_divergence_penalty_sweep_helpers.h`
- `aphi_a_divergence_penalty_ck.*` (unless advisor asks about divA)

---

## 6. Follow-Up Level B Cleanup (After Afternoon)

- [ ] Merge `diagnostics/aphi_contact_*_helpers.h` → single `aphi_contact_diagnostic_helpers.h`
- [ ] Remove top-level forwarding headers; unify include paths
- [ ] Mark obsolete solver paths (BiCGStab main path)
- [ ] Common GMRES driver abstraction
