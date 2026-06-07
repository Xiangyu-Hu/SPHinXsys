# Stage 10 A-phi Test Case Index

> **~121** `test_3d_aphi_ck_*` directories total (2026-06); legacy `test_3d_em_*` moved to `legacy_aphi_archive/` (not built).
>
> Categories: **demo** / **production gate** / **diagnostic** / **historical scaffold**.
> Code deletion/merging pending GPT discussion; this index is for navigation.

---

## Legend

| Label | Meaning |
|---|---|
| 🎯 **Demo** | Main line for afternoon discussion with advisor |
| ✅ **gate** | Production regression / operator equivalence / convergence acceptance |
| 🔬 **diagnostic** | Informational output; no mandatory passed=1 or sweep-only |
| 📦 **scaffold** | Engineering prototype; metrics not fully closed |
| 🧪 **Stage10.7** | Inner A-divergence penalty research (not demo main line) |

---

## A0. TEAM7 Native Main Line (🎯 current benchmark, Stage 10.17+)

| Case | Purpose | Status |
|---|---|---|
| `particle_generation_em` | STL relax + reload + MR | ✅ |
| `test_3d_aphi_ck_team7_native_geometry_audit` | P0 geometry audit | ✅ passed |
| `test_3d_aphi_ck_team7_native_source_rhs_audit` | P1 source RHS | ✅ passed |
| `test_3d_aphi_ck_team7_native_vacuum_source_b_sanity` | P2 vacuum Bz | ✅ passed |
| `test_3d_aphi_ck_team7_native_geometry_bz_reference_probe` | **P3 Bz vs muFEM** | ✅ passed |
| `test_3d_aphi_ck_team7_native_geometry_source_driven_smoke` | Solve smoke | ✅ |
| `test_3d_aphi_ck_team7_native_joule_post_smoke` | Joule post (deferred) | scaffold |

Records: `CURSOR_APHI_STAGE10_17_TEAM7_NATIVE_RECORD.md`; MR: `CURSOR_APHI_TEAM7_SMALL_AIR_MULTIRES_PLAN.md`.

One-click smoke: `run_team7_native_si_smoke.sh`.

---

## A. Contact / Demo Main Line (🎯 operator and Contact validation)

| Case | Purpose | Key Metrics |
|---|---|---|
| `test_3d_aphi_ck_contact_apply_vs_monolithic` | Contact vs monolithic operator equivalence | `lhs_max_diff`, `passed` |
| `test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured` | Two-body interface MMS + coupled GMRES | `contact_converged`, `interface_band_rel` |
| `test_3d_aphi_ck_contact_phi_gauge_penalty_sweep_diagnostic` | φ gauge λ_φ controls raw drift | sweep table |
| `test_3d_aphi_ck_gmres_contact_team7_three_body_scaffold` | Three-body box TEAM7 + E/J/Joule VTP | `passed`, plate Joule |

**Backup (operator equivalence stdout-only)**: `test_3d_aphi_ck_contact_fused_apply_equivalence`

---

## B. Production Gate (✅ core regression)

### B1. Operator / fused vs debug

| Case | Description |
|---|---|
| `test_3d_aphi_ck_fused_apply_vs_debug_lhs` | Single-body fused vs debug LHS |
| `test_3d_aphi_ck_fused_apply_phi_gauge_penalty_vs_debug_lhs` | φ penalty fused vs debug |
| `test_3d_aphi_ck_inner_a_divergence_penalty_fused_vs_debug` | 🧪 A-penalty fused vs debug |
| `test_3d_aphi_ck_contact_fused_apply_equivalence` | Contact fused apply |
| `test_3d_aphi_ck_contact_pairwise_laplace_equivalence` | Contact Laplace equivalence |
| `test_3d_aphi_ck_contact_block_jacobi_diagonal_equivalence` | Contact Jacobi diagonal |
| `test_3d_aphi_ck_contact_team7_three_body_apply_vs_monolithic` | Three-body apply vs mono |

### B2. Stage 0–3 LHS sub-terms

| Case | Description |
|---|---|
| `test_3d_aphi_ck_variable_registration` | Variable registration |
| `test_3d_aphi_ck_lhs_full_zero` | All zero |
| `test_3d_aphi_ck_lhs_grad_only` | grad coupling |
| `test_3d_aphi_ck_lhs_div_only` | div coupling |
| `test_3d_aphi_ck_lhs_reaction_only` | reaction |
| `test_3d_aphi_ck_lhs_discrete_rhs_residual` | discrete RHS |
| `test_3d_aphi_ck_lhs_polynomial_manufactured` | polynomial MMS |
| `test_3d_aphi_ck_pairwise_laplace_manufactured` | pairwise Laplace MMS |
| `test_3d_aphi_ck_gradient_divergence_manufactured` | grad/div MMS |
| `test_3d_aphi_ck_apply_search_and_s_blocks` | search/s blocks |
| `test_3d_aphi_ck_block_vector_ops` | block vector ops |

### B3. Krylov convergence (single-body)

| Case | Description |
|---|---|
| `test_3d_aphi_ck_gmres_scalar_phi_laplace_penalty_convergence` | scalar φ + penalty GMRES |
| `test_3d_aphi_ck_gmres_vector_helmholtz_convergence` | vector Helmholtz |
| `test_3d_aphi_ck_gmres_full_aphi_penalty_convergence` | full A-phi + penalty |
| `test_3d_aphi_ck_gmres_interface_flux_matched_manufactured` | single-body interface MMS |
| `test_3d_aphi_ck_gmres_two_material_interface_manufactured` | two-material interface |
| `test_3d_aphi_ck_pcg_scalar_phi_laplace_penalty_diagnostic` | PCG scalar sanity |
| `test_3d_aphi_ck_block_jacobi_diagonal_diagnostic` | Jacobi diagonal positive |

### B4. Joule / observable

| Case | Description |
|---|---|
| `test_3d_aphi_ck_joule_uniform_field_analytic_verification` | uniform-field Joule analytic |
| `test_3d_aphi_ck_joule_heating_plumbing_diagnostic` | Joule plumbing |

---

## C. Diagnostic / Sweep (🔬 informational, not demo main line)

### C1. Contact diagnostics

| Case | Description |
|---|---|
| `test_3d_aphi_ck_contact_two_body_polish_sweep_diagnostic` | polish count sweep |
| `test_3d_aphi_ck_contact_left_field_error_diagnostic` | left-body field error |
| `test_3d_aphi_ck_contact_interface_rhs_consistency_diagnostic` | interface RHS |
| `test_3d_aphi_ck_contact_phi_gauge_lambda100_ab_diagnostic` | λ=100 A/B comparison |
| `test_3d_aphi_ck_contact_team7_phi_gauge_sweep_diagnostic` | TEAM7 φ sweep |
| `test_3d_aphi_ck_gmres_workspace_contamination` | workspace contamination |

### C2. GMRES / PC sensitivity

| Case | Description |
|---|---|
| `test_3d_aphi_ck_gmres_restart_sensitivity_diagnostic` | restart sensitivity |
| `test_3d_aphi_ck_gmres_preconditioner_comparison_diagnostic` | PC comparison |
| `test_3d_aphi_ck_gmres_manufactured_robustness_sweep_diagnostic` | robustness sweep |
| `test_3d_aphi_ck_block_jacobi_8x8_pc_diagnostic` | 8×8 block PC |
| `test_3d_aphi_ck_laplace_adjointness_diagnostic` | Laplace self-adjoint |
| `test_3d_aphi_ck_laplace_constant_nullspace_diagnostic` | constant nullspace |

### C3. Scalar φ diagnostics

| Case | Description |
|---|---|
| `test_3d_aphi_ck_scalar_phi_tiny_host_matrix_diagnostic` | small matrix comparison |
| `test_3d_aphi_ck_scalar_phi_constant_mode_rhs_diagnostic` | constant-mode RHS |
| `test_3d_aphi_ck_scalar_phi_random_adjointness_diagnostic` | random self-adjoint |
| `test_3d_aphi_ck_scalar_phi_random_curvature_diagnostic` | random curvature |

### C4. Joule sensitivity

| Case | Description |
|---|---|
| `test_3d_aphi_ck_joule_local_distribution_sensitivity_diagnostic` | local distribution |
| `test_3d_aphi_ck_joule_em_tolerance_sensitivity_diagnostic` | EM tol sensitivity |

### C5. BiCGStab legacy path

| Case | Description |
|---|---|
| `test_3d_aphi_ck_bicgstab_discrete_rhs_solve` | discrete RHS |
| `test_3d_aphi_ck_bicgstab_exact_initial_guess` | exact initial guess |
| `test_3d_aphi_ck_bicgstab_scalar_phi_laplace_penalty_diagnostic` | scalar φ |
| `test_3d_aphi_ck_bicgstab_vector_helmholtz_convergence` | vector |
| `test_3d_aphi_ck_bicgstab_full_aphi_penalty_convergence` | full A-phi |
| `test_3d_aphi_ck_bicgstab_vector_laplace_singular_smoke` | smoke |

---

## D. TEAM7 / Engineering Scaffold (📦)

| Case | Description | Status |
|---|---|---|
| `test_3d_aphi_ck_gmres_team7_like_coil_plate_simplified` | simplified TEAM7 | engineering |
| `test_3d_aphi_ck_gmres_team7_like_coil_plate_tight_tol_diagnostic` | tight tol | diagnostic |
| `test_3d_aphi_ck_gmres_team7_like_solver_study_diagnostic` | solver study | diagnostic |
| `test_3d_aphi_ck_gmres_team7_dp_convergence_diagnostic` | dp convergence | diagnostic |
| `test_3d_aphi_ck_gmres_team7_physical_dimensions_smoke` | physical dimensions smoke | smoke |
| `test_3d_aphi_ck_gmres_team7_physical_dp_observable_diagnostic` | dp observable | diagnostic |
| `test_3d_aphi_ck_gmres_team7_quantitative_reference_comparison` | quantitative comparison | in research |
| `test_3d_aphi_ck_gmres_coupled_pc_team7_convergence` | coupled PC | in research |
| `test_3d_aphi_ck_gmres_impressed_current_homogeneous_box` | impressed J box | scaffold |
| `test_3d_aphi_ck_gmres_cold_crucible_scaffold` | cold crucible GMRES scaffold | ✅ |
| `test_3d_aphi_ck_cold_crucible_demo` | cold crucible impressed flow demo + VTP (**non-quantitative**) | ✅ |

---

## E. Stage 10.7 Additions — Inner A-divergence Penalty (🧪 all diagnostic/research)

**10 total**, for A-penalty prototype validation; **not the afternoon demo main line**.

| # | Case | Type | Current Conclusion |
|---|---|---|---|
| 1 | `test_3d_aphi_ck_inner_a_divergence_penalty_fused_vs_debug` | ✅ gate | fused vs debug match, `passed=1` |
| 2 | `test_3d_aphi_ck_inner_a_divergence_penalty_sign_energy_diagnostic` | 🔬 | sign/energy diagnostic pass |
| 3 | `test_3d_aphi_ck_impressed_current_divergence_diagnostic` | 🔬 | divJ ≈ 4.88 (non-solenoidal source) |
| 4 | `test_3d_aphi_ck_solenoidal_current_divergence_diagnostic` | 🔬 | divJ ≈ 2.67 (solenoidal J=curl C) |
| 5 | `test_3d_aphi_ck_curl_a_manufactured_diagnostic` | 🔬 | curl(A) MMS pass |
| 6 | `test_3d_aphi_ck_graddiv_block_diagonal_diagnostic` | 🔬 | 3×3 block PC, FD diff ≈ 0.17 |
| 7 | `test_3d_aphi_ck_inner_a_divergence_penalty_sweep_diagnostic` | 🔬 | impressed λ_A sweep |
| 8 | `test_3d_aphi_ck_inner_solenoidal_a_divergence_penalty_sweep_diagnostic` | 🔬 | solenoidal λ_A sweep |
| 9 | `test_3d_aphi_ck_inner_a_divergence_penalty_gmres_sensitivity_diagnostic` | 🔬 | **0/30 GMRES converged** |
| 10 | `test_3d_aphi_ck_inner_a_divergence_penalty_tol_closure_diagnostic` | 🔬 | **0/13 strict tol converged** |

### Stage 10.7 Key Numerical Summary

- **Operator**: sign corrected to `LhsA -= λ_A grad(div A)`; fused ≡ debug.
- **3×3 JacobiGradDivABlock PC**: at λ=100 true_rel 0.6→**0.07**; at λ=1000 div_A≈**0.0027**.
- **GMRES bottleneck**: stalls ~4% (impressed) / ~8% (solenoidal); increasing outer/restart/tol **ineffective**.
- **strict tol=1e-4**: penalty path **not converged**; baseline without penalty can converge.
- **Recommendation**: keep all 10 cases as research record; for afternoon discussion only mention verbally that "divA control is still in progress".

Detailed record: `CURSOR_APHI_STAGE10_7_A_DIVERGENCE_PENALTY_RECORD.md`

---

## F. Quick Run Commands

```bash
cd /home/yongchuan/sphinxsys/build
ninja <target_name>
./tests/extra_source_and_tests/3d_examples/<case_dir>/bin/<case_dir>
```

Demo quartet:

```bash
ninja test_3d_aphi_ck_contact_apply_vs_monolithic \
      test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured \
      test_3d_aphi_ck_contact_phi_gauge_penalty_sweep_diagnostic \
      test_3d_aphi_ck_gmres_contact_team7_three_body_scaffold
```

---

## G. Follow-Up Case Cleanup (Level B)

- [ ] Add CMake comment `LEGACY` to obsolete BiCGStab cases
- [ ] Consider moving Stage 10.7 directories into `3d_examples/diagnostics/stage10_7/` (directory only, target names unchanged)
- [ ] Merge duplicate TEAM7 diagnostic cases
- [ ] Create `stage10_teacher_demo/` archive for stdout/summary/VTP
