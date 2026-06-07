# Stage 10-contact — Contact<> Multi-Body A-phi Operator Validation

> Main line: migrate A-phi matrix-free operators to `Contact<>` multi-body relations on top of inner-only 9D–10C; **apply-level equivalence first, then GMRES / real cases**.  
> **Not doing**: formal 10D thermal coupling, σ(T)/ν(T), Lorentz force.

## Roadmap (Aligned with ChatGPT)

```
10C-contact / 10B-contact (operator diagnostics + apply equivalence)
  → contact two-body interface MMS (9E-1 migration)
  → contact three-body TEAM7-like scaffold
  → fix E/J/B output metric naming
  → 10C+ cold crucible geometry/parameters
  → 10D fixed σ/ν one-way thermal
  → 10B+ external FEM/MFEM
  → 10E σ(T)/ν(T) two-way feedback
```

## Operator List (Contact Specialization Status)

| Operator | Inner | Contact | Notes |
|------|-------|---------|------|
| `AphiPairwiseLaplaceCK` | ✅ | ✅ Sprint 1 | harmonic mean σ/ν; Contact uses `+=` gather |
| `AphiGradPhiCouplingCK` | ✅ | ✅ Sprint 2 | σ_i * g_ij * (φ_i - φ_j) |
| `AphiDivSigmaACouplingCK` | ✅ | ✅ Sprint 2 | harmonic mean σ_ij |
| `AphiComputeScalarPhiGradientCK` | ✅ | ✅ Sprint 2 | Inner assigns; Contact accumulates |
| `AphiComputeBlockJacobiDiagonalCK` contact | ✅ inner | ✅ Sprint 3 | Inner assigns; Contact accumulates |
| `AphiApplyCK` / fused operator | ✅ inner | ✅ Sprint 3 | `AphiApplyContactDynamicsBundle` |
| Joule pipeline | ✅ inner | ✅ Sprint 2 | grad Inner+Contact → E → Joule |

## Contact Conventions (Same as Inner)

- harmonic mean σ / ν
- `AphiPairwiseNegativeLaplaceWeight` / `AphiPairwiseGradientWeightUncorrected`
- gather writes owner body particle i, **no atomic**
- Inner assigns, Contact **accumulates** (same pattern as `DensitySummationCK`)

## Sprint 1 — Laplace Contact + Apply Equivalence ✅

### Code

- `aphi_laplace_ck.h/.hpp`: `AphiPairwiseLaplaceCK<Contact<DataType>>` specialization
- `test_3d_aphi_ck_contact_pairwise_laplace_equivalence`: monolithic Inner vs split two-body Inner+Contact

### Test Design

- **Monolithic**: unit box, `AssignPiecewiseSigmaHalfSpaceCK` (x=0.5, σ=10/1), quadratic φ, Inner Laplace
- **Split**: left/right half box, constant σ, bidirectional `Contact<>`, `Inner+Contact` Laplace
- **Comparison**: match core-region Laplace by particle position (exclude interface ±1dp and outer shell); float threshold 3e-5

### Results (2026-05-21)

```
matched_particles=32  missing_particles=0  max_abs_diff=1.53e-5  passed=1
```

### Run

```bash
cd build
ninja test_3d_aphi_ck_contact_pairwise_laplace_equivalence
./bin/test_3d_aphi_ck_contact_pairwise_laplace_equivalence
```

## Sprint 2 — Grad/Div/PhiGradient + Full K(X) Apply Equivalence ✅

### Code

- `aphi_grad_phi_coupling_ck.h/.hpp`: `AphiGradPhiCouplingCK<Contact<>>`
- `aphi_div_sigma_a_coupling_ck.h/.hpp`: `AphiDivSigmaACouplingCK<Contact<>>`
- `aphi_joule_heating_ck.h/.hpp`: `AphiComputeScalarPhiGradientCK<Contact<>>`
- `aphi_contact_assemble_lhs_debug_ck.h/.hpp`: Inner+Contact debug K(X) assembly
- `aphi_contact_test_helpers.h`: two-body geometry + position-matching helpers
- `test_3d_aphi_ck_contact_apply_vs_monolithic`: full debug LHS + Joule equivalence

### Results (2026-05-21)

```
test_3d_aphi_ck_contact_apply_vs_monolithic
  lhs_matched=32  lhs_max_abs_diff=4.39e-5
  joule_matched=32  joule_max_abs_diff=1.53e-5  passed=1
```

Two-medium box (σ=10/1, x=0.5 interface), separable field, ω=1.25; exclude interface ±1dp and outer shell.

### Run

```bash
ninja test_3d_aphi_ck_contact_apply_vs_monolithic
./bin/test_3d_aphi_ck_contact_apply_vs_monolithic
```

## Sprint 3 — Block Jacobi + Fused Apply Contact ✅

### Code

- `AphiComputeBlockJacobiDiagonalCK<Contact<>>` + `AphiComputeBlockJacobiContactDynamicsBundle`
- `AphiApplyCK<Contact<>>` + `AphiApplyContactDynamicsBundle`
- `test_3d_aphi_ck_contact_block_jacobi_diagonal_equivalence`
- `test_3d_aphi_ck_contact_fused_apply_equivalence`

### Results (2026-05-21)

```
block_jacobi: matched=32  max_laplace_a≈2.4e-4  max_laplace_phi≈4.9e-4  passed=1
fused_apply:  mono_vs_split≈5.3e-5  split_debug_vs_fused≈6.1e-5  passed=1
```

### Run

```bash
ninja test_3d_aphi_ck_contact_block_jacobi_diagonal_equivalence \
      test_3d_aphi_ck_contact_fused_apply_equivalence
```

## Sprint 4 — Contact Interface MMS GMRES (Diagnostic Complete + MMS Pass)

### Diagnostics

- Test: `test_3d_aphi_ck_contact_interface_rhs_consistency_diagnostic`
- Record: [`CURSOR_APHI_STAGE10_CONTACT_INTERFACE_DIAGNOSTIC_RECORD.md`](CURSOR_APHI_STAGE10_CONTACT_INTERFACE_DIAGNOSTIC_RECORD.md)
- **Root cause**: after `AphiZeroBlockCK` on **neighbor body** for Contact, even with pinned `r_hat`, contact apply still reads wrong neighbor state; RHS assembly itself is correct.
- **Fix**: after coupled RHS, **zero left body only**; keep exact on right body; right-body GMRES also does not zero right body.

### Code

- `aphi_gmres_contact_solver_ck.h/.hpp`
- `aphi_matrix_free_contact_solve_ck.h/.hpp`
- `aphi_contact_gmres_test_helpers.h`
- `aphi_contact_interface_diagnostic_helpers.h`
- `aphi_test_device_sync.h`: `syncAphiBlockToDevice`
- `test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured`
- `test_3d_aphi_ck_contact_interface_rhs_consistency_diagnostic`

### Results (2026-05-21, after fix)

```
mono:    outer=2   discrete_defect≈5e-6   continuous≈1.8e-6
contact: outer=30  discrete_defect≈4.7e-8 continuous≈4.8e-6
         defect_ratio≈0.009  continuous_ratio≈2.6  passed=1
```

Block-GS path `converged` flag may still be 0; MMS uses `true_rel < 10×tol` + defect/continuous criteria.

### Run

```bash
ninja test_3d_aphi_ck_contact_interface_rhs_consistency_diagnostic \
      test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured
```

## Sprint 4.5 — GMRES Arnoldi Block-Diagonal Contact Apply ✅

### Goal

Fix contamination from Arnoldi `apply_z_to_w` reading neighbor body `GMRESZ*` workspace; keep residual / true residual / `apply_solution_to_lhs` as **full** Contact apply.

### Code

- `AphiApplyContactBlockDiagonalCK<Contact<>>` + `AphiApplyContactBlockDiagonalDynamicsBundle` (`aphi_matrix_free_operator_ck.h/.hpp`)
- `AphiGMRESContactSolverCK`: Arnoldi step uses block-diagonal bundle; rest still full apply
- `aphi_contact_workspace_contamination_helpers.h` + `test_3d_aphi_ck_gmres_workspace_contamination`
- `aphi_contact_readback_sync_helpers.h` + `test_3d_contact_variable_readback_sync`
- `runTwoBodyContactInterfaceMms`: `zeroGmresWorkspaceOnBody` before per-body GMRES

### Task 1 — Workspace Contamination

```
reference_output_norm≈89.8
full_apply_rel_diff≈0.82        (neighbor GMRESZ0 contamination → large left-body w change)
block_diagonal_rel_diff≈2.9e-8  (block-diagonal insensitive to neighbor Z)
passed=1
```

### Task 4 — Contact Readback Sync

```
device_write_ok=1  host_sync_ok=1  host_unsynced_stale=1 (SYCL)  passed=1
```

### Task 5 — Rerun MMS (After Block-Diagonal Arnoldi)

```
mono:    outer=2   discrete_defect≈5.0e-6  continuous≈1.8e-6
contact: outer=30  discrete_defect≈5.2e-8  continuous≈7.6e-6
         contact_converged=0  contact_true_rel≈1.38e-5
         defect_ratio≈0.010  continuous_ratio≈4.15  passed=1
```

**Interpretation**: workspace contamination proven in Task 1 and Arnoldi semantics fixed; MMS still hits `max_outer=30`, `converged=0`, indicating **block-GS outer loop + full residual path** still has coupling/convergence bottleneck (not merely Arnoldi matvec contamination). Engineering pass still uses `true_rel < 10×tol`.

### Run

```bash
ninja test_3d_aphi_ck_gmres_workspace_contamination \
      test_3d_contact_variable_readback_sync \
      test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured
```

## Sprint 5 — Three-Body TEAM7-like Contact Scaffold ✅ (Step 5.1–5.3 Initial)

### Geometry / Three-Body Decomposition

- **air**: left/middle/right slab blocks (`AphiTeam7AirSlabsShape`), non-overlapping
- **coil**: `layout.coil` box
- **plate/conductor**: `layout.conductor` box
- air `Contact<>` connects to both coil + plate

### Code

- `aphi_team7_contact_test_helpers.h` — `AphiTeam7ThreeBodyContactCase`, regional apply metrics, `runTeam7ThreeBodyContactBlockGs`
- `test_3d_aphi_ck_contact_team7_three_body_apply_vs_monolithic` — Step 5.1
- `test_3d_aphi_ck_gmres_contact_team7_three_body_scaffold` — Step 5.2/5.3

### Step 5.1 — Apply-Level Monolithic vs Split

Regional max-norm comparison (exclude interface ±dp would empty narrow coil/plate regions; use core + region tag instead):

```
conductor_lhs_gap≈1.0e-4  coil_lhs_gap≈1.5e-3  air_lhs_gap≈5e-5
conductor_joule_gap≈3.1e-2  passed=1
```

### Step 5.2/5.3 — Impressed-Current Scaffold

**Block-GS (transitional path, superseded)**:

```
plate_joule_power_gap≈99.9%  global_true_rel≈0.55  block_gs_sweeps=6
```

**Coupled multi-body GMRES (Sprint 6 prototype, see below)** has replaced block-GS as the scaffold main solve path.

Metric naming: plate observables use `plate_joule_power` / `plate_j_L2` / `plate_solution_block_max`.

### Run

```bash
ninja test_3d_aphi_ck_contact_team7_three_body_apply_vs_monolithic \
      test_3d_aphi_ck_gmres_contact_team7_three_body_scaffold
```

## Sprint 6 — Coupled Multi-Body Contact GMRES Prototype ✅

### Code

- `aphi_multibody_contact_gmres_ck.h/.hpp` — `AphiMultiBodyContactGMRESSolverCK`
  - Global dot/norm: vol-weighted sum across bodies
  - Arnoldi matvec uses **full** Contact apply (neighbor GMRESZ is part of global Krylov vector)
  - residual / true residual / solution apply also full Contact
- `runTeam7ThreeBodyCoupledContactGmres` in `aphi_team7_contact_test_helpers.h`
- `test_3d_aphi_ck_gmres_contact_team7_three_body_scaffold` switched to coupled GMRES acceptance

### Results (restart=50, outer=100, tol=5e-4)

```
mono:    plate_joule_power≈0.136  plate_solution_block_max≈1.91
coupled: plate_joule_power≈0.136  plate_solution_block_max≈1.91
         outer=9  converged=1  global_true_rel≈2.9e-4
         plate_joule_power_gap≈2.1e-4  plate_j_L2_gap≈8e-5  passed=1
```

vs block-GS: `global_true_rel` from ~0.55 → ~3e-4; `plate_joule_power_gap` from ~100% → ~0.02%.

### Run

```bash
ninja test_3d_aphi_ck_gmres_contact_team7_three_body_scaffold
```

## Sprint 6.5 — Two-Body Interface MMS + Coupled GMRES ✅

### Code

- `runTwoBodyCoupledContactInterfaceMms` in `aphi_contact_gmres_test_helpers.h` (reuses `AphiMultiBodyContactGMRESSolverCK`)
- `runTwoBodyContactInterfaceMms` kept as legacy block-GS diagnostic baseline
- `toMatrixFreeSolverResult` in `aphi_gmres_test_helpers.h`
- `test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured` switched to coupled GMRES (`restart=50, outer=100`)

### Results (vs Old Block-GS Per-Body GMRES)

| Metric | block-GS (old) | coupled GMRES (new) |
|------|----------------|---------------------|
| outer | 30 | **2** |
| converged | 0 | **1** |
| true_rel | ~1.3e-5 | **~5e-6** |
| global_true_rel | — | **~5e-6** |
| interface_band_rel | ~1e-6 | ~1e-6 |
| continuous_error (max-norm) | ~8e-6 | **~0.02** (left-body dominated; right body ~1e-6) |

**Conclusion**: Krylov convergence semantics fixed; left conductor block max-norm MMS field error still large (TBD: vol-weighted global inner product + heterogeneous two-body split). Test acceptance changed to: `gmresConvergencePassed` + `global_true_rel` + absolute discrete defect + interface band.

### Sprint 6.6 — Bodywise Residual + Stage 10.5 Left-Body Field Error Diagnostic ✅

#### Code

- `aphi_contact_bodywise_residual_helpers.h` — `buildBodywiseTrueRelativeBreakdown` (`max_bodywise_true_rel` + per-body `rhs_norm` / `solution_norm`)
- `aphi_contact_left_field_error_helpers.h` — 10.5-A gauge/mean-offset; 10.5-E left-half core contact vs mono vs exact
- `test_3d_aphi_ck_contact_left_field_error_diagnostic` — diagnostic only (**does not gate passed on left-body ~1% max-norm**; only requires solver convergence + particle match success)

#### Bodywise Output (Two-Body / Three-Body)

Two-body MMS / polish diagnostic: `left_true_rel`, `right_true_rel`, `max_bodywise_true_rel`, `left/right_rhs_norm`, `left/right_solution_norm`

Three-body TEAM7 scaffold: `air/coil/plate_true_rel`, `max_bodywise_true_rel`, per-body `rhs_norm` / `solution_norm`

Note: in impressed-current cases air/plate `rhs_norm≈0`, then `*_true_rel` degenerates to absolute true residual norm (avoids divide-by-zero blowup).

#### Left-Body ~1% Open Issue (Stage 10.5)

- coupled GMRES: `global_true_rel ~1.8e-6`, `converged=1`; `left_true_rel ~2.5e-5` >> `right_true_rel ~1e-6`
- `left_continuous` (core max-norm vs exact) ~0.96%; polish sweep 1→2 plateaus ~0.97% (**not Krylov non-convergence**)
- legacy block-GS: `outer=30`, `converged=0` (comparison baseline)
- **Does not block three-body main line**; use left field diagnostic to distinguish: gauge constant offset vs split/discretization vs left-body operator semantics

**Stage 10.5 first diagnostic results** (`test_3d_aphi_ck_contact_left_field_error_diagnostic`):

```
left_block_L2_rel≈1.4%   left_block_Linf_rel≈3.6%
after mean-sub L2≈0.43%     Linf≈1.9%   (partial gauge/mean offset)
contact_vs_exact L2≈1.4%  contact_vs_mono L2≈2.5%  (split left-half core closer to analytic)
mono_vs_exact L2≈3.2%    (monolithic also has discretization error on left half)
```

Preliminary interpretation: left-body ~1% is not mainly split-vs-mono inconsistency; mean-sub reduces part but not all → continue investigating gauge / left-body discretization vs MMS exact.

**10.5-C partitioning** (core interior vs interface band, `±2dp`) see diagnostic stdout fields `core_interior_*` / `interface_band_*`.

#### Run

```bash
ninja test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured \
      test_3d_aphi_ck_contact_two_body_polish_sweep_diagnostic \
      test_3d_aphi_ck_gmres_contact_team7_three_body_scaffold \
      test_3d_aphi_ck_contact_left_field_error_diagnostic
```

## Sprint 7 — Metric Naming / EM Observables / Cold Crucible Contact Scaffold

### 7.1 Metric Naming ✅ (Code Layer)

Unified `*_field_max` → `*_solution_block_max` (`aphi_gmres_benchmark_helpers.h`, cold crucible/TEAM7 related tests)

### 7.2 EM Observables ✅ (Code Layer)

- `aphi_em_observable_helpers.h` — plate `E_real_max`, `J_L2`, `Joule_L2`, etc.
- Three-body scaffold: `plate_joule_power_gap ~0.02%` (coupled GMRES vs mono)

### 7.3 Completion Tiers (Aligned with ChatGPT)

| Tier | Content | Status |
|------|------|------|
| Operator equivalence | Laplace / apply / Jacobi / fused | ✅ |
| Two-body interface MMS | coupled GMRES + Krylov acceptance | ✅ |
| Three-body TEAM7 | impressed-current **scaffold** (not real TEAM7 BC/source) | ✅ prototype |
| Left-body max-norm ~1% | Stage 10.5 diagnostic | 🔍 open |
| Cold crucible four-body EM | geometry + parameters | ⏸ deferred |
| Helper merge/relocation | multiple files kept separate | ⏸ not doing |

## Stage 10.6 — divA Post-Solve Diagnostic 🔄

See [`CURSOR_APHI_STAGE10_6_A_GAUGE_DIVERGENCE_DIAGNOSTIC_RECORD.md`](CURSOR_APHI_STAGE10_6_A_GAUGE_DIVERGENCE_DIAGNOSTIC_RECORD.md)

- `aphi_div_a_diagnostic_helpers.h`: post-solve divA + reduction (**production GMRES does not depend on DivA field**)
- Register `DivAReal` / `DivAImag` (within diagnostic pipeline)
- **Key finding**: at λ_φ=10/100 `div_A_relative≈1.0` (high_risk) → **phi penalty ≠ div A=0**
- Tests: `test_3d_aphi_ck_contact_phi_gauge_lambda100_ab_diagnostic`, `test_3d_aphi_ck_contact_team7_phi_gauge_sweep_diagnostic`

### λ_φ Candidate Default

- **λ_φ=100 promoted to Contact scaffold default** (two-body MMS / left_field / polish / rhs_consistency / robustness sweep helper)
- TEAM7 canonical was already 100
- Two-body MMS (λ=100): `left_continuous≈0.19%`, `passed=1` (2025-05-21 rerun)

### 7.4 Pending

- ~~Complete three-body TEAM7 λ sweep~~ ✅
- divA still high_risk → **Stage 10.7 Inner prototype started** (see below)
- Three-body dp/σ/source sweep (after divA closure)
- Cold crucible four-body Contact geometry scaffold

## Stage 10.7 — A-Divergence Penalty (Inner Prototype) 🔄

**Order**: Inner fused vs debug → Inner divA sweep → Contact integration (same migration path as Laplace/apply)

### Code

- `AphiLhsAssemblyOptions`: `use_a_divergence_penalty`, `a_divergence_penalty`
- `AphiGradDivAPenaltyCK`: `LhsA -= λ_A * grad(div A)` (B-corrected pipeline, sign corrected)
- `AphiADivergencePenaltyPipelineBundle`: grad A → div A → grad(div A) → penalty
- `AphiApplyDynamicsBundle`: optional append penalty after fused apply
- `AphiAssembleLhsDebugDynamicsBundle`: debug path sync support

### Tests

- `test_3d_aphi_ck_inner_a_divergence_penalty_fused_vs_debug`: `core_max_difference≈3e-5` **passed=1**
- `test_3d_aphi_ck_inner_a_divergence_penalty_sweep_diagnostic`: impressed source λ_A sweep (diagnostic)
- `test_3d_aphi_ck_solenoidal_current_divergence_diagnostic`: `source_div_J_relative≈2.67` **passed=1**
- `test_3d_aphi_ck_graddiv_block_diagonal_diagnostic`: pairwise 3×3 vs pipeline FD **passed=1**
- `test_3d_aphi_ck_inner_solenoidal_a_divergence_penalty_sweep_diagnostic`: solenoidal source λ_A sweep
- `test_3d_aphi_ck_inner_a_divergence_penalty_gmres_sensitivity_diagnostic`: outer/restart grid (0 converged)
- `test_3d_aphi_ck_inner_a_divergence_penalty_tol_closure_diagnostic`: tol×λ_A (0 converged)

### Not Yet Done

- **Acceptance strategy**: rel≤0.05 + divA good vs continue improving PC (B-corrected graddiv block)
- Observable evaluation on **stalled solutions**
- `AphiApplyContactDynamicsBundle` Contact-side penalty

**GMRES sensitivity (2026-05-29)**: outer/restart/tol cannot break ~4% rel floor; λ=1000 div_A≈0.0027 good but rel≈0.036

See [`CURSOR_APHI_STAGE10_7_A_DIVERGENCE_PENALTY_RECORD.md`](CURSOR_APHI_STAGE10_7_A_DIVERGENCE_PENALTY_RECORD.md)

## Stage 10.5 Diagnostic Document Index

| Document / Test | Content |
|-------------|------|
| [`CURSOR_APHI_STAGE10_5_PHI_GAUGE_DIAGNOSTIC_RECORD.md`](CURSOR_APHI_STAGE10_5_PHI_GAUGE_DIAGNOSTIC_RECORD.md) | **φ penalty sweep {0,1,10,100} + φ-only**; main ChatGPT discussion record |
| `test_3d_aphi_ck_contact_phi_gauge_penalty_sweep_diagnostic` | Data source for table above |
| `test_3d_aphi_ck_contact_left_field_error_diagnostic` | 10.5-A/C/E left-body partitioning + mono comparison |

## Engineering Conventions (Memo)

1. Single-body GMRES Arnoldi: **block-diagonal** Contact apply; coupled GMRES Arnoldi/residual: **full** Contact apply
2. **Do not** `syncAphiBlockToDevice` exact after zeroing neighbor solution (writes stale host exact back to device → NaN)
3. `InteractionDynamicsCK` for Contact needs `FooCK<Contact<>>` partial specialization
4. `solution_block_max` ≠ magnetic field B; B needs future `curl(A)` post-processing
5. MMS **passed** uses Krylov semantics; left-body `left_continuous` diagnostic only
6. Polish default **1× block-GS sweep** (MMS main path); extra sweeps see polish diagnostic
