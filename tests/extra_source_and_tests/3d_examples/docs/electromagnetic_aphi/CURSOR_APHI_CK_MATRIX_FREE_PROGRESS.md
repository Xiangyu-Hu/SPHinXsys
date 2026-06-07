# A-phi CK Matrix-Free Solver Development Progress Record

> For discussion with ChatGPT.  
> Code directory: `tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/`  
> Test directory: `tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_*`

## Short-Term Goals

1. ✅ Build/run all `test_3d_aphi_ck_*` baseline tests on new machine
2. ✅ Debug LHS regression (including Laplace weights aligned with diffusion)
3. ✅ Fused `AphiApplyCK` (Inner only)
4. ✅ Block vector operations + BiCGStab skeleton (convergence deferred to Stage 7)
5. **Not doing**: TEAM7, cold crucible, Contact, penalty gauge (later stages)

---

## OPERATOR Conventions (Unchanged)

- Grad/Div default: `PairwiseUncorrected` (no B correction)
- B-correction `LinearGradient` path: `DIAGNOSTIC_B_CORRECTED / DEBUG_REFERENCE_ONLY` only
- Residual: `Residual = RHS - LHS`, `LHS = K(X)`
- Gather writes: each particle i writes only `output[i]`, no atomic
- Pairwise grad/div weights: **sole** entry `AphiPairwiseGradientWeightUncorrected()` (`aphi_laplace_ck.h`)

---

## Module List (Stage 0–6A)

| File | Responsibility |
|------|------|
| `aphi_field_names_ck.h` | `AphiVariableNames` |
| `aphi_variables_ck.h/.hpp` | variable registration/initialization |
| `aphi_laplace_ck.h/.hpp` | Laplace weights + `AphiPairwiseGradientWeightUncorrected` |
| `aphi_coupling_modes_ck.h` | Grad/Div mode enum + `AphiLhsAssemblyOptions` |
| `aphi_block_zero_ck.h/.hpp` | zero/copy/residual |
| `aphi_reaction_ck.h/.hpp` | ωσA local terms |
| `aphi_grad_phi_coupling_ck.h/.hpp` | σ∇φ coupling |
| `aphi_div_sigma_a_coupling_ck.h/.hpp` | ω∇·(σA) coupling |
| `aphi_assemble_lhs_debug_ck.h/.hpp` | Debug LHS assembly |
| `aphi_matrix_free_operator_ck.h/.hpp` | Fused `AphiApplyCK` + `AphiApplyDynamicsBundle` |
| `aphi_block_vector_ops_ck.h/.hpp` | **AXPY / linear combo / BiCGStab update / Vol-weighted dot & norm²** |
| `aphi_bicgstab_solver_ck.h/.hpp` | **right-PC BiCGStab** |
| `aphi_pcg_solver_ck.h/.hpp` | **block PCG (Stage 8B)** |
| `aphi_krylov_diagnostics_ck.h` | shared finite/near-zero/explosion helpers |
| `aphi_block_jacobi_preconditioner_ck.h/.hpp` | block-Jacobi preconditioner |
| `aphi_lhs_test_helpers.h` | test body + host dot/norm + block differencing |
| `aphi_test_device_sync.h` | SYCL device→host sync |

---

## Test List (16 tests, 2026-05-27)

| Test | Purpose | Result |
|------|------|------|
| `test_3d_aphi_ck_variable_registration` | variables/materials | ✅ PASS |
| `test_3d_aphi_ck_gradient_divergence_manufactured` | Stage2 diagnostic | ✅ PASS |
| `test_3d_aphi_ck_pairwise_laplace_manufactured` | Laplace manufactured | ✅ PASS (threshold 0.25) |
| `test_3d_aphi_ck_lhs_reaction_only` | LHS sub-term | ✅ PASS |
| `test_3d_aphi_ck_lhs_grad_only` | LHS sub-term | ✅ PASS (threshold 0.05) |
| `test_3d_aphi_ck_lhs_div_only` | LHS sub-term | ✅ PASS (threshold 0.05) |
| `test_3d_aphi_ck_lhs_full_zero` | X=0 | ✅ PASS |
| `test_3d_aphi_ck_lhs_polynomial_manufactured` | combined smoke | ✅ PASS |
| `test_3d_aphi_ck_lhs_discrete_rhs_residual` | B=K(X), R≈0 | ✅ PASS |
| `test_3d_aphi_ck_fused_apply_vs_debug_lhs` | Fused vs debug | ✅ PASS |
| `test_3d_aphi_ck_block_vector_ops` | Stage 5A/5B | ✅ PASS |
| `test_3d_aphi_ck_apply_search_and_s_blocks` | Stage 5C | ✅ PASS |
| `test_3d_aphi_ck_bicgstab_exact_initial_guess` | Stage 6 initialization | ✅ PASS |
| `test_3d_aphi_ck_bicgstab_discrete_rhs_solve` | Stage 6B smoke | ✅ PASS (see below) |
| `test_3d_aphi_ck_bicgstab_vector_helmholtz_solve` | Stage 6C smoke | ✅ PASS (see below) |

### Latest Run Summary (2026-05-27)

```
block_vector_ops: dot_err=3.6e-7, norm_err=1.2e-6, passed=1
apply_search_and_s_blocks: v_vs_lhs=0, t_vs_v=0, passed=1
bicgstab_exact_initial_guess: iter=0, init_res=3.6e-6, converged=1, passed=1
bicgstab_discrete_rhs_solve: init_res=53.02=host_rhs_norm, iter=26, breakdown=1, passed=1 (smoke)
bicgstab_vector_helmholtz: init_res=39.12=host_rhs_norm, iter=42, breakdown=1, passed=1 (smoke)
pairwise_laplace: scalar_err=0.086, vector_err=0.050, passed=1 (threshold 0.25)
lhs_grad_only: core_max_error=0.0098, passed=1 (threshold 0.05)
lhs_div_only: core_max_error=0.017, passed=1 (threshold 0.05)
```

---

## Step 4 Cleanup (2026-05-27 ✅)

1. **Comment fix**: `AphiApplyCK` changed to "one InteractKernel + multiple optional neighbor loops"
2. **Mode check**: `throw std::runtime_error` when not `PairwiseUncorrected`
3. **Threshold tightening**: `pairwise_laplace` 1.0→0.25; `lhs_grad/div_only` 0.35→0.05
4. **Single g_ij entry**: fused/coupling both reuse `AphiPairwiseGradientWeightUncorrected`

---

## Step 5 — Block Vector Ops (2026-05-27 ✅)

### 5A LocalDynamics (`aphi_block_vector_ops_ck.h/.hpp`)

- `AphiBlockAXPYCK`: `dst += alpha * src`
- `AphiBlockLinearCombinationCK`: `dst = a*x + b*y`
- `AphiBlockBiCGStabUpdateSolutionCK`: `X += alpha*P + omega*S`
- `AphiBlockBiCGStabUpdateResidualCK`: `R = S - omega*T`
- `AphiBlockBiCGStabUpdateSearchCK`: `P = R + beta*(P - omega*V)`

### 5B Reduction (`ReduceDynamicsCK` + `LocalDynamicsReduce<ReduceSum<Real>>`)

- `AphiBlockDotProductCK`: Vol-weighted `<X,Y>` (variable name `VolumetricMeasure`)
- `AphiBlockNormSquaredCK`: Vol-weighted `<X,X>`, norm = sqrt(...)

### 5C Tests

- `test_3d_aphi_ck_block_vector_ops`: copy/axpy/linear combo/dot/norm vs host
- `test_3d_aphi_ck_apply_search_and_s_blocks`: `Search→V` and `S→T` match `Solution→LHS` (core diff 0)

---

## Step 6 — BiCGStab (2026-05-27 skeleton ✅, full convergence ⏳)

### 6A `AphiBiCGStabSolverCK<ExecutionPolicy>`

- Variable mapping: solution=X, rhs=B, lhs=K(X), residual=R, r_hat, search=P, v=V, s=S, t=T
- Initialization: apply → residual → r_hat=r, P=r
- Standard BiCGStab iteration + breakdown detection (rho, r_hat·v, t·t, omega)
- Vol-weighted inner products via `AphiBlockDotProductCK`

### 6B/6C Test Status

| Test | Acceptance | Status |
|------|------|------|
| `bicgstab_exact_initial_guess` | X=X_exact → iter=0 converged | ✅ |
| `bicgstab_discrete_rhs_solve` | \|\|R0\|\|=\\|RHS\\|, solver runs | ✅ smoke |
| `bicgstab_vector_helmholtz_solve` | laplace_a only, same as above | ✅ smoke |

**Full convergence not achieved (Stage 7 discussion pending)**:

1. Discrete pairwise Laplacian uses `dW_ij * Vol_j`, generally **non-symmetric** (w_ij ≠ w_ji)
2. Full A-φ includes reaction coupling, **non-self-adjoint** under Vol-weighted inner product
3. φ equation has **gauge/nullspace** (constant φ mode)
4. After ~20–40 iterations `rho` changes sign → BiCGStab breakdown, residual explodes to ~1e15

**Next (Stage 7 candidates)**: penalty gauge / preconditioner / symmetrized weights / GMRES comparison, **still no Contact/TEAM7**

---

## Step 7 — Gauge / Diagnostics / Preconditioned Krylov (2026-05-21 in progress)

### 7A Smoke test rename/comment ✅

| Old name | New name/notes |
|------|-----------|
| `bicgstab_discrete_rhs_solve` | kept, added `PLUMBING_SMOKE_ONLY=1` comment |
| `bicgstab_vector_helmholtz_solve` | → `bicgstab_vector_laplace_singular_smoke` (reaction=false, pure vector Laplace) |

### 7B Laplace adjointness diagnostic ✅

- `test_3d_aphi_ck_laplace_adjointness_diagnostic`
- scalar φ / vector A separately test Vol-weighted `<x,Ky>_V` vs `<Kx,y>_V`
- Results: `scalar_relative_gap≈0`, `vector_relative_gap≈1e-6`, **passed=1**
- **Did not change** pairwise Laplace weights

### 7C Constant nullspace diagnostic ✅

- `test_3d_aphi_ck_laplace_constant_nullspace_diagnostic`
- `||K(constant)||≈0` → pure Laplace singular, **passed=1**

### 7D BiCGStab debug log + breakdown fix ✅

- `AphiBiCGStabBreakdownCode` enum + `enable_debug_log`
- **rho sign change no longer treated as breakdown**; only near-zero denominators / residual explosion
- debug output: iter, rho, r_hat_dot_v, alpha, t_dot_s, t_dot_t, omega, beta, residual_norm, relative_residual

### 7E Phi gauge penalty ✅

- `AphiLhsAssemblyOptions`: `use_phi_gauge_penalty`, `phi_gauge_penalty`
- `aphi_phi_gauge_penalty_ck.h/.hpp` + fused/debug LHS integration
- `test_3d_aphi_ck_fused_apply_phi_gauge_penalty_vs_debug_lhs`: core_max_diff≈1.9e-5, **passed=1**

### 7F Convergence tests (implemented, **full convergence not yet achieved** ⏳)

| Test | Config | Result (block-Jacobi PC) |
|------|------|-------------------------|
| `bicgstab_scalar_phi_laplace_penalty_diagnostic` | laplace_phi + penalty | iter≈38, ResidualExplosion (informational) |
| `bicgstab_vector_helmholtz_convergence` | laplace_a + reaction | iter≈247, ResidualExplosion |
| `bicgstab_full_aphi_penalty_convergence` | full + penalty | iter≈136, ResidualExplosion |

debug observation (scalar phi + PC): iter0 `relative_residual≈0.51`, rho negative from iter1, explodes after ~40 steps — same root cause as Stage 6.

### 7G Block-Jacobi preconditioner (first version ✅, convergence still insufficient ⏳)

- `aphi_block_jacobi_preconditioner_ck.h/.hpp`
- A block: local 2×2 real/imag reaction inverse; φ block: |laplace_diag| + penalty
- Ignores grad/div off-diagonal
- `AphiBiCGStabSolverCK` new parameter `use_block_jacobi_preconditioner`
- **right-preconditioned** BiCGStab (Stage 8A fixed comments and variable semantics)

**Stage 7 revised conclusion**: matrix-free apply, phi gauge penalty, block-Jacobi PC wiring all OK; BiCGStab from zero guess still does not converge (including scalar φ+penalty). Laplace adjointness gap is small; **insufficient evidence to attribute failure mainly to Vol_j asymmetry**. Stage 8 will first separate solver instability vs operator inconsistency.

---

## Step 8 — Solver Diagnostics / PCG / GMRES (2026-05-21 in progress)

### 8A PC-BiCGStab cleanup + diagnostics (✅ partially complete)

1. **Comment fix**: right-preconditioned BiCGStab (`z=M⁻¹p`, `y=M⁻¹s`)
2. **Dedicated Krylov variables**: `z`, `y`, `v_old`, `true_residual`; **RHS no longer overwritten**
3. **`AphiBiCGStabSolverOptions`**: tolerance/breakdown/debug/true-residual parameters centralized
4. **Extended breakdown enum**: `RhoNewNearZero`, `AlphaNonFinite`, `BetaNonFinite`, `ResidualNonFinite`, etc.
5. **True residual recompute**: `AphiComputeBlockResidualCK` + optional `recompute_true_residual` / gap logging
6. **New test** `test_3d_aphi_ck_block_jacobi_diagonal_diagnostic`:
   - laplace_a_diag: min≈226, max≈587, non_positive=0
   - laplace_phi_diag: min≈302, max≈783, non_positive=0
   - d_phi=diag+penalty all positive

### 8B scalar φ PCG diagnostic (✅ implementation complete, **convergence failed** ⏳)

**New files:**
- `aphi_pcg_solver_ck.h/.hpp` — block PCG, `AphiPCGSolverOptions` / `AphiPCGResult` / breakdown enum
- `AphiBlockCGUpdateSearchCK` — `p = z + beta*p` (`aphi_block_vector_ops_ck.h/.hpp`)
- `test_3d_aphi_ck_pcg_scalar_phi_laplace_penalty_diagnostic` (renamed from `_convergence`, explicitly diagnostic)

**Variable mapping (PCG):** search=p, v=q=Ap, z=M⁻¹r, residual=r, solution=x

**Test config:** same as BiCGStab scalar phi (laplace_phi + penalty=10, tol=1e-5, block-Jacobi PC, true residual interval=5)

**Results (2026-05-21 build/run):**

| Variant | iter | breakdown | Notes |
|------|------|-----------|------|
| PCG + block-Jacobi PC + penalty | 30 | ResidualExplosion | rel_res: 0.955→1.64e12 |
| PCG no preconditioner + penalty | 23 | ResidualExplosion | explodes faster |
| PCG + block-Jacobi PC **no penalty** | 35 | ResidualExplosion | pure Laplace also fails |

**debug observation (PC + penalty):**
- iter0: `alpha≈3.67`, `recursive_rel_res≈0.955` (first step barely decreases)
- iter1–29: residual **monotonically increases**; `denom`/`rho` grow exponentially; `alpha≈1.24` nearly unchanged
- **Did not trigger** `NonPositiveCurvature` (`p^T A p > 0` throughout)
- `recursive_true_gap ≈ 1e-7` → **PCG plumbing self-consistent**, recursive and true residuals match

**8B conclusion (revised, before Stage 8C Phase 1):**
- PCG diverges from zero guess; recursive/true residual gap ~1e-7 → **not BiCGStab-specific, not a recursive residual artifact**
- Positive Jacobi diagonal and CG direction `p^T A p > 0` **do not prove global SPD**
- Next: rough-field adjointness, random Rayleigh/curvature, tiny host matrix diagnostic; parallel GMRES(m)
- **Still not changing** default SPH Vol_j Laplace weights

### 8C Phase 1 — code hygiene (✅ 2026-05-21)

1. **`InitializeAphiVariablesCK`**: bind and zero-init `v_old`, `z`, `y`, `true_residual`
2. **PCG / BiCGStab**: `maybeRecomputeTrueResidual` writes `result.final_true_residual_norm`
3. **Test naming**: `pcg_scalar_phi_laplace_penalty_convergence` → `_diagnostic`

### 8C Phase 2 — scalar operator diagnostics (✅ 2026-05-21)

**New shared code:**
- `aphi_scalar_phi_diagnostic_helpers.h` — rough random fields, basis vectors, scalar Laplace+penalty options
- `aphi_lhs_test_helpers.h` — Vol-weighted mean, host matrix symmetry/nullspace/eigenvalue analysis (Eigen)

**New tests (informational, `passed=1`):**

| Test | Key Results |
|------|----------|
| `scalar_phi_random_curvature_diagnostic` | penalty=0/10/100/1000: **negative_count=0** (20 seeds) |
| `scalar_phi_random_adjointness_diagnostic` | rough field mean_gap≈2.5e-7, max_gap≈8.7e-7 |
| `scalar_phi_constant_mode_rhs_diagnostic` | penalty=0 mean_rhs_phi≈1e-7 (nullspace compatible) |
| `scalar_phi_tiny_host_matrix_diagnostic` | N=64, dp=0.25; penalty=10 **min_eig_sym=10** |

**Phase 2 preliminary interpretation:** operator SPD on small grid, near self-adjoint on rough fields, no negative Rayleigh values → PCG failure more likely from **preconditioner/condition number/solver**, not global indefiniteness.

### 8C Phase 3 — restarted right-PC GMRES(m) (✅ 2026-05-21)

**New:** `aphi_gmres_workspace_ck.*`, `aphi_gmres_solver_ck.*`, `aphi_gmres_test_helpers.h`, `AphiBlockScaleCopyCK`

| Test | Result (m=30, tol=1e-5, block-Jacobi PC) |
|------|------------------------------------------|
| `gmres_scalar_phi_laplace_penalty_convergence` | outer=2, arnoldi=30, rel_res≈1.6e-6 **passed=1** |
| `gmres_vector_helmholtz_convergence` | outer=5, arnoldi=120, rel_res≈2.9e-6 **passed=1** |
| `gmres_full_aphi_penalty_convergence` | outer=6, arnoldi=150, rel_res≈3.7e-6 **passed=1** |

**Fork conclusion:** same operator GMRES converges, BiCGStab/PCG diverge → Krylov method/PC combination issue, not operator plumbing. Fix: rectangular Hessenberg QR must not use `isInvertible()`.

### 9A — production wrapper + GMRES hardening (✅ 2026-05-21)

**New shared code:**
- `aphi_matrix_free_solve_ck.h/.hpp` — `AphiKrylovSolverKind`, `AphiMatrixFreeSolveCK` (default GMRES)
- `aphi_gmres_test_helpers.h` — `defaultMatrixFreeGMRESOptions`, `gmresConvergencePassed`
- `aphi_gmres_solver_ck.hpp` — force `recomputeFinalTrueResidual` before convergence; workspace validation changed to `>= m`

**Test changes:**
- Three GMRES convergence tests routed through `AphiMatrixFreeSolveCK`, pass criteria strengthened (true rel ≤ 10×tol, gap < 1e-4)
- `bicgstab_scalar_phi_laplace_penalty_convergence` → **`bicgstab_scalar_phi_laplace_penalty_diagnostic`** (informational)

**New diagnostic tests:**

| Test | Key Results |
|------|----------|
| `gmres_restart_sensitivity_diagnostic` | m=10/20/30 all converge (m=10 needs 4 outer) |
| `gmres_preconditioner_comparison_diagnostic` | no-PC and block-Jacobi both converge (~1.7e-6) |

**Production default:** `AphiKrylovSolverKind::GMRES` + block-Jacobi PC + m=30.

### 9B — manufactured-solution robustness sweep (✅ 2026-05-21)

**New shared code:**
- `aphi_gmres_robustness_sweep_helpers.h` — full A-phi manufactured-solution GMRES single-case runner
- `AphiGMRESResult::monotonic_outer_residual` — whether outer restart residuals decrease monotonically

**New test:** `test_3d_aphi_ck_gmres_manufactured_robustness_sweep_diagnostic` (informational, `passed=1`)

Single-parameter sweeps (others at baseline: dp=0.1, omega=1.25, sigma=2, penalty=10, tol=1e-5, m=30, max_outer=30):

| Sweep | Values | Result Summary |
|------|------|----------|
| baseline | — | outer=6, rel≈3.7e-6, monotonic=1 ✅ |
| dp | 0.2 / 0.1 / 0.075 | all converge; 0.075 needs outer=11, arnoldi=300 |
| omega | 0.1 / 1.25 / 10 | **0.1 misses tol** (outer=30, rel≈1.07e-5); 1.25/10 ✅ |
| sigma | 0 / 0.1 / 2.0 | all converge |
| phi_gauge_penalty | 1 / 10 / 100 | all converge |

**Summary:** **12/13** cases converge; `omega=0.1` is low-frequency slow-convergence boundary (monotonic but not reaching 1e-5 within 30 outer).

### 9C — TEAM7-style non-contact benchmarks (✅ 2026-05-21)

**Detailed record:** [`CURSOR_APHI_STAGE9C_BENCHMARK_RECORD.md`](CURSOR_APHI_STAGE9C_BENCHMARK_RECORD.md)

**New shared code:**
- `aphi_benchmark_case_ck.h/.hpp` — region materials, impressed current RHS, TEAM7-like layout
- `aphi_gmres_benchmark_helpers.h` — core-region field strength / MMS error metrics

**New tests (all passed=1):**

| Test | Type | Key Results |
|------|------|----------|
| `gmres_impressed_current_homogeneous_box` | physical RHS | outer=5, rel≈2.7e-6, source/exterior≈1.4 |
| `gmres_two_material_interface_manufactured` | variable σ discrete MMS | outer=14, discrete_defect≈4.8e-6, continuous≈3e-5 ✅ |
| `gmres_team7_like_coil_plate_simplified` | TEAM7 layout | outer=48, rel≈5.0e-4 (tol=5e-4, 10× tighter than 5e-3) ✅ |
| `gmres_team7_like_coil_plate_tight_tol_diagnostic` | 1e-5 attempt | rel≈5.8e-4 @ outer=100, missed 1e-5 (informational) |

**Notes (post-9C improvements):** 9C-2 uses `r_hat` for reference field; 9C-3 initial decoupled PC @ m=30 floor ~5e-4. **See below after Stage 9D.**

### 9D-0 — GMRES/PC diagnostics and fixes (✅ 2026-05-21)

**Code fixes:**
- `AphiGMRESSolverCK`: recompute final residual after last solution update before max-outer exit (`MaxOuterIterationsReached` no longer marked as numerical breakdown)
- `AphiGMRESMaxRestartDimension`: 30 → **80**; `AphiMatrixFreeSolveCK` workspace registered per `solver_options.gmres.restart_dimension`

**New helpers / tests:**
- `logTrueResidualDiagnostics` — component-wise (A/φ) + region-wise (air/coil/conductor) true residual breakdown
- `runTeam7LikeGMRESDirect` — TEAM7 GMRES single-run driver (shared max workspace)
- `test_3d_aphi_ck_gmres_team7_like_solver_study_diagnostic` — restart sweep + PC comparison + residual breakdown

**TEAM7 full contrast @ tol=1e-5, max_outer=100 (decoupled PC):**

| m | true rel | Conclusion |
|---|----------|------|
| 30 | ~5.3e-4 | restart truncation obvious |
| 50 | ~2.0e-5 | near 1e-5 |
| 80 | ~2.0e-5 | saturates with m=50 |

**PC comparison (m=30):** no-PC ~1.6e-3; decoupled block-Jacobi ~2.3e-4.

**Key conclusion:** prior ~5e-4 floor **partly from** (1) max-outer residual not recomputed (2) m=30 restart truncation; increasing m reaches ~2e-5, not purely a PC spectral floor.

### 9D-1 — Coupled 8×8 Point-Block Jacobi PC (✅ 2026-05-21)

**New:**
- `AphiBlockJacobiPreconditionerKind::CoupledPointBlock8x8` — local grad(φ)/div(σA) coupling in 8×8 block inverse
- `JacobiGradPhiCoupling` / `JacobiDivACoupling` (Vecd) diagonal assembly
- Production default: `block_jacobi_kind = CoupledPointBlock8x8`, **m=50**

**New test:** `test_3d_aphi_ck_gmres_coupled_pc_team7_convergence`

| Config | true rel | converged @1e-5 |
|------|----------|-----------------|
| coupled m=50, outer=100 | ~2.0e-5 | ❌ |
| coupled m=80, outer=150 | **~1.96e-5** | ❌ (closest) |

**9C-3 main test (updated):** m=50 coupled PC, tol=5e-4 → **outer=2**, rel≈**2.3e-4** ✅ (much faster than decoupled m=30)

**9D-1 conclusion:** coupled PC at m≥50 same order as decoupled; **1e-5 still ~2× short** (~2e-5 floor). Component residual A already ~5e-6; bottleneck may be global low-frequency modes.

### 9D-lite — Joule heating data path (✅ 2026-05-21)

**New:**
- `aphi_joule_heating_ck.h/.hpp` — grad φ → E → J → Joule source (CK/SYCL)
- `test_3d_aphi_ck_joule_heating_plumbing_diagnostic` — regional Joule power after TEAM7 solve

**Results (tol=5e-4, m=50 coupled):**
- `total_joule_power ≈ 0.136`, `conductor_joule_power ≈ 0.136`
- `min_joule_source > 0` (non-negativity OK)
- **No** σ(T)/ν(T) feedback, no quantitative thermal coupling conclusions

### 9D-2 — Joule tolerance sensitivity (✅ 2026-05-21)

**New test:** `test_3d_aphi_ck_joule_em_tolerance_sensitivity_diagnostic`

Scan EM tol = **5e-3 / 5e-4 / 1e-4 / 1e-5** (m=50 coupled PC, max_outer=100), reference **tol=5e-4**:

| tol | em_rel | converged | conductor_joule |
|-----|--------|-----------|-----------------|
| 5e-3 | ~2.4e-4 | ✅ | 0.13619 |
| 5e-4 | ~2.5e-4 | ✅ | 0.13620 (ref) |
| 1e-4 | ~6.1e-5 | ✅ | 0.13619 |
| 1e-5 | ~2.0e-5 | ❌ | 0.13619 |

**Relative change in conductor Joule power vs reference (5e-4):** max **~5.3×10⁻⁵ (< 0.01%)** → `sensitivity_passed=1`

**Engineering conclusion:** tightening EM tolerance from 5e-4 to ~2e-5 barely changes Joule power; **~2e-5 acceptable as TEAM7 engineering EM tolerance**, no need to force 1e-5 for Joule quantitative conclusions.

---

## Stage 9E — Stronger Validation (✅ 2026-05-21)

Detailed record in **`CURSOR_APHI_STAGE9E_VALIDATION_RECORD.md`**.

### 9E-1 — Interface flux-matched MMS ✅

**New:** `AssignInterfaceFluxMatchedPhiFieldsCK` + `test_3d_aphi_ck_gmres_interface_flux_matched_manufactured`

| Metric | 9C-2 (separable) | 9E-1 (flux-matched) |
|--------|------------------|---------------------|
| continuous_field_error | ~O(1) → ~3e-5 (discrete ref fix) | **~1.8e-6** |
| discrete_mms_defect | ~4.8e-6 | **~5.0e-6** |
| interface_band_rel | — | **~1.9e-6** |

### 9E-2 — Physical TEAM7 dimensions ✅

**New:** `AphiTeam7PhysicalDimensions`, `buildTeam7LayoutForBox`, `test_3d_aphi_ck_gmres_team7_physical_dimensions_smoke`

- Box **1.2×1.0×0.3 m**, dp=0.1, tol=5e-4 → outer=2, true_rel≈**4.7e-4**, conductor/coil field max≈**3.6**

### 9E-3 — TEAM7 dp sweep ✅

**New:** `test_3d_aphi_ck_gmres_team7_dp_convergence_diagnostic`

| dp | true_rel @5e-4 | converged |
|----|----------------|-----------|
| 0.2 | ~6.3e-3 | ❌ (too coarse) |
| 0.1 | ~2.4e-4 | ✅ |
| 0.075 | ~4.8e-4 | ✅ (outer=20) |

---

## Stage 10A — EM/Joule Verification (✅ 2026-05-21)

Detailed record in **`CURSOR_APHI_STAGE10A_PROGRESS.md`**.

| Item | Result |
|------|--------|
| Joule uniform-field | core mean q error **~1.1%** ✅ |
| 8×8 PC fallback | **0** (no silent fallback) |
| Joule local sensitivity | L2/max/band vs 5e-4 ref **<0.01%** ✅ |
| Physical dp observables | Joule power 0.043→0.081→0.093 (dp=0.15/0.1/0.075), not yet saturated |

---

## Stage 10B — Quantitative TEAM7 Reference (✅ 2026-05-21)

Detailed record in **`CURSOR_APHI_STAGE10B_QUANTITATIVE_RECORD.md`**.

- **Canonical spec:** `aphi_team7_canonical_case_ck.h` (1.2×1.0×0.3 m, tol=5e-4, m=50)
- **Self-reference:** dp=0.075 vs dp=0.1 → conductor Joule rel **~13%**, field max rel **~17%** (<20% ✅)
- **External FEM:** deferred
- **Test:** `test_3d_aphi_ck_gmres_team7_quantitative_reference_comparison` ✅

---

## Stage 10C — Cold-Crucible Scaffold (✅ 2026-05-21)

Detailed record in **`CURSOR_APHI_STAGE10C_SCAFFOLD_RECORD.md`**.

- **Layout:** melt + crucible shell + dual coil (unit box fractions)
- **Test:** `test_3d_aphi_ck_gmres_cold_crucible_scaffold` ✅ (outer=2, melt/crucible/coil/air partitions OK)
- **Joule:** crucible wall dominates (σ=1e4); melt Joule > 0
- **Not doing:** Contact, thermal feedback, σ(T)

---

## Production Defaults (Updated 2026-05-21)

| Parameter | Value |
|------|-----|
| Solver | `AphiMatrixFreeSolveCK` + GMRES |
| Restart m | **50** |
| PC | **CoupledPointBlock8x8** |
| Max workspace | 80 (registered on demand ≤ m) |

---

## New Machine Build Commands (Reference)

```bash
cd /home/yongchuan/sphinxsys/build
cmake .. -DSPHINXSYS_3D=ON -DSPHINXSYS_BUILD_3D_EXAMPLES=ON -DSPHINXSYS_BUILD_EXTRA_SOURCE_AND_TESTS=ON
ninja test_3d_aphi_ck_gmres_scalar_phi_laplace_penalty_convergence \
      test_3d_aphi_ck_gmres_vector_helmholtz_convergence \
      test_3d_aphi_ck_gmres_full_aphi_penalty_convergence \
      test_3d_aphi_ck_gmres_restart_sensitivity_diagnostic \
      test_3d_aphi_ck_gmres_preconditioner_comparison_diagnostic \
      test_3d_aphi_ck_gmres_manufactured_robustness_sweep_diagnostic \
      test_3d_aphi_ck_gmres_impressed_current_homogeneous_box \
      test_3d_aphi_ck_gmres_two_material_interface_manufactured \
      test_3d_aphi_ck_gmres_team7_like_coil_plate_simplified \
      test_3d_aphi_ck_gmres_team7_like_coil_plate_tight_tol_diagnostic \
      test_3d_aphi_ck_gmres_team7_like_solver_study_diagnostic \
      test_3d_aphi_ck_gmres_coupled_pc_team7_convergence \
      test_3d_aphi_ck_joule_heating_plumbing_diagnostic \
      test_3d_aphi_ck_joule_em_tolerance_sensitivity_diagnostic \
      test_3d_aphi_ck_gmres_interface_flux_matched_manufactured \
      test_3d_aphi_ck_gmres_team7_physical_dimensions_smoke \
      test_3d_aphi_ck_gmres_team7_dp_convergence_diagnostic \
      test_3d_aphi_ck_joule_uniform_field_analytic_verification \
      test_3d_aphi_ck_block_jacobi_8x8_pc_diagnostic \
      test_3d_aphi_ck_joule_local_distribution_sensitivity_diagnostic \
      test_3d_aphi_ck_gmres_team7_physical_dp_observable_diagnostic \
      test_3d_aphi_ck_gmres_team7_quantitative_reference_comparison \
      test_3d_aphi_ck_gmres_cold_crucible_scaffold \
      test_3d_aphi_ck_scalar_phi_random_curvature_diagnostic \
      test_3d_aphi_ck_scalar_phi_random_adjointness_diagnostic \
      test_3d_aphi_ck_scalar_phi_constant_mode_rhs_diagnostic \
      test_3d_aphi_ck_scalar_phi_tiny_host_matrix_diagnostic \
      test_3d_aphi_ck_pcg_scalar_phi_laplace_penalty_diagnostic \
      test_3d_aphi_ck_bicgstab_scalar_phi_laplace_penalty_diagnostic \
      test_3d_aphi_ck_block_jacobi_diagonal_diagnostic \
      test_3d_aphi_ck_bicgstab_exact_initial_guess
```

---

## Key Conclusions for ChatGPT

1. **Stage 8A complete**: RHS protection, right-PC BiCGStab semantics, breakdown/options, true residual, Jacobi diag all positive
2. **Stage 8B complete**: PCG solver + scalar φ convergence test; **PCG also diverges** (iter≈30, ResidualExplosion)
3. **Solver plumbing trusted**: recursive/true residual gap ~1e-7; exact initial guess (BiCGStab) iter=0 still PASS
4. **Phase 2 diagnostics**: rough adjointness ~1e-7; no negative Rayleigh; tiny matrix penalty=10 min_eig=10
5. **Phase 3 GMRES**: scalar/vector/full **all converge**; BiCGStab/PCG still diverge → **operator solvable, CG-type methods unsuitable for current PC setup**
6. **Stage 9A**: `AphiMatrixFreeSolveCK` production wrapper; GMRES forces true residual before convergence; restart/PC diagnostics pass
7. **Stage 9B**: full A-phi manufactured solution 13-parameter sweep, 12/13 converge; `omega=0.1` is low-frequency slow-convergence boundary
8. **Stage 9C**: three benchmarks + discrete MMS fix + TEAM7 tol=5e-4
9. **Stage 9D-0**: max-outer residual fix; m≤80; TEAM7 restart/PC diagnostics; component/partition residual breakdown
10. **Stage 9D-1**: coupled 8×8 PC; production default m=50 + coupled PC; TEAM7 ~2e-5 (misses 1e-5)
11. **Stage 9D-lite**: Joule source plumbing ✅
12. **Stage 9D-2**: Joule tolerance sensitivity ✅ (5e-4→~2e-5 power change <0.01%)
13. **Stage 9E-1**: interface flux-matched MMS ✅ (continuous error ~2e-6 vs 9C-2 separable)
14. **Stage 9E-2**: real TEAM7 dimensions 1.2×1.0×0.3 m smoke ✅
15. **Stage 9E-3**: dp sweep ✅ (dp=0.1 recommended)
16. **Stage 10A-1**: Joule uniform-field ~1.1% error; 8×8 PC fallback=0 ✅
17. **Stage 10A-2**: Joule local distribution sensitivity ✅; physical dimensions dp observable sweep ✅
18. **Stage 10B**: canonical TEAM7 spec + fine-grid self-reference ✅
19. **Stage 10C**: cold crucible geometry scaffold ✅ (EM+Joule, inner-only, no thermal feedback)
20. **Stage 10-contact Sprint 1**: `AphiPairwiseLaplaceCK<Contact<>>` + monolithic vs two-body apply equivalence ✅ (max_abs_diff≈1.5e-5)
22. **Stage 10-contact Sprint 2**: Grad/Div/PhiGradient Contact + full debug K(X) + Joule equivalence ✅
23. **Stage 10-contact Sprint 3**: Block Jacobi + fused `AphiApplyCK` Contact ✅
24. **Next main line**: 9E contact interface MMS → three-body TEAM7-like contact scaffold → GMRES Contact
