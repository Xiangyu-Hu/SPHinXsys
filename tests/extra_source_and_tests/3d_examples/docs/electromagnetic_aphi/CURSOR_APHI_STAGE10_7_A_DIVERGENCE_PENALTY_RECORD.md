# Stage 10.7 — A-divergence penalty record

> **Date**: 2026-05-28 ~ 2026-05-29  
> **Status**: **Inner prototype phase closed (research prototype, not production)**  
> **Basis**: ChatGPT [`stage10_contact_a_gauge_hourglass_curlA_relaxed_tests_for_cursor.md`](../../../../../../docs/electromagnetic_aphi/stage10_contact_a_gauge_hourglass_curlA_relaxed_tests_for_cursor.md)  
> **Related**: [`CURSOR_APHI_STAGE10_6_A_GAUGE_DIVERGENCE_DIAGNOSTIC_RECORD.md`](CURSOR_APHI_STAGE10_6_A_GAUGE_DIVERGENCE_DIAGNOSTIC_RECORD.md)  
> **Teacher discussion**: see [`stage10_teacher_discussion_demo_and_source_cleanup_plan.md`](../../../../../../docs/electromagnetic_aphi/stage10_teacher_discussion_demo_and_source_cleanup_plan.md) — **not the main demo line for this stage**

---

## 0. Executive summary (closure statement)

### 0.1 Goals

Prototype numerical Coulomb regularization `LhsA -= λ_A · grad(div A)` on Inner single-body; verify operator sign, fused/debug equivalence, source terms, solenoidal alternative, 3×3 graddiv PC, and whether GMRES can converge at tol=1e-4.

### 0.2 Closed items

| Item | Conclusion |
|----|------|
| Operator sign | `LhsA -= λ_A·grad(divA)` correct; sign/energy diagnostic passed |
| fused vs debug | `core_max_diff≈3e-5` passed |
| impressed source | `source_div_J_relative≈4.88`, non-solenoidal, not suitable as primary acceptance |
| solenoidal source | `AssignSolenoidalCurlCurrentRhsCK`; discrete divJ≈2.67 |
| scalar PC | `λ_A·laplace_a_diag` insufficient |
| **3×3 graddiv PC** | `JacobiGradDivABlock`; true_rel@λ=100: 0.6→**0.07** |
| GMRES outer/restart/tol | **cannot break ~4% rel floor** (impressed) |
| divA (stagnated solution) | λ=1000 impressed: div_A≈**0.0027 good** |

### 0.3 Open items (explicitly deferred)

- penalty path **converged=0** under strict tol=1e-4 (except baseline)
- Contact-side A-penalty integration
- B-corrected graddiv PC (current pairwise vs pipeline FD diff≈0.17)
- hourglass / relaxed sphere / projection / cold crucible

### 0.4 Code list (Stage 10.7 new/modified)

**Production-adjacent (in Block Jacobi, but penalty still prototype)**

| File | Role |
|------|------|
| `aphi_coupling_modes_ck.h` | `use_a_divergence_penalty`, `a_divergence_penalty` |
| `aphi_a_divergence_penalty_ck.h/.hpp` | `AphiGradDivAPenaltyCK` |
| `aphi_a_divergence_penalty_pipeline.h` | B-corrected grad→div→grad pipeline |
| `aphi_matrix_free_operator_ck.*` | fused apply append penalty |
| `aphi_assemble_lhs_debug_ck.*` | debug path penalty |
| `aphi_block_jacobi_preconditioner_ck.*` | **`JacobiGradDivABlock` 3×3** + 8×8 assembly |
| `aphi_benchmark_case_ck.*` | `AssignSolenoidalCurlCurrentRhsCK` |
| `aphi_inner_a_divergence_penalty_sweep_helpers.h` | sweep + restart/outer params |
| `aphi_graddiv_block_diagnostic_helpers.h` | FD comparison helper |

**Diagnostic-only**

| File | Role |
|------|------|
| `aphi_a_gauge_diagnostic_helpers.h` | sign energy, source divJ |
| `aphi_curl_a_diagnostic_ck.*` | B=curlA manufactured |
| `aphi_div_a_diagnostic_helpers.h` | post-solve divA |

### 0.5 Test list (Stage 10.7, 11 total)

| Test | Type | passed / purpose |
|------|------|----------------|
| `test_3d_aphi_ck_inner_a_divergence_penalty_fused_vs_debug` | **gate** | passed=1 |
| `test_3d_aphi_ck_inner_a_divergence_penalty_sign_energy_diagnostic` | gate | passed=1 |
| `test_3d_aphi_ck_impressed_current_divergence_diagnostic` | diagnostic | passed=1 |
| `test_3d_aphi_ck_solenoidal_current_divergence_diagnostic` | diagnostic | passed=1 |
| `test_3d_aphi_ck_curl_a_manufactured_diagnostic` | diagnostic | passed=1 |
| `test_3d_aphi_ck_graddiv_block_diagonal_diagnostic` | diagnostic | passed=1 |
| `test_3d_aphi_ck_inner_a_divergence_penalty_sweep_diagnostic` | sweep | diagnostic |
| `test_3d_aphi_ck_inner_solenoidal_a_divergence_penalty_sweep_diagnostic` | sweep | diagnostic |
| `test_3d_aphi_ck_inner_a_divergence_penalty_gmres_sensitivity_diagnostic` | sweep | 0/30 converged |
| `test_3d_aphi_ck_inner_a_divergence_penalty_tol_closure_diagnostic` | sweep | 0/13 converged |

**Teacher discussion**: above **none are afternoon demo cases**; if asked, explain Inner research prototype, Contact mainline already usable.

---

## 1. Mechanism

Numerical Coulomb regularization (non-physical term):

```text
LhsA -= λ_A · grad(div A)    ← corrected (2026-05-28)
```

---

## 2. ChatGPT first batch diagnostics (hourglass / relaxed deferred)

| Test | Status | Key results |
|------|------|----------|
| `test_3d_aphi_ck_inner_a_divergence_penalty_sign_energy_diagnostic` | **passed=1** | after correction `inner_A·(-grad(divA)) > 0` |
| `test_3d_aphi_ck_impressed_current_divergence_diagnostic` | **passed=1** | `source_div_J_relative≈4.88` → box source non-solenoidal |
| `test_3d_aphi_ck_curl_a_manufactured_diagnostic` | **passed=1** | `b_rel_err≈1.5e-7`, `divA≈1.1e-8` |
| `test_3d_aphi_ck_inner_a_divergence_penalty_fused_vs_debug` | **passed=1** | after sign fix `core_max_diff≈3e-5` |

---

## 3. Inner λ_A sweep after sign fix (impressed-current box, λ_φ=100)

| λ_A | converged | outer | true_rel | div_A_relative | source_div_A_relative |
|-----|-----------|-------|----------|----------------|----------------------|
| off | ✓ | 2 | 3.8e-6 | **0.578** | 0.577 |
| 10 | ✗ | 100 | **0.063** | 0.165 | 0.193 |
| 50 | ✗ | 100 | 0.354 | **0.072** | 0.098 |
| 100 | ✗ | 100 | 0.599 | 0.089 | 0.128 |
| 300 | ✗ | 100 | 0.776 | 0.152 | 0.196 |
| 1000 | ✗ | 100 | 0.811 | 0.216 | 0.261 |

**Compared to wrong sign (+grad(divA))**: at λ=10 true_rel~0.37 → after correction **~0.063** (much improved but still not converged).

**Interpretation**:
- Sign correction direction correct, Krylov condition significantly improved
- **Still no converged rows** (except baseline) → PC / source still bottleneck
- divA at λ=50–100 can drop to ~0.07–0.09 (warn), but **unconverged solution not trustworthy**
- impressed-current box source div J large, not suitable as A-penalty primary acceptance

---

## 4. Solenoidal source (`J=curl C`, 2026-05-28)

**Implementation**: `AssignSolenoidalCurlCurrentRhsCK` — `C=(0,0,0.5·g·(x²+y²))`, `g=smoothBoxProfile`, outside box `g=∇g=0`.

| Test | Status | Key results |
|------|------|----------|
| `test_3d_aphi_ck_solenoidal_current_divergence_diagnostic` | **passed=1** | `source_div_J_relative≈2.67` (impressed ~4.88) |
| `test_3d_aphi_ck_inner_solenoidal_a_divergence_penalty_sweep_diagnostic` | diagnostic | see table below |

### Inner λ_A sweep (solenoidal, λ_φ=100)

| λ_A | converged | outer | true_rel | div_A_relative | source_div_A_relative |
|-----|-----------|-------|----------|----------------|----------------------|
| off | ✓ | 2 | 1.0e-5 | **0.401** | 0.256 |
| 10 | ✗ | 100 | **0.075** | 0.120 | 0.075 |
| 50 | ✗ | 100 | 0.468 | **0.087** | 0.058 |
| 100 | ✗ | 100 | 0.686 | 0.094 | 0.065 |
| 300 | ✗ | 100 | 0.879 | 0.120 | 0.081 |

**Interpretation**:
- solenoidal source baseline **can converge**; SPH discrete `div J` still ~2.7 (analytically solenoidal, discrete metric limited)
- after adding penalty Krylov behavior similar to impressed (λ=10 true_rel~0.075), **still not converged**
- on unconverged rows divA can drop to ~0.07–0.12, but not final acceptance
- **Bottleneck confirmed in PC** (scalar `λ_A·laplace_a_diag` insufficient), not impressed source alone

---

## 5. graddiv 3×3 block-diagonal PC (2026-05-28)

**Implementation**: `JacobiGradDivABlock` (3×3 Matd), pairwise formula  
`B[c,α]=Σ_j g_ij,c (g_ji,α + S_α)`, incorporated into 8×8 CoupledPointBlock8x8 and decoupled 3×3 fallback.

| Test | Status | Key results |
|------|------|----------|
| `test_3d_aphi_ck_graddiv_block_diagonal_diagnostic` | **passed=1** | `core_max_diff≈0.17` vs B-corrected pipeline FD; `penalty_to_laplace_diag_ratio≈1.48` |

### Inner λ_A sweep (after 3×3 PC, impressed source, λ_φ=100)

| λ_A | converged | true_rel (old scalar PC) | true_rel (**3×3 PC**) | div_A_relative |
|-----|-----------|--------------------------|------------------------|----------------|
| 10 | ✗ | 0.063 | **0.027** | 0.155 |
| 50 | ✗ | 0.355 | **0.093** | **0.049** |
| 100 | ✗ | 0.600 | **0.070** | **0.024** |
| 300 | ✗ | 0.776 | **0.052** | **0.0087** (good) |

### solenoidal source (3×3 PC)

| λ_A | true_rel (old) | true_rel (**3×3 PC**) | div_A_relative |
|-----|----------------|------------------------|----------------|
| 10 | 0.075 | **0.033** | 0.109 |
| 50 | 0.468 | **0.101** | **0.035** |
| 100 | 0.686 | **0.091** | **0.018** |
| 300 | 0.879 | **0.094** | **0.0063** (good) |

**Interpretation**:
- 3×3 PC reduces true_rel at λ=100 from ~0.6 to ~0.07 (impressed), direction correct
- divA at λ=100–300 can reach warn/good, but **still not converged at tol=1e-4** (true_rel~0.05–0.09)
- Next step: slightly increase outer/restart or tol sensitivity; evaluate divA observable after convergence

---

## 6. GMRES outer/restart sensitivity (3×3 PC, tol=1e-4, 2026-05-29)

**Test**: `test_3d_aphi_ck_inner_a_divergence_penalty_gmres_sensitivity_diagnostic`

**Results**: **converged_rows=0** (30 rows; restart=100 breakdown due to `AphiGMRESMaxRestartDimension=80`)

| Source | Best config | rel | div_A |
|----|----------|-----|-------|
| impressed | λ=300 restart=80 outer=200 | **0.042** | 0.0085 good |
| solenoidal | λ=300 restart=80 outer=200 | 0.083 | 0.0064 good |

Increasing outer 100→300 **cannot** break rel floor.

---

## 7. tolerance / high λ_A closure (restart=80 outer=200, 2026-05-29)

**Test**: `test_3d_aphi_ck_inner_a_divergence_penalty_tol_closure_diagnostic` → **0/13 converged**

| Source | λ_A | rel (same across tol) | div_A |
|----|-----|-------------------|-------|
| impressed | 1000 | **0.036** | **0.0027 good** |
| impressed | 500 | 0.040 | 0.0052 good |
| impressed | 300 | 0.042 | 0.0085 good |
| solenoidal | 300 | 0.082 | 0.0063 good |

**Conclusion**: GMRES **stagnates** at ~4% (impressed) / ~8% (solenoidal); not a tol/outer issue. High λ can push divA to good, but strict tol=1e-4 not closed.

---

## 8. Timeline

| Date | Milestone |
|------|--------|
| 2026-05-28 | Sign fix; sign/curl/impressed-divJ; fused vs debug |
| 2026-05-28 | solenoidal source; scalar PC sweep → bottleneck in PC |
| 2026-05-28 | `JacobiGradDivABlock` 3×3 PC; true_rel greatly improved |
| 2026-05-29 | GMRES outer/restart/tol sensitivity → rel floor ~4%, **Inner prototype closed** |

---

## 10. Stage 10.8 — divA discretization unification (2026-05-21)

### 10.1 Parameter clarification

| Parameter | Equation | Role | Relation to other parameter |
|------|------|------|----------------|
| `λ_φ` (`phi_gauge_penalty`) | `LhsPhi += λ_φ·phi` | Lock φ constant gauge | **Independent** |
| `λ_A` (`a_divergence_penalty`) | `LhsA -= λ_A·grad(divA)` | Penalize divA | **Independent** |

**PC (preconditioner)**: block Jacobi 8×8, local block-diagonal inverse approximation for A/φ; A-penalty needs `JacobiGradDivABlock` 3×3 block.

### 10.2 Root cause: apply and PC discretization inconsistent (Stage 10.7)

| Component | Stage 10.7 | Stage 10.8 default |
|------|------------|-----------------|
| Fused apply | B-corrected pipeline | **PairwiseUncorrected** |
| Debug assemble | B-corrected | **PairwiseUncorrected** (consistent with fused) |
| PC | Pairwise 3×3 | unchanged |

`AphiADivergencePenaltyMode`: `BCorrectedTrace` (diagnostic) / `PairwiseUncorrected` (operator default).

### 10.3 New code / tests

| File | Role |
|------|------|
| `aphi_pairwise_div_a_ck.*` | pairwise divA + grad(divA) |
| `aphi_pairwise_a_divergence_penalty_pipeline.h` | penalty apply same family as PC |
| `diagnostics/aphi_div_a_discretization_comparison_helpers.h` | two-route L2 / energy sign comparison |
| `diagnostics/aphi_gmres_residual_decomposition_helpers.h` | res_A / res_phi block decomposition |
| `test_3d_aphi_ck_div_a_discretization_comparison_diagnostic` | Route B vs Route P gate |

GMRES / sweep output adds: `res_A_frac`, `res_phi_frac`.

### 10.4 Research defaults + PC gate (Stage 10.8 continued)

- Code constants: `AphiADivergencePenaltyResearchDefaults` (`eta_a_min=0.1`, `eta_a_max=0.3`)
- PC gate: `graddiv_block_diagonal` pairwise FD vs pairwise PC, `core_max_diff≈2.9e-6`
- Inner gate: `test_3d_aphi_ck_inner_a_divergence_penalty_gate_diagnostic` **passed=1** (solenoidal η=0.1/0.3)

See [`CURSOR_APHI_STAGE10_8_PAIRWISE_DIVA_RECORD.md`](CURSOR_APHI_STAGE10_8_PAIRWISE_DIVA_RECORD.md) §13.
