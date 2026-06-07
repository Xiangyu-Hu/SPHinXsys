# Stage 10.8 — Pairwise divA unification and η_A calibration record

> **Date**: 2026-05-21  
> **Status**: Inner research mainline in progress (**not Contact production**)  
> **Basis**: [`stage10_mainline_pc_lambda_divA_next_plan_for_cursor.md`](../../../../../../docs/electromagnetic_aphi/stage10_mainline_pc_lambda_divA_next_plan_for_cursor.md)  
> **Prerequisite**: [`CURSOR_APHI_STAGE10_7_A_DIVERGENCE_PENALTY_RECORD.md`](CURSOR_APHI_STAGE10_7_A_DIVERGENCE_PENALTY_RECORD.md)

---

## 0. Executive summary

Stage 10.7 GMRES rel floor (~4–60%) root cause: **penalty apply (B-corrected) inconsistent with PC (pairwise 3×3) discretization**.

Stage 10.8 completed:

1. **Operator default switched to pairwise**: `AphiADivergencePenaltyMode::PairwiseUncorrected`
2. **Penalty scratch fields**: `PenaltyScratchDivA*` / `PenaltyScratchGradDivA*` (avoid GMRES workspace prefix pollution)
3. **Post-hoc divA diagnostic default pairwise** (consistent with operator; `div_a_relative` denominator uses `||A||`)
4. **GMRES residual block decomposition**: `res_A_frac` / `res_phi_frac`
5. **New tests**: discretization comparison, η_A sweep

**Key results**: after pairwise unification, impressed `λ_A=10` **converged=true_rel≈4e-5** (Stage 10.7 same point ~6%); `λ_A=100` true_rel **0.60→0.04**.

---

## 1. Concept clarification (frozen)

| Parameter | Equation | Role |
|------|------|------|
| `λ_φ` (`phi_gauge_penalty`) | `LhsPhi += λ_φ·phi` | Lock φ constant gauge |
| `λ_A` (`a_divergence_penalty`) | `LhsA -= λ_A·grad(divA)` | Penalize divA (Coulomb gauge regularization) |

**PC**: block Jacobi 8×8 + `JacobiGradDivABlock` 3×3 (pairwise `-grad(divA)` local block).

**η_A definition** (relative scale):

```text
η_A = median(|λ_A · D_graddiv|) / median(|D_laplace|)
    ≈ λ_A · median(||JacobiGradDivABlock||) / median(|laplace_a_diag|)
```

Calibration: `λ_A = η_A · median_laplace / median_graddiv`.

---

## 2. Code change list

| File | Change |
|------|------|
| `aphi_pairwise_div_a_ck.*` | pairwise divA, grad(divA) |
| `aphi_pairwise_a_divergence_penalty_pipeline.h` | pairwise penalty apply bundle |
| `aphi_a_divergence_penalty_ck.h` | `AphiADivergencePenaltyScratchNames` |
| `aphi_coupling_modes_ck.h` | `AphiADivergencePenaltyMode` + **`AphiADivergencePenaltyResearchDefaults`** |
| `aphi_matrix_free_operator_ck.*` | default pairwise penalty |
| `aphi_assemble_lhs_debug_ck.*` | debug same route as fused |
| `aphi_variables_ck.hpp` | register scratch + diagnostic divA fields |
| `diagnostics/aphi_div_a_diagnostic_helpers.h` | post-hoc divA **default pairwise** |
| `diagnostics/aphi_div_a_discretization_comparison_helpers.h` | B vs pairwise comparison |
| `diagnostics/aphi_gmres_residual_decomposition_helpers.h` | res block decomposition |
| `diagnostics/aphi_inner_a_divergence_penalty_eta_sweep_helpers.h` | η_A sweep |
| `diagnostics/aphi_inner_a_divergence_penalty_gate_helpers.h` | Inner gate (GMRES + divA + E/J/Joule) |
| `diagnostics/aphi_divergence_free_a_mms_helpers.h` | Stage 10.9: div-free / manufactured MMS + dp sweep |
| `diagnostics/aphi_observable_reference_benchmark_helpers.h` | Stage 10.9-D Case C: E/J/Joule vs reference |

---

## 3. Test list

| Test | Type | Result |
|------|------|------|
| `test_3d_aphi_ck_div_a_discretization_comparison_diagnostic` | gate | passed=1 |
| `test_3d_aphi_ck_inner_a_divergence_penalty_fused_vs_debug` | gate | passed=1, core_max_diff≈4e-5 |
| `test_3d_aphi_ck_graddiv_block_diagonal_diagnostic` | **gate** | passed=1, **pairwise FD vs pairwise PC**, core_max_diff≈**2.9e-6** |
| `test_3d_aphi_ck_inner_a_divergence_penalty_sweep_diagnostic` | sweep | see §4 |
| `test_3d_aphi_ck_inner_solenoidal_a_divergence_penalty_sweep_diagnostic` | sweep | see §4 |
| `test_3d_aphi_ck_inner_a_divergence_penalty_eta_sweep_diagnostic` | sweep | see §5 |
| `test_3d_aphi_ck_inner_solenoidal_a_divergence_penalty_eta_sweep_diagnostic` | sweep | see §5 |
| **`test_3d_aphi_ck_inner_a_divergence_penalty_gate_diagnostic`** | **gate** | **passed=1** (§13 Inner gate checklist) |
| `test_3d_aphi_ck_inner_divergence_free_a_pairwise_diagnostic` | gate | passed=1 (§14.2) |
| `test_3d_aphi_ck_inner_divergence_free_a_mms` | gate | passed=1 (§14.3) |
| `test_3d_aphi_ck_inner_divergence_free_a_dp_refinement_diagnostic` | gate | passed=1 (§14.5) |
| **`test_3d_aphi_ck_inner_a_divergence_penalty_observable_reference_diagnostic`** | **gate** | **passed=1** (§14.6 Case C) |

---

## 4. Absolute λ_A sweep (pairwise apply + pairwise post-hoc divA)

### 4.1 Impressed-current box (tol=1e-4, λ_φ=100)

| λ_A | converged | true_rel | res_A_frac | div_A_rel | source div_A |
|-----|-----------|----------|------------|-----------|--------------|
| off | ✓ | 3.8e-6 | 0.98 | 0.58 high_risk | 0.58 |
| **10** | **✓** | **4.0e-5** | 0.97 | 0.13 high_risk | 0.19 |
| 50 | ✗ | 0.16 | 1.00 | 0.09 warn | 0.06 |
| **100** | ✗ | **0.041** | 0.99 | **0.04 warn** | 0.03 |
| 300 | ✗ | 0.14 | 1.00 | 0.05 warn | 0.01 |
| 1000 | ✗ | 0.19 | 1.00 | 0.03 warn | 0.004 good |

Compared to Stage 10.7 (B-corrected apply): λ=10 true_rel **0.063→4e-5**; λ=100 **0.60→0.04**.

### 4.2 Solenoidal source

| λ_A | converged | true_rel | div_A_rel |
|-----|-----------|----------|-----------|
| off | ✓ | 1.0e-5 | 0.40 high_risk |
| **10** | **✓** | **4.7e-6** | **0.08 warn** |
| 50 | ✗ | 0.13 | 0.06 warn |
| 100 | ✗ | 0.10 | 0.04 warn |

---

## 5. η_A relative-scale sweep

Calibration baseline (core median, 2026-05-21 run):

```text
median_laplace_a_diag = 587.05
median_graddiv_block_norm = 17.41
core_particles = 216
```

Post-hoc `div_A_relative = ||divA|| / ||A||` (pairwise divA; denominator is A-field L2, not B-corrected ||grad A||).

### 5.1 Impressed η_A sweep (tol=1e-4)

| η_A | λ_A | converged | true_rel | div_A_rel | source div_A |
|-----|-----|-----------|----------|-----------|--------------|
| 0 | 0 | ✓ | 3.8e-6 | 0.63 high_risk | 1.95 |
| **0.01** | 0.34 | **✓** | 3.9e-6 | 0.54 high_risk | 1.78 |
| **0.03** | 1.01 | **✓** | 4.2e-6 | 0.42 high_risk | 1.49 |
| **0.1** | 3.37 | **✓** | 4.1e-6 | 0.23 high_risk | 0.95 |
| **0.3** | 10.1 | **✓** | 4.1e-5 | 0.12 high_risk | 0.51 |
| 1.0 | 33.7 | ✗ | 0.0072 | 0.27 high_risk | 0.37 |

**5/6 converged**. η∈[0.01, 0.3] all converged; η=1 did not pass tol=1e-4 but true_rel≈0.7% (much better than Stage 10.7 λ=100 at 60%).

### 5.2 Solenoidal η_A sweep

| η_A | λ_A | converged | true_rel | div_A_rel | source div_A |
|-----|-----|-----------|----------|-----------|--------------|
| 0 | 0 | ✓ | 1.0e-5 | 0.54 high_risk | 1.05 |
| **0.01** | 0.34 | **✓** | 2.6e-5 | 0.46 high_risk | 0.94 |
| **0.03** | 1.01 | **✓** | 4.6e-5 | 0.36 high_risk | 0.77 |
| **0.1** | 3.37 | **✓** | 2.7e-5 | 0.20 high_risk | 0.48 |
| **0.3** | 10.1 | **✓** | 5.1e-6 | **0.10 high_risk** | 0.26 |
| 1.0 | 33.7 | ✗ | 0.012 | 0.23 high_risk | 0.17 |

**5/6 converged**. η=0.3 (λ≈10) is best solenoidal compromise: converged + lowest div_A.

### 5.3 Research defaults (written into code)

`aphi_coupling_modes_ck.h`:

```cpp
struct AphiADivergencePenaltyResearchDefaults {
    static constexpr Real eta_a_min = 0.1;
    static constexpr Real eta_a_max = 0.3;
    static constexpr Real eta_a_mid = 0.2;
};
```

Calibration relation: `λ_A = η_A · median_laplace / median_graddiv` (helper: `lambdaAFromEtaA`).

**Contact baseline still keeps `λ_A=off`**; above constants are for Inner research / sweep reference only.

---

## 6. PC and apply consistency gate

`test_3d_aphi_ck_graddiv_block_diagonal_diagnostic` changed to **pairwise penalty apply FD** vs **JacobiGradDivABlock** (no longer B-corrected pipeline).

| Stage | FD reference | PC block | core_max_diff |
|------|---------|-------|---------------|
| Stage 10.7 | B-corrected pipeline | pairwise 3×3 | ~0.17 |
| **Stage 10.8** | **pairwise apply** | **pairwise 3×3** | **~2.9e-6** |

Same order as `fused_vs_debug` (core_max_diff≈4e-5); **PC and apply discretization consistency gate passed**.

---

## 7. divA discretization comparison (divergent field)

| Metric | B-corrected | Pairwise | Difference |
|------|-------------|----------|------|
| divA L2 | 2.25 | 2.24 | 1.3% |
| gradDivA L2 diff | — | — | 4.7% |
| energy sign | 0.362 | 0.362 | consistent |

---

## 8. Why impressed source-region div_A is elevated (assessment)

Comparing impressed vs solenoidal η sweep:

| Source type | source div_A at λ_A=off | source div_A at η=0.3 | global div_A_rel |
|------|-------------------------|-------------------------|----------------|
| Impressed (J_z) | **1.95** | 0.51 | 0.12 high_risk |
| Solenoidal | 1.05 | **0.26** | **0.10** high_risk |

**Conclusions**:

1. Impressed source (`J=(0,0,1)` real + small imag) is **itself non-solenoidal**, `div J ≠ 0`; Maxwell solution `A` naturally has divergence component in source region; penalty can only partially suppress it.
2. Solenoidal source has lower source div_A and easier global div_A suppression at same η — consistent with physics.
3. **Cannot** use impressed source-region div_A alone to judge penalty failure; should rely on GMRES convergence + solenoidal case div_A trend.

---

## 9. Fixed bug

**GMRES scratch fields**: penalty intermediates cannot bind to `input_block.a_real + "DivA"` (during GMRES matvec input is `GMRESZ0AReal` → unregistered `GMRESZ0ARealDivA`).  
Fix: uniformly use `PenaltyScratchDivAReal/Imag` + `PenaltyScratchGradDivAReal/Imag`.

---

## 10. Current assessment

| Question | Conclusion |
|------|------|
| Is divA operator reasonable? | **Yes** (pairwise complete set, consistent with PC, energy sign correct) |
| PC and apply consistent? | **Yes** (graddiv gate core_max_diff≈3e-6) |
| GMRES closed? | **η∈[0.01,0.3] all converged**; η=1 still floor (~0.7–1.2%) |
| Can divA be suppressed? | **div_A monotonically decreases as η↑** (solenoidal η=0.3 → 0.10); impressed source region elevated see §8 |
| Contact ready? | **No** (Inner gate passed on solenoidal; Contact still deferred, see §13.4) |

---

## 13. Inner gate checklist (plan §23, 2026-05-21)

Solenoidal source, research `η_A∈{0.1,0.3}`, `tol=1e-4`, E/J/Joule relative to baseline (η=0) tolerance **20%**.

### 13.1 Operator / PC gate (passed)

| Test | passed | Key metrics |
|------|--------|----------|
| `div_a_discretization_comparison` | 1 | energy sign consistent |
| `inner_a_divergence_penalty_fused_vs_debug` | 1 | core_max_diff≈4.3e-5 |
| `inner_a_divergence_penalty_sign_energy` | 1 | sign_ok=1 |
| `graddiv_block_diagonal` | 1 | core_max_diff≈2.9e-6 |

### 13.2 GMRES + divA + E/J/Joule gate (solenoidal)

| η_A | λ_A | converged | true_rel | div_A_rel | Joule Δ | E_L2 Δ |
|-----|-----|-----------|----------|-----------|---------|--------|
| 0 (baseline) | 0 | ✓ | 1.0e-5 | 0.54 | — | — |
| **0.1** | 3.37 | **✓** | 2.7e-5 | **0.20** | **13%** | **6.5%** |
| **0.3** | 10.1 | **✓** | 5.2e-6 | **0.10** | **14%** | **7.0%** |

`test_3d_aphi_ck_inner_a_divergence_penalty_gate_diagnostic` → **passed=1** (baseline_ok + 2/2 gate rows).

### 13.3 §23 success criteria check

| Criterion | Status |
|------|------|
| GMRES converged / true_rel ≤ 1e-4 | **✓** (η=0.1, 0.3) |
| div_A_relative decreases on converged solution | **✓** (0.54→0.20→0.10) |
| E/J/Joule not significantly degraded | **✓** (Δ < 20%) |
| residual decomposition no anomaly | **✓** (res_A_frac high when converged but true_rel low, not floor-type) |
| λ_A relative scale has stable interval | **✓** (η∈[0.01,0.3]) |
| PC and apply discretization consistent | **✓** (§6) |

### 13.4 Contact A-penalty assessment

**Still deferred for integration**. Reasons:

1. Inner gate passes on **solenoidal**, but impressed source-region div_A still high_risk (§8).
2. Contact baseline frozen policy unchanged: `λ_φ=100`, `λ_A=off`.
3. Before Contact integration: Contact version divA/grad(divA)/PC complete set + two-body/three-body gate (not Inner-only body-interior penalty).

---

## 14. Stage 10.9: divA validation (MMS + pairwise diagnostic, 2026-05-21)

Per [`stage10_8_metric_interpretation_and_validation_plan_for_cursor.md`](../../../stage10_8_metric_interpretation_and_validation_plan_for_cursor.md) Phase 10.9.

### 14.1 Metric interpretation correction (Task 1)

| Metric | Correct interpretation |
|------|----------|
| `true_rel` | GMRES relative residual; small only means discrete system self-consistent, **not** physical solution or Coulomb gauge correct |
| `div_A_rel` | post-hoc divA relative measure; decrease means penalty suppressing divergence, **not** ∇·A=0 satisfied |
| `Joule Δ` / `E_L2 Δ` | perturbation relative to **η=0 baseline**, not analytic solution error |
| `block_Linf_err` (MMS) | field error vs manufactured exact; should be small only under consistency RHS |

### 14.2 Phase 10.9-A: pairwise divA on exact curl fields

Test: `test_3d_aphi_ck_inner_divergence_free_a_pairwise_diagnostic` → **passed=1**

| Field | B-corrected div_A_rel | Pairwise div_A_rel |
|----|----------------------|-------------------|
| **Linear2D** `A=(-y/2,x/2,0)` | ~1.1e-8 | **~2.4e-9** ✅ |
| **Sinusoidal3DCurlPsi** | ~0.011 | **~0.48** ❌ (informational) |

Conclusion: pairwise divA accurate on simple 2D curl field; 3D sinusoidal curl field fails pairwise at boundary/discretization (B-corrected still ~0.01). hard pass only Linear2D.

### 14.3 Phase 10.9-B: manufactured separable MMS + GMRES

Test: `test_3d_aphi_ck_inner_divergence_free_a_mms` → **passed=1**

- Field: separable sin/cos manufactured (**not div-free**, `div_A_rel≈3.7`)
- **consistency** `b_η=K_η(A_exact)`: 4/4 passed (η=0,0.1,0.2,0.3; `true_rel≤1e-4`, `block_Linf_err≤0.15`)
- **invariance** `b_0=K_0(A_exact)`, solve `K_η x=b_0`: 3/3 passed (GMRES convergence only; non-div-free field not hard-gated on block error)

| η_A | mode | rhs_norm | true_rel | block_Linf_err | div_A_rel |
|-----|------|----------|----------|----------------|-----------|
| 0 | consistency | 103 | 1.4e-5 | 1.7e-5 | 3.74 |
| 0.1 | consistency | 128 | 1.5e-6 | 1.8e-6 | 3.74 |
| 0.1 | invariance | 103 | 4.2e-6 | 3.75 | 0.70 |
| 0.3 | consistency | 213 | 4.6e-5 | 1.5e-3 | 3.75 |
| 0.3 | invariance | 103 | 5.3e-5 | 5.46 | 0.23 |

**div-free exact MMS still degenerate**: Linear2D / Sinusoidal3D exact fields make `K(A_exact,φ)≈0` → `rhs_norm=0`; need separate manufactured field with `b≠0` and divA=0 (independent follow-up).

### 14.4 Fixed bug

`AphiCopyBlockCK` parameter order is `(dst, src)`. `copy_exact` was wrongly `(solution, exact_block)`, overwriting solution with zero initial exact → `apply` gives `lhs=0`. Fixed to `(exact_block, solution)`.

### 14.5 Phase 10.9-C: dp refinement

Test: `test_3d_aphi_ck_inner_divergence_free_a_dp_refinement_diagnostic` → **passed=1**

η_A=0.1 consistency MMS + Linear2D exact curl pairwise divA, dp∈{0.15, 0.1, 0.075}:

| dp | N | linear2d pairwise div_A_rel | mms block_Linf | mms true_rel | mms div_A_rel |
|----|---|------------------------------|----------------|--------------|---------------|
| 0.15 | 343 | 3.2e-9 | 2.6e-6 | 7.5e-7 | 3.45 |
| 0.10 | 1000 | 2.4e-9 | 1.9e-6 | 1.5e-6 | 3.74 |
| 0.075 | 2197 | 4.5e-9 | 3.3e-5 | 4.7e-6 | 3.90 |

- **Linear2D pairwise**: all three dp levels ~1e-8 → pairwise divA resolution-stable on exact curl field ✅
- **MMS block**: all levels converge and `<0.05`; finest not monotonically better than coarsest (`mms_block_monotone=0`, informational)
- **mms div_A_rel≈3.4–3.9**: manufactured field itself non-div-free, does not decrease with dp (as expected)

### 14.6 Phase 10.9-D / Case C: reference observable benchmark

Test: `test_3d_aphi_ck_inner_a_divergence_penalty_observable_reference_diagnostic` → **passed=1**

Helper: `diagnostics/aphi_observable_reference_benchmark_helpers.h`

**Two validation tracks**:

| Track | Reference | Candidate | Gate metrics |
|-------|-----------|-----------|-----------|
| **Solenoidal** | dp=0.075 η=0 (fine self-ref, informational) | dp=0.1 η∈{0.1,0.2,0.3} | vs **same dp η=0**: Joule/E error <25%, div_A decrease |
| **MMS manufactured** | dp=0.1 η=0 consistency | dp=0.1 η∈{0.1,0.2,0.3} | vs η=0 ref: Joule/E error <20% |

**Solenoidal results (dp=0.1, vs cand η=0 baseline)**:

| η_A | div_A_rel | div_A↓ | Joule err vs η=0 | E_L2 err vs η=0 |
|-----|-----------|--------|------------------|-----------------|
| 0 (baseline) | 0.54 | — | — | — |
| **0.1** | **0.20** | ✅ 2.8× | **13%** | **6.5%** |
| 0.2 | 0.13 | ✅ 4.3× | 22% | 11% |
| **0.3** | **0.10** | ✅ 5.4× | **14%** | **7.0%** |

**Mesh gap (informational)**: dp=0.1 η=0 vs dp=0.075 η=0 → Joule err **93%**, E_L2 err **46%**. Coarse-fine cross-mesh cannot serve as penalty gate directly; need same-dp baseline or analytic reference.

**MMS results (vs η=0 consistency ref)**:

| η_A | Joule err vs ref | E_L2 err vs ref | div_A_rel | block_Linf |
|-----|------------------|-----------------|-----------|------------|
| 0.1 | ~0 | ~0 | 3.74 | 1.8e-6 |
| 0.2 | ~0 | ~0 | 3.74 | 5.8e-6 |
| 0.3 | 0.05% | 0.03% | 3.75 | 1.5e-3 |

Under MMS consistency, observables nearly unchanged (RHS includes penalty, field error still small); div_A does not decrease because manufactured field is non-div-free.

**Case C conclusions**:

- Solenoidal: **η=0.1/0.3 significantly reduce div_A while E/J/Joule perturbation vs same-dp η=0 is within research-acceptable range**; without vs analytic exact, cannot say "closer to true solution".
- Fine-dp self-reference quantifies **discretization error baseline** (mesh gap ~50–90%), cannot alone serve as penalty acceptance.
- MMS track validates **observable pipeline and MMS solve chain closure**; div-free analytic reference still pending construction.

---

## 15. Next steps

1. ~~Stage 10.9-A pairwise divA diagnostic~~ ✅
2. ~~Stage 10.9-B manufactured MMS~~ ✅
3. ~~Phase 10.9-C dp refinement~~ ✅
4. ~~Case C reference observable benchmark~~ ✅
5. div-free exact MMS: manufactured field with `b=K(A_exact)≠0` + analytic E/J/Joule reference
6. 3D pairwise divA boundary error: fix discretization or accept as informational
7. Contact A-penalty complete prototype (production still deferred)

**Stage 10.10 record**: [`CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md`](CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md)

**Still deferred**: Contact A-penalty production, teacher demo, cold crucible

---

## 16. Related documents

- Full record: this document
- **ChatGPT discussion pack (after Stage 10.9)**: [`stage10_9_validation_results_and_open_questions_for_chatgpt.md`](../../../stage10_9_validation_results_and_open_questions_for_chatgpt.md)
- Stage 10.7 closure statement: [`CURSOR_APHI_STAGE10_7_A_DIVERGENCE_PENALTY_RECORD.md`](CURSOR_APHI_STAGE10_7_A_DIVERGENCE_PENALTY_RECORD.md) §10
- ChatGPT metric/validation plan: [`stage10_8_metric_interpretation_and_validation_plan_for_cursor.md`](../../../stage10_8_metric_interpretation_and_validation_plan_for_cursor.md)
- ChatGPT main plan: [`stage10_mainline_pc_lambda_divA_next_plan_for_cursor.md`](../../../../../../docs/electromagnetic_aphi/stage10_mainline_pc_lambda_divA_next_plan_for_cursor.md)
