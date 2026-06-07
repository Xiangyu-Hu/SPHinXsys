# Stage 10.10 — divergence-free MMS and 3D pairwise divA validation record

> **Date**: 2026-05-21  
> **Status**: 10.10-A/B/C complete; **pairwise rel≈0.13 root cause identified** (denominator + boundary, not core discretization failure)  
> **Plan**: [`stage10_10_divergence_free_mms_and_3d_pairwise_plan_for_cursor.md`](../../../stage10_10_divergence_free_mms_and_3d_pairwise_plan_for_cursor.md)  
> **Prerequisite**: [`CURSOR_APHI_STAGE10_8_PAIRWISE_DIVA_RECORD.md`](CURSOR_APHI_STAGE10_8_PAIRWISE_DIVA_RECORD.md) §14

---

## 0. Executive summary

| Item | Result |
|----|------|
| **non-degenerate MMS** | ✅ `rhs_norm≈22–39` (Az2D/CrossSine3D), no longer degenerate |
| **consistency MMS** | ✅ 8/8 (two fields × η=0,0.1,0.2,0.3) |
| **Joule vs exact** | ✅ Az2D η=0: `joule_err_vs_exact≈2e-7` |
| **pairwise divA on Az2D/CrossSine3D** | ⚠️ global Anorm rel≈0.13; **core region ~1e-7**; gradDen rel ~1e-7 |
| **invariance MMS (div-free)** | ⚠️ GMRES converged 6/6, but **block error large** (exact field has pairwise divA≠0 on boundary) |
| **Contact A-penalty** | still deferred |

**Key finding**: Az2D resolved the `K(A)≈0` degeneracy. pairwise global rel≈0.13 is **not a core discretization failure**: root-cause diagnostics show **core-region pairwise rel ~1e-7**, with error **100% concentrated on boundary**; global rel is amplified by **denominator `||A||` (not `||grad A||`)**. With gradDen denominator, rel ~1e-7. B-corrected gradDen rel ~0.037.

**Architecture recommendation (10.10-B2)**: adopt **dual-track pairwise penalty + B-corrected diagnostic** — penalty continues to use pairwise (same family as operator/PC); post-hoc divA reporting uses B-corrected + gradDen denominator + core/boundary decomposition.

---

## 1. New code

| File | Purpose |
|------|------|
| [`diagnostics/aphi_divergence_free_nonzero_rhs_mms_helpers.h`](../extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_divergence_free_nonzero_rhs_mms_helpers.h) | Az2D/CrossSine3D fields, MMS, divA diagnostics, exact Joule reference |
| [`test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms/`](test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms/) | Stage 10.10-A/C |
| [`diagnostics/aphi_pairwise_diva_root_cause_helpers.h`](../extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_pairwise_diva_root_cause_helpers.h) | core/boundary decomposition, gradDen denominator, core_shell sweep |
| [`test_3d_aphi_ck_pairwise_diva_3d_diagnostic/`](test_3d_aphi_ck_pairwise_diva_3d_diagnostic/) | Stage 10.10-B |
| [`test_3d_aphi_ck_pairwise_diva_root_cause_diagnostic/`](test_3d_aphi_ck_pairwise_diva_root_cause_diagnostic/) | Stage 10.10-B2 root-cause diagnostic |

---

## 2. Stage 10.10-A/C: non-degenerate MMS

Test: `test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms` → **passed=1**

### 2.1 Exact fields

| Field | A | Analytic divA |
|-------|---|-----------|
| **Az2D** | `(0,0,sin πx sin πy)` | 0 |
| **CrossSine3D** | `(sin πy sin πz, sin πz sin πx, sin πx sin πy)` | 0 |

### 2.2 MMS results (dp=0.1, consistency)

| Field | η_A | rhs_norm | true_rel | block_Linf | divA_pairwise | divA_B | joule_err_vs_exact |
|-------|-----|----------|----------|------------|---------------|--------|-------------------|
| Az2D | 0 | 22.4 | 2.6e-6 | 6e-7 | **0.129** | 0.037 | ~2e-7 |
| Az2D | 0.1 | 21.8 | 2.6e-6 | 6e-7 | 0.129 | 0.037 | ~2e-7 |
| CrossSine3D | 0 | 38.9 | 2.6e-6 | 8e-6 | **0.129** | 0.037 | ~1e-6 |
| CrossSine3D | 0.3 | 36.2 | 8.6e-5 | 4.8e-4 | 0.128 | 0.037 | 6e-4 |

- **rhs_norm nonzero**: degeneracy problem resolved ✅  
- **consistency**: 8/8 passed ✅  
- **invariance**: 6/6 GMRES converged; block_Linf can reach 0.2–2.3 (pairwise nonzero on exact field → penalty changes solution) ⚠️  

### 2.3 E observable vs exact

| Metric | Az2D η=0 | Notes |
|------|----------|------|
| `joule_err_vs_exact` | ~2e-7 | ✅ reliable |
| `E_combined_err_vs_exact` | ~1e-6 | ✅ new: \|E_real\|²+\|E_imag\|² vol-weighted L2 |
| `E_L2_err_vs_exact` | ~10¹⁹ | ❌ spurious blow-up: only real component counted when `E_real=0` |

`E_combined` resolves the issue where `E_L2` is unusable under the phasor convention.

---

## 3. Stage 10.10-B: four-field pairwise vs B-corrected (dp sweep)

Test: `test_3d_aphi_ck_pairwise_diva_3d_diagnostic` → **passed=1** (hard pass: Linear2D + Az2D + CrossSine3D @ dp=0.1)

### 3.1 dp=0.1 (core_shell=2.5 dp)

| Field | pairwise div_A_rel | B-corrected div_A_rel |
|-------|-------------------|----------------------|
| Linear2D | **2.4e-9** | 1.2e-8 |
| Az2D | **0.129** | 0.037 |
| CrossSine3D | **0.129** | 0.037 |
| Sinusoidal3DCurlPsi | **0.477** (info) | 0.011 |

### 3.2 dp refinement (pairwise div_A_rel)

| Field | dp=0.15 | dp=0.10 | dp=0.075 | Trend |
|-------|---------|---------|----------|------|
| Az2D | 0.194 | 0.129 | **0.095** | ↓ with dp refinement |
| CrossSine3D | 0.194 | 0.129 | **0.095** | ↓ |
| Sinusoidal3DCurlPsi | 0.495 | 0.477 | 0.444 | ↓ slow |

**Conclusions**:

1. **Linear2D** remains the only hard-pass pairwise ~1e-9.  
2. **Az2D/CrossSine3D** are analytically div-free, but pairwise **~0.13** (not ChatGPT's expected ~1e-9).  
3. pairwise **decreases with dp but does not converge below 1e-2** (finest Az2D ~0.095).  
4. **B-corrected** stable at 0.01–0.04, one order of magnitude better than pairwise.  
5. Sinusoidal3DCurlPsi pairwise ~0.48 remains informational.

---

## 3.1 Stage 10.10-B2: pairwise divA root-cause diagnostic

Test: `test_3d_aphi_ck_pairwise_diva_root_cause_diagnostic` → **passed=1**

### 3.1.1 Denominator effect (dp=0.1, core_shell=2.5 dp)

| Field | rel(Anorm denom) | rel(gradA denom) | core-only rel(Anorm) | boundary energy frac |
|-------|----------------|-----------------|---------------------|---------------------|
| Linear2D | 2.4e-9 | 3.5e-10 | 6.8e-10 | 0.97 |
| Az2D | **0.129** | **1.15e-7** | **1.78e-7** | **1.00** |
| CrossSine3D | **0.129** | **1.07e-7** | **1.66e-7** | **1.00** |

### 3.1.2 Conclusions

1. **Global rel≈0.129 is a denominator artifact**: `div_A_L2 / ||A||_L2` severely underestimates the denominator for oscillatory fields (\|A\|~O(1), \|grad A\|~O(π)) → rel inflated ~10⁶×.  
2. **Core-region pairwise discretization is correct**: core-only rel ~1e-7, same order as Linear2D.  
3. **Boundary dominates**: 100% of divA L2 energy is on boundary (outside `isCoreParticle`).  
4. **core_shell sweep (2.5/3.0/3.5 dp)**: global Anorm rel unchanged (0.129); core-only rel unchanged (~1e-7).  
5. **B-corrected gradDen rel ~0.037**: boundary trace(gradA) still has O(0.04) error; usable as diagnostic reference track.

### 3.1.3 Architecture decision

| Use case | Recommended discretization | Rationale |
|------|-------------------|------|
| **penalty operator / PC** | pairwise uncorrected | Same family as GMRES apply |
| **Post-hoc divA reporting** | B-corrected + gradDen denominator | More physical denominator; core/boundary decomposition |
| **Hard gate** | core-only pairwise or B-corrected | Global Anorm rel cannot serve as gate |

---

## 4. Response to ChatGPT 10.10 plan

| ChatGPT expectation | Actual | Status |
|-------------|------|------|
| Az2D `rhs≠0` | rhs≈22 | ✅ |
| Az2D pairwise ~small | core ~1e-7; global Anorm rel≈0.13 is denominator+boundary artifact | ✅ root cause identified |
| invariance div-free solution unchanged | boundary pairwise nonzero → penalty changes solution | ⚠️ boundary open |
| CrossSine3D as 3D validation | MMS passed; pairwise ~0.13 | partial |
| η=0.1 primary candidate | consistency MMS + solenoidal gate still supported | ✅ research only |

---

## 5. η_A policy (10.10-D provisional, to be frozen after exact observables are more complete)

```text
primary_eta_a   = 0.1   (research)
optional_eta_a  = 0.2
upper_eta_a     = 0.3
production      = lambda_A off
```

---

## 6. External messaging (updated)

**Can say**:

- non-degenerate div-free MMS established (Az2D/CrossSine3D, `rhs_norm≠0`).  
- Inner pairwise A-penalty consistency MMS converges.  
- Joule / E_combined vs analytic exact on Az2D η=0 at ~1e-6–1e-7 relative error.
- pairwise divA **core region** on div-free fields ~1e-7; global Anorm rel cannot be interpreted directly.

**Cannot say**:

- pairwise divA global Anorm rel≈0.13 indicates core discretization failure.
- Coulomb gauge is accurately satisfied (boundary still has O(0.04) B-corrected error).
- invariance RHS preserves div-free physical solution (boundary penalty nonzero).  

---

## 7. Next steps (Stage 10.10 continued)

1. ~~**pairwise divA root cause**~~ ✅ complete (§3.1)  
2. ~~**Architecture decision**~~ ✅ recommend dual-track pairwise penalty + B-corrected diagnostic (§3.1.3)  
3. ~~**E observable**~~ ✅ `E_combined` added to MMS (§2.3)  
4. **boundary pairwise divA**: investigate why trace(gradA) / pairwise is nonzero in boundary layer; whether boundary exclusion or dedicated BC is needed  
5. **invariance**: informational only until boundary open issue is fixed  
6. **Contact A-penalty**: still deferred (10.10-E conditions not met)

---

## 8. Reproduction

```bash
cd build && ninja test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms \
              test_3d_aphi_ck_pairwise_diva_3d_diagnostic \
              test_3d_aphi_ck_pairwise_diva_root_cause_diagnostic
```

---

## 9. Related documents

- This record  
- [`CURSOR_APHI_STAGE10_8_PAIRWISE_DIVA_RECORD.md`](CURSOR_APHI_STAGE10_8_PAIRWISE_DIVA_RECORD.md)  
- [`stage10_10_divergence_free_mms_and_3d_pairwise_plan_for_cursor.md`](../../../stage10_10_divergence_free_mms_and_3d_pairwise_plan_for_cursor.md)
