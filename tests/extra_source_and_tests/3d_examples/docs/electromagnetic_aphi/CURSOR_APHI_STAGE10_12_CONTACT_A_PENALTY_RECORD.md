# Stage 10.12-B — Contact A-penalty two-body research gate record

> **Date**: 2026-05-21  
> **Status**: C1–C5 + **C3+** + **Stage 10.13 P0 gate hardening** **passed=1**; Contact two-body research gate closed  
> **Plan**: [`stage10_12_contact_a_penalty_and_boundary_b_validation_plan_for_cursor.md`](../../../stage10_12_contact_a_penalty_and_boundary_b_validation_plan_for_cursor.md)  
> **Prerequisite**: [`CURSOR_APHI_STAGE10_11_BOUNDARY_VALIDATION_RECORD.md`](CURSOR_APHI_STAGE10_11_BOUNDARY_VALIDATION_RECORD.md) (includes 10.12-A B dp refinement)

---

## 0. One-line conclusion

Contact two-body **A-penalty research gate closed** for **η=0 / η=0.1** with **InnerOnly** penalty stencil (main K still Inner+Contact); **η=0.2** usable with **`polish_sweeps=0`** (`global_true_rel` ~8e-6). Original "η>0 GMRES non-convergence" root cause: **InnerContact penalty cross split-body contamination** — fixed in C4. Three-body A-penalty **deferred**.

---

## 1. Architecture freeze (aligned with ChatGPT / 10.11)

```text
penalty / PC / apply     → pairwise uncorrected
Contact penalty stencil  → InnerOnly (C4 validated; C1/C2 still explicitly test InnerContact pipeline)
divA posterior diagnostic  → B-corrected + gradDen
global ||divA||/||A||    → not a hard gate
core divA                → validated gate
boundary / interface divA→ open diagnostic
η_A=0.1                  → research primary (under Contact InnerOnly should match Inner order of magnitude)
production λ_A           → off
projection               → defer
Contact A-penalty        → two-body research only; three-body deferred
```

---

## 2. Code change summary

| File | Purpose |
|------|---------|
| `aphi_matrix_free_operator_ck.h/.hpp` | `AphiApplyContactDynamicsBundle` adds `AphiContactPairwiseADivergencePenaltyPipelineBundle` |
| `diagnostics/aphi_contact_a_divergence_penalty_diagnostic_helpers.h` | C1 operator equivalence |
| `diagnostics/aphi_contact_graddiv_pc_diagnostic_helpers.h` | C2 PC consistency |
| `diagnostics/aphi_contact_a_divergence_penalty_two_body_mms_helpers.h` | C3 MMS runner + Joule post |
| `aphi_coupling_modes_ck.h` | C4: `AphiContactADivergencePenaltyStencilMode` + `contact_a_divergence_penalty_stencil` |
| `aphi_matrix_free_operator_ck.h/.hpp` | C4: Contact apply branch InnerContact / InnerOnly penalty pipeline |
| `aphi_block_jacobi_preconditioner_ck.h/.hpp` | C4: skip graddiv accumulation in Contact Jacobi when InnerOnly |
| `diagnostics/aphi_contact_inner_only_penalty_stencil_diagnostic_helpers.h` | C4 comparison runner |
| `test_3d_aphi_ck_contact_a_divergence_penalty_inner_only_stencil_diagnostic/` | C4 |

---

## 3. C1: Contact pairwise apply equivalence — passed=1

Test: `test_3d_aphi_ck_contact_a_divergence_penalty_pairwise_operator_diagnostic`

| Region | Metric | max diff |
|------|------|----------|
| stencil-safe core | divA | ~1.4e-6 |
| stencil-safe core | grad(divA) | ~2.9e-6 |
| stencil-safe core | penalty lhs | ~6.8e-5 |
| interface band (report only) | grad(divA) | ~0.28 |

**Conclusion**: Monolithic Inner vs split Inner+Contact are operator-equivalent in **stencil-safe region**; interface band deviation comes from duplicate particle layers at x=0.5, not an apply implementation bug.

---

## 4. C2: Contact PC consistency — passed=1

Test: `test_3d_aphi_ck_contact_graddiv_pc_consistency_diagnostic`

| Region | max PC diff |
|------|-------------|
| core (8 probe) | ~6.7e-6 |
| interface band (8 probe, report only) | ~43.4 |

**Conclusion**: `JacobiGradDivABlock` (Inner+Contact accumulation) matches FD penalty apply in core.

---

## 5. C3: two-body Contact A-penalty MMS

Test: `test_3d_aphi_ck_contact_a_divergence_penalty_two_body_mms`  
Field: **Az2D**, `lambda_phi=100`, `core_shell=dp`, `eta_A ∈ {0, 0.1, 0.2}`

### 5.1 C3 prototype (InnerContact, 2026-05-21) — superseded by C3+

#### η_A = 0 — strict gate ✅

| Metric | Value |
|------|-----|
| exact_consistency_defect | ~5e-8 |
| converged | 1 |
| global_true_rel | ~2.6e-6 |
| A_error | ~7e-7 |
| E/J/Joule vs exact | ~10⁻⁷ |
| core_divA | ~1.3e-5 |

#### η_A = 0.1 / 0.2 — InnerContact failure ⚠️ (explained by C4)

| Metric | η=0.1 | η=0.2 |
|------|-------|-------|
| exact_consistency_defect | ~1.0% | ~2.0% |
| converged | 0 | 0 |
| global_true_rel | inf | inf |

Root cause: cross-body contamination when penalty uses Inner+Contact divA/gradDivA (C4 fix).

### 5.2 C3+ gate definition and results (InnerOnly, 2026-05-21)

```text
penalty_stencil = InnerOnly
primary strict (η ≤ 0.1): converged + finite true_rel + A/E/J/Joule
upper η=0.2: exact_consistency_defect ≤ 5% (informational)
```

| η_A | true_rel | Result |
|-----|----------|------|
| 0 | 2.5e-6 | ✅ |
| 0.1 | 2.7e-6 | ✅ |
| 0.2 | 7.8e-6 (`polish_sweeps=0`) | ✅ strict (after C5 fix) |

### 5.3 Comparison with Inner monolithic

Inner P4 Az2D η=0.1/0.2 consistency both converged (true_rel ~1e-6).  
C4 proves Contact **InnerOnly** aligns with Inner monolithic; C3+ elevates this to strict gate.

---

## 6. Current status: known vs unknown vs blocked

### 6.1 Resolved / no longer under discussion

| Issue | Status |
|------|------|
| Does Contact matvec include A-penalty | ✅ integrated |
| stencil-safe divA/gradDivA/PC | ✅ C1/C2 |
| η=0 two-body MMS + E/J/Joule | ✅ C3 |
| boundary divA ~27% | ✅ known open (10.11) |
| B=curlA ~3.5% @ dp=0.1 | ✅ 10.12-A converges with dp, informational |

### 6.2 Closed / diagnostic only

| Former open item | Status |
|--------------|------|
| Contact η>0 GMRES non-convergence (InnerContact) | ✅ closed by C4 InnerOnly |
| η=0.2 false convergence / inf true_rel | ✅ C5: `polish_sweeps=0` |
| interface band gradDivA / PC large deviation | diagnostic only (duplicate particle layers) |
| η>0 production | **no**, `lambda_A` production off |

### 6.3 Only escalate to ChatGPT when architecture decisions are needed (**not blocking now**)

| Topic | Why deferred |
|------|----------------|
| Three-body TEAM7 + A-penalty | Plan explicit: do not proceed until C1–C3 pass η>0 converged |
| projection implementation | deferred |
| Interface duplicate layer **geometry** fix vs accept as research limitation | Decide after reviewing 10.12-B+ spatial diagnostic; more efficient to discuss with data |
| cold crucible / 10D thermal | Plan Priority 6, deferred |

**Conclusion: no need to open a ChatGPT round now** — direction is frozen; blockers are **measurable engineering open items**, not unclear concepts.

---

## 7. Next steps (Stage 10.13)

1. ~~Contact two-body research gate~~ ✅ — see §12–14, [`CURSOR_APHI_STAGE10_13_CHATGPT_HANDOFF.md`](CURSOR_APHI_STAGE10_13_CHATGPT_HANDOFF.md)  
2. **P1** cold crucible demo (pipeline smoke, not quantitative validation)  
3. **P4** InnerOnly graddiv PC supplementary test (recommended)  
4. ~~**P4 InnerOnly graddiv PC**~~ ✅ — `test_3d_aphi_ck_contact_inner_only_graddiv_pc_consistency_diagnostic`  
5. **Do not**: three-body A-penalty, production default λ_A, projection implementation

---

## 9. Stage 10.12-B+: interface MMS consistency spatial decomposition — passed=1

Test: `test_3d_aphi_ck_contact_a_divergence_penalty_interface_mms_consistency_diagnostic`  
Method: two Contact applies on exact Az2D; regional stats of ||Ku−b||/||b|| (vol-weighted L2).

| η_A | global | **stencil-safe core** | interface band | boundary |
|-----|--------|----------------------|----------------|----------|
| 0 | ~7e-8 | **~2.6e-7** | ~9e-8 | ~4e-8 |
| 0.1 | ~1.04% | **~2.5×10⁻⁷** | ~0.69% | ~1.14% |
| 0.2 | ~2.13% | **~3.1×10⁻⁷** | ~1.40% | ~2.35% |

(User reproduction 2026-05-21: η=0 global ~3.6e-8; η=0.1 global ~1.04%, interface_max ~0.73; η=0.2 global ~2.13%, interface_max ~1.46)

**Key conclusions**:
1. **stencil-safe core (64 points) still ~1e-7 at η=0.1/0.2** — consistent with C1; penalty in core away from interface **does not** have MMS defect.  
2. **~1–2% global defect almost entirely from interface band + boundary** (interface max single point up to ~0.73).  
3. GMRES η>0 fix direction should focus on **interface/boundary geometry and penalty RHS**, not changing core pairwise operator.  
4. **Still no ChatGPT discussion needed** — data sufficient to guide next engineering steps.

---

## 10. Stage 10.12-B++: interface base vs penalty decomposition — passed=1

Test: `test_3d_aphi_ck_contact_a_divergence_penalty_interface_penalty_breakdown_diagnostic`  
Method: Az2D exact field; `rhs` from **penalty-free** base apply; compare base lhs (`t`) vs full lhs (`v`, with penalty).

| η_A | interface base | interface full | **penalty only** (interface) | penalty only (boundary) |
|-----|----------------|----------------|------------------------------|-------------------------|
| 0 | ~8e-8 | ~9e-8 | ~9e-8 | ~3e-8 |
| 0.1 | **~7e-8** | ~4.6% | **~4.6%** (max ~2.7) | ~8.4% (max ~5.4) |
| 0.2 | ~6e-8 | ~9.3% | **~9.3%** (max ~5.4) | ~17% (max ~11) |

**Interpretation**:
1. **Penalty-free base Contact operator at interface ~1e-7** — again confirms not a Laplace/Contact main operator issue.  
2. **~1% global MMS defect (B+) vs ~4.6% penalty-only (B++)**: B+ uses self-consistent RHS (with penalty); B++ uses η=0 RHS as reference, so penalty increment appears larger.  
3. **η>0 GMRES fix should focus on penalty discrete consistency at interface/boundary**, not base operator.

---

## 11. Stage 10.12-P4: Inner η observable summary table — passed=1

Test: `test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms` (added `eta_observable_matrix` summary row)

Az2D **consistency** (MMS recovers exact, η sweep):

| η_A | conv | true_rel | A_err | B_err | E_comb | Joule | core_divA | boundary_divA |
|-----|------|----------|-------|-------|--------|-------|-----------|---------------|
| 0 | 1 | ~3e-6 | ~3e-6 | ~0.035 | ~6e-7 | ~1e-6 | ~4e-7 | ~0.037 |
| 0.1 | 1 | ~3e-6 | ~5e-6 | ~0.035 | ~1e-6 | ~2e-6 | ~3e-7 | ~0.037 |
| 0.2 | 1 | ~3e-6 | ~2e-6 | ~0.035 | ~2e-6 | ~5e-6 | ~5e-7 | ~0.037 |
| 0.3 | 1 | ~9e-5 | ~3e-4 | ~0.035 | ~3e-4 | ~6e-4 | ~6e-6 | ~0.037 |

CrossSine3D consistency: same trend as Az2D (B ~0.035, E/J ~1e-6, core_divA ~1e-7).

**invariance** (η>0, RHS without penalty): A/B/Joule errors large — informational only (consistent with 10.11).

**Conclusion**: Inner monolithic **η=0.1/0.2/0.3 consistency all converge**; Contact η>0 issue is **split-interface specific**.

---

## 12. Stage 10.12-C4: Inner-only penalty stencil — passed=1 ✅

**Hypothesis**: Contact cross-body neighbor contamination of penalty divA/grad(divA) causes interface ~1–2% MMS defect and η>0 GMRES non-convergence.

**Conclusion**: **Hypothesis fully confirmed**. InnerOnly fixes η=0.1 from ~1% defect + GMRES divergence to ~1e-7 MMS + GMRES convergence (same order as Inner monolithic).

Test: `test_3d_aphi_ck_contact_a_divergence_penalty_inner_only_stencil_diagnostic`

| Metric | η=0 InnerOnly | η=0.1 InnerContact | η=0.1 InnerOnly |
|------|---------------|--------------------|--------------------|
| mms_global_rel_l2 | 7.1e-8 | **1.04%** | **5.2e-8** |
| mms_core_safe_rel_l2 | 2.0e-7 | 3.2e-7 | 2.3e-7 |
| mms_interface_rel_l2 | — | **0.69%** | **6.0e-8** |
| mms_boundary_rel_l2 | — | 1.14% | 3.9e-8 |
| gmres_converged | — | **0** | **1** |
| gmres_exact_consistency_defect | — | **~1.0%** | **4.9e-8** |
| gmres_global_true_rel | — | inf | **2.9e-6** |

gate (Stage 10.13 P0 hardening): `eta0_regression_ok` + `research_improved` (global ratio ≥1e3 + GMRES converged) + `inner_only_eta01_ok` (mms_global ≤1e-6, gmres true_rel ≤1e-5) → `passed=1`

**Root cause localization** (consistent with B++, C4 provides fix):
- base Contact operator at interface ~1e-7 (C1/C2 proven)
- **when penalty uses Inner+Contact divA/gradDivA at split interface**, cross-body neighbors introduce ~1% discrete inconsistency
- **InnerOnly penalty stencil** eliminates cross-body contamination; Contact η=0.1 behavior aligns with Inner monolithic

**Engineering implications**:
- Contact research path: **A-penalty should default to `InnerOnly` stencil** (main K(X) still Inner+Contact)
- **No ChatGPT needed** — data sufficient to decide research default
- Next: C3+ strict gate user run confirmation; C1/C2 still explicitly InnerContact

---

## 13. Stage 10.12-C3+: InnerOnly strict gate — passed=1 (primary η)

**User reproduction** (2026-05-21, `penalty_stencil=InnerOnly`):

| η_A | true_rel | A_err | E/J/Joule | gate |
|-----|----------|-------|-----------|------|
| 0 | 2.7e-6 | 3.3e-6 | ~10⁻⁷ | **strict ✅** |
| 0.1 | 2.7e-6 | 3.8e-6 | ~10⁻⁶ | **strict ✅** |
| 0.2 | **7.8e-6** | 4.9e-6 | ~10⁻⁶ | **strict ✅** (after C5 polish fix) |

After C5 fix rerun (2026-05-21): η=0.2 `global_true_rel=7.79e-6` (was inf), consistent with C5 `polish=0` prediction ~7.8e-6.

`primary_rows_passed=2/2 upper_eta_reported=1/1 passed=1`

**C3+ gate (final)**:

```text
primary strict (η ≤ 0.1): converged + finite true_rel + observables
upper η=0.2: finite exact_consistency_defect and ≤ 5% (informational)
```

User result: `primary_rows_passed=2/2`, `upper_eta_reported=1/1` → **passed=1** (reconfirmed 2026-05-21)

**C1/C2 regression**: pairwise_operator `passed=1`; graddiv_pc `passed=1` (explicit InnerContact, unaffected by default change).

**Conclusion**: Contact two-body **research primary gate (η=0 + η=0.1, InnerOnly) closed**.

**Closed**: η=0.2 requires `polish_sweeps=0` (C5); not a penalty stencil issue

---

## 14. Stage 10.12-C5: upper-η solver diagnostic — passed=1 ✅

**ninja: no work to do** = binaries already built in prior session, **output format correct**, this is the C5 case.

Test: `test_3d_aphi_ck_contact_a_divergence_penalty_upper_eta_solver_diagnostic`

| Config | GMRES recursive/true | global_true_rel | mms_finite |
|------|---------------------|-----------------|------------|
| η=0.1 polish=1 | ~2.8e-6 / gap=0 | 2.8e-6 | ✅ |
| η=0.2 polish=**0** | ~7.8e-6 / gap=0 | **7.8e-6** | ✅ |
| η=0.2 polish=**1** | ~7.8e-6 / gap=0 | **inf** | ❌ |

**Polish curve (η=0.2)**:

| sweep | global_true_rel |
|-------|-----------------|
| 0 | 7.8e-6 ✅ |
| 1 | 9.97e+5 |
| 2 | 1.36e+17 |
| 3–5 | inf |

**Conclusion (root cause closed)**:
- **Not GMRES false convergence** — η=0.2 at `polish_sweeps=0` has true_rel ~8e-6, same order as Inner monolithic
- **Block-GS polish sweep destroys solution at upper η** — C3+ η=0.2 failure due to default `polish_sweeps=1`
- η=0.1 + polish=1 still stable; η≥0.2 should **disable polish** (fixed in MMS helper)

**Engineering fix** (`runContactDivFreeTwoBodyMmsRow`):
```text
η < 0.2  → polish_sweeps = 1
η ≥ 0.2  → polish_sweeps = 0
```

**Next step**: user reruns C3, expect η=0.2 strict also passes. → **Confirmed** (true_rel 7.8e-6)

---

## 8. Reproduction

```bash
cd build
ninja test_3d_aphi_ck_contact_a_divergence_penalty_pairwise_operator_diagnostic \
     test_3d_aphi_ck_contact_graddiv_pc_consistency_diagnostic \
     test_3d_aphi_ck_contact_a_divergence_penalty_two_body_mms \
     test_3d_aphi_ck_contact_a_divergence_penalty_interface_mms_consistency_diagnostic \
     test_3d_aphi_ck_contact_a_divergence_penalty_interface_penalty_breakdown_diagnostic \
     test_3d_aphi_ck_contact_a_divergence_penalty_inner_only_stencil_diagnostic \
     test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms

./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_contact_a_divergence_penalty_two_body_mms/bin/test_3d_aphi_ck_contact_a_divergence_penalty_two_body_mms
```
