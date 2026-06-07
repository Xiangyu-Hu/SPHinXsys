# Stage 10.10 Validation Results + Open Questions (ChatGPT Discussion Pack)

> **Date**: 2026-05-21  
> **Status**: Stage 10.10-A/B/B2/C implemented and run on Cursor side; **pairwise global rel≈0.13 root cause located**; align with ChatGPT on boundary open item, dual-track architecture, η_A freeze, and Stage 10.11 priority  
> **Prior plan**: [`stage10_10_divergence_free_mms_and_3d_pairwise_plan_for_cursor.md`](stage10_10_divergence_free_mms_and_3d_pairwise_plan_for_cursor.md)  
> **Cursor main record**: [`tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md)  
> **10.9 discussion pack (background)**: [`stage10_9_validation_results_and_open_questions_for_chatgpt.md`](stage10_9_validation_results_and_open_questions_for_chatgpt.md)

---

## 0. One sentence for ChatGPT

Stage 10.10 completed per your plan: **non-degenerate div-free MMS** (Az2D/CrossSine3D), **four-field pairwise divA sweep**, **root-cause diagnostic** (core/boundary + gradDen denominator).  
**Key reversal**: Az2D/CrossSine3D global pairwise rel≈0.13 is **not core discretization failure**—core region ~1e-7, error 100% at boundary; global rel amplified ~10⁶× by `||A||` denominator.  
Please help decide: **(1) formally adopt pairwise penalty + B-corrected diagnostic dual-track? (2) must boundary divA O(0.04) be fixed? (3) can η_A=0.1 be frozen? (4) Stage 10.11 priority?**

---

## 1. Completed vs 10.10 Plan Task Mapping

| 10.10 Task | Cursor Status | Test / Record |
|------------|------------|-------------|
| 10.10-A non-degenerate MMS | ✅ | `test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms` |
| 10.10-B four-field pairwise sweep | ✅ | `test_3d_aphi_ck_pairwise_diva_3d_diagnostic` |
| 10.10-B2 root-cause diagnostic | ✅ new | `test_3d_aphi_ck_pairwise_diva_root_cause_diagnostic` |
| 10.10-C observable exact | ✅ partial | Joule + **E_combined** vs exact; `E_L2` still false explosion |
| 10.10-D η_A policy freeze | ⏸ tentative 0.1 primary | discuss Q4 |
| 10.10-E Contact unlock | ❌ not met | still deferred |

---

## 2. Key Numerical Summary (for discussion)

### 2.1 non-degenerate MMS (10.10-A, dp=0.1)

| Field | η_A | rhs_norm | consistency | block_Linf | joule_err_vs_exact | E_combined_err_vs_exact |
|-------|-----|----------|-------------|------------|-------------------|------------------------|
| Az2D | 0 | 22.4 | ✅ | 6e-7 | ~2e-7 | ~1e-6 |
| Az2D | 0.1 | 21.8 | ✅ | 6e-7 | ~2e-7 | ~1e-6 |
| CrossSine3D | 0 | 38.9 | ✅ | 8e-6 | ~1e-6 | — |
| CrossSine3D | 0.3 | 36.2 | ✅ | 4.8e-4 | 6e-4 | — |

- consistency: **8/8** (two fields × η=0,0.1,0.2,0.3)  
- invariance GMRES: **6/6** converged; but block_Linf can reach 0.2–2.3 (boundary penalty nonzero)  
- **10.9 Q1 resolved**: Az2D/CrossSine3D `rhs_norm≠0`, no longer degenerate  

### 2.2 four-field pairwise vs B-corrected (10.10-B, dp=0.1)

| Field | analytic divA | pairwise rel(Anorm) | B-corrected rel |
|-------|-----------|---------------------|-----------------|
| Linear2D | 0 | **2.4e-9** | 1.2e-8 |
| Az2D | 0 | **0.129** | 0.037 |
| CrossSine3D | 0 | **0.129** | 0.037 |
| Sinusoidal3DCurlPsi | 0 | **0.477** (info) | 0.011 |

dp refinement (pairwise Anorm rel): Az2D 0.194→0.129→0.095 (dp=0.15/0.10/0.075), **did not converge below 1e-2**.

### 2.3 root-cause diagnostic (10.10-B2, dp=0.1, core_shell=2.5 dp)— **most important**

| Field | rel(Anorm denom) | rel(gradA denom) | core-only rel | boundary energy frac |
|-------|----------------|-----------------|---------------|---------------------|
| Linear2D | 2.4e-9 | 3.5e-10 | 6.8e-10 | 0.97 |
| Az2D | **0.129** | **1.15e-7** | **1.78e-7** | **1.00** |
| CrossSine3D | **0.129** | **1.07e-7** | **1.66e-7** | **1.00** |

- core_shell sweep (2.5/3.0/3.5 dp): **global Anorm rel constant 0.129**; core-only rel constant ~1e-7  
- B-corrected gradDen rel ~**0.037** (boundary trace(gradA) still O(0.04) error)

### 2.4 10.9 Case C background (still valid)

- Solenoidal η=0.1: div_A↓, Joule/E vs **same dp η=0** perturbation 7–14%  
- Fine-dp self-ref vs dp=0.1: **Joule gap ~93%** → cannot serve as penalty hard gate  
- MMS track: observable vs η=0 ref <0.05%, but div_A≈3.7 (non div-free field)

### 2.5 Cursor-side architecture recommendation (pending ChatGPT confirmation)

| Purpose | Recommendation |
|------|------|
| penalty operator / PC | pairwise uncorrected (same family as GMRES apply) |
| a posteriori divA reporting | B-corrected + gradDen denominator + core/boundary decomposition |
| hard gate | **do not use** global Anorm rel; use core-only pairwise or B-corrected |

---

## 3. Open Questions (please advise line by line)

### Q1. 10.9 Q1 resolved — need other div-free fields?

Az2D `A=(0,0,sinπx sinπy)`, CrossSine3D give `rhs_norm≈22–39`, MMS 8/8 passed.  
**Question**: Are these two fields sufficient as **absolute validation benchmarks**? Or must Sinusoidal3DCurlPsi-type curl-ψ fields also reach pairwise ~1e-9 to close?

### Q2. pairwise global rel≈0.13 — root cause located, how to freeze architecture?

Cursor found:
- core region pairwise ~1e-7 (discretization correct)
- global rel dominated by `||A||` denominator + boundary 100% energy

**Questions**:
1. **Formally adopt** "pairwise penalty + B-corrected diagnostic dual-track"?  
2. In gate / records / external messaging, **deprecate** global `div_A_rel = div_L2 / ||A||_L2` as hard metric?  
3. Is Sinusoidal3DCurlPsi pairwise ~0.48 still informational, or needs same root-cause decomposition?

### Q3. boundary divA O(0.04) — must it be fixed?

- B-corrected gradDen rel ~0.037 stable on Az2D/CrossSine3D  
- invariance MMS: GMRES converges but block error large (exact field divA≠0 at boundary → penalty changes solution)  
- boundary divA L2 energy fraction **100%**

**Questions**:
1. Is boundary error expected SPH free-surface/kernel truncation behavior, or should it be fixed (dedicated BC, boundary exclusion, corrected kernel)?  
2. Before boundary fix, should **invariance MMS** be permanently downgraded to informational?  
3. If only diagnostic fixed not penalty, can research prototype claim "convergent"?

### Q4. Can η_A=0.1 be formally frozen?

Combined 10.9 Case C + 10.10 MMS exact observable:

| Evidence | Supports η=0.1 |
|------|-----------|
| Solenoidal gate div_A↓ | ✅ |
| vs same dp η=0 perturbation 7–14% | ✅ research acceptable |
| Az2D Joule/E_combined vs exact ~1e-6–1e-7 | ✅ |
| invariance block error large | ⚠️ boundary open |
| CrossSine3D η=0.3 joule_err_vs_exact ~6e-4 | ⚠️ |

**Question**: Freeze `primary_eta_a=0.1, optional=0.2, upper=0.3` and write to `AphiADivergencePenaltyResearchDefaults`? Or need boundary fix or full η-sweep vs exact E/J/Joule matrix first?

### Q5. E observable — is E_combined sufficient?

- `E_L2` (real components only) false explosion (~10¹⁹) when `A_imag=0, φ=0`  
- new `E_combined = sqrt(∫ |E_real|²+|E_imag|² dV)` → Az2D η=0 err ~1e-6  
- Joule vs exact ~2e-7

**Question**: Should Case C / MMS reports uniformly use **E_combined + Joule** as observable pair? Need B_exact / J_exact analytic reference too?

### Q6. Stage 10.11 priority (please rank)

1. **boundary divA dedicated study** (why B-corrected ~0.04; exclusion / BC?)  
2. **documentation/external messaging freeze** (dual-track architecture + η_A policy + Part V update)  
3. **η-sweep vs exact observable full matrix** (Az2D/CrossSine3D × η=0,0.1,0.2,0.3)  
4. **Sinusoidal3DCurlPsi root-cause decomposition** (same Az2D pattern?)  
5. **Contact A-penalty prototype** (10.10-E condition: boundary + absolute validation)  
6. **impressed case / cold crucible** informational gate  

### Q7. External messaging (Part V update)

Current Cursor understanding:

**Can say**:
- non-degenerate div-free MMS established (Az2D/CrossSine3D, `rhs_norm≠0`)  
- Inner pairwise A-penalty consistency MMS converges  
- Joule / E_combined vs analytic exact ~1e-6–1e-7 (Az2D η=0)  
- pairwise divA **core region** ~1e-7; global Anorm rel is denominator+boundary artifact  
- η_A=0.1 primary **research** candidate  

**Cannot say**:
- Coulomb gauge enforced  
- pairwise divA globally validated  
- invariance div-free physical solution unchanged  
- production ready / Contact A-penalty can open  

**Question**: After 10.10 completion, do these boundaries need revision? Can we say "3D div-free pairwise discretization validated in core region"?

---

## 4. Recommended Upload File List (by priority)

### 4.1 Must read (discussion decisions)— minimum 3

| # | File | Description |
|---|------|------|
| 1 | [`stage10_10_validation_results_and_open_questions_for_chatgpt.md`](stage10_10_validation_results_and_open_questions_for_chatgpt.md) | **This document** |
| 2 | [`tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md) | **Cursor 10.10 main record** |
| 3 | [`stage10_10_divergence_free_mms_and_3d_pairwise_plan_for_cursor.md`](stage10_10_divergence_free_mms_and_3d_pairwise_plan_for_cursor.md) | Your original 10.10 plan (completion mapping) |

### 4.2 Strongly recommended (background closure)

| # | File | Description |
|---|------|------|
| 4 | [`stage10_9_validation_results_and_open_questions_for_chatgpt.md`](stage10_9_validation_results_and_open_questions_for_chatgpt.md) | 10.9 discussion pack (Q1–Q6 background) |
| 5 | [`tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_8_PAIRWISE_DIVA_RECORD.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_8_PAIRWISE_DIVA_RECORD.md) | Stage 10.8/10.9 main record §13–§15 |
| 6 | [`stage10_8_metric_interpretation_and_validation_plan_for_cursor.md`](stage10_8_metric_interpretation_and_validation_plan_for_cursor.md) | Original validation plan |

### 4.3 New code (10.10 implementation)— optional but helps understand metrics

| # | File | Description |
|---|------|------|
| 7 | [`tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_divergence_free_nonzero_rhs_mms_helpers.h`](tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_divergence_free_nonzero_rhs_mms_helpers.h) | Az2D/CrossSine3D MMS + E_combined |
| 8 | [`tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_pairwise_diva_root_cause_helpers.h`](tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_pairwise_diva_root_cause_helpers.h) | Root-cause diagnostic (core/boundary + gradDen) |
| 9 | [`tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_div_a_diagnostic_helpers.h`](tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_div_a_diagnostic_helpers.h) | divA metric definitions (Anorm vs gradDen denominator) |

### 4.4 New tests (all passed=1)— usually no need to upload source

| Test | Phase | Description |
|------|-------|------|
| [`test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms/`](tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms/) | 10.10-A/C | MMS + exact observable |
| [`test_3d_aphi_ck_pairwise_diva_3d_diagnostic/`](tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_pairwise_diva_3d_diagnostic/) | 10.10-B | four-field + dp sweep |
| [`test_3d_aphi_ck_pairwise_diva_root_cause_diagnostic/`](tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_pairwise_diva_root_cause_diagnostic/) | 10.10-B2 | root-cause diagnostic |

### 4.5 Minimal upload strategy

| Scenario | Upload |
|------|------|
| **Minimal pack (recommended)** | #1 + #2 + #3 |
| **Full discussion** | #1–#6 |
| **If ChatGPT needs metric implementation** | add #7–#9 |

---

## 5. Cursor-side Work Log (for verbal summary to ChatGPT)

### 5.1 New since 10.9 discussion pack

1. **10.10-A**: Az2D / CrossSine3D non-degenerate div-free MMS; `rhs_norm≈22–39`; consistency 8/8  
2. **10.10-B**: four-field pairwise vs B-corrected + dp sweep  
3. **10.10-B2**: root-cause diagnostic — found global rel≈0.13 is **denominator artifact + boundary dominated**, core ~1e-7  
4. **10.10-C**: `E_combined` observable; Joule/E_combined vs exact ~1e-6–1e-7  
5. **Record**: [`CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md)

### 5.2 Status of 10.9 open questions

| 10.9 Q | 10.10 Status |
|--------|-----------|
| Q1 div-free MMS degeneracy | ✅ Az2D/CrossSine3D resolved |
| Q2 3D pairwise ~0.48 | ⚠️ partial: Az2D/CrossSine3D root cause located; Sinusoidal not decomposed |
| Q3 Case C boundaries | ⏸ still open; 10.10 added exact observable |
| Q4 η=0.1 freeze | ⏸ pending ChatGPT confirmation |
| Q5 priority | 10.10 executed; 10.11 pending discussion |
| Q6 external messaging | needs Q7 update |

---

## 6. Reproduction Commands

```bash
cd /home/yongchuan/sphinxsys/build
cmake .. -DCMAKE_BUILD_TYPE=Release
ninja test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms \
      test_3d_aphi_ck_pairwise_diva_3d_diagnostic \
      test_3d_aphi_ck_pairwise_diva_root_cause_diagnostic

BIN=tests/extra_source_and_tests/3d_examples
for t in inner_divergence_free_a_nonzero_rhs_mms \
         pairwise_diva_3d_diagnostic \
         pairwise_diva_root_cause_diagnostic; do
  echo "=== test_3d_aphi_ck_${t} ==="
  "$BIN/test_3d_aphi_ck_${t}/bin/test_3d_aphi_ck_${t}" 2>&1 | grep -E 'passed=|az2d_|grad_den'
done
```

---

## 7. Suggested Prompt to Paste for ChatGPT (copy directly)

```text
You are the technical advisor for A-phi SPHinXsys Stage 10. Cursor completed Stage 10.10-A/B/B2/C per stage10_10_divergence_free_mms_and_3d_pairwise_plan (see stage10_10_validation_results_and_open_questions_for_chatgpt.md and CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md).

10.10 key new findings (please read §2.3):
- Az2D/CrossSine3D global pairwise divA rel≈0.13 is not core discretization failure
- core region pairwise rel ~1e-7; after gradA denominator global rel ~1e-7
- divA L2 energy 100% at boundary; B-corrected gradDen rel ~0.037
- non-degenerate MMS: rhs_norm≈22–39, consistency 8/8, Joule/E_combined vs exact ~1e-6–1e-7

Based on attachments:
1. Give clear advice on Q1–Q7 line by line (including priority ranking);
2. Confirm or reject "pairwise penalty + B-corrected diagnostic dual-track" architecture;
3. Output Stage 10.11 execution plan (task names, acceptance criteria, suggested test names);
4. Update Part V external messaging (can say / cannot say);
5. Agree or not to formally freeze η_A=0.1 as research primary default;
6. Is boundary divA O(0.04) blocking, or acceptable as SPH boundary artifact.

10.9 resolved: div-free MMS degeneracy (Az2D/CrossSine3D).
10.10 still open: boundary divA, invariance MMS block error, Sinusoidal3DCurlPsi ~0.48 not decomposed.
```

---

## 8. Quick Reference Links

| Document | Path |
|------|------|
| **This discussion pack** | [`stage10_10_validation_results_and_open_questions_for_chatgpt.md`](stage10_10_validation_results_and_open_questions_for_chatgpt.md) |
| **10.10 Cursor main record** | [`tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md) |
| **10.10 original plan** | [`stage10_10_divergence_free_mms_and_3d_pairwise_plan_for_cursor.md`](stage10_10_divergence_free_mms_and_3d_pairwise_plan_for_cursor.md) |
| **10.9 discussion pack** | [`stage10_9_validation_results_and_open_questions_for_chatgpt.md`](stage10_9_validation_results_and_open_questions_for_chatgpt.md) |
| **10.8/10.9 main record** | [`tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_8_PAIRWISE_DIVA_RECORD.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_8_PAIRWISE_DIVA_RECORD.md) |
| **Source file map** | [`tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/STAGE10_APHI_FILE_MAP_FOR_DISCUSSION.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/STAGE10_APHI_FILE_MAP_FOR_DISCUSSION.md) |
