# Stage 10.11 Validation Results + Open Questions (ChatGPT Discussion Pack)

> **Date**: 2026-05-21  
> **Status**: Stage 10.11-A/B/C/D/E completed and run on Cursor side; **boundary divA localized**; **B curl ~3.5% attributed**; align with ChatGPT on Contact unlock conditions, projection timing, and Stage 10.12 priority  
> **Prior plan**: [`stage10_11_boundary_divA_B_E_H_projection_eta_plan_for_cursor.md`](stage10_11_boundary_divA_B_E_H_projection_eta_plan_for_cursor.md)  
> **Cursor main record**: [`tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_11_BOUNDARY_VALIDATION_RECORD.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_11_BOUNDARY_VALIDATION_RECORD.md)  
> **10.10 discussion pack (background)**: [`stage10_10_validation_results_and_open_questions_for_chatgpt.md`](stage10_10_validation_results_and_open_questions_for_chatgpt.md)

---

## 0. One sentence for ChatGPT

Stage 10.11 completed per your plan: **boundary divA dedicated study**, **B/E/J/Joule/H exact observable**, **ghost/buffer boundary experiments**, **projection design doc**, **B curl dual-track diagnostic**.  
**Key conclusions**:
1. core pairwise divA ~1e-7 validated; boundary divA is support deficiency (100% energy at boundary)  
2. Joule / E_combined / J_combined vs exact ~1e-6–1e-7 (Az2D η=0 MMS)  
3. B=curlA ~3.5% is **B-corrected grad + oscillatory exact A field** discretization limit, **not** MMS failure, **not** pairwise vs B-corrected operator-family choice (pairwise grad curl error ~200%)  
4. ghost only significantly helps Sinusoidal3DCurlPsi; Az2D/CrossSine lattice buffer sufficient  
5. η_A=0.1 written to research defaults  

Please help decide: **(1) Contact A-penalty unlock conditions met? (2) projection still defer? (3) B ~3.5% blocking? (4) Stage 10.12 mainline?**

---

## 1. Completed vs 10.11 Plan Task Mapping

| 10.11 Task | Cursor Status | Test / Record |
|------------|------------|-------------|
| 10.11-A boundary divA dedicated study | ✅ passed=1 | `test_3d_aphi_ck_pairwise_diva_boundary_diagnostic` |
| 10.11-B B/E/J/Joule/H exact | ✅ passed=1 | `test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms` (extended B/J gate) |
| 10.11-C ghost/buffer experiments | ✅ passed=1 | `test_3d_aphi_ck_ghost_buffer_diva_diagnostic` |
| 10.11-D projection design | ✅ design only | `CURSOR_APHI_STAGE10_PROJECTION_DESIGN_NOTES.md` |
| 10.11-E B curl dual-track (Cursor added) | ✅ passed=1 | `test_3d_aphi_ck_curl_b_dual_track_diagnostic` |
| η_A=0.1 freeze | ✅ | `aphi_coupling_modes_ck.h` |
| Contact A-penalty | ❌ still deferred | see §3 Q4 |

---

## 2. Key Numerical Summary (for discussion)

### 2.1 boundary divA (10.11-A, dp=0.1, core_shell=2.5 dp)

| Region | pairwise gradDen rel | B-corrected gradDen rel |
|------|---------------------|-------------------------|
| core | **~1e-7** | **~1e-7** |
| boundary | **~0.14** | **~0.037** |
| divA energy fraction | — | **100% boundary** |

boundary_width sweep (2/3/4 dp): core metrics unchanged.

### 2.2 B/E/J/Joule exact observable (10.11-B, Az2D η=0 MMS)

| Metric | vs exact | Notes |
|------|----------|------|
| Joule | ~1e-6 | ✅ |
| E_combined | ~1e-7 | ✅ |
| J_combined | ~1e-7 | ✅ |
| **B=curlA** | **~0.035** | ⚠️ oscillatory field curl discretization limit |
| H (=νB) | ~0.035 | constant ν |
| divA core gradDen | ~5e-7 | ✅ |
| divA boundary gradDen | ~0.037 | ⚠️ open |

invariance (η>0, b_0): GMRES converges but block/B/Joule error large — **informational only**.

### 2.3 ghost/buffer (10.11-C, supplemental exterior analytic ghost)

| Field | baseline boundary gradDen | + ghost | Conclusion |
|-------|--------------------------|---------|------|
| Az2D | ~0.031 | ~0.52 (worse) | lattice sufficient |
| CrossSine3D | ~0.031 | ~0.52 (worse) | same |
| Sinusoidal3DCurlPsi | ~0.098 | **~0.045 (-54%)** | support deficiency confirmed |

### 2.4 B curl dual-track (10.11-E)

| Field | B-corrected core rel | Pairwise uncorrected core rel |
|-------|---------------------|------------------------------|
| Linear2D | **~1e-7** | ~2.0 |
| Az2D | **~0.035** | ~2.0 |
| CrossSine3D | **~0.035** | ~2.0 |

**Interpretation**: B post-processing must use B-corrected grad; ~3.5% is oscillatory field discretization error, not operator-family choice issue.

### 2.5 η_A policy (written to code)

```text
primary_eta_a   = 0.1
optional_eta_a  = 0.2
upper_eta_a     = 0.3
production      = lambda_A off
```

### 2.6 Dual-track architecture (executed on Cursor side, pending ChatGPT formal freeze)

| Purpose | Discretization |
|------|--------|
| penalty operator / PC / GMRES apply | pairwise uncorrected |
| a posteriori divA reporting | B-corrected + gradDen denominator + core/boundary decomposition |
| a posteriori B=curlA reporting | **B-corrected grad only** (10.11-E confirmed) |
| hard gate | core-only; **forbid** global Anorm rel |

---

## 3. Open Questions (please advise line by line)

### Q1. How to freeze 10.10 Q2/Q3 after 10.11?

10.11 confirms:
- boundary divA 100% energy; core ~1e-7  
- B-corrected boundary gradDen ~0.037 stable  
- ghost no benefit for Az2D/CrossSine (over-correction)

**Questions**:
1. **Formally close** gates like "global divA must < 1e-2"?  
2. Is boundary ~0.04 **permanently accepted** as SPH support artifact (research prototype)?  
3. Should external messaging adopt: "divA validated in core; boundary is open diagnostic issue"?

### Q2. B=curlA ~3.5% — blocking or informational?

- exact A field + B-corrected curl → ~3.5% vs analytic B (Az2D/CrossSine3D)  
- Linear2D solenoidal → ~1e-7  
- Joule/E/J vs exact ~1e-6 shows **physical solution (E/J path) accurate**  
- pairwise grad curl → ~200% (operator-family issue ruled out)

**Questions**:
1. Is ~3.5% **acceptable** as SPH curl post-processing limit?  
2. Need dp refinement / higher-order gradient experiments before claiming "B validated"?  
3. How should B acceptance criteria be defined in cold crucible / impressed case (vs exact? vs η=0?)?

### Q3. ghost/buffer strategy — enter production diagnostic?

**Questions**:
1. Enable supplemental ghost only for **support-deficient** fields (e.g. Sinusoidal3DCurlPsi)?  
2. **Forbid** extra ghost for Az2D/CrossSine (over-correction)?  
3. Need dummy particle prototype in Contact / wall case, or continue defer?

### Q4. Contact A-penalty unlock — conditions met?

Your 10.11 §20 conditions:

| Condition | 10.11 Status |
|------|-----------|
| core divA validated | ✅ ~1e-7 |
| boundary divA explained or has strategy | ✅ support artifact + ghost experiments |
| B/E/J/Joule exact validation passed | ⚠️ B ~3.5%; E/J/Joule ~1e-6 |
| η_A=0.1 no unacceptable deviation on exact observable | ✅ Az2D η=0.1 MMS |
| projection not needed yet or design exists | ✅ design only |
| production baseline still stable | ✅ λ_A off |

**Question**: **Still not met** for Contact unlock? Which item missing? If unlock, what should first test be?

### Q5. projection — still defer or Stage 10.12 post-solve prototype?

Design doc [`CURSOR_APHI_STAGE10_PROJECTION_DESIGN_NOTES.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_PROJECTION_DESIGN_NOTES.md) answered 7 open questions (favor pairwise Laplace + post-solve).

**Questions**:
1. With core divA already ~1e-7, is projection **still unnecessary**?  
2. If prototype, test **only** on Sinusoidal / boundary-deficient cases?  
3. Is φ synchronized update (keeping E unchanged) blocking?

### Q6. η_A=0.1 — can it be formally frozen as research primary?

Cursor wrote `AphiADivergencePenaltyResearchDefaults::primary_eta_a = 0.1`.

**Questions**:
1. Need **full matrix** η-sweep × (Az2D, CrossSine3D) × (Joule, E, B, J) before freeze?  
2. Does CrossSine3D η=0.3 joule_err ~6e-4 limit upper bound?  
3. Is production **permanently** λ_A off until Contact unlock?

### Q7. invariance MMS — permanently informational?

η>0 + exact div-free A → boundary penalty changes solution → block/B error large, but GMRES converges.

**Question**: **Permanently** not a hard gate? Rename to "penalty sensitivity diagnostic"?

### Q8. Stage 10.12 priority (please rank)

Cursor candidates:

1. **Contact A-penalty prototype** (if Q4 unlock)  
2. **dp refinement B curl** (does ~3.5% decrease with dp?)  
3. **curlH vs J** exact validation (10.11 plan optional)  
4. **impressed / cold crucible informational gate**  
5. **projection post-solve prototype** (if Q5 agrees)  
6. **documentation/Part V external messaging formal freeze**  
7. **10D thermal coupling** — still defer?

---

## 4. Recommended Upload File List (by priority)

### 4.1 Must read (discussion decisions)— minimum 4

| # | File | Description |
|---|------|------|
| 1 | [`stage10_11_validation_results_and_open_questions_for_chatgpt.md`](stage10_11_validation_results_and_open_questions_for_chatgpt.md) | **This document** |
| 2 | [`tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_11_BOUNDARY_VALIDATION_RECORD.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_11_BOUNDARY_VALIDATION_RECORD.md) | **Cursor 10.11 main record** |
| 3 | [`stage10_11_boundary_divA_B_E_H_projection_eta_plan_for_cursor.md`](stage10_11_boundary_divA_B_E_H_projection_eta_plan_for_cursor.md) | Your original 10.11 plan (completion mapping) |
| 4 | [`tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_PROJECTION_DESIGN_NOTES.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_PROJECTION_DESIGN_NOTES.md) | 10.11-D projection design |

### 4.2 Strongly recommended (10.10→10.11 closure)

| # | File | Description |
|---|------|------|
| 5 | [`stage10_10_validation_results_and_open_questions_for_chatgpt.md`](stage10_10_validation_results_and_open_questions_for_chatgpt.md) | 10.10 discussion pack (Q2/Q3 background) |
| 6 | [`tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md) | 10.10 main record |
| 7 | [`stage10_10_divergence_free_mms_and_3d_pairwise_plan_for_cursor.md`](stage10_10_divergence_free_mms_and_3d_pairwise_plan_for_cursor.md) | 10.10 original plan |

### 4.3 New/modified source files (10.11 implementation)

| # | File | Description |
|---|------|------|
| 8 | [`tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_pairwise_diva_boundary_helpers.h`](tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_pairwise_diva_boundary_helpers.h) | 10.11-A boundary/core decomposition |
| 9 | [`tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_ghost_buffer_diva_diagnostic_helpers.h`](tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_ghost_buffer_diva_diagnostic_helpers.h) | 10.11-C supplemental ghost |
| 10 | [`tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_curl_b_dual_track_diagnostic_helpers.h`](tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_curl_b_dual_track_diagnostic_helpers.h) | 10.11-E B curl dual-track pipeline |
| 11 | [`tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_divergence_free_nonzero_rhs_mms_helpers.h`](tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_divergence_free_nonzero_rhs_mms_helpers.h) | 10.11-B extended B/J/divA core-boundary |
| 12 | [`tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_pairwise_diva_root_cause_helpers.h`](tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_pairwise_diva_root_cause_helpers.h) | 10.10-B2 root cause (10.11 reused) |
| 13 | [`tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_div_a_diagnostic_helpers.h`](tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_div_a_diagnostic_helpers.h) | divA dual-track mode definitions |
| 14 | [`tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/aphi_pairwise_div_a_ck.h`](tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/aphi_pairwise_div_a_ck.h) | + `AphiPairwiseVectorGradientCK` |
| 15 | [`tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/aphi_pairwise_div_a_ck.hpp`](tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/aphi_pairwise_div_a_ck.hpp) | pairwise vector grad implementation |
| 16 | [`tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/aphi_coupling_modes_ck.h`](tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/aphi_coupling_modes_ck.h) | η_A=0.1 research defaults |

### 4.4 New test cases (all passed=1)

| # | Test directory | Phase | Description |
|---|----------|-------|------|
| 17 | [`test_3d_aphi_ck_pairwise_diva_boundary_diagnostic/`](tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_pairwise_diva_boundary_diagnostic/) | 10.11-A | boundary divA dedicated study |
| 18 | [`test_3d_aphi_ck_ghost_buffer_diva_diagnostic/`](tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_ghost_buffer_diva_diagnostic/) | 10.11-C | ghost/buffer experiments |
| 19 | [`test_3d_aphi_ck_curl_b_dual_track_diagnostic/`](tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_curl_b_dual_track_diagnostic/) | 10.11-E | B curl dual-track |
| 20 | [`test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms/`](tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms/) | 10.11-B | MMS + B/E/J exact (extended) |

### 4.5 10.10 background tests (optional)

| # | Test directory | Description |
|---|----------|------|
| 21 | [`test_3d_aphi_ck_pairwise_diva_root_cause_diagnostic/`](tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_pairwise_diva_root_cause_diagnostic/) | root-cause diagnostic |
| 22 | [`test_3d_aphi_ck_pairwise_diva_3d_diagnostic/`](tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_pairwise_diva_3d_diagnostic/) | four-field sweep |
| 23 | [`test_3d_aphi_ck_curl_a_manufactured_diagnostic/`](tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_curl_a_manufactured_diagnostic/) | Linear2D solenoidal B |

### 4.6 Minimal upload strategy

| Scenario | Upload |
|------|------|
| **Minimal pack (recommended)** | #1 + #2 + #3 + #4 |
| **Full 10.11 discussion** | #1–#4 + #5–#7 |
| **If ChatGPT needs implementation** | add #8–#16 |
| **If ChatGPT needs test gates** | add #17–#20 |

---

## 5. Cursor-side Work Log (for verbal summary to ChatGPT)

### 5.1 New since 10.10 discussion pack (10.11)

1. **10.11-A**: boundary/core divA decomposition; boundary 100% energy; width sweep does not affect core  
2. **10.11-B**: MMS extended B/H/J_combined vs exact; confirmed Joule/E ~1e-6, B ~3.5%  
3. **10.11-C**: supplemental ghost — Sinusoidal improved 54%; Az2D/CrossSine over-correction  
4. **10.11-D**: projection design doc (no implementation)  
5. **10.11-E**: B curl dual-track — B-corrected ~3.5%; pairwise ~200%  
6. **η_A=0.1** written to `AphiADivergencePenaltyResearchDefaults`  
7. **Record**: [`CURSOR_APHI_STAGE10_11_BOUNDARY_VALIDATION_RECORD.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_11_BOUNDARY_VALIDATION_RECORD.md)

### 5.2 Status of 10.10 open questions

| 10.10 Q | 10.11 Status |
|---------|-----------|
| Q2 dual-track architecture freeze | ✅ executed + B curl extension; pending formal wording |
| Q3 boundary O(0.04) | ✅ explained as support artifact; ghost experiments done |
| Q4 η=0.1 freeze | ✅ written to code; pending ChatGPT confirmation |
| Q5 E_combined | ✅ reused; + J_combined |
| Q6 10.11 priority | ✅ A/B/C/D/E completed |
| Sinusoidal ~0.48 | ⚠️ ghost improved boundary; core not separately reported |

---

## 6. External Messaging (current Cursor understanding, pending ChatGPT confirmation)

**Can say**:
- core region pairwise divA validated for div-free MMS fields (gradDen ~1e-7)  
- boundary divA is support deficiency dominant issue  
- Joule / E_combined / J_combined vs exact ~1e-6–1e-7 (Az2D η=0 MMS)  
- η_A=0.1 primary research candidate  
- B=curlA post-processing should use B-corrected grad  

**Cannot say**:
- Coulomb gauge globally enforced  
- B=curlA validated to machine precision (oscillatory fields ~3.5%)  
- boundary divA solved  
- Contact A-penalty ready / projection implemented  

**Most accurate one-liner**:
> Stage 10.10/10.11 establishes core-region pairwise divA validity and localizes boundary divA to support deficiency; E/J/Joule paths match exact to ~1e-6 while B=curlA on oscillatory fields shows ~3.5% discrete curl gap under B-corrected gradient.

---

## 7. Reproduction Commands

```bash
cd /home/yongchuan/sphinxsys/build
cmake .. -DCMAKE_BUILD_TYPE=Release
ninja test_3d_aphi_ck_pairwise_diva_boundary_diagnostic \
      test_3d_aphi_ck_ghost_buffer_diva_diagnostic \
      test_3d_aphi_ck_curl_b_dual_track_diagnostic \
      test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms

BIN=tests/extra_source_and_tests/3d_examples
for t in pairwise_diva_boundary_diagnostic \
         ghost_buffer_diva_diagnostic \
         curl_b_dual_track_diagnostic \
         inner_divergence_free_a_nonzero_rhs_mms; do
  echo "=== test_3d_aphi_ck_${t} ==="
  "$BIN/test_3d_aphi_ck_${t}/bin/test_3d_aphi_ck_${t}" 2>&1 | grep -E 'passed=|field=|b_corrected|boundary'
done
```

---

## 8. Suggested Prompt to Paste for ChatGPT (copy directly)

```text
You are the technical advisor for A-phi SPHinXsys Stage 10. Cursor completed Stage 10.11-A/B/C/D/E per stage10_11_boundary_divA_B_E_H_projection_eta_plan (see stage10_11_validation_results_and_open_questions_for_chatgpt.md and CURSOR_APHI_STAGE10_11_BOUNDARY_VALIDATION_RECORD.md).

10.11 key conclusions (please read §2):
- core pairwise divA ~1e-7; boundary 100% energy; B-corrected boundary gradDen ~0.037
- Joule/E_combined/J_combined vs exact ~1e-6–1e-7 (Az2D η=0 MMS)
- B=curlA ~3.5% is B-corrected grad + oscillatory exact A discretization limit; pairwise grad curl ~200% (unusable)
- ghost only significantly helps Sinusoidal3DCurlPsi; Az2D/CrossSine lattice sufficient
- η_A=0.1 written to research defaults; projection design only

Based on attachments:
1. Give clear advice on Q1–Q8 line by line (including priority ranking);
2. Confirm or revise dual-track architecture + external messaging (§6);
3. Contact A-penalty still not meeting unlock? What's missing?
4. projection still defer? If done, what are acceptance criteria?
5. Is B ~3.5% blocking cold crucible / impressed case?
6. Output Stage 10.12 execution plan (task names, acceptance criteria, suggested test names).

10.10 resolved: global divA rel≈0.13 root cause (denominator+boundary); non-degenerate MMS.
10.11 still open: Contact unlock, B curl dp refinement, invariance MMS positioning, curlH vs J.
```

---

## 9. Quick Reference Links

| Document | Path |
|------|------|
| **This discussion pack** | [`stage10_11_validation_results_and_open_questions_for_chatgpt.md`](stage10_11_validation_results_and_open_questions_for_chatgpt.md) |
| **10.11 Cursor main record** | [`tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_11_BOUNDARY_VALIDATION_RECORD.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_11_BOUNDARY_VALIDATION_RECORD.md) |
| **10.11 original plan** | [`stage10_11_boundary_divA_B_E_H_projection_eta_plan_for_cursor.md`](stage10_11_boundary_divA_B_E_H_projection_eta_plan_for_cursor.md) |
| **projection design** | [`tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_PROJECTION_DESIGN_NOTES.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_PROJECTION_DESIGN_NOTES.md) |
| **10.10 discussion pack** | [`stage10_10_validation_results_and_open_questions_for_chatgpt.md`](stage10_10_validation_results_and_open_questions_for_chatgpt.md) |
| **10.10 main record** | [`tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md) |
| **Source file map** | [`tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/STAGE10_APHI_FILE_MAP_FOR_DISCUSSION.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/STAGE10_APHI_FILE_MAP_FOR_DISCUSSION.md) |
