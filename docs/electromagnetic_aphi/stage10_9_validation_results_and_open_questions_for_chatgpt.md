# Stage 10.9 Validation Results + Open Questions (ChatGPT Discussion Pack)

> **Date**: 2026-05-21  
> **Status**: Stage 10.9-A/B/C/D implemented and run on Cursor side; **recommend aligning with ChatGPT before writing next batch of code**  
> **Prior plan**: [`stage10_8_metric_interpretation_and_validation_plan_for_cursor.md`](stage10_8_metric_interpretation_and_validation_plan_for_cursor.md)  
> **Cursor full record**: [`tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_8_PAIRWISE_DIVA_RECORD.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_8_PAIRWISE_DIVA_RECORD.md) §14

---

## 0. One sentence for ChatGPT

Stage 10.8 pairwise divA unification + Stage 10.9 validation (pairwise diagnostic, MMS, dp sweep, Case C observable reference) **closed on Inner research stack**. Please help decide: **next priority is div-free MMS (b≠0), 3D pairwise boundary fix, η_A default freeze, or Contact A-penalty prototype?** And whether current results suffice to claim externally "research prototype convergent."

---

## 1. Completed vs Original Plan Task Mapping

| Original Task | Cursor Status | Test / Record |
|-------------|------------|-------------|
| Task 1 metric interpretation | ✅ §14.1 written to record | see record §14.1 |
| Task 2 divergence-free MMS | ⚠️ partial: 10.9-A/B; **exact div-free MMS degenerate** | see §3 open Q1 |
| Task 3 dp refinement | ✅ 10.9-C | `test_3d_aphi_ck_inner_divergence_free_a_dp_refinement_diagnostic` |
| Task 4 observable gate split | ✅ 10.9-D Case C | `test_3d_aphi_ck_inner_a_divergence_penalty_observable_reference_diagnostic` |
| Task 5 η_A default candidates | ⏸ code still 0.1–0.3 research defaults | discuss Q4 |
| Task 6 Contact deferred | ✅ still deferred | — |

---

## 2. Key Numerical Summary (for discussion)

### 2.1 Inner gate (Stage 10.8, solenoidal, dp=0.1)

| η_A | true_rel | div_A_rel | Joule Δ vs η=0 | E_L2 Δ vs η=0 |
|-----|----------|-----------|----------------|---------------|
| 0 | 1.0e-5 | 0.54 | — | — |
| 0.1 | 2.7e-5 | 0.20 | 13% | 6.5% |
| 0.3 | 5.2e-6 | 0.10 | 14% | 7.0% |

### 2.2 Pairwise divA on exact curl field (10.9-A)

| Field | B-corrected div_A_rel | Pairwise div_A_rel |
|----|----------------------|-------------------|
| Linear2D `A=(-y/2,x/2,0)` | ~1.1e-8 | **~2.4e-9** |
| Sinusoidal3DCurlPsi | ~0.011 | **~0.48** (informational) |

### 2.3 MMS manufactured separable (10.9-B, non div-free)

- consistency η=0,0.1,0.2,0.3: **4/4 passed** (`block_Linf_err≤0.15`, `true_rel≤1e-4`)
- invariance η=0.1,0.2,0.3: **3/3 passed** (GMRES convergence gate only)
- div-free exact field: `rhs_norm=0` → MMS degenerate

### 2.4 dp refinement (10.9-C, Linear2D pairwise + MMS η=0.1)

| dp | linear2d pairwise div_A_rel | mms block_Linf |
|----|----------------------------|----------------|
| 0.15 | 3.2e-9 | 2.6e-6 |
| 0.10 | 2.4e-9 | 1.9e-6 |
| 0.075 | 4.5e-9 | 3.3e-5 |

Linear2D pairwise resolution stable; MMS block **did not** monotonically decrease coarse→fine (informational).

### 2.5 Case C observable reference (10.9-D)

**Solenoidal** (dp=0.1, vs same dp η=0 baseline):

| η_A | div_A_rel | Joule err vs η=0 | E_L2 err vs η=0 |
|-----|-----------|------------------|-----------------|
| 0.1 | 0.20 | 13% | 6.5% |
| 0.3 | 0.10 | 14% | 7.0% |

**Mesh gap** (informational): dp=0.1 η=0 vs dp=0.075 η=0 → Joule **93%**, E **46%**.

**MMS** (vs η=0 consistency ref): observable error **<0.05%**; div_A_rel≈3.7 unchanged (field non div-free).

---

## 3. Open Questions (please advise line by line)

### Q1. div-free exact MMS degeneracy: `K(A_exact,φ)≈0`

- Linear2D `A=(-y/2,x/2,0)`, Sinusoidal3DCurlPsi, adding nonzero φ all show `rhs_norm≈0`.
- **Question**: How to construct manufactured field so (1) analytic divA=0, (2) discrete `b=K(A_exact)≠0`, (3) analytic E/J/Joule computable?
- **Alternative**: Use **source-driven + high-resolution η=0 self-reference** instead of div-free MMS for absolute validation?

### Q2. 3D pairwise divA on Sinusoidal3DCurlPsi rel≈0.48

- B-corrected still ~0.01; Linear2D ~1e-9.
- **Question**: Must fix pairwise discretization (boundary/correction terms) to claim gauge diagnostic valid? Or 2D gate + B-corrected 3D sufficient, pairwise 3D marked informational?

### Q3. Case C conclusion boundaries

- Solenoidal: div_A down, observable vs **same dp η=0** perturbation 7–14%, but **no analytic exact**.
- Fine-dp self-ref vs cand-dp gap 50–90% (mesh gap).
- **Question**: Is this enough to say "η=0.1 research candidate"? What reference needed to say "E/J/Joule closer to true solution"?

### Q4. η_A default value

- Current `AphiADivergencePenaltyResearchDefaults`: η∈[0.1, 0.3].
- Case C + gate: η=0.1 perturbation smaller than η=0.2 (Joule 22%).
- **Question**: Freeze primary=**0.1**, optional=**0.2**, 0.3 upper research bound only? Should code constants change?

### Q5. Next implementation priority (pick one or rank)

1. div-free MMS + analytic E/J/Joule reference  
2. 3D pairwise divA boundary fix  
3. Contact A-penalty full stack (Inner+Contact divA/PC)  
4. impressed case informational gate + documentation/external messaging freeze  

### Q6. External messaging (Part V)

Current Cursor understanding:

- **Can say**: Inner pairwise A-penalty research prototype convergent; divA reducible by several times; observable perturbation research acceptable.  
- **Cannot say**: Coulomb gauge enforced; E/J/Joule validated; production ready.  

**Question**: After Stage 10.9 completion, do these boundaries need revision?

---

## 4. Recommended Upload File List (by priority)

### 4.1 Must read (discussion decisions)

| File | Description |
|------|------|
| [`stage10_9_validation_results_and_open_questions_for_chatgpt.md`](stage10_9_validation_results_and_open_questions_for_chatgpt.md) | This document |
| [`tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_8_PAIRWISE_DIVA_RECORD.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_8_PAIRWISE_DIVA_RECORD.md) | **Cursor main record** §13–§15 |
| [`stage10_8_metric_interpretation_and_validation_plan_for_cursor.md`](stage10_8_metric_interpretation_and_validation_plan_for_cursor.md) | Original validation plan |

### 4.2 Plan context (optional)

| File | Description |
|------|------|
| [`stage10_mainline_pc_lambda_divA_next_plan_for_cursor.md`](stage10_mainline_pc_lambda_divA_next_plan_for_cursor.md) | Stage 10.8 main plan |
| [`tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_7_A_DIVERGENCE_PENALTY_RECORD.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_7_A_DIVERGENCE_PENALTY_RECORD.md) | Stage 10.7 closure background |
| [`tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/STAGE10_APHI_FILE_MAP_FOR_DISCUSSION.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/STAGE10_APHI_FILE_MAP_FOR_DISCUSSION.md) | Source file map |

### 4.3 New code (10.9 implementation)

| File | Description |
|------|------|
| [`tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_divergence_free_a_mms_helpers.h`](tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_divergence_free_a_mms_helpers.h) | MMS + dp sweep helper |
| [`tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_observable_reference_benchmark_helpers.h`](tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_observable_reference_benchmark_helpers.h) | Case C helper |

### 4.4 New tests (run passed=1)

| Test directory | Phase |
|----------|-------|
| [`tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_inner_divergence_free_a_pairwise_diagnostic/`](tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_inner_divergence_free_a_pairwise_diagnostic/) | 10.9-A |
| [`tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_inner_divergence_free_a_mms/`](tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_inner_divergence_free_a_mms/) | 10.9-B |
| [`tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_inner_divergence_free_a_dp_refinement_diagnostic/`](tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_inner_divergence_free_a_dp_refinement_diagnostic/) | 10.9-C |
| [`tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_inner_a_divergence_penalty_observable_reference_diagnostic/`](tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_inner_a_divergence_penalty_observable_reference_diagnostic/) | 10.9-D |

### 4.5 Existing gate (reference)

| Test | Description |
|------|------|
| [`tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_inner_a_divergence_penalty_gate_diagnostic/`](tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_inner_a_divergence_penalty_gate_diagnostic/) | Stage 10.8 Inner gate |
| [`tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_inner_a_divergence_penalty_gate_helpers.h`](tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_inner_a_divergence_penalty_gate_helpers.h) | gate helper |

---

## 5. Reproduction Commands

```bash
cd /home/yongchuan/sphinxsys/build
cmake .. -DCMAKE_BUILD_TYPE=Release
ninja test_3d_aphi_ck_inner_a_divergence_penalty_gate_diagnostic \
      test_3d_aphi_ck_inner_divergence_free_a_pairwise_diagnostic \
      test_3d_aphi_ck_inner_divergence_free_a_mms \
      test_3d_aphi_ck_inner_divergence_free_a_dp_refinement_diagnostic \
      test_3d_aphi_ck_inner_a_divergence_penalty_observable_reference_diagnostic

BIN=tests/extra_source_and_tests/3d_examples
for t in inner_a_divergence_penalty_gate_diagnostic \
         inner_divergence_free_a_pairwise_diagnostic \
         inner_divergence_free_a_mms \
         inner_divergence_free_a_dp_refinement_diagnostic \
         inner_a_divergence_penalty_observable_reference_diagnostic; do
  echo "=== test_3d_aphi_ck_${t} ==="
  "$BIN/test_3d_aphi_ck_${t}/bin/test_3d_aphi_ck_${t}" 2>&1 | grep -E 'passed=|row_passed='
done
```

---

## 6. Suggested Prompt to Paste for ChatGPT (copy directly)

```text
You are the technical advisor for A-phi SPHinXsys Stage 10. Cursor completed Stage 10.9-A/B/C/D per stage10_8_metric_interpretation_and_validation_plan (see stage10_9_validation_results_and_open_questions_for_chatgpt.md and CURSOR_APHI_STAGE10_8_PAIRWISE_DIVA_RECORD.md §14).

Based on attachments:
1. Give clear advice on Q1–Q6 line by line (including priority ranking);
2. Output Stage 10.10 execution plan (if needed): task names, acceptance criteria, suggested test names;
3. Update Part V external messaging (can say / cannot say);
4. Agree or not: "Do not touch Contact A-penalty until div-free MMS or equivalent absolute validation closes."

Key blockers:
- div-free exact MMS: K(A_exact)≈0, rhs_norm=0;
- 3D Sinusoidal3DCurlPsi pairwise div_A_rel≈0.48;
- Case C only has same dp η=0 baseline and fine-dp mesh gap, no analytic E/J/Joule.
```

---

## 7. Cursor-side Suggestions (for user reference, not ChatGPT input)

| Option | When to choose |
|------|--------|
| **Discuss first** ✅ recommended | 10.9 plan items done; next fork has large impact; avoid div-free MMS dead ends |
| Continue coding | After ChatGPT clarifies Q1/Q5, then tackle div-free MMS or 3D pairwise fix |

First batch of code Cursor can execute after discussion (depends on ChatGPT answers):

- **If Q1 has construction formula** → new manufactured div-free field + MMS + analytic E/J/Joule test  
- **If Q2 requires 3D pairwise fix** → boundary correction or use B-corrected as gate  
- **If Q4 freezes η=0.1** → change `AphiADivergencePenaltyResearchDefaults` + update gate docs  

---

## 8. Quick Reference Links

- Main record: [`CURSOR_APHI_STAGE10_8_PAIRWISE_DIVA_RECORD.md`](tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_8_PAIRWISE_DIVA_RECORD.md)
- Original validation plan: [`stage10_8_metric_interpretation_and_validation_plan_for_cursor.md`](stage10_8_metric_interpretation_and_validation_plan_for_cursor.md)
- Stage 10.8 main plan: [`stage10_mainline_pc_lambda_divA_next_plan_for_cursor.md`](stage10_mainline_pc_lambda_divA_next_plan_for_cursor.md)
- Discussion pack (this doc): [`stage10_9_validation_results_and_open_questions_for_chatgpt.md`](stage10_9_validation_results_and_open_questions_for_chatgpt.md)
