# Stage 10.12 Progress Roadmap

> **Last updated**: 2026-05-21 (Stage 10.13 P0)  
> **Handoff**: [`CURSOR_APHI_STAGE10_13_CHATGPT_HANDOFF.md`](CURSOR_APHI_STAGE10_13_CHATGPT_HANDOFF.md)  
> **Main record**: [`CURSOR_APHI_STAGE10_12_CONTACT_A_PENALTY_RECORD.md`](CURSOR_APHI_STAGE10_12_CONTACT_A_PENALTY_RECORD.md)

---

## Completed

| ID | Content | Test | Status |
|----|------|------|------|
| 10.12-A | B=curlA dp refinement | `test_3d_aphi_ck_curl_b_dp_refinement_diagnostic` | ✅ |
| 10.12-C1 | Contact pairwise apply equivalence | `..._pairwise_operator_diagnostic` | ✅ |
| 10.12-C2 | Contact PC consistency | `..._graddiv_pc_consistency_diagnostic` | ✅ |
| 10.12-B+ | Interface MMS spatial decomposition | `..._interface_mms_consistency_diagnostic` | ✅ |
| 10.12-B++ | base vs penalty decomposition | `..._interface_penalty_breakdown_diagnostic` | ✅ |
| 10.12-P4 | Inner η observable summary table | `test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms` | ✅ |
| 10.12-C3+ | two-body MMS η∈{0,0.1,0.2} (InnerOnly + polish fix) | `..._two_body_mms` | ✅ |
| 10.12-C4 | Inner-only penalty stencil comparison | `..._inner_only_stencil_diagnostic` | ✅ |
| 10.12-C5 | upper-η solver / polish root cause | `..._upper_eta_solver_diagnostic` | ✅ |
| 10.12-P5 | H=nuB + curlH vs J | `..._h_nu_b_curl_h_j_observable_diagnostic` | ✅ |

---

## Open

| ID | Issue | Known | Next Step |
|----|------|------|--------|
| boundary divA ~27% | support deficiency | confirmed in 10.11 | open diagnostic |
| B ~3.5% @ dp=0.1 | discrete curl | converges with dp | informational |

---

## Stage 10.13 (Current)

| ID | Content | Status |
|----|------|------|
| P0 | handoff + record revision + C4/C5 gate hardening | ✅ |
| P1 | `test_3d_aphi_ck_cold_crucible_demo` | ✅ |
| P2 | demo VTP / magnitude scalars | ✅ |
| P4 | `test_3d_aphi_ck_contact_inner_only_graddiv_pc_consistency_diagnostic` | ✅ |
| P3 | `test_3d_aphi_ck_curl_b_dp_refinement_diagnostic` + dp=0.05 | ✅ |

## Deferred

- Three-body TEAM7 + A-penalty  
- production default λ_A  
- projection implementation  
- 10D thermal (cold crucible **demo** in P1, not quantitative validation)

---

## This Round

Stage 10.12 gates **all complete**. Stage 10.13 **P0–P5** complete (P5 design doc only). Next: start Maxwell-Ampere MMS **implementation** on demand (Stage 10.14, not scheduled).

---

## History (C3+ Primary Gate — Passed)

```text
primary_rows_passed=2/2 upper_eta_reported=1/1 passed=1
```

---

## Historical Reproduction (C4, Completed)

```bash
cd ~/sphinxsys/build
ninja test_3d_aphi_ck_contact_a_divergence_penalty_inner_only_stencil_diagnostic
./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_contact_a_divergence_penalty_inner_only_stencil_diagnostic/bin/test_3d_aphi_ck_contact_a_divergence_penalty_inner_only_stencil_diagnostic
```

C4 result: `research_improved=1` `passed=1`; InnerOnly η=0.1 GMRES converged, true_rel ~3e-6.

---

## Historical Reproduction (P4 / B++, Completed)

```bash
cd ~/sphinxsys/build
cmake .. -DCMAKE_BUILD_TYPE=Release

# P4: Inner η summary table (existing test, new summary rows)
ninja test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms
./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms/bin/test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms

# B++: interface base vs penalty decomposition (new)
ninja test_3d_aphi_ck_contact_a_divergence_penalty_interface_penalty_breakdown_diagnostic
./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_contact_a_divergence_penalty_interface_penalty_breakdown_diagnostic/bin/test_3d_aphi_ck_contact_a_divergence_penalty_interface_penalty_breakdown_diagnostic
```

Watch output:
- P4: `eta_observable_matrix` block (Az2D/CrossSine × η × consistency/invariance)
- B++: `interface_base_rel` vs `interface_penalty_only_rel`

---

## Decision Record

1. **Do P4 + B++ first** (independent diagnostics, do not block main line) — code complete  
2. **Contact η>0 fix** wait for B++ data before engineering (no architecture discussion first)  
3. **ChatGPT** only when a decision is needed on "whether to change split geometry" — **not needed now**
