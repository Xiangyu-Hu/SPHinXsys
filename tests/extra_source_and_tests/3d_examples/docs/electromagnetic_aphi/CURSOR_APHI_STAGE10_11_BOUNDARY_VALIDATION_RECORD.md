# Stage 10.11 — boundary divA, B/E/J/H exact validation record

> **Date**: 2026-05-21  
> **Status**: 10.11-A/B/C/D/E + 10.12-A/B (Contact C1–C5 + C3+) complete; **Stage 10.13 handoff** see [`CURSOR_APHI_STAGE10_13_CHATGPT_HANDOFF.md`](CURSOR_APHI_STAGE10_13_CHATGPT_HANDOFF.md)  
> **Plan**: [`stage10_11_boundary_divA_B_E_H_projection_eta_plan_for_cursor.md`](../../../stage10_11_boundary_divA_B_E_H_projection_eta_plan_for_cursor.md)  
> **Prerequisite**: [`CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md`](CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md)

---

## 0. Alignment summary with ChatGPT 10.11 plan

| ChatGPT conclusion | Cursor opinion | Action |
|-------------|------------|------|
| boundary divA is support artifact, not global A error | ✅ consistent with 10.10 | 10.11-A boundary focus |
| Dual track: pairwise penalty + B-corrected diagnostic | ✅ adopted | metrics/report retained |
| Must validate B/E/J/Joule/H | ✅ | 10.11-B extended MMS |
| projection not implemented yet | ✅ | design document only |
| η_A=0.1 primary research default | ✅ | written to `AphiADivergencePenaltyResearchDefaults` |
| ghost/buffer boundary treatment | ✅ 10.11-C | see §6 |

---

## 1. New code

| File | Purpose |
|------|---------|
| [`diagnostics/aphi_pairwise_diva_boundary_helpers.h`](../extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_pairwise_diva_boundary_helpers.h) | boundary/core gradDen decomposition + boundary_width sweep |
| [`diagnostics/aphi_ghost_buffer_diva_diagnostic_helpers.h`](../extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_ghost_buffer_diva_diagnostic_helpers.h) | supplemental exterior ghost divA |
| [`test_3d_aphi_ck_ghost_buffer_diva_diagnostic/`](test_3d_aphi_ck_ghost_buffer_diva_diagnostic/) | Stage 10.11-C |
| MMS helper extension | B/H/J_combined vs exact; divA core/boundary gradDen |
| [`CURSOR_APHI_STAGE10_PROJECTION_DESIGN_NOTES.md`](CURSOR_APHI_STAGE10_PROJECTION_DESIGN_NOTES.md) | Stage 10.11-D design only |
| `aphi_coupling_modes_ck.h` | `primary_eta_a=0.1` frozen |
| [`test_3d_aphi_ck_curl_b_dual_track_diagnostic/`](test_3d_aphi_ck_curl_b_dual_track_diagnostic/) | Stage 10.11-E B curl dual track |
| [`test_3d_aphi_ck_curl_b_dp_refinement_diagnostic/`](test_3d_aphi_ck_curl_b_dp_refinement_diagnostic/) | Stage 10.12-A B dp sweep |
| [`diagnostics/aphi_curl_b_dp_refinement_diagnostic_helpers.h`](../extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_curl_b_dp_refinement_diagnostic_helpers.h) | B/E/divA dp refinement |

---

## 2. Stage 10.11-A: boundary divA focus

Test: `test_3d_aphi_ck_pairwise_diva_boundary_diagnostic` → **passed=1**

Az2D @ dp=0.1, boundary_width=3dp, core_shell=2.5dp:

| Metric | core | boundary |
|------|------|----------|
| pairwise gradDen rel | **~1e-7** | **~0.14** |
| B-corrected gradDen rel | **~1e-7** | **~0.037** |
| divA energy fraction | — | **100%** |

**Conclusion**: Consistent with 10.10-B2; boundary_width sweep (2/3/4 dp) does not change core metrics.

---

## 3. Stage 10.11-B: B/E/J/H exact observable

Test: `test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms` → **passed=1**

Az2D η=0 consistency (MMS recovers exact field):

| Metric | vs exact | Notes |
|------|----------|------|
| Joule | ~1e-6 | ✅ |
| E_combined | ~5e-7 | ✅ |
| J_combined | ~5e-7 | ✅ |
| **B=curlA** | **~0.035** | ⚠️ discrete curl vs analytic B ~3.5% |
| H (=νB) | ~0.035 | constant ν, same as B |
| divA core gradDen | ~5e-7 | ✅ |
| divA boundary gradDen | ~0.037 | ⚠️ boundary open |

**New finding**: Even with exact A field, **SPH curl(grad A) vs Az2D/CrossSine3D analytic B still has ~3.5% core error** (oscillatory field gradient discretization error). This is not MMS solve failure, but **B post-processing discretization limit**.

invariance (η>0, b_0 mode): block/B/Joule errors still large — boundary penalty changes solution, informational only.

---

## 4. η_A policy (formally frozen research defaults)

```text
primary_eta_a   = 0.1
optional_eta_a  = 0.2
upper_eta_a     = 0.3
production      = lambda_A off
```

---

## 5. External messaging (10.11 update)

**Can say**:
- core region pairwise divA validated for div-free MMS field (gradDen ~1e-7)
- boundary divA is support deficiency dominant issue
- Joule / E_combined / J_combined vs exact ~1e-6–1e-7 (Az2D η=0 MMS)
- η_A=0.1 primary research candidate

**Cannot say**:
- Coulomb gauge globally enforced
- B=curlA validated to machine precision (oscillatory field ~3.5% discrete gap)
- boundary divA solved
- **Contact A-penalty production-ready** (two-body research: η=0/0.1 strict ✅; η=0.2 requires `polish_sweeps=0`; see 10.12-B + 10.13 handoff)
- projection implemented

---

## 6. Stage 10.11-C: ghost/buffer analytic layer

Test: `test_3d_aphi_ck_ghost_buffer_diva_diagnostic` → **passed=1**

Method: device pairwise divA baseline + host-side **supplemental** exterior analytic ghost points (does not duplicate existing lattice particles).

| Field | baseline boundary gradDen | + supplemental ghost | Conclusion |
|-------|--------------------------|----------------------|------|
| Az2D | ~0.031 | ~0.52 (worse) | lattice buffer sufficient; extra ghost over-corrects |
| CrossSine3D | ~0.031 | ~0.52 (worse) | same as above |
| Sinusoidal3DCurlPsi | ~0.098 | **~0.045** | ghost reduces ~54%, support deficiency confirmed |

**Conclusion**:
1. Az2D/CrossSine under existing `boundary_width=3dp` lattice already have boundary gradDen ~0.03.  
2. Sinusoidal3DCurlPsi still lacks exterior support; analytic ghost significantly improves.  
3. Should not blindly add ghost for all fields; decide per field/boundary diagnostic.

---

## 7. Stage 10.11-E: B curl dual-track diagnostic

Test: `test_3d_aphi_ck_curl_b_dual_track_diagnostic` → **passed=1**

Compare **B-corrected LinearGradient** vs **pairwise uncorrected grad** then curl:

| Field | B-corrected core rel | Pairwise uncorrected core rel | Conclusion |
|-------|---------------------|------------------------------|------|
| Linear2D | **~1e-7** | ~2.0 | B-corrected correct; pairwise unsuitable for curl B |
| Az2D | **~0.035** | ~2.0 | ~3.5% is B-corrected oscillatory field discretization limit |
| CrossSine3D | **~0.035** | ~2.0 | same as above |

**Conclusion**:
1. B=curlA post-processing should use **B-corrected grad** (consistent with divA diagnostic dual track).  
2. Az2D/CrossSine ~3.5% **is not** a pairwise vs B-corrected choice issue; pairwise grad is far worse for curl.  
3. ~3.5% is **oscillatory exact A field + B-corrected curl discretization error**, not MMS solve failure (Joule/E ~1e-6 already validated).

---

## 8. Next steps

1. ~~**ghost/buffer**~~ ✅ 10.11-C  
2. ~~**B curl dual track**~~ ✅ 10.11-E  
3. ~~**B dp refinement**~~ ✅ 10.12-A (see §9)  
4. ~~**Contact A-penalty C1–C3**~~ ✅ 10.12-B (see [`CURSOR_APHI_STAGE10_12_CONTACT_A_PENALTY_RECORD.md`](CURSOR_APHI_STAGE10_12_CONTACT_A_PENALTY_RECORD.md))  
5. **Contact η>0 interface MMS spatial decomposition** — ✅ 10.12-B+  
6. ~~**Inner η observable summary table**~~ — ✅ 10.12-P4  
7. **Contact interface penalty engineering fix** — next step (see roadmap)

---

## 9. Stage 10.12-A: B=curlA dp refinement

Test: `test_3d_aphi_ck_curl_b_dp_refinement_diagnostic` → **passed=1**

Method: exact div-free A field + B-corrected grad curl; dp = 0.15 / 0.10 / 0.075.

| Field | dp=0.15 | dp=0.10 | dp=0.075 | Trend |
|-------|---------|---------|----------|------|
| Linear2D B rel | ~2e-7 | ~1.6e-7 | ~1.4e-7 | ✅ stable ~1e-7 |
| Az2D B rel | **~0.077** | **~0.035** | **~0.020** | ✅ decreases with dp |
| CrossSine3D B rel | ~2e11 (under-resolved) | **~0.035** | **~0.020** | ✅ dp≥0.10 same as Az2D |

E_combined / Joule vs exact: all dp **~1e-7 or 0** (φ=0 exact field).

**Conclusion**:
1. Az2D/CrossSine (dp≥0.10) ~3.5% **converges with dp** (0.10→0.075: 0.035→0.020); **discrete curl error**, not implementation bug.  
2. CrossSine3D @ dp=0.15 (N=343) **under-resolved**, B metric unusable; 3D oscillatory field needs dp≤0.10.  
3. **Not blocking** E/J/Joule mainline; claim "B validated" must note **dp dependence** (dp=0.075 still ~2%).  
4. boundary divA gradDen slightly decreases with dp (~0.037→~0.027), still open.

---

## 10. Stage 10.12-B: Contact A-penalty (summary)

See [`CURSOR_APHI_STAGE10_12_CONTACT_A_PENALTY_RECORD.md`](CURSOR_APHI_STAGE10_12_CONTACT_A_PENALTY_RECORD.md).

| Gate | Result |
|------|------|
| C1 apply equivalence | passed=1 (stencil-safe) |
| C2 PC consistency | passed=1 (core) |
| C3 MMS η=0 | passed=1 |
| C3+ MMS η=0.1/0.2 | **closed** (InnerOnly stencil; η=0.2 + `polish_sweeps=0`) — see 10.12 Contact record |

---

## 11. Reproduction

```bash
cd build && ninja test_3d_aphi_ck_curl_b_dp_refinement_diagnostic \
              test_3d_aphi_ck_ghost_buffer_diva_diagnostic \
              test_3d_aphi_ck_curl_b_dual_track_diagnostic \
              test_3d_aphi_ck_pairwise_diva_boundary_diagnostic \
              test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms \
              test_3d_aphi_ck_contact_a_divergence_penalty_two_body_mms
```
