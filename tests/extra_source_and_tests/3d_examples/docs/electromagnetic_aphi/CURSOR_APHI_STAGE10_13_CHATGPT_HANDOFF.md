# Stage 10.13 — ChatGPT / Cursor Handoff (Frozen State)

> **Date**: 2026-05-21  
> **Branch**: `feature/electromagnetic`  
> **Detailed plan**: [`CURSOR_STAGE10_13_APHI_NEXT_WORK_DETAILED.md`](../../../../../../docs/electromagnetic_aphi/CURSOR_STAGE10_13_APHI_NEXT_WORK_DETAILED.md) (see docs/electromagnetic_aphi/)  
> **Prior records**: 10.11 boundary, 10.12 Contact A-penalty, 10.12 B/H observable

---

## Frozen Conclusions (Auto / Subsequent Cursor Must Follow)

```text
production lambda_A = off
research eta_A = 0.1
Contact A-penalty stencil = InnerOnly (default; main K(X) still Inner + Contact)
eta_A >= 0.2 => polish_sweeps = 0 (block-GS polish destroys solution)
projection = deferred (no post-solve A/phi projection implementation)
boundary divA = open diagnostic / support deficiency (not core operator bug)
B=curlA = corrected LinearGradient + LinearCorrectionMatrix route as diagnostic/reference
uncorrected pairwise curl = comparison only, not a hard gate
E/J/Joule = MMS validated, usable for demo post-processing pipeline
curlH ~= J = informational only (div-free MMS does not guarantee Maxwell-Ampere; see P5 design)
cold crucible demo = pipeline smoke test, not quantitative validation
three-body Contact A-penalty = deferred
```

---

## Contact Two-Body A-Penalty (Research Gate Status)

```text
Contact two-body A-penalty research gate closed for eta=0 and eta=0.1 with InnerOnly stencil;
eta=0.2 is usable with polish_sweeps=0 (global_true_rel ~ O(1e-5));
three-body A-penalty remains deferred.
```

**Do not** write "Contact η>0 open" or "η>0 GMRES does not converge" as current blockers — C4 InnerOnly + C5 polish strategy is closed.

C1/C2 diagnostics still **explicitly** use `InnerContact` penalty pipeline to test legacy path; **default** production/research assembly is `InnerOnly`.

---

## Solve vs Post-Processing (External Wording)

```text
The matrix-free A-phi solve still uses the pairwise Laplace-form operator. The magnetic field B=curlA is an observable reconstructed from the solved vector potential. For this first-derivative reconstruction, the corrected SPHinXsys LinearGradient route is used as the diagnostic/reference path. The uncorrected pairwise curl is retained only as a diagnostic comparison and must not be used as a hard accuracy gate.
```

---

## Contact Test Semantics

- **Two independent `SolidBody` instances**, finding each other's particles via `Contact<>`; **not** two regions within one body.
- Typical: `left_body` + `right_body`, each with `Inner<>` + `Contact<>` pointing to the other.

---

## Stage 10.13 Work Order (Current)

| Priority | Content | Status |
|--------|------|------|
| P0 | this document + old record revision + C4/C5 gate hardening | ✅ |
| P1 | `test_3d_aphi_ck_cold_crucible_demo` (MVP-A impressed A) | ✅ |
| P2 | demo VTP + magnitude scalars | ✅ |
| P3 | B dp=0.05 (Az2D gate) | ✅ |
| P4 | InnerOnly graddiv PC supplemental test | ✅ |
| P5 | [`CURSOR_APHI_STAGE10_13_FULL_MAXWELL_MMS_DESIGN.md`](CURSOR_APHI_STAGE10_13_FULL_MAXWELL_MMS_DESIGN.md) | ✅ design only |
| P6–P7 | three-body A-penalty, projection | **forbidden** |

---

## Key Tests and Expectations

| Test | Expectation |
|------|------|
| `..._two_body_mms` | `passed=1`, primary η≤0.1 strict |
| `..._inner_only_stencil_diagnostic` | `passed=1`, includes `research_improved` |
| `..._upper_eta_solver_diagnostic` | `passed=1`, η=0.2 polish=0 true_rel ≤ 1e-4 |
| `..._h_nu_b_curl_h_j_observable_diagnostic` | `passed=1`, curlH/J not a hard gate |

Regression list see Stage 10.13 plan §8.

---

## Record Index

| File | Content |
|------|------|
| [`CURSOR_APHI_STAGE10_11_BOUNDARY_VALIDATION_RECORD.md`](CURSOR_APHI_STAGE10_11_BOUNDARY_VALIDATION_RECORD.md) | boundary divA, B/E/J, ghost |
| [`CURSOR_APHI_STAGE10_12_CONTACT_A_PENALTY_RECORD.md`](CURSOR_APHI_STAGE10_12_CONTACT_A_PENALTY_RECORD.md) | C1–C5, C3+, InnerOnly |
| [`CURSOR_APHI_STAGE10_12_B_H_OBSERVABLE_RECORD.md`](CURSOR_APHI_STAGE10_12_B_H_OBSERVABLE_RECORD.md) | P5 H=nuB |
| [`CURSOR_APHI_STAGE10_12_ROADMAP.md`](CURSOR_APHI_STAGE10_12_ROADMAP.md) | 10.12 completed items |
