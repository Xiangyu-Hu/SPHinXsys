# Stage 10.13-P5 — Full Maxwell MMS design document (design only, no implementation)

> **Date**: 2026-05-21  
> **Status**: design only — **forbidden** to write production/test implementation code in this stage  
> **Handoff**: [`CURSOR_APHI_STAGE10_13_CHATGPT_HANDOFF.md`](CURSOR_APHI_STAGE10_13_CHATGPT_HANDOFF.md)  
> **Prerequisite**: div-free MMS (10.10–10.11), B/H observable (10.12-P5), B dp refinement (10.13-P3)

---

## 1. What MMS is (and is not)

**MMS = Method of Manufactured Solutions**.

It is not a new full Maxwell solver, nor is it changing the current matrix-free A-phi into a full Maxwell FEM/FD system.

The approach: artificially construct analytic fields `(A_exact, phi_exact, …)`, substitute into **the discrete/continuous equations actually solved by current code**, derive RHS/source, then have the numerical solution recover the manufactured field and verify derived quantities.

```text
MMS is not a new full Maxwell solver.
The purpose is to construct exact A, phi, B, H, E, J fields that are mutually consistent
so that curlH ≈ J can be tested as a hard gate when desired.
```

---

## 2. What current div-free MMS already validates

Existing manufactured fields (`AphiDivFreeValidationFieldKind`: Linear2D, Az2D, CrossSine3D, etc.) mainly guarantee:

| Relation | Current gate |
|------|-----------|
| `div A_exact = 0` (by construction) | core divA ~1e-7 |
| `B = curl A` (post-processing) | B-corrected, Az2D ~0.88% @ dp=0.05 (informational convergence) |
| `H = nu B` | P5 `h_nu_b_curl_h_j`: **hard gate** |
| `E = -iωA - grad phi` (when phi=0, `E = -iωA`) | E/J/Joule MMS ~1e-6–1e-7 |
| `J = sigma E`, Joule | MMS ~1e-6–1e-7 |
| `curl H = J` (Maxwell-Ampere) | **informational only** — div-free A does not guarantee |

**User reproduction (2026-05-21, P3)**: Az2D `B_err`: 3.49% @ dp=0.1 → 1.98% @ 0.075 → **0.88% @ 0.05**; `az2d_stage1013_ok=1`.

Therefore: **E/J/Joule chain can serve as production demo and MMS gate; `curlH≈J` cannot be elevated to hard gate directly from existing div-free fields.**

---

## 3. Full Maxwell MMS target relations

To make `curlH ≈ J` a **hard gate**, manufactured fields should simultaneously satisfy (frequency domain, linear medium):

```text
B_exact  = curl A_exact
H_exact  = nu * B_exact
J_exact  = curl H_exact
E_exact  = J_exact / sigma(x)        (isotropic sigma)
E_exact  = -i*omega*A_exact - grad phi_exact
Joule_exact = 0.5 * sigma * |E_exact|^2
```

Optional (gauge/constraints):

```text
div A_exact = 0          (Coulomb-like gauge, can coexist with penalty research)
div H_exact = 0          (inherited from B when mu is constant)
```

**At discrete verification**: both `curl H` and `J` use the same **B-corrected gradient + curl** post-processing path as P5, avoiding uncorrected pairwise curl false negatives.

---

## 4. Why you cannot "add one term to existing Az2D" to satisfy curlH=J

Current Az2D-like fields: `A` is div-free scalar-potential construction, `J_exact = sigma * (-iωA)` derived from potential equation, **not** constrained by `curl H = J`.

Therefore P5 diagnostic shows:

- `H = nu B` holds (ratio≈1)
- `curlH` and `J` magnitudes reportable, but relative error unstable against zero reference or each other

This is not an implementation bug, but **insufficient manufactured field degrees of freedom**.

---

## 5. Recommended construction strategies (reference for implementation phase)

### 5.1 Strategy A: derive from `J_exact` (recommended for first prototype)

1. Choose divergence-free `J_exact` (e.g. polynomial/sine combination, automatically satisfies `div J = 0` in simple domain).
2. `H_exact = curl Psi` or use known vector potential so `curl H = J` (3D must verify `div J = 0`).
3. `B_exact = H_exact / nu`.
4. `A_exact` such that `curl A = B` (can reuse existing div-free `A` family, but must refit coefficients to satisfy curl relation).
5. `E_exact = J/sigma`, `phi_exact` from `E + iωA = -grad phi` integration/construction (or set `phi=0` and adjust `A` for compatibility).

**Pros**: directly locks `curlH=J`.  
**Cons**: `A` and `E` dual-path consistency (Section 3 two equations) requires algebraic verification.

### 5.2 Strategy B: forward construction from `A_exact, phi_exact` (consistent with current framework)

1. First construct `A, phi` satisfying `div A = 0`.
2. Compute `B, H, E`.
3. **Explicitly compute** `J_amp = curl H`, check if equals `sigma E`; if not, record `source_J` or abandon hard gate.

**Pros**: closest to current MMS helper.  
**Cons**: generally **cannot** simultaneously satisfy `curlH=J` and `E=J/sigma`, only diagnostic.

### 5.3 Strategy C: two-tier gate (pragmatic compromise)

| Tier | Content | gate |
|------|------|------|
| Tier-1 (retain) | div-free + E/J/Joule + H=nuB | existing strict |
| Tier-2 (new) | dedicated `MaxwellAmpereMMS` field family | `curlH vs J` **hard** |

When implementing, create new `AphiDivFreeValidationFieldKind::MaxwellAmpere*` or independent namespace; **do not** break existing Az2D gate.

---

## 6. Recommended test structure (future implementation, not P5)

```text
test_3d_aphi_ck_full_maxwell_mms_manufactured/
  - single body or two-body Inner (no Contact first, reduce interface degrees of freedom)
  - assign exact A, phi
  - post: B, H, E, J, curlH, Joule
  - compare: all relations + split core/boundary
```

**Hard gate recommendations (at implementation)**:

```text
core: |curlH - J| / |J|  <= tol_ampere   (e.g. 1e-2 @ dp=0.1, tighten with dp)
core: |E - J/sigma| / |E|  <= tol_ej
core: |H - nu B| / |H|     <= tol_hb
```

**Informational only** (consistent with current policy):

- boundary `curlH` / `J`
- uncorrected pairwise curl vs corrected curl

---

## 7. Integration points with existing code

| Module | Purpose |
|------|------|
| `aphi_divergence_free_a_mms_helpers.h` | current div-free fields; can add parallel `aphi_full_maxwell_mms_fields.h` |
| `aphi_h_curl_h_j_observable_diagnostic_helpers.h` | `execBodyCurlFromVectorFieldDiagnostic`, Maxwell gate reuse |
| `aphi_curl_b_dp_refinement_diagnostic_helpers.h` | B convergence policy; Maxwell fields also need dp study |
| `test_3d_aphi_ck_h_nu_b_curl_h_j_observable_diagnostic` | template; Maxwell fields should make `curlH_ok` hard |

**Do not** modify `AphiLhsAssemblyOptions` defaults or Contact penalty defaults in this stage.

---

## 8. Explicitly out of P5 / next-phase initial scope

- No projection implementation (`A_new = A - grad chi`)
- Do not impose `curlH≈J` on existing Az2D/Linear2D div-free tests
- No three-body Contact + Maxwell MMS joint gate
- Do not upgrade cold crucible demo to quantitative Maxwell validation
- Do not require B_err → 1e-6 (discrete curl limit see P3)

---

## 9. Prerequisites before starting implementation (recommended)

Before ChatGPT/user explicitly requests **Stage 10.14 implementation**, should have:

1. Stage 10.13 **P0–P4** all `passed=1` (including user local P3/P4 reproduction).
2. This design document reviewed (choose one of construction strategies A/B/C).
3. At least first manufactured field decided: **LinearMaxwell** (polynomial) or **AzMaxwell** (sinusoidal, compare with Az2D).

---

## 10. External messaging template

```text
We use Method of Manufactured Solutions (MMS) to verify the A-phi matrix-free pipeline.
The existing div-free MMS validates A/phi, E, J, and Joule; H=nuB is verified separately.
A future Maxwell-Ampere-consistent MMS will be added to validate curlH≈J as a hard gate.
This is not a full Maxwell eigensolver—it is a manufactured-field consistency test only.
```

English summary:

```text
MMS is used to verify the A-phi matrix-free pipeline, not a full Maxwell solver.
Current div-free MMS already covers E/J/Joule; curlH≈J will be a hard gate in a self-consistent Maxwell-Ampere MMS field family.
```

---

## 11. Record index

| Document | Content |
|------|------|
| this document | P5 design |
| [`CURSOR_APHI_STAGE10_12_B_H_OBSERVABLE_RECORD.md`](CURSOR_APHI_STAGE10_12_B_H_OBSERVABLE_RECORD.md) | P5 curlH/J informational |
| [`CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md`](CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md) | div-free MMS |
| [`CURSOR_STAGE10_13_APHI_NEXT_WORK_DETAILED.md`](../../../../../../docs/electromagnetic_aphi/CURSOR_STAGE10_13_APHI_NEXT_WORK_DETAILED.md) | 10.13 master plan |
