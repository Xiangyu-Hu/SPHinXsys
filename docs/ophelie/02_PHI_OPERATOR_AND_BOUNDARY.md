# φ Operator and Boundary Development

> Consolidated from: `OPHELIE_PHI_DEVELOPMENT_LOG.md`, `OPHELIE_PHI_BOUNDARY_AND_NEXT_DEVELOPMENT_PLAN.md`, `OPHELIE_PHI_DISCRETIZATION_NEXT_SPRINT_PLAN.md`, `OPHELIE_PHI_OPERATOR_AUDIT.md`, `OPHELIE_PHI_TEST_UPDATE_AFTER_GRADIENT_RESULT.md`, `OPHELIE_DIVERGENCE_OPERATOR_VALIDATION_AND_EDGE_FLUX_PLAN.md`, `docs/ophelie/PHI_GPT_DISCUSSION_PACKAGE.md`.  
> Archive: [`archive/plans/OPHELIE_PHI_*`](archive/plans/), [`archive/legacy_packages/PHI_GPT_DISCUSSION_PACKAGE.md`](archive/legacy_packages/)

---

## 1. Governing equation (production)

In conductor domain, imag chain drives phasor EM:

```text
E = -∇φ - iωA_active
J = σE
```

Discrete production uses **DivSigmaGrad** LHS and **DivSigmaA** RHS (not legacy flux RHS for French acceptance).

---

## 2. Sprint status

| Sprint | Status | Notes |
|--------|--------|-------|
| Sprint 1 — CLI, RHS fingerprint, linear MMS | ✅ Done | |
| Sprint 2 — discrete self-consistency | ✅ Done | rel ≈ 3.3e-7 |
| Sprint 3 — paired G_c/D_c compatible operator | ✅ MMS / ❌ French Biot | ChatGPT review closed sign issue |
| A — P0 diagnostics (RHS integral, Jn, CSV) | ✅ Done | zero-mean projection inconclusive pre-RHS fix |
| B — P1 one-sided Neumann RHS | ✅ v1 | slab + cylinder MMS |
| B+ — LHS grad Neumann | ❌ Deferred | operator projection ineffective |
| C — P2 corrected gradient alone | ❌ Deferred | worsens phi_eq_res in isolation |
| D/E — virtual shell, A_ind, thermal | ⏸ Deferred | |

---

## 3. Boundary CLI (literature default unchanged)

| Option | Default | Purpose |
|--------|---------|---------|
| `--phi-boundary-mode=none\|one-sided-neumann\|virtual-shell-diagnostic` | `none` | |
| `--phi-boundary-normal-source=analytic-box\|analytic-cylinder\|level-set` | `analytic-cylinder` | |
| `--phi-compatible-correction=0\|1` | `0` | Full G_c/D_c pair |
| `--phi-gradient-correction=0\|1` | `0` | Diagnostic only |
| `--phi-boundary-grad-neumann=0\|1` | `0` | Experimental MMS |

---

## 4. Key numerical findings (2026-06)

### GradPhi sign

Excluded by `test_3d_ophelie_phi_gradient_linear` — do **not** flip `(phi_i - phi_j)`.

### RHS sign audit

`test_3d_ophelie_phi_rhs_flux_sign_audit`: legacy flux RHS ≈ −DivSigmaA on Biot field (cosine = −1, machine ε).

### Solvability reload

`test_3d_ophelie_phi_biot_rhs_solvability` @ dp=0.02, n=16066:

```text
eq_res(div/legacy) ≈ 0.574
divJ_red ≈ 1.74
sign_flip_ok = 1
passed = 1
```

`cross_eq_res ≈ 1.79` is expected when comparing opposite-sign RHS forms, not independent physics error.

---

## 5. Divergence operator validation

MMS and scan outputs archived in `discussion_bundle/ophelie_vector_divergence_mms_scan_v2.csv`.

Plan detail: [`archive/plans/OPHELIE_DIVERGENCE_OPERATOR_VALIDATION_AND_EDGE_FLUX_PLAN.md`](archive/plans/OPHELIE_DIVERGENCE_OPERATOR_VALIDATION_AND_EDGE_FLUX_PLAN.md).

---

## 6. Related tests

```text
test_3d_ophelie_phi_laplace_constant
test_3d_ophelie_phi_gradient_linear
test_3d_ophelie_phi_compatible_operator_mms
test_3d_ophelie_phi_discrete_self_consistency
test_3d_ophelie_phi_neumann_slab
test_3d_ophelie_phi_neumann_cylinder
test_3d_ophelie_phi_biot_rhs_solvability
test_3d_ophelie_phi_rhs_flux_sign_audit
test_3d_ophelie_vector_diverance_mms
```

---

## 7. Next-stage notes (device GMRES, A_ind)

See [`archive/plans/OPHELIE_NEXT_STAGE_DEVICE_GMRES_AIND_PLAN.md`](archive/plans/OPHELIE_NEXT_STAGE_DEVICE_GMRES_AIND_PLAN.md) for device GMRES and self-induction diagnostic roadmap.
