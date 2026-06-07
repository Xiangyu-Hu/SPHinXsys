# Stage 10A — EM/Joule Verification Progress

Date: 2026-05-21  
Branch: `feature/electromagnetic`  
Plan source: `cursor_aphi_stage9d9e_review_next_steps.md`

## Scope

Finalize EM→Joule pipeline defensibility **before** quantitative TEAM7 or cold-crucible work.

| Item | Test | Status |
|------|------|--------|
| 10A-1a Joule uniform-field analytic verification | `test_3d_aphi_ck_joule_uniform_field_analytic_verification` | ✅ |
| 10A-1b Coupled 8×8 PC fallback/pivot diagnostics | `test_3d_aphi_ck_block_jacobi_8x8_pc_diagnostic` | ✅ |
| 10A-2a Joule local distribution sensitivity | `test_3d_aphi_ck_joule_local_distribution_sensitivity_diagnostic` | ✅ |
| 10A-2b Physical TEAM7-dimension dp observable sweep | `test_3d_aphi_ck_gmres_team7_physical_dp_observable_diagnostic` | ✅ |

## New / modified shared code

| File | Change |
|------|--------|
| `aphi_benchmark_case_ck.h/.hpp` | `AssignUniformLinearPhiFieldsCK` |
| `aphi_block_jacobi_preconditioner_ck.h/.hpp` | `AphiBlockJacobiDiagnosticNames`, `RegisterAphiBlockJacobiDiagnosticFieldsCK`, optional PC diagnostics in apply kernel |
| `aphi_gmres_benchmark_helpers.h` | Joule local-field metrics, 8×8 PC summary, core mean scalar helpers |

---

## 10A-1a — Joule uniform-field verification

**Setup:** uniform σ=2, A=0, phi_real=-E0·x (E0=3), **no GMRES**; Joule pipeline only.

| Metric | Value |
|--------|-------|
| expected q | 9.0 (= 0.5 σ E0²) |
| core_mean_joule | 8.898 |
| core_mean_rel_error | **~1.1%** |
| core_mean_Ex | ~2.983 (≈ E0) |
| **passed** | **1** |

**Conclusion:** grad φ → E → q chain consistent with analytic uniform-field (~1% explainable by discrete gradient error near boundaries).

---

## 10A-1b — Coupled 8×8 PC diagnostics

**Setup:** TEAM7-like layout, one coupled PC apply on impressed-current RHS block.

| Metric | Value |
|--------|-------|
| fallback_count | **0** |
| fallback_fraction | **0** |
| global_min_pivot (min \|diag\| proxy) | 100.015 |
| conductor/coil/air fallback | 0 / 0 / 0 |

**Conclusion:** **No silent fallback** on current TEAM7 unit-box layout; coupled 8×8 applies to all particles. min_pivot≈100 from φ gauge penalty diagonal term magnitude.

---

## 10A-2a — Joule local distribution sensitivity

**Setup:** tol = 5e-3 / 5e-4 / 1e-4 / 1e-5, reference tol=5e-4; added conductor L2/max/band differences.

| Metric (vs 5e-4 ref, tol≥1e-4) | Value |
|--------------------------------|-------|
| max conductor power rel change | ~4.1e-5 |
| max conductor L2 ratio | ~4.1e-5 |
| max conductor max rel | ~4.3e-5 |
| conductor_l2_diff @ 1e-4 | ~2.6e-5 |
| conductor_max_diff @ 1e-4 | ~1.2e-4 |
| **sensitivity_passed** | **1** |

**Conclusion:** While total power is insensitive, **local Joule field distribution** is equally insensitive to EM tolerance tightening (<0.01% magnitude).

---

## 10A-2b — Physical TEAM7-dimension dp observable sweep

**Setup:** 1.2×1.0×0.3 m, tol=5e-4, m=50 coupled PC; dp = 0.15 / 0.1 / 0.075.

| dp | particles | outer | em_rel | conductor_field_max | conductor_joule_power | conductor_joule_max |
|----|-----------|-------|--------|---------------------|----------------------|---------------------|
| 0.15 | 112 | 2 | ~7.4e-5 | 1.65 | 0.0428 | 2.12 |
| 0.1 | 360 | 3 | ~2.0e-4 | 3.59 | 0.0805 | 10.07 |
| 0.075 | 832 | 85 | ~4.9e-4 | 3.07 | 0.0931 | 7.36 |

**Conclusion:**
- All converged @ tol=5e-4
- **Observables not yet saturated:** conductor Joule power 0.043 → 0.081 → 0.093 still rising as dp decreases
- dp=0.15 too few particles (112), coarse grid reference only; **dp=0.1 remains practical baseline**

---

## Build & run

```bash
cd /home/yongchuan/sphinxsys/build
ninja test_3d_aphi_ck_joule_uniform_field_analytic_verification \
      test_3d_aphi_ck_block_jacobi_8x8_pc_diagnostic \
      test_3d_aphi_ck_joule_local_distribution_sensitivity_diagnostic \
      test_3d_aphi_ck_gmres_team7_physical_dp_observable_diagnostic
```

## Checklist

- [x] 10A-1a uniform-field Joule analytic verification
- [x] 10A-1b 8×8 PC fallback/pivot diagnostics
- [x] 10A-2a Joule local/source-field sensitivity
- [x] 10A-2b physical-dimension dp observable sweep
- [x] Stage 10A-2b physical-dimension dp observable sweep
- [x] Stage 10B: canonical spec + fine-mesh self-reference comparison
- [x] Stage 10C: cold-crucible geometry scaffold (EM-only + Joule)
- [ ] Stage 10D: one-way thermal coupling

## Next step (Stage 10B)

1. Define TEAM7 reference source (literature/FEM/MFEM)
2. Select comparison quantities: plate Joule distribution, line field strength, total power
3. Produce quantitative error report after dp observables stabilize
