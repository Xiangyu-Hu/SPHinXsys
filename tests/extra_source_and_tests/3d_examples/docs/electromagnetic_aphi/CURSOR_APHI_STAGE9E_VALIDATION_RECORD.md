# Stage 9E — Stronger Validation Record

Date: 2026-05-21  
Branch: `feature/electromagnetic`  
Solver default: GMRES **m=50** + **CoupledPointBlock8x8** PC

## 1. Scope & status

| Sub-stage | Content | Test | Status |
|-----------|---------|------|--------|
| **9E-1** | Interface-aware flux-matched φ MMS (continuous σ∂φ/∂n) | `test_3d_aphi_ck_gmres_interface_flux_matched_manufactured` | ✅ |
| **9E-2** | Real TEAM7 dimensions 1.2×1.0×0.3 m + fraction layout | `test_3d_aphi_ck_gmres_team7_physical_dimensions_smoke` | ✅ |
| **9E-3** | TEAM7 unit-box **dp** convergence sweep | `test_3d_aphi_ck_gmres_team7_dp_convergence_diagnostic` | ✅ (informational) |

**Explicitly out of scope (this stage):** change Laplace weights, Contact, cold crucible, full thermal feedback.

## 2. New / modified shared code

| File | Change |
|------|------|
| `aphi_benchmark_case_ck.h/.hpp` | `AphiTeam7PhysicalDimensions`, `buildTeam7LayoutForBox`, `AssignInterfaceFluxMatchedPhiFieldsCK` |
| `aphi_gmres_benchmark_helpers.h` | `hostInterfaceBandTrueResidualNorm` |

### 9E-1 manufactured field notes

- `phi_real`: piecewise linear in x, φ continuous at interface x=0.5, and `sigma_left * dphi/dx = sigma_right * dphi/dx`
- `A = 0`, `phi_imag = 0`
- Compared to 9C-2 separable sin/cos field: **continuous-field error reduced from O(1) to ~2e-6**

## 3. Run results (2026-05-21)

### 9E-1 Interface flux-matched MMS

**Setup:** unit box, dp=0.1, σ=2 / 1e-4 half-space @ x=0.5, tol=1e-5, m=50 (default)

| Metric | Value |
|--------|-------|
| particles | 1000 |
| outer / arnoldi | 2 / 50 |
| true_rel_res | **4.99e-6** |
| discrete_mms_defect | **4.99e-6** |
| continuous_field_relative_error | **1.84e-6** |
| interface_band_rel (±2dp @ x=0.5) | **1.89e-6** |
| **passed** | **1** |

**Conclusion:** Interface-aware manufactured solution significantly better than 9C-2 non-interface-aware field; interface-band residual same order as global, no abnormal localization.

**Retroactive update:** 9C-2 record §4 "future: interface-aware MMS" → **completed in 9E-1**.

---

### 9E-2 TEAM7 physical dimensions smoke

**Setup:** 1.2×1.0×0.3 m, dp=0.1, tol=5e-4, m=50 coupled PC, impressed coil amp=8

| Metric | Value |
|--------|-------|
| particles | 360 |
| outer / arnoldi | 2 / 50 |
| true_rel_res | **4.73e-4** |
| conductor_field_max | 3.59 |
| coil_field_max | 3.59 |
| **passed** | **1** |

**Implementation note:** With thin box z=0.3 m, default `core_shell=2*dp` empties z-direction core particles; this test uses  
`core_shell = min(2*dp, 0.125*min(L,H,W))` for regional field statistics.

**Conclusion:** Physical dimensions + fraction layout runs successfully; converges in 2 outer iterations at engineering tolerance 5e-4.

**Retroactive update:** STAGE9C record §8 item 5 "real dimensions deferred" → **smoke completed**; quantitative TEAM7 benchmark still pending.

---

### 9E-3 TEAM7 dp convergence sweep (informational)

**Setup:** unit box TEAM7 layout, tol=5e-4, m=50 coupled PC

| dp | particles | outer | true_rel_res | conductor_field_max | converged |
|----|-----------|-------|--------------|---------------------|-----------|
| 0.2 | 125 | 100 | 6.27e-3 | 0 (core too coarse) | ❌ |
| 0.1 | 1000 | 2 | **2.44e-4** | 1.91 | ✅ |
| 0.075 | 2197 | 20 | 4.80e-4 | 3.13 | ✅ |

**Conclusion:**

- dp=0.2 **too coarse**, too few particles and no effective conductor/coil sampling in core region → does not converge
- dp=0.1 is current **recommended baseline** (consistent with 9C-3, outer=2)
- dp=0.075 finer but needs more outer (20) at tol=5e-4; true_rel not monotonically better than dp=0.1 → **finest grid not mandatory under engineering tolerance**

## 4. Update vs prior "not done" items

| Prior record | Current status |
|----------|--------|
| 9C-2 interface-aware MMS | ✅ 9E-1 |
| Geometry: real TEAM7 dimensions | ✅ 9E-2 smoke (not quantitative benchmark) |
| dp convergence study | ✅ 9E-3 diagnostic |
| Stage 9E: 1e-5 closure | ⚠️ **does not block 9E**; 9D-2 already supports engineering acceptance **~2e-5** EM tolerance |
| Joule sensitivity | ✅ completed in 9D-2 (record retained) |
| Full thermal / cold crucible | ❌ still deferred |

## 5. Build & run

```bash
cd /home/yongchuan/sphinxsys/build
cmake .. -DSPHINXSYS_3D=ON -DSPHINXSYS_BUILD_3D_EXAMPLES=ON -DSPHINXSYS_BUILD_EXTRA_SOURCE_AND_TESTS=ON
ninja test_3d_aphi_ck_gmres_interface_flux_matched_manufactured \
      test_3d_aphi_ck_gmres_team7_physical_dimensions_smoke \
      test_3d_aphi_ck_gmres_team7_dp_convergence_diagnostic
```

## 6. Acceptance checklist (Stage 9E)

- [x] 9E-1: flux-matched interface MMS converges @ 1e-5
- [x] 9E-1: discrete defect & interface-band residual logged
- [x] 9E-2: physical 1.2×1.0×0.3 m TEAM7 layout smoke @ tol=5e-4
- [x] 9E-3: dp sweep 0.2 / 0.1 / 0.075 logged
- [x] PROGRESS / STAGE9C record / plan doc updated
- [ ] Optional: strict 1e-5 @ full σ contrast TEAM7 (best ~2e-5, 9D-1)
- [ ] Optional: quantitative TEAM7 reference comparison
- [ ] Cold crucible / thermal feedback (post-9E)

## 7. Suggested next steps

1. For publication-level TEAM7: introduce reference solution or experimental data for quantitative error
2. Physical-dimension case can sweep dp (similar to 9E-3) to confirm convergence behavior under thin box
3. Cold crucible geometry + Contact explicitly deferred
