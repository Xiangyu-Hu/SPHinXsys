# Stage 10B — Quantitative TEAM7 Reference Comparison

Date: 2026-05-21  
Branch: `feature/electromagnetic`

## Scope

Prepare quantitative TEAM7-like validation path before cold-crucible work.

| Item | Status |
|------|--------|
| Canonical case specification | ✅ `aphi_team7_canonical_case_ck.h` |
| Unified physical-box runner + observables | ✅ `runTeam7PhysicalCanonicalCase` |
| Conductor centerline field profile | ✅ (informational) |
| Fine vs coarse self-reference test | ✅ |
| External FEM/MFEM comparison | ❌ deferred (no reference data in repo) |

## Reference strategy

**Primary (Stage 10B):** fine-mesh **self-reference**
- Reference: `dp=0.075` on physical box 1.2×1.0×0.3 m
- Candidate: `dp=0.1` (engineering baseline)
- Compare conductor Joule power, conductor field max, centerline profile (informational)

**Regression anchors:** recorded observables from reference_dp run stored in `AphiTeam7CanonicalCaseSpec`.

**External FEM:** not in repository; import path deferred to future Stage 10B+.

## Canonical case parameters

| Parameter | Value |
|-----------|-------|
| Box | 1.2 × 1.0 × 0.3 m |
| Layout | unit-box fractions → physical regions |
| σ_air / σ_conductor / σ_coil | 1e-4 / 1.0 / 1e-4 |
| ω | 1.25 |
| φ gauge penalty | 100 |
| Impressed coil amplitude | 8 (z-directed) |
| GMRES | m=50, coupled 8×8 PC, tol=5e-4, max_outer=150 |
| reference_dp / candidate_dp | 0.075 / 0.1 |

## Test: `test_3d_aphi_ck_gmres_team7_quantitative_reference_comparison`

### Run results (2026-05-21)

| Run | dp | particles | outer | em_rel | conductor_joule | conductor_field_max |
|-----|-----|-----------|-------|--------|-----------------|---------------------|
| reference | 0.075 | 832 | 134 | ~4.9e-4 | **0.0930** | **3.07** |
| candidate | 0.1 | 360 | 3 | ~2.1e-4 | **0.0805** | **3.59** |

### Relative differences (candidate vs reference)

| Metric | Value | Threshold | Pass |
|--------|-------|-----------|------|
| conductor_joule_rel | **~13.4%** | 20% | ✅ |
| conductor_field_rel | **~17.0%** | 20% | ✅ |
| centerline_profile_rel | ~1.38 | 2.0 (informational) | ✅ logged |
| reference vs recorded joule | ~0.13% | 5% | ✅ |
| reference vs recorded field | ~0.11% | 5% | ✅ |
| **passed** | | | **1** |

### Interpretation

1. **Observable gap dp=0.1 vs 0.075 ~13–17%** — mesh not fully converged; dp=0.1 acceptable for engineering but not publication-grade quantitative claim.
2. **Centerline profile L2 ~138%** — coarse mesh sparsely populates x-bins; kept **informational** until finer mesh or binning fix.
3. **Regression anchors stable** — reference_dp reproducible within ~0.1% of recorded constants.

## New shared code

| File | Purpose |
|------|---------|
| `aphi_team7_canonical_case_ck.h` | Canonical spec + thresholds + regression anchors |
| `aphi_gmres_benchmark_helpers.h` | `runTeam7PhysicalCanonicalCase`, centerline profile, comparison helpers |

## Build & run

```bash
cd /home/yongchuan/sphinxsys/build
ninja test_3d_aphi_ck_gmres_team7_quantitative_reference_comparison
```

## Checklist

- [x] Canonical TEAM7-like spec documented in code
- [x] Coil/source/material/frequency/solver settings unified
- [x] Observable metrics: Joule power, field max, centerline profile
- [x] Fine-mesh self-reference comparison test
- [x] Regression anchors for reference_dp
- [ ] Import external FEM/MFEM reference data
- [ ] Tighten candidate_vs_reference thresholds after dp convergence study
- [ ] Import external FEM/MFEM reference data
- [x] Stage 10C: cold-crucible geometry scaffold
- [ ] Stage 10D: one-way thermal coupling (next)

## Next steps

1. **Optional:** dp=0.06 or 0.05 reference to reduce observable gap
2. **10B+:** import TEAM7 FEM baseline (plate Joule, B/H along lines)
3. **10C:** cold-crucible geometry scaffold (EM-only + Joule)
