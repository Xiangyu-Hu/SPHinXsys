# Stage 10C — Cold-Crucible Geometry Scaffold Record

Date: 2026-05-21  
Branch: `feature/electromagnetic`

## Scope

First cold-crucible **geometry scaffold** per ChatGPT Stage 10C plan:

| Item | Status |
|------|--------|
| Geometry + region tagging | ✅ |
| Dual-coil impressed-current source | ✅ |
| Melt / crucible wall / coil / air materials | ✅ |
| EM-only GMRES solve | ✅ |
| Joule source output (regional) | ✅ |
| σ(T)/ν(T) thermal feedback | ❌ deferred (10D/10E) |
| Contact / CAD geometry | ❌ deferred |

**This is a simplified box scaffold**, not a publication-grade cold-crucible replica.

## Layout (unit box fractions → scalable physical box)

Schematic (x–y mid-plane):

```
[ coil_L ]  [ crucible shell + melt ]  [ coil_R ]
```

| Region | Unit-box fractions (x, y, z) | σ | ν |
|--------|------------------------------|---|---|
| **Melt** | [0.38,0.62]³ center, z∈[0.25,0.75] | 1.0 | 1.0 |
| **Crucible wall** | shell [0.32,0.68]³ \ melt | **1e4** | 1.0 |
| **Coil L/R** | x∈[0.10,0.26] & [0.74,0.90], y,z∈[0.20,0.80] | 1e-4 | 1.0 |
| **Air** | remainder | 1e-4 | 1.0 |

Material priority: **melt > crucible_wall > coil > air**.

Default physical box: **1.0×1.0×1.0 m** (`AphiColdCruciblePhysicalDimensions`); scale via `buildColdCrucibleLayoutForBox`.

## New shared code

| File | Purpose |
|------|---------|
| `aphi_cold_crucible_case_ck.h/.hpp` | Layout, case spec, `AssignColdCrucibleRegionMaterialsCK`, `AssignColdCrucibleCoilSourceCK` |
| `aphi_gmres_benchmark_helpers.h` | `runColdCrucibleScaffoldCase`, region counts, regional Joule helpers |

## Solver settings (canonical)

| Parameter | Value |
|-----------|-------|
| dp | 0.1 |
| ω | 1.25 |
| φ penalty | 100 |
| Coil impressed amplitude | 8 (z-directed) |
| GMRES | m=50, coupled 8×8 PC, tol=**5e-4**, max_outer=150 |

## Test: `test_3d_aphi_ck_gmres_cold_crucible_scaffold`

### Results (2026-05-21)

| Metric | Value |
|--------|-------|
| particles | 1000 |
| em_rel | **~5.4e-6** |
| outer | 2 |
| converged | ✅ |
| region counts (melt / crucible / coil / air) | 24 / 72 / 144 / 760 |
| melt_joule_power | ~1.1e-7 |
| crucible_joule_power | ~**4.8e-5** |
| coil_joule_power | ~2.2e-8 |
| total_joule_power | ~4.8e-5 |
| melt_field_max | ~0.0040 |
| crucible_field_max | ~0.0040 |
| coil_field_max | ~0.040 |
| min_joule_source | > 0 |
| **passed** | **1** |

### Interpretation

1. **EM solve stable** — converges in 2 outer iterations @ tol=5e-4 (better than TEAM7 high-contrast case on first scaffold try).
2. **Region tagging OK** — all four regions populated.
3. **Joule dominated by crucible wall** — expected with σ_wall=1e4 ≫ σ_melt=1; melt Joule is small but **> 0**.
4. **Not yet engineering validation** — box regions, no curved crucible, no skin depth mesh study, no thermal loop.

## Explicit non-goals (Stage 10C)

- No Contact between melt/crucible
- No temperature-dependent σ/ν
- No one-way / two-way thermal coupling (Stage 10D/10E)
- No claim of quantitative cold-crucible melting physics

## Build & run

```bash
cd /home/yongchuan/sphinxsys/build
ninja test_3d_aphi_ck_gmres_cold_crucible_scaffold
```

## Checklist

- [x] 10C-1 geometry + region tagging
- [x] 10C-2 dual-coil source tagging
- [x] 10C-3 melt/crucible/coil/air materials
- [x] 10C-4 EM-only solve
- [x] 10C-5 Joule source regional output
- [ ] 10C+ curved geometry / physical crucible dimensions study
- [ ] Stage 10D one-way thermal coupling

## Next steps

1. **10C+:** refine crucible σ_wall vs σ_melt contrast; verify melt Joule fraction
2. **10D:** map Joule → thermal equation (fixed σ/ν)
3. **External FEM:** optional comparison for coil–melt coupling

## Files to upload for ChatGPT discussion

**Docs:**
- `tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10C_SCAFFOLD_RECORD.md` (this file)
- `tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10A_PROGRESS.md`
- `tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10B_QUANTITATIVE_RECORD.md`
- `tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_CK_MATRIX_FREE_PROGRESS.md`

**New code:**
- `tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/aphi_cold_crucible_case_ck.h`
- `tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/aphi_cold_crucible_case_ck.hpp`

**New test:**
- `tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_gmres_cold_crucible_scaffold/test_3d_aphi_ck_gmres_cold_crucible_scaffold.cpp`
