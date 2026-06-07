# Stage 10.17 — Native TEAM7 SI geometry + coil-path source + Bz reference

> **Date**: 2026-06-03 (closed loop); **Update**: 2026-06-05 (post-10.17 MR / rerun record)  
> **Status**: **P0–P4 passed=1** on SI reload (`stl_scale=0.001`, `dp_0=0.006 m`)  
> **Branch**: `feature/electromagnetic`

---

## 0. One-Line Conclusion

The native SPHinXsys path (STL → relax → reload → three-body Contact A-φ GMRES → Bz probe) is **order-of-magnitude aligned** with the TEAM7 muFEM Bz reference under **SI meters**: A1-B1 peak **~8.05 mT vs ref ~7.81 mT**, profile L2 **~19% / ~24%** (frozen sign `(+real, -imag)`).

---

## 1. Frozen Decisions (Implemented)

| Item | Decision |
|----|------|
| Geometry/particles | STL `scale=0.001`, rerun relax/reload; do not multiply by 0.001 after reload |
| Reluctance | native path `ν = 1/μ₀` |
| Coil source | `rhs = J₀ × tangent`, `J₀ = NI/A_eff`, **do not multiply** by μ₀/ω/V |
| Centerline | `R_out=50 mm`, `R_in=25 mm` fixed; bbox used only for positioning |
| Bz sign P4 | `phase0 = +Bz_real×1000 mT`, `phase90 = -Bz_imag×1000 mT` → **(1,-1)** |
| GMRES | No tuning; `final_true_rel ~ 0.005` |

---

## 2. Test Matrix

| ID | Executable | Purpose | Typical Result |
|----|------------|------|----------|
| P0 | `test_3d_aphi_ck_team7_native_geometry_audit` | SI bbox, particle count | plate ~291 mm |
| P1 | `test_3d_aphi_ck_team7_native_source_rhs_audit` | NI=2742 A, `J₀~1.02e6 A/m²` | `i_eff_ratio≈1` |
| P2 | `test_3d_aphi_ck_team7_native_vacuum_source_b_sanity` | plate σ=0, mT only | ~16 mT, profile gate off |
| P3 | `test_3d_aphi_ck_team7_native_geometry_bz_reference_probe` | Full TEAM7 50 Hz Bz | peak ~8 mT, profile ≤30% |
| P4 | (merged into P3) | Frozen sign + scan consistency | `frozen_matches_scan=1` |

Particle generation: `particle_generation_em` (`tests/.../particle_generation_em/particle_generation_em.cpp`).

---

## 3. Key Code

| Path | Content |
|------|------|
| `diagnostics/aphi_team7_native_reload_geometry_helpers.h` | Three-body case, GMRES, ν_SI |
| `diagnostics/aphi_team7_native_coil_source_helpers.h` | Centerline, tangent, `AssignTeam7CoilPathImpressedCurrentRhsCK` |
| `diagnostics/aphi_team7_native_bz_reference_probe_helpers.h` | Bz probe, P4 frozen sign |
| `cmake/team7_native_reload_sync.cmake` | Copy reload at configure time |
| `sync_team7_native_si_reload.sh` | Manual sync to each test `bin/` |

---

## 4. User Validation Data (2026-06-03)

```text
P1: j0_used=1.019e+06  i_eff_ratio≈1.0  rhs_l2_team7=41566
P2: max_abs_Bz_A1_B1_mT=15.73  final_passed=1 (magnitude_only)
P3: max_abs_Bz_mT=8.05  profile_real=0.189  imag=0.239
    sign=(1,-1)  frozen_matches_scan=1  final_true_rel=0.005
```

---

## 5. Known Limitations / Follow-Up (Non-Blocking)

- Profile L2 ~20%: overall passes 25% gate; peak at x≈126 mm differs from ref by ~3%.
- Pointwise error: near reference zero-crossing (e.g. **x=90 mm**, ref Bz≈0.036 mT) inflates relative point error; diagnostics skip |ref|<0.1 mT.
- Reload sync: `sync_team7_native_si_reload.sh` + `cmake` auto-copies from `particle_generation_em/bin/reload` via `team7_native_reload_sync.cmake` (verified 2026-06-03).
- One-click smoke: `run_team7_native_si_smoke.sh` (default P0+P1+P3; `RUN_TEAM7_P2=1` includes vacuum case).
## Next (validation focus — Bz vs muFEM/COMSOL, not Joule)

1. **High-res small air**: `--team7-dp=0.003 --team7-air=small` → relax → P3 Bz (`run_team7_highres_bz.sh`).
2. **Big air** (after high-res baseline): `--team7-air=big` with chosen dp.
3. Profile gate tuning only after domain + dp sensitivity.

Geometry CLI / env: `aphi_team7_native_geometry_config.h`; metadata `./reload/team7_native_geometry.txt` auto-synced to tests.

Stage 10.18 Joule post: **deferred** (TEAM7 official benchmark has no thermal target).

---

## 6. Quick Run Reference

```bash
# 1) Generate SI reload (once)
cd ~/sphinxsys/build/tests/extra_source_and_tests/3d_examples/particle_generation_em/bin
./particle_generation_em

# 2) Sync + configure (cmake auto-copies reload; see configure log TEAM7 native reload: copied from ...)
~/sphinxsys/tests/extra_source_and_tests/3d_examples/sync_team7_native_si_reload.sh
cd ~/sphinxsys/build && cmake ..

# 3) One-click smoke (P0 + P1 + P3)
~/sphinxsys/tests/extra_source_and_tests/3d_examples/run_team7_native_si_smoke.sh
# Optional with P2: RUN_TEAM7_P2=1 ~/sphinxsys/tests/extra_source_and_tests/3d_examples/run_team7_native_si_smoke.sh
```

Report directory: `team7_native_bz_reference_report/` (`summary.csv`, `comparison_A1_B1_vs_reference.csv`).

---

## 7. Closed-Loop Sign-Off (2026-06-03)

```text
sync + cmake reload copy     OK
P0 geometry_audit            passed=1
P1 source_rhs_audit          passed=1
P2 vacuum (optional)         passed=1  profile_gate=off
P3 bz_reference + P4 sign    passed=1  peak_Bz≈8.05 mT  profile≈19/24%
```

---

## 8. Work After 10.17 (2026-06-04 — 06-05)

### 8.1 Multi-Resolution (AdaptiveBody)

- `particle_generation_em` + P3 solve: `AdaptiveBody`, `SmoothingLengthRatio` written to reload XML.
- `aphi_team7_native_geometry_config.h`: `--team7-case` / `--team7-reload-case` / air preset CLI.
- `aphi_team7_native_mr_probe_helpers.h`: MR probe sampling by adaptive kernel support.
- Scripts: `run_team7_uniform_legacy_bz_vtp.sh`, `run_team7_highres_bz.sh`.

### 8.2 Bz Comparison (P3, `--team7-reload-case=...`)

| Reload case | Air particles | L2 real | Peak Bz | P3 | Notes |
|-------------|----------|---------|---------|-----|------|
| small uniform 6 mm | 173133 | 20.4% | 7.89 mT | pass | 6/4 |
| small MR levels=1 | 233004 | 24.1% | 7.90 mT | pass | 6/4; coarsened 12 mm |
| legacy uniform (`uniform_legacy_dp6`) | 1122763 | **13.1%** | 7.11 mT | pass | 6/5 rerun |
| legacy MR `multires_small_dp6_l2` | — | 100% | 0 | **fail** | GMRES nan |

Default P3 without CLI uses `bin/reload/` (same as `uniform_legacy_dp6`).

### 8.3 ParaView

- Air `NativeProbeBReal` Z component shows dipole-like Bz distribution (consistent with probe CSV).
- VTP includes three-body A/φ/B; Plate also has Jy; Coil has `Team7CoilSourceTangent`.
- Air `ElectricFieldAReal` is computed but not written to VTP (to be added during cleanup).

### 8.4 Still Incomplete

- High-res `dp=3 mm` (`run_team7_highres_bz.sh`) not finished.
- Big air + MR not started.
- TEAM7 quantitative profile target <15% not closed (legacy uniform already ~13% real).
- Stage 10.18 Joule post still deferred.

See: `CURSOR_APHI_TEAM7_SMALL_AIR_MULTIRES_PLAN.md` (aligned with code, revised 2026-06-05).
