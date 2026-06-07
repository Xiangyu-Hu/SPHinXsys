# Stage 10.17 Native TEAM7 Reload + Bz Reference — Handoff for GPT Discussion

> Repo root (SSH): `/home/yongchuan/sphinxsys`  
> Branch: `feature/electromagnetic`  
> Date: 2026-06-03 (initial draft); **historical document** — Q1–Q3 closed on 6/3  
> Purpose: Previously used for GPT discussion on units/source/sign; **current progress see** `CURSOR_APHI_STAGE10_17_TEAM7_NATIVE_RECORD.md` §8.

---

## 1. What we are trying to do

- Use **native SPHinXsys** STL geometry + particle relax + **reload XML** (not hand-made box TEAM7).
- Run **SYCL CK** frequency-domain **A-φ** solve (three-body **Contact** + matrix-free **GMRES**).
- Compare **Bz** on muFEM/COMSOL probe lines (A1-B1, A2-B2) against CSVs under `reference_data/team7/`.

**Status snapshot**

| Date | Conclusion |
|------|------|
| 2026-06-03 | P0–P4 closed: SI + coil-path source + sign (1,-1); small air peak ~8 mT, L2 ~19–24% |
| 2026-06-05 | legacy uniform rerun: peak 7.11 mT, L2 real **13.1%**; small MR L1 pass; legacy MR L2 **fail** |

**Open (not unit/RHS):** profile <15%, high-res dp=3 mm, big air, legacy MR L2 debug. See `CURSOR_APHI_TEAM7_SMALL_AIR_MULTIRES_PLAN.md`.

---

## 2. Case setup (native reload path)

### 2.1 Particle geometry (user-validated relax)

| Item | Value |
|------|--------|
| Program | `particle_generation_em` |
| Bodies | `Coil`, `Plate`, `Air` (names must match reload XML) |
| STL | `./input/coil.stl`, `./input/plate.stl` |
| Small air box (active) | lower `(-50,-50,-50)`, upper `(350,350,200)` |
| Big air box (commented) | lower `(-1353,-1353,-300)`, upper `(1647,1647,449)` |
| Global resolution | `dp_0 = 6.0` (same units as particle positions — **mm-like**, see §6) |
| Relax | 1000 steps CK level-set particle relaxation |
| Material in relax | **None** (only `Solid` placeholder; no σ/ν) |
| Reload output | `./reload/Coil_rld.xml`, `Plate_rld.xml`, `Air_rld.xml` |

**Source:**  
`/home/yongchuan/sphinxsys/tests/extra_source_and_tests/3d_examples/particle_generation_em/particle_generation_em.cpp`

**Typical particle counts (user small-air run):**  
Air ~172414, Plate ~6302, Coil ~7712.

### 2.2 EM solve case (test cases)

| Item | Value |
|------|--------|
| Domain | Same small air bounding box as relax |
| Bodies | `RealBody Air`, `SolidBody Coil`, `SolidBody Plate` |
| Particles | `setReloadParticles(true)`, `generateParticles<Reload>(body_name)` |
| Discretization | SPH matrix-free CK; **Contact** between air↔coil, air↔plate, coil↔air, plate↔air |
| Execution | `MainExecutionPolicy` / `par_ck` (SYCL device CK) |

**σ / ν (assigned in solve, not in relax):**

| Body | σ (S/m) | ν | Role |
|------|---------|---|------|
| Plate | `3.526e7` | `1.0` | Aluminum conductor (TEAM7) |
| Coil | `0` | `1.0` | Source-only region (no eddy in coil) |
| Air | `0` | `1.0` | Non-conducting air |

Dynamics: `SetAphiMaterialPropertiesCK` per body after `InitializeAphiVariablesCK`.

### 2.3 Frequency and solver

| Parameter | Default in native tests |
|-----------|-------------------------|
| `omega` | `2π × 50` Hz |
| `phi_gauge_penalty` | `100` |
| GMRES tolerance | `5e-3` |
| Restart dimension | `30` |
| Max outer iterations | `80` |
| Preconditioner | Coupled contact GMRES (`defaultCoupledContactGMRESOptions`) |
| LHS options | `use_phi_gauge_penalty = true`; no A-divergence penalty in this path |

**Diagnostic pass (smoke / Bz):** also accept `final_true_rel ≤ 0.1` if strict GMRES flag fails.

---

## 3. How excitation current is applied (critical)

> **2026-06-05 errata**: §3 below describes the **pre-10.17-fix** placeholder source. The current production path uses `AssignTeam7CoilPathImpressedCurrentRhsCK` (`aphi_team7_native_coil_source_helpers.h`), `J₀ = NI/A_eff × tangent`, P1 `i_eff_ratio≈1`. The following paragraphs are for historical reference only.

We do **not** implement TEAM7 2742 turns × 1 A/turn yet. We use a **placeholder** copied from simplified TEAM7 benchmarks.

### 3.1 RHS assignment

- Class: `AssignUniformImpressedCurrentRhsCK` (native-only helper).
- Applied on: **`Coil` body only** (all coil particles).
- Air / Plate RHS: **zeroed** (`AphiZeroBlockCK` on `rhs` block).

Per particle `i` on coil:

```text
rhs_a_real[i] = impressed_current_amplitude * current_real
rhs_a_imag[i] = impressed_current_amplitude * current_imag
```

Defaults:

```text
current_real = (0, 0, 1)   // unit vector, z-direction
current_imag = (0, 0, 0)
impressed_current_amplitude = 8.0   // NOT 2742 AT
```

This is **not** `AssignImpressedCurrentRhsCK` (smooth box profile in coil region) used by unit-box TEAM7 cases. It is a **uniform** RHS on every coil particle.

### 3.2 What is solved

Unknowns: complex A-φ block per particle (`a_real`, `a_imag`, `phi_real`, `phi_imag`).

Solver: `AphiMultiBodyContactGMRESSolverCK<MainExecutionPolicy>` — applies discrete LHS (harmonic grad/div/reaction + φ gauge) matrix-free across three Contact-linked bodies.

After solve:

1. `runTeam7NativeContactJoulePipeline` — grad φ, E, J, Joule on each body.
2. (Bz test only) `execBodyCurlBFromADiagnostic` — `B = curl A` with corrected gradient on each body; fields `NativeProbeBReal`, `NativeProbeBImag`.

### 3.3 Post-processing for Bz probe

- Probe lines in **mm coordinates** (same as STL), 17 points: `x = 0, 18, …, 288` mm.
- A1-B1: `y=72`, `z=34` mm.
- A2-B2: `y=144`, `z=34` mm.
- Sample: **nearest particle** among Air + Coil + Plate particles.
- Convert to mT: `our_Bz_*_mT = sign * 1000 * B_*[2]` (Tesla→mT assumption on curl A z-component).
- Sign audit: try `(sign_real, sign_imag) ∈ {(±1,±1)}`; pick min L2 profile error vs reference `Bz_50Hz_phase0_mT` / `phase90_mT`.

muFEM doc convention (from Stage 10.17 plan, **not yet enforced as default**):

```text
Bz_complex_mT = 1e3 * (-B_real.z + i * B_imag.z)
```

Last run picked **(+1, +1)** by error minimization.

---

## 4. Reference data

**Repo path:**

`/home/yongchuan/sphinxsys/tests/extra_source_and_tests/3d_examples/reference_data/team7/`

| File | Content |
|------|---------|
| `TEAM7_Bz_A1_B1_reference_mT.csv` | x=0..288 mm, Bz 50/200 Hz phase0/90 (mT) |
| `TEAM7_Bz_A2_B2_reference_mT.csv` | A2-B2 line |
| `TEAM7_Bz_*_reference_Gauss.csv` | Raw Gauss |
| `TEAM7_probe_definitions.csv` | Probe metadata |

Reference comparison uses **50 Hz** columns: `Bz_50Hz_phase0_mT`, `Bz_50Hz_phase90_mT`.

---

## 5. Test cases and build outputs

### 5.1 Tests (source)

| Test | Path |
|------|------|
| Particle relax only | `tests/extra_source_and_tests/3d_examples/particle_generation_em/` |
| Geometry audit (no solve) | `tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_team7_native_geometry_audit/` |
| Source-driven smoke | `tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_team7_native_geometry_source_driven_smoke/` |
| **Bz vs reference** | `tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_team7_native_geometry_bz_reference_probe/` |

### 5.2 Shared helpers (source)

| File | Role |
|------|------|
| `extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_team7_native_reload_geometry_helpers.h` | Geometry, reload, materials, GMRES, smoke |
| `extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_team7_native_bz_reference_probe_helpers.h` | Bz probe, CSV, sign audit |
| `extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_curl_b_dual_track_diagnostic_helpers.h` | `B = curl A` |
| `extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_probe_metric_helpers.h` | Nearest-particle line sampling |

### 5.3 Build binaries (after `ninja`)

| Target | Binary directory |
|--------|------------------|
| `particle_generation_em` | `build/tests/extra_source_and_tests/3d_examples/particle_generation_em/bin/` |
| `test_3d_aphi_ck_team7_native_geometry_audit` | `build/.../test_3d_aphi_ck_team7_native_geometry_audit/bin/` |
| `test_3d_aphi_ck_team7_native_geometry_source_driven_smoke` | `build/.../test_3d_aphi_ck_team7_native_geometry_source_driven_smoke/bin/` |
| `test_3d_aphi_ck_team7_native_geometry_bz_reference_probe` | `build/.../test_3d_aphi_ck_team7_native_geometry_bz_reference_probe/bin/` |

**Run cwd:** always `.../bin/` (needs `./input/*.stl`, `./reload/*_rld.xml`; Bz test also `./reference_data/team7/*.csv` — CMake copies reference into bin on build).

### 5.4 Runtime inputs user must have in `bin/`

```text
bin/input/coil.stl
bin/input/plate.stl
bin/reload/Coil_rld.xml
bin/reload/Plate_rld.xml
bin/reload/Air_rld.xml
bin/reference_data/team7/TEAM7_Bz_A1_B1_reference_mT.csv   # Bz test only
bin/reference_data/team7/TEAM7_Bz_A2_B2_reference_mT.csv
```

Generate reload by running `particle_generation_em` from its `bin/`, then copy `reload/` + `input/` to other test `bin/` folders.

---

## 6. Recorded outputs

### 6.1 `test_3d_aphi_ck_team7_native_geometry_audit`

Stdout only, e.g.:

- `reload_files_present`, `air/plate/coil_particles`
- Bounding boxes (printed with `_mm` suffix — coordinates are STL/mm scale)
- `probe_A1B1_y_mm=72`, `probe_A1B1_z_mm=34`, etc.
- VTP: `bin/output/` (particle positions)

### 6.2 `test_3d_aphi_ck_team7_native_geometry_source_driven_smoke`

Stdout:

- `passed`, `converged`, `diagnostic_converged`, `final_true_rel`
- `coil_rhs_norm`, `plate/air/coil_joule_integral`, `air_to_conductor_joule`
- `finite_fields`, particle counts

Optional VTP under `bin/output/` (plate A, φ, E, J, σ).

### 6.3 `test_3d_aphi_ck_team7_native_geometry_bz_reference_probe`

Stdout:

- `solver_ok`, `final_true_rel`
- `A1_B1_sign=(sign_real,sign_imag)`
- `A1_B1_profile_rel_err_real/imag`, `A1_B1_max_nearest_distance_mm`
- A2-B2 profile errors (informative; gate focuses on A1-B1)

**Report directory:** `bin/team7_native_bz_reference_report/`

| File | Columns / content |
|------|-------------------|
| `probe_A1_B1_Bz.csv` | x_mm, y_mm, z_mm, nearest_distance_mm, Bz in T, our_Bz_*_mT, our_minus_* |
| `comparison_A1_B1_vs_reference.csv` | x_mm, our_Bz_real/imag_mT, ref phase0/90 mT, nearest_distance_mm |
| `probe_A2_B2_Bz.csv` | same for A2-B2 |
| `comparison_A2_B2_vs_reference.csv` | same |
| `summary.csv` | final_true_rel, solver_ok, signs, profile rel errs, max nearest distance |

**Observed last run (illustrative):**

- `final_true_rel ≈ 0.047`
- `A1_B1_profile_rel_err_* ≈ 1.0` (at 100% gate)
- `our_Bz_*_mT ~ 1e-9`, `ref ~ 0.5–8 mT`

---

## 7. Plans / records to upload

| Document | Path |
|----------|------|
| Stage 10.17 plan (TEAM7 native + Bz) | `tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_17_TEAM7_NATIVE_GEOMETRY_PLAN.md` |
| **This handoff** | `tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_17_NATIVE_RELOAD_BZ_HANDOFF_FOR_GPT.md` |
| Stage 10.16 review (context) | `/home/yongchuan/sphinxsys/CURSOR_APHI_STAGE10_16_NEXT_PLAN_AFTER_1014C_1015_REVIEW.md` |
| Legacy full TEAM7 (reference implementation) | `tests/extra_source_and_tests/3d_examples/legacy_aphi_archive/particle_generation_em/particle_generation_em.cpp` |

---

## 8. Questions for GPT (priority order)

### Q1 — Unit system (blocker)

- Particle positions from STL reload appear to be in **mm** (e.g. plate extent ~O(10³), TEAM7 plate 294 mm).
- Simplified CK TEAM7 benchmarks use **meters** (`AphiTeam7PhysicalDimensions`: 1.2 m × 1.0 m × 0.3 m) and SI `sigma`, `omega=2π×50`.
- Discrete grad/div use particle coordinates directly.

**Ask GPT:** Should we (a) scale all positions to meters after reload, (b) rescale EM equations for mm, or (c) rescale STL to meters at import? What is consistent with TEAM7 reference Bz in mT?

### Q2 — Source strength / 2742 AT

- Current RHS: uniform `J_rhs ∝ (0,0,1) × 8` on coil particles — **not** 2742 A·turn, not solenoidal, not legacy `TEAM7_B_CURVE_*` path.

**Ask GPT:** Map TEAM7 **2742 turns × 1 A/turn** to our `AssignUniformImpressedCurrentRhsCK` or `AssignImpressedCurrentRhsCK` / solenoidal curl-C formulation. What amplitude and direction give COMSOL/muFEM-comparable Bz (mT)?

### Q3 — Sign convention

- Plan: `Bz_complex_mT = 1e3 * (-B_real.z + i * B_imag.z)`.
- Audit picked `(+1,+1)` on near-zero fields.

**Ask GPT:** Correct time-harmonic convention for comparing phase0/phase90 columns to our complex A-φ solve at `omega=2π×50`.

### Q4 — Probe vs physics

- Nearest-particle sampling in mm grid; max distance ~4–8 mm (< `3*dp_0=18` mm gate).
- Is A1-B1 (y=72,z=34) in **air** between coil and plate in our STL frame? Confirm line is inside relaxed air particles.

### Q5 — What “passed” means

- Bz test `passed=1` does **not** mean reference match; only pipeline + loose L2 ≤ 100%.
- Propose staged gates: finite Bz magnitude, order-of-magnitude vs ref, then profile shape, then <50% / <20% error.

---

## 9. Suggested upload bundle for ChatGPT

**Minimum:**

1. This file (`CURSOR_APHI_STAGE10_17_NATIVE_RELOAD_BZ_HANDOFF_FOR_GPT.md`)
2. `CURSOR_APHI_STAGE10_17_TEAM7_NATIVE_GEOMETRY_PLAN.md`
3. `aphi_team7_native_reload_geometry_helpers.h` (excitation + materials + solve)
4. `aphi_team7_native_bz_reference_probe_helpers.h` (probe + comparison)
5. `particle_generation_em.cpp`
6. `team7_native_bz_reference_report/comparison_A1_B1_vs_reference.csv` (from last run)
7. `team7_native_bz_reference_report/summary.csv`
8. `reference_data/team7/TEAM7_Bz_A1_B1_reference_mT.csv`
9. `reference_data/team7/TEAM7_probe_definitions.csv`

**Optional:** legacy `particle_generation_em.cpp` (11772 lines) for how TEAM7 used to set `TEAM7_B_CURVE_*` and env-based current.

---

## 10. Contrast: simplified TEAM7 box cases (not this path)

These **do not** use STL reload; they use a **1.2×1.0×0.3 m** box and region fractions:

- `test_3d_aphi_ck_simplified_team7_source_driven`
- `AssignTeam7LikeRegionMaterialsCK` + `AssignImpressedCurrentRhsCK` with smooth box in coil region

Native path is intentionally separate; do not mix numbers between paths without unit conversion.
