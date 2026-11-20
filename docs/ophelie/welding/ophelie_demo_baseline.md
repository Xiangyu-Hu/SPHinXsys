# OPHELIE Induction Welding Demo — Gate 0 Baseline

**Date:** 2026-10-10  
**Plan:** [OPHELIE_Induction_Welding_Demo_Cursor_Plan.md](OPHELIE_Induction_Welding_Demo_Cursor_Plan.md)  
**Status:** Gate 0–4 **PASS** — Case A + Case B VTP.  
**Nature:** demonstration baseline — **not** validated welding-process prediction.

Smoke: `gate0_smoke/`, `gate12_smoke/`, `gate3_smoke/`, `gate4_smoke/`.  
**VTP:** [case_a_vtp/](case_a_vtp/), [case_b_vtp/](case_b_vtp/) (`eta_joint≈0.49` at coil `z=0` vs `≈0.21` at `z=0.025`).

Also mirrored at: `docs/ophelie_demo_baseline.md`.

---

## 1. Git / worktree state

| Item | Value |
|---|---|
| Branch | `feature/EM-Biot-Savart` @ `67fef44ef` (`origin/feature/EM-Biot-Savart`) |
| Dirty tree | **Yes** — many uncommitted Doc105/French hydro + EM edits (~180 paths). Do **not** overwrite. |
| Suggested demo branch | `feature/ophelie-induction-welding-demo` (not created yet; create from current branch after Gate 0) |
| `nvidia-smi` | Failed (driver not visible); SYCL OPHELIE smokes still ran successfully on this machine |

**Discipline:** New demo case must not modify French stirring/natural convection parameters or their output directories. Prefer a new `test_3d_ophelie_induction_demo/` target.

---

## 2. Module map (reuse, do not reinvent)

### 2.1 Shared OPHELIE library

Root include: `tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie.h`

| Concern | Path / entry |
|---|---|
| Parameters / CLI | `electromagnetic_ophelie_parameters.h`, `electromagnetic_ophelie_cli.h` |
| Field names | `electromagnetic_ophelie_field_names.h` (`ACoil*`, `BCoil*`, `ASrc*`, `BSrc*`, `Phi*`, `E*`, `J*`, `JouleHeatEdgeReconComplex`, …) |
| Multiloop filament Biot–Savart | `electromagnetic_ophelie_multiloop_source.h` → `OphelieMultiloopCoilSpec`, `applyMultiloopFilamentBiotToGlass` |
| Particle-coil Biot (legacy) | `electromagnetic_ophelie_biot_savart.h` → `ComputeOphelieCoilToGlassBiotSavartCK` |
| Current moments | `electromagnetic_ophelie_current_moment_source.h` |
| φ / edge-flux production | `electromagnetic_ophelie_phi.h`, `electromagnetic_ophelie_edge_flux.h` |
| Self-induction Picard (optional) | `stage2/electromagnetic_ophelie_self_induction.h` |
| French geometry helpers | `electromagnetic_ophelie_french_reduced_geometry.h` (CLI coil/glass; **glass melt**, not Al welding) |
| EM→thermal one-way | `diagnostics/electromagnetic_ophelie_joule_to_heat_one_way.h` |

### 2.2 Existing cases relevant to the demo

| Case | Role for welding demo |
|---|---|
| `test_3d_ophelie_biot_savart_circular_loop_axis` | Axis `Bz` vs analytic — **Biot gate** |
| `test_3d_ophelie_biot_savart_multiloop_axis` | Multiloop axis — **Biot gate** |
| `test_3d_ophelie_analytic_induction_simple` | Prescribed `A` + edge-flux φ/J/Q smoke (no coil BS) |
| `test_3d_ophelie_french_aind_diagnostic` | Multiloop Biot → edge-flux → Joule on French reload |
| `test_3d_ophelie_french_stirring_em` / `french_natural_convection_frozen_q` | Full production CLI + VTP writers — **template**, do not overwrite |
| `test_3d_ophelie_analytic_cylinder_skin` | Analytic skin/self-ind — validation, not demo geometry |

**No JSON config framework** for OPHELIE demos today → use C++/CLI (as the plan requires).

---

## 3. Production EM path (real, not assumed)

```text
A_coil, B_coil ← multiloop filament Biot–Savart (peak phasor, real A_coil typical)
A_src ← A_coil   (+ optional A_ind Picard / one-way feedback)
φ_r, φ_i ← PCG under --ophelie-current-form=edge-flux
E/J ← edge reconstruction (EEdgeRecon / JEdgeRecon)
q_J ← JouleHeatEdgeReconComplex   (peak-phasor convention in code comments)
P_joule ← volume integral of recon Q (calibration uses P_recon)
```

**Minimum CLI for production EM:**

```text
--ophelie-current-form=edge-flux
--ophelie-edge-flux-complex=1
```

Common geometry / drive (French helpers; demo should copy pattern with Al parameters):

```text
--dp= --glass-radius= --glass-height=
--frequency= --sigma=
--coil-radius= --coil-num-loops= --coil-segments-per-loop=
--coil-z-min= --coil-z-max= --coil-current=
--reload-dir=   # optional particle reload
--self-induction   # optional; demo Gate 2/3 default = coil-only source field
```

Other useful knobs: `--phi-pcg-tol=`, `--phi-pcg-max-iter=`, `--target-power=`, `--no-power-scaling`, `--ophelie-aind-one-way-feedback=`.

**Axis convention:** French / multiloop stacks use **z** as coil axis. Keep Case A/B on **z**; do not mix y-axis modules.

**Known physics limit for Al demo:** plan `f=10 kHz`, `σ≈3.5e7` → skin depth ~0.85 mm ≪ `dp=4 mm`. Coarse mesh can show workflow / |B| and qualitative |J|/q_J localization only — **must label as demo / not resolved skin**.

---

## 4. Gate 0 smoke results (2026-10-10)

| Test | Command (cwd notes) | Result |
|---|---|---|
| Analytic induction | `test_3d_ophelie_analytic_induction_simple --dp=0.05` | **passed=1**; `Q_o/c≈24.81` (analytic 25); `phi_eq_res_vol≈3.7e-4`; ~1.5 s |
| Biot circular loop | `test_3d_ophelie_biot_savart_circular_loop_axis` | **passed=1**; `max_rel_err≈1.5e-5` |
| Biot multiloop | `test_3d_ophelie_biot_savart_multiloop_axis` | **passed=1**; `max_rel_err≈4.4e-5` |
| French A_ind diag | `test_3d_ophelie_french_aind_diagnostic --reload-dir=<natural_glass_relax/bin/reload> --ophelie-current-form=edge-flux --ophelie-edge-flux-complex=1` | **passed=1**; n=37428; `phi_eq_res_vol≈1.89e-4`; finite `P_joule≈0.86 W` (glass σ=16, I~1 A demo scale); ~4 s |

Build targets used:

```bash
cd /home/ubuntu/SPHinXsys/build
ninja test_3d_ophelie_analytic_induction_simple \
      test_3d_ophelie_biot_savart_circular_loop_axis \
      test_3d_ophelie_biot_savart_multiloop_axis \
      test_3d_ophelie_french_aind_diagnostic
```

**Which φ path produces finite J/Q?**  
`--ophelie-current-form=edge-flux` + `--ophelie-edge-flux-complex=1` with PCG φ — confirmed on analytic + French reload smokes.

**VTK note:** `french_aind_diagnostic` smoke wrote only `ShapeSPHSystemDomain.vtp` under its `output/` (no field body VTP). Full field VTPs exist on French stirring / frozen_q paths via `BodyStatesRecording` — Case A must add writers for `BCoil*`, `|B|`, and later `J*` / `JouleHeatEdgeReconComplex`.

---

## 5. Known issues / risks (for demo planning)

1. **Dirty worktree** — create a dedicated demo branch; isolate commits from Doc105 hydro WIP.  
2. **Source-field approximation** — coil-only Biot omits full workpiece reaction `B` unless Picard/`A_ind` enabled; README must say so.  
3. **Al skin under-resolved** at plan `dp=4 mm` / `f=10 kHz` — either lower `f` for feasibility, or accept qualitative Joule maps.  
4. **French σ≈16 S/m** ≠ aluminum; do not reuse glass material defaults for welding demos.  
5. **Power scaling / literature_mode** defaults may retarget `I` toward 50 kW glass — disable or set explicit `I_peak` for welding demo.  
6. **φ residual floor** often ~1e-4–2e-4 on French meshes; use documented engineering gates, do not claim “fully converged Krylov”.  
7. GPU driver may be unavailable; still verify demo scripts on a machine with CUDA if ParaView/runtime targets GPU.

---

## 6. Gate 1–2 Case A (done 2026-10-10)

Implemented: `tests/extra_source_and_tests/3d_examples/test_3d_ophelie_induction_demo/`

| Item | Detail |
|---|---|
| Entry | `test_3d_ophelie_induction_demo.cpp` |
| CLI | `--case=a` `--stage=geometry\|biot` `--dp=` `--i-peak=` |
| Geometry | Al cylinder `R=0.02`, `L=0.10`; 5-loop coil `R=0.035`, pitch `0.012`, z-axis |
| Biot | `applyMultiloopFilamentBiotToGlass` + `Combine…` → `BSrc*` / `BSrcMag` |
| Checks | axis filament vs analytic; `I→2I` ratio = 2; finite `|B|` |
| Smoke | `n=1971`, `max\|B\|≈0.037 T`, `passed=1` |

**Gate 3 (next, only when instructed):** edge-flux complex φ → recon `J` / `JouleHeatEdgeReconComplex` with residual gates; else keep B-only for external slides. Do **not** start Case B until Gate 3 decision.

---

## 7. Gate 0 acceptance checklist

| Criterion | Met? |
|---|---|
| Existing OPHELIE targets compile | Yes (ninja targets above) |
| Identified φ path with finite J/Q | Yes — edge-flux complex |
| Repeatable smoke commands + logs | Yes — `gate0_smoke/*.log` |
| No large demo coding yet | Yes |
| French cases not overwritten | Yes |

**Stop here per plan §8.** Next instruction: Gate 1 + Gate 2 only (Case A geometry + Biot–Savart `B_src`).
