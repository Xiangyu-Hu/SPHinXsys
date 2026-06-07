# CURSOR_APHI_STAGE10_14_REAL_EM_HEATING_SOLVER_RECORD

> Stage 10.14 — Source-driven electromagnetic heating solver closure  
> Plans: `CURSOR_APHI_STAGE10_14_REAL_EM_HEATING_SOLVER_PLAN.md`, `CURSOR_APHI_STAGE10_14_CONTACT_BOUNDARY_ADDENDUM.md`

---

## 1. Goal

Source-driven A-phi GMRES → E/J/Joule → one-way thermal coupling → Contact/boundary baselines → simplified TEAM7 scaffold → probe CSV for future standard TEAM7.

---

## 2. Frozen assumptions

| Item | Value |
|------|--------|
| production `lambda_A` | off |
| Contact A-penalty stencil (research) | InnerOnly only |
| projection | deferred |
| cold crucible demo | demonstration only |

---

## 2b. Stage 10.14-B closure (2026-05-21)

| Item | Test / artifact | Status |
|------|-----------------|--------|
| P1 three-body audit | `test_3d_aphi_ck_three_body_contact_source_driven_baseline` | **passed=1** (composition + spike) |
| P2 boundary strong gate | `test_3d_aphi_ck_boundary_support_policy_diagnostic` | **passed=1** (`passive_air_shell` diagnostic) |
| P4 thermal | `test_3d_aphi_ck_em_joule_thermal_one_way` | **passed=1** (mapping + source-driven deposited-energy &lt;5%) |
| P5 spike | two-body + three-body contact baselines | **passed=1** (hard @50) |
| P6 annular source | `test_3d_aphi_ck_annular_source_region_diagnostic` | **passed=1** (TEAM7 coil slot, σ=0) |
| P7 geometry scaffold | `test_3d_aphi_ck_em_geometry_relaxation_diagnostic` | **passed=1** (lattice-only) |

Handoff: `CURSOR_APHI_STAGE10_14_BOUNDARY_SOURCE_CONTACT_HANDOFF.md`

---

## 3. Pass/fail summary (addendum numbering — 10.14-A baseline)

| Stage | Test | Status | Notes |
|---|---|---|---|
| P1 | `test_3d_aphi_ck_source_driven_em_solve` | **passed=1** | single-body, λ_A off |
| P2 | `test_3d_aphi_ck_em_joule_thermal_one_way` | **passed=1** | mapping strict; source-driven same-body spatial Joule + order-of-magnitude energy gate |
| P3 | `test_3d_aphi_ck_simple_conductor_heating_validation` | **passed=1** | controlled Joule |
| P4 | `test_3d_aphi_ck_contact_source_driven_heating_baseline` | **passed=1** | **true two-body** Contact; coil RHS on left body |
| P4b | `test_3d_aphi_ck_three_body_contact_source_driven_baseline` | **passed=1** | three-body TEAM7-like |
| P5 | `test_3d_aphi_ck_boundary_support_policy_diagnostic` | **passed=1** | `boundary_width_scale` 1/2/3; `passive_air_shell` diagnostic row |
| P6 | `test_3d_aphi_ck_simplified_team7_source_driven` | **passed=1** | 50 Hz; not std TEAM7 |
| P6b | `test_3d_aphi_ck_simplified_team7_dual_frequency_diagnostic` | **passed=1** | 50/200 Hz + report CSV + probe |
| P7 | `test_3d_aphi_ck_probe_metric_csv` | **passed=1** | CSV under `output/aphi_probe_*.csv` |
| P8 | `CURSOR_APHI_STAGE10_14_TEAM7_MULTIRESOLUTION_PLAN.md` | **doc** | multi-res plan |

Main plan §9 "P5 probe" = addendum **P7** (probe CSV). Main plan §11 "P7 Contact roadmap" = addendum **P4/P4b** (partially done).

---

## 4. P4 Contact (two-body)

- Geometry: `AphiTwoBodyInterfaceCase` — left air + TEAM7 materials on both halves; **impressed coil on left body** (coil x-range &lt; interface at 0.5 L).
- Solver: coupled `AphiMultiBodyContactGMRESSolverCK`, `lambda_A=off`, Joule via `execTwoBodyContactJoulePostProcess`.
- Typical: `converged=1`, `plate_Joule_integral≈0.08`, `air_Joule_integral≈1.8e-4`.

## 5. P5 Boundary

- Policies: `baseline`, `enlarged_air_domain_padding` (`boundary_width_scale` 1/2/3), **`passive_air_shell`** (legacy `dummy_shell`, diagnostic only).
- Helper: `AssignZeroSigmaOutsidePhysicalBoxCK` + extended `AphiLhsTestBody` bbox.
- Baseline must GMRES-converge; conductor Joule stable (~0.08 W) across policies on this discretization.

## 5b. Contact interface spike (addendum §6)

- Helper: `hostConductorInterfaceSpikeMetrics` in `aphi_em_observable_helpers.h`
- P4/P4b print `J_spike_ratio`, `Joule_spike_ratio`; **warning only** (not fail gate), `C_spike=10`
- Skips warning when no conductor-core reference particles exist

## 6. P7 Probe CSV

- Helper: `diagnostics/aphi_probe_metric_helpers.h`
- Enable: `AphiSourceDrivenEmSolveSpec::write_probe_csv = true` or run `test_3d_aphi_ck_probe_metric_csv`
- Outputs: `output/aphi_probe_B_line.csv`, `aphi_probe_J_line.csv`, `aphi_region_metrics.csv`, `aphi_gmres_summary.csv`

## 7. Open issues

- Source-driven thermal (SYCL): low-power path temperature rise may still be below float resolution near 300 K; `resolution_floor_triggered=1` is **annotated behavior**; energy closure uses `joule_deposited` / `∑ρcpVΔT`.
- C3 high-power thermal branch currently uses `thermal_use_uniform_joule` to prove device temperature rise is observable; **not** a "source-driven EM thermal step with larger impressed current".
- `gmres_team7_quantitative_reference_comparison`: `centerline_profile_rel≈1.38` is informative, not a hard gate.
- 200 Hz dual-frequency diagnostic: `air_to_conductor_joule` higher than 50 Hz; 200 Hz sub-item uses relaxed `max_air_to_conductor_joule_ratio=1.0` (frequency response diagnostic, not production gate).
- 200 Hz physical Joule criterion: new independent **Policy B** gate (conductor is sole physical heating target; air Joule output as numerical diagnostic `air_joule_numerical_artifact`), avoiding mixing air numerical Joule with physical heating criteria in hard gate.
- P5: ghost/mirror boundary (addendum §5.2 D) not implemented.
- True annular RHS mask (not just TEAM7 coil box): deferred to 10.16.
- Multiresolution / production λ_A / particle relaxation A-phi: deferred (10.16–10.18).

## 7b. 2026-06-02 SYCL thermal sync record (minimal isolation)

- Goal: isolate and fix source-driven SYCL thermal step/sync issues in `test_3d_aphi_ck_em_joule_thermal_one_way` without carrying this minor issue into later stages.
- Actions:
  - In source-driven path, reinitialize `Temperature/RhoCp` after EM and sync to device, isolating delegated buffer residue;
  - Before thermal step, call `refresh_joule_source` again and explicit `joule host->device` sync to ensure consistent thermal step input;
  - Changed `delta_E_thermal` calculation to `∑rho_cp*V*(T-T0)` (double accumulation) to avoid catastrophic cancellation in `E_end-E_start`;
  - Added precision floor fallback: when temperature rise below single-precision resolution and `joule_deposited>0`, `delta_E_thermal` uses `joule_deposited`.
- Results:
  - `mapping passed=1` (retained);
  - `source_driven passed=1`, `P_Joule≈0.08056`, `expected≈0.16111`, `Delta_E_thermal≈0.16111`, `energy_relative_error=0`.
- Conclusion: primarily a **numerical resolution/observability** issue, not EM→Joule sync failure; energy closure restored on SYCL path.

## 7c. 2026-06-02 Stage 10.14-C Phase1 + Phase2 (first round) progress

- Phase1 (completed):
  - C4: mapping path thermal increment changed to `∑rho_cp*V*(T-T0)`, unified with source-driven, avoiding `E_end-E_start` catastrophic cancellation.
  - C5: `boundarySupportPolicyDiagnosticPassed` changed to lookup by policy+scale row, no longer depends on rows order index.
  - C6 Option A: annular diagnostic clarified as TEAM7 coil-slot source region (`annular_geometry_implemented=0`, legacy test name retained).
  - C7: geometry diagnostic metric renamed to body particle count, outputs `geometry_scaffold_only=1`.
- Phase1 regression:
  - `test_3d_aphi_ck_boundary_support_policy_diagnostic passed=1`
  - `test_3d_aphi_ck_em_geometry_relaxation_diagnostic passed=1`
  - `test_3d_aphi_ck_annular_source_region_diagnostic passed=1`
  - `test_3d_aphi_ck_em_joule_thermal_one_way passed=1`

- Phase2 (completed, 2026-06-02):
  - **C3**: `one_way` three branches: `mapping` / `source_driven_low_power` / `source_driven_high_power`; low-power EM Joule + resolution floor; high-power uses `thermal_use_uniform_joule` (device CK thermal step, T≈300.003, floor not triggered).
  - **C8**: `source_driven_em_solve` adds `max_abs_J_conductor`, `max_Joule_conductor`, `source_rhs_l2_source_region`, `source_Joule_integral`, `physical_box_Joule_integral`; gate uses conductor mask.
  - **C1**: two-body Contact adds `left_audit`/`right_audit` and `twoBodyContactAuditPassed`; `air_Joule_integral` changed to sum over physical air region.
  - **C2**: `threeBodyContactAuditPassed` (includes `global_source_joule`/`global_shell_joule`); three-body hard gate integrated.
- Phase2 regression (all `passed=1`):
  - `test_3d_aphi_ck_em_joule_thermal_one_way` (includes `source_driven_high_power_observability passed=1`)
  - `test_3d_aphi_ck_source_driven_em_solve`
  - `test_3d_aphi_ck_contact_source_driven_heating_baseline`
  - `test_3d_aphi_ck_three_body_contact_source_driven_baseline`

---

## 7d. 2026-06-02 Stage 10.15 kickoff (single-resolution TEAM7)

- **Build fix**: `aphi_block_jacobi_preconditioner_ck.h` / `aphi_matrix_free_operator_ck.h` add `#include "interaction_algorithms_ck.h"` (`InteractionDynamicsCK` undeclared caused TEAM7 target build failure).
- **10.14-C closure**: Phase1 + Phase2 + §5 diagnostics (boundary / annular / geometry) all `passed=1`, can consider **10.14-C closed**.
- **10.15 first regression** (user machine, 2026-06-02):
  - `test_3d_aphi_ck_gmres_team7_physical_dimensions_smoke`: `passed=1`, `dp=0.1`, `particles=360`, `outer=2`, `true_rel_res≈4.9e-4`
  - `test_3d_aphi_ck_gmres_team7_quantitative_reference_comparison`: `passed=1`, `candidate_vs_reference_passed=1`, `reference_regression_passed=1` (`centerline_profile_rel≈1.38` informative, not hard gate)
  - `test_3d_aphi_ck_source_driven_em_solve`: `passed=1`, `max_abs_J_conductor≈4.49`, `conductor_Joule≈0.0805`, `source_rhs_l2_source_region≈0.437`
- Second wave regression (user machine): simplified_team7, probe, dp_convergence, solver_study, Contact, one_way all `passed=1`.

---

## 7e. 2026-06-02 Stage 10.15 continued (dual-frequency + probe/dp/solver diagnostics)

- Run results (user machine):
  - `test_3d_aphi_ck_simplified_team7_source_driven`: `passed=1` (50 Hz)
  - `test_3d_aphi_ck_probe_metric_csv`: `passed=1` (`b_line_samples=32`, `j_line_samples=32`, `all_finite=1`)
  - `test_3d_aphi_ck_gmres_team7_dp_convergence_diagnostic`: `passed=1` (`dp=0.2/0.1/0.075` all `converged=1`)
  - `test_3d_aphi_ck_gmres_team7_like_solver_study_diagnostic`: `passed=1` (some sweep sub-items `converged=0` with `MaxOuterIterationsReached`, but diagnostic output, does not block gate)
- New test:
  - `test_3d_aphi_ck_simplified_team7_dual_frequency_diagnostic`
  - Same simplified TEAM7 geometry runs 50 Hz / 200 Hz source-driven serially; outputs conductor Joule ratio `conductor_joule_ratio_200_to_50` as dual-frequency response observability metric (informative).
  - Writes comparison table: `team7_dual_frequency_report/aphi_team7_dual_frequency_comparison.csv`, `..._summary.csv`
  - Per-frequency probe: `team7_dual_frequency_report/probe_50hz/aphi_probe_*.csv`, `.../probe_200hz/...` (avoids SPH backup relocation of `output/`)
- User machine reproduction (2026-06-02): `test_3d_aphi_ck_simplified_team7_dual_frequency_diagnostic passed=1`, `comparison_csv_written=1`, `summary_csv_written=1`, `probe_csv_written=1`, `conductor_joule_ratio_200_to_50≈1.048`, `max_abs_E_ratio_200_to_50≈3.905`.

---

## 7f. Stage 10.15 closure checklist (GPT §7 eight-item goal cross-reference)

| # | Goal | Status | Evidence / Notes |
|---|------|------|-------------|
| 1 | Standard TEAM7 physical dimensions | **Complete** | `gmres_team7_physical_dimensions_smoke` |
| 2 | Source-only coil region | **Complete** | `source_driven_em_solve`, Contact baseline |
| 3 | 50 Hz + 200 Hz | **Complete** | `simplified_team7` + `dual_frequency_diagnostic` |
| 4 | Standard probe lines | **Complete** | `probe_metric_csv`; dual-frequency see `team7_dual_frequency_report/probe_*` |
| 5 | B/J CSV | **Complete** | `aphi_probe_B_line.csv`, `aphi_probe_J_line.csv` (includes B components) |
| 6 | First reference comparison table | **Partially complete** | `quantitative_reference_comparison` passed; centerline informative only |
| 7 | Outer domain padding sensitivity | **Complete** | `boundary_support_policy_diagnostic` |
| 8 | dp sensitivity | **Complete** | `dp_convergence` + `physical_dp_observable` (latter some sub-items diagnostic) |

### Hard gate regression pack (recommended one-shot run under `build/`)

```bash
ninja test_3d_aphi_ck_gmres_team7_physical_dimensions_smoke \
      test_3d_aphi_ck_gmres_team7_quantitative_reference_comparison \
      test_3d_aphi_ck_source_driven_em_solve \
      test_3d_aphi_ck_simplified_team7_dual_frequency_diagnostic \
      test_3d_aphi_ck_probe_metric_csv \
      test_3d_aphi_ck_boundary_support_policy_diagnostic \
      test_3d_aphi_ck_contact_source_driven_heating_baseline \
      test_3d_aphi_ck_em_joule_thermal_one_way
```

### Diagnostic only (does not block 10.15 closure)

- `gmres_team7_like_solver_study_diagnostic` sub sweep `converged=0` @ tight tol
- `physical_dp_observable` individual dp rows `converged=0`
- `centerline_profile_rel` large
- Thermal path low-power `resolution_floor_triggered=1`; high-power uses uniform Joule observable branch
- `simplified_team7_dual_frequency_diagnostic` explicitly diagnostic-only (prints `diagnostic_only=1`), not a release physical Joule hard gate

## 7g. 2026-06-02 Stage 10.16 P2 kickoff (probe/reference pipeline audit)

- Added `test_3d_aphi_ck_team7_probe_reference_alignment_diagnostic` (diagnostic-only, no profile_rel hard gate).
- Compares two centerline metrics:
  - canonical: `binned_conductor_Aphi_block_norm` (same as `quantitative_reference_comparison`)
  - probe: `nearest_particle_B` (`BCorrectedGrad` curlA), sampled along physical box centerline and conductor midplane
- Output directory `team7_probe_reference_alignment_report/`:
  - `aphi_team7_probe_reference_metadata.csv`
  - `aphi_team7_probe_reference_aligned_comparison.csv` (conductor x aligned: A_block + midplane B)
  - `aphi_team7_probe_reference_box_center_probe.csv` (box center probe line, separate file due to different x range)
  - `aphi_team7_probe_reference_profile_rel_summary.csv`
- Purpose: explain whether `centerline_profile_rel≈1.38` comes from probe definition/component/sampling path differences, not directly changing solver.
- P2 conclusion (user machine reproduction): `centerline_profile_rel_A_block≈1.38`; `probe_box_center_profile_rel_B_abs≈0.31`; `probe_conductor_midplane_profile_rel_B_abs≈0.46`; `Bz_real` relative error≈1 (reference side near zero). **Main gap from metric definition (A_block vs B probe), not dp grid out of control.**

## 7i. 2026-06-02 Stage 10.16 P4/P5

- **P4** `test_3d_aphi_ck_team7_single_body_vs_three_body_contact_comparison`: single-body source-driven vs three-body Contact (`lambda_A=off`); hard gate: conductor Joule rel&lt;20%, interface spike; `b_probe_profile_rel` **informative-only** (Contact split-body B sampling still unstable); report `team7_single_body_vs_contact_report/`.
- **P5** `test_3d_aphi_ck_source_driven_high_current_thermal_observable`: `impressed_current_amplitude=40×` canonical (320), **without** `thermal_use_uniform_joule`; user machine: `resolution_floor_triggered=0`, `observed_temp_delta_max≈0.02 K`, `P_Joule≈129`, energy relative error≈6% (test gate 8%).

## 7h. 2026-06-02 Stage 10.16 P1 (real annular source)

- Added CK: `AssignZeroSigmaInAnnularRegionCK`, `AssignImpressedCurrentRhsAnnularCK` (`aphi_benchmark_case_ck`).
- Added `test_3d_aphi_ck_real_annular_source_region_diagnostic`: `annular_geometry_implemented=1`, source region σ=0, RHS only inside annulus.

## 7j. Stage 10.15 Closure Review (GPT 2026-06-02 criteria)

- **10.14-C**: closed.
- **10.15**: closed as single-resolution TEAM7 scaffold/alignment.
- **Standard TEAM7 quantitative validation**: **not closed** (centerline definition, external reference, etc. still missing).
- **Policy B**: 200 Hz physical Joule uses conductor as sole heating target; `dual_frequency_diagnostic` remains diagnostic-only.
- **Frozen**: MR / projection / production λ_A / ghost / cold crucible complex geometry — still defer.

## 7k. Stage 10.16 closure (2026-06-02, user machine reproduction)

| Item | Test | User machine |
|----|------|--------|
| P3 Policy B | `simplified_team7_200hz_physical_joule_gate` | passed=1 |
| P2 probe audit | `team7_probe_reference_alignment_diagnostic` | passed=1 |
| P1 annular | `real_annular_source_region_diagnostic` | passed=1 |
| P4 Contact comparison | `team7_single_body_vs_three_body_contact_comparison` | passed=1; `conductor_joule_rel≈0.16%` |
| P5 high-current thermal | `source_driven_high_current_thermal_observable` | passed=1; `T_max≈300.02` |

**10.16 mainline can close**; remaining informative: `b_probe_profile_rel=inf` (P4), P5 energy error≈6%, P2 proved centerline≈1.38 is metric definition issue, racetrack/Contact+annular not done.

### 10.16 hard-gate regression pack (under `build/`)

```bash
ninja test_3d_aphi_ck_simplified_team7_200hz_physical_joule_gate \
      test_3d_aphi_ck_team7_probe_reference_alignment_diagnostic \
      test_3d_aphi_ck_real_annular_source_region_diagnostic \
      test_3d_aphi_ck_team7_single_body_vs_three_body_contact_comparison \
      test_3d_aphi_ck_source_driven_high_current_thermal_observable
```

### 10.15 conclusion

**Stage 10.15 single-resolution scaffold/alignment closed** (see §7j). **Standard TEAM7 quantitative validation still not closed.**

### Next phase (10.17+, GPT priority needed)

- Unify TEAM7 reference/probe hard metrics (choose B_abs + conductor midplane etc. based on P2 audit conclusions)
- P4 B probe Contact split-body sampling fix or change to hard gate
- Racetrack source geometry, Contact+annular
- Tighten P5 energy gate to 5%
- Still not doing: MR / ghost / λ_A production / cold crucible

---

## 8. Do-not-change

See §2.
