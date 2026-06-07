# Stage 10.14 Round 2 — Cursor Handoff for ChatGPT (Optional)

> **Conclusion**: P1–P7 (supplement numbering) baselines have run; Cursor can continue P2 energy closure, boundary dummy_shell, Contact spike warnings, **without** pausing for ChatGPT on those items.  
> Upload this file + listed files **only** when aligning with ChatGPT on **standard TEAM7 reference data / production acceptance thresholds**.

---

## 1. Two P Numbering Schemes

| Topic | Main plan `CURSOR_APHI_STAGE10_14_REAL_EM_HEATING_SOLVER_PLAN.md` | Addendum `CURSOR_APHI_STAGE10_14_CONTACT_BOUNDARY_ADDENDUM.md` |
|------|----------------------------------------------------------------|-------------------------------------------------------------|
| Contact baseline | §11 P7 roadmap | **P4 / P4b** (implemented) |
| Boundary diagnostics | (not listed separately) | **P5** (implemented padding scale) |
| Probe CSV | §9 **P5** | **P7** (implemented) |
| Multi-resolution docs | §10 **P6** | **P8** (existing `CURSOR_APHI_STAGE10_14_TEAM7_MULTIRESOLUTION_PLAN.md`) |

---

## 2. Current passed=1 (Local)

| Test | Purpose |
|------|------|
| `test_3d_aphi_ck_source_driven_em_solve` | P1 single-body source-driven |
| `test_3d_aphi_ck_em_joule_thermal_one_way` | P2 thermal coupling (see limitations) |
| `test_3d_aphi_ck_simple_conductor_heating_validation` | P3 |
| `test_3d_aphi_ck_contact_source_driven_heating_baseline` | P4 **true two-body** Contact |
| `test_3d_aphi_ck_three_body_contact_source_driven_baseline` | P4b three-body |
| `test_3d_aphi_ck_boundary_support_policy_diagnostic` | P5 boundary padding |
| `test_3d_aphi_ck_simplified_team7_source_driven` | P6 simplified TEAM7 @ 50 Hz |
| `test_3d_aphi_ck_probe_metric_csv` | P7 CSV probe |

Build directory: `/home/yongchuan/sphinxsys/build`  
Run example: `./tests/extra_source_and_tests/3d_examples/<test_name>/bin/<test_name>`

---

## 3. Cursor Can Proceed Independently (No ChatGPT Needed)

1. **P2**: Strict energy closure after EM/thermal split by body on SYCL (`aphi_em_joule_thermal_coupling_helpers.h`).
2. **P5 extension**: `dummy_shell` / ghost boundary strategy (`aphi_boundary_support_policy_diagnostic_helpers.h`).
3. **Contact addendum §6**: interface/core J, Joule spike **warnings** (not fail gate), wired into P4/P4b.
4. **Probes**: attach `write_probe_csv=true` on P6 or Contact cases.

---

## 4. Questions for ChatGPT Discussion (Needs Reference Data / Product Judgment)

1. **Standard TEAM7**: Which literature/table provides **reference coordinates and reference B, J curves** for probe lines/points? We only have simplified geometry + CSV infrastructure, no quantitative comparison gate.
2. **P2 acceptance**: For source-driven thermal coupling, is "conductor lumped uniform Joule power + ΔE&gt;0" acceptable, or must we enforce **&lt;5% energy_relative_error** and spatial Joule distribution?
3. **Boundary**: `boundary_width_scale` is stable; is **dummy_shell** still required to claim TEAM7-ready? How to set shell thickness/σ?
4. **Contact**: Two-body baseline converges; should three-body **A-penalty research** still rank before standard TEAM7? Should `C_spike=10` be a hard gate?
5. **Stage 10.14 closure criteria**: Is minimum set still addendum §10 `P1+P2+P3+P4+P5`, or must we add P7 probe comparison + strict P2?

---

## 5. Suggested Upload File List (For ChatGPT)

### 5.1 Plans and Records (Required Reading)

| File | Path |
|------|------|
| Main plan | [`CURSOR_APHI_STAGE10_14_REAL_EM_HEATING_SOLVER_PLAN.md`](../../../../CURSOR_APHI_STAGE10_14_REAL_EM_HEATING_SOLVER_PLAN.md) |
| Contact/boundary addendum | [`CURSOR_APHI_STAGE10_14_CONTACT_BOUNDARY_ADDENDUM.md`](../../../../CURSOR_APHI_STAGE10_14_CONTACT_BOUNDARY_ADDENDUM.md) |
| Stage record | [`CURSOR_APHI_STAGE10_14_REAL_EM_HEATING_SOLVER_RECORD.md`](CURSOR_APHI_STAGE10_14_REAL_EM_HEATING_SOLVER_RECORD.md) |
| Multi-resolution plan | [`CURSOR_APHI_STAGE10_14_TEAM7_MULTIRESOLUTION_PLAN.md`](CURSOR_APHI_STAGE10_14_TEAM7_MULTIRESOLUTION_PLAN.md) |
| This handoff | [`CURSOR_APHI_STAGE10_14_CHATGPT_HANDOFF_ROUND2.md`](CURSOR_APHI_STAGE10_14_CHATGPT_HANDOFF_ROUND2.md) |

### 5.2 Core Helpers (As Needed)

| File | Path |
|------|------|
| Source-driven solve | [`extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_source_driven_em_solve_helpers.h`](../extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_source_driven_em_solve_helpers.h) |
| Thermal coupling | [`extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_em_joule_thermal_coupling_helpers.h`](../extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_em_joule_thermal_coupling_helpers.h) |
| Contact source-driven | [`extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_contact_source_driven_heating_baseline_helpers.h`](../extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_contact_source_driven_heating_baseline_helpers.h) |
| Boundary diagnostics | [`extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_boundary_support_policy_diagnostic_helpers.h`](../extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_boundary_support_policy_diagnostic_helpers.h) |
| Probe CSV | [`extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_probe_metric_helpers.h`](../extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_probe_metric_helpers.h) |
| Two-body Contact infrastructure | [`extra_src/shared/electromagnetic_dynamics/test_helpers/aphi_contact_gmres_test_helpers.h`](../extra_src/shared/electromagnetic_dynamics/test_helpers/aphi_contact_gmres_test_helpers.h) |
| Three-body TEAM7 Contact | [`extra_src/shared/electromagnetic_dynamics/test_helpers/aphi_team7_contact_test_helpers.h`](../extra_src/shared/electromagnetic_dynamics/test_helpers/aphi_team7_contact_test_helpers.h) |

### 5.3 Test Entry Points (Each Directory Has `CMakeLists.txt`)

| Test | Path |
|------|------|
| P1 | [`test_3d_aphi_ck_source_driven_em_solve/`](test_3d_aphi_ck_source_driven_em_solve/) |
| P2 | [`test_3d_aphi_ck_em_joule_thermal_one_way/`](test_3d_aphi_ck_em_joule_thermal_one_way/) |
| P4 | [`test_3d_aphi_ck_contact_source_driven_heating_baseline/`](test_3d_aphi_ck_contact_source_driven_heating_baseline/) |
| P4b | [`test_3d_aphi_ck_three_body_contact_source_driven_baseline/`](test_3d_aphi_ck_three_body_contact_source_driven_baseline/) |
| P5 | [`test_3d_aphi_ck_boundary_support_policy_diagnostic/`](test_3d_aphi_ck_boundary_support_policy_diagnostic/) |
| P6 | [`test_3d_aphi_ck_simplified_team7_source_driven/`](test_3d_aphi_ck_simplified_team7_source_driven/) |
| P7 | [`test_3d_aphi_ck_probe_metric_csv/`](test_3d_aphi_ck_probe_metric_csv/) |

### 5.4 Sample Output (P7)

After running `test_3d_aphi_ck_probe_metric_csv`, under `output/`:

- `aphi_probe_B_line.csv`
- `aphi_probe_J_line.csv`
- `aphi_region_metrics.csv`
- `aphi_gmres_summary.csv`

---

## 6. Known Limitations (Honest)

- **P4**: Coil in left half-domain, excitation on `left_body` (geometry constraint), not "right-half conductor + coil in same body".
- **P2 source-driven**: EM and thermal are separate `SPHBody` on SYCL; energy gate is temporarily `ΔE&gt;0`, not 5% closure.
- **P5**: `boundary_width_scale` only, no `dummy_shell`.
- **projection / production λ_A**: not enabled (frozen).

---

## 7. One-Line Instruction for ChatGPT

```text
Stage 10.14 addendum P1–P7 baselines passed on feature/electromagnetic.
Please answer: (1) standard TEAM7 reference probe specs and pass/fail thresholds,
(2) whether P2 may use lumped conductor Joule for thermal validation,
(3) whether dummy_shell is required before TEAM7 scaffold sign-off,
(4) recommended order for three-body A-penalty vs standard TEAM7.
Do not request full test source dumps; use RECORD + helpers + this handoff.
```
