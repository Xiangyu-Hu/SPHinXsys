# CURSOR_APHI_STAGE10_16_PLAN — Status (2026-06-02)

> Parent: `CURSOR_APHI_STAGE10_16_NEXT_PLAN_AFTER_1014C_1015_REVIEW.md`  
> Record: `CURSOR_APHI_STAGE10_14_REAL_EM_HEATING_SOLVER_RECORD.md` §7g–7k

## Executive status

**Stage 10.16 mainline: CLOSED** (user machine regression confirmed 2026-06-02).

Do **not** claim standard TEAM7 quantitative validation complete.

## Completed tasks

| ID | Task | Test / artifact | Status |
|----|------|-----------------|--------|
| P3 | 200 Hz Policy B (conductor-only physical Joule) | `test_3d_aphi_ck_simplified_team7_200hz_physical_joule_gate` | passed=1 |
| P2 | probe/reference pipeline audit | `test_3d_aphi_ck_team7_probe_reference_alignment_diagnostic` | passed=1 |
| P1 | real annular source (Option A) | `test_3d_aphi_ck_real_annular_source_region_diagnostic` | passed=1 |
| P4 | single-body vs three-body Contact | `test_3d_aphi_ck_team7_single_body_vs_three_body_contact_comparison` | passed=1 |
| P5 | high-current source-driven thermal | `test_3d_aphi_ck_source_driven_high_current_thermal_observable` | passed=1 |

## Key user-machine numbers

- P2: `centerline_profile_rel_A_block≈1.38`; `probe_box_center_B_abs≈0.31` → metric definition gap, not dp blow-up.
- P4: `conductor_joule_rel≈0.16%`; `b_probe_profile_rel=inf` (informative only).
- P5: `impressed_current=320`, `T_max≈300.02 K`, `resolution_floor_triggered=0`, `energy_rel≈6%` (gate 8%).

## Known gaps (informative / follow-up)

1. **centerline / TEAM7 validation metric** — adopt B-based probe after P2; do not use A-block norm alone.
2. **P4 B probe** — Contact split-body sampling returns `inf`; conductor Joule agreement is good.
3. **P5 energy** — 6% vs 5% target; tighten when SYCL thermal path stable.
4. **racetrack geometry** — not implemented (annular only).
5. **Contact + annular** — not implemented.

## Frozen (unchanged)

- No multiresolution, projection, production lambda_A, ghost/mirror EM BC, full cold-crucible geometry.

## Hard-gate regression (build/)

```bash
ninja test_3d_aphi_ck_simplified_team7_200hz_physical_joule_gate \
      test_3d_aphi_ck_team7_probe_reference_alignment_diagnostic \
      test_3d_aphi_ck_real_annular_source_region_diagnostic \
      test_3d_aphi_ck_team7_single_body_vs_three_body_contact_comparison \
      test_3d_aphi_ck_source_driven_high_current_thermal_observable
```

## Recommended GPT discussion topics (next milestone)

1. Formally close 10.16 and name **10.17** scope.
2. Choose **official TEAM7 comparison metric** (probe line, component, normalization) from P2 findings.
3. Whether P4 `b_probe` should become hard gate after Contact B sampling fix.
4. Whether P5 energy tolerance stays 8% or must reach 5% before closing thermal story.
5. Priority: racetrack vs Contact+annular vs external FEM reference import.
