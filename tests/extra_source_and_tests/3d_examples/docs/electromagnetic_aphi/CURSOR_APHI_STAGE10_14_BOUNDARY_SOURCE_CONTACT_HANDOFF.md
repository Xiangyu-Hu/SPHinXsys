# Stage 10.14-B handoff — boundary / source / contact closure

> User confirmed ChatGPT action list (2026-05-21). Baseline **10.14-A** frozen; this tracks **10.14-B** closure.

## Frozen (do not change)

- Production `lambda_A = off`
- Contact A-penalty research: **InnerOnly** only; projection deferred
- Default outer boundary: **`enlarged_air_domain_padding`** (`boundary_width_scale`), not ghost/mirror
- `passive_air_shell` (legacy `dummy_shell`): diagnostic only; not SphinxSys ghost BC
- `aphi_ghost_buffer_diva`: divA research only — **not** wired into heating solver
- TEAM7 coil: prescribed current RHS, **coil σ=0**, coil Joule excluded from conductor stats

## 10.14-B deliverables

| Item | Status | Notes |
|------|--------|--------|
| P1 three-body bodywise audit | done | `aphi_physical_region_audit_helpers.h`, composition gate |
| P2 boundary strong gate | done | `PassiveAirShell`, padding sensitivity, `passed=1` |
| P3 prescribed source metrics | partial | `source_rhs_l2`, particle counts, coil σ=0 |
| P4 thermal closure | in progress | mapping &lt;1%; source-driven spatial &lt;5% (host adiabatic step for SYCL) |
| P5 interface spike hard @50 | done | bulk reference uses `x > x_interface + band + dp` |
| P6 annular source diagnostic | done | `test_3d_aphi_ck_annular_source_region_diagnostic` (TEAM7 coil slot, σ=0) |
| P7 geometry relaxation scaffold | done | `test_3d_aphi_ck_em_geometry_relaxation_diagnostic` (lattice only) |

## Key files

- `extra_src/shared/electromagnetic_dynamics/test_helpers/aphi_physical_region_audit_helpers.h`
- `extra_src/shared/electromagnetic_dynamics/test_helpers/aphi_em_observable_helpers.h` — physical region + spike
- `extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_source_driven_em_solve_helpers.h` — `runSourceDrivenEmSolveWithLayout`
- `extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_boundary_support_policy_diagnostic_helpers.h`
- `extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_contact_source_driven_heating_baseline_helpers.h`
- `extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_em_joule_thermal_coupling_helpers.h`
- `extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_annular_source_region_diagnostic_helpers.h`
- `extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_em_geometry_relaxation_diagnostic_helpers.h`

## Tests to run (build dir)

```bash
ninja test_3d_aphi_ck_source_driven_em_solve \
  test_3d_aphi_ck_em_joule_thermal_one_way \
  test_3d_aphi_ck_contact_source_driven_heating_baseline \
  test_3d_aphi_ck_three_body_contact_source_driven_baseline \
  test_3d_aphi_ck_boundary_support_policy_diagnostic \
  test_3d_aphi_ck_annular_source_region_diagnostic \
  test_3d_aphi_ck_em_geometry_relaxation_diagnostic
```

## Not in scope (10.14-B)

- Ghost/mirror outer BC
- Maxwell-Ampere MMS (P8) — plan only
- Standard TEAM7 reference comparison
