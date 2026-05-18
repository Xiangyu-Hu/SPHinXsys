# Matrix-Free A-phi Physical Response

This case is the simple physical-response companion to the manufactured `staggered_smoke` baselines.

## Purpose
- keep the same thin-band matrix-free `A-phi` solver route
- use an explicit forcing in the `A_y` equation instead of a back-computed manufactured exact solution
- check physically interpretable response quantities before moving to heavier geometry or heating cases

## Current mode
- `coulomb_variable_sigma_forced_response`

## Forcing profiles
- default: `EM_APHI_FORCED_PROFILE=sin`
- localized option: `EM_APHI_FORCED_PROFILE=gaussian` with
  - `EM_APHI_FORCED_CENTER_X`
  - `EM_APHI_FORCED_WIDTH`

## Meaning of the name
- `coulomb`: the setup stays on the Coulomb-compatible `div(A)`-controlled solver route
- `variable_sigma`: conductivity varies across the band as `sigma(y)`
- `forced_response`: the response is driven by an explicit forcing term, not by a manufactured exact field

## What to check
Because this mode has no strong manufactured reference, the primary diagnostics are:
- `converged`
- `residual_ay_l2`
- `electric_field_l2`
- `current_density_l2`
- `joule_density_l2`
- `joule_density_max`
- `ay_mirror_symmetry_l2`
- `phi_mirror_antisymmetry_l2`
- `phi_centered_mirror_symmetry_l2`
- `phi_centered_mirror_antisymmetry_l2`
- `joule_mirror_symmetry_l2`

The summary prints `has_discrete_reference=0` to make it explicit that reference-based L2 errors are not the governing checks here. For `phi`, prefer the centered parity diagnostics first so the reference-index offset does not contaminate the interpretation.

## Exploratory command
```bash
EM_APHI_STAGGERED_SMOKE_CASE=coulomb_variable_sigma_forced_response \
EM_APHI_STAGGERED_SMOKE_DP=0.05 \
EM_APHI_STAGGERED_SMOKE_ENABLE_GAUGE=0 \
EM_APHI_ENABLE_GAUGE_PENALTY=1 \
EM_APHI_STAGGERED_SMOKE_OUTER_ITERS=600 \
/home/yongchuan/sphinxsys/build/tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_matrix_free_physical_response/bin/test_3d_em_aphi_matrix_free_physical_response
```

## Symmetry-oriented gate variables
- `EM_APHI_STAGGERED_SMOKE_MAX_AY_MIRROR`
- `EM_APHI_STAGGERED_SMOKE_MAX_PHI_MIRROR`
- `EM_APHI_STAGGERED_SMOKE_MAX_PHI_ANTIMIRROR`
- `EM_APHI_STAGGERED_SMOKE_MAX_PHI_CENTERED_MIRROR`
- `EM_APHI_STAGGERED_SMOKE_MAX_PHI_CENTERED_ANTIMIRROR`
- `EM_APHI_STAGGERED_SMOKE_MAX_JOULE_MIRROR`


## Recommended first validation gate
```bash
EM_APHI_STAGGERED_SMOKE_CASE=coulomb_variable_sigma_forced_response \
EM_APHI_STAGGERED_SMOKE_DP=0.05 \
EM_APHI_STAGGERED_SMOKE_ENABLE_GAUGE=0 \
EM_APHI_ENABLE_GAUGE_PENALTY=1 \
EM_APHI_STAGGERED_SMOKE_OUTER_ITERS=600 \
EM_APHI_STAGGERED_SMOKE_REQUIRE_CONVERGED=1 \
EM_APHI_STAGGERED_SMOKE_MAX_RESIDUAL_AY=1.2e-5 \
EM_APHI_STAGGERED_SMOKE_MAX_AY_MIRROR=1e-6 \
EM_APHI_STAGGERED_SMOKE_MAX_JOULE_MIRROR=1e-6 \
EM_APHI_STAGGERED_SMOKE_MAX_PHI_CENTERED_ANTIMIRROR=1e-6 \
/home/yongchuan/sphinxsys/build/tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_matrix_free_physical_response/bin/test_3d_em_aphi_matrix_free_physical_response
```

This gate intentionally avoids raw `phi` parity checks and instead uses the centered anti-mirror diagnostic, because the scalar-potential reference constraint introduces a removable constant offset.
