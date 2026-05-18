# Matrix-Free A-phi Joule Heating Baseline

This case extends the standalone physical-response baseline by taking the computed Joule heat density and using it as a one-way thermal source.

## Purpose
- reuse the penalty-stabilized matrix-free `A-phi` forced-response route
- compute `E`, `J`, and Joule heat on the thin variable-conductivity band
- feed Joule heat into a simple thermal diffusion update
- check that the resulting temperature field is symmetric and that the mean temperature rise matches the injected average heating rate

## Current mode
- `coulomb_variable_sigma_forced_response`

## Forcing profiles
- default: `EM_APHI_FORCED_PROFILE=sin`
- localized option: `EM_APHI_FORCED_PROFILE=gaussian` with
  - `EM_APHI_FORCED_CENTER_X`
  - `EM_APHI_FORCED_WIDTH`

## Main diagnostics
- `converged`
- `residual_ay_l2`
- `joule_density_l2`
- `joule_density_max`
- `thermal_time`
- `thermal_steps`
- `temperature_average`
- `temperature_min`
- `temperature_max`
- `temperature_delta_average`
- `expected_temperature_delta_average`
- `temperature_delta_average_error`
- `temperature_mirror_symmetry_l2`
- `temperature_delta_mirror_symmetry_l2`

## Exploratory command
```bash
EM_APHI_STAGGERED_SMOKE_CASE=coulomb_variable_sigma_forced_response \
EM_APHI_STAGGERED_SMOKE_DP=0.05 \
EM_APHI_STAGGERED_SMOKE_ENABLE_GAUGE=0 \
EM_APHI_ENABLE_GAUGE_PENALTY=1 \
EM_APHI_STAGGERED_SMOKE_OUTER_ITERS=600 \
/home/yongchuan/sphinxsys/build/tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_matrix_free_joule_heating_baseline/bin/test_3d_em_aphi_matrix_free_joule_heating_baseline
```

## Thermal knobs
- `EM_APHI_FORCED_PROFILE`
- `EM_APHI_FORCED_CENTER_X`
- `EM_APHI_FORCED_WIDTH`
- `EM_APHI_JOULE_THERMAL_DIFFUSIVITY`
- `EM_APHI_JOULE_RHO_CP`
- `EM_APHI_JOULE_INITIAL_TEMPERATURE`
- `EM_APHI_JOULE_END_TIME`
- `EM_APHI_JOULE_MAX_TEMPERATURE_MIRROR`
- `EM_APHI_JOULE_MAX_TEMPERATURE_DELTA_MIRROR`
- `EM_APHI_JOULE_MAX_TEMPERATURE_DELTA_AVG_ERROR`


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
EM_APHI_JOULE_MAX_TEMPERATURE_MIRROR=1e-6 \
EM_APHI_JOULE_MAX_TEMPERATURE_DELTA_MIRROR=1e-6 \
EM_APHI_JOULE_MAX_TEMPERATURE_DELTA_AVG_ERROR=1e-5 \
/home/yongchuan/sphinxsys/build/tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_matrix_free_joule_heating_baseline/bin/test_3d_em_aphi_matrix_free_joule_heating_baseline
```


## Localized hot-spot note
A localized variant is now available with:
```bash
EM_APHI_FORCED_PROFILE=gaussian \
EM_APHI_FORCED_CENTER_X=0.5 \
EM_APHI_FORCED_WIDTH=0.12
```

The first `dp=0.05` localized run stayed converged and produced a smaller but sharper hot spot, with representative diagnostics:
- `joule_density_max=1.78e-02`
- `temperature_delta_average=1.31e-03`
- `temperature_delta_average_error=4.64e-06`
- `temperature_delta_mirror_symmetry_l2=0`
