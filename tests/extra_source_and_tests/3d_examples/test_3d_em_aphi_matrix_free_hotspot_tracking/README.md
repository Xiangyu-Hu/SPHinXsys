# Matrix-Free A-phi Hotspot Tracking

This case extends the localized Gaussian forced-response Joule-heating baseline by checking whether the Joule and temperature hot spots stay near the prescribed forcing center and respond symmetrically when the forcing center is mirrored.

## Purpose
- keep the same thin-band matrix-free electromagnetic and one-way thermal pipeline
- drive the system with a localized Gaussian forcing profile
- verify that the response hot spot is spatially where we expect it to be
- verify that off-center hot-spot drift mirrors correctly between left and right forcing locations

## Main diagnostics
- `forced_profile`
- `forced_center_x`
- `forced_width`
- `source_centroid_x`
- `joule_centroid_x`
- `temperature_delta_centroid_x`
- `source_center_error`
- `joule_center_error`
- `temperature_delta_center_error`

## Exploratory commands
```bash
EM_APHI_STAGGERED_SMOKE_CASE=coulomb_variable_sigma_forced_response \
EM_APHI_FORCED_PROFILE=gaussian \
EM_APHI_FORCED_CENTER_X=0.35 \
EM_APHI_FORCED_WIDTH=0.12 \
EM_APHI_STAGGERED_SMOKE_DP=0.05 \
EM_APHI_STAGGERED_SMOKE_ENABLE_GAUGE=0 \
EM_APHI_ENABLE_GAUGE_PENALTY=1 \
EM_APHI_STAGGERED_SMOKE_OUTER_ITERS=600 \
/home/yongchuan/sphinxsys/build/tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_matrix_free_hotspot_tracking/bin/test_3d_em_aphi_matrix_free_hotspot_tracking
```

```bash
EM_APHI_STAGGERED_SMOKE_CASE=coulomb_variable_sigma_forced_response \
EM_APHI_FORCED_PROFILE=gaussian \
EM_APHI_FORCED_CENTER_X=0.65 \
EM_APHI_FORCED_WIDTH=0.12 \
EM_APHI_STAGGERED_SMOKE_DP=0.05 \
EM_APHI_STAGGERED_SMOKE_ENABLE_GAUGE=0 \
EM_APHI_ENABLE_GAUGE_PENALTY=1 \
EM_APHI_STAGGERED_SMOKE_OUTER_ITERS=600 \
/home/yongchuan/sphinxsys/build/tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_matrix_free_hotspot_tracking/bin/test_3d_em_aphi_matrix_free_hotspot_tracking
```

## Current interpretation
- the left and right source centers are both placed accurately
- Joule and temperature hot spots shift by nearly equal magnitudes in mirrored directions
- this indicates a stable mirrored physical response rather than random positional drift

## Hotspot gate knobs
- `EM_APHI_HOTSPOT_MAX_SOURCE_CENTER_ERROR`
- `EM_APHI_HOTSPOT_MAX_JOULE_CENTER_ERROR`
- `EM_APHI_HOTSPOT_MAX_TEMPERATURE_CENTER_ERROR`
