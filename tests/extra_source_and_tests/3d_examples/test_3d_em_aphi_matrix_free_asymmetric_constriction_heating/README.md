# Matrix-Free A-phi Asymmetric Constriction Heating

This case extends the one-way Joule-heating baseline to a one-sided narrowed conductor geometry.

## Purpose
- keep the same matrix-free electromagnetic and one-way thermal chain
- use a symmetric forcing profile together with an asymmetric conductor shape
- check whether Joule heating and temperature shift toward the narrower side purely because of geometry
- distinguish absolute centroid bias from response-relative drift by reporting source-to-response shifts directly

## Main diagnostics
- `asym_split_x`
- `asym_narrow_height`
- `source_centroid_x`
- `joule_centroid_x`
- `temperature_delta_centroid_x`
- `source_centroid_shift`
- `joule_centroid_shift`
- `temperature_delta_centroid_shift`
- `joule_relative_to_source_shift`
- `temperature_delta_relative_to_source_shift`
- `narrow_joule_average`
- `wide_joule_average`
- `narrow_temperature_delta_average`
- `wide_temperature_delta_average`
- `joule_narrow_to_wide_ratio`
- `temperature_delta_narrow_to_wide_ratio`

## Exploratory command
```bash
EM_APHI_STAGGERED_SMOKE_CASE=coulomb_variable_sigma_forced_response \
EM_APHI_FORCED_PROFILE=sin \
EM_APHI_STAGGERED_SMOKE_DP=0.05 \
EM_APHI_STAGGERED_SMOKE_ENABLE_GAUGE=0 \
EM_APHI_ENABLE_GAUGE_PENALTY=1 \
EM_APHI_STAGGERED_SMOKE_OUTER_ITERS=600 \
/home/yongchuan/sphinxsys/build/tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_matrix_free_asymmetric_constriction_heating/bin/test_3d_em_aphi_matrix_free_asymmetric_constriction_heating
```

## Centered-Gaussian follow-up
Use a centered Gaussian to reduce source-shape bias and focus more directly on geometry-induced response drift:

```bash
EM_APHI_STAGGERED_SMOKE_CASE=coulomb_variable_sigma_forced_response \
EM_APHI_FORCED_PROFILE=gaussian \
EM_APHI_FORCED_CENTER_X=0.5 \
EM_APHI_FORCED_WIDTH=0.12 \
EM_APHI_ASYM_NARROW_HEIGHT=0.08 \
EM_APHI_STAGGERED_SMOKE_DP=0.025 \
EM_APHI_STAGGERED_SMOKE_ENABLE_GAUGE=0 \
EM_APHI_ENABLE_GAUGE_PENALTY=1 \
EM_APHI_STAGGERED_SMOKE_OUTER_ITERS=800 \
/home/yongchuan/sphinxsys/build/tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_matrix_free_asymmetric_constriction_heating/bin/test_3d_em_aphi_matrix_free_asymmetric_constriction_heating
```

Interpretation guide:
- `source_centroid_shift` close to zero means the forcing itself is centered
- `joule_relative_to_source_shift` and `temperature_delta_relative_to_source_shift` then measure geometry-induced drift more cleanly
- `joule_narrow_to_wide_ratio > 1` and `temperature_delta_narrow_to_wide_ratio > 1` indicate local enhancement in the narrow side

## Optional geometry knobs
- `EM_APHI_ASYM_SPLIT_X`
- `EM_APHI_ASYM_NARROW_HEIGHT`

## Optional validation knobs
- `EM_APHI_ASYM_MIN_JOULE_RATIO`
- `EM_APHI_ASYM_MIN_TEMPERATURE_RATIO`
- `EM_APHI_ASYM_MIN_JOULE_SHIFT`
- `EM_APHI_ASYM_MIN_TEMPERATURE_SHIFT`
