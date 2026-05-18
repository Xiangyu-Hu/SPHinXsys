# Matrix-Free A-phi Constriction Heating

This case extends the one-way Joule-heating baseline to a simple necked conductor geometry.

## Purpose
- keep the same matrix-free electromagnetic and one-way thermal chain
- replace the uniform thin band with a symmetric constricted conductor
- check whether Joule heating and temperature rise become stronger in the neck than in the wider shoulders

## Main diagnostics
- `constriction_length`
- `constriction_height`
- `neck_joule_average`
- `shoulder_joule_average`
- `neck_temperature_delta_average`
- `shoulder_temperature_delta_average`
- `joule_neck_to_shoulder_ratio`
- `temperature_delta_neck_to_shoulder_ratio`

## First validated result
At `dp=0.05` with the default neck geometry (`constriction_length=0.20`, `constriction_height=0.16`) and smooth `sin` forcing:
- `converged=1`
- `joule_neck_to_shoulder_ratio=1.209`
- `temperature_delta_neck_to_shoulder_ratio=1.198`

Interpretation:
- the neck already concentrates Joule heating by about twenty percent relative to the shoulders
- the thermal response follows the same trend, which makes this a clean first geometry-induced heating concentration baseline

## First height scan
At `dp=0.05`:
- `constriction_height=0.20`
  - `neck_particle_count=16`
  - `joule_neck_to_shoulder_ratio=1.209033`
  - `temperature_delta_neck_to_shoulder_ratio=1.198151`
- `constriction_height=0.16`
  - `neck_particle_count=16`
  - `joule_neck_to_shoulder_ratio=1.209033`
  - `temperature_delta_neck_to_shoulder_ratio=1.198151`
- `constriction_height=0.10`
  - `neck_particle_count=8`
  - `joule_neck_to_shoulder_ratio=1.316765`
  - `temperature_delta_neck_to_shoulder_ratio=1.301144`

Interpretation:
- the narrowing trend is physically consistent: a smaller neck produces stronger heating concentration
- `0.20` and `0.16` are not meaningfully distinguished at `dp=0.05`, which is why they produce nearly identical ratios
- the `0.10` case is the first clearly resolved stronger constriction at this grid spacing

## Exploratory command
```bash
EM_APHI_STAGGERED_SMOKE_CASE=coulomb_variable_sigma_forced_response \
EM_APHI_FORCED_PROFILE=sin \
EM_APHI_STAGGERED_SMOKE_DP=0.05 \
EM_APHI_STAGGERED_SMOKE_ENABLE_GAUGE=0 \
EM_APHI_ENABLE_GAUGE_PENALTY=1 \
EM_APHI_STAGGERED_SMOKE_OUTER_ITERS=600 \
/home/yongchuan/sphinxsys/build/tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_matrix_free_constriction_heating/bin/test_3d_em_aphi_matrix_free_constriction_heating
```

## Recommended next scan
To separate geometry effects from coarse geometric resolution, prefer one of these next:
- repeat at finer resolution, e.g. `EM_APHI_STAGGERED_SMOKE_DP=0.025`
- or keep `dp=0.05` but choose heights separated by more than one particle spacing

## Optional geometry knobs
- `EM_APHI_CONSTRICTION_LENGTH`
- `EM_APHI_CONSTRICTION_HEIGHT`

## Optional validation knobs
- `EM_APHI_CONSTRICTION_MIN_JOULE_NECK_RATIO`
- `EM_APHI_CONSTRICTION_MIN_TEMPERATURE_NECK_RATIO`
