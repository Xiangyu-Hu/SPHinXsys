# Matrix-Free A-phi TEAM7-Like Three-Region Baseline

This target is the first matrix-free step from the simplified strip validations toward a TEAM7-like layout.

## Purpose
- keep the same matrix-free electromagnetic solver and one-way Joule-heating chain
- replace the strip-only geometry with a single box body partitioned into `air`, `coil`, and `conductor` subregions
- drive the system from a localized source inside the coil block
- report region-wise `A / phi / E / J / Joule / temperature-delta` summaries for the conductor, coil, air, and source regions

## Main diagnostics
- `conductor_particles`
- `coil_particles`
- `air_particles`
- `source_particles`
- `source_centroid_x`
- `conductor_joule_centroid_x`
- `conductor_temperature_delta_centroid_x`
- `conductor_joule_relative_to_source_shift`
- `conductor_temperature_delta_relative_to_source_shift`
- `conductor_to_air_joule_ratio`
- `conductor_to_air_temperature_delta_ratio`
- `conductor_mean_abs_a`
- `conductor_mean_abs_phi`
- `conductor_mean_abs_e`
- `conductor_mean_abs_j`
- `conductor_mean_joule`
- `conductor_max_joule`
- `conductor_mean_temperature_delta`
- `conductor_max_temperature_delta`
- `coil_mean_joule`
- `air_mean_joule`
- `source_mean_joule`

## Exploratory command
```bash
EM_APHI_STAGGERED_SMOKE_CASE=coulomb_variable_sigma_forced_response \
EM_APHI_FORCED_PROFILE=sin \
EM_APHI_STAGGERED_SMOKE_DP=0.05 \
EM_APHI_STAGGERED_SMOKE_ENABLE_GAUGE=0 \
EM_APHI_ENABLE_GAUGE_PENALTY=1 \
EM_APHI_STAGGERED_SMOKE_OUTER_ITERS=600 \
/home/yongchuan/sphinxsys/build/tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_matrix_free_team7_like/bin/test_3d_em_aphi_matrix_free_team7_like
```

## Centered-Gaussian follow-up
Use a centered Gaussian source inside the coil to reduce source-shape bias and inspect how the conductor response shifts relative to the source:

```bash
EM_APHI_STAGGERED_SMOKE_CASE=coulomb_variable_sigma_forced_response \
EM_APHI_FORCED_PROFILE=gaussian \
EM_APHI_FORCED_CENTER_X=0.35 \
EM_APHI_FORCED_WIDTH=0.10 \
EM_APHI_STAGGERED_SMOKE_DP=0.05 \
EM_APHI_STAGGERED_SMOKE_ENABLE_GAUGE=0 \
EM_APHI_ENABLE_GAUGE_PENALTY=1 \
EM_APHI_STAGGERED_SMOKE_OUTER_ITERS=600 \
/home/yongchuan/sphinxsys/build/tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_matrix_free_team7_like/bin/test_3d_em_aphi_matrix_free_team7_like
```

## Optional geometry knobs
- `EM_APHI_TEAM7LIKE_LENGTH`
- `EM_APHI_TEAM7LIKE_HEIGHT`
- `EM_APHI_TEAM7LIKE_WIDTH`
- `EM_APHI_TEAM7LIKE_BOUNDARY_WIDTH`
- `EM_APHI_TEAM7LIKE_CONDUCTOR_*_FRACTION`
- `EM_APHI_TEAM7LIKE_COIL_*_FRACTION`

## Optional material and source knobs
- `EM_APHI_TEAM7LIKE_SIGMA_AIR`
- `EM_APHI_TEAM7LIKE_SIGMA_CONDUCTOR`
- `EM_APHI_TEAM7LIKE_SIGMA_COIL`
- `EM_APHI_TEAM7LIKE_NU_AIR`
- `EM_APHI_TEAM7LIKE_NU_CONDUCTOR`
- `EM_APHI_TEAM7LIKE_NU_COIL`
- `EM_APHI_TEAM7LIKE_COIL_J0`
- `EM_APHI_FORCED_PROFILE`
- `EM_APHI_FORCED_CENTER_X`
- `EM_APHI_FORCED_WIDTH`

## ParaView output (matrix-free TEAM7-like)

After a run with state recording enabled (default in many builds), inspect:

- `output/MatrixFreeTeam7LikeBody_ite_0000000000.vtp`

Environment:

- `EM_APHI_TEAM7LIKE_WRITE_VTP=1` (default) — write VTP after EM + thermal post-process
- `EM_APHI_TEAM7LIKE_WRITE_VTP=0` — skip VTP

Suggested ParaView scalars: `RegionId` (0=air, 1=conductor, 2=coil), `AbsPhi`, `AbsE`, `AbsJ`, `JouleDensity`, `TemperatureDelta`.

Implementation checklist: `../../extra_src/shared/aphi_sphinxsys/records/CURSOR_APHI_TEAM7LIKE_IMPLEMENTATION_CHECKLIST.md`.

## Sparse vs matrix-free comparison (planned)

Target binary: `test_3d_em_aphi_team7_like_sparse_vs_matrix_free` (see checklist). Cases: `EM_APHI_COMPARE_CASE=source_free|driven`.

## Optional validation knobs
- `EM_APHI_TEAM7LIKE_MIN_CONDUCTOR_TO_AIR_JOULE_RATIO`
- `EM_APHI_TEAM7LIKE_MIN_CONDUCTOR_TO_AIR_TEMPERATURE_RATIO`
- `EM_APHI_TEAM7LIKE_MIN_CONDUCTOR_JOULE_SHIFT`
- `EM_APHI_TEAM7LIKE_MIN_CONDUCTOR_TEMPERATURE_SHIFT`

## SYCL：图 USM 与标量场 USM（可选）

- 本目录下的 `run_matrix_free_aphi_sycl_accuracy_smoke.sh`：CPU 基线、SYCL-all、policy selftest、图拓扑 USM、图+场 USM 等组合。
- 上级目录的 `matrix_free_aphi_sycl_usm_smoke.sh`：仅对比 SYCL-all 与「双 USM」两档，便于套用到其它 `test_3d_em_aphi_matrix_free_*` 可执行文件；环境变量说明见 `../../extra_src/shared/aphi_sphinxsys/records/MATRIX_FREE_APHI_SYCL_USM.md`。
