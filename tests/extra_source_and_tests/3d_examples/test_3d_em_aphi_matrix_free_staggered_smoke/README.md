# Matrix-Free Staggered A-phi Smoke

This case contains relation-based manufactured baselines for the matrix-free staggered `A-step / phi-step` solver and is now the main validation entry point for the current Laplace-structured matrix-free `A-phi` branch.

## Case modes

Current status note:
- `transverse_source_free`: verified, including the full gauge-enabled outer iteration
- `coupled_source_free`: verified as a nonzero-`phi` source-free baseline without gauge in quasi-1D
- `coulomb_variable_sigma_source_free`: current source-free manufactured baseline for the penalty-stabilized route
- `coulomb_variable_sigma_driven`: current source-driven manufactured physical baseline for the penalty-stabilized route
- `coulomb_variable_sigma_forced_response`: transitional physical-response mode; the dedicated standalone target is now preferred

### `transverse_source_free`
- thin quasi-1D body
- manufactured field uses only transverse `A_y(x)`
- `phi = 0`
- `source_phi = 0`
- this is the cleanest first source-free current-continuity baseline
- it is compatible with the gauge-projection diagnostic route

### `coupled_source_free`
- thin quasi-1D body
- manufactured field uses nonzero longitudinal `A_x(x)`
- a discrete nonzero `phi_exact` is constructed from the same scalar Laplace operator so that `source_phi = 0`
- this is the first nontrivial coupled source-free baseline with nonzero `phi`
- in the present quasi-1D setting, it does **not** match a Coulomb-gauge representation, so gauge projection is intentionally rejected for this mode

### `coulomb_variable_sigma_source_free`
- thin 2D-like band
- divergence-compatible transverse reference `A_y(x)`
- spatially varying conductivity `sigma(y)`
- nonzero discrete `phi_exact` with `source_phi = 0`
- this is the current mainline manufactured baseline for the penalty-stabilized matrix-free solver

### `coulomb_variable_sigma_driven`
- thin 2D-like band
- nonzero manufactured source in both the vector-potential and scalar-potential equations
- divergence-compatible transverse reference `A_y(x)` plus directly prescribed nonzero `phi_exact`
- derived physical quantities `E`, `J`, and Joule heat are compared against the same discrete manufactured reference
- this is the current source-driven physical manufactured baseline for the branch

### `coulomb_variable_sigma_forced_response`
- thin 2D-like band
- explicit imposed forcing in the `A_y` equation, without a back-computed manufactured exact solution
- symmetric variable conductivity `sigma(y)` and symmetric forcing are used to produce a simple physically interpretable response
- the main diagnostics are field norms, Joule-heat level, and mirror-symmetry errors in `A_y`, `phi`, and Joule density
- this is the current simple physical-response baseline; it is not a strong manufactured-reference case

## Naming note
- `coulomb_variable_sigma_*` does **not** mean a Coulomb source term is injected directly
- `coulomb` here means the reference setup is intended to stay compatible with the Coulomb-style `div(A)` control route we have been studying
- `variable_sigma` means the conductivity varies across the band as `sigma(y)`, which redistributes `phi`, `E`, `J`, and Joule heating even when the imposed forcing is simple

## What the case does
- builds a thin quasi-1D or thin-band `SolidBody`
- builds an `InnerRelation` and `MatrixFreePairwiseGraph`
- assembles source terms from the same discrete graph operators used by the solver
- runs the staggered matrix-free outer iteration
- reports residuals, discrete-reference field errors, and divergence diagnostics

## Current mainline route
- use `EM_APHI_STAGGERED_SMOKE_CASE=coulomb_variable_sigma_source_free`
- keep `EM_APHI_STAGGERED_SMOKE_ENABLE_GAUGE=0`
- enable `EM_APHI_ENABLE_GAUGE_PENALTY=1`
- rely on the default baseline `EM_APHI_GAUGE_PENALTY_COEFF=1` unless running an explicit penalty experiment
- use `divergence_a_error_l2` rather than absolute `divergence_a_l2` as the primary divergence metric

## Resolution-scan entry
The smoke case now supports direct manufactured convergence-style runs through:
- `EM_APHI_STAGGERED_SMOKE_DP`

Useful current settings are:
- `dp=0.1`
- `dp=0.05`
- `dp=0.025`

For finer resolutions, increase `EM_APHI_STAGGERED_SMOKE_OUTER_ITERS` beyond the coarse-grid default when needed.

## Gauge diagnostics
- `EM_APHI_STAGGERED_SMOKE_ENABLE_GAUGE=1` keeps the explicit gauge route available as a diagnostic path
- `EM_APHI_STAGGERED_SMOKE_GAUGE_ONLY=1` runs isolated gauge-step analysis
- the dense reference gauge solve remains useful for operator diagnostics
- however, hard post-gauge projection is **not** the current production-facing route for this branch

## Useful environment variables
- `EM_APHI_STAGGERED_SMOKE_CASE`
- `EM_APHI_STAGGERED_SMOKE_DP`
- `EM_APHI_STAGGERED_SMOKE_OUTER_ITERS`
- `EM_APHI_STAGGERED_SMOKE_A_ITERS`
- `EM_APHI_STAGGERED_SMOKE_PHI_ITERS`
- `EM_APHI_STAGGERED_SMOKE_A_TOL`
- `EM_APHI_STAGGERED_SMOKE_PHI_TOL`
- `EM_APHI_STAGGERED_SMOKE_RESIDUAL_TOL`
- `EM_APHI_STAGGERED_SMOKE_ENABLE_GAUGE`
- `EM_APHI_STAGGERED_SMOKE_GAUGE_ONLY`
- `EM_APHI_STAGGERED_SMOKE_GAUGE_OPERATOR`
- `EM_APHI_STAGGERED_SMOKE_GAUGE_ITERS`
- `EM_APHI_STAGGERED_SMOKE_GAUGE_TOL`
- `EM_APHI_STAGGERED_SMOKE_GAUGE_SEED`
- `EM_APHI_ENABLE_GAUGE_PENALTY`
- `EM_APHI_GAUGE_PENALTY_COEFF`
- `EM_APHI_GAUGE_PENALTY_RAMP_ITERS`
- `EM_APHI_GAUGE_PENALTY_INITIAL_RATIO`
- `EM_APHI_STAGGERED_SMOKE_REQUIRE_CONVERGED`
- `EM_APHI_STAGGERED_SMOKE_MAX_RESIDUAL_AY`
- `EM_APHI_STAGGERED_SMOKE_MAX_DIVA_ERROR`
- `EM_APHI_STAGGERED_SMOKE_MAX_AY_ERROR`
- `EM_APHI_STAGGERED_SMOKE_MAX_PHI_ERROR`
- `EM_APHI_STAGGERED_SMOKE_MAX_E_ERROR`
- `EM_APHI_STAGGERED_SMOKE_MAX_J_ERROR`
- `EM_APHI_STAGGERED_SMOKE_MAX_JOULE_ERROR`
- `EM_APHI_STAGGERED_SMOKE_MAX_AY_MIRROR`
- `EM_APHI_STAGGERED_SMOKE_MAX_PHI_MIRROR`
- `EM_APHI_STAGGERED_SMOKE_MAX_JOULE_MIRROR`

## Validation gate
The smoke executable can now act as a lightweight regression gate.

If any of these are set, the program still prints the usual summary line, but it will return a nonzero exit code when a requested validation condition is not satisfied:
- `EM_APHI_STAGGERED_SMOKE_REQUIRE_CONVERGED=1`
- `EM_APHI_STAGGERED_SMOKE_MAX_RESIDUAL_AY=...`
- `EM_APHI_STAGGERED_SMOKE_MAX_DIVA_ERROR=...`
- `EM_APHI_STAGGERED_SMOKE_MAX_AY_ERROR=...`
- `EM_APHI_STAGGERED_SMOKE_MAX_PHI_ERROR=...`

This is intended for future regression-style checks of the penalty-stabilized mainline route.

Current recommended medium-resolution regression baseline:
```bash
EM_APHI_STAGGERED_SMOKE_CASE=coulomb_variable_sigma_source_free \
EM_APHI_STAGGERED_SMOKE_DP=0.05 \
EM_APHI_STAGGERED_SMOKE_ENABLE_GAUGE=0 \
EM_APHI_ENABLE_GAUGE_PENALTY=1 \
EM_APHI_STAGGERED_SMOKE_OUTER_ITERS=600 \
EM_APHI_STAGGERED_SMOKE_REQUIRE_CONVERGED=1 \
EM_APHI_STAGGERED_SMOKE_MAX_DIVA_ERROR=1e-6 \
EM_APHI_STAGGERED_SMOKE_MAX_AY_ERROR=2e-6 \
EM_APHI_STAGGERED_SMOKE_MAX_PHI_ERROR=2e-7 \
EM_APHI_STAGGERED_SMOKE_MAX_E_ERROR=5e-6 \
EM_APHI_STAGGERED_SMOKE_MAX_J_ERROR=6e-6 \
EM_APHI_STAGGERED_SMOKE_MAX_JOULE_ERROR=3e-6 \
/home/yongchuan/sphinxsys/build/tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_matrix_free_staggered_smoke/bin/test_3d_em_aphi_matrix_free_staggered_smoke
```

Current recommended medium-resolution source-driven validation baseline:
```bash
EM_APHI_STAGGERED_SMOKE_CASE=coulomb_variable_sigma_driven \
EM_APHI_STAGGERED_SMOKE_DP=0.05 \
EM_APHI_STAGGERED_SMOKE_ENABLE_GAUGE=0 \
EM_APHI_ENABLE_GAUGE_PENALTY=1 \
EM_APHI_STAGGERED_SMOKE_OUTER_ITERS=600 \
EM_APHI_STAGGERED_SMOKE_REQUIRE_CONVERGED=1 \
EM_APHI_STAGGERED_SMOKE_MAX_DIVA_ERROR=1e-6 \
EM_APHI_STAGGERED_SMOKE_MAX_AY_ERROR=1e-6 \
EM_APHI_STAGGERED_SMOKE_MAX_PHI_ERROR=2e-7 \
EM_APHI_STAGGERED_SMOKE_MAX_E_ERROR=4e-6 \
EM_APHI_STAGGERED_SMOKE_MAX_J_ERROR=5e-6 \
EM_APHI_STAGGERED_SMOKE_MAX_JOULE_ERROR=1e-5 \
/home/yongchuan/sphinxsys/build/tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_matrix_free_staggered_smoke/bin/test_3d_em_aphi_matrix_free_staggered_smoke
```

## Typical summary fields
- `residual_ax_l2`, `residual_ay_l2`, `residual_phi_l2`: staggered-solver residual quality
- `ax_l2_error`, `ay_l2_error`, `phi_l2_error`: error to the discrete manufactured reference
- `divergence_a_l2`: absolute discrete divergence of the numerical solution
- `exact_divergence_a_l2`: absolute discrete divergence of the manufactured reference field
- `divergence_a_error_l2`: difference between numerical and reference divergence; this is the preferred mainline divergence metric
- `divergence_j_l2`: current-continuity diagnostic
- `field_update_l2`, `relative_field_update_l2`: outer fixed-point update size diagnostics; useful when judging fine-grid solver maturity
- `reference_phi_build_residual`: relevant when the nonzero `phi_exact` reference is built internally
- `effective_gauge_penalty_coeff`: actual penalty coefficient used during the run
- `validation_require_converged`, `validation_max_*`: the active validation thresholds for the run
- `validation_pass`: whether the current run satisfied the requested validation gate
- `electric_l2_error`, `current_l2_error`, `joule_l2_error`: physical-field error diagnostics against the discrete reference


Current exploratory simple physical-response command:
```bash
EM_APHI_STAGGERED_SMOKE_CASE=coulomb_variable_sigma_forced_response \
EM_APHI_STAGGERED_SMOKE_DP=0.05 \
EM_APHI_STAGGERED_SMOKE_ENABLE_GAUGE=0 \
EM_APHI_ENABLE_GAUGE_PENALTY=1 \
EM_APHI_STAGGERED_SMOKE_OUTER_ITERS=600 \
/home/yongchuan/sphinxsys/build/tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_matrix_free_staggered_smoke/bin/test_3d_em_aphi_matrix_free_staggered_smoke
```

For this mode, treat `ay_mirror_symmetry_l2`, `phi_mirror_symmetry_l2`, and `joule_mirror_symmetry_l2` as first-pass signals only. A dedicated standalone target `test_3d_em_aphi_matrix_free_physical_response` is now preferred for further physical-response work.
