# Laplace A-phi Matrix-Free Worklog

## Current Direction

We are extending the existing SPHinXsys electromagnetic work toward a matrix-free,
Laplace-structured frequency-domain `A-phi` solver.

The first implementation target is:

- staggered single-field style iteration
- matrix-free pairwise Laplace / gradient / divergence operators
- explicit gauge projection
- standalone verification cases before large physical heating runs

The first target is not the final locally coupled `4x4` solve and not the final TEAM7 heating workflow.

## Key Decisions

1. Start with staggered matrix-free updates.
2. Keep local coupled correction as a later upgrade path.
3. Keep global Eigen sparse solve only as a reference path.
4. Reuse SPHinXsys diffusion / relaxation patterns whenever possible.
5. Keep the new work in `tests/extra_source_and_tests/extra_src/shared/` first.

## Test Placement Strategy

- First-phase verification should live as standalone test cases under `tests/extra_source_and_tests/3d_examples/`.
- Later, the smallest stable invariants can be moved into tighter regression or unit-like tests.
- Full physical or geometry-heavy validation should stay as standalone cases.

## File Roadmap

Planned shared files:

- `electromagnetic_aphi_matrix_free_types.h`
- `electromagnetic_aphi_matrix_free_fields.h`
- `electromagnetic_aphi_matrix_free_operators.h`
- `electromagnetic_aphi_matrix_free_operators.hpp`
- `electromagnetic_aphi_matrix_free_residuals.h`
- `electromagnetic_aphi_matrix_free_residuals.hpp`
- `electromagnetic_aphi_matrix_free_solver.h`
- `electromagnetic_aphi_matrix_free_solver.hpp`
- `electromagnetic_aphi_matrix_free_gauge_projection.h`
- `electromagnetic_aphi_matrix_free_gauge_projection.hpp`

Planned first verification cases:

- `test_3d_em_aphi_matrix_free_complex_helmholtz`
- `test_3d_em_aphi_matrix_free_gauge_projection`

## Extraction Map

### From `electromagnetic_aphi_laplace_eigen.*`

Extract or adapt:

- diagnostics layout
- harmonic mean helper
- pairwise diffusion weight helper
- pairwise gradient accumulation logic
- pairwise Laplace accumulation logic
- `divA`, `divJ`, and Joule post-processing definitions

### From `electromagnetic_team7_aphi_frequency_dynamics.*`

Reuse the update style from:

- `ElectricPotentialSourceFromVectorFieldInner`
- `ElectricPotentialSourceFromVectorFieldContact`
- `ScalarRelaxationInnerByName`
- `ScalarRelaxationContactByName`
- `ScalarRelaxationComplexByName`
- `UpdateScalarByRelaxationRateByName`
- `VectorPotentialFrequencyCoupledEquationComplex`
- `VectorPotentialFrequencyCoupledBlockEquationComplex`

### From diffusion splitting

Reuse the algorithmic pattern:

- build local residuals from pairwise interactions
- estimate local diagonal-like scaling
- apply local correction
- iterate and report residual history

## Code-Level Task Breakdown

### Task A

Create shared complex types and lightweight containers.

Status:

- completed

Deliverables:

- `electromagnetic_aphi_matrix_free_types.h`
- `electromagnetic_aphi_matrix_free_fields.h`

### Task B

Create reusable matrix-free operator helpers.

Status:

- completed for the first helper layer

Deliverables:

- `electromagnetic_aphi_matrix_free_operators.h`
- `electromagnetic_aphi_matrix_free_operators.hpp`

Current scope:

- harmonic mean
- pairwise diffusion weight
- scalar Laplace contribution
- scalar gradient contribution
- squared-magnitude helpers

### Task C

Create residual containers and evaluators for the scalar complex Helmholtz prototype.

Status:

- next

Deliverables:

- `electromagnetic_aphi_matrix_free_residuals.h`
- `electromagnetic_aphi_matrix_free_residuals.hpp`

### Task D

Create the scalar matrix-free relaxation solver and the first runnable Helmholtz case.

Deliverables:

- `electromagnetic_aphi_matrix_free_solver.h`
- `electromagnetic_aphi_matrix_free_solver.hpp`
- runnable `test_3d_em_aphi_matrix_free_complex_helmholtz`

### Task E

Create the gauge projection layer and the second runnable verification case.

Deliverables:

- `electromagnetic_aphi_matrix_free_gauge_projection.h`
- `electromagnetic_aphi_matrix_free_gauge_projection.hpp`
- runnable `test_3d_em_aphi_matrix_free_gauge_projection`

## Work Entries

### 2026-05-10

- confirmed that the first route should be the staggered single-field matrix-free route
- confirmed that the local coupled route should be postponed until the single-field route is stable
- confirmed that first-phase tests should be standalone verification cases, not only small unit tests
- created the new matrix-free shared type and field headers
- created the first operator helper headers for matrix-free reuse
- created the landing directories and README specs for:
  - `test_3d_em_aphi_matrix_free_complex_helmholtz`
  - `test_3d_em_aphi_matrix_free_gauge_projection`
- completed the first residual helper layer for the scalar complex Helmholtz prototype
- added `electromagnetic_aphi_matrix_free_residuals.h/.hpp`
- extended the residual layer with a simple edge-based pairwise Laplace builder for smoke-level verification
- completed the first scalar matrix-free solver skeleton
- added `electromagnetic_aphi_matrix_free_solver.h/.hpp`
- added a first runnable smoke target in `test_3d_em_aphi_matrix_free_complex_helmholtz/`
- upgraded `test_3d_em_aphi_matrix_free_complex_helmholtz` from a smoke target to a first runnable manufactured complex Helmholtz case
- the current Helmholtz case uses a 1D pairwise chain with switchable discrete or continuum right-hand side
- updated the test CMake file so the target builds the real case source only
- reviewed the first manufactured Helmholtz runs: discrete and continuum modes matched closely, which supports the residual assembly path
- identified that the current issue is solver speed rather than a wrong residual form: the plain Jacobi-like iteration was accurate but too slow under the earlier default parameters
- local parameter checks showed that increasing relaxation and iteration count reduces the error to the 1e-5 to 1e-6 range for the 1D prototype
- updated the Helmholtz test defaults to a more appropriate verification-oriented configuration
- started Task E by adding shared matrix-free gauge-projection helpers in `electromagnetic_aphi_matrix_free_gauge_projection.h/.hpp`
- added the first runnable gauge-projection prototype case in `test_3d_em_aphi_matrix_free_gauge_projection/`
- the current gauge-projection case uses a 1D chain and focuses on chi solve, gauge transform, divergence reduction, and electric-field invariance before the later full relation-based version
- refined the 1D gauge-projection prototype to use an operator-consistent discrete negative-Laplace source and projection-residual diagnostic
- kept the earlier gradient-based divergence metric as a comparison diagnostic instead of the primary pass/fail criterion
- added a new relation-based shared helper layer: `electromagnetic_aphi_matrix_free_pairwise_graph.h/.hpp`
- this helper builds a matrix-free pairwise graph directly from `MatrixFreeAPhiDiscreteView` and exposes reusable scalar-Laplace and gradient application paths
- this is the first shared step from the 1D prototypes toward real relation-based assembly for Helmholtz, gauge projection, and later full A-phi updates
- upgraded `test_3d_em_aphi_matrix_free_complex_helmholtz` from the earlier synthetic 1D edge-list version to a real body-and-inner-relation version using `MatrixFreePairwiseGraph`
- the current Helmholtz prototype now builds its matrix-free scalar Laplace residual from an SPH body, lattice particles, and an actual inner relation
- relation-based Helmholtz testing showed that the earlier default box thickness made the baseline unnecessarily 3D and slowed plain Jacobi-like convergence
- switched the default relation-based Helmholtz baseline to a thin quasi-1D body so the first real inner-relation validation remains clean and interpretable
- upgraded `test_3d_em_aphi_matrix_free_gauge_projection` from the earlier synthetic prototype to a real body-and-inner-relation version using `MatrixFreePairwiseGraph`
- the current gauge-projection prototype now builds its chi solve, projection residual, divergence-like diagnostic, and electric-field comparison from an SPH body and actual inner relation
- added a new shared relation-based residual layer for the full matrix-free Laplace A-phi formulation: `electromagnetic_aphi_matrix_free_aphi_residuals.h/.hpp`
- this new layer evaluates residuals for `A_x`, `A_y`, `A_z`, `phi`, and also reports `div(A)` and `div(J)` from the shared pairwise graph
- added the first shared relation-based staggered A-phi solver skeleton in `electromagnetic_aphi_matrix_free_aphi_solver.h/.hpp`
- added a new smoke target `test_3d_em_aphi_matrix_free_staggered_smoke` so the new outer-iteration skeleton is compiled and exercised instead of remaining header-only
- upgraded `test_3d_em_aphi_matrix_free_staggered_smoke` from a bare skeleton run to a first relation-based discrete manufactured baseline for the staggered `A-step / phi-step` path
- the current staggered smoke case now assembles `source_ax` and `source_phi` from the same discrete graph operators used by the solver, so it can report solution errors against a discrete manufactured reference instead of only residuals
- made the staggered smoke case tunable through environment variables for outer-iteration count, scalar subsolver iterations, tolerances, and optional gauge enablement
- identified and fixed a convergence-logic issue in `electromagnetic_aphi_matrix_free_aphi_solver.hpp`: when gauge projection is disabled, `div(A)` should remain a diagnostic rather than a hard convergence blocker
- local staggered manufactured runs show that the current outer fixed-point iteration does converge on the relation-based quasi-1D baseline, but it needs substantially more outer iterations than the scalar Helmholtz and gauge prototypes
- updated the staggered smoke defaults to a verification-oriented configuration so the baseline converges without manual parameter overrides
- added `test_3d_em_aphi_matrix_free_staggered_smoke/README.md` to document the current manufactured baseline, its intended scope, and the environment-variable knobs for future scans
- replaced the earlier mixed `A_x + phi` staggered manufactured smoke with a cleaner source-free current-continuity baseline using only transverse `A_y(x)` and `phi = 0`
- in this new quasi-1D manufactured baseline, `source_phi` is kept exactly zero and the discrete current-continuity diagnostic `div(J)` collapses to zero on the default run
- local default run of the updated staggered source-free baseline now gives: `outer_iterations=79`, `converged=1`, `residual_ay_l2=7.06e-06`, `residual_phi_l2=0`, `divergence_a_l2=0`, `divergence_j_l2=0`, `ay_l2_error=1.26e-06`, `phi_l2_error=0`
- this is the first relation-based coupled matrix-free `A-step / phi-step` verification case in the current branch that simultaneously shows correct field recovery and source-free current continuity on the same discrete baseline
- extended the source-free staggered smoke baseline so that the gauge-enabled mode starts from a deliberately contaminated longitudinal `A_x` and `phi` seed while keeping the manufactured target equal to the divergence-free transverse `A_y(x)` solution
- local rebuilt comparison now shows:
  - non-gauge baseline: `outer_iterations=79`, `converged=1`, `residual_ay_l2=6.00e-06`, `divergence_a_l2=0`, `divergence_j_l2=0`, `ay_l2_error=1.17e-06`
  - full `A-step / phi-step / gauge` baseline: `outer_iterations=60`, `converged=1`, `residual_ax_l2=9.23e-06`, `residual_ay_l2=7.54e-08`, `residual_phi_l2=5.28e-07`, `divergence_a_l2=1.27e-07`, `divergence_j_l2=5.90e-07`, `ax_l2_error=1.83e-06`, `ay_l2_error=0`, `phi_l2_error=5.54e-08`
- this is the first relation-based source-free manufactured verification in the branch that exercises the full matrix-free staggered `A-step / phi-step / gauge` outer iteration and shows that the gauge stage removes a seeded longitudinal mode without damaging the transverse physical solution
- added a second manufactured mode `coupled_source_free` to `test_3d_em_aphi_matrix_free_staggered_smoke` so the branch now has a source-free baseline with nonzero `phi`
- in `coupled_source_free`, a discrete `phi_exact` is built from the same scalar Laplace operator using a nonzero longitudinal `A_x(x)` while keeping `source_phi = 0`
- local default run of `coupled_source_free` now gives: `outer_iterations=63`, `converged=1`, `residual_ax_l2=8.56e-06`, `residual_phi_l2=9.55e-06`, `ax_l2_error=1.54e-06`, `phi_l2_error=7.44e-07`, `reference_phi_build_residual=2.48e-06`
- this confirms that the current matrix-free staggered solver can recover a nontrivial source-free coupled `A_x + phi` discrete baseline in the quasi-1D relation-based setting
- also identified an important formulation constraint: the quasi-1D `coupled_source_free` baseline does not match a Coulomb-gauge representation because nonzero `phi` in this setting requires a longitudinal `A_x` contribution
- therefore `coupled_source_free + gauge` is now rejected intentionally with the message `coupled_source_free_does_not_match_coulomb_gauge_in_quasi_1d` instead of producing misleading diagnostics
- added a third manufactured mode `coulomb_variable_sigma_source_free` to `test_3d_em_aphi_matrix_free_staggered_smoke`
- this new mode uses a thin 2D-like band, a divergence-free transverse `A_y(x)` field, and spatially varying `sigma(y)` so that `source_phi = 0` still produces a nonzero discrete `phi_exact`
- local non-gauge run of `coulomb_variable_sigma_source_free` now gives: `points=120`, `outer_iterations=78`, `converged=1`, `residual_ax_l2=4.91e-07`, `residual_ay_l2=9.88e-06`, `residual_phi_l2=1.13e-06`, `ax_l2_error=9.53e-09`, `ay_l2_error=8.30e-07`, `phi_l2_error=7.41e-08`, `reference_phi_build_residual=9.91e-07`
- this is the first relation-based baseline in the branch that combines all of the following at once: nonzero `phi`, `source_phi = 0`, spatially varying conductivity, and a Coulomb-compatible divergence-free reference `A`
- however, the gauge-enabled run of `coulomb_variable_sigma_source_free` does not converge under the current staggered outer iteration, even when the outer and scalar subsolver iteration budgets are increased substantially
- representative gauge-enabled run: `outer_iterations=320`, `converged=0`, `residual_ax_l2=1.74e-01`, `residual_ay_l2=1.92e-01`, `residual_phi_l2=1.95e-01`, `phi_l2_error=4.64e-03`
- reducing the seeded gauge contamination amplitude did not materially change this stalled state, so the present issue looks structural rather than a simple lack of iteration budget
- current interpretation: the 2D-like variable-conductivity source-free baseline itself is valid and useful, but the full `A-step / phi-step / gauge` fixed-point iteration still needs a more careful treatment in this regime
- added shared gauge-step diagnostics to `electromagnetic_aphi_matrix_free_aphi_solver.h/.hpp` without changing the outer-iteration structure
- the staggered solver state now records, for the most recent gauge application: phi-reference offset after the phi solve, after the raw gauge update, and after the final reference enforcement
- it also records `div(A)` and `div(J)` before and after the gauge step, plus `E` and `J` norms before and after the gauge step together with their L2 changes
- exposed these new diagnostics in the summary line printed by `test_3d_em_aphi_matrix_free_staggered_smoke`
- purpose of this instrumentation: isolate whether the current stalled `coulomb_variable_sigma_source_free + gauge` run is caused primarily by gauge-invariance breakage in `E`/`J`, by phi-reference enforcement side effects, or by the outer fixed-point coupling itself
- extended the shared staggered A-phi solver diagnostics to separate raw gauge-update effects from the final reference-enforced state
- new gauge-step diagnostics now include `chi` and `grad(chi)` norms, `div(A)` and `div(J)` before / after raw update / after final reference enforcement, and `E` / `J` norm changes both before and after the final reference projection
- added a new `EM_APHI_STAGGERED_SMOKE_GAUGE_ONLY=1` analysis mode to `test_3d_em_aphi_matrix_free_staggered_smoke`
- this mode first runs a non-gauge staggered solve for a configurable number of outer iterations and then applies exactly one isolated shared gauge step, so the variable-sigma gauge issue can be diagnosed outside the full outer fixed-point loop
- purpose of the new mode: distinguish three possibilities cleanly on the next run
  - the gauge update itself breaks `E` / `J`
  - the final phi-reference enforcement perturbs the projected state
  - the isolated gauge step is acceptable, but the full `A-step / phi-step / gauge` outer coupling is what stalls

## 2026-05-11: Gauge operator-consistency wiring
- Added a switchable gauge-projection backend in the staggered matrix-free A-phi solver.
- Set `use_operator_consistent_gauge_projection = true` as the default solver behavior.
- Wired `applyMatrixFreeAPhiGaugeProjectionStep(...)` and the full staggered outer iteration to pass the selected gauge operator through to the scalar `chi` solve.
- Extended the staggered smoke test with `EM_APHI_STAGGERED_SMOKE_GAUGE_OPERATOR=consistent|laplace` so the old Laplace-based gauge solve can be compared directly against the new operator-consistent variant.
- Extended both the full and isolated gauge-step smoke summaries to print `gauge_operator_mode`, making result comparison easier when diagnosing `div(A)` behavior.
- No cases were run in this step; the next comparison should be done by building `test_3d_em_aphi_matrix_free_staggered_smoke` and running the `coulomb_variable_sigma_source_free` baseline in both `consistent` and `laplace` modes.

## 2026-05-11: Consistent gauge stabilization pass
- Adjusted the operator-consistent `D(G chi)` gauge assembly so its diagonal proxy uses the sign convention expected by the scalar Jacobi update.
- Added scalar solver diagonal diagnostics: `current_min_diagonal_abs_`, `current_max_diagonal_abs_`, and `current_nonfinite_diagonal_count_`.
- Exposed these gauge-solver diagonal diagnostics in both the full staggered smoke summary and the isolated gauge-step summary.
- Tightened the default consistent-gauge smoke parameters to safer values:
  - `EM_APHI_STAGGERED_SMOKE_GAUGE_ITERS` default: `2000`
  - `EM_APHI_STAGGERED_SMOKE_GAUGE_RELAX` default: `0.05`
  - `EM_APHI_STAGGERED_SMOKE_GAUGE_DIAG_REG` default: `1e-6`
- No cases were run in this step; the next check should rerun the `coulomb_variable_sigma_source_free` case in `consistent` mode and compare the new diagonal diagnostics against the previous `inf/nan` behavior.

## 2026-05-11: Consistent gauge RHS sign adjustment
- After stabilizing the consistent `D(G chi)` iteration, the gauge step remained finite but moved `div(A)` in the wrong direction.
- Updated the consistent-gauge branch to solve against `-div(A)` instead of `div(A)` so the RHS sign matches the current `D(G)` implementation and the scalar negative-Laplace residual convention used by the solver.
- No cases were run in this step; the next check should rerun the `coulomb_variable_sigma_source_free` baseline in `consistent` full-gauge and isolated-gauge modes and compare `gauge_div_a_before_l2` vs `gauge_div_a_after_raw_l2`.

## 2026-05-11: Static seeded-gauge operator diagnostic
- Added a static consistency diagnostic to the staggered smoke case for gauge-enabled runs.
- The smoke summary now reports four seeded-gauge metrics:
  - `gauge_seed_div_delta_l2`: L2 norm of `div(A_seeded) - div(A_exact)`
  - `gauge_seed_operator_response_l2`: L2 norm of `applyScalarDivergenceOfGradientFromGraph(chi_seed)`
  - `gauge_seed_match_l2`: L2 mismatch between the seeded divergence delta and the operator response
  - `gauge_seed_match_neg_l2`: L2 mismatch between the seeded divergence delta and the negated operator response
- These metrics are intended to distinguish a pure sign error from a deeper operator inconsistency in the current repeated-gradient `D(G)` construction.
- No cases were run in this step; the next check should rerun the isolated `consistent` gauge case and inspect `gauge_seed_match_l2` vs `gauge_seed_match_neg_l2`.

## 2026-05-11: Gauge zero-space compatibility pass
- Reverted the temporary consistent-gauge sign hacks after the static seeded-gauge diagnostic showed that `div(A_seeded)-div(A_exact)` and `D(G chi_seed)` already match with the same sign.
- Reverted the temporary negative diagonal proxy in the consistent `D(G)` assembly.
- Added gauge RHS zero-mean projection inside `applyMatrixFreeAPhiGaugeProjectionStep(...)` and recorded:
  - `gauge_rhs_mean_abs_before_projection`
  - `gauge_rhs_mean_abs_after_projection`
- Replaced the gauge scalar solve with a dedicated zero-mean-aware iteration loop that removes the mean of `chi` after every Jacobi update.
- Added `gauge_chi_mean_abs_after_solve` to confirm whether the nullspace handling is effective at the end of the gauge solve.
- No cases were run in this step; the next check should rerun the pure seeded-gauge case with `PREGAUGE_OUTER_ITERS=0` and inspect the new RHS/chi zero-space diagnostics together with `div(A)` before/after projection.

## 2026-05-11: Guarded Consistent-Gauge Iteration

- Rechecked the relationship between SPHinXsys diffusion-style local projection and our consistent gauge operator.
- Conclusion stayed the same: existing SPHinXsys splitting/projection patterns are useful inspiration, but they do not directly provide a stable solver for the repeated-gradient `D(G chi)` projection operator.
- Updated `solveGaugeComponent(...)` so the consistent-gauge branch now uses a guarded normalized residual update instead of the previous raw explicit step.
- New protection logic:
  - scale the update by the current maximum residual magnitude;
  - keep a `last_stable_chi` snapshot;
  - stop and roll back if the residual becomes non-finite;
  - stop and roll back if the residual grows too aggressively relative to the previous accepted state;
  - reject non-finite candidate `chi` values before accepting the update.
- Immediate validation target remains the one-shot pure-gauge isolation case with `PREGAUGE_OUTER_ITERS=0`.

## 2026-05-11: Gauge-Only Reference Projection Path

- Added a small Eigen-based reference projection path inside `test_3d_em_aphi_matrix_free_staggered_smoke.cpp` for `gauge_only` diagnostics.
- This path explicitly assembles the current graph-based consistent `D(G)` operator by probing basis vectors with `applyScalarDivergenceOfGradientFromGraph(...)`.
- The reference solve is meant only for diagnosis, not as the final production route.
- New comparison fields are intended to tell us whether the projection equation itself can reduce `div(A)` when the operator is solved directly, or whether the difficulty is deeper than the current iterative update.

## 2026-05-11: First Gauge-Penalty Path

- Added a first explicit-lagged gauge-penalty path to the staggered A-step.
- New solver controls:
  - `enable_gauge_penalty`
  - `gauge_penalty_coefficient`
- First formulation is conservative and easy to diagnose:
  - compute `div(A)` from the current outer-iteration state;
  - compute `grad(divA)` from the same graph operator;
  - subtract `gamma * grad(divA)` from the A-step RHS;
  - add `gamma * grad(divA)` into residual evaluation so diagnostics match the actual stabilized equation.
- This is intentionally an explicit/lagged stabilization, not a fully coupled reformulation yet.

## 2026-05-11: Manufactured Source Aligned With Gauge Penalty

- The first gauge-penalty scans showed a meaningful tradeoff: `div(A)` and `Ay` residuals improved, but the reported field errors also grew quickly.
- That comparison was not fully fair yet, because the manufactured source still matched the unpenalized equation.
- Updated the smoke test so that when `EM_APHI_ENABLE_GAUGE_PENALTY=1`, the manufactured source is also augmented by the exact-field penalty contribution `gamma * grad(divA_exact)`.
- This keeps the discrete reference field aligned with the stabilized equation and makes the next penalty scans much more interpretable.

## 2026-05-11: Exact Divergence Baseline Diagnostics

- The updated manufactured-source-aligned penalty scans suggest that absolute `divergence_a_l2` is not a reliable standalone metric in the current raw-operator setting.
- Reason: the discrete reference field itself may already have nonzero `D(A_exact)` under the current graph gradient/divergence pair.
- Added diagnostics to report:
  - `exact_divergence_a_l2`
  - `divergence_a_error_l2`
  - `pregauge_divergence_a_error_l2`
  - `post_gauge_divergence_a_error_l2`
  - `reference_gauge_div_a_error_l2`
- The next comparisons should focus on whether the numerical solution approaches the discrete reference divergence, not only whether absolute `divA` goes to zero.

## 2026-05-11: Penalty Ramp And Outer-Iteration Guard

- Added a first continuation-style ramp for the gauge-penalty coefficient.
- New controls:
  - `gauge_penalty_ramp_iterations`
  - `gauge_penalty_initial_ratio`
- During the outer iteration, the effective penalty now ramps from `initial_ratio * gamma_target` to `gamma_target`.
- Added a last-stable-state rollback guard in the staggered solver: if field residuals become non-finite or grow too aggressively, the solver restores the previous accepted state instead of leaving the fields contaminated.
- Smoke summaries now print `effective_gauge_penalty_coeff` so ramped runs can be interpreted correctly.
## 2026-05-11: Penalty baseline freeze and convergence-entry cleanup
- Stopped treating gauge-penalty coefficient scanning as the current bottleneck for the main route.
- Froze the smoke-test penalty baseline around the currently validated stable setting: when `EM_APHI_ENABLE_GAUGE_PENALTY=1`, the default `EM_APHI_GAUGE_PENALTY_COEFF` is now `1.0` unless overridden explicitly.
- Simplified the default penalty continuation knobs for the mainline validation path: `EM_APHI_GAUGE_PENALTY_RAMP_ITERS` now defaults to `0` and `EM_APHI_GAUGE_PENALTY_INITIAL_RATIO` defaults to `1.0` unless an experiment requests otherwise.
- Added `EM_APHI_STAGGERED_SMOKE_DP` so the manufactured staggered smoke case can be rerun directly at multiple resolutions without editing the source.
- Extended both smoke summaries to print `dp=...`, making it easier to compare convergence-style runs across resolutions.
- Current interpretation: the next meaningful validation step is resolution-based manufactured convergence under the penalty-stabilized main route rather than more penalty-coefficient tuning.
## 2026-05-11: Penalty-Stabilized Mainline Baseline Validated

- Ran the penalty-stabilized `coulomb_variable_sigma_source_free` manufactured baseline at three resolutions using the current mainline settings:
  - `EM_APHI_STAGGERED_SMOKE_ENABLE_GAUGE=0`
  - `EM_APHI_ENABLE_GAUGE_PENALTY=1`
  - default `EM_APHI_GAUGE_PENALTY_COEFF=1`
  - default no-ramp penalty application
- Resolution scan results:
  - `dp=0.1`: `points=30`, `outer_iterations=14`, `converged=1`, `residual_ay_l2=5.05e-06`, `divergence_a_error_l2=5.41e-08`, `ay_l2_error=6.83e-07`, `phi_l2_error=2.64e-08`
  - `dp=0.05`: default `outer_iterations=160` was not enough to satisfy the convergence criterion, but increasing to `EM_APHI_STAGGERED_SMOKE_OUTER_ITERS=600` gave `outer_iterations=175`, `converged=1`, `residual_ay_l2=8.95e-06`, `divergence_a_error_l2=5.07e-08`, `ay_l2_error=1.30e-06`, `phi_l2_error=1.40e-07`
  - `dp=0.025`: default `outer_iterations=160` was also insufficient, but increasing to `EM_APHI_STAGGERED_SMOKE_OUTER_ITERS=800` gave `outer_iterations=135`, `converged=1`, `residual_ay_l2=9.89e-06`, `divergence_a_error_l2=4.11e-08`, `ay_l2_error=6.92e-07`, `phi_l2_error=7.72e-08`
- Interpretation:
  - the penalty-stabilized route is now the validated mainline baseline for the current Laplace-structured matrix-free `A-phi` formulation;
  - the previous non-converged `dp=0.05` and `dp=0.025` runs were caused by outer-iteration budget, not by a formulation failure;
  - `divergence_a_error_l2` is the correct main divergence metric for this route because the discrete reference field itself carries a nonzero absolute `div(A)` baseline under the current raw graph operator pair.
- Current practical baseline for further development:
  - use `coulomb_variable_sigma_source_free`
  - keep post-gauge projection disabled in the main solve path
  - enable penalty stabilization with `gamma = 1`
  - increase outer-iteration budget for finer resolutions when doing convergence-style runs.
- Current development conclusion:
  - hard post-projection should remain a diagnostic/reference tool;
  - penalty stabilization should remain the production-facing route for the present matrix-free staggered Laplace `A-phi` branch.
## 2026-05-11: Smoke Validation Gate

- Extended `test_3d_em_aphi_matrix_free_staggered_smoke` so it can act as a lightweight regression gate for the current penalty-stabilized mainline baseline.
- Added optional validation controls:
  - `EM_APHI_STAGGERED_SMOKE_REQUIRE_CONVERGED`
  - `EM_APHI_STAGGERED_SMOKE_MAX_RESIDUAL_AY`
  - `EM_APHI_STAGGERED_SMOKE_MAX_DIVA_ERROR`
  - `EM_APHI_STAGGERED_SMOKE_MAX_AY_ERROR`
  - `EM_APHI_STAGGERED_SMOKE_MAX_PHI_ERROR`
- The smoke summary now reports the active validation thresholds plus `validation_pass=0|1`.
- If any requested validation condition is violated, the executable now returns a nonzero exit code instead of only printing diagnostics.
- Purpose of this step: turn the manually validated penalty-stabilized baseline into a reusable regression-style entry point before moving on to larger physical baselines.
## 2026-05-11: First Passing Regression-Style Gate

- Confirmed that the penalty-stabilized mainline baseline now supports an explicit pass/fail regression-style gate.
- Passing validated command (current recommended medium-resolution regression baseline):
  - `EM_APHI_STAGGERED_SMOKE_CASE=coulomb_variable_sigma_source_free`
  - `EM_APHI_STAGGERED_SMOKE_DP=0.05`
  - `EM_APHI_STAGGERED_SMOKE_ENABLE_GAUGE=0`
  - `EM_APHI_ENABLE_GAUGE_PENALTY=1`
  - `EM_APHI_STAGGERED_SMOKE_OUTER_ITERS=600`
  - `EM_APHI_STAGGERED_SMOKE_REQUIRE_CONVERGED=1`
  - `EM_APHI_STAGGERED_SMOKE_MAX_RESIDUAL_AY=1e-5`
  - `EM_APHI_STAGGERED_SMOKE_MAX_DIVA_ERROR=1e-6`
  - `EM_APHI_STAGGERED_SMOKE_MAX_AY_ERROR=2e-6`
  - `EM_APHI_STAGGERED_SMOKE_MAX_PHI_ERROR=2e-7`
- Observed passing result:
  - `outer_iterations=176`
  - `converged=1`
  - `residual_ay_l2=8.70e-06`
  - `divergence_a_error_l2=5.13e-08`
  - `ay_l2_error=1.23e-06`
  - `phi_l2_error=1.29e-07`
  - `validation_pass=1`
- This is the first run in the current branch that upgrades the penalty-stabilized staggered manufactured baseline from a manually interpreted diagnostic to an explicit regression-style validation point.
## 2026-05-11: Outer Relaxation Improves Fine-Grid Solver Behavior

- Added outer-field under-relaxation and update-based convergence diagnostics to the staggered matrix-free `A-phi` solver.
- New solver controls:
  - `outer_relaxation_factor`
  - `update_tolerance`
- New reported diagnostics:
  - `field_update_l2`
  - `relative_field_update_l2`
- Fine-grid recheck on the current penalty-stabilized mainline baseline (`coulomb_variable_sigma_source_free`, `dp=0.025`, penalty on, gauge off) showed a clear solver-behavior improvement:
  - previous run without the new outer stabilization reached the outer-iteration cap and did not satisfy the convergence gate;
  - with `outer_relaxation_factor=0.8`, the same case reached `converged=1` in `outer_iterations=61`.
- Observed fine-grid result with the new outer stabilization:
  - `field_update_l2=1.17e-08`
  - `relative_field_update_l2=1.64e-08`
  - `divergence_a_error_l2=4.83e-08`
  - `ay_l2_error=8.96e-07`
  - `phi_l2_error=1.16e-07`
  - but `residual_ay_l2=2.22e-05`, which remained above the currently chosen fine-grid regression threshold.
- Current interpretation:
  - the outer fixed-point solver itself became materially more stable and mature under the new under-relaxed update strategy;
  - the remaining fine-grid regression miss is now dominated by the chosen absolute `residual_ay_l2` gate rather than by obviously poor field error or divergence behavior;
  - for the current mainline route, absolute residual should be treated as a secondary diagnostic, while discrete-reference field error and `divergence_a_error_l2` remain the more trustworthy primary validation targets.
- Development implication:
  - do not treat the remaining fine-grid residual gate miss as a formulation failure;
  - keep the current outer stabilization;
  - continue the mainline validation path toward more physical/source-driven test content rather than pausing again for penalty-coefficient tuning.
## 2026-05-11: Physical-Field Validation Added To Mainline Smoke

- Extended the penalty-stabilized mainline smoke baseline to compare not only `A` and `phi`, but also the derived physical quantities:
  - electric field `E`
  - current density `J`
  - Joule-heating density
- New reported diagnostics:
  - `electric_l2_error`
  - `current_l2_error`
  - `joule_l2_error`
- New optional validation-gate controls:
  - `EM_APHI_STAGGERED_SMOKE_MAX_E_ERROR`
  - `EM_APHI_STAGGERED_SMOKE_MAX_J_ERROR`
  - `EM_APHI_STAGGERED_SMOKE_MAX_JOULE_ERROR`
- First measured medium-resolution mainline result (`coulomb_variable_sigma_source_free`, `dp=0.05`, penalty on, gauge off, converged baseline):
  - `electric_l2_error=2.99e-06`
  - `current_l2_error=3.65e-06`
  - `joule_l2_error=1.62e-06`
- Recommended first conservative physical-field gates for this baseline:
  - `EM_APHI_STAGGERED_SMOKE_MAX_E_ERROR=5e-6`
  - `EM_APHI_STAGGERED_SMOKE_MAX_J_ERROR=6e-6`
  - `EM_APHI_STAGGERED_SMOKE_MAX_JOULE_ERROR=3e-6`
- Development meaning: the current manufactured baseline now exercises the chain from `A/phi` to `E/J/Joule`, which is a more physical validation layer than field-only checks while still retaining a strong discrete reference.
## 2026-05-11: Source-Driven Manufactured Physical Baseline Validated

- Added and exercised a new source-driven manufactured baseline in `test_3d_em_aphi_matrix_free_staggered_smoke`:
  - `EM_APHI_STAGGERED_SMOKE_CASE=coulomb_variable_sigma_driven`
- Purpose of this new mode:
  - keep the geometry simple and the reference strong;
  - move one step beyond the earlier source-free consistency checks;
  - validate the full chain `source -> A/phi -> E/J -> Joule heat` while still using a discrete manufactured reference built from the same graph operators.
- Baseline run (`dp=0.05`, gauge off, penalty on, outer stabilization active) gave:
  - `outer_iterations=117`, `converged=1`
  - `residual_ay_l2=9.91e-06`
  - `divergence_a_error_l2=6.84e-08`
  - `ax_l2_error=1.43e-07`
  - `ay_l2_error=3.80e-07`
  - `phi_l2_error=5.12e-08`
  - `electric_l2_error=2.00e-06`
  - `current_l2_error=2.45e-06`
  - `joule_l2_error=7.11e-06`
  - `source_ax_l2=2.72e-01`
  - `source_ay_l2=1.31e+01`
  - `source_phi_l2=8.85e+00`
- An earlier large `phi_l2_error` in this mode was traced to a reference-offset mismatch rather than a solver failure; the reference `phi` field is now aligned with the solver's enforced phi reference before comparison.
- Current interpretation:
  - this new mode is the first validated source-driven physical manufactured baseline in the current branch;
  - it confirms that the stabilized solver can reproduce not only `A` and `phi`, but also derived `E`, `J`, and Joule-heating response under a nontrivial manufactured source.
- Recommended first conservative validation thresholds for this mode at `dp=0.05`:
  - `EM_APHI_STAGGERED_SMOKE_MAX_DIVA_ERROR=1e-6`
  - `EM_APHI_STAGGERED_SMOKE_MAX_AY_ERROR=1e-6`
  - `EM_APHI_STAGGERED_SMOKE_MAX_PHI_ERROR=2e-7`
  - `EM_APHI_STAGGERED_SMOKE_MAX_E_ERROR=4e-6`
  - `EM_APHI_STAGGERED_SMOKE_MAX_J_ERROR=5e-6`
  - `EM_APHI_STAGGERED_SMOKE_MAX_JOULE_ERROR=1e-5`


## 2026-05-12: Simple physical-response baseline
- clarified the naming of the `coulomb_variable_sigma_*` family: `variable_sigma` refers to the imposed conductivity variation `sigma(y)`, not to a direct Coulomb source term
- the existing `source_free` and `driven` modes remain strong manufactured-reference checks; they use a Coulomb-compatible transverse reference field together with varying conductivity
- added a new exploratory mode `coulomb_variable_sigma_forced_response`
- this new mode applies a simple explicit forcing in the `A_y` equation without back-computing a manufactured exact solution
- the purpose of this mode is to start validating physically interpretable response patterns before moving to heavier geometry or heating cases
- added physical-response diagnostics to the smoke summary:
  - `electric_field_l2`
  - `current_density_l2`
  - `joule_density_l2`
  - `joule_density_max`
  - `ay_mirror_symmetry_l2`
  - `phi_mirror_symmetry_l2`
  - `joule_mirror_symmetry_l2`
- updated the smoke summary with `has_discrete_reference=0/1` so reference-based L2 errors are only treated as manufactured-solution checks when a true discrete reference exists
- updated the validation logic so `forced_response` can be gated by mirror-symmetry thresholds instead of manufactured-reference error thresholds
- new symmetry-oriented gate variables:
  - `EM_APHI_STAGGERED_SMOKE_MAX_AY_MIRROR`
  - `EM_APHI_STAGGERED_SMOKE_MAX_PHI_MIRROR`
  - `EM_APHI_STAGGERED_SMOKE_MAX_JOULE_MIRROR`
- first `dp=0.05` exploratory run of `coulomb_variable_sigma_forced_response` converged and produced stable field / Joule norms, which is enough to keep this mode as the next-step physical plausibility baseline

## 2026-05-12: Physical-response case split
- split the exploratory `coulomb_variable_sigma_forced_response` path out of the overloaded `staggered_smoke` workflow into a dedicated standalone target
- new target: `test_3d_em_aphi_matrix_free_physical_response`
- purpose of the split: keep manufactured/regression baselines separate from simpler physical-plausibility checks
- the new physical-response target keeps the same thin-band geometry and penalty-stabilized matrix-free solver route, but it is documented and reported as a non-manufactured validation mode
- added a `phi_mirror_antisymmetry_l2` diagnostic in the standalone physical-response target so `phi` parity can be checked with a physically better-matched quantity
- added `EM_APHI_STAGGERED_SMOKE_MAX_PHI_ANTIMIRROR` to the standalone physical-response target for symmetry-oriented gating

- refined the standalone physical-response diagnostics for `phi`: the target now also reports mean-centered mirror / anti-mirror parity errors so scalar-potential symmetry is less sensitive to reference-offset bias
- added two more optional physical-response gate variables:
  - `EM_APHI_STAGGERED_SMOKE_MAX_PHI_CENTERED_MIRROR`
  - `EM_APHI_STAGGERED_SMOKE_MAX_PHI_CENTERED_ANTIMIRROR`
- current role split:
  - `test_3d_em_aphi_matrix_free_staggered_smoke`: manufactured/reference-driven validation and regression gates
  - `test_3d_em_aphi_matrix_free_physical_response`: simple physical-response, symmetry, and Joule-pattern validation

- first standalone `physical_response` run at `dp=0.05` showed:
  - `converged=1`, `outer_iterations=97`, `residual_ay_l2=9.73e-06`
  - `ay_mirror_symmetry_l2=1.65e-09`
  - `joule_mirror_symmetry_l2=9.67e-09`
  - raw `phi` parity looked poor, but the refined centered metric gave `phi_centered_mirror_antisymmetry_l2=3.61e-09`
- interpretation: the scalar potential in this simple forced-response setup is anti-symmetric up to the expected removable constant offset introduced by the reference constraint
- adopted first conservative physical-response validation gate at `dp=0.05`:
  - `EM_APHI_STAGGERED_SMOKE_REQUIRE_CONVERGED=1`
  - `EM_APHI_STAGGERED_SMOKE_MAX_RESIDUAL_AY=1.2e-5`
  - `EM_APHI_STAGGERED_SMOKE_MAX_AY_MIRROR=1e-6`
  - `EM_APHI_STAGGERED_SMOKE_MAX_JOULE_MIRROR=1e-6`
  - `EM_APHI_STAGGERED_SMOKE_MAX_PHI_CENTERED_ANTIMIRROR=1e-6`
- current validation split is now:
  - manufactured/reference-driven: `test_3d_em_aphi_matrix_free_staggered_smoke`
  - physical-response/symmetry-driven: `test_3d_em_aphi_matrix_free_physical_response`

## 2026-05-12: One-way Joule-heating baseline
- added a new standalone target `test_3d_em_aphi_matrix_free_joule_heating_baseline`
- this target reuses the converged penalty-stabilized forced-response electromagnetic solve, then feeds the resulting `JouleHeatSource` into a simple one-way thermal diffusion update
- reused the existing SPHinXsys diffusion closure route for `Temperature` together with the already established `JouleHeatSource` / `TemperatureChangeRateByJoule` variable names
- first `dp=0.05` exploratory run showed:
  - `converged=1`, `outer_iterations=97`, `residual_ay_l2=9.76e-06`
  - `joule_density_l2=1.09e-01`, `joule_density_max=1.20e-01`
  - `thermal_time=1.0e-01`, `thermal_steps=4`
  - `temperature_average=300.010864`, `temperature_delta_average=1.08968e-02`
  - `expected_temperature_delta_average=1.08995e-02`, `temperature_delta_average_error=2.65e-06`
  - `temperature_mirror_symmetry_l2=0`, `temperature_delta_mirror_symmetry_l2=0`
- interpretation: the one-way thermal response is symmetric and its mean temperature rise matches the injected average Joule heating to a tight tolerance, which is enough to keep this case as the current simple EM-to-thermal validation baseline
- adopted first conservative Joule-heating validation gate at `dp=0.05`:
  - `EM_APHI_STAGGERED_SMOKE_REQUIRE_CONVERGED=1`
  - `EM_APHI_STAGGERED_SMOKE_MAX_RESIDUAL_AY=1.2e-5`
  - `EM_APHI_STAGGERED_SMOKE_MAX_AY_MIRROR=1e-6`
  - `EM_APHI_STAGGERED_SMOKE_MAX_JOULE_MIRROR=1e-6`
  - `EM_APHI_STAGGERED_SMOKE_MAX_PHI_CENTERED_ANTIMIRROR=1e-6`
  - `EM_APHI_JOULE_MAX_TEMPERATURE_MIRROR=1e-6`
  - `EM_APHI_JOULE_MAX_TEMPERATURE_DELTA_MIRROR=1e-6`
  - `EM_APHI_JOULE_MAX_TEMPERATURE_DELTA_AVG_ERROR=1e-5`
- validation layering is now:
  - manufactured/reference-driven: `test_3d_em_aphi_matrix_free_staggered_smoke`
  - physical-response/symmetry-driven: `test_3d_em_aphi_matrix_free_physical_response`
  - one-way EM-to-thermal: `test_3d_em_aphi_matrix_free_joule_heating_baseline`

- added a switchable forcing-profile knob to the standalone physical-response and Joule-heating baselines:
  - default `EM_APHI_FORCED_PROFILE=sin` keeps the original smooth distributed drive
  - optional `EM_APHI_FORCED_PROFILE=gaussian` adds a simple localized drive controlled by `EM_APHI_FORCED_CENTER_X` and `EM_APHI_FORCED_WIDTH`
- purpose: let the branch move from globally smooth response checks toward localized hot-spot style response without introducing a heavier geometry yet

- added a localized Gaussian forcing option to the physical-response and Joule-heating baselines and verified the first `dp=0.05` Joule-heating hot-spot response:
  - `forced_profile=gaussian`, `forced_center_x=0.5`, `forced_width=0.12`
  - `converged=1`, `outer_iterations=87`, `residual_ay_l2=9.81e-06`
  - `joule_density_max=1.78e-02`, `temperature_delta_average=1.31e-03`
  - `temperature_delta_average_error=4.64e-06`
  - mirror-symmetry diagnostics for `Ay`, Joule, and temperature remained essentially machine-zero
- interpretation: the branch now resolves both smooth distributed heating and localized hot-spot heating on the same simple one-way EM-to-thermal baseline

## 2026-05-12: Mirrored localized hotspot tracking
- added a standalone target `test_3d_em_aphi_matrix_free_hotspot_tracking`
- purpose: verify not only that a localized Gaussian forcing produces a hot spot, but that the Joule and thermal hot-spot locations move in a physically consistent way when the forcing center is shifted
- first off-center run used `EM_APHI_FORCED_CENTER_X=0.35` with `EM_APHI_FORCED_WIDTH=0.12` at `dp=0.05`
- observed:
  - `source_centroid_x=0.350005`
  - `joule_centroid_x=0.377319`
  - `temperature_delta_centroid_x=0.376314`
  - `joule_center_error=2.73e-02`
  - `temperature_delta_center_error=2.63e-02`
- mirrored run used `EM_APHI_FORCED_CENTER_X=0.65` with the same width and resolution
- observed:
  - `source_centroid_x=0.649995`
  - `joule_centroid_x=0.622681`
  - `temperature_delta_centroid_x=0.623686`
  - `joule_center_error=2.73e-02`
  - `temperature_delta_center_error=2.63e-02`
- interpretation:
  - forcing placement is accurate on both sides
  - Joule and temperature hot spots drift by nearly identical magnitudes in opposite directions
  - the drift is therefore a mirrored physical-response feature of the present EM setup, not random numerical bias

## 2026-05-12: Next geometry-response validation target
- added a new standalone target `test_3d_em_aphi_matrix_free_constriction_heating`
- purpose: move from localized hot-spot tracking in a uniform band to a simple geometry-driven heating response in a necked conductor
- geometry design:
  - left wide shoulder
  - right wide shoulder
  - narrow central constriction
- first diagnostics are ratio-based rather than exact-reference based:
  - `neck_joule_average`
  - `shoulder_joule_average`
  - `neck_temperature_delta_average`
  - `shoulder_temperature_delta_average`
  - `joule_neck_to_shoulder_ratio`
  - `temperature_delta_neck_to_shoulder_ratio`
- intended validation meaning:
  - if the constriction geometry is physically active, Joule heating and temperature rise should be stronger in the neck than in the shoulders under the same broad forcing profile
- optional gate knobs prepared:
  - `EM_APHI_CONSTRICTION_MIN_JOULE_NECK_RATIO`
  - `EM_APHI_CONSTRICTION_MIN_TEMPERATURE_NECK_RATIO`
- first `dp=0.05` constriction-heating run with `EM_APHI_CONSTRICTION_LENGTH=0.20` and `EM_APHI_CONSTRICTION_HEIGHT=0.16` gave:
  - `converged=1`, `outer_iterations=72`, `residual_ay_l2=9.51e-06`
  - `joule_density_max=1.22056e-01`
  - `temperature_delta_average=1.01820e-02`
  - `temperature_delta_average_error=7.02e-06`
  - `neck_joule_average=1.19616e-01`
  - `shoulder_joule_average=9.89355e-02`
  - `neck_temperature_delta_average=1.18637e-02`
  - `shoulder_temperature_delta_average=9.90168e-03`
  - `joule_neck_to_shoulder_ratio=1.209`
  - `temperature_delta_neck_to_shoulder_ratio=1.198`
- interpretation:
  - the constricted geometry already produces a clear and stable heating concentration in the neck region;
  - both electromagnetic power deposition and one-way thermal response are stronger in the neck than in the shoulders by about twenty percent;
  - this is the first geometry-induced heating-concentration validation in the current branch.
- next scan priority:
  - vary `EM_APHI_CONSTRICTION_HEIGHT` while holding forcing and length fixed;
  - confirm that neck-to-shoulder Joule and temperature ratios increase as the neck gets narrower.
- performed a first neck-height scan for `test_3d_em_aphi_matrix_free_constriction_heating` at `dp=0.05` with smooth `sin` forcing:
  - `EM_APHI_CONSTRICTION_HEIGHT=0.20`
    - `points=112`, `neck_particle_count=16`, `shoulder_particle_count=96`
    - `joule_neck_to_shoulder_ratio=1.209033`
    - `temperature_delta_neck_to_shoulder_ratio=1.198151`
  - `EM_APHI_CONSTRICTION_HEIGHT=0.16`
    - `points=112`, `neck_particle_count=16`, `shoulder_particle_count=96`
    - `joule_neck_to_shoulder_ratio=1.209033`
    - `temperature_delta_neck_to_shoulder_ratio=1.198151`
  - `EM_APHI_CONSTRICTION_HEIGHT=0.10`
    - `points=104`, `neck_particle_count=8`, `shoulder_particle_count=96`
    - `joule_neck_to_shoulder_ratio=1.316765`
    - `temperature_delta_neck_to_shoulder_ratio=1.301144`
- interpretation:
  - the coarse-resolution geometry trend is physically consistent: a narrower neck increases both Joule concentration and thermal concentration;
  - however, `0.20` and `0.16` are effectively unresolved as distinct geometries at `dp=0.05`, which is why their ratios are nearly identical and their particle counts match exactly;
  - the first cleanly resolved narrowing appears at `constriction_height=0.10`, where both the particle count in the neck and the neck-to-shoulder heating ratios change materially.
- development implication:
  - the basic geometry-induced heating trend is now supported;
  - the next most valuable check is either a finer-grid repeat of the constriction-height scan or a coarser set of height values separated by more than one particle spacing.

## 2026-05-12: Asymmetric constriction interpretation and centered-Gaussian follow-up
- first asymmetric constriction runs showed a consistent local narrow-side enhancement but did not yet move the global Joule / thermal centroid to the narrow side:
  - at `dp=0.05`, `asym_narrow_height=0.12`, the case produced `joule_narrow_to_wide_ratio ≈ 1.051` and `temperature_delta_narrow_to_wide_ratio ≈ 1.048`, while the global centroids remained left-shifted because the wide side still carries more total weight
  - at `dp=0.025`, `asym_narrow_height=0.08`, the local enhancement strengthened to `joule_narrow_to_wide_ratio ≈ 1.084` and `temperature_delta_narrow_to_wide_ratio ≈ 1.077`, which confirms the narrow-side local concentration trend more clearly under stronger asymmetry and finer resolution
- interpretation:
  - the current asymmetric case is validating local geometry sensitivity, not yet a full global hotspot migration to the narrow side
  - with broad `sin` forcing, the source itself carries a left-biased geometric weighting, so absolute centroid shifts alone are not a clean geometry-only diagnostic
- upgraded `test_3d_em_aphi_matrix_free_asymmetric_constriction_heating` to report:
  - `joule_relative_to_source_shift`
  - `temperature_delta_relative_to_source_shift`
- next validation step for this case:
  - rerun the asymmetric geometry with a centered Gaussian forcing profile so source bias is reduced and source-relative response drift becomes easier to interpret
- development implication:
  - this is the last lightweight geometry-bias check before moving to a simplified matrix-free TEAM7-like three-region case with explicit air / conductor / coil partitioning

## 2026-05-12: First matrix-free TEAM7-like three-region target
- added a new standalone target `test_3d_em_aphi_matrix_free_team7_like`
- purpose of this step:
  - stop iterating only on strip-style geometries
  - carry the stabilized matrix-free `A-phi -> E/J -> Joule -> one-way thermal` chain into a first TEAM7-like spatial layout
- current scope of the new target:
  - a single box body discretized with the existing matrix-free route
  - runtime partition into `air`, `coil`, and `conductor` regions using box fractions inspired by the existing eigen-based TEAM7-like case
  - explicit source applied only inside the coil region
  - region-wise summaries for `A`, `phi`, `E`, `J`, Joule heat, and temperature rise
  - source-centroid and conductor-response-centroid diagnostics so conductor heating drift can be interpreted relative to the coil source
- this is intentionally still lighter than the final TEAM7 workflow:
  - no STL geometry
  - no full coil/plate/air multi-body particle-generation pipeline
  - no direct benchmark curve comparison yet
- development role:
  - bridge from the simplified strip/constriction validation family to a first matrix-free three-region TEAM7-like layout before touching the full `particle_generation_em` scaffold


## 2026-05-13: TEAM7-like matrix-free instability localized to conductivity-jump handling
- aligned the matrix-free and eigen TEAM7-like comparisons enough to separate geometry/source issues from matrix-free coupling issues:
  - the eigen single-body three-region TEAM7-like reference solves cleanly on the same geometry, region partition, frequency, and coil-profile shape;
  - the matrix-free TEAM7-like route diverges on the corresponding forced-response case, so the problem is not the TEAM7-like spatial layout by itself.
- enabled additional case modes in `test_3d_em_aphi_matrix_free_team7_like` so the same three-region geometry could be exercised with built-in manufactured/reference modes:
  - `coulomb_variable_sigma_source_free`
  - `coulomb_variable_sigma_driven`
- first observation from those manufactured/reference runs:
  - with heterogeneous conductivity (`sigma_air=1e-4`, `sigma_conductor=1`, `sigma_coil=1e-4` or `1`), the matrix-free route also becomes unstable even without the explicit TEAM7-like impressed-coil forcing;
  - therefore the instability is not limited to the external coil forcing path.
- continuation was ruled out as the primary cause:
  - disabling source ramp and staged load stepping (`effective_source_scale=1`, `load_steps=1`) did not restore stability for the heterogeneous three-region reference cases.
- constant-conductivity control runs on the same TEAM7-like geometry were then used to isolate the jump effect:
  - with `sigma_air=sigma_conductor=sigma_coil=1`, both `coulomb_variable_sigma_source_free` and `coulomb_variable_sigma_driven` remain bounded and produce small residuals / reasonable field errors after long outer iteration runs;
  - this shows the geometry itself is not enough to trigger the blow-up.
- conductivity-contrast scan for `coulomb_variable_sigma_source_free` with `sigma_conductor=sigma_coil=1` and varying `sigma_air`:
  - `sigma_air=1e-1` (`10:1` jump): bounded and near-stationary
    - `outer_iterations=400`
    - `residual_ay_l2=1.51e-05`
    - `phi_l2_error=7.97e-06`
    - `current_l2_error=3.46e-06`
  - `sigma_air=5e-2` (`20:1` jump): still bounded but noticeably degraded
    - `outer_iterations=400`
    - `residual_ay_l2=2.47e-03`
    - `phi_l2_error=8.87e-02`
    - `current_l2_error=4.47e-02`
  - `sigma_air=4e-2` (`25:1` jump): no longer trustworthy / effectively unstable
    - `outer_iterations=400`
    - `residual_ay_l2=3.96e+06`
    - `phi_l2_error=2.96e+08`
    - `current_l2_error=2.35e+08`
  - `sigma_air=3.5e-2` (`28.6:1` jump): catastrophic growth
    - `outer_iterations=400`
    - huge field norms and thermal blow-up
  - `sigma_air=3e-2` (`33:1` jump): catastrophic growth by `outer_iterations=297`
  - `sigma_air=2e-2` (`50:1` jump): catastrophic growth by `outer_iterations=146`
  - `sigma_air=1e-2`, `1e-3`, `1e-4`: progressively earlier catastrophic growth
- interpretation:
  - the current matrix-free staggered TEAM7-like route has a conductivity-contrast stability boundary in the neighborhood of roughly `20:1` to `25:1` on this geometry and resolution;
  - that behavior is consistent with an interface-coupling inconsistency rather than a generic failure of the matrix-free geometry handling.
- current hypothesis from code inspection:
  - the diffusion/Laplace part is assembled through pairwise edge contributions using harmonic averaging of the diffusion coefficient, which is interface-aware;
  - the `A-phi` coupling terms are still evaluated in a nodewise form (`sigma[i] * grad(phi)` and `grad(sigma*A)` built from pointwise `sigma*A` fields), which is not obviously flux-consistent with the harmonic-mean edge treatment used in the Laplace part;
  - this mismatch is likely harmless for constant or mild-contrast conductivity but becomes unstable once the conductivity jump is large enough.
- next implementation direction:
  - add interface diagnostics that isolate large-jump edges and quantify their contribution to `A`/`phi` residuals;
  - prototype a face-consistent / edge-consistent treatment for the conductivity-coupling terms so the interface discretization is aligned with the harmonic-mean diffusion operator.
