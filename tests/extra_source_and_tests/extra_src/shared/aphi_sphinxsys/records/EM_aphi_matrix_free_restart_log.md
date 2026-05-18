# Laplace A-phi Matrix-Free Restart Log

## Goal

This log restarts the Laplace-structured frequency-domain `A-phi` work from a clean baseline.
The first implementation target is not the final coupled electromagnetic heating case.
It is a matrix-free prototype that stays close to existing SPHinXsys diffusion / relaxation patterns
while keeping the current Eigen-based Laplace prototype as a reference.

## Working Decisions

1. Phase 1 uses a staggered matrix-free solve, not a local `4x4` coupled solve.
2. Phase 1 keeps the current dual-real-field storage style for particle states.
3. Eigen complex types remain allowed for local algebra and reference solves.
4. Global `Eigen::SparseMatrix<std::complex<Real>>` is kept only as a reference path.
5. Gauge projection is treated as a required correction step, not an optional post-process.

## Test Strategy

We should not start with tiny unit tests only, and we should not jump directly to TEAM7 heating.

Recommended split:

- Use standalone example-style verification cases first under `tests/extra_source_and_tests/3d_examples`.
- After the operators stabilize, move the smallest invariant checks into tighter regression or unit-style tests.
- Keep geometry-heavy, interface-heavy, and multi-resolution electromagnetic verification as standalone cases.

This means the first phase is driven by dedicated test cases, not by trying to force everything into the smallest unit-test format.

## Phase Order

### Phase 1

- Complex scalar matrix-free Helmholtz
- Complex scalar matrix-free Poisson / diffusion
- Pairwise gradient / divergence helpers
- Gauge projection
- Staggered `A-step -> phi-step -> gauge projection`

### Phase 2

- Optional local `4x4` complex block correction
- Better preconditioning / acceleration
- Stronger diagnostics and operator comparisons

### Phase 3

- Reconnect to larger TEAM-like electromagnetic heating validation
- Revisit where small Eigen dense blocks are helpful
- Keep global sparse Eigen only as a reference or debugging route

## File Plan

The first-stage code should stay in `tests/extra_source_and_tests/extra_src/shared/`.

Planned files:

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

## Extraction Map

### From the current Laplace Eigen prototype

Source:

- `extra_src/shared/legacy_aphi_archive/baselines/electromagnetic_aphi_laplace_eigen.h`
- `extra_src/shared/legacy_aphi_archive/baselines/electromagnetic_aphi_laplace_eigen.hpp`

Extract or adapt first:

- field and diagnostics layout
- harmonic mean helper
- pairwise diffusion weight helper
- pairwise gradient accumulation logic
- pairwise Laplace accumulation logic
- `divA`, `divJ`, and Joule diagnostics definitions

### From the current frequency-domain local relaxation code

Source:

- `extra_src/shared/aphi_case_support/electromagnetic_team7_aphi_frequency_dynamics.h`
- `extra_src/shared/aphi_case_support/electromagnetic_team7_aphi_frequency_dynamics.hpp`

Reuse the implementation style from:

- `ElectricPotentialSourceFromVectorFieldInner`
- `ElectricPotentialSourceFromVectorFieldContact`
- `ScalarRelaxationInnerByName`
- `ScalarRelaxationContactByName`
- `ScalarRelaxationComplexByName`
- `UpdateScalarByRelaxationRateByName`
- `VectorPotentialFrequencyCoupledEquationComplex`
- `VectorPotentialFrequencyCoupledBlockEquationComplex`

### From the SPHinXsys diffusion-splitting pattern

Source:

- `src/shared/particle_dynamics/diffusion_optimization_dynamics/diffusion_splitting_base.h`
- `src/shared/particle_dynamics/diffusion_optimization_dynamics/diffusion_splitting_state.h`
- `src/shared/particle_dynamics/diffusion_optimization_dynamics/diffusion_splitting_state.hpp`

Reuse the algorithmic pattern:

- assemble local residual from pairwise interactions
- estimate a local diagonal-like scaling
- apply local correction
- iterate and report residual

## First Two Test Cases

### Test 1

Name:

- `test_3d_em_aphi_matrix_free_complex_helmholtz`

Purpose:

- verify complex pairwise Laplace plus reaction solve
- verify matrix-free residual correction for a single complex scalar field
- compare against a manufactured analytic field and the current Eigen reference path

### Test 2

Name:

- `test_3d_em_aphi_matrix_free_gauge_projection`

Purpose:

- construct `A = A0 + grad(psi)`
- solve `Laplace(chi) = div(A)`
- update `A <- A - grad(chi)` and `phi <- phi + i omega chi`
- verify `||divA||` decreases while `E` remains nearly unchanged

## Code-Level Task Breakdown

### Task A

Create shared complex types and lightweight containers.

Deliverables:

- `electromagnetic_aphi_matrix_free_types.h`
- `electromagnetic_aphi_matrix_free_fields.h`

Status:

- completed

### Task B

Create reusable operator helpers that do not depend on global Eigen sparse assembly.

Deliverables:

- `electromagnetic_aphi_matrix_free_operators.h`
- `electromagnetic_aphi_matrix_free_operators.hpp`

Scope:

- harmonic mean
- pairwise diffusion weight
- pairwise scalar Laplace contribution
- pairwise scalar gradient contribution
- simple complex norm helpers

Status:

- in progress

### Task C

Create residual containers and evaluators for the scalar Helmholtz prototype.

Deliverables:

- `electromagnetic_aphi_matrix_free_residuals.h`
- `electromagnetic_aphi_matrix_free_residuals.hpp`

### Task D

Create a scalar complex matrix-free relaxation solver and use it for the Helmholtz verification case.

Deliverables:

- `electromagnetic_aphi_matrix_free_solver.h`
- `electromagnetic_aphi_matrix_free_solver.hpp`
- first runnable Helmholtz test case

### Task E

Create the gauge projection module and build the second verification case.

Deliverables:

- `electromagnetic_aphi_matrix_free_gauge_projection.h`
- `electromagnetic_aphi_matrix_free_gauge_projection.hpp`
- first runnable gauge projection test case

## Work Log

### 2026-05-10

- restarted the Laplace `A-phi` matrix-free branch with a separate implementation log
- fixed the earlier adaptive cell-linked-list workaround so it stays in `extra_src/shared` instead of `src/`
- created the first shared matrix-free files:
  - `electromagnetic_aphi_matrix_free_types.h`
  - `electromagnetic_aphi_matrix_free_fields.h`
- clarified that the first implementation step is the staggered single-field route, not the local coupled route
- clarified that the first validation step should be standalone example-style verification cases
- started Task B by introducing operator-helper scaffolding for matrix-free reuse
