# EM A-phi Matrix-Free To SphinxSys SYCL Refactor Plan

## Goal

Refactor the current electromagnetic `A-phi` matrix-free prototype into a `SphinxSys`-style solver that:

- follows the existing framework organization instead of living mainly as test-only prototype code
- supports frequency-domain electromagnetic iteration in the same spirit as other `SphinxSys` dynamics
- can be migrated to `SYCL` and run on GPU using the existing `src/src_sycl` execution model
- scales toward a real three-body TEAM7-style setup instead of a single box body with region flags
- exposes iteration history, residuals, pseudo-step progress, and diagnostics in a way that matches `SphinxSys` usage patterns

This document is a migration and redesign plan, not a statement that the current code is already in the desired final form.

## Current State

### What is working now

- The matrix-free electromagnetic prototype has already passed an important numerical checkpoint:
  - the interface handling for `sigma grad(phi)` and `div(sigma A)` was corrected to use a harmonic-weighted edge-consistent form
  - the manufactured/reference `phi` source construction in the TEAM7-like test driver was synchronized with the solver operator
- With those fixes, the following high-contrast cases are now numerically healthy:
  - `coulomb_variable_sigma_source_free`
  - `coulomb_variable_sigma_forced_response` with `EM_APHI_FORCED_PROFILE=eigen_z_sin`
  - `coulomb_variable_sigma_forced_response` with `EM_APHI_FORCED_PROFILE=sin`
- The `source_free` case now converges with very small field and Joule errors when under-relaxation and load stepping are enabled.

### What is not acceptable yet

- The `driven` case is still too slow and not robust enough under aggressive conductivity contrast.
- The current implementation is still a prototype in `tests/extra_source_and_tests/extra_src/shared/`.
- The current control flow is host-side `StdVec` based and not organized around `SphinxSys` `Dynamics` / `Relation` / `execution policy` patterns.
- The current matrix-free loop is not `SYCL` ready:
  - it uses `std::complex`
  - it uses `Eigen::Matrix<Complex, Dimensions, 1>` as `Vec3c`
  - it relies on host vectors and repeated temporary allocations inside iterative loops
- The current solver behaves like a standalone fixed-point experiment, not like a formal `SphinxSys` solver object with pseudo-step execution and residual history output.

### Why the current iterative version can be slower than Eigen

Right now, the matrix-free approach is not yet benefiting from the things that usually make iterative GPU solvers worthwhile:

- no device execution
- no device-friendly data layout
- repeated host temporary vector creation
- serial outer fixed-point orchestration
- no serious preconditioner or Krylov accelerator
- multiple load steps and outer under-relaxation added only to recover stability

So the current state is:

- numerically promising
- architecturally immature
- not yet performance-competitive

That is normal for a prototype, but it means we should stop treating the current test driver as the final foundation.

## Existing Prototype Files

### Matrix-free prototype core

Current files under `tests/extra_source_and_tests/extra_src/shared/`:

- `electromagnetic_aphi_matrix_free_types.h`
- `electromagnetic_aphi_matrix_free_fields.h`
- `electromagnetic_aphi_matrix_free_pairwise_graph.h/.hpp`
- `electromagnetic_aphi_matrix_free_operators.h/.hpp`
- `electromagnetic_aphi_matrix_free_residuals.h/.hpp`
- `electromagnetic_aphi_matrix_free_solver.h/.hpp`
- `electromagnetic_aphi_matrix_free_aphi_residuals.h/.hpp`
- `electromagnetic_aphi_matrix_free_aphi_solver.h/.hpp`
- `electromagnetic_aphi_matrix_free_gauge_projection.h/.hpp`

### Existing `SphinxSys`-style electromagnetic dynamics already worth reusing

- `electromagnetic_team7_aphi_dynamics.h/.hpp`
- `electromagnetic_team7_aphi_frequency_dynamics.h/.hpp`

These files are much closer to the framework style we want long term.

### Existing `SYCL` execution templates already in the repository

Under `src/src_sycl/shared/` there are already reusable patterns:

- `particle_dynamics/particle_iterators_sycl.h`
- `particle_dynamics/implementation_sycl.h`
- `common/sphinxsys_variable_sycl.hpp`
- `common/sphinxsys_variable_array_sycl.hpp`
- `common/sphinxsys_buffer_array_sycl.hpp`

These should be treated as the target style, not bypassed.

## Main Diagnosis

The current prototype has three separate problems:

### Problem 1: Numerical kernel and architecture are mixed together

The current code mixes:

- graph construction
- operator application
- residual evaluation
- fixed-point orchestration
- manufactured reference logic
- exploratory diagnostics

This makes it hard to:

- profile performance
- switch execution backend
- expose per-iteration monitoring
- reuse the solver for a real multi-body case

### Problem 2: The current math types are not good GPU boundary types

Current types:

- `Complex = std::complex<Real>`
- `Vec3c = Eigen::Matrix<Complex, Dimensions, 1>`

These are convenient on CPU, but they are not good long-term `SYCL` kernel payload types. Even if some compilers can digest them, they are the wrong foundation for maintainable device code.

### Problem 3: The solver is not expressed in SphinxSys dynamics style

The current solver entry point:

- `solveMatrixFreeAPhiStaggered(...)`

is a monolithic host function that runs until convergence inside one call.

That is unlike typical `SphinxSys` usage, where we prefer:

- explicit state variables on bodies
- clear local/contact dynamics
- policy-driven iteration
- reusable `exec()` calls
- observable iteration progress

## Target End State

The desired end state should look like this:

### Layer 1: Physics state on bodies

Body variables should carry the electromagnetic state in a framework-native way:

- vector potential real/imag
- scalar potential real/imag
- current density real/imag
- electric field real/imag
- conductivity
- magnetic reluctivity
- source current density real/imag
- Joule heat source

This should reuse the philosophy already present in:

- `electromagnetic_team7_aphi_frequency_dynamics.h/.hpp`

### Layer 2: Solver state object

Introduce a dedicated solver-state structure with fields like:

- current pseudo-iteration index
- current pseudo-step size or relaxation factor
- max residuals
- divergence norms
- relative update norm
- convergence flags
- load-step index
- source-ramp scale

This should be persistent across calls instead of living only inside one monolithic function.

### Layer 3: Dynamics-based pseudo-iteration

The electromagnetic frequency solve should be decomposed into `SphinxSys`-style operations such as:

- build or refresh graph connectivity
- compute `sigma grad(phi)`
- compute `div(sigma A)`
- update `A_x`
- update `A_y`
- update `A_z`
- update `phi`
- apply reference / gauge / boundary conditions
- compute residuals
- check convergence

These can be composed in a solver driver class or a small orchestration module, but each step should be an identifiable operation.

### Layer 4: SYCL backend

The same operator layer should have:

- a CPU implementation first
- a `SYCL` execution path second

The `SYCL` path should not be a separate mathematical algorithm. It should be the same operator sequence on a different execution backend.

### Layer 5: Three-body TEAM7 setup

The final TEAM7-style case should be built from:

- air body
- coil body
- conductor body

with proper:

- body relations
- contact relations
- source assignment
- Joule coupling
- diagnostics per physical body

instead of one box with region flags.

## Recommended Refactor Strategy

Do this in stages. Do not jump directly from the current prototype to a full GPU implementation.

## Stage 0: Freeze The Prototype As Reference

### Objective

Preserve a known-good numerical reference before major refactor work.

### Actions

- Keep the current test-only matrix-free implementation untouched except for critical bug fixes.
- Record known-good benchmark outputs for:
  - `coulomb_variable_sigma_source_free`
  - `coulomb_variable_sigma_forced_response` with `eigen_z_sin`
  - `coulomb_variable_sigma_forced_response` with `sin`
- Record one known-problem benchmark:
  - `coulomb_variable_sigma_driven`

### Deliverables

- a small benchmark log table in markdown
- a statement of which prototype version is the numeric baseline

### Why this matters

Once the architecture changes, we need a stable truth set to confirm we did not break the corrected operator logic.

## Stage 1: Split Math Kernels From Test Driver Logic

### Objective

Reduce coupling between the TEAM7-like example and the matrix-free solver internals.

### Actions

Move example-specific responsibilities out of the solver layer:

- manufactured solution building
- region summaries
- thermal-only exploratory reporting
- exploratory load-step scripting

Keep in the solver layer only:

- field update kernels
- residual kernels
- graph operators
- gauge/reference enforcement
- solver-state updates

### Suggested file split

Keep the test example:

- `test_3d_em_aphi_matrix_free_team7_like.cpp`

focused on:

- geometry
- material setup
- source setup
- run configuration
- reporting

Keep solver internals focused on:

- apply operator
- update field
- compute residual
- convergence check

### Deliverables

- a thinner test driver
- a cleaner operator/solver separation

## Stage 2: Replace Device-Unfriendly Complex Types

### Objective

Prepare the code for `SYCL` kernels by removing fragile host-centric complex container types from kernel-facing code.

### Current problem

The current type layer is:

- `Complex = std::complex<Real>`
- `Vec3c = Eigen::Matrix<Complex, Dimensions, 1>`

### Recommendation

Introduce explicit POD-friendly types for device-facing kernels, for example:

- `ComplexValue`
- `ComplexVec3`

Suggested shape:

```cpp
struct ComplexValue
{
    Real real_;
    Real imag_;
};

struct ComplexVec3
{
    ComplexValue x_, y_, z_;
};
```

Then provide:

- helper arithmetic
- conjugate
- norm
- fused multiply-add style helpers

### Important rule

Do not replace all host-side `std::complex` and Eigen usage at once.

Instead:

- keep host reference and testing utilities simple if needed
- migrate the operator layer first to device-friendly types

### Deliverables

- a new device-facing complex type layer
- a narrow conversion boundary between host reference math and solver kernel math

## Stage 3: Introduce A Formal Electromagnetic Solver State

### Objective

Turn the monolithic host solve into an observable iterative process.

### Add a new state structure

Create a persistent state object with at least:

- `outer_iteration`
- `load_step`
- `effective_source_scale`
- `effective_load_scale`
- `field_update_l2`
- `relative_field_update_l2`
- `residual_ax_l2`
- `residual_ay_l2`
- `residual_az_l2`
- `residual_phi_l2`
- `divergence_a_l2`
- `divergence_j_l2`
- `converged`
- `stagnated`
- `blow_up_detected`

### Add iteration history

Do not only print final diagnostics.

Maintain a history buffer or sampled history with entries like:

- iteration index
- residual summary
- update norm
- divergence summary
- source scale
- load scale

### Why this matters

This is the missing equivalent of “time-step style” introspection for a frequency-domain iterative solve.

Even though the physics is steady-state, the solver still has a pseudo-time or pseudo-iteration evolution that should be visible.

### Deliverables

- per-iteration solver monitoring
- a reusable convergence/debugging record

## Stage 4: Recast The Solve As SphinxSys-Style Dynamics

### Objective

Make the electromagnetic solve feel like a `SphinxSys` solver instead of a standalone algorithm blob.

### Proposed design

Introduce classes along these lines:

- `InitializeMatrixFreeAphiVariables`
- `BuildMatrixFreeAphiGraph`
- `ComputeMatrixFreeAphiResiduals`
- `UpdateMatrixFreeAphiVectorPotential`
- `UpdateMatrixFreeAphiScalarPotential`
- `ApplyMatrixFreeAphiGaugeOrReference`
- `CheckMatrixFreeAphiConvergence`
- `MatrixFreeAphiPseudoTimeStep`
- `MatrixFreeAphiFrequencySolver`

### Style target

The top-level solve should become something like:

```cpp
while (!solver_state.converged && solver_state.iteration < max_iterations)
{
    compute_rhs.exec();
    update_a.exec();
    update_phi.exec();
    apply_constraints.exec();
    evaluate_residuals.exec();
    record_history.exec();
}
```

This is much closer to the way `SphinxSys` users expect to read solver logic.

### Deliverables

- a formal solver driver
- identifiable sub-dynamics
- clearer profiling boundaries

## Stage 5: Reduce Host-Side Temporary Allocation

### Objective

Remove the current per-iteration host memory churn.

### Current issue

Inside the current `solveMatrixFreeAPhiStaggered(...)`, each outer iteration builds multiple temporary vectors:

- `effective_sources`
- `grad_phi`
- `sigma_grad_phi`
- `gauge_penalty_gradient`
- `rhs_ax`
- `rhs_ay`
- `rhs_az`
- `divergence_sigma_a`
- `rhs_phi`

This is expensive even before any GPU discussion.

### Recommendation

Introduce a workspace/cache object:

- `MatrixFreeAphiWorkspace`

with reusable buffers:

- temporary RHS arrays
- temporary gradients
- temporary divergence arrays
- temporary residual arrays

### Deliverables

- stable scratch storage
- much better CPU performance
- an easier path to device memory ownership

## Stage 6: Introduce Multi-Body TEAM7 Structure

### Objective

Move from a single box with region tags to a real physical decomposition.

### Target bodies

- air body
- coil body
- conductor body

### What changes

Instead of:

- one particle set with `is_coil`, `is_conductor`, `sigma[i]`, `nu[i]`

Move to:

- separate `SPHBody` instances
- body-specific material assignment
- contact relations for electromagnetic coupling
- body-local source assignment for the coil

### Why this matters

The single-body prototype was useful for debugging operator consistency.
It is not the right shape for the final TEAM7-like physical workflow.

### Deliverables

- three-body case layout
- clearer boundary and contact handling
- easier extension toward realistic TEAM7 workflows

## Stage 7: CPU SphinxSys-Style Solver Before SYCL

### Objective

First stabilize the architecture in framework style on CPU.

### Why this stage is mandatory

If we try to jump directly from:

- prototype host vectors

to:

- production `SYCL`

we will conflate:

- math bugs
- architecture bugs
- device bugs
- data movement bugs

That will slow development down more than it speeds it up.

### Exit criterion for this stage

The CPU framework-style solver should:

- reproduce the corrected prototype results
- expose iteration history
- run the three-body case layout
- support stable convergence tuning

Only then should we freeze the CPU architecture and port the hot kernels to `SYCL`.

## Stage 8: Introduce SYCL Data Ownership

### Objective

Map the electromagnetic state and graph data onto the existing repository `SYCL` memory model.

### Use existing repository templates

Reuse the patterns already present in:

- `src/src_sycl/shared/particle_dynamics/implementation_sycl.h`
- `src/src_sycl/shared/particle_dynamics/particle_iterators_sycl.h`
- `src/src_sycl/shared/common/sphinxsys_variable_sycl.hpp`
- `src/src_sycl/shared/common/sphinxsys_variable_array_sycl.hpp`
- `src/src_sycl/shared/common/sphinxsys_buffer_array_sycl.hpp`

### New device-side data that must be represented explicitly

- device field arrays for `A_real`, `A_imag`, `phi_real`, `phi_imag`
- device conductivity and reluctivity arrays
- device source arrays
- device graph edge arrays
- device workspace buffers
- device residual buffers

### Graph layout recommendation

Do not port the graph as pointer-heavy CPU objects.

Use a flat structure-of-arrays representation such as:

- `edge_i`
- `edge_j`
- `pair_weight`
- optional harmonic coefficients

or a CSR-like neighborhood form if that matches other neighbor structures better.

### Deliverables

- device-resident electromagnetic state
- device-resident graph
- explicit host-device synchronization boundaries

## Stage 9: Port Hot Kernels To SYCL

### Objective

Port only the expensive inner kernels, not the entire application logic at once.

### First kernels to port

1. graph gradient application
2. harmonic-weighted gradient application
3. divergence of `sigma A`
4. A-component Jacobi or relaxed update
5. phi update
6. residual evaluation
7. reductions for norms

### Style target

Each hot operator should use `particle_for`, `particle_reduce`, or a graph-edge-specific analogous pattern.

### Important recommendation

Do not keep a host-only outer iteration with repeated host-device transfers for every kernel.

Instead:

- keep full solver state on device as much as possible
- reduce only the scalar convergence diagnostics needed to decide whether to continue

### Deliverables

- a device-enabled matrix-free electromagnetic operator layer
- minimal synchronization overhead

## Stage 10: Replace Fixed-Point With A Faster Iterative Strategy

### Objective

Make the iterative approach actually competitive.

### Current limitation

The present staggered fixed-point scheme is numerically useful for prototyping, but it is not obviously the fastest final algorithm.

### Recommended options

Evaluate one or more of these after the architecture is cleaned up:

- block Gauss-Seidel style A-phi iteration
- Anderson acceleration on the outer fixed-point map
- Krylov solve on the coupled linearized system
- better diagonal or block preconditioning
- adaptive pseudo-step control

### Practical recommendation

Do not introduce a sophisticated accelerator before the architecture cleanup.

First make the operator framework clean.
Then add acceleration in a measurable way.

## Stage 11: Add Solver-Style Logging And Output

### Objective

Make the solve feel like a real `SphinxSys` iterative process.

### Add outputs

- per-iteration residual history
- sampled convergence CSV or simple text record
- load-step history
- source-ramp history
- final solver summary

### Optional but recommended

Add a verbosity control with levels:

- `quiet`
- `summary`
- `iteration`
- `debug`

### Why this matters

The user concern is correct:

The current solver does not feel like a `SphinxSys` iterative process because it does not expose its progress the way a pseudo-time solver should.

## Stage 12: Final TEAM7 Three-Body Production Case

### Objective

Use the cleaned CPU/SYCL framework to build the actual three-body TEAM7-style example.

### Final example should include

- three bodies with explicit relations
- coil source specification
- conductor Joule heating
- body-wise diagnostics
- optional thermal coupling
- optional gauge handling
- optional regression thresholds

### Acceptance criteria

- numerically stable at target contrast
- converges in reasonable wall time
- CPU and SYCL results agree within tolerances
- solver history shows predictable convergence behavior
- code layout matches framework conventions

## Recommended New File Organization

### Keep prototype code temporarily

Keep current prototype in `tests/extra_source_and_tests/extra_src/shared/` as a numeric reference during refactor.

### Introduce formal CPU module

Add a new formal module under `src/shared/electromagnetics/` or the closest repository-consistent location, for example:

- `matrix_free_aphi_base.h/.hpp`
- `matrix_free_aphi_graph.h/.hpp`
- `matrix_free_aphi_operator.h/.hpp`
- `matrix_free_aphi_residual.h/.hpp`
- `matrix_free_aphi_solver.h/.hpp`
- `matrix_free_aphi_dynamics.h/.hpp`

### Add SYCL counterparts

Under `src/src_sycl/shared/electromagnetics/`, add:

- `matrix_free_aphi_operator_sycl.h/.hpp`
- `matrix_free_aphi_residual_sycl.h/.hpp`
- `matrix_free_aphi_solver_sycl.h/.hpp`
- optional shared device type headers

### Keep tests focused on integration

The example under:

- `tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_matrix_free_team7_like/`

should eventually become only:

- case setup
- driver logic
- reporting

not the primary location of solver experimentation.

## Short-Term Priorities

If implementation begins now, the most efficient order is:

1. Freeze current numeric baseline.
2. Split example logic from solver core.
3. Introduce device-friendly complex/vector types.
4. Add persistent solver state and iteration history.
5. Reduce per-iteration host allocations with a workspace object.
6. Re-express the CPU solve in `SphinxSys` dynamics style.
7. Build the three-body TEAM7 case on CPU.
8. Port hot kernels to `SYCL`.
9. Add acceleration only after the `SYCL` backend is structurally sound.

## Phase 1 Implementation Checklist

This section turns the refactor strategy into the first concrete work package. The intent of Phase 1 is not to finish the final solver. The intent is to stop further architectural drift and move the current prototype into a shape that can be cleanly extended.

### Phase 1 Goal

At the end of Phase 1, we should have:

- the same corrected numerical operators as the current prototype
- a thinner TEAM7-like test driver
- a reusable solver core with explicit state and workspace
- per-iteration monitoring
- no major behavior change yet in physics or convergence logic

### Phase 1 Non-Goals

Phase 1 should not try to do these yet:

- full three-body decomposition
- full `SYCL` port
- solver acceleration redesign
- full replacement of every host-side complex helper
- broad migration into final `src/shared` and `src/src_sycl` locations

### Phase 1 Suggested Order

1. freeze the current benchmark outputs in the plan file or a sibling markdown file
2. introduce explicit solver state and iteration history types
3. introduce a reusable workspace object for temporary arrays
4. move reporting and manufactured-reference helpers out of the solver path
5. refactor the main staggered solve into smaller internal steps
6. update the TEAM7-like driver to use the thinner solver API

## Phase 1 File-Level Change Plan

### 1. Benchmark Baseline

Before touching solver internals, record the known-good outputs for these cases:

- `coulomb_variable_sigma_source_free`
- `coulomb_variable_sigma_driven`
- `coulomb_variable_sigma_forced_response` with `EM_APHI_FORCED_PROFILE=eigen_z_sin`
- `coulomb_variable_sigma_forced_response` with `EM_APHI_FORCED_PROFILE=sin`

Record at least:

- `converged`
- `completed_load_steps`
- `outer_iterations`
- `residual_ax_l2`
- `residual_ay_l2`
- `residual_phi_l2`
- `divergence_a_error_l2`
- `ax_l2_error`
- `ay_l2_error`
- `phi_l2_error`
- `electric_l2_error`
- `current_l2_error`
- `joule_l2_error`

This can stay in markdown at first. It does not need a new parser or regression harness yet.

### 2. `electromagnetic_aphi_matrix_free_types.h`

This file should become the home for solver-facing data carriers, not only scalar aliases.

Add or prepare:

- `MatrixFreeAPhiIterationRecord`
- `MatrixFreeAPhiIterationHistory`
- `MatrixFreeAPhiWorkspaceSizes` or equivalent sizing helper

The immediate purpose is not device migration yet. The immediate purpose is to stop scattering iteration bookkeeping and scratch-buffer assumptions across solver code.

Recommended additions:

```cpp
struct MatrixFreeAPhiIterationRecord
{
    size_t outer_iteration_ = 0;
    size_t load_step_ = 0;
    Real effective_source_scale_ = 1.0;
    Real effective_load_scale_ = 1.0;
    Real field_update_l2_ = 0.0;
    Real relative_field_update_l2_ = 0.0;
    Real residual_ax_l2_ = 0.0;
    Real residual_ay_l2_ = 0.0;
    Real residual_az_l2_ = 0.0;
    Real residual_phi_l2_ = 0.0;
    Real divergence_a_l2_ = 0.0;
    Real divergence_j_l2_ = 0.0;
};
```

Do not over-design this yet. A flat POD-style history record is enough.

### 3. `electromagnetic_aphi_matrix_free_solver.h`

This header should stop being only a function declaration surface and start describing the persistent solver-side objects.

Add:

- `MatrixFreeAPhiWorkspace`
- `MatrixFreeAPhiMonitorOptions`
- optional `MatrixFreeAPhiConvergenceStatus`

`MatrixFreeAPhiWorkspace` should own reusable arrays such as:

- effective source buffers
- `grad_phi`
- `sigma_grad_phi`
- gauge penalty gradient
- RHS arrays for `ax`, `ay`, `az`, `phi`
- divergence buffers
- previous field copies if update norms need them

The point is to allocate once per solve, not once per iteration.

### 4. `electromagnetic_aphi_matrix_free_solver.hpp`

This file needs the biggest structural cleanup.

Refactor the monolithic solve into internal helper stages such as:

- `initializeMatrixFreeAPhiWorkspace(...)`
- `buildEffectiveSources(...)`
- `buildVectorPotentialRhs(...)`
- `buildScalarPotentialRhs(...)`
- `updateVectorPotentialComponents(...)`
- `updateScalarPotential(...)`
- `evaluateAndStoreResiduals(...)`
- `appendIterationHistory(...)`
- `checkMatrixFreeAPhiConvergence(...)`

Important guideline:

Keep the public solve entry point for now if that minimizes disruption, but make it a thin orchestrator instead of a giant all-in-one body.

The first internal cleanup target should be the current outer-iteration loop in `solveMatrixFreeAPhiStaggered(...)`. That loop should become readable enough that a future `Dynamics` wrapper could call the same stages one pseudo-step at a time.

### 5. `electromagnetic_aphi_matrix_free_residuals.h/.hpp`

These files should become the single place responsible for residual summaries and convergence metrics.

Move or centralize here:

- relative update norm logic
- divergence norm logic
- high-contrast edge diagnostics
- any per-iteration residual aggregation

The solver should ask the residual layer for a structured summary instead of manually stitching many scalar diagnostics together in several locations.

### 6. `electromagnetic_aphi_matrix_free_aphi_solver.h/.hpp`

Use these files as the staging ground for the future `A-phi`-specific orchestration layer, but do not mix back in test-only reporting.

This layer should eventually expose:

- one-step update functions
- one full pseudo-iteration function
- a top-level solve function that uses persistent state and workspace

If there are still utility routines here that belong to the TEAM7-like example rather than the solver, move them out.

### 7. `test_3d_em_aphi_matrix_free_team7_like.cpp`

This file is currently doing too much.

Phase 1 should reduce it to four categories:

- case setup
- material/source assignment
- solver execution
- postprocess/reporting

Move out of the main flow where possible:

- manufactured-reference construction helpers that can live in a separate example helper
- repeated metric helper logic that is not specific to this exact case
- exploratory solver internals that belong in solver history rather than final reporting

The test should stop knowing too much about transient scratch fields inside the solver.

## Phase 1 Acceptance Criteria

Phase 1 should be considered complete only if all of these are true:

- the known-good benchmark cases still match the current prototype within small tolerance
- the solver allocates its large scratch arrays outside the inner outer-iteration loop
- iteration history can be enabled without editing the solver core
- the TEAM7-like driver no longer contains solver-internal bookkeeping logic
- the main solve path is easier to map onto future `Dynamics` objects

Recommended tolerance policy for Phase 1:

- exact bitwise identity is not required
- relative behavior and final errors should remain in the same regime as the current prototype

## Phase 2 Preview: CPU SphinxSys-Style Solver Shape

Once Phase 1 is done, the next step should be to make the CPU solver look like `SphinxSys`, not merely be cleaner C++.

### Proposed CPU-side classes

- `InitializeMatrixFreeAPhiFields`
- `PrepareMatrixFreeAPhiSources`
- `MatrixFreeAPhiSingleIteration`
- `MatrixFreeAPhiIterationDiagnostics`
- `MatrixFreeAPhiConvergenceCheck`
- `MatrixFreeAPhiLoadStepDriver`
- `MatrixFreeAPhiFrequencySolver`

### Proposed CPU-side execution model

The desired shape is:

```cpp
for (size_t load_step = 0; load_step != load_step_count; ++load_step)
{
    while (!solver_state.converged_ && solver_state.outer_iteration_ < max_outer_iterations)
    {
        single_iteration.exec();
        diagnostics.exec();
        convergence_check.exec();
    }
}
```

This is still a steady-state frequency solver, but it starts to behave like a pseudo-time iterative process with observable state transitions.

## Phase 3 Preview: SYCL Migration Checklist

Only start this after Phase 2 is stable.

### Data migration targets

First migrate data that is naturally array-based:

- `ax`, `ay`, `az`, `phi`
- source arrays
- conductivity and reluctivity
- graph edge lists
- residual buffers
- workspace buffers

### Kernel migration order

Migrate in this order:

1. graph gradient kernels
2. harmonic-weighted gradient kernels
3. residual reductions
4. vector-potential update kernels
5. scalar-potential update kernels
6. gauge/reference postprocessing

### SYCL-specific design rule

Do not let the first `SYCL` version depend on `Eigen::Matrix<Complex, Dimensions, 1>` in device code. That boundary should already be cleaned before the migration begins.

## Risks And Watch Items

### Risk 1: Cleaning architecture while accidentally changing numerics

Mitigation:

- freeze benchmark baselines first
- refactor in small patches
- compare after each structural step

### Risk 2: Trying to finish CPU cleanup and SYCL port in one patch series

Mitigation:

- keep CPU architecture cleanup as an explicit milestone
- do not mix device-memory redesign into the first solver-state cleanup

### Risk 3: Keeping too much case-specific logic inside solver files

Mitigation:

- if a helper talks about TEAM7 geometry, region summaries, or thermal reporting, it probably does not belong in the long-term solver core

### Risk 4: Overfitting the architecture to the current single-box test

Mitigation:

- every new solver-side object should be checked against the question:
  can this still make sense after we split the setup into air, coil, and conductor bodies?

## Recommended Immediate Next Patch Series

If coding starts right away, the lowest-risk patch sequence is:

1. add iteration history/state structs in `electromagnetic_aphi_matrix_free_types.h`
2. add workspace ownership in `electromagnetic_aphi_matrix_free_solver.h/.hpp`
3. move scratch-buffer allocation out of the inner solve loop
4. split the main staggered solve into internal helper functions
5. trim the TEAM7-like test driver so it only orchestrates the case

## What Is Still Missing Before The Final Three-Body TEAM7 Goal

Relative to the final objective, the current code is still missing:

- proper three-body decomposition
- formal solver object design
- `SphinxSys`-style dynamics orchestration
- residual history output
- device-friendly data types
- device memory ownership
- `SYCL` kernels
- performance-oriented iteration strategy

So the project is not at the final integration stage yet. It is at the transition point between:

- successful numerical prototype

and

- formal framework integration

## Final Recommendation

Do not keep pushing the current test-driver-centered prototype toward production by tuning more environment variables.

Use the current prototype as:

- a numerical correctness reference
- a source of operator formulas
- a source of useful diagnostics

Then start a structured refactor toward:

- framework-style CPU solver first
- `SYCL` backend second
- three-body TEAM7 production case last

That path is much more likely to produce code that is:

- maintainable
- consistent with `SphinxSys`
- portable to GPU
- actually faster than the current prototype and potentially faster than the current Eigen reference in the intended large-scale regime

---

## 实现进度台账（外部维护）

落地代码变更与 SYCL 试点说明见同目录：

- `EM_APHI_MATRIX_FREE_SYCL_REFACTOR_PROGRESS.md`

