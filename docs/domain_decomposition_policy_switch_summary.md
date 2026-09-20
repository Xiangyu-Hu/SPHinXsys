# Domain decomposition by execution policy: summary of the changes

Date: 2026-09-10. Branch `shuang_and_ming_1`, on top of commit `6993ba6c8` ("13.55").
Twenty files changed, one added, two removed. Nothing committed yet.

## Goal

Make the domain decomposition a property of the execution policy alone. A case is
written once, and the same source runs plain or decomposed depending on whether
`MainExecutionPolicy` is a `DecomposedExecution<>`. Compared with the original library
the only added preprocessor symbol and CMake option is `SPHINXSYS_DECOMPOSITION`.

## 1. One macro

`SPHINXSYS_MULTI_DEVICE` and `SPHINXSYS_MULTI_SUBDOMAIN_HOST` are replaced by
`SPHINXSYS_DECOMPOSITION`. `SPHINXSYS_USE_SYCL` picks the backend policy,
`SYCLDevicePolicy` or `ParallelPolicy`; the new option wraps it in
`DecomposedExecution<>` to form `MainExecutionPolicy`:

| `SPHINXSYS_USE_SYCL` | `SPHINXSYS_DECOMPOSITION` | `MainExecutionPolicy` |
|---|---|---|
| OFF | OFF | `ParallelPolicy` |
| OFF | ON  | `DecomposedExecution<ParallelPolicy>` |
| ON  | OFF | `SYCLDevicePolicy` |
| ON  | ON  | `DecomposedExecution<SYCLDevicePolicy>` |

Files: `CMakeLists.txt`, `execution_policy.h`, `sphinxsys.h`. (Committed as "13.55".)

## 2. No environment variables, no per-case helpers

Everything the bring-up steered through the environment is now an argument:

| Former | Now |
|---|---|
| `SPHINXSYS_SUBDOMAINS` | `--subdomains=N` on the command line, or `SPHSystem::setNumberOfSubdomains(N)`; both must precede the first body |
| `SPHINXSYS_RUNNER_THREADED` | `execution::subdomain_runner.initialize(N, Mode::Threaded)` |
| `SPHINXSYS_THREADS` | removed; all threads are used |
| `SPHINXSYS_VTP_INTERVAL`, `SPHINXSYS_DEBUG_STEPS` | were case-local; gone with the bring-up case |

More than one subdomain in a build without the option prints a warning and runs plain.
Files: `sph_system.{h,cpp}`.

## 3. The decomposition glue lives in the library

New `src/shared/domain_decomposition/body_decomposition.h`: `BodyDecomposition<Policy>`,
the decomposition of one body over the subdomains. The primary template is the
non-decomposed no-op; the partial specialization on `DecomposedExecution<>` owns the
slab cut planes (balanced on the initial particle positions), the halo and migration
exchange, and an optional owner tag for the vtp output. The exchange set is fixed at
construction: the evolving variables, the variables registered for output, and
`addExchangeVariable<T>(name)` additions made before `scatterFromHost()`.

The loop dynamics `UpdateHaloCK`, `MigrateParticlesCK`, `SyncHaloStateCK` and
`RebalanceSubdomainsCK` take the decomposition object and are no-ops with the plain
one. The method container offers `addDecomposition(body)`; the dynamics are added with
the existing `addGeneralDynamics<UpdateHaloCK>(decomposition)`.

`SubdomainExchange` gained `addExchangeVariable<T>()` for a single late addition, and
its implementation file includes `cell_linked_list.h` before `particle_iterators.h`,
which used to rely on include order for `ConcurrentVec`.

## 4. The I/O handles a decomposed body on its own

The particle-level interface, formerly the halo refresher, is now
`SubdomainExchangeInterface` with `refreshHalo()`, `gatherToHost()`,
`finishHostAccess()` and `subdomainOf(position)`. `BodyDecomposition` implements it and
installs itself on the body's particles; all four are no-ops on a plain body.

- The vtp recorder, the restart output and the reload output bracket their host access
  with gather and finish. A case places no gather any more.
- `ObservedQuantityRecording` under a decomposed policy: the observer body is
  replicated, so each subdomain interpolates every observer point from its own owned
  and halo particles; the recorder takes each value from the subdomain owning the
  observation point, through the subdomain copy helper, so the device path shares the
  code.
- Restart works: the output gathers, the read precedes `scatterFromHost()`.

Files: `base_particles.h`, `io_base_ck.{h,hpp}`, `io_observation_ck.h`.

## 5. Two defects fixed

Both surface only when one subdomain grows its neighbor list after the others have
already built theirs in the same step. In the dam break that happens at about three
seconds of physical time, which the two-second bring-up case never reached; the run
then went to NaN and looked like a hang.

- **Replica reallocation discarded contents.** `HostOnlyDiscreteVariable` and
  `DeviceOnlyDiscreteVariable` now keep their contents on `reallocateData()`.
- **The relation update did not register its kernel with the relation.** The growth
  invalidated the interaction kernels of every subdomain but not the update kernels, so
  the other subdomains kept pointers into the freed replica and counted garbage
  neighbors. `UpdateRelation` registers its kernel now, for inner and contact relations.

Files: `sphinxsys_variable.h`, `sphinxsys_variable_sycl.hpp`, `update_body_relation.hpp`.

## 6. Cases and tests

- `test_2d_dambreak_sycl` runs decomposed by policy. It gained the decomposition object,
  two exchange variables (`Pressure`, `LinearCorrectionMatrix`), the four loop dynamics
  and one `scatterFromHost()` after the restart read. No conditional compilation.
- In a decomposed build ctest additionally runs it with `--subdomains=2` and a restart
  from those files (`tests/tests_sycl/2d_examples/test_2d_dambreak_sycl/CMakeLists.txt`).
- `test_2d_dambreak_decomp`, the bring-up case, is removed. Its job is done by the
  original case.

## 7. Verification

Two build directories of this checkout, `build` (option OFF) and `cmake-build-decomp`
(option ON), both without SYCL, all threads:

| Run | Result |
|---|---|
| OFF | both DTW regression tests pass |
| ON, `--subdomains=1` | both DTW regression tests pass |
| ON, `--subdomains=2` | both DTW regression tests pass, 20 s wall time (12 s with one subdomain) |
| restart at step 4000, OFF and `--subdomains=2` | both complete |
| energy and observer records, ON vs. OFF | identical through 1.8 s of physical time |

Beyond about two seconds two plain runs of this case differ from each other, on this
branch and without any decomposition: the one-sided inner relation appends neighbors
with atomic counters, so the summation order and then the trajectory vary run to run.
Bitwise comparisons are therefore limited to the first records, and the pressure
regression test occasionally trips its tolerance for the same reason. Not changed.

Not verified: SYCL, no compiler on this machine. The device-side changes use the same
helpers as the existing SYCL code.

## 8. Open

- Threaded host runner: still restricted by the reallocation issue listed in the design
  document (§10); step 4 of its bring-up order.
- Multi-GPU: step 5 of the bring-up order.
- The two `SyncHaloStateCK` lines of the case (volume, correction matrix) still encode
  case knowledge; a "produced variable" declaration on the writing dynamics would let
  them go.
- Run-to-run reproducibility of the inner relation, if bitwise gates over a whole run
  are wanted.

## Related documents

`docs/domain_decomposition.md` (design, §8 shows the case-level API, §11–12 the test
state), `docs/domain_decomposition_build_option.md` (the step-by-step log with the
file lists).
