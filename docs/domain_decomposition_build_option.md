# One build option for domain decomposition

Date: 2026-09-10. Branch `shuang_and_ming_1`.

## Motivation

The decomposed execution policy is a template, `DecomposedExecution<PolicyType>`,
and the backend it wraps is already decided by `SPHINXSYS_USE_SYCL`. Two separate
build options, `SPHINXSYS_MULTI_DEVICE` (GPU) and `SPHINXSYS_MULTI_SUBDOMAIN_HOST`
(CPU), therefore encoded a choice the compiler could make on its own. They also had
to be kept mutually exclusive by hand, and every `#if` that wanted "any decomposed
run" had to test both.

This change replaces them with a single option, `SPHINXSYS_DECOMPOSITION`. It only
says whether the run is decomposed. Which decomposition is used follows from the
backend.

## Selection rule

`src/shared/particle_dynamics/execution/execution_policy.h`:

```cpp
#if SPHINXSYS_USE_SYCL
using BackendExecutionPolicy = SYCLDevicePolicy;
#else
using BackendExecutionPolicy = ParallelPolicy;
#endif

#if SPHINXSYS_DECOMPOSITION
using MainExecutionPolicy = DecomposedExecution<BackendExecutionPolicy>;
#else
using MainExecutionPolicy = BackendExecutionPolicy;
#endif
inline constexpr auto par_ck = MainExecutionPolicy{};
```

| `SPHINXSYS_USE_SYCL` | `SPHINXSYS_DECOMPOSITION` | `MainExecutionPolicy`                          |
|----------------------|---------------------------|------------------------------------------------|
| OFF                  | OFF                       | `ParallelPolicy`                               |
| OFF                  | ON                        | `DecomposedExecution<ParallelPolicy>` (`MultiHostPolicy`)   |
| ON                   | OFF                       | `SYCLDevicePolicy`                             |
| ON                   | ON                        | `DecomposedExecution<SYCLDevicePolicy>` (`MultiDevicePolicy`) |

The named aliases `MultiDevicePolicy`, `MultiHostPolicy` and
`SequencedMultiHostPolicy`, and the objects `multi_device`, `multi_host` and
`seq_multi_host`, are unchanged. `DecomposedExecution<SequencedPolicy>` has no build
option: it is used explicitly where a sequenced decomposed loop is wanted, and the
sequential visiting of subdomains on one thread is a run time choice of
`SubdomainRunner::Mode`, not a compile time one.

## Files touched

- `CMakeLists.txt`: one `option()` and one `target_compile_definitions()` replace the
  former two of each.
- `src/shared/particle_dynamics/execution/execution_policy.h`: the nested four-way
  `#if` becomes backend selection followed by the decomposition wrap.
- `src/shared/include/sphinxsys.h`: `domain_decomposition_dynamics.h` and
  `subdomain_exchange.hpp` are included under `#if SPHINXSYS_DECOMPOSITION`.
- `tests/tests_sycl/2d_examples/test_2d_dambreak_decomp/dambreak_decomp.cpp`: all
  guards renamed; the OFF build still follows the original single-domain path.
- `docs/domain_decomposition.md`: build instructions and test plan updated.

`docs/domain_decomposition_commit summary.md` and
`docs/domain_decomposition_review_after_merge.md` keep the old names, as a record of
the state they describe.

## Building

```
cmake -DSPHINXSYS_DECOMPOSITION=ON ...                          # host subdomains
cmake -DSPHINXSYS_USE_SYCL=ON -DSPHINXSYS_DECOMPOSITION=ON ...  # one subdomain per GPU
```

An existing build directory configured with the old options should be reconfigured
with the new one. CMake ignores unknown `-D` values silently, so a stale
`-DSPHINXSYS_MULTI_SUBDOMAIN_HOST=ON` now produces a non-decomposed build without
any warning. The stale cache entries can be dropped with
`cmake -U SPHINXSYS_MULTI_DEVICE -U SPHINXSYS_MULTI_SUBDOMAIN_HOST`.

## Verification

Both build directories of this checkout (`build`, OFF; `cmake-build-decomp`, ON) were
reconfigured with the old cache entries removed and the target
`test_2d_dambreak_decomp` rebuilt without errors or warnings. The case was then run
with `SPHINXSYS_THREADS=1` and `WaterBody_TotalMechanicalEnergy.dat` compared with
`diff`:

| Comparison                                                | Result    |
|-----------------------------------------------------------|-----------|
| new OFF build vs. OFF reference recorded before the change | identical |
| ON, 1 subdomain vs. OFF                                    | identical |
| ON, 2 subdomains vs. OFF                                   | identical |
| ON, 2 subdomains vs. 2-subdomain reference before the change | identical |

The SYCL path was not compiled: this machine has no IntelLLVM compiler. Only the alias
selection changes on that path.

## Step 1: the decomposition glue moved into the library

Date: 2026-09-10, same branch. Goal: a case drives the decomposition without any
`#if`, so that the original `test_2d_dambreak_sycl` can later take it over by policy
alone. The decomposition case `test_2d_dambreak_decomp` is kept as the bitwise gate.

New: `src/shared/domain_decomposition/body_decomposition.h`. `BodyDecomposition<Policy>`
is the decomposition of one body over the subdomains. The primary template is the
non-decomposed no-op; the partial specialization on `DecomposedExecution<>` owns the
`SlabDecomposition` (balanced on the initial particle positions), the
`SubdomainExchange`, and the optional `SubdomainID` output tag, which `gatherToHost()`
fills on the host. The exchange set is fixed at construction: evolving variables, the
variables registered for output, and `addExchangeVariable<T>(name)` additions made
before `scatterFromHost()`.

Changed:

- `domain_decomposition_dynamics.h`: `UpdateHaloCK`, `MigrateParticlesCK`,
  `SyncHaloStateCK`, `RebalanceSubdomainsCK` take a `BodyDecomposition<Policy> &` and
  are no-ops with it under a plain policy. They are added through the existing
  `main_methods.addGeneralDynamics<UpdateHaloCK>(decomposition)`.
- `particle_method_container.h`: `addDecomposition(body, split_axis)`.
- `subdomain_exchange.{h,hpp}`: `addExchangeVariable<T>(variable)` creates the send
  buffers of one more variable after construction. The `.hpp` now includes
  `cell_linked_list.h` before `particle_iterators.h`, which uses `ConcurrentVec` without
  declaring it.
- `sph_system.{h,cpp}`: `setNumberOfSubdomains(N)` and the command line option
  `--subdomains=N` initialize the subdomain runner; both must precede the first body.
  A value larger than 1 in a build without `SPHINXSYS_DECOMPOSITION` prints a warning
  and stays unused. (A first version read environment variables instead; they are gone,
  see below.)
- `sphinxsys.h`: the decomposition headers are included unconditionally; the macro now
  only selects `MainExecutionPolicy`.
- `dambreak_decomp.cpp`: every `#if SPHINXSYS_DECOMPOSITION` is gone. The runner
  initialization, the balancing loop, the exchange construction and the output tag left
  the case; what remains is the placement of the calls in the loop, which is the part
  that still needs a case's knowledge.

Verification, energy record and every vtp frame compared with `diff`:

| Comparison                                        | Result    |
|---------------------------------------------------|-----------|
| OFF vs. OFF reference before this step             | identical (energy and all vtp) |
| ON, 1 subdomain vs. OFF                            | identical |
| ON, 2 subdomains vs. OFF                           | identical |
| ON, 2 subdomains vs. 2-subdomain reference before this step | identical (energy and all vtp, including `SubdomainID`) |

Not covered by this step: the gather around output still sits in the case (step 2, an
output hook), the observer path, restart, and the automatic refresh of variables such as
the volume and the correction matrix, which the case still names in `SyncHaloStateCK`.

The case's `SPHINXSYS_THREADS` control and its `environmentInt` helper were removed
afterwards: the case uses all available threads like every other case. The gates were
rerun that way and stayed bitwise identical, in two consecutive OFF runs as well as ON
with 1 and 2 subdomains against OFF, so single-threaded runs are not needed for them.

## Only one macro

Compared with the original library, `SPHINXSYS_DECOMPOSITION` is the only added
preprocessor symbol and the only added CMake option. Everything that was steered by
environment variables in the bring-up (`SPHINXSYS_THREADS`, `SPHINXSYS_SUBDOMAINS`,
`SPHINXSYS_RUNNER_THREADED`, `SPHINXSYS_VTP_INTERVAL`, `SPHINXSYS_DEBUG_STEPS`) is now a
plain argument:

| Former environment variable   | Now                                                        |
|-------------------------------|------------------------------------------------------------|
| `SPHINXSYS_THREADS`           | removed; the case uses all threads, runs are deterministic |
| `SPHINXSYS_SUBDOMAINS`        | `--subdomains=N`, or `SPHSystem::setNumberOfSubdomains(N)` |
| `SPHINXSYS_RUNNER_THREADED`   | `execution::subdomain_runner.initialize(N, Mode::Threaded)` |
| `SPHINXSYS_VTP_INTERVAL`      | the constant `state_recording_interval` of the case        |
| `SPHINXSYS_DEBUG_STEPS`       | the constructor argument of the case's `TimeStepTracer`    |

The gates are run as `./test_2d_dambreak_decomp` and `./test_2d_dambreak_decomp --subdomains=2`.
`SYCL_EXT_ONEAPI_PEER_ACCESS`, tested in `device_environment_sycl.cpp`, is a feature
macro of the SYCL implementation, not one of ours.

## Step 2: the original case runs decomposed, the bring-up case is retired

Date: 2026-09-10, same branch.

Library:

- `SubdomainExchangeInterface` (formerly the halo refresher) in `base_particles.h` now
  also offers `gatherToHost()`, `finishHostAccess()` and `subdomainOf(position)`;
  `BodyDecomposition` implements it and installs itself on the body's particles. The
  vtp recorder, the restart output and the reload output bracket their host access with
  the two calls themselves, so a case places no gather any more.
- `ObservedQuantityRecording` under a decomposed policy: the observer body is
  replicated, every subdomain interpolates every observer point from its own owned and
  halo particles, and the recorder takes each value from the subdomain owning the
  observation point (through the subdomain copy helper, so the multi-GPU path uses the
  same code).
- Two fixes for a subdomain growing its neighbor list after the others built theirs in
  the same step: `HostOnlyDiscreteVariable` and `DeviceOnlyDiscreteVariable` keep their
  contents on `reallocateData()`, and `UpdateRelation` registers its kernel with the
  relation so that the growth invalidates the update kernels of the other subdomains
  too. Without them the run went to NaN at the second growth, around three seconds of
  physical time, which the short bring-up case never reached.

Cases:

- `test_2d_dambreak_sycl` gained the decomposition object, two exchange variables, the
  four loop dynamics and one `scatterFromHost()` after the restart read. In a decomposed
  build ctest additionally runs it with `--subdomains=2` and a restart from those files.
- `test_2d_dambreak_decomp` is removed; the original case is the gate now.

Verification (this machine, `SPHINXSYS_DECOMPOSITION` OFF and ON builds):

| Run                                        | Result                                   |
|--------------------------------------------|------------------------------------------|
| OFF                                        | both regression tests pass               |
| ON, `--subdomains=1`                       | both regression tests pass               |
| ON, `--subdomains=2`                       | both regression tests pass, 20 s wall    |
| restart at step 4000, OFF and `--subdomains=2` | both complete                        |
| energy record, ON-1 and ON-2 vs. OFF       | identical through 1.8 s physical time; from 2.5 s the plain run itself differs between two of its own runs (atomic neighbor append, see the design document §11) |

Not verified: SYCL. The device replica fix and the observer copy on the device path
compile against the same helpers the existing SYCL code uses, but no SYCL compiler is
available here.
