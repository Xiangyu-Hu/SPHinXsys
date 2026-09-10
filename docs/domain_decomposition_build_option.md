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
