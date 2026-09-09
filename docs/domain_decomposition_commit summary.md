# Domain decomposition bring-up: the two-subdomain host path on the 2D dambreak

Date: 2026-09-09. Branch: `hackathon/domain-decomp`, on top of 7bf4b8882 "Domain decomposition second attempt".

## 1. What this change does

The domain decomposition draft described in `docs/domain_decomposition.md` had never been compiled against the real dependencies. This change brings it up on the **host path** (`SPHINXSYS_MULTI_SUBDOMAIN_HOST=ON`, `SubdomainRunner::Mode::Sequential`) and validates it on a new case, `tests/tests_sycl/2d_examples/test_2d_dambreak_decomp/`.

State reached:

| Configuration | Compared with the OFF build, single TBB thread |
|---|---|
| ON, 1 subdomain | total mechanical energy identical bit for bit (23 rows); vtp identical per particle (20 frames) |
| ON, 2 subdomains, 1 thread | energy identical bit for bit; per-particle position difference in vtp ≤ 1e-9 (Float32 print precision) |
| ON, 2 subdomains, 8 threads | energy identical bit for bit |
| ON, 2 subdomains, AddressSanitizer | no out-of-bounds access, no use-after-free |
| Consistency | `checkConsistency()` empty at every advection step; particle total stays 3200; owned counts 1600/1600 → 1347/1853 (migration is happening) |
| Original case `test_2d_dambreak_sycl` | still builds and runs in the OFF build; DTW regression test passes |

Size: 18 existing files (+382 / −108), one new case directory (`CMakeLists.txt` and `dambreak_decomp.cpp`).

## 2. Changes and why

### 2.1 The draft did not compile (not even the OFF build), 6 items

| File | Problem | Fix |
|---|---|---|
| `src/shared/common/sphinxsys_variable.h` | both `DiscreteVariable` constructors initialized `device_only_variable_`, now a `std::array`, with `nullptr` | removed the initializer; the member is value-initialized at its declaration |
| `src/shared/domain_decomposition/domain_decomposition.h` | included a non-existent `base_data_package.h` | `base_data_type_package.h` |
| `src/shared/domain_decomposition/domain_decomposition.cpp` | `BoundingBox` has no `first_/second_` | `lower_/upper_` |
| `src/shared/include/sphinxsys.h` | the decomposition headers were included first, before `ConcurrentVec` used by `particle_iterators.h` is declared | moved after all other headers |
| `src/shared/shared_ck/particle_dynamics/particle_iterators_ck.h` | no `particle_for` / `particle_reduce` overload for `LoopRangeCK<ParallelMultiHostPolicy, …>` (template arguments do not undergo derived-to-base conversion) | loop bodies factored into 4 generic functions; overloads added for `SequencedMultiHostPolicy` and `ParallelMultiHostPolicy` |
| `src/shared/particle_dynamics/particle_iterators.h`, `src/shared/common/algorithm_primitive.h`, `src/shared/meshes/mesh_iterators.h` | immediate exit at run time: `particle_for(ParallelMultiHostPolicy{}, IndexRange, f)` is an exact match for the generic catch-all template, so the `ParallelPolicy&` overload, which needs a derived-to-base conversion, is never chosen | forwarding overloads on `MultiHostExecution<PolicyType>` for `particle_for`, `particle_reduce`, `generic_for`, `exclusive_scan`, `mesh_for`, `package_for`; each `static_cast`s to the base policy and calls again |

### 2.2 Logic bugs in the draft, 5 items

| Id | File | Problem | Fix |
|---|---|---|---|
| B1 | `subdomain_exchange.hpp` constructor | subdomain counts were written with `setValue`, but replicas are created lazily; before a replica exists `SingleVariable::Data()` addresses the host value, so the last subdomain overwrote the others | the constructor pre-touches the replicas of `TotalRealParticles` and `TotalLocalParticles` for every subdomain (`SubdomainScope` + `DelegatedData()`) |
| B2 | `subdomain_exchange.{h,hpp}` `gatherToHost` | the final `setValue(total)` runs on the host thread, i.e. it writes the replica of subdomain 0, whose next loop then covers the global count | writes `setValue(0, total)` explicitly; new `finishHostAccess()` restores the owned count of subdomain 0; the case does gather → write vtp → finish |
| B3 | `shared_ck/.../update_body_relation.hpp` | the inner relation is one-sided: a pair of an owned `i` and a halo `j` writes `particle_offset_[j]`, but the scan only covered `n_owned + 1` | the four loops and the `exclusive_scan` cover `TotalLocalParticles()`; equal to the owned count outside a decomposed run |
| B4 | `domain_decomposition_dynamics.h`, `subdomain_exchange.h`, `docs/domain_decomposition.md` §8 | the documented order was cell linked list → UpdateHalo, but the list is built over `TotalLocalParticles`, which UpdateHalo publishes | order is migrate → sort → **UpdateHalo → cell linked list** → relation |
| B7 | `subdomain_exchange.hpp` `CompactVariables` | migration compacted the whole owned range (≈3200 entries) into a scratch buffer sized for the send capacity (1024): heap overflow, found by ASan | eliminated by the migration rewrite of 2.4; the scratch buffer is gone |

### 2.3 Halo refresh inside an acoustic half step: split fan-out plus a hook

`src/shared/shared_ck/particle_dynamics/interaction_algorithms_ck.{h,hpp}`

- `InteractionDynamicsCK<…, OneLevel, …>::exec()` fans out three times instead of once: Initialize → hooks → Interact (with its pre- and post-processes) → Update. A fan-out ends with a barrier over the subdomains, and the hooks run on the host thread when every subdomain has finished Initialize.
- `InteractionDynamicsCK<OneLevel>` gets `addPostInitialization(BaseDynamics<void>&)`.
- Under a non-decomposed policy the three fan-outs collapse into plain calls, equivalent to the previous `runAllSteps()`.

Why this is necessary: the Interact step of the first half reads the neighbors' `Pressure`, which the Initialize step of the same `exec()` has just written. A refresh can only sit between two fan-outs; with one fan-out around the whole `exec()` there was no such point. Calling `SyncHaloStateCK` from inside a fan-out degenerates to a single-subdomain call through `insideFanOut()`: no barrier, and no copy for the other subdomains.

What is actually transferred per acoustic step: `{Pressure}` in the first half, `{Velocity}` in the second. In the advection step, `{VolumetricMeasure}` after `AdvectionStepSetup` and `{LinearCorrectionMatrix}` after `LinearCorrectionMatrix`; `UpdateHaloCK` after migrate/sort refreshes everything. The sets come from checking, kernel by kernel, which quantities are read at the neighbor index `j`.

### 2.4 Migration: pack first, delete by filling from the tail, arrivals are the generation

`src/shared/domain_decomposition/subdomain_exchange.{h,hpp}`

- Removed the `Scratch` buffer, `CompactVariables` and `dv_keep_index_`; added `dv_hole_index_`, `dv_donor_index_`, `fill_count_` and `FillHolesFromTail`.
- The first fan-out uses the existing flag + `exclusive_scan` to get the number of departing particles, hence the new owned count `n_new`; the departing slots below `n_new` (holes, ascending); the staying particles at or above `n_new` (donors, ascending). The two counts are necessarily equal; a mismatch is an `exit(1)`.
- The second fan-out packs the departing particles into the send buffers first (their slots are about to be overwritten), then does `data[hole[k]] = data[donor[k]]`. O(departing) copies, no atomics, deterministic.
- The third fan-out is unchanged: pull from the neighbors' send buffers into the slots after `n_new`, `owned += count`. This is the generation side: the state of a new particle comes from the neighbor's buffer, and `OriginalID` travels with it.

Why `RemoveRealParticle` is not called directly: it claims the tail with `fetch_sub`; its own comment restricts it to a sequenced policy, and an interleaving that loses a particle can be constructed under parallel execution. Its `LifeStatus` bit lives in `ParticleGroups`, an evolving variable, and would be packed and shipped to the neighbor.

### 2.5 A latent bug in the CK library: `Force` was not an evolving variable

`src/shared/shared_ck/particle_dynamics/fluid_dynamics/acoustic_step_1st_half.hpp`

The second half step assigns `Force` with `=`; the next first half step accumulates onto it with `+=`; an advection step with a particle sort sits in between. The sort permutes the evolving variables only, and `Force` was not among them (the classic `fluid_integration.hpp:39` does add it). The first half step therefore accumulated onto another particle's dissipative force. In a single domain the error is deterministic; with two subdomains each sorts its own particles and the error differs, so the runs diverged over the whole field, at the 1e-4 level, right after the first sort. With `Force` evolving, two subdomains match the single domain bit for bit; the original case still passes its regression test.

### 2.6 Threaded runner: partial fix, not finished

- `src/shared/particle_dynamics/execution/subdomain_scope.h` adds `replicaCreationMutex()`; the lazy replica creation in `sphinxsys_variable.h` and `src_sycl/.../sphinxsys_variable_sycl.hpp` is now double-checked under that lock. This removes the use-after-free in `UniquePtrsKeeper` when two subdomain threads create replicas of the same variable at once.
- Still crashes: `DiscreteVariable::reallocateData` grows the shared capacity and every subdomain's replica from whichever subdomain thread triggers it, while the other threads still hold the old pointers in their computing kernels. Per-replica capacities are needed; recorded in `docs/domain_decomposition.md` §10. **This change guarantees the sequential runner only.**

### 2.7 New case `tests/tests_sycl/2d_examples/test_2d_dambreak_decomp/`

Copied from `test_2d_dambreak_sycl` without the observer, the DTW regression tests and RestartIO; `end_time = 2.0`; the energy is recorded every 20 advection steps. Everything related to the decomposition is inside `#if SPHINXSYS_MULTI_SUBDOMAIN_HOST`, so the OFF build follows the same physics path as the original case.

Environment variables:

| Variable | Meaning | Default |
|---|---|---|
| `SPHINXSYS_THREADS` | TBB threads; 1 makes a run bit reproducible | all cores |
| `SPHINXSYS_SUBDOMAINS` | number of subdomains (ON build only) | 1 |
| `SPHINXSYS_RUNNER_THREADED` | 1 = one host thread per subdomain (known to crash) | 0 |
| `SPHINXSYS_VTP_INTERVAL` | state recording interval in seconds | 0.1 |
| `SPHINXSYS_DEBUG_STEPS` | print the time step of the first N acoustic steps | 0 |

Additional work in the ON build: the initial cut plane is moved from the middle of the tank to the middle of the water column with `rebalance` (otherwise subdomain 1 starts empty); the exchange set is `EvolvingVariables()` plus `Pressure` and `LinearCorrectionMatrix` (needed on the halo) and `Density` (needed for the vtp output); `checkConsistency()` runs every advection step and exits on a non-empty report; owned counts and the halo load factor are printed every 100 steps; the vtp carries an extra `SubdomainID` field, filled on the host after `gatherToHost()`, which lays the owned particles out subdomain by subdomain.

The comparisons below were done with small throw-away Python scripts (energy compared on matching time stamps and by interpolation; vtp frames compared per particle, aligned on `OriginalID`). They are not part of the repository.

### 2.8 Documentation

`docs/domain_decomposition.md`: §7 describes migration as tail filling; §8 puts UpdateHalo before the cell linked list; §9 marks the one-sided relation as fixed and adds the notes "state that survives an advection step must be evolving" and "halo refresh points"; §10 drops the one-sided relation entry and adds the threaded-runner limitation.

## 3. How it was tested

Two build directories from the same sources: `build` (OFF, pre-existing) and `cmake-build-decomp` (ON: `-DSPHINXSYS_MULTI_SUBDOMAIN_HOST=ON`, otherwise the same as `build`). A third, `cmake-build-asan`, is ON with `-fsanitize=address`. All three match `cmake-build-*/` in `.gitignore`.

Every run is started in the case's `bin/` directory, since the case uses the relative `output/` path.

**Level 0, calibrating the noise.** OFF with `SPHINXSYS_THREADS=1`, run twice: the energy is identical bit for bit, so a deterministic reference exists. OFF with 8 threads against 1 thread: relative difference about 1e-9, the size of floating-point association noise. Multi-threaded TBB runs of the same binary are not bit reproducible from one run to the next, which is why every bitwise statement below is made with a single thread.

**Level 1, ON with 1 subdomain against OFF.** Identical bit for bit. Proves that the replicas, scatter/gather, the three fan-outs and the empty hooks are transparent to the physics in the degenerate case.

**Level 2, ON with 2 subdomains against OFF.** Energy identical bit for bit, vtp per particle ≤ 1e-9; on screen the owned counts always sum to 3200 and shift as the water advances; no `Decomposition inconsistent` line. Proves that the four refresh points transfer the right quantities at the right time, that migration neither loses nor misplaces particles, and that building the relation over the local range is correct.

**Locating a difference (used to find the `Force` bug).** `SPHINXSYS_DEBUG_STEPS=30` prints the first 30 acoustic time steps in OFF and ON; `paste` the two lists, and the first differing step is where the states first differ (the time step is a global maximum, so a change in it means some particle already differed in the previous step). `SPHINXSYS_VTP_INTERVAL=0.001` writes a vtp at every advection step; comparing the frames particle by particle (aligned on `OriginalID`) and looking at where the differing particles sit tells the two cases apart: concentrated near the cut plane means a halo problem, everywhere at once means a global quantity or the sort.

**ASan.** Two subdomains run to completion with no report. It had earlier pointed directly at the heap overflow of B7 and at the use-after-free of the threaded mode.

**Regression.** The original `test_2d_dambreak_sycl` re-run in the OFF build passes its DTW regression test (the difference introduced by the `Force` fix is within tolerance).

## 4. Reproducing

```bash
cd ~/vcpkg/sphinxsys
cmake -S . -B cmake-build-decomp -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_TOOLCHAIN_FILE=$HOME/vcpkg/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=arm64-osx \
  -DSPHINXSYS_MULTI_SUBDOMAIN_HOST=ON -DSPHINXSYS_BUILD_3D_EXAMPLES=OFF -DSPHINXSYS_BUILD_PYTHON_INTERFACE=OFF
cmake --build build --target test_2d_dambreak_decomp
cmake --build cmake-build-decomp --target test_2d_dambreak_decomp

# OFF reference
cd build/tests/tests_sycl/2d_examples/test_2d_dambreak_decomp/bin
SPHINXSYS_THREADS=1 ./test_2d_dambreak_decomp && cp -r output baseline

# ON, two subdomains
cd ../../../../../cmake-build-decomp/tests/tests_sycl/2d_examples/test_2d_dambreak_decomp/bin
SPHINXSYS_THREADS=1 SPHINXSYS_SUBDOMAINS=2 ./test_2d_dambreak_decomp
B=../../../../../../build/tests/tests_sycl/2d_examples/test_2d_dambreak_decomp/bin/baseline
diff $B/WaterBody_TotalMechanicalEnergy.dat output/WaterBody_TotalMechanicalEnergy.dat
```

Expected: `diff` prints nothing (the energy files are identical byte for byte); the vtp frames, compared per particle after aligning on `OriginalID`, differ by at most 1e-9 in position. On screen, the owned counts sum to 3200 at every report and no `Decomposition inconsistent` line appears.

## 5. Explicitly not done

- The threaded runner (see 2.6).
- The SYCL multi-device path: the device side of `copyBetweenSubdomains` and `device_environment_sycl` were not compiled on this machine (no SYCL toolchain); the locking change in `sphinxsys_variable_sycl.hpp` is likewise uncompiled.
- Load balancing with `RebalanceSubdomainsCK` (`rebalance` is only called once, at setup), per-subdomain cell meshes, deeper halos, restart.
- `setValue` on a `SingleVariable` from the host thread only reaches the replica of subdomain 0 (`PhysicalTime`, for instance). The constant gravity of the dambreak does not read the time, so this case is unaffected; time-dependent cases will need it addressed.
