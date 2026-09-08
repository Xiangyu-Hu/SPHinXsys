# Domain decomposition in SPHinXsys: multi-GPU, and a CPU path to debug it

**Status: draft.** The parts that can be built without the project's dependencies are
built and tested here (see §11); the rest — anything touching particles, SYCL or Eigen —
is a reviewed design sketch with the plumbing written out, not tested code. This machine
has neither a SYCL toolchain nor Simbody/TBB/Eigen. §10 lists what is stubbed and what
is likely to break first.

## 1. Two backends, one implementation

There are two decomposed execution policies, and they share essentially all of their
code:

| | `ParallelMultiDevicePolicy` | `ParallelMultiHostPolicy` |
| --- | --- | --- |
| a subdomain is | one GPU | one host thread, or one pass of a loop |
| replicas live in | SYCL USM, shared context | ordinary host memory |
| a copy between them is | `queue.memcpy` | `std::copy` |
| purpose | the actual speedup | **debugging vehicle for the left column** |

Everything else — the decomposition geometry, the halo plan, the pack/pull protocol,
migration, load balancing, the fan-out points in every algorithm's `exec()` — is one
piece of code, instantiated on both. The only backend-specific call in the exchange is
`execution::copyBetweenSubdomains`.

That is deliberate. A decomposition bug reproduced on the CPU path is *the same bug*,
reachable in a debugger, on a laptop, with no device toolchain, and with deterministic
ordering. Chasing it on eight GPUs is not a good use of anyone's time.

### Why queues per device rather than MPI (for the GPU path)

For a single node, one process driving N devices avoids everything that makes MPI
expensive to introduce here: no separate ranks, no duplicated `SPHSystem`, no
serialization layer, no launcher change. The decisive point is that all devices can be
placed in **one `sycl::context`**, so a USM pointer allocated for device A is a legal
argument on device B's queue: the halo exchange is a plain `memcpy`, using the peer link
when there is one and staging through the host when there is not, with no second code
path. The cost is that it does not scale past one node.

## 2. The central idea: thread-local subdomain scope

Every subdomain is driven by a host thread that announces itself with an
`execution::SubdomainScope`, setting a thread-local id. Every per-subdomain resource is
then resolved implicitly through `currentSubdomainID()`:

| Resource | Where it is resolved |
| --- | --- |
| `sycl::queue` | `ExecutionInstance::getQueue()` → `DeviceEnvironment::getCurrentQueue()` |
| device replica of a `DiscreteVariable` | `DiscreteVariable::DelegatedOnDevice()` |
| **host replica** of a `DiscreteVariable` | `DiscreteVariable::DelegatedOnHostSubdomain()` |
| replica of a `SingleVariable` | `SingleVariable::DelegatedOn{Device,HostSubdomain}()` |
| computing kernel | `Implementation<...>::getComputingKernel()` |
| freshness flag | `Implementation<Base>::isUpdated()` |
| loop range | `LoopRangeCK<ParallelMultiDevicePolicy, ...>` |

The consequence is that **the physics code does not change at all**. A computing
kernel's constructor calls `DelegatedData()` on the variables it reads; run that
constructor inside `SubdomainScope(2)` and every pointer in it addresses subdomain 2's
replica. `AcousticStep1stHalf`, `LinearCorrectionMatrix` and the rest are untouched.

`HostOnlyDiscreteVariable` is deliberately shaped exactly like
`DeviceOnlyDiscreteVariable`, including the rule that the replicas hold disjoint
particle sets and `data_` is only I/O staging. Symmetry is the point: if the host path
diverged in its aliasing rules, a bug found there would not be the bug you have.

When `numberOfSubdomains() == 1`, `currentSubdomainID()` is always 0, every array has
one live entry, and behavior is identical to the current single-GPU / single-CPU path.

### Where the fan-out lives

Parallelism over subdomains is expressed **at the algorithm level**, not inside
`particle_for`. `StateDynamics::exec()` becomes:

```cpp
execution::fanOutOverSubdomains(ExecutionPolicy{}, [&]() {
    UpdateKernel *update_kernel = kernel_implementation_.getComputingKernel();
    particle_for(LoopRangeCK<ExecutionPolicy, RangeIdentifier>(*this->identifier_),
                 [=](size_t i) { update_kernel->update(i, dt); });
});
```

`fanOutOverSubdomains` is a no-op wrapper for every non-decomposed policy, so one
`exec()` body serves all of them. This placement is what makes the kernel lookup and the
loop-range construction happen *inside* the subdomain thread. Putting the fan-out inside
`particle_for` would not work: the kernel pointer is captured before the loop is
entered, so every subdomain would get subdomain 0's kernel.

Reductions take the symmetric route: `reduceOverSubdomains<Operation>` runs the body per
subdomain and combines the partials in subdomain order. A global reduction therefore
differs from the non-decomposed result by floating-point association — the same caveat
as an MPI reduction. **Measure that difference on the CPU path first**; it tells you how
much of a multi-GPU regression failure is expected before you go looking for a bug.

Dynamics compose (an interaction runs its pre- and post-processes, themselves dynamics
with their own `exec()`), so a nested fan-out must degenerate to a plain call.
`insideFanOut()` handles that; it is covered by a unit test, because getting it wrong is
either a deadlock or a silent N× duplication of work.

## 3. The CPU path's two modes

`SubdomainRunner` runs the fan-out in one of two modes, and running both is the intended
workflow:

- **Sequential** (default) — subdomains visited one after another on the calling thread.
  Deterministic, single-threaded, steppable in a debugger, bit-reproducible. Barriers are
  trivially satisfied, so a data race cannot mask a logic error.
- **Threaded** — one persistent host thread per subdomain, with the same barrier
  structure as the multi-GPU path. Meant to be run under ThreadSanitizer.

The diagnostic value is in the difference:

> A difference between **1 and N subdomains** is a decomposition bug.
> A difference between **sequential and threaded** is a synchronization bug.

Separating those two questions is the single biggest reason this path is worth having.

## 4. Particle data layout

Each subdomain holds a **local particle set**, not a copy of the global one:

```
[0, n_owned)                 particles owned by this subdomain
[n_owned, n_local)           halo copies, authored by a neighbor
[n_local, particles_bound)   spare capacity
```

Indices are local. There is no global index; the original id already carried by the
particles is what survives migration and what output uses to reassemble a global order.

Two counters are needed, and `BaseParticles` carries both:

- `TotalRealParticles()` — the **owned** particles. All physics, all reductions,
  particle sorting and the relation build loop over this.
- `TotalLocalParticles()` — owned **plus** halo. Only `UpdateCellLinkedList` uses it. If
  the cell linked list saw only the owned particles, a particle next to a cut plane
  would lose part of its support and the interaction would be silently wrong there.

Outside a decomposed run the two are equal, which is why the change is invisible to
existing code paths.

## 5. What had to change in shared code

Each of these is a no-op when one subdomain is used.

| File | Change |
| --- | --- |
| `execution/subdomain_scope.h` | new: `MaxSubdomains`, thread-local id, `SubdomainScope`, fan-out guard |
| `execution/subdomain_runner.{h,cpp}` | new: worker pool and the sequential/threaded runner |
| `execution/subdomain_fan_out.h` | new: `fanOutOverSubdomains` / `reduceOverSubdomains` / `copyBetweenSubdomains` |
| `execution/execution_policy.h` | new `MultiDeviceExecution` and `MultiHostExecution` policies |
| `execution/base_implementation.h` | `is_updated_` becomes one flag per subdomain |
| `execution/implementation.h` | computing kernel and staging keeper become per-subdomain arrays |
| `common/sphinxsys_variable.h` | per-subdomain delegates; new `HostOnlyDiscreteVariable` |
| `particles/base_particles.{h,cpp}` | second counter `TotalLocalParticles` |
| `simple_algorithms_ck.h`, `interaction_algorithms_ck.hpp`, `update_cell_linked_list.hpp`, `update_body_relation.hpp`, `particle_sort_ck.hpp` | `exec()` bodies wrapped in the fan-out |
| `implementation_sycl.h` | `ExecutionInstance` becomes a facade over `DeviceEnvironment` |
| `sphinxsys_variable_sycl.hpp` | per-device allocation and staging |

New files:

```
shared/domain_decomposition/domain_decomposition.{h,cpp}          slab geometry, load balancing
shared/domain_decomposition/subdomain_exchange.{h,hpp}            halo + migration, policy generic
shared/domain_decomposition/domain_decomposition_dynamics.h       time loop entries
src_sycl/shared/common/device_environment_sycl.{h,cpp}            devices, shared context, queues
src_sycl/shared/particle_dynamics/particle_iterators_multi_device_sycl.h
tests/unit_tests_src/for_2D_build/domain_decomposition/...        the tests of §11
```

Build options, all default `OFF`, and mutually exclusive in practice:

```
cmake -DSPHINXSYS_USE_SYCL=ON -DSPHINXSYS_MULTI_DEVICE=ON        # multi-GPU
cmake -DSPHINXSYS_MULTI_SUBDOMAIN_HOST=ON                        # CPU debugging path
```

Each switches `MainExecutionPolicy`. Headers are globbed, so no `CMakeLists.txt` edits
are needed for the new sources.

## 6. Decomposition

A **slab decomposition**: cut planes normal to one axis (by default the longest, which
minimizes interface area), one slab per subdomain, so a subdomain has at most two
neighbors. `SubdomainMap` is a small trivially-copyable struct holding the cut planes,
captured by value into kernels rather than reached through a pointer.

Load balancing (`SlabDecomposition::rebalance`) moves the interior planes towards an
equal particle count, assuming the density is locally uniform along the split axis —
which is what makes a 1-D rebalance cheap: only per-subdomain counts are needed, no
histogram. A relaxation factor damps the oscillation a fully applied correction would
cause. A slab is never allowed to become thinner than twice the halo width, below which
a halo would reach past the adjacent subdomain and break the two-neighbor assumption.

**Known cost:** every subdomain currently allocates the cell-linked-list mesh for the
*whole* domain, most of it empty. This keeps cell indices identical across subdomains
and avoids any mesh remapping — the right trade for a first implementation — but it is
O(total cells) memory per subdomain. `Mesh::setLinearCellIndexOffset()` is the hook for
giving each subdomain only its own cells plus halo. First optimization after
correctness.

## 7. Exchange protocol

Halo exchange and migration use one mechanism, with two barriers over the subdomains
(each `fanOutOverSubdomains` call is itself a barrier):

1. every subdomain packs what it must send into its own send buffers;
2. **barrier**;
3. every subdomain *pulls* from its neighbors' send buffers into its own arrays;
4. **barrier**.

The pull direction is deliberate: a subdomain only ever writes memory it owns, so the
exchange needs no cross-subdomain atomics and no locks. The barriers alone order it —
and under the sequential runner they are satisfied trivially, which is exactly the
configuration to debug the *logic* in, before turning on the threaded runner to debug
the *synchronization*.

Packing is a gather per (variable, side); the send index lists come from a flag kernel
plus `exclusive_scan`, with a sentinel entry at `n_owned` so the scan returns the total.
Buffers are plain `DiscreteVariable`s — per-subdomain replication for free — and are
pre-touched at construction so a `DelegatedData()` from a neighbor's thread cannot race
with a lazy allocation. Halo slots are contiguous per side, so the pull writes directly
into the particle arrays at `n_owned + offset`: no receive buffer, no unpack kernel.

### Plan versus refresh

Separated because their costs differ by an order of magnitude:

- **`updateHaloPlan()`** — recompute which particles fall in a neighbor's halo band, the
  counts, the offsets and `n_local`. A scan plus a capacity check. Once per
  configuration update, right after the cell linked list is rebuilt.
- **`refreshHalo(variables)`** — re-send values along the existing plan. A pack plus a
  copy. After every stage that writes a state variable the next interaction reads.

`SyncHaloStateCK` takes a named subset precisely so a per-stage refresh moves only what
the next interaction needs.

### Migration

Ownership transfer, needed only at the advection step:

1. flag departing particles per destination side, and the complement (staying);
2. pack the departing ones;
3. compact the staying ones — gather into scratch, then copy back; an in-place
   compaction would overwrite entries still to be read;
4. **barrier**;
5. pull arrivals into the slots after the compacted particles, where they become owned.

Packing must precede compaction (packing reads the arrays, compaction rewrites them),
but no barrier is needed between them since compaction touches only local memory.

`SubdomainExchange::checkConsistency()` verifies the invariants — counts within bounds,
and every owned particle actually inside its own slab. On the host path it reads the
replicas directly; that check catching a missed migration, immediately rather than as a
slow physics drift, is most of the value of the CPU path.

## 8. Placement in the time loop

```
advection step:
    water_update_particle_position.exec();
    migrate_particles.exec();              // ownership follows the positions
    particle_sort.exec();                  // optional, per subdomain
    water_cell_linked_list.exec();         // over owned + halo
    update_halo.exec();                    // new plan + full state refresh
    water_block_update_complex_relation.exec();
    ... rebalance every few hundred steps ...

acoustic step:
    fluid_acoustic_step_1st_half.exec(dt);
    sync_halo_state.exec();                // velocity, pressure, density, ...
    fluid_acoustic_step_2nd_half.exec(dt);
    sync_halo_state.exec();
```

Setup, before any particles are generated (this ordering matters — it fixes how many
replicas each variable allocates):

```cpp
execution::subdomain_runner.initialize(4, SubdomainRunner::Mode::Sequential); // CPU path
// or, on the SYCL path:
execution::device_environment.initialize(0);   // 0 = all visible GPUs
```

The per-stage refresh is the dominant new cost and follows from a halo one cut-off deep:
a halo particle has an incomplete neighborhood, so its own update is untrustworthy and
must be replaced by the owner's value before the next stage reads it.

The alternative is a **deeper halo**: with a halo `k` cut-offs deep, `k` stages can run
between exchanges, at the cost of redundant computation on halo particles. Which wins
depends on particles per subdomain and on the interconnect. `halo_width` is a
constructor argument of `SlabDecomposition`, so the experiment is supported — but the
"compute on halo, then discard" variant additionally needs the dynamics loop bound
switched from owned to local, which is *not* wired up.

## 9. Correctness notes

- **Reductions** differ from the non-decomposed result by floating-point association.
  Regression tolerances need revisiting; quantify the effect on the CPU path first.
- **One-sided inner relation.** `UpdateRelation<Inner<...>>::incrementNeighborSize`
  registers the reverse neighbor of each pair, writing `neighbor_index_[tar_index]` for
  a target that may be a halo particle, while offsets are only scanned over
  `[0, n_owned]`. With halos present this writes out of range. The inner relation must
  be built **two-sided** under decomposition, or the scan extended to `n_local`. This is
  the most likely source of silent corruption and is **not fixed**. It should be the
  first thing the CPU path is pointed at — under ASan it is an immediate, localized
  failure rather than a drifting result.
- **Contact relations to a non-decomposed body** (walls, observers) work only if that
  body is fully replicated on every subdomain. Replication is right for small static
  bodies, but the draft does not distinguish replicated from decomposed bodies, and that
  distinction needs adding to `SubdomainExchange`.
- **Observers and I/O** run host-side and need `gatherToHost()` first.
- **Particles leaving the domain** are clamped to the end subdomains, never unowned.
- **A particle may only move to an adjacent subdomain per step.** A larger jump means
  the time step or the slab thickness is wrong; `checkConsistency()` detects the result.

## 10. What is stubbed

| Item | State |
| --- | --- |
| Replicated (non-decomposed) bodies | not distinguished from decomposed ones |
| One-sided inner relation with halos | broken as described above |
| Per-subdomain cell mesh | full-domain mesh per subdomain |
| Restart with a decomposition | untouched |
| NUMA placement on the host path | none; first-touch is incidental, `tbb::task_arena` constraints would be the fix if this path ever needs to be fast |

`scatterFromHost()` and `gatherToHost()` are now implemented (they were stubs in the
first draft): scatter groups the host arrays by owner and stages one contiguous slice
per subdomain; gather reverses it over the owned ranges only, so halo duplicates are not
written out twice.

## 11. What is actually tested

Built and run on this machine, under AddressSanitizer, UndefinedBehaviorSanitizer and
ThreadSanitizer, all clean:

- **Fan-out semantics**, in both runner modes: body runs once per subdomain; nested
  fan-out collapses onto the bound subdomain rather than deadlocking or duplicating;
  reduction combines partials correctly; exceptions propagate out of a worker; the
  thread binding is restored afterwards; a non-decomposed policy runs the body once.
- **Decomposition geometry**: ownership is an exact partition over 8000 sample positions
  (no gaps, no overlaps); positions outside the domain stay owned; halo bands match the
  slab widened by the cut-off; a particle in a neighbor's halo band is still owned by
  us; end subdomains have one neighbor.
- **Load balancing**: from a 3.86× imbalance, `rebalance` converges to 1.0005× in 40
  iterations while keeping the cut planes monotone, the outer planes pinned, and every
  slab at or above the minimum thickness even under a degenerate load.

These live in `tests/unit_tests_src/for_2D_build/domain_decomposition/` as GTest cases.
The geometry tests were additionally run here against a stub of the Eigen-backed types,
since Eigen is not installed on this machine; in the repo they compile against the real
`Vecd`/`BoundingBoxd`.

**Not tested:** anything touching `BaseParticles`, the variable replication, the
exchange itself, or SYCL. That needs the real dependencies.

## 12. Suggested bring-up order

1. Build with both options `OFF`, confirm existing tests unchanged. Every edit in §5 is
   designed to be a no-op there; this is the regression gate for the whole refactor.
2. `SPHINXSYS_MULTI_SUBDOMAIN_HOST=ON` with **1** subdomain. Exercises the fan-out, the
   per-subdomain arrays and the host replicas while the answer must still match step 1.
3. Same, **2 subdomains, sequential**, on `dambreak`. First real decomposition. Turn on
   `checkConsistency()` every step and build with ASan — this is where the one-sided
   inner relation of §9 should surface.
4. Same, **threaded**, under ThreadSanitizer. Any difference from step 3 is a
   synchronization bug, and the barrier structure of §7 is where to look.
5. `SPHINXSYS_MULTI_DEVICE=ON`, one GPU, then two. By this point the decomposition logic
   is already known good, so a failure here is device-specific: USM lifetime, queue
   ordering, or peer access.
6. Then: per-subdomain mesh, deeper halos, and overlapping the exchange with interior
   computation — the barrier structure already isolates where that overlap goes.
