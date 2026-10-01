# Domain decomposition in SPHinXsys: multi-GPU, and a CPU path to debug it

**Status: draft.** The parts that can be built without the project's dependencies are
built and tested here (see §11); the rest — anything touching particles, SYCL or Eigen —
is a reviewed design sketch with the plumbing written out, not tested code. This machine
has neither a SYCL toolchain nor Simbody/TBB/Eigen. §10 lists what is stubbed and what
is likely to break first.

## 1. Two backends, one implementation

There are two decomposed execution policies, and they share essentially all of their
code:

| | `MultiDevicePolicy` = `DecomposedExecution<SYCLDevicePolicy>` | `MultiHostPolicy` = `DecomposedExecution<ParallelPolicy>` |
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
| loop range | `LoopRangeCK<DecomposedExecution<P>, ...>`, built per subdomain inside `particle_for` |

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

Parallelism over subdomains is expressed **inside `particle_for` / `particle_reduce`**, in the
overloads for `LoopRangeCK<DecomposedExecution<P>, ...>` (`particle_iterators_ck.h`). The
algorithms' `exec()` bodies are the plain single-domain ones:

```cpp
// StateDynamics::exec(), unchanged for every policy
particle_for(LoopRangeCK<ExecutionPolicy, RangeIdentifier>(*this->identifier_),
             kernel_implementation_, dt);

// the decomposed overload
template <class PolicyType, class Identifier, class KernelImplementationType>
void particle_for(const LoopRangeCK<DecomposedExecution<PolicyType>, Identifier> &loop_range,
                  KernelImplementationType &implementation, Real dt)
{
    fanOutOverSubdomains(DecomposedExecution<PolicyType>{}, [&]()
                         { particle_for(loop_range.onCurrentSubdomain(), implementation, dt); });
}
```

Two things make this possible. The computing kernel is fetched *inside* the loop function
(`implementation.getComputingKernel()`, resolved on `currentSubdomainID()`), so each subdomain
gets its own kernel. And the loop range of a decomposed policy is deferred
(`DeferredLoopRangeCK` in `loop_range.h`): it only stores the identifier, and
`onCurrentSubdomain()` builds the base policy's range inside the fan-out, where
`DelegatedData()` addresses that subdomain's replica. A range built on the host thread would
bind every subdomain to the replica of subdomain 0.

`DecomposedExecution<P>` means "P, plus a fan-out" and nothing else: every other operation
(`DelegatedData`, kernel allocation, `exclusive_scan`, `particle_for` on an `IndexRange`, the
copy between subdomains) is forwarded to `P` by a one-line overload on
`DecomposedExecution<P>`. Those forwarders are not optional: a catch-all template
`foo(const ExecutionPolicy &)` is an exact match for `DecomposedExecution<SYCLDevicePolicy>` and
would otherwise win over `foo(const SYCLDevicePolicy &)`, silently running the host branch.

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
| `execution/execution_policy.h` | `DecomposedExecution<P>`, with `MultiDevicePolicy`, `MultiHostPolicy`, `SequencedMultiHostPolicy` |
| `execution/base_implementation.h` | `is_updated_` becomes one flag per subdomain |
| `execution/implementation.h` | computing kernel and staging keeper become per-subdomain arrays |
| `common/sphinxsys_variable.h` | per-subdomain delegates; new `HostOnlyDiscreteVariable` |
| `particles/base_particles.{h,cpp}` | second counter `TotalLocalParticles` |
| `loop_range.h`, `particle_iterators_ck.h` | deferred loop range and the loop-level fan-out for `DecomposedExecution<P>` |
| `update_cell_linked_list.hpp`, `update_body_relation.hpp`, `particle_sort_ck.hpp` | `exec()` bodies wrapped in the fan-out (loops on an `IndexRange` with a captured kernel) |
| `particles/base_particles.h`, `interaction_algorithms_ck.hpp` | `HaloRefresher` hook; interactions refresh their interact variables before the interaction step |
| `implementation_sycl.h` | `ExecutionInstance` becomes a facade over `DeviceEnvironment` |
| `sphinxsys_variable_sycl.hpp` | per-device allocation and staging |

New files:

```
shared/domain_decomposition/domain_decomposition.{h,cpp}          slab geometry, load balancing
shared/domain_decomposition/subdomain_exchange.{h,hpp}            halo + migration, policy generic
shared/domain_decomposition/domain_decomposition_dynamics.h       time loop entries
src_sycl/shared/common/device_environment_sycl.{h,cpp}            devices, shared context, queues
tests/unit_tests_src/for_2D_build/domain_decomposition/...        the tests of §11
```

One build option, default `OFF`. Which decomposition it selects follows from the backend:

```
cmake -DSPHINXSYS_USE_SYCL=ON -DSPHINXSYS_DECOMPOSITION=ON       # multi-GPU
cmake -DSPHINXSYS_DECOMPOSITION=ON                               # CPU debugging path
```

`SPHINXSYS_DECOMPOSITION` wraps the backend policy (`SYCLDevicePolicy` with SYCL,
`ParallelPolicy` otherwise) in `DecomposedExecution<>` to form `MainExecutionPolicy`.
Headers are globbed, so no `CMakeLists.txt` edits are needed for the new sources.

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

1. flag departing particles per destination side; from the same flags derive the new
   owned count `n_new`, the list of departing slots below `n_new` (the holes) and the
   list of staying particles at or above `n_new` (the donors), both ascending. There
   are exactly as many donors as holes;
2. pack the departing ones;
3. remove them by copying the k-th donor into the k-th hole. This is the "swap with
   the last real particle" used by particle deletion, driven by scans instead of an
   atomic counter, so it is deterministic and costs O(departing) rather than a full
   compaction; no slot is both read and written;
4. **barrier**;
5. pull arrivals into the slots after `n_new`, where they become owned. This is the
   particle generation side of the exchange: the state of a new particle comes from
   the neighbor's send buffer rather than from another particle of the same array.

Packing must precede the removal (packing reads the departing slots, the removal
overwrites them), but no barrier is needed between them since the removal touches
only local memory.

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
    update_halo.exec();                    // new plan + full state refresh; publishes n_local
    water_cell_linked_list.exec();         // over owned + halo
    water_block_update_complex_relation.exec();
    ... rebalance every few hundred steps ...

acoustic step:
    fluid_acoustic_step_1st_half.exec(dt); // refreshes the halo pressure itself, before
    fluid_acoustic_step_2nd_half.exec(dt); // its interaction step; the velocity likewise
```

Setup. The subdomain runner must be initialized before any particles are generated,
since that fixes how many replicas each variable allocates. The number of subdomains is
a run time choice of the `SPHSystem`: the command line option `--subdomains=N`, or
`sph_system.setNumberOfSubdomains(N)` right after construction; both initialize the
runner in its sequential mode. The threaded host runner is selected by calling
`execution::subdomain_runner.initialize(N, SubdomainRunner::Mode::Threaded)` at the same
point. On the SYCL path the device environment is initialized there as well
(`execution::device_environment.initialize(0)`, 0 = all visible GPUs).

The case then drives everything through one `BodyDecomposition<MainExecutionPolicy>` and
the dynamics built on it (`body_decomposition.h`, `domain_decomposition_dynamics.h`). All
of them are no-ops when the policy is not a `DecomposedExecution<>`, so the same source
serves the decomposed and the plain build:

```cpp
// after every dynamics of the body and after its output variables: fixes the exchange set
auto &decomposition = main_methods.addDecomposition(water_block, 0 /* split axis */);
decomposition.addExchangeVariable<Real>("Pressure");             // read at the neighbors,
decomposition.addExchangeVariable<Matd>("LinearCorrectionMatrix"); // neither evolving nor written
decomposition.addSubdomainIDToWrite(body_state_recorder);        // optional owner tag in the vtp
auto &update_halo = main_methods.addGeneralDynamics<UpdateHaloCK>(decomposition);
auto &migrate_particles = main_methods.addGeneralDynamics<MigrateParticlesCK>(decomposition);
auto &sync_volume = main_methods.addGeneralDynamics<SyncHaloStateCK>(decomposition);
sync_volume.addVariable<Real>("VolumetricMeasure");

decomposition.scatterFromHost();   // once, before the first dynamics on the body
...
decomposition.gatherToHost();      // around every host side access: output, restart
body_state_recorder.writeToFile();
decomposition.finishHostAccess();
```

The exchange set of a `BodyDecomposition` is the evolving variables plus the variables
registered for output when it is constructed, plus `addExchangeVariable()` calls made
before `scatterFromHost()`.

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
  a target that may be a halo particle. The lists, the offsets and the scan therefore
  cover the local range `[0, n_local]` (`UpdateRelation<Inner>::updateOnCurrentDevice`);
  outside a decomposed run the two counts are equal and nothing changes.
- **State that survives an advection step must be an evolving variable.** The particle
  sort permutes the evolving variables only. `Force` is assigned by the second acoustic
  half step and accumulated onto by the next first half step, across the sort, and was
  not evolving in the CK acoustic step; a single domain scrambles it deterministically,
  two subdomains scramble it differently, and the runs diverge at the first sort. It is
  evolving now. The same reasoning applies to any variable a case adds: whatever a
  particle carries from one advection step to the next has to be in the evolving set,
  which is also the set that migrates.
- **Halo refresh points.** Each quantity read at the neighbors is refreshed right after
  the stage that writes it. The pressure and the velocity are refreshed by the acoustic
  steps themselves: an interaction algorithm calls `BaseParticles::refreshHalo()` with
  its `interact_variables_` right before its interaction step (a no-op unless a
  `SubdomainExchange` installed itself as the `HaloRefresher` of that body). The volume
  and the correction matrix change once per advection step and are refreshed there by the
  case (`SyncHaloStateCK`), rather than being listed as interact variables and re-sent
  every acoustic step.
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
| Threaded host runner | lazy replica creation is serialized (`replicaCreationMutex()`), but `DiscreteVariable::reallocateData` grows the shared capacity and every subdomain's replica from whichever subdomain thread triggers it, while the other threads still hold the old pointers in their computing kernels. Neighbor list growth during a threaded run therefore crashes; use the sequential runner until capacities are per replica |
| Per-subdomain cell mesh | full-domain mesh per subdomain |
| Restart with a decomposition | works: the restart output gathers to the host, the restart read precedes `scatterFromHost()` |
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

Since then, with the real dependencies (2026-09-10): `test_2d_dambreak_sycl` itself runs
decomposed, `--subdomains=2` in a `SPHINXSYS_DECOMPOSITION=ON` build, through
`BodyDecomposition` and the loop dynamics of §8, with its observer, its restart output
and both of its dynamic time warping regression tests passing, and a restart from the
files of the decomposed run completing. The corresponding ctest entries exist in that
build. The plain build is unchanged: with one subdomain, or without the option, the
energy record matches the plain run bit for bit for as long as the plain run is itself
reproducible (about two seconds of physical time; the one-sided inner relation appends
neighbors with atomic counters, so the summation order and, from there, the trajectory
varies from run to run, on this branch and independently of the decomposition).

Two defects found on the way, both only visible once a subdomain grows its neighbor
list after the others have built theirs in the same step: replica reallocation used to
discard the contents of every replica, and the relation update dynamics did not
register its kernel with the relation, so the other subdomains kept pointers into the
freed replica. Both are fixed; the host and the device replicas now keep their contents
on growth, and the update kernels are invalidated like the interaction kernels.

## 12. Bring-up order and where it stands

1. Build with the option `OFF`, confirm existing tests unchanged. Every edit in §5 is
   designed to be a no-op there; this is the regression gate for the whole refactor.
   **Done.**
2. `SPHINXSYS_DECOMPOSITION=ON` (without SYCL) with **1** subdomain. Exercises the fan-out, the
   per-subdomain arrays and the host replicas while the answer must still match step 1.
   **Done**, on `test_2d_dambreak_sycl --subdomains=1`.
3. Same, **2 subdomains, sequential**, on `dambreak`. First real decomposition.
   **Done**, `--subdomains=2`, regression tests and restart included (§11).
4. Same, **threaded**, under ThreadSanitizer. Any difference from step 3 is a
   synchronization bug, and the barrier structure of §7 is where to look. Open; the
   threaded runner is still restricted by the reallocation issue of §10.
5. `SPHINXSYS_USE_SYCL=ON -DSPHINXSYS_DECOMPOSITION=ON`, one GPU, then two. By this point the decomposition logic
   is already known good, so a failure here is device-specific: USM lifetime, queue
   ordering, or peer access. Open; no SYCL compiler on the development machine.
6. Then: per-subdomain mesh, deeper halos, and overlapping the exchange with interior
   computation — the barrier structure already isolates where that overlap goes.
