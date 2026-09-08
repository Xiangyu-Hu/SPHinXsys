# Multi-GPU / distributed-memory SPHinXsys — design notes

Status: **draft skeleton, not compiled, not tested.** No MPI code existed in the tree
before this (`grep -ri mpi src/` returns only false positives such as `compile`/`simple`).

## 1. What the CK (computing-kernel) path actually does per step

Verified by reading `shared_ck/particle_dynamics/configuration_dynamics/` and the
SYCL dam-break case (`tests/tests_sycl/3d_examples/test_3d_dambreak_sycl/dambreak.cpp`).

`UpdateCellLinkedList::exec()` runs **every step** and is a counting sort, not a radix sort:

1. `clearAllLists`   over `[0, number_of_cells_)`      — zero offsets/sizes
2. `incrementCellSize` over `[0, total_real_particles)` — atomic histogram of particles per cell
   (the histogram is temporarily stored in `particle_index_` — deliberate buffer reuse)
3. `exclusive_scan(particle_index_ -> cell_offset_, number_of_cells_ + 1)` — CSR row pointers
4. `updateCellList`  over `[0, total_real_particles)`  — atomic fetch-add scatter of particle
   indices into `particle_index_[cell_offset_[cell] + local_rank]`

`ParticleSortCK` (which *is* the radix sort, via `SortMethod<ExecutionPolicy>::type`) is a
**separate, periodic** operation — the dam-break case calls it every 100 iterations. It
permutes all evolving particle variables into Morton order purely for memory locality.

Consequence: the per-step neighbour-search cost is histogram + scan + scatter. The radix
sort is amortised over ~100 steps, so swapping it for CUB is a much smaller win than a
per-step sort would imply. Profile before spending effort there.

## 2. The existing single-rank analogue of a halo exchange

`Ghost<PeriodicAlongAxis>` + `PeriodicConditionUsingGhostParticles`
(`shared/particle_dynamics/general_dynamics/domian_bouding/ghost_bounding.{h,cpp}`)
already implements the exact structure a halo exchange needs, and — importantly — already
splits it along the same seam MPI needs:

| periodic ghost (existing, single rank)                | MPI halo (this work)                          |
|-------------------------------------------------------|-----------------------------------------------|
| `reserveGhostParticles` — allocate ghost index block   | same, sized per neighbour rank                |
| `CreatPeriodicGhostParticles` — find boundary particles, `updateGhostParticle(slot, src)` | build the *send list* per neighbour rank; exchange counts |
| `UpdatePeriodicGhostParticles` — `copyFromAnotherParticle(ghost, sorted_id_[i])` each step | `MPI_Isend`/`MPI_Irecv` the payload each step |

"Creation" = establish the communication pattern (expensive, occasional).
"Update"   = refresh the payload (cheap, every step).
That is exactly the classic halo-exchange decomposition.

Caveat: this code is legacy `shared/` path only — `execution::ParallelPolicy`, `std::mutex`,
and the old `InsertListDataEntry` cell-list API. **There is no periodic BC in `shared_ck/`**,
so there is no CK/SYCL template to copy; it has to be written.

## 3. Migration primitives already exist

`BaseParticles` (`shared/particles/base_particles.cpp`) provides exactly what particle
migration needs:

- `switchToBufferParticle(index)` — remove a real particle (swap with last, decrement count)
- `createRealParticleFrom(index)` — promote a buffer particle to real (increment count)
- `allocateGhostParticles(n)` / `updateGhostParticle(ghost, src)` — ghost block management

## 4. Where the new code goes

`src/shared/shared_ck/domain_decomposition/` (new):

- `mpi_environment.{h,cpp}`   — MPI init/finalize, rank/size, guarded by `SPHINXSYS_USE_MPI`
- `domain_partition.{h,cpp}`  — Cartesian split of the system bounding box, neighbour-rank
                                topology, owner lookup. Pure geometry: unit-testable with no MPI.
- `halo_exchange_ck.{h,hpp}`  — per-step halo exchange dynamics (CK style, policy-templated)
- `particle_migration_ck.{h,hpp}` — ownership transfer (drafted, not compiled — section 8)

## 5. Insertion point in the time loop

The time loop lives in *user* code (see `dambreak.cpp`), so a distributed case wires this in
explicitly. Order matters — the halo must be populated **before** the cell list is rebuilt,
because the cell list is what neighbour search reads:

```
water_update_particle_position.exec();
particle_migration.exec();        // NEW: hand off particles that left this subdomain
halo_exchange.exec();             // NEW: receive neighbours' boundary particles
water_cell_linked_list.exec();    // existing — now also bins halo particles
water_block_update_complex_relation.exec();   // existing, unchanged
```

## 6. The central design decision: how halo particles enter the cell list

`UpdateCellLinkedList::exec()` loops `IndexRange(0, total_real_particles)`, so anything
above `total_real_particles_` is invisible to neighbour search.

**Option A (chosen here): halo particles are appended as real particles.**
Set `total_real_particles_ = n_owned + n_halo` before the cell-list build. Cell list,
relations and neighbour search then include halo particles with *zero* modification.
Physics must not write back to halo slots, so physics dynamics loop over a new
`sv_total_owned_particles_` bound instead. `LoopRangeCK<Policy, SPHBody>` already has a
constructor taking an explicit `SingleVariable<UnsignedInt> *`, which is the hook.

- pro: no changes to the hot neighbour-search path
- con: every physics call site must be given the "owned" identifier; forgetting one is a
  silent correctness bug (halo slots get integrated, then overwritten next exchange)

**Option B: keep `total_real_particles_ = n_owned`, extend the cell-list build.**
Modify `UpdateCellLinkedList::exec()` to also loop the halo block in `incrementCellSize`
and `updateCellList`.

- pro: physics call sites untouched; halo cannot leak into physics by accident
- con: touches the hot path and the shared `UpdateCellLinkedList` used by every case

Option A is drafted below. Option B is the safer choice if silent physics bugs are the
bigger worry — the decision should be made before this goes further.

### Where the halo particles physically sit

Not in the `Ghost<>` block. `allocateGhostParticles()` hands out slots at the *fixed*
offset `particles_bound_`, which sits above the buffer region; the owned particle count
moves at runtime, so a fixed ghost block leaves a gap between the owned particles and the
halo. Setting `total_real_particles_` past that gap would make the cell linked list bin the
stale buffer particles in between, producing wrong neighbours.

So the halo occupies the **buffer region directly above the owned particles**,
`[n_owned, n_owned + n_halo)` — the same region the emitter/inflow boundaries write into via
`createRealParticleFrom`. A distributed case must therefore reserve buffer particles sized
for its expected halo. Each step: reset `total_real_particles_` to `n_owned`, refill the
halo, then set it to `n_owned + n_halo`.

## 7. Overlap opportunity

`ExecutionEvent`-based async dispatch (PR #456) already exists on the SYCL path. The
standard optimisation is: post the halo `MPI_Irecv`/`MPI_Isend`, run the interior
particles' kernels (those farther than one cutoff radius from any subdomain boundary),
then wait on the requests before the boundary particles' kernels. Not attempted in this draft.

## 8. What is actually implemented

| file | state |
|------|-------|
| `domain_partition.{h,cpp}` | **complete and tested** — 24 assertions pass, clean under `-Wall` |
| `mpi_environment.{h,cpp}`  | complete, not compiled (needs an MPI toolchain) |
| `halo_exchange_ck.{h,hpp}` | structure complete; host packing path written, device path is a TODO |
| `particle_migration_ck.{h,hpp}` | structure complete, **not compiled, not run** — see below |
| CMake `SPHINXSYS_USE_MPI`  | option, `find_package(MPI)`, `MPI::MPI_CXX` linkage, compile definition |

`multi_gpu_standalone_test/` compiles the real `domain_partition.cpp` against a minimal
type shim, so the topology logic is verifiable without Simbody/Eigen/TBB/Boost/spdlog.
Nothing else here has been compiled — the full library needs a dependency stack that is
not installed on this machine, and none of it has been run on more than one rank.

### `particle_migration_ck` design notes

Reuses `HaloVariablePacker`/`TypedHaloVariablePacker` from `halo_exchange_ck.h` — packing
N particles' worth of one variable into a buffer is the same operation for a migration
payload as for a halo payload. Case code registers variables via `addToMigrate<T>(name)`,
mirroring `HaloExchangeCK::addToExchange` (deliberately not automatic: the migration
payload should be the same or a superset of what halo exchange sends, and there is no way
to enumerate "everything a case cares about" generically without pulling in the private
`all_state_data_` assemble in `BaseParticles`).

Per-step order inside `exec()`: `buildMigrationLists()` (scan owned particles, sort each
one whose `DomainPartition::ownerRank` differs from this rank into a per-neighbour send
list; errors out if the owner isn't one of the ~26 precomputed neighbours, i.e. dt moved a
particle more than one subdomain width) → `negotiateCounts()` → `packSendBuffers()` (must
happen **before** any particle is removed, while indices are still valid) →
`removeMigratedParticles()` (swap-with-last via the existing `switchToBufferParticle`,
processed in descending index order so a removal's swap partner is never itself pending
removal) → `exchangeAndAppend()` (MPI transfer, then one `createRealParticleFrom(0)` per
arrival to get a correctly-numbered new slot — the copied placeholder state is immediately
overwritten by `unpack()`).

Ends by writing `particles_.TotalRealParticles()` into the *same*
`SingleVariable<UnsignedInt>` that `HaloExchangeCK::svTotalOwnedParticles()` exposes,
passed into the constructor — this is the integration seam with the halo exchange that
runs immediately afterwards each step (see section 5): without it, `HaloExchangeCK`
would keep reading the owned count from setup time and never see migration's effect.

Known gaps: no compile/run yet (same dependency-stack blocker as everything but
`domain_partition`); `createRealParticleFrom(0)` assumes the rank owns at least one real
particle already when the first arrival lands, which is false for the edge case of a rank
losing every particle in one step (memory-safe — `particles_bound_`-sized arrays make slot
0 always addressable — but reads uninitialized data as the throwaway copy source); no
attempt yet at overlapping the migration and halo MPI rounds.

## 9. Next steps, in order

1. ~~`particle_migration_ck`~~ — drafted (see above), not yet compiled or run.
2. Decide option A vs option B (section 6) before more code depends on the choice.
3. `MPI_Allreduce` on the CFL time step, or ranks silently desynchronise (see section 8) —
   `MPIEnvironment::allReduceMin` already exists for this; it just isn't wired into a
   `ReduceDynamicsCK` time-step call site yet.
4. Device-side send-list construction and packing — the host path in `buildSendLists()`/
   `buildMigrationLists()` is a serial scan over every owned particle and would dominate
   on GPU.
5. Only then: the interior/boundary overlap in section 7.

## 10. Not addressed in this draft

- Multi-level `AdaptiveSmoothingLength` cell lists across ranks
- Load balancing / repartitioning as particles migrate
- Global reductions (`ReduceDynamicsCK`, e.g. the CFL time step) need an `MPI_Allreduce`;
  currently each rank would compute its own local `dt` and desynchronise
- Distributed I/O (`BodyStatesRecordingToVtpCK` writes per-rank files today)
- GPU-aware MPI (sending USM device pointers directly vs. staging through host)
