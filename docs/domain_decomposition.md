# Domain decomposition

This document describes the domain-decomposition design and its implementation in
SPHinXsys. It is the current reference for the feature; it intentionally omits the
development and review chronology.

## Purpose and scope

Domain decomposition assigns particles to subdomains so that their physics can be
computed independently, with neighboring particle data exchanged through halos.
The same policy-generic implementation supports a host execution path and a SYCL
device path:

| Backend | Decomposed policy | Subdomain execution |
| --- | --- | --- |
| Host | `DecomposedExecution<ParallelPolicy>` (`MultiHostPolicy`) | Sequential host runner by default; a threaded runner is also available |
| SYCL | `DecomposedExecution<SYCLDevicePolicy>` (`MultiDevicePolicy`) | One device per subdomain |

The host path is useful for debugging the decomposition and exchange logic without
requiring SYCL hardware. The device path uses a shared SYCL context for its device
queues and cross-device copies. This is a single-process, single-node design; it is
not an MPI or multi-node implementation.

The decomposition is selected by the execution policy, not by case-level
preprocessor branches. `SPHINXSYS_USE_SYCL` selects the backend, and
`SPHINXSYS_DECOMPOSITION` determines whether `MainExecutionPolicy` wraps it in
`DecomposedExecution`.

| `SPHINXSYS_USE_SYCL` | `SPHINXSYS_DECOMPOSITION` | `MainExecutionPolicy` |
| --- | --- | --- |
| OFF | OFF | `ParallelPolicy` |
| OFF | ON | `DecomposedExecution<ParallelPolicy>` |
| ON | OFF | `SYCLDevicePolicy` |
| ON | ON | `DecomposedExecution<SYCLDevicePolicy>` |

## Build and run configuration

The decomposition option is disabled by default:

```sh
cmake -DSPHINXSYS_DECOMPOSITION=ON ...
cmake -DSPHINXSYS_USE_SYCL=ON -DSPHINXSYS_DECOMPOSITION=ON ...
```

The first command selects host subdomains; the second selects the SYCL device path.
The runtime subdomain count is limited by `execution::MaxSubdomains` (currently 8).
Set it before creating bodies or particle variables, since the count determines the
number of replicas:

```sh
simulation --subdomains=2
```

Alternatively, call `SPHSystem::setNumberOfSubdomains(N)` before creating bodies.
For host builds, that API and the command-line option select the sequential runner
by default. The threaded host runner can be selected through
`execution::subdomain_runner.initialize(N, SubdomainRunner::Mode::Threaded)`.
Do not change the subdomain count after variables have allocated replicas.

When the library is built without `SPHINXSYS_DECOMPOSITION`, requesting more than
one subdomain emits a warning; the main execution policy remains non-decomposed.

## Execution model

`execution::SubdomainScope` binds a host thread to a subdomain using a thread-local
ID. Particle variables, counters, computing kernels, freshness flags, and device
queues resolve their per-subdomain replica through that ID. Physics kernels can
therefore use the same variable access pattern for decomposed and non-decomposed
policies.

CK loops use deferred `LoopRangeCK` specializations for decomposed policies. The
range and computing kernel are resolved only after a subdomain has been bound.
`particle_for` fans out over subdomains; `particle_reduce` computes per-subdomain
partial results and combines them in subdomain order. Calls nested inside an
existing fan-out execute on the already-bound subdomain instead of starting another
fan-out.

Other operations, including index-range loops, scans, and computing-kernel
allocation, forward to the wrapped backend policy. These forwarders are important:
generic overloads can otherwise match `DecomposedExecution<P>` more closely than
backend-specific overloads.

For the host backend, sequential mode visits subdomains in order on the calling
thread. Threaded mode uses persistent worker threads and barriers at fan-out
boundaries. A difference between one and multiple subdomains points to a
decomposition or exchange issue; a difference between sequential and threaded
execution points to a synchronization issue. Floating-point reductions can differ
from a non-decomposed run because the association order changes.

## Decomposition geometry

The current geometry is a one-dimensional slab decomposition. Cut planes are normal
to one axis; by default, the longest domain extent is selected to reduce interface
area. A `SubdomainMap` stores the cut planes and provides owner, neighbor, and halo
queries. A subdomain has at most two neighbors. Positions outside the system bounds
are assigned to the first or last subdomain rather than becoming unowned.

The initial slabs are equally spaced. `SlabDecomposition::rebalance()` moves interior
cut planes towards equal owned-particle counts using a relaxed correction. This
estimates the position of a balanced cut assuming density is locally uniform along
the split axis; it does not build a particle histogram. The minimum slab thickness
is twice the halo width, preserving the adjacent-neighbor halo assumption. The
constructor warns if the requested subdomain count makes the initial slabs thinner
than this minimum.

## Particle layout and exchange

Each subdomain stores its own particle replicas:

```text
[0, n_owned)                 particles owned by this subdomain
[n_owned, n_local)           read-only halo copies from neighbors
[n_local, particle_bound)    spare capacity
```

`TotalRealParticles()` is the owned count. Physics loops, reductions, and particle
sorting use owned particles. `TotalLocalParticles()` includes halo copies and is
used where neighbor construction must see the halo. Outside a decomposed run, the
two counts are equal.

The exchange protocol is a pack-then-pull operation:

1. Each subdomain packs its outbound particle data into its own send buffers.
2. A fan-out boundary ensures all send buffers are ready.
3. Each subdomain pulls data from neighboring send buffers into its own arrays.
4. A second fan-out boundary completes the exchange.

Each subdomain writes only its own destination arrays. Host copies use ordinary
memory copies; SYCL copies use the device environment. Send-buffer replicas are
created before exchange so cross-subdomain reads do not race with lazy allocation.

The halo plan and halo values are separate:

- `updateHaloPlan()` rebuilds send indices, counts, offsets, and the local-particle
  count after particle positions or ownership have changed. It then refreshes the
  evolving variables.
- `refreshHalo(variables)` copies only the specified variable set through the
  existing plan. The exchange buffers are created from evolving and registered
  interaction variables when the exchange is first created; explicit refreshes can
  select from that exchange set. Interaction dynamics invoke this for their
  registered interaction variables before the interaction loop.

### Migration

Migration transfers ownership when particles cross a cut plane. The implementation
flags departing particles, packs their state, fills the resulting holes from
staying particles at the tail of the owned range, and pulls arrivals from
neighboring send buffers. The state that must follow a particle across an ownership
change must be included in `EvolvingVariables()`, because that is the set migrated
and sorted.

A particle is expected to move by no more than one adjacent subdomain in a step.
Larger jumps violate the current slab-neighbor exchange assumption and should be
prevented by the timestep and slab sizing.

## Particle-variable responsibilities

The current `BodyDecomposition` initializes the exchange-variable set lazily from
the body's evolving variables and registered interaction variables
(`AllInteractVariables()`). Variables needed at neighboring particles must be
registered as interaction variables or included in an explicit halo refresh.

Output-only variables have a different role. They do not need to be scattered or
exchanged merely because they are written to output; they must be included in the
host-side gather when their values are needed there. If a per-particle value must
follow a particle during sorting or migration, it is not output-only in the
ownership sense and must be managed as particle state accordingly.

## Setup and lifecycle

The intended lifecycle is:

1. Select the execution policy through the build configuration.
2. Set the runtime subdomain count before creating bodies or particle variables.
3. Create all bodies and construct the relevant dynamics so their interaction
   variables are registered.
4. Ensure each body has its `BodyDecomposition` and finalize exchange state only
   after the required methods and variables are known.
5. Scatter the initial host particle data before the first decomposed dynamics.
6. Run dynamics, migration, halo-plan updates, and output with the ordering
   appropriate to the case.

The case should not need to select decomposition separately: the main execution
policy is the decomposition choice. In the current solver/container path,
`SPHSolver::getMainMethodContainer()` attaches decompositions to bodies already
registered with the `SPHSystem`; exchange state is created later on first use.
`SPHSolver::getTimeStepper()` scatters initial particles for a non-restart run.
Restart loading uses the I/O helper's scatter after reading host data.

**Lifecycle limitation:** direct construction of a dynamics method can bypass the
main-method container. In that path, the current automatic attachment is not
guaranteed to have created the body's decomposition before the method runs.
Likewise, the exchange set must not be finalized before all relevant methods have
registered their interaction variables. The desired invariant is that all bodies
and relevant methods are registered, and all decomposition state is ready, before
the first dynamics execute. A policy-driven initialization path shared by container
and directly constructed methods is needed to guarantee that invariant generally.

## Dirty-flag and refresh semantics

Halo refresh is controlled by two independent pieces of information:

- **Which particles move** is fixed by the halo plan. `updateHaloPlan()` builds
  per-side send index lists containing only owned particles inside a neighbor's
  halo band. Packing and pulling loop only over these lists.
- **Which variables move** is selected by the per-variable dirty flag. Pack and
  pull skip a variable unless it is dirty and in the requested set.

A refresh therefore transfers halo-band particles of dirty, requested variables,
not whole arrays. No per-particle dirty tracking is needed.

The dirty flag is a single flag per variable object, shared by all subdomains.
Decomposition runs the same operation on different data, so all replicas of a
variable become dirty and clean together. The flag is host-only state:

- Marking happens in the host-side `setupDynamics` of `StateDynamics` and
  `InteractionDynamicsCK` (decomposed policies only), on
  `to_be_interact_variables_`, before any per-subdomain fan-out.
- `updateHaloPlan()` and `migrateParticles()` mark the evolving variables dirty.
- `refreshHalo()` and `migrateParticles()` clean the variables once, after the
  pull barrier. Fan-out bodies only read the flag.

Dirty means "needs publishing"; the registered interaction set means "needed by
the consumer". A refresh handles the intersection and cleans what it published;
a dirty variable outside the requested set stays dirty for a later refresh.

`migrateParticles()` returns early unless the position variable is dirty. The
migration trigger therefore relies on position being marked dirty after
advection; the ordering of migration relative to the cell-list update (which
calls `updateHaloPlan()` and cleans the flags) is called out in
`TimeStepper::incrementIterationStep()` and has not been verified here.

## Comparison with common single-domain decomposition approaches

Sources: LAMMPS documentation (partitioning, communication, neighbor lists,
`comm_modify`, `balance`) and AMReX particle documentation. The OpenFPM primary
paper, DualSPHysics multi-GPU pages and LIGGGHTS documentation could not be
retrieved; no claims are made about them, and any DEM remarks are inferred from
LAMMPS. The periodic-boundary observation covers only this module.

### Partitioning

LAMMPS defaults to a regular 3-D brick grid (up to 6 neighbors), with `tiled`
communication and recursive coordinate bisection (RCB) for non-uniform density.
AMReX uses an arbitrary `BoxArray` plus `DistributionMapping`. The slab here is
the simplest member of the family: a 1-D cut with at most two neighbors, chosen
on purpose behind the `SubdomainMap` interface so that it can be replaced later.
It corresponds to forcing a 1xNx1 processor grid in LAMMPS, which its
documentation calls suboptimal because thin slices increase communication.
Expect worse scaling than brick or tiled schemes as the subdomain count grows or
for domains elongated across the split axis.

### Owned plus halo layout

`[0, n_owned)` owned and `[n_owned, n_local)` halo matches the LAMMPS owned-then-
ghost array. `inHaloBandOf()` plays the role of the cutoff-widened ghost region.
Unlike LAMMPS, ghosts are not used for periodic boundaries; positions outside the
domain are clamped to the end subdomains.

### Exchange protocol

LAMMPS performs forward (owner to ghost) and reverse (ghost to owner, for example
summed forces) communication in staged axis sweeps over reusable send lists. Here
the exchange is pack, barrier, pull, barrier over send lists built by
`updateHaloPlan()`.

- Single exchange axis: no corner double-communication to handle.
- Pull instead of push: a subdomain writes only memory it owns, so no atomics or
  locks are needed. This suits shared memory and one SYCL context, and avoids the
  buffer-ownership races of a naive MPI-style push port.
- Selectivity: LAMMPS restricts communication by field set and, with
  `comm_modify mode multi`, by cutoff collection. Here the halo plan selects the
  particles and the per-variable dirty flag selects the variables (see the
  dirty-flag section). The two are comparable in kind.
- No reverse communication exists. Any dynamics that writes into halo particles
  and expects the owner to receive the result would need one.

### Migration

LAMMPS migrates atoms only on reneighboring steps and applies periodic wrap at the
same time. AMReX `Redistribute()` has a local mode for bounded movement and a
global mode for arbitrary jumps. `migrateParticles()` uses hole-fill from the tail
plus appended arrivals. It is local only: a particle that skips a subdomain is
detected by `checkConsistency()` as a count mismatch, not handled, and there is no
global fallback, for example after a large rebalance move or an aggressive time
step.

### Ordering and load balancing

Migrate, rebuild spatial index, refresh halo, then interact is the standard order
in both systems, with infrequent rebalancing. `rebalance()` moves cut planes from
cumulative counts, assuming uniform density inside each slab, with damping. This
is equivalent to LAMMPS `balance shift`, the weaker of its two dynamic tiers, but
it has no intra-slab histogram and no escalation to RCB. It is weakest for voids
and fronts such as free surfaces and dam-break fronts.

### Spatial index

LAMMPS stores only the neighbor bins overlapping the local subdomain extended by
the cutoff. Here the cell mesh covers the whole domain on every subdomain, so
mesh memory and per-cell update cost do not shrink with subdomain count. This is
the clearest gap relative to established practice.

### Lessons

1. Per-subdomain cell mesh limited to slab plus halo band.
2. A reverse-communication audit, and an explicit rule for halo writes.
3. A safety net for multi-hop migration: an explicit check, or a global mode.
4. Cost-aware balancing, with an intra-slab density profile or RCB if slabs prove
   too coarse.
5. Periodic ghost images if periodic boundaries must work with decomposition.
6. Keep pull-based, barrier-ordered exchange; it is a strength in the
   shared-memory context.

## Typical ordering in a simulation

For a case that changes particle positions during advection, the decomposition
operations must preserve this dependency order:

```text
update particle positions
    -> migrate ownership
    -> sort owned particles, if required
    -> update halo plan and refresh the full evolving state
    -> update cell linked list over owned plus halo particles
    -> rebuild body relations
    -> run interaction stages (refresh their registered interaction variables)
```

Rebalancing is optional and should be done infrequently; if it moves cut planes,
particles must be migrated before subsequent use.

## I/O and restart

Host-side output must gather each requested variable from the owned particle ranges;
halo duplicates must not be emitted as additional particles. The gather publishes
the global owned-particle count for the host reader. `finishHostAccess()` restores
the subdomain-local counter and must be called before the next decomposed fan-out.
The CK VTP and restart writers use `VariablesWriteHelper` to bracket host access.
Reloading restart data scatters host particle state back to the subdomains.

`addSubdomainIDToWrite()` is an optional API intended to write the owner subdomain
for each particle. Its implementation currently registers the `SubdomainID`
variable, but the gather path does not populate that variable. Treat this output
feature as unimplemented and unverified until the owner IDs are assigned in the
gathered particle order and a test confirms the output.

## Current limitations and validation

- The cell-linked-list mesh is currently allocated for the whole domain on each
  subdomain. This avoids mesh remapping but scales mesh memory with the number of
  subdomains.
- Threaded host execution is available, but replica growth/reallocation and its
  synchronization need dedicated validation before relying on it for production.
- The SYCL decomposed path requires validation with a supported IntelLLVM/SYCL
  toolchain and device hardware. Host success alone does not establish device
  correctness.
- Contact interactions, including static walls and observer bodies, must be
  validated for the body's chosen decomposition and replication behavior; there
  is not a general documented replicated-body mode in the current API.
- Reduction association can change numerical results. Independently, the
  one-sided inner relation can append neighbors using atomic counters, so
  run-to-run ordering and trajectories may vary.
- `addSubdomainIDToWrite()` is not yet a verified output path, as noted above.

The repository has focused unit tests for slab geometry, rebalance, and host
fan-out semantics in
[`tests/unit_tests_src/for_2D_build/domain_decomposition/test_2d_domain_decomposition/test_2d_domain_decomposition.cpp`](../tests/unit_tests_src/for_2D_build/domain_decomposition/test_2d_domain_decomposition/test_2d_domain_decomposition.cpp).
When `SPHINXSYS_DECOMPOSITION` is enabled, the SYCL dambreak test configuration
adds two-subdomain and decomposed-restart CTest entries in
[`tests/tests_sycl/2d_examples/test_2d_dambreak_sycl/CMakeLists.txt`](../tests/tests_sycl/2d_examples/test_2d_dambreak_sycl/CMakeLists.txt).
The presence of those tests does not imply they have been run for the current
checkout or with every backend.

## Implementation map

- Policy selection: `src/shared/particle_dynamics/execution/execution_policy.h`
- Subdomain scope and runner: `src/shared/particle_dynamics/execution/`
- Slab geometry: `src/shared/domain_decomposition/domain_decomposition.{h,cpp}`
- Per-body decomposition: `src/shared/domain_decomposition/body_decomposition.h`
- Halo exchange and migration: `src/shared/domain_decomposition/subdomain_exchange.{h,hpp}`
- Loop integration: `src/shared/shared_ck/particle_dynamics/loop_range.h` and
  `particle_iterators_ck.h`
- Automatic interaction refresh: `src/shared/shared_ck/particle_dynamics/interaction_algorithms_ck.hpp`
- Host I/O and restart hooks: `src/shared/shared_ck/io_system/io_base_ck.hpp`
- Main-method setup: `src/shared/shared_ck/particle_dynamics/sph_solver.cpp`
