# Review of the merge with `xiangyu/hackathon_domain_decomposition`, and what it lets us drop

Date: 2026-09-10. Branch `shuang_and_ming_1` at b56d2c095, reviewed against my commit f7a1f80ca
("dambreak 2d test on cpu"). Line numbers are a snapshot of b56d2c095; re-check with `grep -n`.

## 0. Verdict in four lines

1. **Yes, the bring-up implementation can be simplified**, and noticeably: the fan-out moves out
   of every `exec()` into one `particle_for` / `particle_reduce` overload keyed on
   `LoopRangeCK<DecomposedExecution<P>, …>`, and the halo refresh inside the acoustic steps can be
   driven by the new `interact_variables_` instead of the `addPostInitialization()` hook. §3.
2. **HEAD does not compile, not even with both decomposition options OFF.** The rename
   `MultiHostExecution` → `DecomposedExecution` was only applied in `execution_policy.h`; ten other
   files still use the old names, and `particle_iterators_ck.h` carries dead overloads from the
   merge. §1 lists every site. This has to be fixed before anything else.
3. **The new policy shape has an overload trap on the SYCL side.** `MultiDevicePolicy` is now
   `DecomposedExecution<SYCLDevicePolicy>`; every place that pairs a catch-all template
   `foo(const ExecutionPolicy &)` with a non-template `foo(const SYCLDevicePolicy &)` picks the
   catch-all for it, i.e. the host branch. Separately, the SYCL `allocateComputingKernel` overload
   is not viable at all any more. Both are proven with a standalone test in §2.
4. Nothing in the exchange itself (halo plan, refresh, migration, scatter/gather, the B1–B4/B7
   fixes, `Force` evolving, `TotalLocalParticles`, relation over the local range) is affected by
   the merge; that all stays.

## 1. Hard blockers: the tree does not build

Checked with `clang++ -fsyntax-only` on `interaction_algorithms_ck.cpp` using the compile
command from `build/compile_commands.json` (options OFF). First errors:

```
sphinxsys_variable.h:179: error: no template named 'MultiHostExecution'
sphinxsys_variable.h:366: error: no template named 'MultiHostExecution'
sphinxsys_variable.h:431: error: no template named 'MultiHostExecution'
mesh_iterators.h:98:  error: no template named 'MultiHostExecution' in namespace 'SPH::execution'
subdomain_fan_out.h:75: error: no template named 'MultiHostExecution'
```

Complete list of stale names (from `grep -rn`), with the new name each should become:

| File | Lines | Old | New |
|---|---|---|---|
| `src/shared/common/sphinxsys_variable.h` | 179, 366, 431, 483 (comment) | `MultiHostExecution<PolicyType>` | `DecomposedExecution<PolicyType>` (but see §2 for the SYCL case) |
| `src/shared/common/algorithm_primitive.h` | 283, 314 | `MultiHostExecution<PolicyType>` | `DecomposedExecution<PolicyType>` |
| `src/shared/particle_dynamics/particle_iterators.h` | 64, 199 (+ comment 60) | `MultiHostExecution<PolicyType>` | `DecomposedExecution<PolicyType>` |
| `src/shared/meshes/mesh_iterators.h` | 98, 105 | `execution::MultiHostExecution<PolicyType>` | `execution::DecomposedExecution<PolicyType>` |
| `src/shared/particle_dynamics/execution/subdomain_fan_out.h` | 75, 89 | `MultiHostExecution<PolicyType>` | `DecomposedExecution<PolicyType>` (host only, see §2.3) |
| `src/shared/shared_ck/particle_dynamics/particle_iterators_ck.h` | 100–147 | whole block: calls `sequenced_particle_for` etc. that no longer exist, and names `ParallelMultiHostPolicy` | delete the block; it is the pre-merge unary-function API, superseded by the `(loop_range, implementation, dt)` API above it |
| `src/shared/shared_ck/particle_dynamics/configuration_dynamics/base_configuration_dynamics.h` | 77, 83, 89, 95 | `ParallelMultiHostPolicy` / `SequencedMultiHostPolicy` | `MultiHostPolicy` / `SequencedMultiHostPolicy` (the second name still exists) |
| `src/shared/domain_decomposition/subdomain_exchange.hpp` | 779–780 | `ParallelMultiHostPolicy` | `MultiHostPolicy` |
| `src/src_sycl/shared/particle_dynamics/particle_iterators_multi_device_sycl.h` | 65 | `LoopRangeCK<ParallelMultiDevicePolicy, …>` | goes away entirely under §3.1 |
| `src/src_sycl/shared/particle_dynamics/configuration_dynamics/base_configuration_dynamics_sycl.h` | 51 | `ParallelMultiDevicePolicy` | `MultiDevicePolicy` |
| `tests/unit_tests_src/…/test_2d_domain_decomposition.cpp` | 171, 185, 188, 201, 208 | `par_multi_host` | `multi_host` |
| `src/shared/particle_dynamics/execution/subdomain_scope.h` | 35–36 (comment) | old names | cosmetic |

And one macro mismatch, which makes the ON build silently behave like OFF even after the
renames: `execution_policy.h:85` now tests `SPHINXSYS_MULTI_HOST`, while `CMakeLists.txt:19,74`,
`sphinxsys.h:63` and the case file `dambreak_decomp.cpp` (17 places) still use
`SPHINXSYS_MULTI_SUBDOMAIN_HOST`. Pick one name (I would keep the CMake option as it is and
change the one line in `execution_policy.h`, since the case and the docs already use it).

## 2. The overload trap of `DecomposedExecution<SYCLDevicePolicy>`

### 2.1 What happens

Before the merge, `MultiDeviceExecution<P>` derived from `DeviceExecution<P>` and the SYCL
overloads were templates on `DeviceExecution<PolicyType>`. After the merge the SYCL overloads are
plain non-template functions on `const SYCLDevicePolicy &` (e.g. `sphinxsys_variable.h:172`,
`sphinxsys_constant.h:86-91`, `implementation.h:78-95`), next to catch-all templates
`template <class ExecutionPolicy> foo(const ExecutionPolicy &)`.

For an argument of type `MultiDevicePolicy = DecomposedExecution<SYCLDevicePolicy>`:

- the catch-all template deduces `ExecutionPolicy = MultiDevicePolicy`: **exact match**;
- the SYCL overload needs a derived-to-base conversion: **conversion rank**.

Overload resolution ranks the implicit conversion sequence before the template/non-template
tie-break, so the catch-all wins. That is the same rule that forced the forwarding overloads in
§2.1 of the previous commit summary, only now it hits the device side.

Standalone proof (`scratchpad/ovl.cpp`, clang++ -std=c++17), reproducing both patterns exactly:

```
SYCLDevicePolicy  -> DelegatedData: SYCL device | allocate: malloc on host
MultiDevicePolicy -> DelegatedData: host catch-all | allocate: malloc on host
```

Two findings in that output:

1. **`MultiDevicePolicy` resolves to the host branch** of every such pair. Affected, from the
   grep: `DiscreteVariable::DelegatedData/reallocateData/prepareForOutput/finalizeLoadIn`,
   `SingleVariable::DelegatedData`, `ConstantArray::DelegatedData`,
   `ComputingKernelArray::DelegatedData`, the variable-array equivalents, `particle_for` /
   `particle_reduce` on an `IndexRange`, `exclusive_scan`, `generic_for`, `package_for`,
   `mesh_for`.
2. **Even plain `SYCLDevicePolicy` allocates its computing kernel with `malloc`**. The SYCL
   overloads in `implementation.h:78, 83, 91` kept a `template <class ComputingKernelType, class
   PolicyType>` header (from 948f6a525) while their parameter became `const SYCLDevicePolicy &`,
   so `PolicyType` can no longer be deduced and the overload is never viable. The kernel is then
   host memory dereferenced inside the device lambda. This is independent of the decomposition
   and worth fixing first: drop `class PolicyType` from those three templates.

Neither can be compiled here (no SYCL toolchain), which is why the standalone test.

### 2.2 The rule that fixes it, once

`DecomposedExecution<P>` should mean "P, plus a fan-out" and nothing else. So for every operation
that has a `SYCLDevicePolicy` overload, the decomposed policy must land on the **same** overload as
its base. Two ways, pick one and apply it uniformly:

| Option | Shape | Pros / cons |
|---|---|---|
| **A (recommended)**: forwarders on the decomposed policy | `template <class P> auto foo(const DecomposedExecution<P> &) { return foo(P{}); }` (what `particle_iterators.h:64` already does for the host) | one generic line per operation, serves host and device alike; needs the host-replica cases (`DelegatedOnHostSubdomain`) to be expressed as `foo(ParallelPolicy)`-with-subdomain rather than as a separate decomposed overload |
| B: explicit `MultiDevicePolicy` overloads | `DataType *DelegatedData(const MultiDevicePolicy &) { return DelegatedOnDevice(); }` next to each SYCL overload | non-template, exact match, always wins; but one extra overload per operation and it is easy to forget one (the failure is silent: the host branch runs) |

For `DiscreteVariable::DelegatedData` specifically, the host replica and the device replica are
different functions (`DelegatedOnHostSubdomain()` vs `DelegatedOnDevice()`), so this one needs a
`MultiDevicePolicy` overload under either option; the generic `DecomposedExecution<P>` template
then covers only the host policies. Same for `reallocateData` and `SingleVariable::DelegatedData`.

### 2.3 `fanOutOverSubdomains` has the mirror-image problem

`subdomain_fan_out.h:75` (host runner) will become a template on `DecomposedExecution<P>` and
`device_environment_sycl.h:187` is a non-template on `MultiDevicePolicy`. Here the non-template
is the exact match and wins, so this direction is fine as long as the device header is visible at
the point of instantiation. To make it robust rather than lucky, constrain the host template:
`std::enable_if_t<!std::is_base_of_v<SYCLDevicePolicy, P>>`.

## 3. What the merge lets us simplify

### 3.1 The fan-out moves into the loop (your intention, and it works now)

Why it did not work before: the kernel pointer was fetched *before* `particle_for` and captured
into the lambda, so a fan-out inside the loop would have handed subdomain 0's kernel to every
subdomain (`docs/domain_decomposition.md` §2, "Where the fan-out lives"). The new API
`particle_for(loop_range, implementation, dt)` calls `implementation.getComputingKernel()` inside
the loop function (`particle_iterators_ck.h:18, 27`), and `getComputingKernel()` resolves on
`currentSubdomainID()` (`implementation.h:105`). So a per-subdomain call of the *base policy's*
`particle_for` from inside the fan-out gets the right kernel for free.

The one remaining obstacle is the **loop range**: `LoopRangeCK<P, SPHBody>` calls
`DelegatedData(P{})` in its constructor (`loop_range.h:24`), and the algorithms construct it on the
host thread before calling `particle_for`. Under the old design this was solved by constructing it
inside the fan-out; under the new one the range object must *defer* that until the fan-out has
bound a subdomain. Target shape (all in `loop_range.h` / `particle_iterators_ck.h`; replaces
`particle_iterators_multi_device_sycl.h`):

```cpp
// A decomposed range only remembers the identifier; the per-subdomain range is built inside
// the fan-out, where DelegatedData() resolves to that subdomain's replica.
template <class PolicyType, class Identifier>
class LoopRangeCK<DecomposedExecution<PolicyType>, Identifier>
{
  public:
    explicit LoopRangeCK(Identifier &identifier) : identifier_(identifier) {};
    LoopRangeCK<PolicyType, Identifier> onCurrentSubdomain() const
    {
        return LoopRangeCK<PolicyType, Identifier>(identifier_);
    };
  private:
    Identifier &identifier_;
};

template <class PolicyType, class Identifier, class KernelImplementationType>
void particle_for(const LoopRangeCK<DecomposedExecution<PolicyType>, Identifier> &loop_range,
                  KernelImplementationType &implementation, Real dt)
{
    fanOutOverSubdomains(DecomposedExecution<PolicyType>{}, [&]()
                         { particle_for(loop_range.onCurrentSubdomain(), implementation, dt); });
}

template <typename Operation, class PolicyType, class Identifier, class ReturnType,
          class KernelImplementationType>
ReturnType particle_reduce(const LoopRangeCK<DecomposedExecution<PolicyType>, Identifier> &loop_range,
                           ReturnType temp, KernelImplementationType &implementation, Real dt)
{
    return reduceOverSubdomains<Operation>(
        DecomposedExecution<PolicyType>{}, temp, [&]()
        { return particle_reduce<Operation>(loop_range.onCurrentSubdomain(), temp, implementation, dt); });
}
```

Why this is enough: `PolicyType` is `ParallelPolicy`, `SequencedPolicy` or `SYCLDevicePolicy`, so
the inner call lands on the existing TBB or SYCL overload unchanged; the fan-out function is chosen
by the runner (host) or the device environment (SYCL) as today. The `LoopRangeCK<SPHBody>`
constructor taking a `SingleVariable<UnsignedInt> *` (`loop_range.h:26`) needs the same deferred
treatment if any caller uses it with a decomposed policy (none does in the CK algorithms today).

What this removes from the previous commit:

| Removed | Where | Why it is no longer needed |
|---|---|---|
| `fanOutOverSubdomains` around `runAllSteps` | `interaction_algorithms_ck.hpp:154, 193` | each step's `particle_for` fans out on its own |
| the three-fan-out `exec()` of `OneLevel` | `interaction_algorithms_ck.hpp:255-266` | every `particle_for` is already a barrier over the subdomains, so a refresh can sit between any two steps; `exec()` becomes `setUpdated; setupDynamics; runAllSteps` for all three variants, identical to the non-decomposed code |
| `addPostInitialization()` and `post_initialization_` | `interaction_algorithms_ck.h:59-66` | replaced by §3.2 |
| `LoopRangeCK<ParallelMultiDevicePolicy, …>` | `particle_iterators_multi_device_sycl.h` | replaced by the generic specialization above, which serves both backends |
| (already gone in the merge) fan-outs in `StateDynamics` / `ReduceDynamicsCK` | `simple_algorithms_ck.h` | now handled by the overload above; note that at HEAD they are simply *missing* for the decomposed policies, which is why the ON build would be wrong even after the renames |

What stays as explicit fan-outs, on purpose: `SubdomainExchange` (the pack/barrier/pull protocol
*is* the fan-out structure, `subdomain_exchange.hpp:407, 565-583, 596, 712, 727`) and the
configuration dynamics (`update_cell_linked_list.hpp:81`, `update_body_relation.hpp:112, 242`,
`particle_sort_ck.hpp:81`), whose loops run on an `IndexRange` with a captured kernel inside their
`…OnCurrentDevice()` bodies. Those go through the forwarding overloads of `particle_iterators.h`,
which do not fan out, so there is no double fan-out. The `insideFanOut()` guard remains necessary
for exactly that reason.

### 3.2 Halo refresh driven by `interact_variables_`

What the merge added: `Interaction<Inner<…>>::interact_variables_` /
`addInteractVariable()` (`interaction_ck.h:88-91`) and the contact twin; the acoustic steps register
`Pressure` (1st half, `acoustic_step_1st_half.hpp:59`) and `Velocity` (2nd half). That is precisely
the set the previous commit refreshed by hand through `sync_pressure` / `sync_velocity` +
`addPostInitialization()` (commit summary §2.3), so the two designs agree on *what*; the question
is *who triggers it*.

Proposed shape: the interaction algorithm refreshes its interact variables right before its
interaction step, through a hook that is a no-op outside a decomposed run.

```cpp
// base_particles.h: a tiny interface, so the shared code never sees SubdomainExchange.
class HaloRefresher { public: virtual void refreshHalo(DiscreteVariables &variables) = 0; };
// BaseParticles gets: HaloRefresher *halo_refresher_ = nullptr; setHaloRefresher(); refreshHalo(vars)
//   { if (halo_refresher_) halo_refresher_->refreshHalo(vars); }

// subdomain_exchange.hpp constructor: particles.setHaloRefresher(this);  (SubdomainExchange
//   implements HaloRefresher by calling its refreshHalo(variables))

// interaction_algorithms_ck.hpp, both runInteraction():
//   this->particles_->refreshHalo(this->interact_variables_);   // Inner
//   this->contact_particles_->refreshHalo(this->contact_interact_variables_);   // Contact
//   particle_for(...);
```

Why here and not in `runAllSteps`: `runInteraction` is the only place that knows the variable set,
and it sits after the initialization step's `particle_for` (a barrier) and before the interaction's
own `particle_for`. For a `WithUpdate` algorithm it is simply the first thing that happens. The
contact call is a no-op for the wall (never decomposed, hence no refresher), and is already right
for a future decomposed contact body.

What this removes from the case file and the library: the four `SyncHaloStateCK` objects and the
two `addPostInitialization` lines in `dambreak_decomp.cpp:172-184`, and `addPostInitialization`
itself. `SyncHaloStateCK` can stay as a manual escape hatch or go.

**The point that needs a decision (D2 below): the two advection-step refreshes.** The interaction
kernel of the first half reads at the neighbor index `Vol_[j]`, `p_[j]`, `correction_(j)` and the
position (`acoustic_step_1st_half.hpp:104-114`); the merge registers only `Pressure`. That matches
the hand-placed sets because `VolumetricMeasure` and `LinearCorrectionMatrix` change once per
advection step and were refreshed there (`sync_volume`, `sync_correction`,
`dambreak_decomp.cpp:255-259, 366-370`). Two consistent readings of `interact_variables_`:

- **"what changed since the last interaction" (current listing)**: minimal traffic; the advection
  step refreshes stay in the case file as they are; the kernel author must know the time loop.
- **"everything read at `j`"**: `Vol`, `B`, `pos`, `p` for the first half; the case file loses all
  manual refresh points; but `Vol` and `B` are then re-sent every acoustic step (~5-10× the
  necessary traffic for those two), unless a per-variable "written since last refresh" mark is
  added later.

I recommend the first for now, since it keeps the measured configuration, and to revisit when the
device path exists and traffic can be measured.

### 3.3 Things the merge does not change, so they stay

- Migration by tail filling, the halo plan, `scatterFromHost`/`gatherToHost`/`finishHostAccess`
  and the `setValue(0, total)` fix (B2), the replicas pre-touch (B1), the relation over
  `TotalLocalParticles` (B3), the order migrate → sort → UpdateHalo → cell list (B4), the
  `Force` evolving fix, `replicaCreationMutex()`. None of them depend on where the fan-out lives.
- The forwarding overloads for `IndexRange` loops (`particle_iterators.h:64, 199`,
  `algorithm_primitive.h:283, 314`, `mesh_iterators.h:98, 105`): still needed, only renamed, and
  under §2.2 option A they become the model for the whole rule.
- The limitations in `docs/domain_decomposition.md` §10 (threaded runner reallocation, replicated
  bodies, per-subdomain mesh).

## 4. Suggested order

| Step | Content | Gate |
|---|---|---|
| 0 | Renames of §1, the macro name, delete the dead block in `particle_iterators_ck.h`, drop `class PolicyType` in `implementation.h:78-95` | OFF build compiles; `test_2d_dambreak_sycl` regression passes |
| 1 | `LoopRangeCK<DecomposedExecution<P>>` + the two overloads of §3.1; delete the exec-level fan-outs and `particle_iterators_multi_device_sycl.h` | ON build, 1 subdomain, energy identical to OFF bit for bit (the baseline output of f7a1f80ca is in `~/vcpkg/sphinxsys/build/…/test_2d_dambreak_decomp/bin/baseline`) |
| 2 | `HaloRefresher` hook + refresh in `runInteraction`; remove `addPostInitialization` and the acoustic `SyncHaloStateCK` objects from the case | ON build, 2 subdomains, energy identical to OFF (commit summary §3 levels 1–2) |
| 3 | The `MultiDevicePolicy` overload audit of §2.2 on the SYCL headers | cannot be tested on this machine; review only |
| 4 | Decide D2 and update `docs/domain_decomposition.md` §2 ("Where the fan-out lives") and §8 accordingly | docs match code |

Steps 1 and 2 are independent of each other and both depend on 0.

## 5. Decisions to take

| Id | Question | Options | Recommendation |
|---|---|---|---|
| D1 | How to make `DecomposedExecution<SYCLDevicePolicy>` hit the SYCL overloads | A: generic forwarders `DecomposedExecution<P>` → `P`; B: explicit `MultiDevicePolicy` overloads | A, with B for the three variable-replica functions that genuinely differ between host and device |
| D2 | Meaning of `interact_variables_` | "changed since the last interaction" (keep `sync_volume`/`sync_correction` in the case); "everything read at `j`" (no manual refresh, more traffic) | the first, for now |
| D3 | Keep `SyncHaloStateCK` after step 2 | keep as manual escape hatch; delete | keep until D2 is settled |
| D4 | Macro name | `SPHINXSYS_MULTI_SUBDOMAIN_HOST` (CMake, case, docs); `SPHINXSYS_MULTI_HOST` (`execution_policy.h`) | keep the CMake name, fix the one line |

## 6. Done on 2026-09-10 (working tree, not committed)

Decisions taken: D1 = A with explicit `MultiDevicePolicy` overloads for the three replica
functions; D2 = "changed since the last interaction"; D3 = `SyncHaloStateCK` kept; D4 = CMake
name kept, `execution_policy.h` fixed.

| Step | What changed | Gate |
|---|---|---|
| 0 | renames of §1; macro fixed; dead block of `particle_iterators_ck.h` removed; `class PolicyType` dropped from the three SYCL kernel-allocation templates; forwarders `DecomposedExecution<P>` → `P` in `implementation.h`, `subdomain_fan_out.h` (`copyBetweenSubdomains`), `sphinxsys_constant.h`, `sphinxsys_variable_array.h`, `sphinxsys_variable.h` (`prepareForOutput`, `finalizeLoadIn`); explicit `MultiDevicePolicy` overloads of `DelegatedData` (both variable kinds) and `reallocateData`; host fan-out templates constrained to non-SYCL bases | OFF build compiles; fan-out unit test compiles |
| 1 | `DeferredLoopRangeCK` + three `LoopRangeCK<DecomposedExecution<P>, Identifier>` specializations in `loop_range.h` (one per identifier, since a single partial specialization is ambiguous with the per-identifier ones); loop-level `particle_for` / `particle_reduce` in `particle_iterators_ck.h`; all exec-level fan-outs and `addPostInitialization` removed from `interaction_algorithms_ck.{h,hpp}`; `particle_iterators_multi_device_sycl.h` deleted | ON, 1 subdomain: energy and 20 vtp frames identical to OFF bit for bit |
| 2 | `HaloRefresher` interface and hook in `base_particles.h`; `SubdomainExchange` implements it and installs itself; both `runInteraction()` refresh their interact variables first; the acoustic `SyncHaloStateCK` objects and hook calls removed from the case | ON, 2 subdomains, 1 thread: energy identical to OFF, vtp positions ≤ 1e-9 (Float32 print); 8 threads: energy identical; owned counts 1600/1600 → 1347/1853, `checkConsistency()` clean |
| 4 | `docs/domain_decomposition.md` §1, §2 ("Where the fan-out lives"), §5, §8, §9 updated | — |

Not done / not testable here: the SYCL side was only reviewed and syntax-checked through the
shared headers (no SYCL toolchain on this machine); the threaded runner limitation of
`docs/domain_decomposition.md` §10 is unchanged; no ASan run this time (the exchange code
did not change).
