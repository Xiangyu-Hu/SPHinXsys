# SPHinXsys SYCL/GPU Paper — Benchmark Preparation Plan

## 1. Purpose of this document

This document defines the benchmark preparation and data-collection plan for the SPHinXsys SYCL/GPU paper.

The paper is primarily about the heterogeneous/SYCL implementation of SPHinXsys, rather than proposing a new SPH formulation. Therefore, the numerical experiments should be organized around three distinct questions:

1. **Numerical verification / backend consistency**
   - Does the SYCL/GPU implementation reproduce the expected physical/numerical solution?
   - Are CPU and GPU results numerically consistent?

2. **Computational scalability**
   - How does performance change with particle number?
   - How do CPU, RTX 2080 Super, RTX 3090, and RTX 4090 compare for the same controlled problem?

3. **Industrial-scale demonstrations**
   - Can the same heterogeneous framework execute genuinely complex multiphysics / industrial workloads?
   - What speedup and throughput are obtained in realistic applications?

The first two parts should be prepared in the **main SPHinXsys repository**.

The industrial applications currently live in **SPHinXsys-EinSIMO** and should be benchmarked separately later.

---

# 2. Repository and branch strategy

## 2.1 Main SPHinXsys repository

Create a dedicated branch for the GPU-paper benchmark work.

Suggested branch name:

```text
paper/sycl-gpu-benchmarks
```

or, if the project prefers feature-style naming:

```text
feature/gpu-paper-benchmarks
```

The purpose of this branch is to:

- collect all selected verification and scalability cases;
- standardize runtime statistics;
- standardize benchmark output files;
- allow CPU/GPU runs from the same benchmark definitions where possible;
- avoid modifying unrelated examples;
- freeze a reproducible code state for the paper.

Do **not** mix industrial EinSIMO applications into this branch.

---

## 2.2 SPHinXsys-EinSIMO repository

Industrial cases will be handled later in a separate benchmark branch.

Suggested future branch:

```text
paper/sycl-gpu-industrial-benchmarks
```

Potential industrial applications include:

- gearbox;
- nut–screw solid contact/friction;
- stirring;
- two-phase plate heat-exchanger filling / air displacement.

The current task is **not** to prepare or benchmark these industrial cases yet.

For now, only ensure that the data-collection format designed for the main repository can later be reused in the EinSIMO repository.

---

# 3. Recommended benchmark structure in the main repository

Create a common GPU-paper benchmark location, preferably without destroying the existing regression examples.

Suggested structure:

```text
tests/
└── gpu_paper_benchmarks/
    ├── README.md
    ├── common/
    │   ├── benchmark_recorder.h
    │   ├── benchmark_config.h
    │   └── benchmark_output/
    │
    ├── verification/
    │   ├── poiseuille_2d/
    │   ├── dambreak_3d/
    │   ├── fish_fsi_2d/
    │   ├── oscillating_beam_2d/
    │   ├── twisting_column_3d/
    │   └── diffusion_neumann_2d/
    │
    └── scalability/
        └── dambreak_3d/
```

Preferred implementation approach:

- reuse/call existing benchmark code as much as possible;
- avoid maintaining a second independent physical implementation;
- if possible, factor benchmark configuration into reusable headers/functions;
- preserve the original existing tests;
- the GPU-paper version may add runtime flags, resolution controls, statistics, and output suppression.

The goal is reproducibility, not creating duplicate physics implementations.

---

# 4. General benchmark philosophy

## 4.1 CPU versus GPU terminology

For this paper, do not make the distinction between:

- native host/TBB CPU;
- SYCL CPU backend;

a central scientific topic.

This is a computational-physics / multiphysics software paper, not a compiler-runtime comparison paper.

Use the simpler terminology:

```text
CPU execution
GPU execution
```

However, the exact backend should still be stored in the benchmark metadata.

If a case can run through the normal host execution path, use that as the preferred CPU baseline.

If a specific case is only available through a SYCL CPU device, that is acceptable, provided that the backend is explicitly recorded.

Do not spend substantial development effort creating a separate native CPU implementation solely for the paper if the existing SYCL implementation already runs correctly on CPU.

---

## 4.2 Same physical problem

For all CPU/GPU comparisons, keep the following identical:

- geometry;
- particle spacing;
- particle generation method;
- material properties;
- physical models;
- kernel formulation;
- smoothing length relation;
- CFL/time-step criteria;
- boundary conditions;
- initial conditions;
- end physical time;
- floating-point precision;
- code commit;
- output settings used for timing.

The paper should compare the **same physical problem**, not simply the same number of iterations.

---

## 4.3 Timing rule

Performance runs should use the same physical end time:

```text
t = 0 -> t_end
```

Do not artificially force the same number of time steps at different resolutions.

Because the CFL-controlled time step changes with resolution, record the actual number of steps.

For each run report both:

```text
Total compute wall time
Average compute time per step
```

The primary practical metric is total compute wall time required to reach the same physical end time.

---

## 4.4 Output / I/O

Performance timing should exclude expensive visualization output from the main compute time where possible.

Record separately:

```text
T_compute
T_IO
```

Use the same output policy for all compared hardware.

For pure performance runs, output frequency should be strongly reduced or disabled after initial validation.

---

# 5. Numerical verification cases

Recommended core verification suite:

| ID | Case | Dimension | Main physics | Main purpose |
|---|---|---:|---|---|
| V1 | Mixed Poiseuille flow | 2D | viscous internal flow | analytical accuracy + CPU/GPU consistency |
| V2 | Dam-break | 3D | free-surface fluid | nonlinear 3D fluid verification |
| V3 | Flow around deformable fish | 2D | deformable FSI | FSI / fluid-solid coupling |
| V4 | Oscillating beam | 2D | elastic solid | structural dynamics |
| V5 | Twisting column | 3D | nonlinear hyperelastic solid | large-deformation 3D solid |
| V6 | Diffusion with Neumann BC | 2D | diffusion / thermal analogue | diffusion operator verification |

Optional:

| ID | Case | Purpose |
|---|---|---|
| V7 | 3D Taylor bar | elastoplastic solid + basic contact |

Taylor bar should be prepared only if convenient. It is not required for the main six-case verification set.

---

# 6. Verification case details

## V1 — 2D Mixed Poiseuille Flow

### Scientific purpose

This should be one of the strongest quantitative verification cases because an analytical velocity profile is available.

Use it to demonstrate:

- analytical accuracy;
- CPU/GPU solution consistency;
- representative spatial convergence.

### Runs

Run at approximately three spatial resolutions.

Example labels:

```text
coarse
medium
fine
```

The exact `dp` values should be selected based on the existing test geometry.

Do not change physical parameters between resolutions.

### Required solution data

Record:

- analytical velocity profile;
- CPU velocity profile;
- GPU velocity profile;
- observer/sample coordinates;
- L2 error relative to analytical solution;
- Linf error relative to analytical solution;
- CPU–GPU L2 difference;
- CPU–GPU Linf difference.

If convenient, also retain pressure data.

### Required performance data

For every run also record the full benchmark statistics defined in Section 9, even if only a subset is eventually used in the paper.

### Expected paper figure

Likely:

```text
Velocity profile:
analytical vs CPU vs GPU
```

and possibly:

```text
error or convergence versus dp
```

---

## V2 — 3D Dam-Break

### Scientific purpose

This case has two roles:

1. numerical verification;
2. the dedicated scalability benchmark.

The numerical-verification subsection should use one established/reference resolution.

The scalability subsection should use multiple resolutions.

### Verification run

For one representative resolution record:

- CPU result;
- GPU result;
- pressure history at observer(s);
- mechanical energy history if available;
- representative free-surface configuration;
- CPU–GPU difference.

If an experimental/reference curve is already associated with the benchmark, retain it.

### Performance data

Record all timing and throughput fields in Section 9.

### Important

Do not generate separate visualization figures for every scalability resolution.

One representative dam-break geometry/result image is enough.

---

## V3 — 2D Flow Around a Deformable Fish

### Scientific purpose

Demonstrate that the heterogeneous framework supports actual deformable fluid–structure interaction rather than only fluid-only kernels.

### Resolution

One standard validated resolution is sufficient.

No full resolution-convergence study is required for this paper.

### Required solution data

Record at least:

- representative fish deformation / position;
- tip or characteristic displacement if available;
- mechanical energy if already used by regression;
- characteristic kinematics;
- CPU and GPU histories;
- CPU–GPU difference metrics.

### Performance data

Record all standard benchmark fields.

### Expected paper output

One compact quantitative comparison is enough.

Do not over-expand the FSI physics discussion because the paper is about the implementation framework.

---

## V4 — 2D Oscillating Beam

### Scientific purpose

Simple, clean structural-dynamics verification.

### Resolution

One standard validated resolution.

### Required solution data

Record:

- beam tip displacement versus time;
- CPU result;
- GPU result;
- reference/regression solution if available;
- CPU–GPU L2/Linf or time-series difference.

### Performance data

Record all standard benchmark fields.

---

## V5 — 3D Twisting Column

### Scientific purpose

Demonstrate:

- 3D structural dynamics;
- nonlinear/hyperelastic material response;
- large deformation;
- GPU execution of solid mechanics.

### Resolution

One standard validated resolution.

### Required solution data

Record:

- monitored tip position versus time;
- monitored tip velocity versus time;
- CPU result;
- GPU result;
- existing reference/regression solution;
- CPU–GPU differences.

### Performance data

Record all standard benchmark fields.

### Visualization

One representative deformed-shape snapshot may be retained.

---

## V6 — 2D Diffusion with Neumann Boundary Condition

### Scientific purpose

Demonstrate the diffusion/thermal class of operators.

This should be the second representative convergence case.

### Runs

Use approximately three resolutions.

### Required solution data

Record:

- scalar field;
- observer values;
- reference/analytical/regression solution;
- CPU values;
- GPU values;
- L2 error;
- Linf error;
- CPU–GPU difference;
- convergence trend.

### Performance data

Record all standard benchmark fields.

---

## V7 — Optional 3D Taylor Bar

Use only if preparation cost is small.

Purpose:

- elastoplastic material model;
- large deformation;
- basic contact.

Record:

- final/deformed shape;
- characteristic dimensions/displacement if available;
- CPU/GPU result consistency;
- full timing data.

This case is optional because oscillating beam + twisting column already provide a strong structural verification set.

---

# 7. Dedicated scalability benchmark

## Selected case

Use:

```text
3D dam-break
```

as the dedicated scalability benchmark.

This avoids introducing another physics case solely for performance scaling.

---

## 7.1 Resolution / particle-count levels

Use approximately 5–6 particle-number levels.

Suggested target range:

| Level | Approximate particle count |
|---|---:|
| S1 | 0.1–0.2 M |
| S2 | 0.4–0.5 M |
| S3 | ~1 M |
| S4 | ~2 M |
| S5 | ~4 M |
| S6 | ~6–10 M |

These are targets only.

The actual particle numbers should result naturally from the chosen `dp`.

The largest level should be determined by practical GPU memory and runtime constraints.

If a smaller GPU cannot run the largest problem, record:

```text
OOM
```

rather than reducing the maximum problem size solely to accommodate that device.

---

## 7.2 Hardware matrix

Target platforms:

```text
AMD 128-thread CPU
NVIDIA RTX 2080 Super
NVIDIA RTX 3090
NVIDIA RTX 4090
```

Run the same case and same physical configuration on every available platform.

---

## 7.3 Quantities to record

For each resolution and hardware platform:

- actual particle number;
- actual total neighbor interactions if available;
- physical end time;
- total number of time steps;
- advection steps if applicable;
- acoustic/substeps if applicable;
- total compute wall time;
- average time per time step;
- particle updates per second;
- neighbor interactions per second;
- GPIPS if available;
- peak memory usage;
- detailed kernel/component timings;
- output/I/O time;
- exact backend/device metadata.

---

## 7.4 Expected scalability plots

Prepare raw data sufficient for the following plots:

### Plot S1

```text
Compute wall time vs particle number
```

Curves:

- CPU;
- RTX 2080S;
- RTX 3090;
- RTX 4090.

### Plot S2

```text
Speedup relative to CPU vs particle number
```

Curves:

- RTX 2080S / CPU;
- RTX 3090 / CPU;
- RTX 4090 / CPU.

### Plot S3

```text
MPPS or MNIPS/GPIPS vs particle number
```

### Optional Plot S4

```text
Peak memory usage vs particle number
```

This may be useful for explaining maximum feasible problem sizes.

---

# 8. Industrial cases — future EinSIMO phase

Do not implement this section in the main repository now.

The current candidate industrial benchmarks are:

| Case | Main computational workload |
|---|---|
| Gearbox | complex free-surface flow + moving geometry + thermal coupling |
| Nut–screw | solid mechanics + moving contact + friction |
| Stirring | 3D moving geometry + strong fluid rearrangement/free surface |
| Plate heat exchanger filling | complex internal geometry + gas/liquid filling / multiphase flow |

The final paper will probably use 3–4 industrial cases.

Each industrial case normally needs only one realistic production resolution.

Industrial cases are **not** intended to perform systematic grid-convergence studies in this GPU paper.

The systematic particle-number scalability study is handled separately by the 3D dam-break benchmark.

For industrial cases later:

- record all available detailed statistics;
- preserve one representative physical-result snapshot;
- use only the most useful subset in the final paper.

---

# 9. Unified benchmark data to record for EVERY case

Even if the final paper does not use all of these fields, record the most complete dataset possible now.

The philosophy is:

> collect once, decide later what to publish.

---

## 9.1 Run identification

Record:

```text
case_name
run_id
date
git_repository
git_commit
git_branch
build_type
precision
backend
execution_policy
device_name
device_vendor
```

---

## 9.2 Hardware and software environment

Record:

```text
CPU model
CPU physical cores
CPU logical threads
system memory
GPU model
GPU memory
GPU driver
operating system
compiler
compiler version
SYCL implementation
SYCL backend
CMake configuration
relevant compile flags
SPHinXsys commit
```

If available:

```text
CUDA version
oneAPI version
Codeplay plugin version
```

---

## 9.3 Physical/numerical configuration

Record:

```text
dimension
particle spacing dp
smoothing length ratio
initial particle number
maximum/average active particle number if variable
physical end time
output interval
material/model identifier
time-step/CFL settings
```

Where relevant:

```text
number of solid bodies
number of fluid bodies
number of contact relations
number of observers
```

---

## 9.4 Step counts

Record as many as apply:

```text
total outer steps
advection steps
acoustic steps
solid substeps
contact updates
configuration/relation updates
diffusion steps
```

Do not assume CPU and GPU necessarily take exactly the same number of steps.

Record the actual values.

---

## 9.5 Overall timing

Record:

```text
total wall time
compute wall time excluding I/O
I/O wall time
initialization time
particle generation time if relevant
average compute time per outer step
average compute time per acoustic step
```

The final paper may use only compute time, but raw data should retain all categories.

---

## 9.6 Component/kernel timings

Instrument as many of the following as applicable:

```text
particle sorting
cell-linked-list construction
neighbor/relation/configuration update
density summation
advection-step setup
acoustic time-step calculation
pressure relaxation
density relaxation
viscous force
wall force/contact
solid stress integration
solid first-half integration
solid second-half integration
solid damping
contact search
contact factor
normal contact force
friction/tangential contact force
FSI viscous force
FSI pressure force
diffusion
thermal conduction
convection boundary
heat-source terms
reduction operations
force reduction
torque reduction
observer interpolation
data transfer if explicitly present
```

Not every case uses every field.

Unused fields may be blank/zero.

The objective is to preserve detailed profiling data even if only one representative case is later shown with a full component breakdown.

---

## 9.7 Workload / throughput metrics

Where possible calculate:

```text
particle_updates
particle_updates_per_second
MPPS
neighbor_interactions
neighbor_interactions_per_second
MNIPS
GPIPS
```

Suggested definitions should be kept consistent across hardware.

Do not change the metric definition between cases without documenting it.

---

## 9.8 Memory metrics

Where feasible record:

```text
peak host memory
peak GPU memory
initial GPU memory
maximum particle capacity
OOM / allocation failure
```

GPU memory is useful for interpreting why the RTX 2080S may fail for large cases while the RTX 3090/4090 can continue.

---

## 9.9 Numerical verification metrics

Where a reference solution is available, record:

```text
L2_CPU_reference
Linf_CPU_reference
L2_GPU_reference
Linf_GPU_reference
L2_GPU_CPU
Linf_GPU_CPU
```

For time-series data, also preserve:

```text
time
reference
CPU
GPU
```

For regression-based benchmarks, retain the existing regression metric in addition to these quantities.

Do not require bitwise-identical CPU/GPU results.

Floating-point reduction/order differences are expected.

The goal is numerical equivalence within the benchmark/discretization tolerance.

---

# 10. Recommended machine-readable output

Every benchmark run should generate a summary file automatically.

Recommended:

```text
benchmark_summary.csv
```

or one row appended to a common CSV database.

Suggested columns:

```text
case
run_id
git_commit
backend
device
precision
dp
particle_number
physical_time
outer_steps
advection_steps
acoustic_steps
solid_steps
compute_time
io_time
time_per_step
sorting_time
cll_time
configuration_time
advection_time
acoustic_time
pressure_time
density_time
viscous_time
solid_time
contact_time
fsi_time
diffusion_time
reduction_time
particle_updates
mpps
neighbor_interactions
mnips
gpips
peak_host_memory
peak_gpu_memory
L2_reference
Linf_reference
L2_cpu_gpu
Linf_cpu_gpu
status
```

Also create an environment file, for example:

```text
benchmark_environment.json
```

containing:

```json
{
  "repository": "...",
  "commit": "...",
  "branch": "...",
  "os": "...",
  "compiler": "...",
  "compiler_version": "...",
  "sycl_implementation": "...",
  "sycl_backend": "...",
  "cpu": "...",
  "cpu_threads": "...",
  "gpu": "...",
  "gpu_driver": "...",
  "gpu_memory": "...",
  "precision": "...",
  "cmake_flags": "..."
}
```

The exact schema may be improved by Cursor, but it should be shared by all benchmark cases.

---

# 11. Repeatability

For short/medium performance cases, run each configuration at least 3 times if practical.

Record:

```text
run_1
run_2
run_3
```

and calculate:

```text
mean
standard deviation
minimum/median if useful
```

For extremely expensive large cases, one warm-up plus 1–2 production runs may be acceptable.

Do not include compilation time.

If SYCL JIT/device initialization creates a substantial first-run penalty, identify it explicitly and ensure the timing methodology is consistent across devices.

---

# 12. Visualization policy

## 12.1 Verification cases

Verification cases should contain quantitative physical results.

Examples:

- Poiseuille: analytical/CPU/GPU velocity profile;
- dam-break: pressure/energy history and representative free surface;
- fish FSI: deformation/kinematics history;
- oscillating beam: tip displacement;
- twisting column: tip position/velocity;
- diffusion: reference/CPU/GPU scalar values.

These figures answer:

```text
Does the heterogeneous implementation compute the correct solution?
```

---

## 12.2 Scalability benchmark

The scalability subsection does not need physical-result figures for every resolution.

Use:

- one representative dam-break geometry/snapshot;
- performance curves for all resolutions.

The focus is computational scaling.

---

## 12.3 Industrial cases

Later, each industrial case should retain:

```text
1 geometry/setup view
+
1 representative real computed field/result
```

Examples:

- gearbox: oil/free-surface + velocity or temperature;
- nut–screw: deformation/stress/contact;
- stirring: velocity/free surface;
- heat exchanger: phase distribution/velocity.

These are not intended as new physical validation figures.

They demonstrate that the performance numbers correspond to genuine complex simulations rather than synthetic workloads.

A combined multi-panel figure may be sufficient for all industrial cases.

---

# 13. Proposed paper-results structure

## 4. Numerical verification of heterogeneous kernels

### 4.1 2D mixed Poiseuille flow
- analytical comparison;
- representative convergence;
- CPU/GPU consistency.

### 4.2 3D dam-break
- nonlinear free-surface verification;
- CPU/GPU consistency.

### 4.3 2D deformable FSI
- fish benchmark.

### 4.4 Structural dynamics
- 2D oscillating beam;
- 3D twisting column.

### 4.5 Diffusion
- Neumann-boundary diffusion;
- representative convergence.

Optional Taylor bar may go to supplementary material.

---

## 5. Computational performance

### 5.1 Benchmark methodology
- hardware;
- compiler;
- backend;
- precision;
- timing rules;
- reproducibility;
- workload metrics.

### 5.2 Particle-number scalability
- 3D dam-break;
- 5–6 resolutions;
- CPU + RTX 2080S + RTX 3090 + RTX 4090;
- wall time;
- speedup;
- MPPS/MNIPS/GPIPS;
- memory if useful.

### 5.3 Component-level profiling
Use one representative workload for detailed kernel/component breakdown.

The exact case can be selected after data collection.

Do not require every case to have a full profiling figure in the final paper, although all raw profiling data should be collected.

### 5.4 Industrial-scale applications
Later from SPHinXsys-EinSIMO:

- gearbox;
- nut–screw;
- stirring;
- optional heat exchanger.

---

# 14. Immediate Cursor task — main SPHinXsys repository only

The current implementation task should be limited to the following.

## Task A — create benchmark branch

Create:

```text
paper/sycl-gpu-benchmarks
```

from the current intended GPU-paper code base.

Do not modify `master/main` directly.

---

## Task B — inventory existing cases

Locate and inspect the existing implementations of:

```text
test_2d_mixed_poiseuille_flow_sycl
test_3d_dambreak_sycl
test_2d_flow_stream_around_fish_sycl
test_2d_oscillating_beam_sycl
test_3d_twisting_column_sycl
test_2d_diffusion_NeumannBC_sycl
```

Optional:

```text
test_3d_taylor_bar_sycl
```

For each case determine:

- current source location;
- current execution policies;
- whether CPU execution is already available from the same source;
- current regression/reference output;
- current timer fields;
- current observer outputs;
- current particle-number controls;
- how `dp` can be exposed as a runtime/build parameter;
- what changes are required for standardized benchmark output.

Do not change physics at this stage unless required to expose benchmark controls.

---

## Task C — implement shared benchmark recorder

Add a lightweight common recorder for:

- environment metadata;
- numerical setup;
- step counts;
- total/component timing;
- throughput;
- memory where practical;
- numerical error metrics.

Prefer a reusable implementation rather than duplicating timing code in each case.

---

## Task D — expose benchmark configuration

Each benchmark should support a clean way to choose at least:

```text
resolution / dp
physical end time
output enabled/disabled
output interval
benchmark mode
device/backend through existing build/runtime mechanism
```

Prefer command-line arguments or a compact benchmark configuration layer.

Avoid manually editing source code for every performance point.

---

## Task E — prepare verification modes

Prepare the six core cases so that one command can generate:

- physical/reference output;
- benchmark summary;
- metadata.

Poiseuille and diffusion must support at least three resolution levels.

Other verification cases need only one standard resolution initially.

---

## Task F — prepare 3D dam-break scalability sweep

Add at least 5–6 selectable resolution levels.

The levels should span approximately:

```text
0.1 M -> several million particles
```

with the final maximum chosen after memory testing.

The same executable/configuration should be usable on:

```text
CPU
RTX 2080S
RTX 3090
RTX 4090
```

without physics changes.

---

## Task G — create reproducible run scripts

Prepare scripts such as:

```text
scripts/
├── run_verification_cpu.sh
├── run_verification_gpu.sh
├── run_dambreak_scaling_cpu.sh
├── run_dambreak_scaling_gpu.sh
└── collect_benchmark_results.py
```

Exact naming may differ.

Scripts should:

- log the exact command;
- create a unique result folder;
- store metadata;
- store benchmark CSV/JSON;
- store stdout/stderr;
- never overwrite a previous run accidentally.

Example result structure:

```text
benchmark_results/
└── 3d_dambreak/
    └── 2026-xx-xx_commit-device-resolution/
        ├── summary.csv
        ├── environment.json
        ├── run.log
        └── solution_data/
```

---

# 15. What NOT to do yet

Do not yet:

- benchmark gearbox;
- modify nut–screw;
- modify stirring;
- modify plate heat exchanger;
- move EinSIMO cases into the public/main SPHinXsys repository;
- create a new 3D two-phase dam-break;
- port additional solid-contact functionality solely for this paper;
- perform exhaustive convergence studies for every validation case;
- run every industrial case at many resolutions;
- optimize individual kernels before baseline benchmark data are collected.

First establish a clean, reproducible benchmark infrastructure in the main SPHinXsys repository.

Industrial EinSIMO benchmarking will be a separate second phase.

---

# 16. Working principle for the paper

The benchmark program should ultimately support the following three claims:

### Claim 1 — Numerical consistency

Existing SPHinXsys formulations retain their expected numerical behavior after heterogeneous/SYCL execution.

### Claim 2 — Scalable acceleration

GPU acceleration increases computational throughput for large 3D particle systems, with systematic scaling demonstrated using the controlled 3D dam-break benchmark.

### Claim 3 — Multiphysics industrial applicability

The same framework can execute substantially different industrial workloads, including complex fluid, thermal, FSI, and solid-contact problems.

The benchmark code should therefore collect more information than the final paper will necessarily show.

**Collect first; select the clearest subset for publication later.**
