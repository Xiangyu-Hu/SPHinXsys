# CURSOR_APHI_STAGE10_14C_CLEANUP_AND_NEXT_PLAN.md

> Branch: `feature/electromagnetic`  
> Scope: SPHinXsys SYCL CK A-phi electromagnetic heating solver  
> Stage: 10.14-C cleanup and robustness before standard TEAM7 alignment  
> Purpose: close the remaining structural, semantic, and test-gate issues after Stage 10.14-B; keep the code aligned with SPHinXsys SYCL CK style; prepare for standard TEAM7 single-resolution validation.

---

## 0. Executive summary

After reviewing the missing helper files:

- `aphi_source_driven_em_solve_helpers.h`
- `aphi_team7_contact_test_helpers.h`
- `aphi_physical_region_audit_helpers.h`
- `aphi_em_geometry_relaxation_diagnostic_helpers.h`

the current implementation is substantially better than the previous review suggested.

### 0.1 Confirmed progress

The following are now confirmed:

1. **Source-driven A-phi GMRES pipeline exists and is CK-style.**
   - It assigns TEAM7-like materials.
   - It zeros the coil/source region conductivity for source-only coil.
   - It assigns impressed current RHS in the coil/source region.
   - It solves A-phi using `AphiMatrixFreeSolveCK`.
   - It computes B/H/E/J/Joule and writes observables.

2. **Three-body Contact case is a true three-body geometry scaffold.**
   - `AirBody`: non-overlapping air slabs.
   - `CoilBody`: coil/source box.
   - `PlateBody`: conductor/plate box.
   - Contact relations:
     - air contacts coil + plate;
     - coil contacts air;
     - plate contacts air.
   - Source RHS is assigned only on the coil body in the three-body coupled GMRES path.

3. **Physical-region audit helper exists.**
   - It reports bodywise particle counts.
   - It separates physical box and passive shell.
   - It separates air/conductor/source particles.
   - It reports sigma/nu ranges, RHS norm, and Joule integrals by region.

4. **SYCL CK coding style is mostly respected.**
   - LocalDynamics + UpdateKernel + DelegatedData patterns are used.
   - StateDynamics / InteractionDynamicsCK style is retained.
   - Host-side diagnostics use explicit device-to-host sync helpers.
   - The source-driven helper remains under `tests/extra_source_and_tests`, which is correct for this prototype stage.

### 0.2 Remaining concerns

The following still require cleanup before moving to standard TEAM7:

1. **The two-body Contact test semantics are still potentially ambiguous.**
   - The true three-body case is clean.
   - The two-body Contact source-driven helper still needs explicit bodywise region audit to avoid misleading labels such as `air_joule_integral` when the body may contain non-air regions.

2. **`annular source` is not yet proven to be a real annular geometry.**
   - If the current annular helper still uses TEAM7 box layout, rename it to `team7_coil_source_region_diagnostic`.
   - If the test name stays `annular`, implement a real annular/racetrack mask.

3. **`geometry relaxation diagnostic` is not true relaxation.**
   - Current helper creates a three-body lattice scaffold and reports body particle counts.
   - It does not run particle relaxation dynamics.
   - It currently names particle counts as neighbor counts, which is misleading and must be fixed.

4. **Joule-to-thermal source-driven closure needs a low-power/high-power resolution-floor test matrix.**
   - The current resolution-floor fallback is acceptable for low-temperature-increment observability.
   - But a high-power branch must prove that the SYCL thermal update actually changes temperature without fallback.

5. **Boundary policy gate should not rely on row order.**
   - Use explicit policy tags: Baseline / EnlargedAir / PassiveAirShell.
   - Do not assume rows `[0,1,2,3]`.

---

## 1. Frozen assumptions for this cleanup round

Do not change these unless explicitly instructed:

```text
lambda_A production = off
Contact A-penalty = research only
Contact A-penalty stencil = InnerOnly when enabled
projection = deferred
ghost/mirror EM boundary = not implemented
passive_air_shell = diagnostic only
default EM outer-domain policy = enlarged air-domain padding
thermal boundary conditions = not a focus in this stage
standard TEAM7 reference comparison = next stage, not this cleanup round
multiresolution = deferred
```

---

## 2. Code structure requirements

All new or modified code must continue to follow SPHinXsys SYCL CK conventions.

### 2.1 Required CK style

Use patterns like:

```cpp
class SomeDynamicsCK : public LocalDynamics
{
  public:
    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : variable_(encloser.dv_variable_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            // device-side update
        }
    };

  protected:
    DiscreteVariable<Real> *dv_variable_;
};
```

Allowed dynamics wrappers:

```cpp
StateDynamics<MainExecutionPolicy, SomeLocalDynamicsCK>
InteractionDynamicsCK<MainExecutionPolicy, SomeInteractionDynamicsCK<Inner<> / Contact<>>>
AphiMatrixFreeSolveCK<MainExecutionPolicy>
AphiMultiBodyContactGMRESSolverCK<MainExecutionPolicy>
```

### 2.2 Host diagnostics

Host diagnostics are allowed only if they are explicitly diagnostic/test helpers. Before host access, always call the appropriate sync helper:

```cpp
syncVariableToHost<Real>(particles, variable_name);
syncVariableToHost<Vecd>(particles, variable_name);
syncAphiBlockToHost(particles, block_names);
```

If host values are modified and then used by device dynamics, sync back:

```cpp
syncVariableToDevice<Real>(particles, variable_name);
syncVariableToDevice<Vecd>(particles, variable_name);
```

### 2.3 Do not move prototype diagnostics into production solver yet

Keep the following under test/diagnostic helpers:

```text
aphi_physical_region_audit_helpers.h
aphi_em_geometry_relaxation_diagnostic_helpers.h
aphi_annular_source_region_diagnostic_helpers.h
boundary support policy diagnostics
thermal resolution-floor diagnostics
```

Do not merge them into core production headers until the tests stabilize.

---

## 3. Current file-level review

### 3.1 `aphi_source_driven_em_solve_helpers.h`

#### What is good

The helper correctly implements a source-driven A-phi pipeline:

```text
Initialize A-phi variables
Assign TEAM7-like material regions
Zero sigma in coil/source region for source-only coil
Optionally zero sigma outside physical box for passive shell
Assign impressed current RHS in coil/source region
Solve A-phi by GMRES
Compute corrected B = curl A
Compute H = nu B
Compute grad(phi)
Compute E
Compute J
Compute Joule
Write VTP / probe CSV
Collect physical-region metrics
```

The following implementation choices are good and should be preserved:

```text
coil/source region is source-only: sigma = 0
source RHS still exists in coil region
passive_air_shell is sigma=0 outside physical box
A-divergence penalty remains off
phi gauge penalty remains on
B is computed with corrected gradient/curl diagnostic
VTP includes A/B/H/E/J/Joule and material fields
physical box masks are used for region metrics
```

#### Required cleanup

1. `dummy_shell_width_scale` is still present as a deprecated alias. This is acceptable temporarily, but all new tests and records must use:

```text
passive_air_shell_width_scale
```

2. `sourceDrivenEmSolvePassed()` currently uses:

```cpp
metrics.max_abs_J > spec.min_conductor_J
```

but `metrics.max_abs_J` is the global max over all particles, not conductor-only.

Required change:

Add:

```cpp
Real max_abs_J_conductor = 0.0;
Real max_Joule_conductor = 0.0;
```

to `AphiSourceDrivenEmSolveMetrics`, compute them with the conductor mask, and update the pass gate:

```cpp
metrics.max_abs_J_conductor > spec.min_conductor_J
metrics.max_Joule_conductor > 0
```

Keep global `max_abs_J` only for diagnostic output.

3. If `source_rhs_l2` is computed over all particles, also add:

```cpp
Real source_rhs_l2_source_region = 0.0;
```

This avoids false positives if RHS appears outside the intended source region.

---

### 3.2 `aphi_team7_contact_test_helpers.h`

#### What is good

The three-body Contact scaffold is structurally correct:

```text
AirBody: air slabs
CoilBody: coil/source box
PlateBody: conductor box
AirBody contacts CoilBody and PlateBody
CoilBody contacts AirBody
PlateBody contacts AirBody
```

The source-driven three-body GMRES path is also reasonable:

```text
initialize air/coil/plate material fields
source_only_coil option sets coil sigma to zero
zero RHS on air and plate
assign impressed current RHS on coil
zero solutions on all bodies
update relations
run multi-body contact GMRES
compute bodywise true residual
run B/E/J/Joule postprocess
collect plate observables
```

This is a much cleaner physical scaffold than a single-body region-tagged test.

#### Required cleanup

1. Confirm that `source_only_coil = true` is used for TEAM7-style source tests unless intentionally testing physical coil self-Joule.

2. Add or confirm hard gates for three-body case:

```text
air_body particle count > 0
coil_body particle count > 0
plate_body particle count > 0
coil RHS norm > 0
plate conductor Joule > 0
air Joule / plate conductor Joule < 0.05
coil Joule should be ~0 when source_only_coil = true
all body residuals finite
GMRES converged
interface spike hard gate passed
```

3. Contact relations are currently air-mediated:

```text
coil <-> air
plate <-> air
air <-> coil, plate
```

This is acceptable for TEAM7-like coil/air/plate geometry if coil and plate are not directly adjacent. Add a comment:

```cpp
// Coil and plate are not directly connected by Contact<>.
// Coupling is mediated through AirBody. This is intentional for the current TEAM7-like scaffold.
```

If future geometry has coil and plate within the kernel support or directly adjacent, add a direct coil-plate Contact relation only after a dedicated test.

---

### 3.3 `aphi_physical_region_audit_helpers.h`

#### What is good

The audit helper is useful and should be retained.

It already reports:

```text
body_name
particle_count_total
particle_count_physical_box
particle_count_shell
particle_count_air
particle_count_conductor
particle_count_source
sigma_min/max
nu_min/max
rhs_l2
joule_total
joule_air
joule_conductor
joule_source
max_J_conductor
max_Joule_conductor
```

This is exactly the right direction.

#### Required cleanup

1. `metrics.joule_integral_total` currently uses all particles:

```cpp
[](const Vecd &) { return true; }
```

For passive shell cases, add:

```cpp
joule_integral_physical_box
joule_integral_shell
```

and avoid using total-all-particles for physical claims.

2. `max_J_source` is computed using:

```cpp
layout.coil_material.sigma * E_source
```

But for source-only coil, the actual sigma is zero. This diagnostic can be misleading.

Replace or add:

```cpp
max_E_source
max_J_source_from_actual_sigma
max_impressed_rhs_source
```

Do not use material-layout coil sigma as if the source-only coil were conductive.

3. `AphiThreeBodyContactAuditSummary::global_conductor_joule` currently uses only:

```cpp
summary.plate.joule_integral_conductor
```

This is acceptable for the current source-only coil case, but add:

```cpp
global_source_joule = summary.coil.joule_integral_source
global_shell_joule = air.shell + coil.shell + plate.shell
```

4. Add explicit pass helper:

```cpp
bool threeBodyContactAuditPassed(const AphiThreeBodyContactAuditSummary &summary, Real max_air_to_conductor_joule_ratio);
```

Expected gates:

```text
summary.air.particle_count_air > 0
summary.coil.particle_count_source > 0
summary.plate.particle_count_conductor > 0
summary.total_source_rhs_l2 > 0
summary.global_conductor_joule > 0
summary.global_air_joule / (summary.global_conductor_joule + TinyReal) < max_ratio
source-only coil: global_source_joule <= max(1e-12, 1e-6 * global_conductor_joule)
```

---

### 3.4 `aphi_em_geometry_relaxation_diagnostic_helpers.h`

#### What is good

This helper creates a lightweight three-body lattice scaffold:

```text
AirBody
CoilBody
PlateBody
Inner<> for each
Contact<> relations
lattice particle generation
```

This is useful as a geometry-scaffold smoke test.

#### Problem

It is not particle relaxation.

Current metrics:

```cpp
min_neighbor_count
max_neighbor_count
mean_neighbor_count
```

are actually body particle counts:

```cpp
air_particles
coil_particles
plate_particles
```

This is misleading.

#### Required cleanup

Rename metrics:

```cpp
min_neighbor_count      -> min_body_particle_count
max_neighbor_count      -> max_body_particle_count
mean_neighbor_count     -> mean_body_particle_count
neighbor_count_finite   -> body_particle_count_finite
```

Rename the helper/test documentation from:

```text
geometry relaxation diagnostic
```

to:

```text
geometry lattice scaffold diagnostic
```

unless real relaxation dynamics are implemented.

Recommended names:

```text
test_3d_aphi_ck_em_geometry_lattice_scaffold_diagnostic
aphi_em_geometry_lattice_scaffold_diagnostic_helpers.h
```

If renaming files is too disruptive, at least update all record/README text:

```text
This test is a lattice geometry scaffold only. It does not perform particle relaxation yet.
```

#### Future real relaxation requirements

A real relaxation test must include:

```text
actual particle relaxation dynamics
before/after particle quality
neighbor count distribution
shape boundary quality
volume conservation
A-phi solve after relaxed particles
comparison against lattice baseline
```

Do not claim this stage has real relaxation.

---

## 4. 10.14-C required tasks

### Task C1 — fix Contact two-body semantics and audit

#### Problem

Previous review found potential ambiguity in two-body Contact naming, especially if a body labeled as air also contains TEAM7-like material tags.

#### Required action

Add bodywise audit to the two-body Contact test as well, not only three-body.

For each body output:

```text
body_name
particle_count_total
particle_count_air
particle_count_conductor
particle_count_source
rhs_l2
joule_air
joule_conductor
joule_source
joule_total_physical_box
sigma_min/max
nu_min/max
```

Do not label a whole body as `air_joule_integral` unless it is guaranteed to contain only air particles.

Use names like:

```text
left_body_joule_integral
right_body_joule_integral
left_air_joule
left_conductor_joule
left_source_joule
right_air_joule
right_conductor_joule
right_source_joule
```

#### Hard gate

Two-body Contact test passes only if:

```text
GMRES converged
fields finite
source RHS norm > 0
total conductor Joule > 0
total air Joule / total conductor Joule < 0.05
interface spike hard gate passed
```

---

### Task C2 — strengthen three-body Contact audit gates

Use the new `AphiThreeBodyContactAuditSummary` or an improved version.

Hard gate:

```text
air body has air particles
coil body has source particles
plate body has conductor particles
coil RHS norm > 0
plate conductor Joule > 0
global air Joule / plate conductor Joule < 0.05
source-only coil Joule ~ 0
all bodywise residuals finite
GMRES converged
interface spike hard gate passed
```

Print all audit metrics in test output and record.

---

### Task C3 — thermal resolution-floor matrix

#### Current status

The source-driven thermal path now uses a resolution-floor fallback. This is acceptable for small temperature rises, but the current test must explicitly prove both branches.

#### Add or modify test

Recommended name:

```text
test_3d_aphi_ck_em_joule_thermal_resolution_floor_diagnostic
```

or extend:

```text
test_3d_aphi_ck_em_joule_thermal_one_way
```

with two source-driven subcases.

#### Required subcases

1. Low-power case:

```text
small impressed_current_amplitude or short dt
expected ΔT below Real precision near T0
fallback_triggered = true
energy_relative_error < 5%
```

2. High-power case:

```text
larger impressed_current_amplitude or larger dt
expected ΔT above precision floor
fallback_triggered = false
actual Temperature field changes visibly
energy_relative_error < 5%
```

#### Required metrics

Add to thermal metrics:

```cpp
bool resolution_floor_triggered;
Real expected_temperature_rise_max;
Real temperature_resolution_floor;
Real observed_temperature_delta_max;
Real joule_deposited_from_field;
Real delta_E_from_temperature_increment;
Real delta_E_effective;
```

#### Gate

```text
low-power: resolution_floor_triggered == true
high-power: resolution_floor_triggered == false
high-power: observed_temperature_delta_max > temperature_resolution_floor
both: energy_relative_error < 5%
```

---

### Task C4 — unify thermal energy increment calculation

The source-driven path uses:

```cpp
sum rho_cp * V * (T - T0)
```

The mapping path should use the same method.

Required change:

```text
Do not use E_end - E_start for mapping-path closure.
Use hostRegionThermalEnergyIncrementFromReference(...) or equivalent.
```

Reason:

```text
avoid future large-number subtraction precision errors.
```

---

### Task C5 — boundary support policy gate small refactor

#### Current issue

If `boundarySupportPolicyDiagnosticPassed()` assumes row order, it is fragile.

Required change:

```text
Find rows by policy enum/name, not by index.
```

Hard gate should apply by policy:

```text
Baseline:
    must converge
    conductor Joule > 0
    air/conductor Joule < 0.05

EnlargedAir:
    all enlarged-air rows converge
    conductor Joule relative change vs baseline < 2%

PassiveAirShell:
    diagnostic only
    sigma shell max <= 1e-14
    shell Joule <= max(1e-12, 1e-6 * conductor_joule)
    conductor Joule relative change vs baseline < 2%
```

Final output:

```text
default_boundary_policy = enlarged_air_domain_padding
passive_air_shell_status = diagnostic_only
ghost_mirror_boundary = not_implemented
legacy_ghost_buffer_diva = decoupled
```

---

### Task C6 — annular source naming or implementation

#### Current uncertainty

The uploaded `aphi_annular_source_region_diagnostic_helpers.h` from the previous package appeared to still use TEAM7-like box layout and not a real annular mask.

#### Required action

Choose exactly one:

### Option A — rename to TEAM7 coil-source diagnostic

If no real annular geometry is implemented, rename/document:

```text
test_3d_aphi_ck_team7_coil_source_region_diagnostic
```

and write:

```text
This is a TEAM7 box/slot source-region diagnostic, not a real annular coil.
```

### Option B — implement real annular/racetrack source mask

If keeping the annular name, implement actual geometry:

```cpp
r = sqrt((x - cx)^2 + (y - cy)^2);
inside_annulus = inner_radius <= r && r <= outer_radius && abs(z - cz) <= z_half_height;
```

or racetrack equivalent.

Hard gate:

```text
source_particle_count > 0
conductor_particle_count > 0
source_rhs_l2 > 0
GMRES converged
fields finite
conductor Joule > 0
air/conductor Joule < 0.05
```

Recommendation for now: choose Option A, then implement real annular geometry in Stage 10.16.

---

### Task C7 — geometry scaffold naming cleanup

Current helper is a lattice scaffold, not relaxation.

Required action:

```text
Rename metrics or comments:
neighbor_count -> body_particle_count
relaxation -> lattice scaffold
```

If file/test names cannot be renamed immediately, update README and record clearly:

```text
This is not particle relaxation. It only verifies three-body lattice geometry construction.
```

Do not claim particle relaxation is complete.

---

### Task C8 — source-driven metrics cleanup

In `AphiSourceDrivenEmSolveMetrics`, add conductor-specific and source-specific metrics:

```cpp
Real max_abs_J_conductor;
Real max_Joule_conductor;
Real source_rhs_l2_source_region;
Real source_Joule_integral;
Real physical_box_Joule_integral;
Real shell_Joule_integral;
```

Update pass gate:

```text
max_abs_J_conductor > threshold
max_Joule_conductor > 0
source_rhs_l2_source_region > 0
conductor_Joule_integral > threshold
air/conductor Joule < 0.05
```

Keep global max values for diagnostic output only.

---

## 5. Tests to run after 10.14-C cleanup

Run at least these tests:

```text
test_3d_aphi_ck_source_driven_em_solve
test_3d_aphi_ck_em_joule_thermal_one_way
test_3d_aphi_ck_em_joule_thermal_resolution_floor_diagnostic  # if added
test_3d_aphi_ck_contact_source_driven_heating_baseline
test_3d_aphi_ck_three_body_contact_source_driven_baseline
test_3d_aphi_ck_boundary_support_policy_diagnostic
test_3d_aphi_ck_annular_source_region_diagnostic  # or renamed team7 source diagnostic
test_3d_aphi_ck_em_geometry_relaxation_diagnostic # or renamed lattice scaffold diagnostic
```

All must print:

```text
passed=1
```

For each test, update the record with:

```text
command
passed value
key metrics
changed files
remaining risk
```

---

## 6. Stage 10.14-C closure criteria

Stage 10.14-C may be closed only after:

```text
1. two-body Contact audit semantics are fixed;
2. three-body Contact audit hard gates are present;
3. thermal low/high-power resolution-floor branches are tested;
4. mapping and source-driven thermal energy increments use the same reference-increment method;
5. boundary policy gate no longer depends on row order;
6. annular source name is fixed or real annular mask is implemented;
7. geometry scaffold is no longer described as true relaxation;
8. source-driven metrics include conductor-specific J/Joule and source-region RHS norms;
9. all listed tests pass.
```

Do not close 10.14-C if the following remain:

```text
source-driven thermal only passes via unlabelled fallback;
two-body body-level Joule metrics still use misleading air/conductor labels;
annular source still ignores annulus parameters while retaining annular name;
geometry relaxation still reports body particle counts as neighbor counts;
passive_air_shell particles can enter physical conductor/Joule metrics.
```

---

## 7. Next stage after 10.14-C

After 10.14-C cleanup, proceed to:

```text
Stage 10.15 — Standard TEAM7 single-resolution alignment
```

### 10.15 main goals

```text
1. Standard TEAM7 geometry dimensions;
2. source-only coil region matching TEAM7 simplification;
3. 50 Hz and 200 Hz runs;
4. standard probe line / observation points;
5. B real/imag or magnitude CSV;
6. first reference comparison table;
7. air-domain padding sensitivity;
8. dp sensitivity if affordable.
```

### Do not include yet

```text
multiresolution;
ghost/mirror boundary;
projection;
production lambda_A;
real helical coil;
full cold-crucible geometry;
particle relaxation A-phi solve.
```

---

## 8. Medium-term roadmap

### Stage 10.16 — real annular/racetrack source and non-box geometry

```text
real annular or racetrack coil/source mask
single-body source-driven solve
then Contact version
then cold-crucible-like cylindrical source/melt geometry
```

### Stage 10.17 — particle relaxation coupled to A-phi

```text
actual relaxation dynamics
particle quality metrics
A-phi solve after relaxation
comparison with lattice baseline
```

### Stage 10.18 — multiresolution

```text
single-resolution TEAM7 must pass first
then far-air coarsening
then conductor/coil fine resolution
then Contact + multiresolution
```

### Stage 10.19 — Maxwell-Ampere MMS

This may be done earlier if needed for paper rigor, but it is not blocking 10.14-C cleanup.

---

## 9. Final wording for the record

Use this wording in the updated record:

```text
Stage 10.14-B established the source-driven A-phi electromagnetic heating baseline. Stage 10.14-C focuses on cleanup and robustness before standard TEAM7 alignment. The source-driven GMRES pipeline, three-body Contact scaffold, passive-air-shell diagnostics, and Joule-to-thermal coupling are now present. However, two-body Contact semantics, thermal resolution-floor branch coverage, annular-source naming, geometry-scaffold terminology, and boundary-gate robustness must be cleaned before claiming a complete Stage 10.14 closure. Standard TEAM7 reference validation remains a Stage 10.15 task.
```

