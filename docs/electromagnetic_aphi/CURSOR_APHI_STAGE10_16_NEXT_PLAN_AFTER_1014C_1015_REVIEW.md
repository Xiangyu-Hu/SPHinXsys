# CURSOR_APHI_STAGE10_16_NEXT_PLAN_AFTER_1014C_1015_REVIEW.md

> Branch: `feature/electromagnetic`  
> Scope: SPHinXsys SYCL CK A-phi electromagnetic heating solver  
> Context: Review after Stage 10.14-C cleanup and Stage 10.15 single-resolution TEAM7 alignment  
> Purpose: summarize current status, identify remaining gaps, and define the next development stage for Cursor.  
> Required style: continue following SPHinXsys SYCL CK coding patterns.

---

## 0. Executive conclusion

After reviewing the Stage 10.14-C and Stage 10.15 source files, helper files, records, test cases, and Cursor notes, the current status is:

```text
Stage 10.14-C:
    Can be closed.

Stage 10.15:
    Can be closed as "single-resolution TEAM7 scaffold/alignment".

Standard TEAM7 quantitative validation:
    Not closed yet.

Production A-phi electromagnetic heating solver:
    Not final yet, but the source-driven EM -> E/J/Joule -> one-way thermal baseline is now established.
```

The project has moved beyond "can the solver run?" and is now in the stage of:

```text
benchmark realism,
source geometry realism,
probe/reference consistency,
Contact-vs-single-body comparison,
frequency-response interpretation,
and solver robustness.
```

The next stage should be:

```text
Stage 10.16 — Real source geometry + TEAM7 reference/probe pipeline hardening
```

Do **not** move to multiresolution yet.  
Do **not** implement projection.  
Do **not** enable production lambda_A.  
Do **not** implement ghost/mirror EM boundaries.  
Do **not** build full cold-crucible geometry yet.

---

## 1. Frozen assumptions for the next stage

These assumptions remain frozen unless explicitly changed by the project owner:

```text
production lambda_A = off
Contact A-penalty = research only
Contact A-penalty stencil = InnerOnly when used
projection = deferred
ghost/mirror EM boundary = deferred
passive_air_shell = diagnostic only
default EM outer-domain policy = enlarged air-domain padding
thermal boundary conditions = not the current focus
multiresolution = deferred
real helical coil = deferred
standard TEAM7 validation = not yet completed
```

Current coil/source modeling principle:

```text
TEAM7-style coil/source should be modeled as a prescribed current source region.
Do not model real helical copper wire at this stage.
Do not add circuit coupling at this stage.
Do not add coil self-heating at this stage.
```

---

## 2. Current progress summary

### 2.1 Stage 10.14-C completed items

The following cleanup tasks have been completed and may be marked as closed:

```text
C1: two-body Contact bodywise audit and corrected air/Joule semantics
C2: three-body Contact audit hard gate
C3: Joule->thermal low/high-power resolution-floor tests
C4: mapping/source-driven thermal energy increment unified as sum(rho_cp * V * (T - T0))
C5: boundary policy gate now finds rows by policy name/type
C6: annular-source diagnostic clarified as TEAM7 coil-slot semantic, annular_geometry_implemented=0
C7: geometry relaxation diagnostic clarified as lattice scaffold, not actual relaxation
C8: source-driven metrics now use conductor/source-region masks
```

Key conclusion:

```text
10.14-C can be closed.
```

But do not claim:

```text
full electromagnetic heating solver is final
standard TEAM7 validation is complete
true annular geometry exists
true particle relaxation is complete
```

---

### 2.2 Stage 10.15 completed items

The following Stage 10.15 single-resolution TEAM7 alignment items have been added or run:

```text
TEAM7 physical dimensions smoke
TEAM7 quantitative reference comparison
simplified TEAM7 source-driven 50 Hz
simplified TEAM7 dual-frequency 50/200 Hz diagnostic
probe CSV infrastructure
TEAM7 dp convergence diagnostic
TEAM7 physical dp observable diagnostic
TEAM7-like solver study diagnostic
```

A new dual-frequency path exists with CSV outputs, including 50 Hz and 200 Hz probe reports. Recent reported key values:

```text
conductor_joule_ratio_200_to_50 ≈ 1.048
max_abs_E_ratio_200_to_50 ≈ 3.905
probe_csv_written = 1
passed = 1
```

Key conclusion:

```text
Stage 10.15 can be closed only as "single-resolution TEAM7 scaffold/alignment".
```

Do not write:

```text
standard TEAM7 quantitative validation completed.
```

---

## 3. SPHinXsys SYCL CK style requirements

All future code must continue to follow SPHinXsys SYCL CK conventions.

### 3.1 Required patterns

Use existing CK dynamics patterns:

```cpp
StateDynamics<MainExecutionPolicy, SomeLocalDynamicsCK>
InteractionDynamicsCK<MainExecutionPolicy, SomeInteractionDynamicsCK<Inner<> / Contact<>>>
AphiMatrixFreeSolveCK<MainExecutionPolicy>
AphiMultiBodyContactGMRESSolverCK<MainExecutionPolicy>
```

For new CK local dynamics:

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
            // device-side operation
        }

      protected:
        DeviceOnlyVariableDataType variable_;
    };

  protected:
    DiscreteVariable<Real> *dv_variable_;
};
```

### 3.2 Device/host sync

For host-side diagnostics, CSV writing, probe metrics, or energy/Joule integrals, always sync device data to host first:

```cpp
syncVariableToHost<Real>(particles, "VariableName");
syncVariableToHost<Vecd>(particles, "VectorVariableName");
syncAphiBlockToHost(particles, block_names);
```

If host code modifies variables used later by device dynamics, sync them back:

```cpp
syncVariableToDevice<Real>(particles, "VariableName");
syncVariableToDevice<Vecd>(particles, "VectorVariableName");
```

### 3.3 Keep diagnostics in test/helper layer

The following must remain in `tests/extra_source_and_tests` until stable:

```text
source-driven diagnostics
TEAM7 physical/reference diagnostics
thermal resolution-floor diagnostics
physical-region audit helpers
boundary support policy diagnostics
solver study diagnostics
dp observable diagnostics
probe CSV diagnostics
```

Do not move these into core production SPHinXsys modules yet.

---

## 4. Current known gaps and required decisions

### 4.1 200 Hz air/conductor Joule ratio

Current issue:

```text
50 Hz air/conductor Joule ratio ≈ 0.0355
200 Hz air/conductor Joule ratio ≈ 0.4719
```

The 200 Hz test currently uses a relaxed gate:

```text
max_air_to_conductor_joule_ratio = 1.0
```

This is acceptable only if the test remains diagnostic.

#### Required action

Keep the existing dual-frequency test as:

```text
diagnostic-only frequency response test
```

Do **not** use it as a hard physical Joule gate.

Add or plan a separate strict physical gate:

```text
test_3d_aphi_ck_simplified_team7_200hz_physical_joule_gate
```

This test should explicitly decide one of two acceptable interpretations:

#### Option A: strict nonconductive-air Joule gate

```text
air/conductor Joule ratio < 0.05
```

If this fails, inspect whether the air conductivity used for numerical conditioning is creating artificial Joule heat.

#### Option B: exclude numerical air Joule from physical heating target

If air uses tiny nonzero sigma for conditioning, do not treat air Joule as physical heat. Then the test must clearly report:

```text
air_joule_numerical_artifact = true
physical_heating_target = conductor Joule only
```

Do not mix diagnostic air Joule with physical heating validation.

---

### 4.2 centerline profile mismatch

Current issue:

```text
centerline_profile_rel ≈ 1.38
```

The current quantitative reference comparison still passes, treating this as informative.

#### Required action

Keep `centerline_profile_rel` informative in Stage 10.15.  
Do not make it a hard gate yet.

However, Stage 10.16 must include a dedicated probe/reference audit.

Checklist:

```text
confirm probe line coordinates
confirm z-height and plate-relative position
confirm whether comparing B_real, B_imag, B_magnitude, or component
confirm normalization
confirm reference/self-reference sampling points
confirm interpolation method
confirm corrected B=curlA is used
confirm frequency and current amplitude scaling
```

Expected output:

```text
B_real probe CSV
B_imag probe CSV
B_magnitude probe CSV
component-wise CSV if applicable
reference/candidate aligned CSV
profile_rel by component
```

The goal is to determine whether `centerline_profile_rel≈1.38` comes from:

```text
probe definition mismatch,
normalization mismatch,
component mismatch,
interpolation mismatch,
insufficient dp,
or actual solver error.
```

Do not change the solver until the probe/reference pipeline is audited.

---

### 4.3 SYCL thermal high-power branch

Current issue:

The high-power branch proving observable temperature change uses:

```text
thermal_use_uniform_joule
```

rather than a higher-current source-driven EM Joule field.

This is acceptable for closing Stage 10.14-C because it proves:

```text
SYCL thermal kernel can produce observable temperature increments when the source is above Real precision.
```

But it does not fully prove:

```text
source-driven high-current EM Joule -> thermal observable temperature increment
```

#### Required action

Add a later test:

```text
test_3d_aphi_ck_source_driven_high_current_thermal_observable
```

Requirements:

```text
increase impressed_current_amplitude or thermal window
do not use uniform Joule
source-driven EM solve produces Joule
thermal update uses actual Joule field
resolution_floor_triggered = false
observed_temperature_delta_max > temperature_resolution_floor
energy_relative_error < 5%
```

This is recommended for Stage 10.16 or 10.15-post-cleanup, not a blocker for closing 10.15.

---

### 4.4 Diagnostic tests with converged=0

The following tests may include diagnostic subcases with `converged=0`:

```text
test_3d_aphi_ck_gmres_team7_like_solver_study_diagnostic
test_3d_aphi_ck_gmres_team7_physical_dp_observable_diagnostic
```

This is acceptable if and only if:

```text
the test is explicitly diagnostic-only;
failed/tight/unpreconditioned configurations are reported as study data;
these tests are not part of release hard gates.
```

#### Required action

Update records and README text:

```text
This test is diagnostic-only. Some solver configurations may report converged=0. The test passes if the diagnostic report is generated and expected baseline configurations are evaluated.
```

---

## 5. Current Contact status and next work

### 5.1 Current Contact status

Contact now has:

```text
two-body source-driven Contact baseline
three-body air/coil/plate Contact scaffold
bodywise audit
source RHS audit
air/conductor/source Joule audit
interface spike hard gate
lambda_A off
```

The three-body helper is now confirmed to be a true three-body setup:

```text
AirBody: air slabs / background domain
CoilBody: coil/source box
PlateBody: conductor/plate box
Contacts:
    AirBody contacts CoilBody and PlateBody
    CoilBody contacts AirBody
    PlateBody contacts AirBody
```

This is acceptable as the current TEAM7-like three-body scaffold.

### 5.2 Remaining Contact gaps

Contact is not final production closure because the following are still missing:

```text
Contact A-penalty three-body gate
standard TEAM7 reference comparison under Contact
Contact vs single-body comparison
Contact dp convergence
Contact with real annular/racetrack source geometry
Contact with real non-box geometry
multiresolution Contact
```

Do not enable Contact A-penalty in production. Keep:

```text
lambda_A = off
```

### 5.3 Required next Contact task

Add a comparison test:

```text
test_3d_aphi_ck_team7_single_body_vs_three_body_contact_comparison
```

Purpose:

```text
compare single-body region-tagged TEAM7 solve
against three-body Contact TEAM7 solve
using the same physical setup
```

Metrics:

```text
GMRES convergence
GMRES iterations
final residual
conductor Joule
air Joule
source RHS norm
B probe line
J probe line
Joule probe or region integral
centerline/profile rel difference
interface spike metrics
```

Hard/initial gates:

```text
both solves converge
both fields finite
both conductor Joule > 0
relative conductor Joule difference < 20% initially
B probe profile relative difference < 30% initially
interface spike hard gate passes
```

The initial tolerances are intentionally loose. Tighten them after the comparison is stable.

---

## 6. Boundary strategy going forward

### 6.1 Current boundary strategy

The current EM outer-domain policy is:

```text
default: enlarged air-domain padding
diagnostic: passive_air_shell
deferred: ghost/mirror boundary
```

This remains correct.

### 6.2 Do not implement EM Dirichlet/Neumann/Robin now

Do not implement electromagnetic Dirichlet/Neumann/Robin boundary APIs yet. For A-phi, these are not simple thermal boundary conditions; they affect:

```text
operator rows
RHS
preconditioner
GMRES boundary DOFs
A/phi complex block
B=curlA postprocess near boundary
Contact interface behavior
```

Current strategy remains:

```text
large/enlarged air domain + padding sensitivity
```

### 6.3 Passive air shell

`passive_air_shell` remains diagnostic only.

It must satisfy:

```text
shell sigma max <= 1e-14
shell Joule <= max(1e-12, 1e-6 * conductor_joule)
physical metrics exclude shell particles
```

---

## 7. Stage 10.16 recommended plan

Stage 10.16 should be:

```text
Stage 10.16 — Real source geometry + TEAM7 reference/probe pipeline hardening
```

### 7.1 P1 — real annular/racetrack source geometry

Current "annular" diagnostic was correctly renamed/clarified as TEAM7 coil-slot semantic:

```text
annular_geometry_implemented = 0
```

Next step is to implement real annular or racetrack source geometry.

#### Requirements

Use source-only coil:

```text
sigma_source = 0
source_rhs_l2_source_region > 0
source Joule not counted as physical conductor Joule
```

Geometry options:

#### Option A: annular source

```cpp
r = sqrt((x - cx)^2 + (y - cy)^2);
inside_source =
    inner_radius <= r &&
    r <= outer_radius &&
    abs(z - cz) <= z_half_height;
```

#### Option B: racetrack source

Implement a rectangular-loop/racetrack mask using straight segments + rounded ends.

Start with Option A unless TEAM7 reference geometry requires racetrack-like coil.

#### Test

```text
test_3d_aphi_ck_real_annular_source_region_diagnostic
```

Hard gate:

```text
source_particle_count > 0
conductor_particle_count > 0
source_rhs_l2_source_region > 0
source_Joule_integral <= max(1e-12, 1e-6 * conductor_Joule)
GMRES converged
fields finite
max_abs_J_conductor > threshold
conductor_Joule_integral > threshold
```

Start with single-body. Add Contact version later.

---

### 7.2 P2 — TEAM7 probe/reference pipeline audit

Do not modify the solver first. Audit the measurement pipeline.

Create:

```text
test_3d_aphi_ck_team7_probe_reference_alignment_diagnostic
```

Output:

```text
probe coordinates
probe z-height
probe component choice
B_real
B_imag
B_magnitude
candidate/reference aligned table
profile_rel by component
normalization factor
frequency
current amplitude
dp
boundary_width_scale
```

This test should explain why:

```text
centerline_profile_rel ≈ 1.38
```

is currently large.

No hard gate on profile error yet. The pass condition is:

```text
CSV written
all values finite
alignment metadata written
```

After interpretation, decide which profile metric becomes hard gate for standard TEAM7 validation.

---

### 7.3 P3 — 200 Hz physical Joule policy

Create a separate test or diagnostic mode:

```text
test_3d_aphi_ck_simplified_team7_200hz_physical_joule_gate
```

Decide and document one policy:

#### Policy A: strict nonconductive-air Joule

```text
air/conductor Joule < 0.05
```

If this fails, inspect air sigma and Joule calculation.

#### Policy B: conductor-only physical Joule

If air sigma is nonzero only for numerical conditioning:

```text
air Joule is diagnostic only
physical Joule = conductor Joule
```

Then pass gate should not use air Joule as physical heat. But it must print the numerical air Joule ratio.

Do not silently use `max_air_to_conductor_joule_ratio=1.0` in a hard validation test.

---

### 7.4 P4 — Contact vs single-body TEAM7 comparison

See Section 5.3.

This is required before claiming Contact-based TEAM7 validity.

---

### 7.5 P5 — source-driven high-current thermal observable

Create:

```text
test_3d_aphi_ck_source_driven_high_current_thermal_observable
```

Requirements:

```text
actual source-driven EM Joule field
no uniform Joule replacement
resolution_floor_triggered = false
observed_temperature_delta_max > resolution floor
energy_relative_error < 5%
```

This closes the current high-power thermal engineering compromise.

---

## 8. Things explicitly deferred

Do not implement these in Stage 10.16 unless the project owner explicitly changes priority:

```text
multiresolution
projection
production lambda_A
ghost/mirror EM boundary
full cold-crucible geometry
real helical coil
coil circuit coupling
coil self-heating
Contact A-penalty production
full particle-relaxed A-phi production case
```

---

## 9. Tests to keep as hard gates

The following should remain hard-gate tests:

```text
test_3d_aphi_ck_source_driven_em_solve
test_3d_aphi_ck_em_joule_thermal_one_way
test_3d_aphi_ck_contact_source_driven_heating_baseline
test_3d_aphi_ck_three_body_contact_source_driven_baseline
test_3d_aphi_ck_boundary_support_policy_diagnostic
test_3d_aphi_ck_simplified_team7_source_driven
test_3d_aphi_ck_probe_metric_csv
```

The following may remain diagnostic-only:

```text
test_3d_aphi_ck_simplified_team7_dual_frequency_diagnostic
test_3d_aphi_ck_gmres_team7_like_solver_study_diagnostic
test_3d_aphi_ck_gmres_team7_physical_dp_observable_diagnostic
```

If a diagnostic test contains `converged=0`, it must not be part of release hard gates.

---

## 10. Record updates required

Update:

```text
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_14_REAL_EM_HEATING_SOLVER_RECORD.md
```

Add a new section:

```text
## Stage 10.15 Closure Review
```

Must state:

```text
10.14-C closed.
10.15 single-resolution TEAM7 scaffold/alignment closed.
Standard TEAM7 quantitative validation not closed.
Known gaps:
    200 Hz air Joule ratio high;
    centerline_profile_rel informative gap;
    high-power thermal uses uniform Joule branch;
    solver study contains diagnostic non-converged cases;
    real annular/racetrack source not implemented yet;
    multiresolution deferred.
Next stage:
    10.16 real source geometry + TEAM7 probe/reference hardening.
```

Add a new file:

```text
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_16_PLAN.md
```

It should include this document's Stage 10.16 tasks.

---

## 11. Final wording for project status

Use this wording in summaries:

```text
Stage 10.14-C is closed. Stage 10.15 is closed as a single-resolution TEAM7 scaffold/alignment stage. The source-driven A-phi solver, E/J/Joule postprocessing, one-way Joule-to-thermal coupling, Contact baselines, boundary policy diagnostics, and probe CSV infrastructure are now established in the SPHinXsys SYCL CK style. However, standard TEAM7 quantitative validation is not yet closed. Known gaps include 200 Hz air Joule interpretation, centerline profile mismatch, diagnostic-only solver sweep non-convergence, and the absence of real annular/racetrack source geometry. The next recommended stage is 10.16: real source geometry and TEAM7 reference/probe pipeline hardening.
```

