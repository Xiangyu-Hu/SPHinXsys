# CURSOR_APHI_STAGE10_14_BOUNDARY_SOURCE_CONTACT_NEXT_ACTIONS.md

> Applicable branch: `feature/electromagnetic`  
> Scope: SPHinXsys SYCL CK A-phi electromagnetic heating solver  
> Current goal: Continue from the Stage 10.14 baseline toward a "usable source-driven A-phi electromagnetic heating solver", prioritizing boundary policy, TEAM7-style prescribed source, Contact three-body composition audit, Joule→thermal energy closure, and non-rectangular geometry preparation.  
> Important principle: This stage does not pursue a customer-facing demo; prioritize real solver functional closure and verifiability.

---

## 0. Purpose of this document

The previous round completed several Stage 10.14 baselines:

- `source-driven A-phi GMRES` baseline;
- `E/J/Joule` postprocessing baseline;
- `Joule -> thermal one-way` preliminary baseline;
- two-body / three-body Contact source-driven baseline;
- boundary support policy diagnostic initial version;
- simplified TEAM7 scaffold initial version;
- probe/metric CSV initial version.

However, we cannot yet claim that a complete A-phi electromagnetic heating solver is finished. The following key issues must still be addressed:

1. Definition and standardization of TEAM7-style prescribed current source;
2. Distinction between enlarged air domain and passive air shell boundary policies;
3. Physical composition audit of the Contact three-body case;
4. Gradually upgrading interface spike from warning to hard gate;
5. Raising Joule -> thermal energy closure from "order-of-magnitude pass" to strict closure;
6. Preparing non-rectangular/annular coil/source geometry for standard TEAM7 and cold crucible;
7. Explicitly deferring ghost/mirror, electromagnetic Dirichlet/Neumann/Robin from the current stage.

This document requires Cursor to execute tasks in order, avoiding Auto model misinterpretation or premature implementation of features that should not be built yet.

---

## 1. Terminology freeze: these concepts must be strictly distinguished

### 1.1 prescribed current source / impressed current RHS

In this document, `prescribed current source` means:

```text
Given current density J_source or equivalent RHS in the coil/source region
-> assemble A-phi equation RHS
-> solve A, phi with GMRES
-> postprocess B, H, E, J, Joule
```

It is not directly prescribing the A field.

Must distinguish:

```text
impressed A demo:
    directly specify A
    only do B/E/J/Joule postprocessing

source-driven solve:
    specify J_source / RHS
    solver solves A, phi itself
    then do B/H/E/J/Joule postprocessing
```

The current main line must be `source-driven solve`, not impressed-A demo.

---

### 1.2 TEAM7 coil/source modeling approach

TEAM7 excitation can be understood as prescribed current source. Standard TEAM7 does not require modeling a real helical copper wire shape; many benchmark/software implementations simplify the coil as a coil/source region with cross-section, then apply an AC current source in that region.

Therefore the first-version TEAM7-style source should use:

```text
coil/source region:
    geometrically an annular or rectangular-cross-section source region;
    no real helical wire needed;
    no circuit coupling;
    no coil self-Joule heating;
    apply prescribed sinusoidal current density in the source region;
    source enters A-phi RHS.
```

Forbidden in the current stage:

```text
real helical copper wire;
circuit coupling;
coil resistive heating;
coil internal current redistribution.
```

---

### 1.3 enlarged air domain / enlarged air padding

`enlarged air domain` or `enlarged air padding` means:

```text
Expand the air computational domain participating in electromagnetic solve;
outer expanded region is air;
sigma = 0;
nu = 1 / mu0;
this region participates in GMRES solve;
fields in this region are not the primary observation target;
final metrics only count inside plate / conductor / probe / physical box.
```

Two purposes:

1. Physically simulate a larger open air domain;
2. Numerically push outer boundary kernel support deficiency / truncated far boundary farther away, so plate/coil/probe regions are not strongly affected by the outer boundary.

This is not ghost boundary, nor non-physical particles. It is the default outer-boundary strategy for current TEAM7/cold crucible A-phi.

---

### 1.4 passive air shell / support shell

The strategy formerly called `dummy_shell` in the code must be renamed or at least documented as:

```text
passive_air_shell
```

Its actual meaning is:

```text
Expand BoundingBox;
generate ordinary lattice particles;
set sigma = 0 for particles outside the physical box;
do not use SphinxSys ghost/dummy BC API;
no mirror ghost;
no analytic boundary value imposed.
```

It is essentially:

```text
passive sigma=0 support padding
```

It can be used to diagnose boundary kernel support impact on internal results, but should not be the current default production boundary policy.

Requirements:

```text
passive_air_shell is diagnostic only;
do not call it ghost;
do not call it wall dummy;
do not call it formal Dirichlet/Neumann/Robin boundary conditions;
all formal physical metrics must explicitly exclude shell particles.
```

---

### 1.5 ghost/mirror boundary

`ghost/mirror` means actually constructing mirror particles or ghost values to impose boundary conditions, e.g. Neumann/Dirichlet/Robin.

Forbidden in the current stage:

```text
A-phi ghost/mirror boundary;
analytic ghost boundary;
production ghost particle boundary;
mirror ghost particle boundary.
```

The existing `aphi_ghost_buffer_diva_diagnostic_helpers.h` in the repository may only remain as a divA research diagnostic line; it must not be connected to source-driven heating / TEAM7 / production A-phi solver.

---

### 1.6 Thermal boundary is not the focus of this stage

Thermal boundary SPHinXsys diffusion/thermal module already has a foundation. Do not spend main effort redesigning thermal Dirichlet / Neumann / Robin in this stage.

Current focus is:

```text
A-phi electromagnetic air domain truncation;
Contact material interface;
source-driven A-phi;
Joule -> thermal energy closure.
```

---

## 2. Current solver status: allowed wording

### 2.1 Baselines that can be claimed as complete

May write in records:

```text
Stage 10.14 baseline has been established:
    source-driven A-phi GMRES baseline passed;
    E/J/Joule postprocessing baseline passed;
    two-body and three-body Contact source-driven baselines passed;
    boundary support policy diagnostic initial version passed;
    simplified TEAM7 scaffold initial version passed;
    probe/metric CSV initial version passed.
```

### 2.2 Do not overstate

Forbidden to write:

```text
complete A-phi electromagnetic heating solver is finished;
Contact production is closed;
standard TEAM7 is validated;
boundary conditions fully support Dirichlet/Neumann/Robin;
ghost/mirror boundary is implemented;
Joule -> thermal energy closure strictly passed.
```

Accurate current wording:

```text
solver baseline is working, but production/benchmark closure is still ongoing.
```

---

## 3. Next-round task order

Proceed in the following order; do not skip steps.

```text
P0: documentation and naming corrections
P1: Contact three-body composition audit
P2: Boundary policy remediation and stronger gate
P3: TEAM7-style prescribed current source standardization
P4: Joule -> thermal energy closure remediation
P5: Contact interface spike gate strengthening
P6: Annular / racetrack source geometry minimal case
P7: Particle relaxation / non-rectangular geometry preparation
P8: Maxwell-Ampere MMS implementation plan enters next round
```

P0-P5 are the closure items that must be prioritized now; P6-P7 prepare TEAM7/cold crucible geometry capability; P8 may be design or minimal implementation, but must not affect P0-P5.

---

## 4. P0: Documentation and naming corrections

### 4.1 Update wording in records and documentation

In all Stage 10.14 related records, change:

```text
dummy_shell
```

to or annotate as:

```text
passive_air_shell / passive sigma=0 support shell
```

If code variable names are not changed yet, README/record must clearly state:

```text
dummy_shell in old code = passive_air_shell diagnostic;
it is NOT SphinxSys ghost/dummy boundary condition.
```

### 4.2 Must write into handoff

Add or update:

```text
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_14_BOUNDARY_SOURCE_CONTACT_HANDOFF.md
```

Must include:

```text
Default EM outer-domain policy: enlarged air-domain padding.
Optional diagnostic policy: passive sigma=0 air shell.
No ghost/mirror boundary in Stage 10.14.
Legacy aphi_ghost_buffer_diva is divA diagnostic only.
TEAM7-style source is prescribed current source, not real helical coil.
Thermal boundary is not the focus of Stage 10.14.
```

### 4.3 Acceptance

Documentation grep check:

```bash
grep -R "dummy_shell" tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_14* -n
```

If `dummy_shell` still appears, there must be an adjacent explanation:

```text
legacy name; actually passive_air_shell, not ghost BC.
```

---

## 5. P1: Contact three-body composition audit

### 5.1 Purpose

There is currently a three-body Contact source-driven baseline, but its actual physical composition must be clarified. Do not just say "three-body"; must output particle distribution of air/conductor/source in each body.

### 5.2 Must add bodywise/regionwise metrics

Output in three-body Contact helper:

```text
number_of_bodies
body_name
particle_count_total
particle_count_air
particle_count_conductor
particle_count_source
particle_count_shell_if_any
sigma_min
sigma_max
nu_min
nu_max
source_rhs_integral_real
source_rhs_integral_imag
source_rhs_l2
joule_integral_total
joule_integral_air
joule_integral_conductor
joule_integral_source
max_J_by_region
max_Joule_by_region
```

Suggested implementation locations:

```text
tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_contact_source_driven_heating_baseline_helpers.h
tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_em_observable_helpers.h
```

### 5.3 Three-body physical composition requirements

Three-body scaffold target composition should be as close as possible to:

```text
body 1: air / background magnetic domain
body 2: conductive plate or melt
body 3: source/coil region
```

If each body in actual code uses the same TEAM7 layout tagging, must still output counts proving:

```text
source particles exist;
conductor particles exist;
air particles exist;
Joule is mainly in conductor;
air Joule is near 0;
source RHS mainly from source region.
```

### 5.4 hard gate

Three-body Contact baseline must satisfy:

```text
GMRES converged = true
all field values finite
total source_rhs_l2 > 0
conductor particle count > 0
source particle count > 0
air particle count > 0
conductor Joule integral > 0
air Joule integral / conductor Joule integral < 0.05
```

If coil/source is built as source-only region:

```text
source Joule is not a hard gate;
but source_rhs_l2 must > 0.
```

---

## 6. P2: Boundary policy remediation and stronger gate

### 6.1 Default policy

Current default outer-boundary policy must be:

```text
enlarged air-domain padding
```

Not passive_air_shell, nor ghost/mirror.

Tests must retain:

```text
boundary_width_scale = 1
boundary_width_scale = 2
boundary_width_scale = 3
```

or equivalent air-domain scale sweep.

### 6.2 passive_air_shell is diagnostic only

Retain one optional test policy:

```text
passive_air_shell_width_scale = 2
```

But it is for diagnosis only, not the default TEAM7 scaffold policy.

### 6.3 Explicit physical-box mask

All physical statistics must explicitly use physical box / physical region mask. Passive shell particles must not enter formal conductor/Joule/probe statistics.

Must implement or confirm:

```text
isInsideTeam7PhysicalBox(position)
isInsidePhysicalAirDomain(position)
isInsidePassiveAirShell(position)
team7ParticleInRegion(position, region_id) must not include shell particles
```

If existing `team7ParticleInRegion()` cannot exclude shell, add wrapper:

```cpp
bool team7ParticleInPhysicalRegion(const Vecd& p, int region_id)
{
    return isInsideTeam7PhysicalBox(p) && team7ParticleInRegion(p, region_id);
}
```

### 6.4 Required boundary metrics

In `test_3d_aphi_ck_boundary_support_policy_diagnostic`, output per policy:

```text
policy_name
boundary_width_scale
passive_air_shell_width_scale
gmres_converged
gmres_iterations
final_residual
particle_count_total
particle_count_physical_box
particle_count_shell
conductor_particle_count
source_particle_count
air_particle_count
shell_sigma_max
shell_joule_integral
conductor_joule_integral
air_joule_integral
source_rhs_l2
max_A
max_B
max_E
max_J
max_Joule
probe_B_metric_if_available
probe_J_metric_if_available
relative_change_conductor_joule_vs_baseline
relative_change_probe_B_vs_baseline
relative_change_probe_J_vs_baseline
```

### 6.5 stronger hard gate

`boundarySupportPolicyDiagnosticPassed()` must change from weak gate to strong gate.

Requirements:

```text
baseline GMRES converged
all enlarged_air policies GMRES converged
all fields finite
all conductor Joule integral > 0
air Joule / conductor Joule < 0.05
passive shell sigma max == 0 or <= 1e-14
passive shell Joule integral <= max(1e-12, conductor_joule * 1e-6)
boundary_width_scale=2 conductor Joule relative change < 2%
boundary_width_scale=3 conductor Joule relative change < 2%
passive_air_shell conductor Joule relative change < 2%  # diagnostic hard gate; if fails, mark passive_air_shell rejected
```

Probe B/J relative change: print warning first, no hard fail, unless stable probe framework exists.

### 6.6 Conclusion output

Test end must print:

```text
default_boundary_policy = enlarged_air_domain_padding
passive_air_shell_status = diagnostic_only / accepted_as_diagnostic / rejected
ghost_mirror_boundary = not_implemented
legacy_ghost_buffer_diva = decoupled_from_heating_solver
```

---

## 7. P3: TEAM7-style prescribed current source standardization

### 7.1 Goal

Current source-driven case must align with TEAM7 excitation:

```text
source/coil region = simplified annular or racetrack-like source body/region
source = prescribed sinusoidal current density
frequency = 50 Hz and 200 Hz eventually
```

At minimum in current stage:

```text
source RHS has clear units;
spatial integral of source RHS can be output;
source region particle count can be output;
source orientation can be output;
source-only region does not count its own Joule.
```

### 7.2 Do not build real helical wire

Currently forbidden:

```text
real helical coil geometry
multi-turn wire geometry
circuit coupling
coil self-heating
```

### 7.3 source-only coil strategy

First-version TEAM7-style source uses:

```text
source region:
    sigma = 0 or treated as non-conductive for Joule statistics
    prescribed J_source != 0
    contributes to A-phi RHS
    not counted as conductor Joule
```

### 7.4 source metrics hard gate

Each source-driven test must output:

```text
source_particle_count
source_rhs_l2
source_rhs_integral_real
source_rhs_integral_imag
source_rhs_direction
frequency
omega
```

hard gate:

```text
source_particle_count > 0
source_rhs_l2 > 0
omega > 0
```

---

## 8. P4: Joule -> thermal energy closure remediation

### 8.1 Current problem

Current `Joule -> thermal one-way` is not a CPU/GPU comparison issue, but:

```text
EM solve produces JouleHeatSource
-> thermal update uses JouleHeatSource
-> thermal energy increment inconsistent with ∫Joule dV * dt
```

Current order-of-magnitude gate passes, but strict 5% gate does not. This must be remediated.

### 8.2 Investigation order

Must check item by item:

```text
1. Whether Joule power integration region matches thermal energy integration region exactly;
2. Whether conductor / air / shell / source particles are mixed into energy statistics;
3. Whether VolumetricMeasure is consistent;
4. Whether JouleHeatSource unit is W/m^3;
5. Whether phasor 0.5 factor is used only once;
6. Whether AdiabaticJouleTemperatureStepCK heats only the target region;
7. Whether rho, cp are assigned correctly per region;
8. Whether host-device sync completes before integration;
9. Whether same-body spatial Joule and mapping-only uniform Joule use consistent statistics.
```

### 8.3 Suggested mask-based thermal step

If current thermal update runs on all particles, suggest adding or temporary test:

```cpp
AdiabaticJouleTemperatureStepInRegionCK
```

Or control with mask in helper:

```text
add Joule only to conductor/physical heated region;
thermal energy also counts only the same mask.
```

### 8.4 hard gate

Two levels:

#### P4-a uniform Joule mapping

```text
energy_relative_error < 1%
```

#### P4-b source-driven spatial Joule

First round:

```text
energy_relative_error < 5%
```

Subsequent target:

```text
energy_relative_error < 2%
```

If not yet achievable, must output failed reason; do not close Stage 10.14 with 0.25x~25x order-of-magnitude gate.

---

## 9. P5: Contact interface spike gate strengthening

### 9.1 Current problem

Current interface spike metric is only warning, not hard gate. Needs strengthening.

### 9.2 Must output

For two-body and three-body Contact source-driven tests, output:

```text
interface_particle_count
near_interface_particle_count
bulk_conductor_particle_count
max_J_interface
mean_J_interface
max_J_bulk
mean_J_bulk
J_spike_ratio_max
J_spike_ratio_mean
max_Joule_interface
mean_Joule_interface
max_Joule_bulk
mean_Joule_bulk
Joule_spike_ratio_max
Joule_spike_ratio_mean
```

### 9.3 hard gate

Initial hard gate:

```text
no NaN / Inf
J_spike_ratio_max < 50
Joule_spike_ratio_max < 50
```

warning gate:

```text
J_spike_ratio_max > 10 -> warning
Joule_spike_ratio_max > 10 -> warning
```

If hard gate fails, cannot claim Contact baseline passed.

### 9.4 Keep lambda_A off

Do not open A-divergence penalty to suppress interface spike.

Keep:

```text
production lambda_A = off
Contact A-penalty = research only
Contact penalty stencil = InnerOnly when research enabled
projection = deferred
```

---

## 10. P6: Annular / racetrack source geometry minimal case

### 10.1 Purpose

Prepare non-rectangular coil/source geometry for TEAM7 and cold crucible.

Many current cases still use box/rectangular tagging. Next step should implement a minimal annular or racetrack source region, but keep the case simple.

### 10.2 Recommended test name

```text
test_3d_aphi_ck_annular_source_region_diagnostic
```

### 10.3 First version: single-body only

First version: no Contact, no thermal, no multiresolution.

Geometry:

```text
air box
conductive block or plate
annular / racetrack source region
```

Solve:

```text
source-driven A-phi GMRES
B/H/E/J/Joule postprocess
VTP output
```

### 10.4 hard gate

```text
source_particle_count > 0
conductor_particle_count > 0
GMRES converged
all fields finite
max_A > threshold
max_B > threshold
max_E > threshold
max_J in conductor > threshold
conductor Joule integral > threshold
air Joule / conductor Joule < 0.05
source RHS l2 > threshold
```

### 10.5 No real helical wire

Still forbidden:

```text
real helical wire
circuit coupling
coil self-heating
```

---

## 11. P7: Particle relaxation / non-rectangular geometry preparation

### 11.1 Purpose

Support body-fitted particles for subsequent:

```text
TEAM7 asymmetric plate with hole
annular coil/source
cold crucible cylindrical melt
segmented crucible wall
```

### 11.2 Current stage: geometry quality test only

Add or plan:

```text
test_3d_aphi_ck_em_geometry_relaxation_diagnostic
```

Do not connect EM solve yet, or only simplest source-free field diagnostic.

### 11.3 Must output

```text
particle_count_before_relaxation
particle_count_after_relaxation
min_neighbor_count
mean_neighbor_count
max_neighbor_count
volume_error
shape_boundary_particle_count
source_region_particle_count
conductor_region_particle_count
```

### 11.4 Gate when connecting to A-phi later

When relaxed particles connect to A-phi, must check:

```text
GMRES iteration count not severely worse than lattice baseline
fields finite
B=curlA diagnostic not severely degraded
Joule integral close to lattice baseline for equivalent geometry
interface spike not worse than threshold
```

---

## 12. P8: Maxwell-Ampere MMS positioning

### 12.1 Current status

`Maxwell-Ampere MMS` is currently in design phase, not implementation complete. It is not a technical blocker, but not yet implemented.

### 12.2 Why it is needed

Ordinary div-free MMS cannot guarantee:

```text
curl H = J
```

Maxwell-Ampere MMS needs to construct self-consistent analytic fields:

```text
B_exact = curl A_exact
H_exact = nu * B_exact
J_exact = curl H_exact
E_exact = J_exact / sigma
E_exact = -i omega A_exact - grad phi_exact
```

### 12.3 Current priority

Before P0-P5 complete, do not let Maxwell-Ampere MMS take the main line.

Recommended order:

```text
complete P0-P5 closure
complete P6 annular source minimal diagnostic
then enter Maxwell-Ampere MMS implementation
```

---

## 13. Near-term roadmap for TEAM7

### 13.1 Phase 1: TEAM7-style simplified source-driven scaffold

Must include:

```text
large/enlarged air domain
asymmetric or simplified conductor plate
source/coil region with prescribed current density
frequency parameter
source-driven GMRES
B/H/E/J/Joule output
probe CSV
boundary_width_scale sensitivity
```

### 13.2 Phase 2: Standard TEAM7 geometry and probe

Then add:

```text
standard aluminum plate dimensions
eccentric hole
standard coil/source dimensions
50 Hz / 200 Hz
standard probe line
reference data comparison
```

### 13.3 Phase 3: Multiresolution

Do not do multiresolution yet. Defer until standard TEAM7 single-resolution runs through.

---

## 14. Near-term roadmap for cold crucible

Cold crucible shares with TEAM7:

```text
prescribed source / coil RHS
large air domain
conductive region Joule
B/H/E/J/Joule output
regionwise metrics
Contact multibody
non-rectangular / annular source geometry
```

Do not jump to full cold crucible complex geometry yet. Do first:

```text
annular source + conductive cylinder/block
then annular source + cylindrical melt
then segmented wall / crucible body
```

---

## 15. Prohibited items summary

Cursor must not do in this stage:

```text
1. Do not implement A-phi projection;
2. Do not enable production lambda_A;
3. Do not set Contact A-penalty default to InnerContact;
4. Do not call passive_air_shell ghost boundary;
5. Do not connect legacy aphi_ghost_buffer_diva to TEAM7 heating;
6. Do not implement ghost/mirror boundary;
7. Do not prioritize real helical coil;
8. Do not claim heating solver complete when Joule->thermal strict energy closure not passed;
9. Do not write simplified TEAM7 scaffold as standard TEAM7 validation;
10. Do not let shell/padding particles enter formal conductor Joule statistics.
```

---

## 16. Current recommended closure criteria

Stage 10.14 cannot be fully closed just because baseline passed. Suggested split:

```text
Stage 10.14-A baseline: already passed
Stage 10.14-B closure: pending
```

### 16.1 Stage 10.14-B closure must satisfy

```text
P1 three-body Contact composition metrics complete
P2 boundary policy strong gate passed
P3 prescribed source metrics complete
P4 Joule->thermal energy closure:
    uniform < 1%
    source-driven < 5%
P5 Contact interface spike hard gate passed
P6 annular/racetrack source diagnostic passed or at least implemented and documented
```

### 16.2 May defer to 10.15

```text
standard TEAM7 reference comparison
Maxwell-Ampere MMS full implementation
particle relaxation coupled A-phi solve
multi-resolution
ghost/mirror boundary
A-phi Dirichlet/Neumann/Robin boundary API
full cold crucible geometry
```

---

## 17. Cursor execution hints

Cursor should organize by the following commands/order:

```text
1. Update documentation and naming: dummy_shell -> passive_air_shell explanation
2. Strengthen P5 boundary_support_policy_diagnostic
3. Strengthen P4/P4b Contact bodywise/regionwise metrics
4. Strengthen interface spike hard gate
5. Fix P2 Joule->thermal energy closure
6. Implement annular/racetrack source minimal diagnostic
7. Update STAGE10_APHI_TEST_CASE_INDEX.md
8. Update CURSOR_APHI_STAGE10_14_REAL_EM_HEATING_SOLVER_RECORD.md
9. Add CURSOR_APHI_STAGE10_14_BOUNDARY_SOURCE_CONTACT_HANDOFF.md
```

After each step, must write in record:

```text
what changed
files changed
tests run
passed/failed
remaining risk
next step
```

---

## 18. Final expected wording

When this round completes, may state:

```text
The Stage 10.14 A-phi electromagnetic heating solver has progressed from baseline source-driven solve to a more controlled closure state. The default electromagnetic outer-domain treatment is enlarged air-domain padding, while passive sigma=0 air shell is retained only as a diagnostic support policy. The TEAM7-style excitation is represented by a prescribed current source region rather than a real helical coil. Contact source-driven two-body and three-body baselines are kept with lambda_A off, and interface spike diagnostics are strengthened. Joule-to-thermal one-way coupling is required to pass strict energy closure before claiming a complete electromagnetic heating solver. Standard TEAM7 validation, multiresolution, ghost/mirror boundary, and projection remain deferred.
```
