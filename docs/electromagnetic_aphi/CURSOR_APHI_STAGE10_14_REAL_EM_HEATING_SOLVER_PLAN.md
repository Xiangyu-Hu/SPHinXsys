# CURSOR_APHI_STAGE10_14_REAL_EM_HEATING_SOLVER_PLAN

> Purpose: continue the A-phi electromagnetic module toward a **usable source-driven electromagnetic heating solver**, not merely a visualization demo.  
> Current branch assumed: `feature/electromagnetic`.  
> Main development area: `tests/extra_source_and_tests` and shared A-phi SYCL CK helpers/operators.  
> This file is intended for Cursor/Auto mode. Follow it literally. Do not reinterpret priorities.

---

## 0. Current status summary

We have finished Stage 10.13. The current state is:

```text
A-phi matrix-free core solver: available
GMRES + block-Jacobi path: available
Inner MMS / E/J/Joule chain: mostly stable
B=curlA corrected postprocess: dp-convergent on Az2D, about 0.88% at dp=0.05
Contact two-body baseline: mostly stable
Contact A-divergence penalty: two-body research gate closed with InnerOnly stencil
Cold crucible demo: impressed-A pipeline demo works, but not a solve validation
Projection: design only, not implemented
Three-body Contact A-penalty: deferred
Standard TEAM7: not yet implemented or validated
Joule-to-thermal coupling: not yet closed as a standard solver case
```

The most important decision now:

```text
Stop prioritizing customer-facing demo features.
Prioritize a real electromagnetic heating solver closure:
source-driven A-phi solve -> E/J/Joule -> thermal source -> energy validation -> simplified TEAM7 -> standard TEAM7.
```

---

## 1. Absolute rules / do-not-change list

These rules are mandatory.

### 1.1 Do not implement projection now

Do not implement post-solve projection:

```text
A_new = A - grad(chi)
phi_new = phi + i omega chi
```

Reason: boundary/interface divA is still an open diagnostic dominated by support deficiency. Projection may destroy E/J/Joule consistency if implemented before source, boundary and thermal coupling are stable.

### 1.2 Do not enable production lambda_A

Keep production A-divergence penalty off:

```text
production lambda_A = off
production eta_A = 0
```

The research setting remains:

```text
eta_A = 0.1 primary research value
eta_A = 0.2 optional diagnostic value
eta_A >= 0.2 requires polish_sweeps = 0
```

### 1.3 Do not change Contact penalty default back to InnerContact

The current Contact A-divergence penalty stencil must remain:

```text
Contact main operator K(A,phi): Inner + Contact
A-divergence penalty stencil: InnerOnly
```

Do not use Contact neighbors in grad-div penalty unless a new explicit diagnostic requires it. InnerContact graddiv is known to create large interface pollution.

### 1.4 Do not claim cold crucible demo is validation

The current cold crucible demo is:

```text
impressed-A / analytic-A pipeline smoke test
single-body region tagging
not source-driven GMRES validation
not quantitative cold crucible validation
```

Do not rename it or document it as a verified physical benchmark.

### 1.5 Do not jump directly to standard TEAM7

Standard TEAM7 requires large air domain, coil/source, conductor/air interface, probes, reference data and possibly multi-resolution. First implement simplified source-driven cases and heating closure.

---

## 2. Meaning of “complete electromagnetic heating solver” for this stage

For this project, a usable first-version electromagnetic heating solver means:

```text
Input:
  geometry / region tags / body layout
  material parameters sigma, nu, rho, cp, k
  excitation frequency omega
  coil/source current or source RHS

Solve:
  complex A, complex phi by matrix-free GMRES

Postprocess:
  B = curl(A)
  H = nu B
  E = -i omega A - grad(phi)
  J = sigma E
  q_Joule = 0.5 sigma |E|^2

Thermal coupling:
  q_Joule is mapped into the thermal source term
  rho cp dT/dt = div(k grad T) + q_Joule

Diagnostics:
  GMRES convergence
  finite field checks
  regionwise Joule integral
  energy conservation / temperature-rise check
  VTP fields for ParaView
```

The first complete version may use:

```text
constant sigma
constant nu
fixed frequency
one-way EM -> thermal coupling
no sigma(T)
no nu(T)
no Lorentz force
no fluid feedback
no projection
lambda_A off
```

Later versions can add sigma(T), nu(T), EM-thermal iteration and Lorentz force.

---

## 3. Stage 10.14 goal

Stage 10.14 should be named:

```text
Stage 10.14 — Source-driven electromagnetic heating solver closure
```

Main objective:

```text
Move from impressed-A demonstration to source-driven A-phi solve, then close Joule-to-thermal coupling with energy validation.
```

Deliverables:

1. Source-driven A-phi solve case.
2. JouleHeatSource -> thermal source coupling.
3. Simple conductor heating energy validation.
4. Simplified TEAM7 source-driven scaffold, single-resolution first.
5. TEAM7-oriented probe / region metric infrastructure.
6. Multi-resolution planning document for large TEAM7 air domain.

---

## 4. P0 — preflight cleanup and documentation freeze

### 4.1 Required reading before coding

Read these files first if present:

```text
CURSOR_STAGE10_13_APHI_NEXT_WORK_DETAILED.md
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_13_CHATGPT_HANDOFF.md
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_13_COLD_CRUCIBLE_DEMO_RECORD.md
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_13_FULL_MAXWELL_MMS_DESIGN.md
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_12_CONTACT_A_PENALTY_RECORD.md
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_12_B_H_OBSERVABLE_RECORD.md
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/STAGE10_APHI_TEST_CASE_INDEX.md
```

### 4.2 Create new record file

Create:

```text
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_14_REAL_EM_HEATING_SOLVER_RECORD.md
```

This record must contain sections:

```text
1. Goal
2. Current frozen assumptions
3. P1 source-driven solve
4. P2 Joule-to-thermal coupling
5. P3 simple conductor heating validation
6. P4 simplified TEAM7 scaffold
7. P5 probe / metric infrastructure
8. P6 TEAM7 large-air-domain / multi-resolution plan
9. Test results summary
10. Open issues
11. Do-not-change decisions
```

Update this record after every completed subtask.

### 4.3 Preserve previous freeze decisions

Write these in the new record explicitly:

```text
production lambda_A = off
Contact A-penalty stencil = InnerOnly
eta_A = 0.1 research only
eta_A >= 0.2 uses polish_sweeps = 0
projection = deferred
three-body Contact A-penalty = deferred until source-driven and thermal closure pass
cold crucible demo = demonstration only
```

---

## 5. P1 — source-driven A-phi solve case

### 5.1 Purpose

Build the first standard case that solves A and phi from a source/RHS. It must not rely on directly assigning the final analytic/impressed A field as the solution.

The goal is to move from:

```text
given A -> postprocess
```

to:

```text
given source/current/RHS -> solve A,phi -> postprocess B,H,E,J,Joule
```

### 5.2 Recommended test name

Create a new test case:

```text
tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_source_driven_em_solve/
```

Suggested files:

```text
test_3d_aphi_ck_source_driven_em_solve.cpp
README.md
CMakeLists.txt
```

Suggested helper:

```text
tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_source_driven_em_solve_helpers.h
```

If there is already a suitable source/cold-crucible scaffold helper, reuse it instead of duplicating geometry code. But keep the new test entry clean and documented.

### 5.3 Geometry for first version

Use the simplest stable geometry:

```text
air box
conductive block or conductive plate
source/coil region represented by analytic current density or RHS source
constant sigma and nu
fixed omega
single resolution
```

Do not use complicated cold crucible geometry in P1 unless the simple block case already passes.

### 5.4 Solver settings

Use stable conservative settings:

```text
lambda_A = off
eta_A = 0
phi_gauge_penalty = 100 or the current stable default
GMRES tolerance = 5e-4 first, then optional 1e-4
restart m = current stable value, e.g. 30
polish_sweeps = current safe baseline, but do not use polish to hide divergence
```

Do not turn on Contact A-penalty in this case.

### 5.5 Required outputs

The case must output or print:

```text
passed
converged
num_iterations
final_residual
max_abs_A
max_abs_phi
max_abs_B
max_abs_H
max_abs_E
max_abs_J
max_Joule
air_Joule_integral
conductor_Joule_integral
total_Joule_integral
finite_field_check
```

VTP fields should include at least:

```text
AReal, AImag
PhiReal, PhiImag
BReal, BImag
HReal, HImag
ElectricFieldReal, ElectricFieldImag or ElectricField magnitude
CurrentDensityReal, CurrentDensityImag or CurrentDensity magnitude
JouleHeatSource
Sigma
Nu
MaterialRegionId
```

### 5.6 P1 hard gates

The case passes only if:

```text
passed = 1
converged = 1
all field values finite
max_abs_A > small_threshold
max_abs_E > small_threshold
max_abs_J > small_threshold in conductor
conductor_Joule_integral > 0
air_Joule_integral is small compared with conductor Joule, unless air sigma is intentionally nonzero
VTP file is written and contains required fields
```

### 5.7 P1 documentation

In the Stage 10.14 record, write:

```text
- Geometry description
- Source/RHS definition
- Material parameters
- Frequency
- Solver configuration
- Printed metrics
- Pass/fail status
- Any observed instability
```

---

## 6. P2 — JouleHeatSource to thermal source coupling

### 6.1 Purpose

This is the most important step toward a real electromagnetic heating solver.

The solver must map:

```text
q_Joule = 0.5 sigma |E|^2
```

into the thermal equation:

```text
rho cp dT/dt = div(k grad T) + q_Joule
```

Start with one-way coupling:

```text
solve EM once
compute JouleHeatSource
hold Joule fixed
solve transient heat conduction for a short time
```

Do not add sigma(T) or nu(T) yet.

### 6.2 Recommended test name

Create:

```text
tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_em_joule_thermal_one_way/
```

Suggested helper:

```text
tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_em_joule_thermal_coupling_helpers.h
```

If thermal helpers already exist in SPHinXsys extra tests, reuse them. Do not duplicate the whole thermal solver.

### 6.3 First thermal validation geometry

Use the simplest possible conductive body:

```text
conductive block
constant rho
constant cp
constant k
adiabatic boundary if possible
q_Joule from EM solve or controlled Joule field
```

First run can use the P1 source-driven EM result. If that is unstable, use a controlled fixed Joule field only to validate thermal mapping, but document this clearly.

### 6.4 Energy accounting

Compute:

```text
P_Joule = integral(q_Joule dV)
E_thermal(t) = integral(rho cp T dV)
Delta_E_thermal = E_thermal(t_end) - E_thermal(t_start)
expected_energy = P_Joule * Delta_t
energy_relative_error = |Delta_E_thermal - expected_energy| / max(expected_energy, epsilon)
```

For adiabatic short-time heating, expect:

```text
Delta_E_thermal ≈ P_Joule * Delta_t
```

If there are boundary heat losses, report them separately and do not use strict energy conservation until losses are accounted for.

### 6.5 P2 hard gates

The case passes only if:

```text
thermal solver runs without NaN/Inf
average conductor temperature increases
P_Joule > 0
Delta_E_thermal > 0
energy_relative_error < target threshold for adiabatic case
```

Initial threshold recommendation:

```text
energy_relative_error < 5e-2 for first gate
then improve toward 1e-2 after stable
```

If using fixed controlled Joule field instead of EM-generated Joule, mark:

```text
thermal_mapping_only = 1
source_driven_em_thermal = 0
```

If using P1 EM-generated Joule, mark:

```text
source_driven_em_thermal = 1
```

### 6.6 Required outputs

Print:

```text
passed
thermal_mapping_only
source_driven_em_thermal
P_Joule
Delta_t
expected_energy
Delta_E_thermal
energy_relative_error
T_min_start, T_max_start, T_avg_start
T_min_end, T_max_end, T_avg_end
```

VTP fields:

```text
Temperature
JouleHeatSource
Sigma
MaterialRegionId
```

---

## 7. P3 — simple conductor heating validation

### 7.1 Purpose

After P1 and P2, create a clean validation case for electromagnetic heating in a conductor. This case is not TEAM7. It is a minimal solver-closure validation.

### 7.2 Recommended test name

Create:

```text
tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_simple_conductor_heating_validation/
```

### 7.3 Validation idea

Use a conductor block/slab with either:

```text
A. source-driven EM Joule field from P1
or
B. controlled near-uniform Joule field for thermal validation
```

The case should verify that the average temperature rise is consistent with the heat source:

```text
Delta_T_avg ≈ P_Joule * Delta_t / integral(rho cp dV)
```

### 7.4 P3 hard gates

Pass only if:

```text
all fields finite
Joule source finite and positive in conductor
temperature increases
energy_relative_error < 5e-2 initially
VTP output contains required fields
```

### 7.5 Documentation

In the record file, write:

```text
- Is Joule source EM-generated or controlled?
- What is the expected average temperature rise?
- What is the observed average temperature rise?
- Energy error
- Any boundary loss assumptions
```

---

## 8. P4 — simplified TEAM7 source-driven scaffold

### 8.1 Purpose

Only start this after P1 passes. Prefer starting after P2/P3 have at least thermal mapping validation.

Goal: create a simplified TEAM7-like source-driven electromagnetic case, not full standard TEAM7.

### 8.2 Recommended test name

Create:

```text
tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_simplified_team7_source_driven/
```

### 8.3 Geometry

Start with a simplified geometry:

```text
large air box, but not full-size TEAM7 yet
conductive plate / block
optional hole or slot after first version passes
coil/source region near plate
single resolution first
```

Do not start with full TEAM7 dimensions if it creates too many particles. First create a stable reduced problem.

### 8.4 Physics

Use:

```text
constant sigma in conductor
sigma = 0 or near zero in air
constant nu
frequency = 50 Hz first
optional 200 Hz later
source-driven A-phi solve
lambda_A = off
phi gauge on
```

### 8.5 Outputs

Required regionwise metrics:

```text
plate/conductor Joule integral
air Joule integral
max_B in air
max_B in conductor
max_J in conductor
max_Joule in conductor
GMRES iterations and residual
```

Required VTP fields:

```text
A, phi, B, H, E, J, Joule, Sigma, Nu, MaterialRegionId
```

### 8.6 P4 hard gates

Pass only if:

```text
GMRES converged
fields finite
conductor J nonzero
conductor Joule integral > 0
air Joule integral small
VTP written
eddy current distribution qualitatively concentrated in conductor
```

This is still not standard TEAM7 validation.

Mark explicitly:

```text
simplified_team7 = 1
standard_team7_validation = 0
```

---

## 9. P5 — probe and metric infrastructure for TEAM7

### 9.1 Purpose

Before standard TEAM7, implement reusable probe/metric infrastructure. TEAM7 validation depends on comparing field values at specified lines/points/regions.

### 9.2 Required metrics

Implement or standardize helpers for:

```text
line probe of |B|
line probe of B components
line probe of |J|
point probe at selected coordinates
region integral of Joule
region max/avg/min for |B|, |J|, Joule
GMRES convergence summary
particle count by region
```

### 9.3 Output format

For each test, write probe data to simple text/CSV files, for example:

```text
bin/output/aphi_probe_B_line.csv
bin/output/aphi_probe_J_line.csv
bin/output/aphi_region_metrics.csv
```

CSV columns should be explicit:

```text
x,y,z,Bx_real,By_real,Bz_real,B_abs,Jx_real,Jy_real,Jz_real,J_abs,Joule,Sigma,RegionId
```

For complex fields, either output real/imag columns or magnitude, but document which one is used.

### 9.4 P5 hard gates

Pass if:

```text
probe files written
number of probe samples > 0
all probe values finite
region metrics finite
```

---

## 10. P6 — TEAM7 large-air-domain and multi-resolution planning

### 10.1 Purpose

Standard TEAM7 likely has a very large air domain. A uniform SPH discretization may create too many particles. Multi-resolution will probably be needed, but it should not be mixed into the first source-driven scaffold.

### 10.2 Create planning document

Create:

```text
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_14_TEAM7_MULTIREsolution_PLAN.md
```

Use the exact filename style preferred by the repository. If lowercase is preferred, use:

```text
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_14_TEAM7_MULTIRESOLUTION_PLAN.md
```

### 10.3 The document must discuss

```text
1. Standard TEAM7 geometry and air-domain size issue
2. Why uniform resolution is expensive
3. Single-resolution first scaffold
4. Two-level resolution strategy:
   - conductor/coil fine
   - near air medium/fine
   - far air coarse
5. Transition layer / buffer requirements
6. Impact on pairwise Laplace consistency
7. Impact on corrected curl B
8. Impact on Contact relation if multibody
9. Diagnostics needed before trusting multi-resolution
10. Proposed particle-count estimates
```

### 10.4 Do not implement multi-resolution immediately

Only plan it in Stage 10.14 unless P1-P4 are already stable. Multi-resolution should become a later stage after single-resolution source-driven TEAM7 scaffold passes.

---

## 11. P7 — Contact roadmap after solver closure

Contact is not finished, but do not prioritize three-body A-penalty before source-driven heating closure.

After P1-P4 pass, schedule:

```text
1. three-body Contact baseline source-driven solve
2. air + conductor + wall/source body layout
3. E/J/Joule consistency across bodies
4. regionwise Joule integral consistency
5. interface spike diagnostics
6. optional three-body A-penalty research with InnerOnly stencil
```

Keep:

```text
lambda_A production off
projection deferred
```

Do not start P7 until P1 source-driven solve and P2/P3 thermal coupling validation are stable.

---

## 12. P8 — Maxwell-Ampere MMS positioning

The full Maxwell-Ampere MMS design from Stage 10.13 remains important, but it should not block P1-P3.

Purpose of Maxwell-Ampere MMS:

```text
verify B = curl A
verify H = nu B
verify curl H ≈ J with a manufactured field that is Maxwell-consistent
```

This is not a full Maxwell solver. It is a manufactured-solution validation for the current A-phi eddy-current formulation.

Recommended order:

```text
P1 source-driven solve first
P2/P3 Joule-to-thermal closure second
P4 simplified TEAM7 scaffold third
Maxwell-Ampere MMS can run in parallel only if it does not delay P1-P3
```

Do not make curlH≈J a hard gate using old div-free MMS fields. Only use a Maxwell-consistent MMS for that.

---

## 13. Suggested implementation order

Follow this order exactly unless a blocking compile error forces local fixes:

```text
Step 1: Create Stage 10.14 record file.
Step 2: Implement P1 source-driven A-phi solve case.
Step 3: Run P1 and document solver metrics.
Step 4: Implement P2 JouleHeatSource -> thermal source mapping.
Step 5: Run P2 with controlled Joule or P1-generated Joule.
Step 6: Implement P3 simple conductor heating validation and energy accounting.
Step 7: Implement P5 probe/metric infrastructure if needed for P4.
Step 8: Implement P4 simplified TEAM7 source-driven scaffold, single-resolution.
Step 9: Write P6 TEAM7 multi-resolution planning document.
Step 10: Only after these pass, revisit Contact three-body and Maxwell MMS scheduling.
```

---

## 14. Suggested pass/fail summary table for record file

Add a table like this to the Stage 10.14 record:

```markdown
| Stage | Test | Purpose | Status | Key metrics | Notes |
|---|---|---|---|---|---|
| P1 | test_3d_aphi_ck_source_driven_em_solve | source-driven A-phi solve | TODO | GMRES residual, Joule integral | lambda_A off |
| P2 | test_3d_aphi_ck_em_joule_thermal_one_way | Joule -> thermal source | TODO | energy_relative_error | one-way coupling |
| P3 | test_3d_aphi_ck_simple_conductor_heating_validation | heating energy validation | TODO | Delta_T_avg, energy error | adiabatic preferred |
| P4 | test_3d_aphi_ck_simplified_team7_source_driven | simplified TEAM7 scaffold | TODO | conductor Joule, probes | not standard TEAM7 |
| P5 | probe helpers | TEAM7 probe infrastructure | TODO | finite CSV outputs | reusable |
| P6 | multi-resolution plan | TEAM7 air-domain strategy | TODO | particle estimates | document only |
```

---

## 15. Standard TEAM7 long-term roadmap

After Stage 10.14, standard TEAM7 should proceed as a separate stage.

### Required before standard TEAM7

```text
source-driven solve passes
Joule-to-thermal mapping passes
simple conductor heating energy validation passes
simplified TEAM7 source-driven scaffold passes
probe infrastructure exists
large-air-domain strategy exists
```

### Standard TEAM7 stage should include

```text
1. exact geometry definition
2. exact material parameters
3. exact frequency cases, likely 50 Hz and 200 Hz if following standard references
4. coil/source definition
5. large air domain or equivalent truncation strategy
6. single-resolution baseline if affordable
7. multi-resolution version if necessary
8. probe locations and reference comparisons
9. dp/refinement convergence
10. B/J/Joule absolute and relative error metrics
```

Do not claim standard TEAM7 validation until field probes are compared with reference data.

---

## 16. Code style requirements

Maintain SPHinXsys SYCL CK style. Do not introduce a separate coding style.

Expected style:

```text
CK UpdateKernel / InteractionDynamicsCK style where applicable
DelegatedData / variable registration pattern
existing A-phi helper conventions
clear test-level helpers in diagnostics or benchmark subdirectories
minimal duplication
CMake integration consistent with other 3d_examples
```

Avoid:

```text
large monolithic test cpp files with hidden logic
host-only code paths that bypass the existing CK data structure
new global state
hard-coded magic values without comments
changing existing defaults silently
```

Every new helper should document:

```text
what it computes
which variables it reads/writes
whether it is source-driven or impressed
whether it is validation or diagnostic/demo
```

---

## 17. Required final response from Cursor after implementation

After completing each subtask, Cursor should report in this exact structure:

```text
Completed stage:
Files added:
Files modified:
Tests run:
Pass/fail:
Key printed metrics:
VTP/CSV outputs:
Known limitations:
Next recommended step:
```

For Stage 10.14 closure, Cursor should provide:

```text
- final Stage 10.14 record file
- list of all new/modified source files
- list of all new test cases
- terminal summary for each test
- whether P1/P2/P3/P4/P5/P6 passed or remains partial
- explicit statement that projection and production lambda_A were not changed
```

---

## 18. Final priority reminder

The priority is now:

```text
1. source-driven A-phi solve
2. JouleHeatSource -> thermal source coupling
3. energy validation of electromagnetic heating
4. simplified TEAM7 source-driven scaffold
5. TEAM7 probe infrastructure
6. TEAM7 large-air-domain / multi-resolution planning
7. then multi-body Contact production validation
8. then Contact A-penalty three-body research / projection decision
```

Do not spend time polishing the customer demo unless the source-driven solver closure is blocked or explicitly requested.
