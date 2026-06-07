# CURSOR_APHI_STAGE10_14_CONTACT_BOUNDARY_ADDENDUM

> Purpose: add an explicit Contact-relation and boundary-treatment roadmap to Stage 10.14.  
> This addendum supplements `CURSOR_APHI_STAGE10_14_REAL_EM_HEATING_SOLVER_PLAN.md`.  
> Follow this document literally. It exists because the previous Stage 10.14 plan mentioned Contact/boundary only at roadmap level and did not give enough implementation detail.

---

## 0. Executive decision

Stage 10.14 remains focused on building a usable electromagnetic heating solver:

```text
source-driven A-phi solve
-> E/J/Joule postprocess
-> Joule-to-thermal coupling
-> energy validation
-> simplified TEAM7 scaffold
```

However, **Contact and boundary handling must now be treated as explicit solver-closure work**, not as vague future items.

The corrected Stage 10.14 priority order is:

```text
P0  Preserve current numerical freeze decisions
P1  Source-driven single-body region-tagged solve
P2  Joule-to-thermal one-way heating closure
P3  Simple conductor heating energy validation
P4  Contact baseline production validation, lambda_A off
P5  Boundary/support policy diagnostics
P6  Simplified TEAM7 source-driven scaffold, single resolution
P7  TEAM7 probe/metric infrastructure
P8  TEAM7 large-air-domain / multi-resolution plan
P9  Three-body Contact A-penalty research gate, optional and after P4-P6
```

Do not skip P4/P5 before claiming the solver is suitable for TEAM7-like multibody cases.

---

## 1. Current Contact status

The current Contact status is:

```text
Contact main A-phi operator: available
Two-body Contact baseline: mostly stable
Two-body Contact MMS: passed in previous stages
Variable sigma interface MMS: passed in previous stages
Contact A-divergence penalty: two-body research gate closed using InnerOnly penalty stencil
Contact graddiv PC: InnerOnly skip of Contact graddiv has been verified diagnostically
Three-body Contact production validation: not closed
Three-body Contact A-penalty: deferred
Boundary/interface divA: open diagnostic
Projection: design only, not implemented
```

Important frozen decisions:

```text
production lambda_A = off
production eta_A = 0
Contact main operator K(A,phi) = Inner + Contact
Contact A-divergence penalty stencil = InnerOnly
eta_A = 0.1 only as research primary
eta_A >= 0.2 requires polish_sweeps = 0
projection = defer
```

Do not reinterpret two-body Contact success as complete Contact closure.

---

## 2. What Contact must support for a complete electromagnetic heating solver

A usable electromagnetic heating solver must eventually support at least these body layouts:

### 2.1 Single-body region-tagged layout

Used for early source-driven closure and simple validation.

```text
One SPH body
region tags: air / conductor / coil-source / wall-dummy if needed
material variables assigned by region
no Contact relation needed inside the case
```

Purpose:

```text
fastest path to source-driven A-phi solve and Joule heating closure
lowest debugging complexity
best first target for energy validation
```

This is acceptable for P1-P3.

### 2.2 Two-body Contact layout

Used for interface validation.

Example:

```text
body 1: air or weakly conductive exterior
body 2: conductive plate / melt / solid conductor
relations:
  Inner<body1>
  Inner<body2>
  Contact<body1 -> body2>
  Contact<body2 -> body1>
```

Purpose:

```text
verify that multi-body neighbor search, interface coupling, material jumps, E/J/Joule and power integral remain stable
```

### 2.3 Three-body Contact layout

Used before TEAM7/cold-crucible multibody production.

Example:

```text
body 1: air domain
body 2: conductive plate or melt
body 3: wall / crucible / source-support body
optional region: coil source represented as impressed current/RHS rather than real conductor initially
```

Purpose:

```text
prove production Contact baseline works beyond two-body MMS
```

This must be done with:

```text
lambda_A = off
eta_A = 0
phi gauge penalty on if needed
main operator = Inner + Contact
no projection
```

Only after this passes should Contact A-penalty be reconsidered for three-body research.

---

## 3. P4 — Contact baseline production validation

Add this as a required Stage 10.14 task before simplified TEAM7 is called production-relevant.

### 3.1 Goal

Validate that the source-driven A-phi solver works in multibody Contact mode with realistic material jumps, while keeping A-divergence penalty off.

### 3.2 New test case

Create:

```text
tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_contact_source_driven_heating_baseline/
```

Suggested files:

```text
test_3d_aphi_ck_contact_source_driven_heating_baseline.cpp
README.md
```

Suggested helper:

```text
tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_contact_source_driven_heating_baseline_helpers.h
```

### 3.3 Geometry

Start with a reduced two-body setup:

```text
body A: air box
body B: conductive block or plate embedded inside / adjacent to air
source: analytic coil-like RHS or compact impressed current region
```

Do not start with full TEAM7 geometry.

### 3.4 Relations

Use real Contact relations, not same-body region tags:

```cpp
Inner<> air_inner(air_body);
Inner<> conductor_inner(conductor_body);
Contact<> air_contact(air_body, {&conductor_body});
Contact<> conductor_contact(conductor_body, {&air_body});
```

If the actual SPHinXsys syntax differs, use the local project convention, but the structure must be two distinct bodies and bidirectional Contact.

### 3.5 Numerical settings

Use:

```text
lambda_A = off
eta_A = 0
Contact A-penalty disabled
phi gauge penalty = current stable default, e.g. 100 if already used in Contact solve
GMRES tolerance = start with 5e-4, then try 1e-4 if stable
polish_sweeps = current stable production default, but do not use polish to hide non-convergence
```

### 3.6 Required outputs

The test must print and record:

```text
passed
GMRES converged flag
GMRES iteration count
final relative residual
max |A|
max |phi|
max |B|
max |E|
max |J|
max Joule
air Joule integral
conductor Joule integral
ratio conductor_Joule / max(air_Joule, eps)
finite-field flag
VTP output path
```

### 3.7 Hard gate

The test passes only if:

```text
GMRES converged = true
all fields finite
max |A| > small_threshold
max |E| > small_threshold
max |J| > small_threshold in conductor
conductor Joule integral > positive_threshold
air Joule integral is near zero if sigma_air = 0, or much smaller than conductor if sigma_air is tiny
no NaN or Inf in any body
VTP exists and contains body-wise/region-wise fields
```

### 3.8 Diagnostic but not hard gate

Report but do not fail on:

```text
global divA
boundary divA
interface divA
B error against analytic field, if no analytic source-driven reference exists
```

Do not turn divA into a hard gate in this Contact source-driven baseline case.

---

## 4. P4b — Three-body Contact baseline validation

After P4 two-body passes, add a three-body Contact baseline.

### 4.1 New test case

Create:

```text
tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_three_body_contact_source_driven_baseline/
```

### 4.2 Geometry

Use a small reduced domain, not TEAM7 size:

```text
body 1: air
body 2: conductive plate / melt block
body 3: wall or non-conductive support shell
source: compact analytic coil/source RHS
```

### 4.3 Required Contact graph

Use bidirectional contacts for all physically adjacent pairs:

```text
air <-> conductor
air <-> wall/support
conductor <-> wall/support, only if geometrically adjacent
```

Do not connect bodies that are not within kernel support unless the case geometry requires it.

### 4.4 Settings

Use:

```text
lambda_A = off
eta_A = 0
A-divergence penalty disabled
projection disabled
phi gauge penalty as needed
```

### 4.5 Hard gate

Same as P4, plus:

```text
all three bodies have finite fields
non-conductive wall/support Joule integral is near zero
conductor Joule dominates
Contact neighbor counts are nonzero at expected interfaces
no abnormal interface spike in E/J/Joule compared with conductor interior maximum
```

Interface spike check suggestion:

```text
interface_J_max <= C_spike * conductor_core_J_max
interface_Joule_max <= C_spike * conductor_core_Joule_max
```

Start with:

```text
C_spike = 10
```

Do not make this too strict at first; use it to detect catastrophic interface blow-up.

---

## 5. Boundary/support handling policy

Boundary handling must be explicit. Do not hide boundary issues behind projection.

### 5.1 Current boundary interpretation

Current observations indicate:

```text
core divA can be very small
boundary/global divA can be large due to support deficiency and denominator artifacts
B=curlA is dp-convergent but boundary-sensitive
E/J/Joule are more stable than global divA
```

Therefore:

```text
global divA is not a production hard gate yet
core divA and observables are more meaningful
boundary treatment must be studied with controlled diagnostics
```

### 5.2 Boundary strategies to compare

Implement diagnostics for the following boundary policies, in this order.

#### Strategy A — baseline truncated domain

```text
No extra ghost/dummy particles
regular air/conductor domain only
```

Purpose:

```text
baseline for all comparisons
```

#### Strategy B — enlarged air buffer

```text
same physical conductor/source
larger surrounding air domain
same particle spacing
```

Purpose:

```text
check sensitivity to far-field truncation
estimate how large air domain must be before conductor observables stabilize
```

#### Strategy C — passive dummy support shell

```text
add non-physical dummy/support shell outside the computational air domain
used only to improve kernel support near outer boundary
sigma = 0 or same as air depending on implementation
not included in physical Joule integral
```

Purpose:

```text
reduce boundary support deficiency without pretending projection solved gauge
```

Do not use dummy shell as a physical conductor.

#### Strategy D — ghost / mirror particles

Only implement after Strategy B/C results are understood.

```text
explicit ghost/mirror particles for boundary support
must not contribute to physical material integrals unless intentionally designed
```

Purpose:

```text
research option for support completion
```

Not production until separately validated.

### 5.3 New boundary diagnostic test

Create:

```text
tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_boundary_support_policy_diagnostic/
```

Suggested helper:

```text
tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_boundary_support_policy_diagnostic_helpers.h
```

### 5.4 Required command-line or compile-time modes

Support at least:

```text
--boundary_policy baseline
--boundary_policy enlarged_air
--boundary_policy dummy_shell
```

If command-line parsing is not convenient, implement separate diagnostic functions and print the policy name clearly.

### 5.5 Metrics to print

For each policy, print:

```text
particle counts per region/body
average neighbor count in core
average neighbor count near outer boundary
min neighbor count near outer boundary
GMRES converged flag, if solving
final residual, if solving
max |A|
max |B|
max |E|
max |J|
conductor Joule integral
air Joule integral
core divA rel
interface divA rel
outer-boundary divA rel
global divA rel
B smoothness indicator near conductor
VTP path
```

### 5.6 Hard gate

This diagnostic should not have a strict physics pass/fail initially. It should pass if:

```text
all policies run without NaN/Inf
metrics are printed
VTP files are generated
conductor Joule remains finite and positive
```

### 5.7 Decision criteria after diagnostic

Use the diagnostic to decide:

```text
If enlarged air strongly changes conductor Joule, current air domain is too small.
If enlarged air stabilizes conductor Joule but costs too many particles, plan multi-resolution.
If dummy shell improves boundary divA/B smoothness without changing conductor Joule materially, keep it as a candidate support strategy.
If dummy shell changes conductor Joule significantly, do not use it in production until understood.
If ghost/mirror improves divA but damages E/J/Joule, do not use it.
```

Suggested quantitative thresholds for first pass:

```text
relative change in conductor Joule integral < 2% between boundary policies: acceptable for scaffold
relative change in max conductor J < 5%: acceptable for scaffold
no NaN/Inf: mandatory
```

Do not claim final TEAM7 boundary validation from these thresholds.

---

## 6. Interface diagnostics for Contact and boundary

Add reusable diagnostics to avoid hiding interface problems.

### 6.1 Regions

For any multibody/contact case, classify particles into:

```text
core particles: far from body/contact boundary
interface particles: within one kernel support of another body
outer-boundary particles: near domain truncation boundary
source-near particles: near coil/source region
```

If exact distance-to-boundary is unavailable, approximate with neighbor-count deficiency and contact-neighbor count.

### 6.2 Metrics

For each class, report:

```text
count
avg neighbor count
min neighbor count
max |A|
max |B|
max |E|
max |J|
max Joule
avg Joule
Joule integral
relative divA if available
```

### 6.3 Required warning conditions

Print warnings if:

```text
interface max J > 10 * core max J
interface max Joule > 10 * core max Joule
outer-boundary max B > 10 * conductor-near max B without physical reason
any class has min neighbor count too low
non-conductive body has large Joule integral
```

Warnings do not necessarily fail the test, but they must be visible in terminal output and record md files.

---

## 7. Contact A-divergence penalty policy

Do not use A-divergence penalty to solve production boundary/contact issues yet.

### 7.1 Production

Production remains:

```text
lambda_A = off
eta_A = 0
```

### 7.2 Research

Only after P4/P4b and P5 pass, create optional research tests:

```text
test_3d_aphi_ck_three_body_contact_a_penalty_research/
```

Settings:

```text
Contact main operator = Inner + Contact
A-divergence penalty stencil = InnerOnly
eta_A = 0.1 first
eta_A = 0.2 optional
polish_sweeps = 0 for eta_A >= 0.2
```

Hard gate should focus on:

```text
GMRES convergence
E/J/Joule not degraded relative to eta=0
no interface spike
core divA improvement
```

Do not require global divA to become small.

### 7.3 Explicitly forbidden

Do not re-enable:

```text
Contact A-divergence penalty stencil = InnerContact
```

unless the test name includes `diagnostic` and the expected result is to show it fails or pollutes the interface.

---

## 8. Projection policy

Projection remains deferred.

Do not implement:

```text
A_new = A - grad(chi)
phi_new = phi + i omega chi
```

until all of the following are true:

```text
source-driven A-phi solve passes
Joule-to-thermal coupling passes energy validation
Contact two-body and three-body baseline source-driven cases pass
boundary policy diagnostic identifies a stable support strategy
simplified TEAM7 scaffold passes
```

Even then, projection must be a separate research stage and must preserve E/J/Joule invariance. It cannot be introduced silently into production.

---

## 9. TEAM7 boundary and large air domain

TEAM7 introduces a large air domain. Do not solve it by immediately using full uniform fine resolution.

### 9.1 First TEAM7-like scaffold

Use single resolution and reduced air domain first:

```text
reduced air box
conductive plate with simplified shape
source coil/RHS
lambda_A off
projection off
```

Goal:

```text
prove source-driven solve and Joule output work in a TEAM7-like geometry
```

### 9.2 Air-domain sensitivity sweep

Before multi-resolution, perform an air-box size sweep:

```text
air_scale = 1.0, 1.5, 2.0, maybe 3.0 if particle count allows
same conductor resolution
same source
```

Metrics:

```text
GMRES iterations
final residual
conductor Joule integral
max conductor J
probe B values near conductor
outer boundary divA/B indicators
particle count
runtime
memory if available
```

Decision:

```text
If conductor Joule/probe values keep changing strongly with air_scale, boundary is not far enough.
If large air_scale is too expensive, multi-resolution is required.
```

### 9.3 Multi-resolution planning

Do not implement multi-resolution before the single-resolution scaffold and air-size sweep are understood.

When planning multi-resolution, address:

```text
fine conductor region
fine coil/source region
transition buffer
coarse far air region
neighbor relation across resolution transition
operator consistency across resolution jump
Contact relation if bodies have different dp
B=curlA reconstruction across resolution transition
Joule integral consistency
GMRES conditioning
```

---

## 10. Updated Stage 10.14 acceptance criteria

Stage 10.14 should not be declared closed until these are done:

```text
P1 source-driven single-body solve passed
P2 Joule-to-thermal one-way coupling passed
P3 simple conductor heating energy validation passed
P4 two-body Contact source-driven baseline passed
P4b three-body Contact source-driven baseline passed or explicitly deferred with reason
P5 boundary/support policy diagnostic completed
P6 simplified TEAM7 source-driven scaffold passed
P7 probe/metric infrastructure available
P8 TEAM7 air-domain/multi-resolution plan written
```

If time is limited, minimum closure requires:

```text
P1 + P2 + P3 + P4 + P5
```

Do not claim readiness for TEAM7 without P5/P6/P7.

---

## 11. Recommended record updates

Update:

```text
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_14_REAL_EM_HEATING_SOLVER_RECORD.md
```

Add these sections:

```text
12. Contact baseline production validation
13. Boundary/support policy diagnostics
14. Interface spike diagnostics
15. TEAM7 air-domain implications
16. Contact/projection decisions after Stage 10.14
```

Each completed test must record:

```text
test name
build/run command
passed flag
key terminal output
VTP output path
interpretation
remaining open issues
```

---

## 12. One-sentence instruction for Cursor

Use this exact instruction if needed:

```text
Do not treat Contact and boundary as solved. Add explicit Contact source-driven baseline tests and boundary/support policy diagnostics before claiming the electromagnetic heating solver is TEAM7-ready. Keep lambda_A off, keep projection deferred, and use InnerOnly A-penalty only for optional research diagnostics after baseline Contact and boundary tests pass.
```
