# TEAM7 Reference Data and Native Geometry Plan for Cursor

> Project: SPHinXsys SYCL CK A-phi electromagnetic heating solver  
> Stage target: Native TEAM7 geometry + B-based probe/reference validation  
> Prepared: 2026-06-02  
> Scope: TEAM7 only. Do not start cold-crucible cases in this stage.

---

## 0. Executive summary

Next stage should be:

```text
Stage 10.17 — Native TEAM7 geometry and Bz reference pipeline
```

Immediate target:

```text
SPHinXsys native TEAM7 geometry / XML / reload particles
+ truncated air domain
+ source-driven A-phi solve
+ Bz probe CSV
+ comparison against TEAM7 reference data
```

Do **not** implement a new hand-made particle generator. Use SPHinXsys native geometry, level set, particle generation, relaxation, and reload pipeline as specified by the project owner.

---

## 1. Reference hierarchy

Use the following hierarchy.

### 1.1 Primary benchmark definition

Official TEAM Problem 7 PDF:

```text
https://www.compumag.org/jsite/images/stories/TEAM/problem7.pdf
```

Use it as the official definition for the problem, geometry concept, material, excitation, and validation targets.

### 1.2 Implementation / probe-location guide

COMSOL TEAM7 model documentation:

```text
https://doc.comsol.com/6.4/doc/com.comsol.help.models.acdc.multiturn_coil_asymmetric_conductor/multiturn_coil_asymmetric_conductor.html
https://www.comsol.com/model/download/1175261/models.acdc.multiturn_coil_asymmetric_conductor.pdf
```

Use COMSOL for:

```text
coil turns/current
frequencies
material properties
Bz probe line
Jy probe line
plot expressions
```

### 1.3 Machine-readable Bz reference

muFEM open implementation:

```text
https://github.com/Raiden-Numerics/mufem/tree/main/Electromagnetics/Compumag-Team7-Asymmetrical-Conductor-with-a-Hole
```

This provides machine-readable Bz reference files for:

```text
A1-B1: y = 72 mm, z = 34 mm
A2-B2: y = 144 mm, z = 34 mm
```

The Bz CSVs are included in this bundle.

---

## 2. Physical/model parameters to use

| Item | Value / instruction |
|---|---:|
| Problem | TEAM Problem 7: asymmetrical conductor with a hole |
| Geometry | Thick aluminum plate with eccentric hole below a multi-turn coil |
| Symmetry | No symmetry; solve full 3D geometry |
| Coil turns | 2742 |
| Coil current | 1 A/turn |
| Ampere-turns | 2742 AT |
| Frequencies | 50 Hz and 200 Hz |
| Aluminum conductivity | 3.526e7 S/m |
| Aluminum relative permeability | 1 |
| Plate/conductor | eddy-current domain |
| Coil/source | prescribed-current source region for current solver stage |
| Air | non-eddy-current region; air Joule is not a physical heat target |

Recommended project policy:

```text
Air:
    sigma = 0, or tiny numerical sigma only if required for conditioning.
    Do not count air Joule as physical heating.

Plate:
    sigma = 3.526e7 S/m.
    mu_r = 1.
    Eddy current enabled.
    Conductor Joule is the physical heating target.

Coil/source:
    prescribed current source.
    2742 turns * 1 A/turn equivalent.
    Source-only region for current solver stage.
    Do not model real helical copper wire.
    Do not count coil/source Joule as physical heating target.
```

---

## 3. Probe definitions

Use meters internally, but write CSV outputs in mm for TEAM7 comparison.

### 3.1 COMSOL Bz probe line

```text
quantity: Bz magnetic flux density
x range: 0 <= x <= 288 mm
y = 72 mm
z = 34 mm
frequencies: 50 Hz, 200 Hz
unit: mT
COMSOL expression: sign(real(mf.Bz))*abs(mf.Bz)
```

This line is between the aluminum plate and the coil.

### 3.2 COMSOL Jy probe line

```text
quantity: Jy induced current density
x range: 0 <= x <= 288 mm
y = 72 mm
z = 18.99 mm, approximately z = 19 mm
frequencies: 50 Hz, 200 Hz
unit: A/m^2
COMSOL expression: sign(real(mf.Jy))*abs(mf.Jy)
```

This line is on or slightly below the conductor surface and crosses the gap.

The numerical Jy table was not directly scraped from the webpage. Cursor must not invent Jy reference values. Treat Jy as diagnostic until `multiturn_coil_asymmetric_conductor_table2.txt` or digitized data is supplied.

### 3.3 muFEM Bz reference lines

```text
A1-B1:
    x = 0..288 mm
    y = 72 mm
    z = 34 mm

A2-B2:
    x = 0..288 mm
    y = 144 mm
    z = 34 mm
```

The included reference CSVs contain 17 points at x = 0, 18, ..., 288 mm.

---

## 4. Included CSV files

Put these files in the repository under:

```text
tests/extra_source_and_tests/3d_examples/reference_data/team7/
```

Files:

```text
TEAM7_Bz_A1_B1_reference_Gauss.csv
TEAM7_Bz_A2_B2_reference_Gauss.csv
TEAM7_Bz_A1_B1_reference_mT.csv
TEAM7_Bz_A2_B2_reference_mT.csv
TEAM7_probe_definitions.csv
```

The public raw data are in Gauss. The mT files use:

```text
1 Gauss = 0.1 mT
```

---

## 5. Complex/sign convention warning

muFEM reconstructs the reported complex Bz using:

```text
Bz_complex_mT = 1e3 * (-B_real.z + i * B_imag.z)
```

Therefore Cursor must not blindly compare our `BReal.z` directly against the reference. First output all four sign options:

```text
(+real, +imag)
(-real, +imag)
(+real, -imag)
(-real, -imag)
```

Then compute profile errors and record which convention matches the reference.

Create a sign-convention audit helper and keep the result in the test record.

---

## 6. Native TEAM7 geometry validation

The user will provide or configure SPHinXsys-native TEAM7 geometry / XML / reload setup.

Cursor must not write a new geometry generator.

Cursor must only:

```text
read the XML / reload particles;
build AirBody, PlateBody, Coil/SourceBody;
assign sigma/nu/source;
connect to existing A-phi source-driven solve;
write VTP and probe CSV.
```

### 6.1 Geometry audit required

Create:

```text
test_3d_aphi_ck_team7_native_geometry_audit
```

It must print:

```text
air particle count
plate particle count
source particle count
air bounding box
plate bounding box
source bounding box
plate thickness
hole center if available
hole dimensions if available
source/coil center
source/coil bounding box
air domain bounding box
probe A1-B1 start/end
probe A2-B2 start/end
probe Jy start/end
VTP path
```

Pass gate:

```text
air particles > 0
plate particles > 0
source particles > 0
finite bounding boxes
probe lines inside truncated air domain or valid sampling region
VTP written
```

No A-phi solve in this audit test.

---

## 7. Native TEAM7 source-driven smoke

Create:

```text
test_3d_aphi_ck_team7_native_geometry_source_driven_smoke
```

Purpose:

```text
run source-driven A-phi on native low-resolution TEAM7 geometry.
```

Initial pass gate:

```text
GMRES converged or diagnostic-converged on first low-resolution run
fields finite
source_rhs_l2 > 0
plate/conductor Joule > 0
air/source Joule is not used as physical heating target
VTP written
```

Do not require strict Bz accuracy yet.

---

## 8. Native TEAM7 Bz reference probe

Create:

```text
test_3d_aphi_ck_team7_native_geometry_bz_reference_probe
```

Outputs:

```text
probe_A1_B1_Bz.csv
probe_A2_B2_Bz.csv
comparison_A1_B1_vs_reference.csv
comparison_A2_B2_vs_reference.csv
summary.csv
```

Required output columns:

```text
x_mm
our_Bz_real_mT
our_Bz_imag_mT
our_minus_Bz_real_mT
our_minus_Bz_imag_mT
ref_Bz_50Hz_phase0_mT
ref_Bz_50Hz_phase90_mT
ref_Bz_200Hz_phase0_mT
ref_Bz_200Hz_phase90_mT
```

First implementation may compare only 50 Hz. Add 200 Hz once 50 Hz is stable.

---

## 9. Initial tolerance strategy

Do not set a strict numerical gate immediately.

### Stage A: smoke/probe correctness

```text
CSV written
all probe values finite
probe coordinates correct
source_rhs_l2 > 0
conductor Joule > 0
```

### Stage B: loose reference gate

After sign convention is established:

```text
A1-B1 Bz real relative profile error < 100%
A1-B1 Bz imag relative profile error < 100%
A2-B2 remains informative initially
```

### Stage C: medium reference gate

After geometry and source scaling are stable:

```text
A1-B1 Bz real relative profile error < 50%
A1-B1 Bz imag relative profile error < 50%
```

### Stage D: validation gate

After dp/domain/source calibration:

```text
Bz profile error < 20% or better
```

Do not start with 5-10% tolerance.

---

## 10. Jy reference status

COMSOL defines the Jy probe, but this bundle does not include machine-readable Jy values.

Expected COMSOL application file:

```text
multiturn_coil_asymmetric_conductor_table2.txt
```

Recommended path if COMSOL is installed:

```text
<COMSOL>/Multiphysics/applications/ACDC_Module/Verifications/
    multiturn_coil_asymmetric_conductor_table2.txt
```

Until that file is available:

```text
Jy probe is diagnostic-only.
Do not create fake Jy reference data.
```

---

## 11. Do not do in this stage

Do not implement:

```text
cold crucible case
multiresolution
projection
production lambda_A
ghost/mirror EM boundary
real helical coil
coil self-heating
circuit coupling
strict Bz < 10% gate
hand-written custom geometry/particle generator
invented Jy reference data
```

---

## 12. Stage 10.17 closure criteria — **CLOSED 2026-06-03**

```text
1. TEAM7 reference CSVs under reference_data/team7 — done
2. native geometry audit — passed
3. native source-driven smoke — passed
4. A1-B1 / A2-B2 Bz CSV + comparison — done
5. Bz vs reference — peak ~8 mT, profile L2 ~20%, sign (1,-1) frozen
6. P1 source audit NI=2742 A — passed
7. Jy — diagnostic only (no fake reference)
8. SYCL CK + native geometry pipeline — done
```

See `CURSOR_APHI_STAGE10_17_TEAM7_NATIVE_RECORD.md`. Next: **Stage 10.18** Joule post (`CURSOR_APHI_STAGE10_18_NATIVE_JOULE_POST_PLAN.md`).

---

## 13. Final instruction for Cursor

The immediate validation target is:

```text
Bz along A1-B1 and A2-B2 compared against reference data.
```

The immediate geometry target is:

```text
native SPHinXsys TEAM7 geometry, particle generation, relaxation/reload, and truncated air-domain pipeline.
```

Do not continue with hand-coded box/annulus masks as the main geometry route.
