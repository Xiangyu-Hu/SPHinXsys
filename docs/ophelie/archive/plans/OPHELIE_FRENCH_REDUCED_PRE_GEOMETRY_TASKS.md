# OPHELIE-like French Cold-Crucible Direction: Development Task Manifest Before French Geometry Is Ready

This document guides Cursor to continue the OPHELIE-like / Biot-Savart particle-integral EM heating module before the formal French cold-crucible geometry body is ready. Conclusion: **TEAM7 is no longer the main validation target; the French cold-crucible reduced case is the customer delivery main line.** TEAM7 remains as benchmark / source smoke / total-field reference, and should no longer block the French main line.

---

## 0. Current Route Decision

### 0.1 Current main line

The main line is now:

```text
OPHELIE-like French cold-crucible reduced case
```

That is:

```text
external coil prescribed source
    -> Biot-Savart gives A_src / B_src on glass/probes
        -> conductive glass internal phi/V correction
            -> E, J, JouleHeat
                -> JouleHeat as input to later heat conduction/flow
```

Still out of scope:

```text
full PDE A-phi solver
air-domain far-field boundary
Coulomb gauge / divA projection
contact A-penalty
cold-crucible metal wall EM solve
real skin-depth surface model
sigma(T) feedback
thermal-EM two-way iteration
water-cooled wall thermal boundary
real helical coil / lead wires / impedance solve
```

### 0.2 TEAM7 current role

TEAM7 is not the main validation for the French plan. Existing results show:

```text
volume-etheta / circle loop / rect-loop / filament-racetrack sweep
cannot stably move Bz peak from x≈198 mm to reference x≈126 mm at reasonable RMS.
```

This means:

```text
1. no-phi / no-self coil-only B cannot strictly validate against TEAM7 reference;
2. TEAM7 reference is likely total B = coil B + plate-induced B;
3. simple path fitting / shell fitting is not the main line;
4. if TEAM7 continues, start with total-field probe, not racetrack path tuning.
```

Therefore TEAM7 is kept only as:

```text
1. completed source-field smoke archive;
2. future total-field benchmark;
3. future A-phi main-line source model reference.
```

---

## 1. Overall Goals Before French Geometry Body Is Ready

Before formal French geometry is ready, Cursor should complete infrastructure so the reduced case runs as soon as geometry arrives:

```text
A. Clean and stabilize current OPHELIE module;
B. Solidify Biot-Savart kernel, source model, power scaling, VTP output;
C. Build analytic/semi-analytic operator tests to avoid finding low-level errors via TEAM7;
D. Prepare parametric French reduced cylinder case without real French STL;
E. Prepare minimal reliable phi correction acceptance;
F. Prepare JouleHeat -> thermal interface fields;
G. Archive current TEAM7 results; stop path fitting.
```

---

## 2. Code Cleanup To Do Immediately

### 2.1 Add formal code to git

If OPHELIE files are still `??`, add to git soon:

```bash
git add tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie
git add tests/extra_source_and_tests/3d_examples/test_3d_ophelie
git add tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7
git add docs/ophelie
git status
```

Recommended initial commit:

```bash
git commit -m "Add OPHELIE-like particle induction module and TEAM7 smoke tests"
```

### 2.2 Clean reference bundles and obsolete builds

Reference directories can be deleted after confirming docs moved to `docs/ophelie/`:

```bash
mkdir -p docs/ophelie
cp FRENCH_OPHELIE_PARTICLE_DEVELOPMENT_PLAN.md docs/ophelie/  # if still at repo root
```

Can delete:

```bash
rm -rf biot_reference_selection_for_cursor
rm -rf french_ophelie_particle_plan_and_starter
rm -rf build_ophelie
rm -rf build_ophelie_cpu
```

Keep in use:

```bash
build/
```

---

## 3. OPHELIE Core Module File Structure To Stabilize

Recommended layout:

```text
tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/
    electromagnetic_ophelie.h
    electromagnetic_ophelie_parameters.h
    electromagnetic_ophelie_field_names.h
    electromagnetic_ophelie_register_fields.h

    electromagnetic_ophelie_source.h/.hpp
    electromagnetic_ophelie_biot_savart.h/.hpp
    electromagnetic_ophelie_postprocess.h/.hpp

    electromagnetic_ophelie_laplace.h
    electromagnetic_ophelie_phi.h/.hpp
    electromagnetic_ophelie_phi_gmres.h

    electromagnetic_ophelie_observables.h
    electromagnetic_ophelie_device_sync.h
    electromagnetic_ophelie_cli.h

    electromagnetic_ophelie_racetrack_source.h
    electromagnetic_ophelie_team7_probe.h
    electromagnetic_ophelie_team7_native_geometry.h

    # later add:
    electromagnetic_ophelie_french_reduced_geometry.h
    electromagnetic_ophelie_multiloop_source.h
    electromagnetic_ophelie_probe_sampling.h
    electromagnetic_ophelie_current_moment_source.h
```

---

## 4. Source Model Must Unify as Current Moment

### 4.1 Why abstract

Three source types may coexist:

```text
volume source:
    J [A/m²] * dV [m³] = A·m

surface source:
    K [A/m] * dS [m²] = A·m

filament source:
    I [A] * dl [m] = A·m
```

Biot-Savart needs:

```text
current moment = A·m
```

Recommended unify Biot-Savart source items as:

```cpp
Vecd CurrentMomentReal;
Vecd CurrentMomentImag;
```

Biot-Savart only:

```cpp
A += coeff * CurrentMoment * inv_r;
B += coeff * cross(CurrentMoment, r) * inv_r3;
```

This avoids volume / surface / filament mixing `VolumetricMeasure` unit errors.

### 4.2 Interface for Cursor

Add or consolidate:

```cpp
struct OphelieCurrentMomentSource
{
    Vecd position_;
    Vecd moment_real_;  // A*m
    Vecd moment_imag_;  // A*m
};
```

Or register on particle body:

```text
CurrentMomentReal : Vecd
CurrentMomentImag : Vecd
```

Per source model:

```text
volume:
    CurrentMomentReal = JsrcReal * Vol

filament:
    CurrentMomentReal = Ieff * dl

surface:
    CurrentMomentReal = KsrcReal * Area
```

First stage at least:

```text
volume
filament
```

Surface deferred.

---

## 5. Biot-Savart Analytic/Semi-Analytic Tests To Solidify

Before French geometry is ready, these small tests matter more than more TEAM7 runs.

### 5.1 Single loop on-axis Bz analytic test

Recommended name:

```text
test_3d_ophelie_biot_savart_circular_loop_axis
```

Known analytic:

```math
B_z(z)
=
\frac{\mu_0 I R^2}
{2(R^2+z^2)^{3/2}}
```

Test content:

```text
1. discretize filament circular loop into line quadrature points;
2. sample Bz at multiple on-axis z;
3. compare to analytic formula;
4. validate B direction, units, I scaling.
```

Recommended acceptance:

```text
relative error < 1%   # with enough line segments
```

### 5.2 Multiloop superposition test

Recommended name:

```text
test_3d_ophelie_biot_savart_multiloop_axis
```

Setup:

```text
multiple circular loops along z;
same current per loop;
Bz analytic = sum of single-loop formulas.
```

Purpose:

```text
validate basis for French annular/multi-turn coil source.
```

### 5.3 B cross-product direction test

If `test_3d_ophelie_biot_savart_direction` exists, keep as regression:

```text
source at (R,0,0)
J or dl along +y
observer at origin
expect Bz > 0
```

### 5.4 I scaling / Q scaling test

Recommended name:

```text
test_3d_ophelie_source_scaling
```

Acceptance:

```text
I -> 2I:
    A, B, E, J -> 2x
    raw Joule power -> 4x
```

### 5.5 Level0 uniform A analytic test

Recommended name:

```text
test_3d_ophelie_level0_uniform_A
```

Specify:

```text
ASrcReal = constant vector
ASrcImag = 0
```

Check:

```text
EImag = -omega * ASrcReal
JImag = sigma * EImag
Q = 0.5 * sigma * |EImag|^2
```

---

## 6. Power Scaling Must Be Normalized

### 6.1 Current direction

For customer heat-source demo, use scaled heat source:

```text
JouleHeat = scaled Q for thermal input
```

But preserve raw diagnostics:

```text
P_raw
P_target
power_scale
field_scale
effective_current_amplitude
```

### 6.2 Recommended implementation

For VTP field self-consistency:

```math
power\_scale = P_{target}/P_{raw}
```

```math
field\_scale = \sqrt{power\_scale}
```

Scaling:

```text
A, B, Phi, GradPhi, E, J multiply by field_scale
JouleHeat multiply by power_scale
```

So:

```math
Q = 0.5 J \cdot E
```

still holds.

### 6.3 CLI

Must support:

```text
--target-power=50000
--no-power-scaling
```

Where:

```text
--no-power-scaling
    for reference / analytic comparison;

--target-power
    for customer heat-source demo.
```

---

## 7. French Reduced Geometry Preparation

Before real French geometry, implement parametric reduced case:

```text
test_3d_ophelie_french_reduced
```

### 7.1 Body setup

First version needs:

```text
GlassBody:
    cylinder
    radius = 0.325 m
    height = 0.50 m
    center = configurable

CoilSource:
    multi circular loops or annular source
    radius / z positions configurable
    no coil internal potential solve

Optional CrucibleWallBody:
    geometry / future thermal wall only
    not in EM
```

Cold-crucible metal wall first version not in EM:

```text
crucible wall is passive geometry only
```

### 7.2 CLI parameters

Recommended:

```text
--glass-radius=0.325
--glass-height=0.50
--coil-radius=0.40
--coil-turns=8
--coil-current=1.0
--frequency=300000
--sigma=16
--target-power=50000
--coil-source-model=filament-circular-multiloop
--coil-z-min=...
--coil-z-max=...
--coil-num-loops=...
--state_recording=1
```

### 7.3 Source model

French reduced case priority:

```text
filament-circular-multiloop
```

Multiple circular loop line integral sources:

```math
A(x)=\frac{\mu_0}{4\pi}\sum_{loops} I \oint \frac{dl}{|x-x'|}
```

```math
B(x)=\frac{\mu_0}{4\pi}\sum_{loops} I \oint \frac{dl\times r}{|r|^3}
```

For equivalent N turns:

```text
Ieff per loop = total_ampere_turns / number_of_loops
```

Or each loop is one turn:

```text
number_of_loops = N
I per loop = I0
```

If N is large, do not build 2742 loops; use weighted loops.

### 7.4 Output fields

Minimum:

```text
ASrcReal
ASrcImag
BSrcReal
BSrcImag
PhiImag
GradPhiImag
EImag
JImag
JouleHeat
Sigma
```

With power scaling also:

```text
JouleHeatRaw
JouleHeatScaled or PowerScale log
```

---

## 8. Phi Correction Tasks in French Reduced Case

### 8.1 Goal

French OPHELIE form:

```math
j = -\sigma(\nabla V + i\omega A)
```

Phi/V correction is a main-line capability. In French reduced case, first version in two steps:

```text
Level0:
    E = -i omega A
    J = sigma E
    Q = 0.5 sigma |E|²

Level1:
    solve PhiImag
    EImag = -GradPhiImag - omega ASrcReal
    JImag = sigma EImag
    Q = 0.5 sigma |EImag|²
```

### 8.2 Acceptance metrics

Do not only look at `phi_rel_res`. Must output:

```text
divJ_level0
divJ_phi
divJ_reduction = divJ_level0 / divJ_phi
```

On regular particles / simple cylinder:

```text
divJ_phi < 0.3 * divJ_level0    # initial
divJ_phi < 0.1 * divJ_level0    # stable target
```

### 8.3 Required phi tests

Recommended names:

```text
test_3d_ophelie_phi_laplace_constant
test_3d_ophelie_phi_rhs_constant_A
test_3d_ophelie_phi_reduces_divJ
```

Content:

```text
1. constant PhiImag -> Laplace ≈ 0;
2. constant ASrcReal -> RHS ≈ 0;
3. after phi correction divJ decreases.
```

---

## 9. JouleHeat -> Thermal Interface Preparation

Main task is EM, but for customer demo JouleHeat must feed heat module.

Cursor should reserve:

```text
JouleHeat field:
    scaled W/m^3 heat source

TotalJoulePower:
    sum(Q_i * Vol_i)

Thermal input:
    dT/dt = Q / (rho * cp)
```

First version: no complex heat flow; one-step or dummy thermal test:

```text
test_3d_ophelie_joule_to_heat_one_way
```

Check:

```text
under uniform Q temperature rise satisfies:
DeltaT = Q * dt / (rho * cp)
```

Later French case output:

```text
Q distribution
Total power
mean / max glass heating rate
```

---

## 10. How To Close Out Current TEAM7 Work

### 10.1 Results to archive

Write current TEAM7 source smoke into README:

```text
volume-etheta:
    RMS 66.5%, peak_x=198

filament/racetrack sweep:
    best RMS ~44%, but peak_x still 198, or RMS breaks when peak aligned

Conclusion:
    simple coil-only source path cannot stably reproduce TEAM7 Bz reference;
    TEAM7 next needs total-field / plate-induced B;
    stop path fitting now.
```

### 10.2 If TEAM7 continues later

Only:

```text
probe total B = coil_B + plate_induced_B
phase0 / phase90
50 Hz / 200 Hz
```

Do not:

```text
more inset/z fitting
surface shell + e_theta
blind Jy chasing
```

---

## 11. Cursor Work Not Yet Complete

### P0: Must complete

```text
1. git commit current OPHELIE module;
2. clean obsolete reference bundles and builds;
3. solidify current moment abstraction;
4. add analytic circular loop Bz test;
5. add multi-loop Bz test;
6. add source scaling test;
7. add Level0 uniform A test;
8. confirm --no-power-scaling and field scaling logic;
9. create test_3d_ophelie_french_reduced;
10. implement filament-circular-multiloop source;
11. output French reduced case A/B/E/J/Q/VTP;
12. add divJ diagnostics.
```

### P1: Complete soon

```text
1. PhiImag correction on French reduced cylinder;
2. divJ_phi clearly less than divJ_level0;
3. JouleHeat power normalization to 50 kW;
4. JouleHeat -> dummy thermal one-way test;
5. generate docs/ophelie/FRENCH_REDUCED_CASE_README.md.
```

### P2: Later enhancements

```text
1. probe total B = coil + plate-induced for TEAM7;
2. self-induction complex channel fix;
3. volume-racetrack for TEAM7 full-chain source;
4. AirSampleBody / field sampling body;
5. SurfaceParticles / shell current source;
6. cold crucible passive conducting shell.
```

---

## 12. Recommended Development Order

### Step 1: Close TEAM7 source smoke

```text
write existing TEAM7 Bz / racetrack sweep results to docs;
stop path fitting;
keep test case but not main line.
```

### Step 2: Add Biot analytic tests

Ensure:

```text
circular loop Bz analytic
multi-loop Bz analytic
scaling
uniform A Level0
```

all pass.

### Step 3: Implement French reduced cylinder case

New test:

```text
test_3d_ophelie_french_reduced
```

First without phi:

```text
multi-loop coil -> glass A/B -> E/J/Q -> power normalization -> VTP
```

### Step 4: Add phi correction

```text
PhiImag solver
GradPhiImag
E/J/Q with phi
divJ diagnostics
```

### Step 5: Joule to thermal

```text
Q as heat source;
one-way/dummy thermal first;
real fluid-thermal later.
```

### Step 6: TEAM7 total-field later

If time:

```text
plate-induced B at probe
phase0/phase90
50/200 Hz
```

---

## 13. One-Line Instruction for Cursor

```text
Stop using TEAM7 as main validation. TEAM7 showed simple coil-only source fitting is exhausted; archive smoke and keep total-field benchmark only. Next stage main line: French reduced cold-crucible case—first lock Biot-Savart with analytic/multi-loop circular coil, then cylindrical glass A/B/E/J/Q, phi correction, power normalization, and JouleHeat thermal interface. Do not continue shell or more racetrack fitting unless French reduced main line is running.
```

---

## 14. Final Goal

Before French geometry body is ready, Cursor should deliver:

```text
1. stable OPHELIE-like electromagnetic module;
2. analytic Biot-Savart regression tests;
3. runnable French reduced cylinder induction case;
4. VTP output A/B/E/J/Q/Phi/JouleHeat;
5. JouleHeat normalizable to target power;
6. phi correction with divJ reduction diagnostic;
7. TEAM7 source-smoke results archived but not blocking main line.
```

When real French geometry is ready, only replace:

```text
GlassBody geometry
CoilSource geometry / loop distribution
Crucible wall visualization
```

without rebuilding the EM chain.
