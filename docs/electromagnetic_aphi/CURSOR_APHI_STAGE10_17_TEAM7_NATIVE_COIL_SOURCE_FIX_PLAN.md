# CURSOR_APHI_STAGE10_17_TEAM7_NATIVE_COIL_SOURCE_FIX_PLAN.md

> Branch: `feature/electromagnetic`  
> Scope: SPHinXsys SYCL CK A-phi electromagnetic heating solver + TEAM7 native geometry/reference validation  
> Current stage: Stage 10.17 native TEAM7 geometry and Bz reference pipeline  
> Purpose: Fix the current native TEAM7 reload case by correcting the physical unit system and coil source definition.  
> Priority: Highest. Do **not** tune GMRES/probe gates before fixing units and source.

---

## 0. Executive summary

The native TEAM7 pipeline has made an important engineering step:

```text
TEAM7 STL geometry
-> SPHinXsys native level-set / particle generation / relaxation / reload
-> Air / Coil / Plate three-body Contact
-> source-driven A-phi GMRES
-> B = curl(A)
-> Bz probe CSV
-> comparison with TEAM7 reference data
```

However, the current result is not a physical TEAM7 validation yet.

Current observed facts:

```text
GMRES final_true_rel ≈ 0.04
our_Bz ≈ 1e-9 mT
reference Bz ≈ 0.5–8 mT
profile_rel_err ≈ 1.0
```

Interpretation:

```text
The discrete A-phi system can be solved.
But the physical Bz magnitude is wrong by many orders of magnitude.
This is not primarily a GMRES tolerance problem.
This is not primarily a Contact A-penalty problem.
This is not primarily a sign-convention problem.
The main problems are:
    1. geometry/unit mismatch: native STL/reload is in mm, but EM equations use SI;
    2. magnetic reluctivity nu is likely not SI;
    3. coil source is currently a placeholder uniform z-directed RHS, not TEAM7 2742 AT winding current.
```

Immediate next step:

```text
Fix native TEAM7 units and coil source before any more reference tuning.
```

---

## 1. Frozen decisions for this stage

Do not change these unless explicitly instructed:

```text
production lambda_A = off
Contact A-penalty = research only
projection = deferred
ghost/mirror EM boundary = deferred
multiresolution = deferred
cold crucible = deferred
real helical copper wire = deferred
coil self-heating = deferred
coil circuit coupling = deferred
surface shell current = deferred
```

Current model should remain:

```text
Coil/source region:
    source-only prescribed current region
    sigma_coil = 0
    no eddy current in coil
    no coil Joule as physical heat target

Plate:
    conductive aluminum
    sigma_plate = 3.526e7 S/m
    mu_r = 1

Air:
    nonconductive magnetic domain
    air Joule is not physical heat target
```

---

## 2. What is already correct

### 2.1 Native particle generation / relaxation

The current particle generation pipeline is correct in direction:

```text
Coil STL
Plate STL
Air box minus coil/plate
LevelSetShape
ParticleGeneratorLattice
RandomizeParticlePositionCK
KernelGradientIntegral
LevelsetKernelGradientIntegral
PositionRelaxationCK
LevelsetBounding
ReloadParticleIOCK
BodyStatesRecordingToVtpCK
```

This follows SPHinXsys native geometry and relaxation style.

Do **not** replace it with a custom Python particle generator or hand-written random point sampler.

### 2.2 Three-body Contact topology

The current native TEAM7 solve uses:

```text
AirBody
CoilBody
PlateBody
```

with Contact coupling through air:

```text
Air contacts Coil and Plate
Coil contacts Air
Plate contacts Air
```

This topology is acceptable for TEAM7 because the coil and plate are separated by air.

### 2.3 Bz probe/reference pipeline

The current Bz probe pipeline exists:

```text
A1-B1 line:
    x = 0..288 mm
    y = 72 mm
    z = 34 mm

A2-B2 line:
    x = 0..288 mm
    y = 144 mm
    z = 34 mm
```

It writes:

```text
probe_A1_B1_Bz.csv
comparison_A1_B1_vs_reference.csv
summary.csv
```

This is the right direction. However, the current Bz magnitude is near zero because the source/units are wrong.

---

## 3. Current source implementation and why it is wrong

### 3.1 Current implementation

Current native source helper uses something equivalent to:

```text
AssignUniformImpressedCurrentRhsCK
```

Applied only on `CoilBody`.

For each coil particle:

```text
rhs_a_real = impressed_current_amplitude * (0, 0, 1)
rhs_a_imag = impressed_current_amplitude * (0, 0, 0)
```

Current settings:

```text
impressed_current_amplitude = 8.0
current_real = (0, 0, 1)
current_imag = (0, 0, 0)
```

This is a frequency-domain / phasor RHS, but it is only a placeholder.

### 3.2 Why this is not TEAM7

TEAM7 coil excitation is:

```text
2742 turns
1 A/turn
total ampere-turns NI = 2742 A
frequency = 50 Hz and 200 Hz
```

The coil current must flow along the coil winding direction, not uniformly along z.

Current source is wrong because:

```text
1. direction is wrong:
   current is currently uniform z-direction;
   TEAM7 current should be tangent to the closed coil path.

2. magnitude is wrong:
   amplitude=8 is arbitrary;
   TEAM7 requires NI=2742 A.

3. units are undefined:
   A-phi RHS should correspond to current density J_s [A/m^2] if LHS uses SI reluctivity nu=1/mu.

4. coil is source-only:
   sigma_coil=0 is fine,
   but source RHS must still represent correct ampere-turn excitation.
```

---

## 4. Do not use surface shell current as the first solution

The user asked whether current should be applied on shell/surface particles of the coil.

Decision:

```text
Do not use shell/surface current as the main TEAM7 source model now.
```

Reasons:

```text
1. TEAM7 coil is a multi-turn source, not a closed conducting shell with uniform surface current.
2. Shell particles include inner wall, outer wall, top and bottom surfaces; current direction and weighting become ambiguous.
3. Surface current density K has units A/m, while current A-phi RHS likely expects volume current density J_s [A/m^2].
4. It is easy to double-count current on multiple surfaces.
5. For both A-phi and Ophelie/Biot-Savart versions, a centerline/volume-current source is easier to calibrate against NI=2742 A.
```

Recommended source models:

```text
A-phi solver:
    source-only coil volume with prescribed current density J_s [A/m^2]

Ophelie / Biot-Savart long-range solver:
    first version: centerline line current NI = 2742 A
    optional later: same volume current density J_s for volume integration
```

---

## 5. Required SI unit conversion

### 5.1 Absolute rule

All native TEAM7 EM solve coordinates must be SI meters.

Do not run A-phi EM equations on mm coordinates.

Recommended conversion:

```text
If STL is in mm:
    import STL with scale = 0.001
or:
    transform particle positions to meters immediately after reload
```

Preferred: scale at geometry import/generation time.

Use:

```text
dp_0 = 0.006 m
air box:
    lower = (-0.05, -0.05, -0.05) m
    upper = ( 0.35,  0.35,  0.20) m
probe A1-B1:
    x = 0..0.288 m
    y = 0.072 m
    z = 0.034 m
probe A2-B2:
    x = 0..0.288 m
    y = 0.144 m
    z = 0.034 m
Jy diagnostic:
    x = 0..0.288 m
    y = 0.072 m
    z = 0.01899 m
```

CSV can still report x in mm:

```text
x_mm = x_m * 1000
```

### 5.2 Fix geometry audit print helpers

If current helper prints:

```cpp
printVec3dMm(x) = 1000 * x
```

then it is only correct if `x` is already in meters.

After this stage:

```text
Position must be stored in meters.
All *_mm prints must be derived by multiplying meter coordinates by 1000.
```

Do not mix mm and m in helper names or test outputs.

---

## 6. Required SI magnetic reluctivity

Current native setup reportedly uses:

```text
nu = 1
```

This is only acceptable for nondimensional tests. It is not acceptable for comparison against Tesla/mT reference.

If `nu` in the A-phi operator is magnetic reluctivity:

```text
nu = 1 / mu
```

then for air, aluminum, and nonmagnetic coil region:

```text
mu0 = 4 * pi * 1e-7 H/m
nu0 = 1 / mu0 ≈ 7.957747154594767e5
```

Required:

```text
Air nu = 1/mu0
Plate nu = 1/mu0
Coil/source nu = 1/mu0
```

Before changing, confirm from `AphiMatrixFreeOperator` whether `nu` is indeed reluctivity. If yes, use SI `nu0`.

If previous simplified benchmark uses `nu=1`, keep that for nondimensional tests only. Native TEAM7 reference path must be SI.

---

## 7. TEAM7 coil centerline definition from the official figure

The TEAM7 coil is a rounded-rectangle / racetrack-like coil.

From the official plan-view figure:

```text
Outer corner radius: R_out = 50 mm
Inner corner radius: R_in  = 25 mm
Coil planar width:   w = R_out - R_in = 25 mm
Centerline radius:   R_c = (R_out + R_in) / 2 = 37.5 mm
```

From the cross-section:

```text
Plate thickness = 19 mm
Plate-to-coil gap = 30 mm
Coil height = 100 mm
```

If z=0 is the plate bottom:

```text
Plate z range: 0..19 mm
Coil z range: 49..149 mm
Coil z center: z_c = 99 mm
```

The observed STL/native bbox is approximately:

```text
x_min = 94 mm
x_max = 294 mm
y_min = 0 mm
y_max = 200 mm
z_min = 49 mm
z_max = 149 mm
```

After converting to SI:

```text
x_min = 0.094 m
x_max = 0.294 m
y_min = 0.000 m
y_max = 0.200 m
z_min = 0.049 m
z_max = 0.149 m
```

---

## 8. Correct TEAM7 coil source: volume current density along winding path

### 8.1 Formula

Use a source-only coil volume current density:

```text
J_s(x_i) = J0 * t_i
```

where:

```text
t_i = unit tangent direction of coil winding path at the closest point to particle i
J0 = NI / A_eff
NI = 2742 A
A_eff = effective cross-sectional area of coil bundle
```

For the TEAM7 figure:

```text
coil planar width = 25 mm = 0.025 m
coil height = 100 mm = 0.100 m
A_eff ≈ 0.025 * 0.100 = 0.0025 m^2
J0 = 2742 / 0.0025 = 1.0968e6 A/m^2
```

A more robust particle-based estimate:

```text
L_path = centerline total length
V_coil = sum_i V_i over coil particles
A_eff = V_coil / L_path
J0 = NI / A_eff
```

The two estimates should be close. Output both.

### 8.2 Expected path length

The centerline consists of four straight segments plus four corner arcs.

Using:

```text
R_c = 37.5 mm = 0.0375 m
straight segment length = 100 mm = 0.100 m
```

Approximate total:

```text
L_path = 4 * 0.100 + 2*pi*0.0375
       ≈ 0.400 + 0.235619
       ≈ 0.635619 m
```

If using particle volume:

```text
A_eff = V_coil / 0.635619
```


---

## 9. Centerline construction: exact rounded-rectangle definition

Let the coil bbox be in SI meters:

```cpp
Real x_min = coil_bbox.lower_bound_[0];
Real x_max = coil_bbox.upper_bound_[0];
Real y_min = coil_bbox.lower_bound_[1];
Real y_max = coil_bbox.upper_bound_[1];
Real z_min = coil_bbox.lower_bound_[2];
Real z_max = coil_bbox.upper_bound_[2];

Real R_out = 0.050;   // 50 mm
Real R_in  = 0.025;   // 25 mm
Real R_c   = 0.5 * (R_out + R_in); // 0.0375 m
Real z_c   = 0.5 * (z_min + z_max);
```

Corner centers:

```cpp
Vec2d C_lb(x_min + R_out, y_min + R_out);
Vec2d C_rb(x_max - R_out, y_min + R_out);
Vec2d C_rt(x_max - R_out, y_max - R_out);
Vec2d C_lt(x_min + R_out, y_max - R_out);
```

For the expected TEAM7 bbox:

```text
C_lb = (0.144, 0.050) m
C_rb = (0.244, 0.050) m
C_rt = (0.244, 0.150) m
C_lt = (0.144, 0.150) m
```

Centerline straight sections:

```text
bottom:
    from (0.144, 0.0125, 0.099) to (0.244, 0.0125, 0.099)

right:
    from (0.2815, 0.050, 0.099) to (0.2815, 0.150, 0.099)

top:
    from (0.244, 0.1875, 0.099) to (0.144, 0.1875, 0.099)

left:
    from (0.1065, 0.150, 0.099) to (0.1065, 0.050, 0.099)
```

Centerline arcs:

```text
bottom-right:
    center (0.244, 0.050, 0.099), radius 0.0375
    angle -90° -> 0°

top-right:
    center (0.244, 0.150, 0.099), radius 0.0375
    angle 0° -> 90°

top-left:
    center (0.144, 0.150, 0.099), radius 0.0375
    angle 90° -> 180°

bottom-left:
    center (0.144, 0.050, 0.099), radius 0.0375
    angle 180° -> 270°
```

This direction is counter-clockwise viewed from +z. If reference sign is reversed, multiply all tangents by -1.

---

## 10. Recommended implementation: dense centerline polyline

To avoid complicated segment classification, build a closed polyline and project each coil particle onto it.

### 10.1 C++ helper: build centerline polyline

```cpp
std::vector<Vecd> buildTeam7CoilCenterline(const BoundingBox &coil_bbox)
{
    constexpr Real Pi = 3.14159265358979323846;

    const Real R_out = 0.050; // 50 mm
    const Real R_in  = 0.025; // 25 mm
    const Real R_c   = 0.5 * (R_out + R_in);

    const Real x_min = coil_bbox.lower_bound_[0];
    const Real x_max = coil_bbox.upper_bound_[0];
    const Real y_min = coil_bbox.lower_bound_[1];
    const Real y_max = coil_bbox.upper_bound_[1];
    const Real z_c   = 0.5 * (coil_bbox.lower_bound_[2] + coil_bbox.upper_bound_[2]);

    const Real cxL = x_min + R_out;
    const Real cxR = x_max - R_out;
    const Real cyB = y_min + R_out;
    const Real cyT = y_max - R_out;

    std::vector<Vecd> pts;
    pts.reserve(256);

    auto add_line = [&](const Vecd &a, const Vecd &b, int n)
    {
        for (int k = 0; k < n; ++k)
        {
            const Real s = Real(k) / Real(n);
            pts.push_back(a * (1.0 - s) + b * s);
        }
    };

    auto add_arc = [&](const Vecd &c, Real a0, Real a1, int n)
    {
        for (int k = 0; k < n; ++k)
        {
            const Real s = Real(k) / Real(n);
            const Real a = a0 * (1.0 - s) + a1 * s;
            pts.push_back(c + Vecd(R_c * std::cos(a), R_c * std::sin(a), 0.0));
        }
    };

    // CCW current direction viewed from +z.
    // bottom straight
    add_line(Vecd(cxL, y_min + 0.5 * (R_out - R_in), z_c),
             Vecd(cxR, y_min + 0.5 * (R_out - R_in), z_c), 32);

    // bottom-right arc
    add_arc(Vecd(cxR, cyB, z_c), -0.5 * Pi, 0.0, 32);

    // right straight
    add_line(Vecd(x_max - 0.5 * (R_out - R_in), cyB, z_c),
             Vecd(x_max - 0.5 * (R_out - R_in), cyT, z_c), 32);

    // top-right arc
    add_arc(Vecd(cxR, cyT, z_c), 0.0, 0.5 * Pi, 32);

    // top straight
    add_line(Vecd(cxR, y_max - 0.5 * (R_out - R_in), z_c),
             Vecd(cxL, y_max - 0.5 * (R_out - R_in), z_c), 32);

    // top-left arc
    add_arc(Vecd(cxL, cyT, z_c), 0.5 * Pi, Pi, 32);

    // left straight
    add_line(Vecd(x_min + 0.5 * (R_out - R_in), cyT, z_c),
             Vecd(x_min + 0.5 * (R_out - R_in), cyB, z_c), 32);

    // bottom-left arc
    add_arc(Vecd(cxL, cyB, z_c), Pi, 1.5 * Pi, 32);

    return pts;
}
```

### 10.2 C++ helper: tangent from polyline

```cpp
Vecd tangentFromClosedPolyline(const Vecd &x, const std::vector<Vecd> &path)
{
    Real best_d2 = std::numeric_limits<Real>::max();
    Vecd best_tangent = Vecd::Zero();

    for (size_t k = 0; k < path.size(); ++k)
    {
        const Vecd &a = path[k];
        const Vecd &b = path[(k + 1) % path.size()];
        const Vecd ab = b - a;

        const Real ab2 = dot(ab, ab);
        if (ab2 <= TinyReal) continue;

        Real s = dot(x - a, ab) / ab2;
        s = std::max(Real(0), std::min(Real(1), s));

        const Vecd q = a + s * ab;
        const Real d2 = dot(x - q, x - q);

        if (d2 < best_d2)
        {
            best_d2 = d2;
            best_tangent = ab / std::sqrt(ab2);
        }
    }

    return best_tangent;
}
```

For reverse current direction:

```cpp
best_tangent *= -1.0;
```

---

## 11. Required new CK source class

Implement a new source class. Suggested name:

```text
AssignTeam7CoilPathImpressedCurrentRhsCK
```

Suggested location:

```text
tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/
    aphi_team7_native_coil_source_helpers.h
```

or:

```text
tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/benchmark/
    aphi_team7_native_coil_source_ck.h/.hpp
```

### 11.1 Specification

```cpp
struct Team7CoilPathSourceSpec
{
    Real turns = 2742.0;
    Real current_per_turn = 1.0;
    Real total_ampere_turns = 2742.0;

    Real R_out = 0.050;
    Real R_in = 0.025;
    Real R_center = 0.0375;

    bool reverse_winding = false;

    // If true, use particle volume / path length.
    // If false, use width * height = 0.025 * 0.100.
    bool use_particle_volume_cross_section = true;

    Real fallback_cross_section_area = 0.025 * 0.100; // m^2

    // Optional source scaling for diagnostics.
    Real source_scale = 1.0;
};
```

### 11.2 Source strength

```cpp
Real NI = spec.turns * spec.current_per_turn;
Real A_eff = spec.use_particle_volume_cross_section ? coil_volume / path_length
                                                     : spec.fallback_cross_section_area;
Real J0 = spec.source_scale * NI / A_eff;
```

Expected first value:

```text
A_eff ≈ 0.0025 m^2
J0 ≈ 1.0968e6 A/m^2
```

### 11.3 Device update

For each coil particle:

```cpp
Vecd t_i = tangent_i; // precomputed or computed from path
rhs_a_real[i] = J0 * t_i;
rhs_a_imag[i] = Vecd::Zero();
```

Because closed polyline search over all segments inside a device kernel may be cumbersome, recommended approach:

```text
Precompute tangent_i on host after reload and sync to device;
then CK UpdateKernel only assigns rhs using stored CoilSourceTangent variable.
```

This is acceptable if explicit host/device sync is done.

### 11.4 Variables to register

Add coil particle variables:

```text
CoilSourceTangent        Vecd
CoilSourceCurrentDensity Vecd
CoilSourcePathDistance   Real optional
```

---

## 12. Source audit required

Add a source audit function.

Print:

```text
coil_particle_count
coil_volume
coil_path_length
A_eff_particle = coil_volume / path_length
A_eff_fallback = 0.0025
A_eff_used
turns
current_per_turn
total_ampere_turns
J0
source_scale
max_abs_Js
mean_abs_Js
min_tangent_norm
max_tangent_norm
rhs_l2
estimated_integrated_current = J0 * A_eff_used
```

Hard gate:

```text
coil_particle_count > 0
path_length > 0
A_eff_used > 0
J0 > 1e5 A/m^2
0.9 < estimated_integrated_current / 2742 < 1.1
0.9 < min_tangent_norm <= max_tangent_norm < 1.1
rhs_l2 > 0
```


---

## 13. Relationship with Ophelie / long-range integral solver

Both A-phi and Ophelie versions should use the same coil centerline definition.

### 13.1 A-phi solver

Use volume current:

```text
J_s(x_i) = J0 * t_i [A/m^2]
```

This is assigned to coil particles as the A-equation impressed source / RHS.

### 13.2 Ophelie / Biot-Savart first version

Use line current along the same centerline:

```text
I_equiv = N I = 2742 A
```

For vector potential:

```text
A(r) = mu0 / (4*pi) * integral_over_path [ I_equiv * dl / |r-r'| ]
```

This is the cleanest way to compare the two solvers.

### 13.3 Optional later volume Biot-Savart

Use:

```text
A(r) = mu0 / (4*pi) * integral_over_coil_volume [ J_s(r') / |r-r'| dV' ]
```

This should converge toward the line-current result away from the coil cross-section scale.

### 13.4 Do not use shell current first

Do not use surface shell current as the first TEAM7 source model.

---

## 14. Validation sequence after source fix

Do not directly go to full TEAM7 reference error.

### 14.1 P0 — SI geometry audit

Run native geometry audit and verify:

```text
plate bbox approximately:
    x = 0..0.294 m
    y = 0..0.294 m
    z = 0..0.019 m

coil bbox approximately:
    x = 0.094..0.294 m
    y = 0..0.200 m
    z = 0.049..0.149 m

probe A1-B1:
    x = 0..0.288 m
    y = 0.072 m
    z = 0.034 m
```

### 14.2 P1 — source audit only

Run source assignment without solving.

Check:

```text
J0 ≈ 1.0968e6 A/m^2
tangent norms ≈ 1
rhs_l2 > 0
estimated integrated current ≈ 2742 A
```

### 14.3 P2 — vacuum/source B sanity

Set:

```text
plate sigma = 0
coil source on
air on
nu = 1/mu0
SI geometry
```

Run A-phi or source field solve.

Goal:

```text
Bz should no longer be ~1e-9 mT.
Expected order should move toward mT scale.
```

Do not require exact reference yet.

### 14.4 P3 — full plate eddy TEAM7 50 Hz

Set:

```text
plate sigma = 3.526e7 S/m
frequency = 50 Hz
source = corrected coil-path J_s
nu = 1/mu0
```

Run:

```text
Bz A1-B1
Bz A2-B2
Jy diagnostic
Joule diagnostic
```

### 14.5 P4 — sign convention audit

Only after Bz magnitude is nonzero, apply sign scan:

```text
(+real,+imag)
(-real,+imag)
(+real,-imag)
(-real,-imag)
```

Expected likely convention from muFEM:

```text
phase0 = -B_real.z * 1000 [mT]
phase90 =  B_imag.z * 1000 [mT]
```

But do not hard-code this before the magnitude issue is fixed.

---

## 15. Probe sampling improvement

Current nearest-particle sampling is acceptable for smoke.

For reference comparison, improve sampling:

```text
Prefer AirBody particles for Bz probe because A1-B1/A2-B2 lines are in air.
Print nearest body name.
Print nearest distance.
Reject probe if nearest distance is too large.
```

Initial threshold:

```text
max_nearest_distance < 2 * dp
```

If using `dp=0.006 m`, threshold:

```text
max_nearest_distance < 0.012 m
```

Later implement SPH/kernel interpolation if nearest sampling is too noisy.

---

## 16. What not to tune right now

Do not spend time on:

```text
GMRES tolerance tuning
profile error gate tightening
sign convention
Contact A-penalty
projection
multiresolution
cold crucible
surface shell current
thermal coupling
```

until after:

```text
SI units + nu=1/mu0 + 2742 AT tangential coil source
```

are fixed and source/vacuum B sanity passes.

---

## 17. Required tests to add or update

### Test A — source audit

```text
test_3d_aphi_ck_team7_native_coil_source_audit
```

Checks source path, `J0`, tangent norms, integrated current.

### Test B — SI geometry audit update

Update existing:

```text
test_3d_aphi_ck_team7_native_geometry_audit
```

Must fail if coordinates still look like mm values in EM solve.

Example:

```text
coil bbox x_max should be around 0.294 m, not 294.
```

### Test C — vacuum B sanity

```text
test_3d_aphi_ck_team7_native_vacuum_source_b_sanity
```

Uses:

```text
plate sigma = 0
source on
nu = 1/mu0
SI geometry
```

Gate:

```text
Bz finite
Bz max > small threshold
Bz not ~1e-9 mT
```

Keep threshold loose.

### Test D — corrected Bz reference smoke

Update:

```text
test_3d_aphi_ck_team7_native_geometry_bz_reference_probe
```

Only after source fix. Do not hard gate <20% yet.

---

## 18. Closure criteria for this fix stage

This stage can close when:

```text
1. Native TEAM7 coordinates are SI meters;
2. nu uses SI 1/mu0 in native reference path;
3. coil source uses tangential centerline current, not uniform z RHS;
4. source audit estimates integrated current ≈ 2742 A;
5. vacuum/source B sanity produces non-negligible Bz;
6. full plate 50 Hz run produces Bz in the correct mT order of magnitude;
7. Bz probe sign convention can be meaningfully evaluated.
```

Do not claim TEAM7 validation is complete until Bz profile error is actually reduced against reference data.

---

## 19. Final one-line instruction for Cursor

```text
The native TEAM7 reload pipeline already works, but the physical input is wrong. Do not tune GMRES or gates. Convert native geometry to SI meters, use nu=1/mu0, replace uniform z RHS with a 2742 AT tangential rounded-rectangle coil-path current density, audit the source, run vacuum B sanity, then rerun Bz A1-B1/A2-B2 reference comparison.
```
