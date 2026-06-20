# OPHELIE-like French Reduced Cold-Crucible Geometry Generation and CoilSource Implementation Plan

This document is for Cursor. The goal is to build a French-literature-style reduced cold-crucible induction case using SPHinXsys analytic/code-generated geometry **without STL files** and **without segmented crucible**, reusing the current OPHELIE-like / Biot-Savart / phi / JouleHeat chain.

Current main line:

```text
analytic cylindrical GlassBody
    + code-generated multi circular loop CoilSource
        -> Biot-Savart A_src / B_src on glass particles
            -> optional PhiImag correction
                -> E / J / JouleHeat
                    -> optional power normalization
                        -> later JouleHeat thermal input
```

Explicitly out of scope:

```text
STL geometry import
segmented cold crucible
cold-crucible metal EM solve
surface / skin-depth model
water-cooling boundary
sigma(T)
thermal-EM two-way iteration
real helical coil and lead wires
TEAM7 path fitting
full PDE A-phi solver
```

---

## 1. Why First Version Does Not Use STL

French reduced case first-version geometry is unconventional:

```text
GlassBody:
    cylindrical conductive glass body

CoilSource:
    external multiloop line source; no solid coil particles required

CoilVisualBody:
    optional annular cylinder or multiloop markers for display only

CrucibleWallVisualBody:
    optional full cylindrical shell + bottom for display or future thermal boundary
```

Advantages of no STL:

```text
1. avoid mm/m unit errors;
2. avoid STL bbox, normal, closure, and reload issues;
3. no relax required;
4. parameters controllable for radius/height/coil radius/loop count sweeps;
5. single-loop / multiloop analytic validation of Biot-Savart;
6. better fit for current OPHELIE-like reduced model.
```

So first version should use SPHinXsys analytic geometry / code-generated source, not STL.

---

## 2. Add Test Case

Recommended add:

```text
tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/
    test_3d_ophelie_french_reduced.cpp
    CMakeLists.txt
    README.md
```

CMake wiring: reference existing:

```text
test_3d_ophelie
test_3d_ophelie_team7
```

Include:

```cpp
#include "electromagnetic_ophelie.h"
```

Later optional headers:

```cpp
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_multiloop_source.h"
#include "electromagnetic_ophelie_line_source_biot.h"
```

---

## 3. Default Parameters

Recommended parameter struct:

```cpp
struct OphelieFrenchReducedParameters
{
    // Glass melt
    Real glass_radius_ = 0.325;       // [m], French literature cold-crucible diameter 650 mm
    Real glass_height_ = 0.50;        // [m], first-version simplified height
    Real dp_ = 0.02;                  // [m], first smoke; 0.01 for finer run

    // Coil source
    Real coil_radius_ = 0.40;         // [m], slightly outside glass
    Real coil_z_min_ = 0.05;          // [m]
    Real coil_z_max_ = 0.45;          // [m]
    int coil_num_loops_ = 8;
    int coil_segments_per_loop_ = 256;

    // EM
    Real frequency_ = 300.0e3;        // [Hz]
    Real sigma_glass_ = 16.0;         // [S/m]
    Real current_amplitude_ = 1.0;    // [A]
    Real ampere_turns_ = 8.0;         // [A-turn], default coil_num_loops * current
    Real target_joule_power_ = 50.0e3;// [W]
    bool power_scaling_ = true;
    bool use_phi_ = false;
};
```

CLI support:

```text
--glass-radius=0.325
--glass-height=0.50
--dp=0.02

--coil-radius=0.40
--coil-z-min=0.05
--coil-z-max=0.45
--coil-num-loops=8
--coil-segments-per-loop=256

--frequency=300000
--sigma=16
--current=1
--ampere-turns=8
--target-power=50000
--no-power-scaling
--no-phi
--phi-solver=GMRES
--state_recording=1
```

---

## 4. GlassBody: Cylindrical Conductive Glass Geometry

### 4.1 Geometry

First version uses z-axis cylinder:

```text
radius = glass_radius
z_min = 0
z_max = glass_height
center = (0, 0, glass_height/2)
```

Geometry test:

```cpp
bool isInGlassCylinder(const Vecd &p, Real R, Real z_min, Real z_max)
{
    Real r2 = p[0] * p[0] + p[1] * p[1];
    return r2 <= R * R && p[2] >= z_min && p[2] <= z_max;
}
```

### 4.2 SPHinXsys generation

If current SPHinXsys branch has cylinder shape, use it. Otherwise add analytic shape. Recommended form:

```cpp
class GlassCylinderShape : public ComplexShape
{
public:
    GlassCylinderShape(const std::string &shape_name, Real radius, Real height)
        : ComplexShape(shape_name)
    {
        // Use native cylinder primitive if available.
        // Otherwise implement custom contain check or level-set shape.
    }
};
```

Generate particles:

```cpp
FluidBody glass_body(system, makeShared<GlassCylinderShape>(
    "GlassCylinder", params.glass_radius_, params.glass_height_));

glass_body.defineParticlesAndMaterial<BaseParticles, SomeMaterial>();
glass_body.generateParticles<Lattice>();
```

If real fluid material is not wired yet, generic particle body is fine as long as EM variables and volume are registered.

### 4.3 GlassBody variables to register

```text
Sigma : Real

ASrcReal : Vecd
ASrcImag : Vecd
BSrcReal : Vecd
BSrcImag : Vecd

PhiImag : Real
GradPhiImag : Vecd

EReal : Vecd
EImag : Vecd
JReal : Vecd
JImag : Vecd

JouleHeatRaw : Real   # recommended
JouleHeat : Real      # scaled or raw depending on CLI

DivJImag : Real       # diagnostic
```

First version minimum output:

```text
Sigma
ASrcReal
BSrcReal
EImag
JImag
JouleHeat
```

---

## 5. CoilSource: Multiloop Line Integral Source, No Solid Particles

### 5.1 Basic idea

CoilSource first version is not an SPHBody. It is a set of line quadrature sources:

```cpp
struct OphelieLineSourcePoint
{
    Vecd pos_;          // source point position [m]
    Vecd moment_real_;  // I * dl, unit A*m
    Vecd moment_imag_;  // first zero
};
```

Biot-Savart uses current moment:

```math
A(x_i)=\frac{\mu_0}{4\pi}\sum_k\frac{I_k\Delta l_k}{|x_i-x_k|}
```

```math
B(x_i)=\frac{\mu_0}{4\pi}\sum_k
\frac{I_k\Delta l_k\times(x_i-x_k)}
{|x_i-x_k|^3}
```

Note:

```text
B direction must be CurrentMoment cross r
not r cross CurrentMoment.
```

### 5.2 Multiloop generation helper

Recommended add:

```text
electromagnetic_ophelie_multiloop_source.h
electromagnetic_ophelie_multiloop_source.hpp
```

Interface:

```cpp
class OphelieMultiLoopCircularSource
{
public:
    OphelieMultiLoopCircularSource(
        Real coil_radius,
        Real z_min,
        Real z_max,
        int num_loops,
        int segments_per_loop,
        Real total_ampere_turns);

    const StdVec<Vecd>& positions() const;
    const StdVec<Vecd>& momentsReal() const;
    const StdVec<Vecd>& momentsImag() const;

private:
    void build();
};
```

Build logic:

```cpp
void OphelieMultiLoopCircularSource::build()
{
    const Real I_loop = total_ampere_turns_ / Real(num_loops_);

    for (int l = 0; l < num_loops_; ++l)
    {
        Real z = z_min_ + (Real(l) + 0.5) / Real(num_loops_) * (z_max_ - z_min_);

        for (int k = 0; k < segments_per_loop_; ++k)
        {
            Real theta0 = 2.0 * Pi * Real(k) / Real(segments_per_loop_);
            Real theta1 = 2.0 * Pi * Real(k + 1) / Real(segments_per_loop_);
            Real theta_mid = 0.5 * (theta0 + theta1);

            Vecd pos = Vecd(
                coil_radius_ * std::cos(theta_mid),
                coil_radius_ * std::sin(theta_mid),
                z);

            // dl is a directed segment vector, not unit tangent.
            Vecd dl = Vecd(
                coil_radius_ * (std::cos(theta1) - std::cos(theta0)),
                coil_radius_ * (std::sin(theta1) - std::sin(theta0)),
                0.0);

            source_pos_.push_back(pos);
            moment_real_.push_back(I_loop * dl);
            moment_imag_.push_back(Vecd::Zero());
        }
    }
}
```

### 5.3 total_ampere_turns

First version:

```text
total_ampere_turns = coil_num_loops * current_amplitude
```

If customer provides real ampere-turns:

```text
--ampere-turns=...
```

If using power normalization, absolute current is not sensitive in first version:

```text
ampere_turns = 1 or 8
target_power = 50 kW
```

---

## 6. Line-Source Biot-Savart on GlassBody

### 6.1 New operator

Recommended add:

```text
electromagnetic_ophelie_line_source_biot.h
electromagnetic_ophelie_line_source_biot.hpp
```

Interface:

```cpp
class ComputeOphelieLineSourceBiotSavartOnBody
{
public:
    ComputeOphelieLineSourceBiotSavartOnBody(
        SPHBody &target_body,
        const StdVec<Vecd> &source_pos,
        const StdVec<Vecd> &source_moment_real,
        const StdVec<Vecd> &source_moment_imag,
        const OphelieFrenchReducedParameters &params);

    void exec();
};
```

Implementation choice:

```text
First version can host-side loop into glass arrays;
later convert to SYCL CK device kernel.

If framework already has device buffer/helper, CK directly.
```

### 6.2 Formula

For each glass particle:

```cpp
const Real coeff = params.mu0_ / (4.0 * Pi);
const Real eps2 = params.softening_ * params.softening_;

for each glass particle i:
    Vecd xi = pos[i];
    Vecd a_r = Zero;
    Vecd a_i = Zero;
    Vecd b_r = Zero;
    Vecd b_i = Zero;

    for each source k:
        Vecd r = xi - source_pos[k];
        Real r2 = dot(r, r) + eps2;
        Real inv_r = 1.0 / sqrt(r2);
        Real inv_r3 = inv_r / r2;

        Vecd mr = source_moment_real[k];
        Vecd mi = source_moment_imag[k];

        a_r += coeff * mr * inv_r;
        a_i += coeff * mi * inv_r;

        b_r += coeff * cross(mr, r) * inv_r3;
        b_i += coeff * cross(mi, r) * inv_r3;

    ASrcReal[i] = a_r;
    ASrcImag[i] = a_i;
    BSrcReal[i] = b_r;
    BSrcImag[i] = b_i;
```

Softening:

```text
softening = 0.25 * dp
```

---

## 7. CoilVisualBody: Optional Visualization, Not in EM

### 7.1 Purpose

Show external coil in ParaView; does not participate in Biot-Savart or phi.

### 7.2 Implementation options

Option A: annular cylinder visual body:

```text
inner_radius = coil_radius - 0.02
outer_radius = coil_radius + 0.02
z in [coil_z_min, coil_z_max]
```

Option B: multiloop marker particles:

```text
marker points along each loop, output only.
```

Recommend option A or skip entirely first; output GlassBody fields. Add visual body for customer demo later.

---

## 8. CrucibleWallVisualBody: Optional Full Cylindrical Shell, Not Segmented

User explicit: **no segmented crucible**.

CrucibleWallVisualBody is visualization / future thermal boundary only; not in EM.

Geometry:

```text
inner_radius = glass_radius
outer_radius = glass_radius + wall_thickness
z in [0, glass_height]
bottom optional
```

Default:

```text
wall_thickness = 0.02 m
```

Containment:

```cpp
bool isInCrucibleWall(const Vecd &p)
{
    Real r = sqrt(p[0]*p[0] + p[1]*p[1]);
    return r >= R_inner && r <= R_outer
        && p[2] >= 0.0 && p[2] <= glass_height;
}
```

Do not assign:

```text
Jsrc
Sigma
Phi
E/J/Q
```

---

## 9. EM Pipeline

### 9.1 Level 0

```math
\hat E = -i\omega \hat A
```

Code:

```cpp
EReal =  omega * ASrcImag;
EImag = -omega * ASrcReal;

JReal = sigma * EReal;
JImag = sigma * EImag;

JouleHeatRaw = 0.5 * sigma * (
    dot(EReal, EReal) + dot(EImag, EImag));
```

### 9.2 PhiImag correction

Equation:

```math
\nabla \cdot (\sigma \nabla \phi_i)
=
-\nabla \cdot (\omega \sigma A_r)
```

Update:

```cpp
EImag = -GradPhiImag - omega * ASrcReal;
JImag = sigma * EImag;
JouleHeatRaw = 0.5 * sigma * dot(EImag, EImag);
```

First version support:

```text
--no-phi
--phi-solver=GMRES
--phi-solver=PCG
```

### 9.3 divJ diagnostics

Must output:

```text
divJ_level0
divJ_phi
divJ_reduction = divJ_level0 / divJ_phi
```

Initial acceptance:

```text
divJ_phi < divJ_level0
```

Stable target:

```text
divJ_phi < 0.3 * divJ_level0
```

---

## 10. Power Scaling

### 10.1 Logic

Compute raw power:

```math
P_{raw} = \sum_i Q_i V_i
```

If target power specified:

```math
power\_scale = P_{target}/P_{raw}
```

```math
field\_scale = \sqrt{power\_scale}
```

Scaling:

```text
ASrc, BSrc, Phi, GradPhi, E, J *= field_scale
JouleHeat = JouleHeatRaw * power_scale
```

If `--no-power-scaling`:

```text
field_scale = 1
power_scale = 1
JouleHeat = JouleHeatRaw
```

### 10.2 Log

Must print:

```text
P_raw
P_target
power_scale
field_scale
effective_ampere_turns = ampere_turns * field_scale
min/max/mean JouleHeat
```

---

## 11. VTP Output Fields

GlassBody:

```text
Sigma
ASrcReal
ASrcImag
BSrcReal
BSrcImag
PhiImag
GradPhiImag
EReal
EImag
JReal
JImag
JouleHeatRaw
JouleHeat
DivJImag
```

Visual body:

```text
Position only, or BodyID
```

---

## 12. Analytic Tests

Before French reduced case is used as customer demo, add:

### 12.1 Single loop on-axis Bz

```text
test_3d_ophelie_biot_savart_circular_loop_axis
```

Analytic:

```math
B_z(z)=
\frac{\mu_0 I R^2}{2(R^2+z^2)^{3/2}}
```

Acceptance:

```text
relative error < 1% when segments_per_loop is large enough
```

### 12.2 Multiloop on-axis Bz

```text
test_3d_ophelie_biot_savart_multiloop_axis
```

Analytic value is sum of single-loop formulas.

### 12.3 Source scaling

```text
test_3d_ophelie_source_scaling
```

Acceptance:

```text
I -> 2I:
    A/B/E/J -> 2x
    raw Joule power -> 4x
```

### 12.4 Uniform A Level0

```text
test_3d_ophelie_level0_uniform_A
```

Set constant ASrcReal; check E/J/Q analytic values.

---

## 13. Recommended Run Commands

### 13.1 no-phi smoke

```bash
cd /home/yyc/SPHinXsysSYCL/build

./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced \
  --glass-radius=0.325 \
  --glass-height=0.50 \
  --dp=0.02 \
  --coil-radius=0.40 \
  --coil-z-min=0.05 \
  --coil-z-max=0.45 \
  --coil-num-loops=8 \
  --coil-segments-per-loop=256 \
  --frequency=300000 \
  --sigma=16 \
  --target-power=50000 \
  --no-phi \
  --state_recording=1
```

### 13.2 no-scaling analytic run

```bash
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced \
  --dp=0.02 \
  --coil-num-loops=1 \
  --coil-segments-per-loop=512 \
  --ampere-turns=1 \
  --no-power-scaling \
  --no-phi \
  --state_recording=1
```

### 13.3 phi run

```bash
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced \
  --dp=0.02 \
  --coil-num-loops=8 \
  --coil-segments-per-loop=256 \
  --target-power=50000 \
  --phi-solver=GMRES \
  --state_recording=1
```

---

## 14. README Must State Model Assumptions

`test_3d_ophelie_french_reduced/README.md` should state:

```text
This is a reduced OPHELIE-like French cold-crucible induction case.

It models:
    cylindrical conductive glass;
    prescribed multi-loop circular coil source;
    Biot-Savart A/B;
    optional PhiImag correction;
    E/J/JouleHeat;
    optional power normalization.

It does not model:
    segmented cold crucible EM;
    metal skin depth;
    cooling water;
    sigma(T);
    thermal feedback;
    stirrer;
    full PDE A-phi.
```

Summary:

```text
This case is a reduced geometry/physics version of the French cold-crucible literature:
only external induction coil, conductive glass body, Biot-Savart fields, in-glass scalar potential correction, and Joule heat source are preserved.
Crucible wall first version is visualization or future thermal boundary only, not in EM solve.
```

---

## 15. Cursor Current Task Manifest

### P0: Must implement

```text
1. Create test_3d_ophelie_french_reduced;
2. Generate GlassBody with analytic/code-generated cylinder;
3. Implement multi circular loop line source;
4. Implement line-source Biot on GlassBody;
5. Run no-phi E/J/Q end-to-end;
6. Support --target-power and --no-power-scaling;
7. Output VTP;
8. Write README;
9. Add circular-loop analytic Biot test.
```

### P1: Implement soon

```text
1. PhiImag correction;
2. divJ diagnostics;
3. JouleHeatRaw / JouleHeat scaled logic;
4. multi-loop analytic test;
5. source scaling test;
6. optional CoilVisualBody;
7. optional CrucibleWallVisualBody.
```

### P2: Later

```text
1. JouleHeat -> thermal one-way test;
2. temperature field / heat diffusion;
3. replace with real French geometry;
4. optional skull mask;
5. optional stirrer geometry;
6. optional passive crucible EM.
```

---

## 16. Final Implementation Principles

First-version French reduced geometry:

```text
Must:
    analytic GlassBody cylinder
    code-generated multi circular loop CoilSource
    Biot-Savart
    E/J/Q
    power scaling
    VTP

Should:
    PhiImag correction
    divJ diagnostics
    analytic Biot tests

Optional:
    CoilVisualBody
    CrucibleWallVisualBody

Do not:
    STL
    relax
    segmented crucible
    shell
    water cooling
    sigma(T)
    self-induction
    full A-phi
```

Final goal before real French geometry is ready: lock the EM chain:

```text
source current moment
    -> Biot-Savart A/B
        -> phi correction
            -> E/J
                -> JouleHeat
```

When real French geometry arrives, only replace:

```text
GlassBody shape
CoilSource loop distribution
Visual wall geometry
```

without rebuilding the EM solve chain.
