# OPHELIE-like TEAM7 Coil Source Modeling: Issue Summary and Next-Steps Implementation Plan

This document consolidates the last two rounds of discussion. The goal is to give Cursor a clear development roadmap for the next stage. The core issue is: the current OPHELIE-like / Biot-Savart branch can already run TEAM7 native STL reload, coil source, Biot, no-phi Level0, and Bz probe comparison, but the Bz line still has large field-shape error versus TEAM7 reference. Next steps should not blindly run more TEAM7 cases; coil excitation/source model must be clarified and corrected first.

---

## 1. Current Status Summary

Current native TEAM7 STL:

```text
coil.stl:
    x = 94–294 mm
    y = 0–200 mm
    z = 49–149 mm

plate.stl:
    x = 0–294 mm
    y = 0–294 mm
    z = 0–19 mm
```

A1–B1 reference probe line:

```text
(x, y, z) = (0..288 mm, 72 mm, 34 mm)
```

This line lies:

```text
above the plate top surface;
below the coil bottom surface;
i.e. in the air region between plate and coil.
```

Current test pipeline:

```text
native STL + relax/reload
    -> CoilSourceBody / PlateBody particle loading
        -> initialize Jsrc on coil
            -> coil -> probe / coil -> plate Biot-Savart
                -> no-phi Level0: E = -i omega A, J = sigma E
                    -> optional power scaling
                        -> Bz probe vs TEAM7 reference CSV
```

Current key results:

```text
dp = 6 mm reload:
    n_plate ≈ 7700
    n_coil ≈ 7712

no-phi, no-power-scaling, volume J = J0 e_theta:
    full-line RMS relative error ≈ 66.5%
    reference Bz peak location: x ≈ 126 mm
    simulation Bz peak location: x ≈ 198 mm
    single point at x=126 mm: sim/ref ≈ 8.95/7.81 mT, error ~15%
    x=162–234 mm region: simulation clearly overshoots
    dp=6 -> dp=3: RMS basically unchanged
```

Conclusion:

```text
1. Not a "can the code run" issue; smoke already passed;
2. Not mainly a particle-resolution issue; finer dp does not help;
3. Not a total ampere-turns off-by-order-of-magnitude issue; single-point amplitude at x=126 is already close;
4. Main cause is more likely wrong spatial distribution / current-path model of coil source.
```

---

## 2. Issues with the Current Coil Source Model

Current approach:

```math
\mathbf{J}_{src} = J_0 \mathbf{e}_\theta
```

Where:

```math
\mathbf{e}_\theta =
\left(
-\frac{y-y_c}{r},
\frac{x-x_c}{r},
0
\right)
```

This is equivalent to assuming the coil is a ring/axisymmetric current source rotated about some z-axis.  
But TEAM7 coil is closer to a closed multi-turn / stranded coil topology, i.e. approximately rectangular or racetrack closed winding, not circular current about the z-axis.

Current simulation peak lies at:

```text
x ≈ 198 mm
```

Very close to the coil geometry centroid side, which matches expectation for a uniform volume e_theta source, but does not match reference peak near:

```text
x ≈ 126 mm
```

Therefore, next steps should not simply keep tuning `sigma`, `dp`, or power scaling; new coil source models should be tested.

---

## 3. Whether Surface Particles / Shell Are Needed

### 3.1 TEAM7 Should Not Be Prioritized as Skin-Depth Surface Current

TEAM7 is a low-frequency stranded / multi-turn coil benchmark. Its core is:

```text
closed winding path;
2742 turns;
1 A/turn;
50 Hz / 200 Hz;
coil current direction determined by closed coil topology.
```

This is not the same issue as high-frequency metal skin-depth surface treatment in French cold-crucible literature.

So:

```text
TEAM7:
    primary issue is coil current path / topology;
    should not first pivot to surface shell + e_theta.

French cold crucible:
    high-frequency metal parts may use surface / shell / skin-depth;
    but current client simplified version can first use prescribed coil source + glass volume.
```

### 3.2 Surface Particles Can Be Diagnostic, but Not the Main Next-Step Line

Can add later:

```text
--coil-source-model=surface-etheta
--coil-source-model=surface-racetrack
```

But not recommended as next-step main line: "whole coil outer surface 6 faces + e_theta".

Reasons:

```text
1. e_theta direction itself may be wrong;
2. only moves error from volume to surface, not necessarily improving peak location;
3. surface particle measure may be area, but internally may still be called VolumetricMeasure, units easy to get wrong;
4. already tried solid radial outer-shell truncation, RMS worse; "moving current outward" is not a sufficient fix.
```

---

## 4. Surface Treatment in French OPHELIE Literature vs. This Project Goal

In French OPHELIE:

```text
coil / crucible / base metal parts:
    at high frequency skin depth is very small;
    can use surface mesh + skin-depth exponential decay;

molten glass:
    thicker skin depth;
    needs volume mesh.
```

This does not mean "all coil excitation should only be applied on outer surface".

More accurately:

```text
French high-frequency cold crucible:
    metal conductors can use surface reduction;
    glass uses volume integral;
    crucible is passive induced conductor, not active excitation source.

TEAM7 stranded coil:
    should prioritize building closed coil path;
    whether surface is a later reduction choice, not the first issue.
```

---

## 5. Cold-Crucible Wall Current Should Not Be Applied Directly

Must clearly distinguish:

```text
coil:
    active source;
    can apply prescribed current / ampere-turns.

cold crucible wall:
    passive conductor;
    if EM is included later, induced current should be solved;
    should not apply excitation current directly.
```

Current client delivery first version can continue:

```text
coil = prescribed source
glass = conducting volume
crucible wall = geometry / thermal boundary only
```

Later enhanced version:

```text
crucible segments = passive conducting shell
induced j solved by OPHELIE-like integral equation
```

Do not apply current directly to cold crucible wall as source; physics would be wrong.

---

## 6. Full A–phi Solver Faces the Same Issue

This issue is not OPHELIE/Biot-specific. Full A–phi solver must also handle coil excitation correctly:

```text
J_src must be a closed loop;
direction must match coil topology;
total ampere-turns must be correct;
source should preferably be approximately divergence-free;
cannot arbitrarily replace complex stranded coil with J0 e_theta.
```

OPHELIE/Biot branch exposes this issue more easily because Bz probe directly reflects source error.  
Clarifying TEAM7 source model on OPHELIE branch first is valuable for later A–phi main line.

---

# 7. Core Implementation Plan: From e_theta to Racetrack Source Model

Next-step recommendation: do not start with shell; introduce multiple coil source models and quickly screen with reference Bz line.

Recommended CLI additions:

```text
--coil-source-model=volume-etheta
--coil-source-model=filament-racetrack
--coil-source-model=volume-racetrack
--coil-source-model=surface-etheta      # later diagnostic
--coil-source-model=surface-racetrack   # later diagnostic
```

Priority:

```text
P0: filament-racetrack
P1: volume-racetrack
P2: surface-racetrack
P3: surface-etheta diagnostic only
```

---

## 8. filament-racetrack: Minimal Diagnostic Source Model

### 8.1 Why filament-racetrack First

Advantages:

```text
1. no coil particles needed;
2. no J0 computation needed;
3. no surface/volume measure handling;
4. only needs one closed polyline;
5. fastest way to test whether current-path direction is main cause of Bz peak shift.
```

### 8.2 How to Generate Rectangular Closed Segments

Use coil STL bbox to generate rectangular closed centerline. Current bbox:

```text
xmin = 94 mm
xmax = 294 mm
ymin = 0 mm
ymax = 200 mm
zmin = 49 mm
zmax = 149 mm
```

Default center plane:

```text
z_mid = 99 mm
```

Four corners:

```text
p0 = (94, 0, 99) mm
p1 = (294, 0, 99) mm
p2 = (294, 200, 99) mm
p3 = (94, 200, 99) mm
```

Convert to SI:

```text
p *= 1e-3
```

Segments:

```text
p0 -> p1
p1 -> p2
p2 -> p3
p3 -> p0
```

Supported parameters:

```text
--racetrack-inset-mm=0/10/20/30
--racetrack-z-mm=49/99/149
--racetrack-ds-mm=1/2/3/6
```

With inset:

```text
xmin = 94 + inset
xmax = 294 - inset
ymin = 0 + inset
ymax = 200 - inset
```

### 8.3 Whether Particles Are Needed

Not needed.

Directly subdivide each segment into small integration segments, generate source quadrature points:

```cpp
struct LineSourcePoint
{
    Vecd pos;
    Vecd dl;   // directed, unit m
};
```

For each segment `[a,b]`:

```cpp
Vecd ab = b - a;
Real L = norm(ab);
int n = ceil(L / ds);
Vecd dl = ab / n;

for (int m = 0; m < n; ++m)
{
    pos = a + (m + 0.5) * dl;
    source_pos.push_back(pos);
    source_dl.push_back(dl);
}
```

### 8.4 Biot-Savart Formula

Total ampere-turns:

```math
I_{\mathrm{eff}} = N I
```

For TEAM7:

```text
N = 2742
I = 1 A/turn
I_eff = 2742 A
```

Magnetic vector potential:

```math
\mathbf{A}(x_i)
=
\frac{\mu_0}{4\pi}
I_{\mathrm{eff}}
\sum_k
\frac{\Delta \mathbf{l}_k}
{|x_i-x_k|}
```

Magnetic field:

```math
\mathbf{B}(x_i)
=
\frac{\mu_0}{4\pi}
I_{\mathrm{eff}}
\sum_k
\frac{
\Delta \mathbf{l}_k \times (x_i-x_k)
}
{|x_i-x_k|^3}
```

Code pseudocode:

```cpp
const Real Ieff = number_of_turns * current_amplitude;
const Real coeff = mu0 / (4.0 * Pi);

for each probe or plate particle i:
    Vecd A = Zero;
    Vecd B = Zero;

    for each line source point k:
        Vecd r = x_i - source_pos[k];
        Real r2 = dot(r, r) + eps2;
        Real inv_r = 1.0 / sqrt(r2);
        Real inv_r3 = inv_r / r2;

        Vecd Idl = Ieff * source_dl[k];

        A += coeff * Idl * inv_r;
        B += coeff * cross(Idl, r) * inv_r3;
```

### 8.5 filament-racetrack Acceptance Criteria

Run:

```text
--coil-source-model=filament-racetrack
--no-phi
--no-power-scaling
--compare-team7-bz
```

Check:

```text
1. does Bz peak move from x≈198 toward x≈126;
2. does RMS drop clearly from 66%;
3. is overshoot in x=162–234 relieved;
4. is negative lobe in x=0–72 improved.
```

If filament-racetrack shows clear improvement, then current path / topology is likely the main cause.

If peak still at x≈198, then possible causes:

```text
1. reference includes plate-induced total B;
2. racetrack size/position wrong;
3. coil winding path differs from simple rectangle;
4. coordinate or reference interpretation still has issues.
```

---

# 9. volume-racetrack: Assign Path Direction to Coil Particles

If filament-racetrack proves path is main cause, then do volume-racetrack.

## 9.1 Basic Idea

Still use reloaded coil solid particles, but no longer use:

```math
\mathbf{J} = J_0 \mathbf{e}_\theta
```

Change to:

```math
\mathbf{J} = J_0 \mathbf{t}_{coil}(x)
```

Where \(t_{coil}(x)\) is tangent direction of nearest racetrack segment.

## 9.2 How to Get Tangent from Particle Position

Define segment:

```cpp
struct RacetrackSegment
{
    Vecd a;
    Vecd b;
    Vecd tangent;  // normalize(b-a)
    Real length;
};
```

For each coil particle `x`:

```cpp
Real best_d2 = Infinity;
Vecd best_t = Zero;

for each segment s:
{
    Vecd ab = s.b - s.a;
    Real L2 = dot(ab, ab);

    Real u = dot(x - s.a, ab) / L2;
    u = clamp(u, 0.0, 1.0);

    Vecd q = s.a + u * ab;       // nearest point on segment
    Real d2 = dot(x - q, x - q);

    if (d2 < best_d2)
    {
        best_d2 = d2;
        best_t = normalize(ab);
    }
}

Jsrc[x] = J0 * best_t;
```

This is "assign closed-winding tangent direction based on particle position relative to racetrack coil path".

## 9.3 volume-racetrack J0

For closed-path homogenized volume source:

```math
A_{cross} \approx \frac{V_{coil}}{L_{path}}
```

Therefore:

```math
J_0
=
\frac{NI}{A_{cross}}
=
\frac{NI L_{path}}{V_{coil}}
```

Where:

```text
V_coil = sum_j Vol_j
L_path = sum segment lengths
```

This definition fits rectangular/racetrack coil better than:

```math
A_{cross} = \sum_j \frac{V_j}{2\pi r_j}
```

The latter only suits ring/near-ring \(e_\theta\) source model.

## 9.4 volume-racetrack Acceptance Criteria

Same run:

```text
--coil-source-model=volume-racetrack
--no-phi
--no-power-scaling
--compare-team7-bz
```

Compare:

```text
volume-etheta
filament-racetrack
volume-racetrack
```

Focus on:

```text
peak_x
RMS rel error
x=126
x=198
x=162–234 over-shoot
x=0–72 negative lobe
```

---

# 10. Surface Particles / Shell Plan Comes Later

If shell is needed, note two things.

## 10.1 Source Weight Units

Biot-Savart source moment must be unified as:

```text
current moment = A·m
```

Weights for three source types:

```text
volume current:
    J [A/m²] * dV [m³] = A·m

surface current:
    K [A/m] * dS [m²] = A·m

filament current:
    I [A] * dl [m] = A·m
```

If using surface particles, must confirm SPHinXsys particle measure is `dS` or `dS*thickness`. Cannot directly reuse volume source `J * Vol` logic.

Recommend abstract unified field:

```cpp
CurrentMomentReal = source_direction * source_weight;
```

Biot-Savart only sums:

```cpp
A += coeff * CurrentMomentReal / r;
B += coeff * cross(CurrentMomentReal, r) / r^3;
```

## 10.2 Shell Direction Should Also Prioritize Racetrack

If doing shell, recommend first:

```text
surface-racetrack
```

Not:

```text
surface-etheta
```

Because current main issue may be current path, not volume vs surface.

---

# 11. Cursor Next-Step Development Tasks

Can send the following tasks directly to Cursor.

## Task A: Add coil source model enum / CLI

Add:

```cpp
enum class OphelieCoilSourceModel
{
    VolumeETheta,
    FilamentRacetrack,
    VolumeRacetrack,
    SurfaceETheta,      // later
    SurfaceRacetrack    // later
};
```

CLI:

```text
--coil-source-model=volume-etheta
--coil-source-model=filament-racetrack
--coil-source-model=volume-racetrack
--racetrack-inset-mm=0
--racetrack-z-mm=99
--racetrack-ds-mm=2
```

Default keep:

```text
volume-etheta
```

for regression comparison.

## Task B: Implement racetrack geometry helper

Recommended new files:

```text
electromagnetic_ophelie_racetrack_source.h
electromagnetic_ophelie_racetrack_source.hpp
```

Core functions:

```cpp
struct RacetrackSegment
{
    Vecd a;
    Vecd b;
    Vecd tangent;
    Real length;
};

StdVec<RacetrackSegment> makeRectangularRacetrack(
    Real xmin, Real xmax,
    Real ymin, Real ymax,
    Real z);

Vecd nearestRacetrackTangent(
    const Vecd& x,
    const StdVec<RacetrackSegment>& segments);

void buildFilamentQuadrature(
    const StdVec<RacetrackSegment>& segments,
    Real ds,
    StdVec<Vecd>& source_pos,
    StdVec<Vecd>& source_dl);
```

## Task C: Implement filament-racetrack probe Biot

No need to change plate field first; only support A1–B1 probe comparison.

Implement:

```text
--coil-source-model=filament-racetrack
--compare-team7-bz
```

On host side for probe points:

```cpp
B += coeff * Ieff * cross(dl, r) / r^3;
```

Output:

```text
x_mm
Bz_sim_mT
Bz_ref_mT
abs_err
rel_err
RMS_rel
peak_x_sim
peak_x_ref
```

## Task D: Implement volume-racetrack Jsrc

At coil particle initialization:

```cpp
if model == volume-racetrack:
    tangent = nearestRacetrackTangent(pos_i, segments);
    JsrcReal_i = J0 * tangent;
```

Where:

```cpp
J0 = N * I * L_path / V_coil;
```

And print:

```text
source_model
V_coil
L_path
A_cross = V_coil / L_path
J0
```

## Task E: Add comparison table

Each `--compare-team7-bz` output:

```text
source_model
RMS_rel
peak_x_sim
peak_x_ref
sim/ref at x=126
sim/ref at x=198
```

For quick comparison:

```text
volume-etheta
filament-racetrack
volume-racetrack
```

---

# 12. Recommended Run Matrix

First round:

```bash
./test_3d_ophelie_team7 \
  --native-stl --reload=1 \
  --no-phi --no-power-scaling \
  --compare-team7-bz \
  --coil-source-model=volume-etheta \
  --state_recording=0
```

Second round:

```bash
./test_3d_ophelie_team7 \
  --native-stl --reload=1 \
  --no-phi --no-power-scaling \
  --compare-team7-bz \
  --coil-source-model=filament-racetrack \
  --racetrack-inset-mm=0 \
  --racetrack-z-mm=99 \
  --racetrack-ds-mm=2 \
  --state_recording=0
```

Parameter sweep:

```text
--racetrack-inset-mm=0,10,20,30
--racetrack-z-mm=49,99,149
```

Third round:

```bash
./test_3d_ophelie_team7 \
  --native-stl --reload=1 \
  --no-phi --no-power-scaling \
  --compare-team7-bz \
  --coil-source-model=volume-racetrack \
  --racetrack-inset-mm=0 \
  --racetrack-z-mm=99 \
  --state_recording=0
```

Only after all of these, consider:

```text
surface-racetrack
phi / Jy
self-induction / phase90
200 Hz
```

---

# 13. Acceptance Criteria

Current P0 is not full TEAM7 validation, only source model selection.

Recommendation:

```text
RMS_rel > 60%:
    source model still clearly wrong.

RMS_rel 40–60%:
    acceptable as smoke, but not good enough.

RMS_rel 20–40%:
    source model basically iterable; enter volume-racetrack / phi / self-induction.

RMS_rel < 20%:
    coil source field already good.
```

Also check:

```text
peak_x_sim near 126 mm;
x=162–234 still severely overshooting;
x=0–72 negative lobe improved.
```

Do not look only at single point \(x=126\).

---

# 14. Final Recommendation

Should not make "shell surface current" the main line now.  
Minimum and highest-information next-step plan:

```text
1. Add filament-racetrack;
2. no particles, direct line-segment Biot;
3. compare with A1–B1 Bz reference;
4. if peak improves, then do volume-racetrack;
5. only then decide whether to do surface-racetrack.
```

In one sentence:

```text
First generate closed racetrack centerline, no shell particles;
use line integration first to test whether current path explains Bz peak shift;
then assign same path tangent to coil volume particles;
only then consider surface particles.
```
