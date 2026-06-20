# OPHELIE/RH200 Single-Phase Glass EM–Joule–Stirring Midterm Delivery Task Plan

## 0. Direction Decision for This Round

The main line no longer continues TEAM7 boundary patch. TEAM7 is kept as a benchmark side track; midterm delivery main line switches to:

```text
RH200 stirring geometry + OPHELIE/French-style analytic equivalent coil + single-phase glass flow + Joule heat + thermal diffusion
```

This case is a **midterm coupled prototype**, not strict reproduction of French literature original geometry. In midterm defense, state clearly:

- RH200 STL provides stirring vessel, wall, and rotor geometry;
- OPHELIE/French reduced solver provides complex edge-flux EM solve, Joule heat, and target-power scaling;
- Coil is not STL but analytic multiloop/equivalent coil auto-generated from RH200 fluid domain size;
- Water cooling, external heat loss, segmented metal crucible, TEAM7 strict validation, and self-induction Picard main line are out of scope for now.

---

## 1. Geometry and Coil Strategy

### 1.1 Must use RH200 geometry

This stage uses RH200 STL as computational geometry:

```text
rotor.stl
oil/fluid domain STL, e.g. oil.stl or full_oil.stl
wall/container STL
```

Rename `oil` to `glass` in the new case. Do not keep grain/two-phase/injection logic.

### 1.2 Cursor must measure RH200 fluid domain size first

Before coil setup, audit RH200 glass/fluid domain dimensions. At minimum output:

```text
bbox_min = (xmin, ymin, zmin)
bbox_max = (xmax, ymax, zmax)
bbox_size = (Lx, Ly, Lz)
bbox_center = (cx, cy, cz)
radial_extent_xy = max sqrt((x-cx)^2 + (y-cy)^2)
fluid_particle_count
geometry_scale_used
```

Recommended two bboxes:

1. **STL mesh bbox**: check STL units and geometry size;
2. **relaxed/reloaded glass particle bbox**: used for actual coil generation.

Both must go to log and CSV, e.g.:

```text
rh200_geometry_audit.csv
rh200_coil_geometry_log.csv
```

### 1.3 Unit check is very important

Cursor must check STL units:

```text
If bbox max dimension is O(0.1) or O(1), likely already in m;
If bbox max dimension is O(100) or O(1000), likely mm; need geometry_scale=0.001.
```

Add/preserve CLI:

```bash
--geometry-scale=1.0
--geometry-scale=0.001
```

Log must state:

```text
raw STL bbox
applied geometry scale
physical bbox after scaling
```

### 1.4 Coil size auto-generated from RH200 fluid domain

RH200 has no coil STL. Midterm version uses analytic multiloop/equivalent coil. Coil geometry from glass/fluid bbox, default z-axis:

```text
coil_center = glass_bbox_center
coil_axis   = z-axis
R_fluid     = radial_extent_xy
H_fluid     = zmax - zmin
R_coil      = coil_radius_factor * R_fluid
z_coil_min  = zmin + z_margin_factor_low  * H_fluid
z_coil_max  = zmax - z_margin_factor_high * H_fluid
n_turns     = default 6 or 8
frequency   = 300 kHz
target_power = 50000 W
sigma_initial = 16 S/m
```

Recommended defaults:

```text
coil_radius_factor = 1.15
z_margin_factor_low = 0.05
z_margin_factor_high = 0.05
n_turns = 6 or 8
frequency = 300000 Hz
sigma_initial = 16 S/m
target_power = 50000 W
```

Must provide CLI override:

```bash
--coil-radius-factor=1.15
--coil-turns=8
--coil-z-margin-low=0.05
--coil-z-margin-high=0.05
--frequency=300000
--sigma0=16
--target-power=50000
```

Target-power scaling calibrates total power; first version need not match real coil current exactly. But coil position, radius, height, and turn count must be fully logged or Joule heat distribution cannot be explained.

---

## 2. Thermal Boundary Strategy: No Heat Loss Out of Domain

This stage excludes water cooling and heat leaving the computational domain. Default thermal boundary: adiabatic / no external heat loss.

First-version thermal model only needs:

```text
Glass internal thermal diffusion;
Joule heat temperature source term (later step);
Glass particles advect temperature with flow;
Walls: flow boundary only, **no** glass–wall thermal contact (adiabatic walls);
Rotor: keep **glass–rotor** thermal contact (user confirmed).
```

Default do not:

```text
wall fixed temperature;
wall Newton cooling;
glass-wall heat loss;
glass-rotor heat loss;
water cooling;
external heat sink.
```

If existing RH200 thermal diffusion template has oil-wall thermal contact, **disable wall thermal contact** in new case; **keep glass–rotor thermal contact** (user confirmed).

First-version energy budget:

\[
\Delta E_{thermal} \approx \int P_{joule}\,dt
\]

No external loss term. CSV may keep `P_wall_loss=0` or `external_heat_loss=0`; do not introduce actual heat loss.

---

## 3. New Case Recommended Structure

New standalone case; do not modify original RH200 two-phase case:

```text
tests/extra_source_and_tests/3d_examples/test_3d_ophelie_rh200_glass_em_stirring/
  test_3d_ophelie_rh200_glass_em_stirring.cpp
  CMakeLists.txt
  README.md
  data/
    rotor.stl
    oil.stl or glass.stl
    full_oil.stl or full_glass.stl
    wall/container STL if needed
```

CMake must copy `data/` to executable run directory `input/`:

```text
source data/*.stl -> binary input/*.stl
```

Recommended single-program three-stage CLI:

```bash
# Stage 0: particle relax
./test_3d_ophelie_rh200_glass_em_stirring --relax=1 --dp=0.005 --geometry-scale=1.0

# Stage 1: reload, one EM solve, write JouleHeat
./test_3d_ophelie_rh200_glass_em_stirring --reload=1 --em-solve=1 --target-power=50000

# Stage 2: fixed/periodic JouleHeat + single-phase stirring flow + thermal diffusion
./test_3d_ophelie_rh200_glass_em_stirring --reload=1 --run=1 --end-time=60
```

Later split relax target and run target; midterm demo uses single three-stage program to reduce body name/reload/variable registration mismatch risk.

---

## 4. Implementation Steps and Acceptance

### Step 0: Single-phase RH200 glass stirring + thermal diffusion

Borrow from `test_3d_stirring_injecting_RH200.cpp`:

```text
STL geometry loading;
wall / rotor / fluid body;
fluid acoustic step;
fluid-wall and fluid-rotor contact;
Simbody rotor rotation or direct imposed rotor motion;
Temperature registration;
thermal diffusion framework;
VTP output.
```

Delete / do not borrow:

```text
GrainFluidGeometry;
grain_f body;
grain_f emitter;
grain_f inflow condition;
oil-grain contact;
grain-wall contact;
grain thermal relaxation;
all two-phase injection logic.
```

Rename:

```text
oil -> glass
OilGeometry -> GlassGeometry
rho0_f -> rho0_glass
mu_f / nu_f -> mu_glass / nu_glass
cp_oil -> cp_glass
k_oil -> k_glass
initial_temperature_oil -> initial_temperature_glass
```

First-version glass parameters: stable runnable values; do not use extreme viscosity that crashes flow step. Tune toward real molten glass later.

Acceptance:

```text
1. --relax=1 generates reload;
2. --reload=1 --run=1 runs;
3. rotor rotates;
4. glass velocity/pressure/temperature in VTP;
5. no GrainFluid/inlet/two-phase output;
6. temperature field stable with particle motion output.
```

---

### Step 1: Fake JouleHeat; validate heat source and energy closure first

Before OPHELIE EM, fake heat source to validate source term and energy budget.

Add CLI:

```bash
--joule-mode=off|uniform|analytic-position
--fake-joule-power=50000
```

Temperature source:

\[
T_i^{n+1}=T_i^n+\Delta t\frac{Q_i}{\rho_i c_p}
\]

Where `Q_i` is volumetric power density W/m³.

Two fake sources:

1. `uniform`: same power density per glass particle for strict energy closure;
2. `analytic-position`: fixed spatial distribution from current particle position, e.g. hotter at center or near coil, for stirring transport demo.

Note: real EM heat is spatial source fixed by coil field. Non-uniform fake heat that moves permanently with particles becomes material-bound source; less physical than `Q(x_current)`. Therefore `analytic-position` should update from current particle position.

Acceptance:

```text
1. P_joule_total near fake-joule-power;
2. thermal_energy increment ≈ ∫P_joule dt;
3. no external heat loss;
4. T_mean/T_max rise over time;
5. velocity field transports temperature field.
```

---

### Step 2: Connect OPHELIE EM one-shot solve

Register OPHELIE variables on `FluidBody glass`. If existing interface is SolidBody-only, generalize to `SPHBody&` or `BaseParticles&`.

Minimum glass particle variables:

```text
Sigma
PhiReal, PhiImag
ACoilReal, ACoilImag
AIndReal, AIndImag optional
ATotalReal, ATotalImag
EReal, EImag
JReal, JImag
JouleHeat
JouleHeatRaw
JouleHeatScaled
```

First-version EM mode:

```bash
--ophelie-em-mode=off|fixed-joule|solve-once|update-sigma
```

Midterm default path:

```text
solve-once or fixed-joule
```

Flow:

```text
1. reload glass particles;
2. measure glass bbox;
3. generate equivalent multiloop coil from bbox;
4. compute A_coil on glass particles;
5. solve complex edge-flux phi;
6. reconstruct E/J;
7. compute JouleHeat = 0.5*sigma*(|EReal|^2 + |EImag|^2);
8. scale to target_power;
9. write JouleHeat to glass particles;
10. output VTP and CSV.
```

Acceptance:

```text
1. VTP outputs |E|, |J|, JouleHeat, Sigma;
2. JouleHeat >= 0, no NaN/Inf;
3. P_joule_total = target_power within tolerance;
4. coil geometry log complete;
5. glass bbox log complete;
6. JouleHeat spatial distribution visually reasonable.
```

---

### Step 3: JouleHeat + stirring flow + thermal diffusion main loop

Recommended main loop order:

```text
1. fluid acoustic / transport step;
2. rotor motion / Simbody constraint update;
3. apply Joule heat temperature source;
4. glass inner thermal diffusion;
5. update cell linked list / configurations as required;
6. output VTP/CSV.
```

First version may use fixed JouleHeat from initial EM solve. More physical: periodic update:

```bash
--em-update-interval=0.1
```

Every interval:

```text
1. recompute A_coil/E/J/JouleHeat from current particle positions;
2. if sigma(T) off, update spatial field only;
3. if sigma(T) on, update sigma too;
4. target-power scaling again.
```

Midterm first version may disable sigma(T); preserve CLI and code structure.

Acceptance:

```text
1. T_mean/T_max rise over time;
2. velocity/stirring changes temperature spatial distribution;
3. energy budget closed;
4. VTP shows Temperature, Velocity, JouleHeat together;
5. no external heat loss default holds.
```

---

### Step 4: sigma(T) weak coupling, optional later

Not blocking midterm defense. Optional:

```bash
--sigma-temperature-coupling=0|1
--sigma-relaxation-alpha=0.2
```

Weak coupling:

\[
\sigma_{new}=(1-\alpha)\sigma_{old}+\alpha\sigma(T)
\]

Each EM update output:

```text
sigma_min
sigma_max
sigma_mean
P_joule_total_before_scaling
P_joule_total_after_scaling
T_mean
T_max
```

Do not make Picard self-induction midterm main path.

---

## 5. Output Requirements

### 5.1 VTP output

Minimum:

```text
Temperature
Velocity
Pressure
JouleHeat
Sigma
|E|
|J|
```

If possible, components:

```text
EReal, EImag
JReal, JImag
ACoilReal, ACoilImag
PhiReal, PhiImag
```

### 5.2 CSV output

Recommended files:

```text
rh200_geometry_audit.csv
rh200_coil_geometry_log.csv
rh200_em_joule_summary.csv
rh200_coupled_energy_history.csv
```

`rh200_coupled_energy_history.csv` at minimum:

```text
time_s
P_joule_total_W
T_mean_K
T_min_K
T_max_K
thermal_energy_J
thermal_energy_increment_J
external_heat_loss_W
sigma_min_S_per_m
sigma_max_S_per_m
sigma_mean_S_per_m
```

With no external loss:

```text
external_heat_loss_W = 0
```

### 5.3 README must state

```text
1. This is a midterm coupled prototype.
2. RH200 STL provides the stirring vessel and rotor geometry.
3. The induction coil is an equivalent analytic multiloop coil generated from the RH200 glass-domain bbox.
4. No water cooling or external thermal loss is included in the default setup.
5. TEAM7 validation and full French literature geometry reproduction are not part of this case.
6. Default thermal boundary is adiabatic.
```

---

## 6. Explicitly Out of Scope This Stage

```text
1. TEAM7 boundary/no-flux patch;
2. TEAM7 strict validation;
3. segmented cold-crucible metal wall EM;
4. water cooling;
5. external heat loss;
6. two-phase grain injection;
7. strong sigma(T)-EM Picard coupling;
8. self-induction Picard as default;
9. empirical J scaling / Cij scaling / RHS scaling;
10. a_sign hack.
```

---

## 7. Midterm Defense Recommended Presentation

Prepare at least four figure groups:

```text
1. RH200 geometry + equivalent coil layout;
2. EM results: |E|, |J|, JouleHeat;
3. Coupled thermal-flow: Temperature at several times + Velocity/glyph;
4. Energy history: P_joule_total, T_mean/T_max, thermal_energy.
```

Narrative focus:

```text
1. OPHELIE-inspired complex edge-flux EM solver integrated on RH200 stirring geometry;
2. Target-power Joule heat deposition implemented;
3. Single-phase glass stirring flow and thermal diffusion implemented;
4. Default adiabatic boundary; no water cooling or external heat loss;
5. Later: sigma(T), more realistic coil parameters, full French/cold-crucible geometry.
```

---

## 8. Cursor Execution Priority

Strict order:

```text
Priority 1: Step 0 single-phase RH200 stirring + thermal diffusion running.
Priority 2: Step 1 fake JouleHeat + energy closure.
Priority 3: RH200 bbox measurement + equivalent coil auto-generation and logging.
Priority 4: Step 2 OPHELIE EM solve once + target power JouleHeat.
Priority 5: Step 3 fixed/periodic JouleHeat + flow/thermal coupling.
Priority 6: sigma(T) weak coupling optional.
```

Only after Step 0 and Step 1 pass, connect real OPHELIE EM. Do not debug RH200 stirring, EM, sigma(T), and energy budget all at once from the start.
