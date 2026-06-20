# OPHELIE French Glass EM--Thermal--Stirring Midterm Delivery Plan

Date: 2026-06-18  
Purpose: redirect the near-term work from TEAM7 strict validation to a midterm customer-deliverable French-reduced glass-bath coupled case.

---

## 0. Executive decision

The project direction should be adjusted immediately for the midterm customer presentation.

The current TEAM7 line has been useful for diagnosing the high-conductivity metal benchmark, but the recent P5/P5.6-lite no-flux boundary work shows that continuing to stack boundary patches on the TEAM7 perforated-plate case has low short-term return. For the customer-facing midterm goal, the main line should switch to:

```text
French-reduced molten-glass geometry
+ OPHELIE-style complex edge-flux EM solve
+ Joule heat deposition in the glass region
+ single-phase stirring flow
+ thermal diffusion / heat transport
+ optional weak sigma(T) coupling
```

TEAM7 should remain a benchmark side track, not the delivery blocker.

The midterm target should not be framed as a full cold-crucible metallic-structure solver. It should be framed as a Phase-1 reduced OPHELIE-like glass-bath induction-heating solver.

---

## 1. What the recent TEAM7 no-flux work tells us

### 1.1 Correct physical boundary interpretation

The heat-transfer no-flux analogy is valid only if it is translated correctly.

For heat conduction, an adiabatic wall means:

\[
q_n=-k\nabla T\cdot n=0 .
\]

For the EM conductor boundary, the analogous condition is not \(B=0\), and not \(A=0\). It is the no-conduction-current boundary:

\[
n\cdot J = 0 .
\]

For an isotropic conductor, with

\[
J=\sigma E,
\]

this is equivalent to

\[
n\cdot E = 0,
\]

where the total electric field is

\[
E=-\nabla\phi-i\omega A .
\]

In the current real/imag edge-flux implementation, this translates into a zero boundary edge EMF condition, not a magnetic insulation condition.

Important clarification:

```text
A and B are not zero outside the conductor.
The coil and the conductor are magnetically coupled through A/B in air.
Only the conduction current J is no-flux across the conductor boundary.
```

---

### 1.2 What P5/P5.6-lite achieved

The P5 no-flux line did not fail conceptually. It narrowed down what does and does not work.

Completed work:

```text
P5.0  surface normal / signed distance reload: successful
P5.1  post-projection no-flux diagnostic: runs, but worsens EMF mismatch
P5.2  tangent-LS diagnostic: distance_t correction improves shell metrics but not TEAM7 core
P5.3  EMF component diagnostics: useful, clarifies previous interpretation
P5.4  ghost-edge no-flux reconstruction: small Jn/Jt improvement, no Bind/B improvement
P5.5  missing moment / phi-RHS ghost: essentially baseline, no useful gain
P5.6-lite box benchmark / L1 audit: diagnostic value, no TEAM7 L2 closure
```

Core TEAM7 L2 numbers remain essentially unchanged:

```text
boundary mode      Bind/Bcoil     hole e_edge_em     phase90 RMS
none               ~3.913         ~0.338             ~9.80
no-flux ghost      ~3.911         ~0.351             ~9.80
no-flux full       ~3.911         ~0.351             ~9.80
```

Interpretation:

```text
E/J reconstruction-level no-flux patch is not enough.
RHS-only ghost is not enough.
Current missing-moment implementation is not a true level-set/static-confinement correction.
A full phi-LHS / Neumann / static-confinement implementation may still be physically meaningful, but it is not worth doing directly on the TEAM7 perforated geometry before simpler benchmarks and L1 reference closure.
```

---

### 1.3 Current TEAM7 status

TEAM7 should now be downgraded to a benchmark side track.

Reasons:

```text
1. L1 source/reference/probe is still not closed.
   Current phase0 reference is not a pure source-only reference.
   Peak location mismatch remains large: x_sim roughly 198--216 mm vs x_ref roughly 126 mm.

2. L2 one-way remains diagnostic-only.
   It does not include self-consistent A_total = A_coil + A_ind feedback.

3. Boundary no-flux E/J patching does not improve the core metrics.

4. Continuing to patch TEAM7 hole geometry is unlikely to produce a short-term customer deliverable.
```

TEAM7 should continue only through narrow items:

```text
TEAM7 side-track task 1: L1 source/reference/probe geometry closure.
TEAM7 side-track task 2: box/slab phi-Neumann or static-confinement manufactured benchmark.
```

Do not continue broad TEAM7 patching now.

---

## 2. Midterm customer deliverable: correct scope

The customer-facing midterm result should be scoped as:

```text
Phase-1 OPHELIE-inspired SPH/SYCL electromagnetic--thermal--flow solver
for the molten-glass region of the French reduced geometry.
```

### 2.1 Include in the midterm demo

```text
1. French reduced glass-bath geometry.
2. Given coil source / impressed source model.
3. Complex edge-flux solve for phi/E/J in the glass body.
4. Joule heat density:
   Q = 0.5 sigma ( |E_real|^2 + |E_imag|^2 ).
5. Target-power scaling, e.g. 50 kW or specified customer power.
6. Single-phase glass stirring flow.
7. Thermal diffusion and Lagrangian thermal transport.
8. Joule heat source term coupled into temperature.
9. VTK output of |E|, |J|, JouleHeat, Temperature, Velocity, Pressure.
10. Energy-budget CSV.
11. Optional weak sigma(T) coupling.
```

### 2.2 Explicitly do not include in this midterm scope

```text
1. TEAM7 strict validation pass.
2. Metallic segmented crucible eddy currents.
3. Water-cooled wall EM response.
4. Full high-conductivity skin-effect benchmark.
5. Fully quantitative reproduction of all French natural-convection data.
6. Guaranteed converged self-induction Picard for every geometry.
```

This scope is realistic and defensible for a midterm answer to the customer.

---

## 3. Current distance to the target

### 3.1 For a midterm demonstration

Target:

```text
French geometry
+ EM Joule heat
+ single-phase stirring
+ thermal evolution
+ credible output and energy budget
```

Distance: medium. The main effort is engineering integration, not new theory.

Available pieces:

```text
1. Complex edge-flux solver exists.
2. JouleHeatEdgeReconComplex exists.
3. Solver-local normalization exists.
4. Target-power scaling has been used.
5. Thermal source / DeltaT energy closure has been tested before.
6. SPHinXsys single-/multi-body fluid and thermal modules exist.
7. Uploaded stirring case provides rotor/wall/thermal/contact structure.
```

Still needed:

```text
1. Clean single-phase glass stirring case.
2. EM solver call integrated into this case.
3. Joule heat source term on the glass particles.
4. Energy-budget output.
5. Optional sigma(T) update loop.
6. Stable run and VTK outputs for the midterm presentation.
```

Estimated development effort:

```text
First coupled demo: 3--5 focused working days if code integration is smooth.
More credible midterm version with energy budget and optional sigma(T): about 1--2 weeks.
Full French quantitative validation: still much farther away.
```

---

## 4. Uploaded stirring case: how it should be reused

The uploaded file is:

```text
test_3d_stirring_injecting_RH200.cpp
```

It is useful as a structural reference, but not as the final target case.

### 4.1 What the uploaded case currently contains

```text
1. Rotor solid body.
2. Oil fluid body.
3. GrainFluid second fluid body.
4. Grain inflow emitter.
5. Oil--Grain two-phase contact.
6. Wall boundary.
7. Rotor--fluid interaction.
8. Simbody-controlled rotation.
9. Temperature fields.
10. Contact heat diffusion.
11. Fluid acoustic/advection time stepping.
```

### 4.2 What should be removed for the midterm case

Remove or disable all two-phase injection logic:

```text
GrainFluidGeometry
grain_f body
grain_f material
grain_f inner/contact relations
grain_f emitter / buffer / inflow condition
grain_f_emitter_injection
oil_grain_f_contact
grain_f_oil_contact
grain_f_wall_contact
grain_f_thermal_relaxation
grain_f VTK / temperature output
```

### 4.3 What should be kept

Keep and adapt:

```text
Rotor / stirrer body
WallBoundary body
Single fluid body, renamed glass
Fluid-wall contact
Rotor-glass contact
Rotor imposed rotation or Simbody pin rotation
Single-phase SPH acoustic/advection stepping
Density regularization if stable
Thermal diffusion on glass
Body state output
Temperature / pressure / velocity output
```

---

## 5. New case Cursor should create

Do not directly mutate the old two-phase case. Create a new case:

```text
test_3d_ophelie_french_glass_stirring.cpp
```

Suggested case purpose:

```text
Single-phase molten-glass stirring case with OPHELIE Joule heating.
```

### 5.1 Body structure

Recommended bodies:

```cpp
FluidBody glass(system, makeShared<GlassGeometry>(...));
SolidBody wall_boundary(system, makeShared<WallBoundaryGeometry>(...));
SolidBody stirrer(system, makeShared<StirrerGeometry>(...));
```

Relations:

```cpp
InnerRelation glass_inner(glass);
ContactRelation glass_wall_contact(glass, {&wall_boundary, &stirrer});
ContactRelation stirrer_glass_contact(stirrer, {&glass});
```

If the French reduced geometry does not yet include an actual stirring impeller, use the rotor/stirrer from the uploaded case as a temporary mechanical stirring device for the midterm demo, and document that it is a demonstration stirring mechanism.

---

## 6. Single-phase glass material setup

Rename oil variables to glass variables.

Replace:

```text
oil -> glass
rho0_f -> rho0_glass
mu_f / nu_f -> mu_glass / nu_glass
cp_oil -> cp_glass
k_oil -> k_glass
T0_oil -> T0_glass
```

Minimum material fields needed:

```text
rho_glass
mu_glass or nu_glass
cp_glass
k_glass
sigma_glass
Temperature
JouleHeat
```

For the midterm demo, use stable effective values first. Do not immediately chase exact high-temperature viscosity if it destabilizes the flow. Keep all values in SI units.

---

## 7. EM Joule heat integration

### 7.1 EM solve modes

Add CLI or parameters:

```text
--ophelie-em-mode=off|fixed-joule|update-sigma
--target-power=50000
--ophelie-em-update-interval=<value>
--sigma-temperature-coupling=0|1
--ophelie-self-induction=off|one-way|picard-experimental
```

Recommended default for midterm:

```text
--ophelie-em-mode=fixed-joule
--target-power=50000
--ophelie-self-induction=off or one-way diagnostic
```

Do not make Picard self-induction the default midterm path.

### 7.2 First-stage coupling: fixed Joule heat

At initialization:

```text
1. Build French reduced glass geometry.
2. Assign coil source A_coil.
3. Solve complex edge-flux phi/E/J in the glass body.
4. Compute JouleHeatEdgeReconComplex.
5. Compute total Joule power:
   P = sum_i JouleHeat_i * Vol_i.
6. Scale source/current or JouleHeat to target power.
7. Store JouleHeat on glass particles.
```

Then run flow/thermal with fixed JouleHeat.

This is the fastest robust coupled demo.

### 7.3 Temperature source update

Add a dynamics class similar to:

```cpp
class ApplyJouleHeatingToTemperature : public LocalDynamics
{
  // For each glass particle:
  // Temperature[i] += dt * JouleHeat[i] / (rho_i * cp_i);
};
```

Use:

\[
\frac{dT}{dt}=\frac{Q_{Joule}}{\rho c_p}.
\]

If the code evolves internal energy instead of temperature, use the corresponding energy update, but for the midterm demo a direct temperature source is acceptable if documented.

Suggested time-step order:

```text
1. Fluid advection/acoustic update.
2. Apply Joule heat source to Temperature.
3. Thermal diffusion relaxation.
4. Rotor / wall contact updates as required.
5. Output.
```

Because SPHinXsys is Lagrangian, thermal advection is naturally carried by particles. The explicit thermal equation only needs source + diffusion.

---

## 8. Optional sigma(T) weak coupling

This is valuable but should be second-stage.

Weak coupling loop:

```text
At every EM update interval:
    sigma_raw_i = sigma(Temperature_i)
    sigma_i = (1 - alpha) * sigma_i + alpha * sigma_raw_i
    solve EM again
    update JouleHeat
    optionally rescale to target power
```

Recommended:

```text
alpha = 0.1--0.3
EM update interval >> acoustic timestep
```

Output every EM update:

```text
sigma_min
sigma_max
sigma_mean
P_joule_total
T_min
T_max
T_mean
EM_update_index
```

For the midterm, fixed JouleHeat is enough for the first demo. sigma(T) can be shown as an implemented or near-implemented extension.

---

## 9. Output requirements

### 9.1 VTK / VTP fields

At minimum output:

```text
Temperature
Velocity
Pressure
JouleHeat
|E|
|J|
sigma
```

If possible also output:

```text
EReal, EImag
JReal, JImag
ACoilReal, ACoilImag
AIndReal, AIndImag, if available
ATotalReal, ATotalImag, if available
```

### 9.2 CSV energy budget

Create a file such as:

```text
french_glass_em_thermal_budget.csv
```

Columns:

```text
time
P_joule_total
E_joule_cumulative
T_mean
T_min
T_max
thermal_energy
delta_thermal_energy
sigma_min
sigma_max
sigma_mean
EM_update_count
```

If wall heat loss is not included yet, state explicitly:

```text
wall_loss = not included in this phase-1 demo
```

Energy-budget consistency is more persuasive for the customer than field snapshots alone.

---

## 10. Suggested code construction steps for Cursor

### Step 1: create the new case skeleton

```text
Create test_3d_ophelie_french_glass_stirring.cpp.
Start from test_3d_stirring_injecting_RH200.cpp.
Remove all grain/two-phase/inflow logic.
Rename oil to glass.
Keep wall, rotor/stirrer, single-phase fluid stepping, and thermal diffusion.
```

Acceptance:

```text
Case compiles and runs as a pure single-phase stirring/thermal case without EM.
Outputs Temperature, Velocity, Pressure.
```

### Step 2: add a constant/fake Joule heat source first

Before integrating EM, add a controlled source:

```text
JouleHeat = constant or analytic spatial function.
```

Acceptance:

```text
Temperature rises according to Q/(rho cp).
Energy budget matches imposed Joule power.
No EM dependency yet.
```

This isolates thermal coupling bugs from EM bugs.

### Step 3: connect OPHELIE French EM JouleHeat

Use existing OPHELIE French/reduced solver functions to compute:

```text
EReal/EImag
JReal/JImag
JouleHeatEdgeReconComplex
```

Then map/use JouleHeat on the same glass particles.

Acceptance:

```text
P_joule_total computed from particles.
Target-power scaling works.
VTK shows nonuniform JouleHeat in glass.
Thermal run uses this JouleHeat.
```

### Step 4: fixed JouleHeat coupled run

Run:

```text
EM solve once at t=0.
JouleHeat fixed during flow/thermal simulation.
```

Acceptance:

```text
T_mean/T_max increase.
Spatial temperature distribution follows JouleHeat and stirring transport.
Energy CSV is reasonable.
```

### Step 5: sigma(T) update

Add weak update:

```text
update sigma(T) every N steps
re-solve EM
update JouleHeat
under-relax sigma
```

Acceptance:

```text
No numerical blow-up.
P_joule_total, sigma_min/max, T_min/max logged.
Can switch sigma(T) on/off.
```

---

## 11. What not to do now

Do not spend midterm-preparation time on:

```text
1. TEAM7 hole boundary patch.
2. TEAM7 strict validation_passed.
3. no-flux ghost-edge variants.
4. Cij empirical scaling.
5. J empirical scaling.
6. a_sign=-1 phase matching.
7. TEAM7 50 kW scaling.
8. Full metallic segmented crucible EM.
9. Water-cooled wall EM.
10. Large Picard relaxation sweeps.
```

TEAM7 remains a side benchmark only.

---

## 12. Midterm presentation plan

Recommended figures and results:

### 12.1 EM field plots

```text
|E|
|J|
JouleHeat
```

Show that the OPHELIE EM solver creates a spatial Joule heat distribution in the French glass body.

### 12.2 Thermal evolution plots

```text
Temperature at t0
Temperature at t1
Temperature at t2
```

Show heat accumulation and redistribution.

### 12.3 Stirring flow plots

```text
Velocity magnitude
Velocity vectors/glyphs
Temperature + velocity overlay
```

Show that stirring transports heat.

### 12.4 Energy budget curves

```text
P_joule_total vs time
T_mean/T_max vs time
Thermal energy vs time
```

These curves should be used to explain numerical consistency.

### 12.5 Scope slide

Explicitly state:

```text
Included in current phase:
  molten-glass EM, Joule heat, stirring flow, thermal transport.

Not included yet:
  metallic segmented crucible eddy currents, water-cooling wall EM, TEAM7 strict validation.
```

This prevents over-commitment.

---

## 13. Final recommendation

The correct near-term decision is:

```text
Stop treating TEAM7 as the main blocker.
Use TEAM7 only as a side benchmark.
Build the French reduced glass EM--thermal--stirring coupled case immediately.
```

The new case should be built from the uploaded stirring case, but simplified to single-phase glass and extended with OPHELIE Joule heat.

The first successful deliverable should be:

```text
fixed JouleHeat from one EM solve
+ single-phase stirring flow
+ thermal diffusion
+ energy budget
+ VTK visualization
```

The second deliverable should add:

```text
sigma(T) weak coupling
+ periodic EM update
```

This path is the best match for the customer midterm requirement and for the August delivery risk.
