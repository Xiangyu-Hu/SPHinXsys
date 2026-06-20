# OPHELIE Complex Edge-Flux Validation Roadmap, TEAM7 Standard Case, and Cursor Follow-On Development Plan

Date: 2026-06-09  
Purpose: Guide Cursor / Codex to continue the SPHinXsys-SYCL OPHELIE complex edge-flux solver.  
Current mainline: **complex edge-flux EM solver → self-induction Picard → 50 kW reduced French scale → TEAM7 standard eddy-current validation → one-way thermal demonstration**.

---

## 0. This Round Discussion Conclusion Summary

This discussion mainly answers two questions:

1. What does “Picard @ 50 kW literature calibration definition / scope” mean?
2. Since strict quantitative comparison with French OPHELIE literature is not possible now, should we validate with TEAM7 or another standard case?

Conclusions:

- **Picard @ 50 kW literature** refers to the total Joule power scale anchor from the French Jacoutot / OPHELIE cold-crucible reduced operating point, not pointwise reproduction of French literature temperature fields or thermocouple data.
- The current reduced French case lacks real cold-crucible structure, segmented walls, skull layer, natural convection, \(\sigma(T)\), cooling-water exchange, and thermo-flow coupling from full French literature; therefore we cannot claim strict quantitative validation against French literature.
- We should still run **Picard @ 50 kW** because it checks whether complex edge-flux + \(A_{ind}\) self-induction remains stable, calibratable, and convergent at realistic cold-crucible power scale.
- Standard EM validation should introduce **TEAM7**. TEAM7 validates coil/plate eddy currents, electromagnetic induction, \(B\) real/imag phase, and self-induction feedback—not cold-crucible thermal processes—but fits the current complex edge-flux EM solver well.
- Besides TEAM7, prioritize a **conducting cylinder harmonic skin-depth benchmark** because it is closer to “cylindrical conductor induction heating” physics and offers analytic or semi-analytic radial distributions for error localization.

Recommended validation hierarchy:

```text
1. Internal consistency: complex edge-flux uniform field / sign / scaling / Picard @ 50 kW
2. Analytic or semi-analytic heating: conducting cylinder harmonic skin-depth benchmark
3. Standard eddy-current benchmark: TEAM7 Level 1/2/3
4. French reduced cold-crucible scale: 50 kW + Picard + one-way heat
5. Later real cold crucible: real geometry, segmented walls, σ(T), thermo-flow coupling
```

---

## 1. Accurate Meaning of “Picard @ 50 kW Literature Calibration Definition / Scope”

“50 kW literature calibration definition / scope” means:

```text
French OPHELIE cold-crucible literature reduced / literature-inspired operating point with total induction heating power about 50 kW.
```

It does not mean we can already compare pointwise with French literature:

```text
Not: pointwise temperature field comparison
Not: quantitative thermocouple curve reproduction
Not: natural convection velocity field reproduction
Not: full cold crucible / skull layer / cooling-water exchange reproduction
```

Its role is a **power-scale acceptance**:

```text
If we say the current case is French reduced literature mode,
at minimum total Joule power should be calibrated to 50 kW,
and complex edge-flux + A_ind Picard should remain convergent, stable, and physically reasonable at that power scale.
```

Current Picard results included:

```text
coil-only complex one-way:       P ≈ 7.28 W
one-way feedback after A_ind:    P ≈ 8.64 W
Picard self-consistent:          P ≈ 6.39 W
```

These show the complex Picard prototype runs, but they are not yet 50 kW literature scale. With fixed EM model material parameters and geometry, theoretically:

\[
P \propto I^2 .
\]

Therefore use:

\[
I_{scale} = \sqrt{\frac{50000}{P_{Picard}}}
\]

to estimate coil current scaling, rerun Picard @ 50 kW, and validate:

```text
P_picard_scaled ≈ 50000 W
A_ind/A_coil essentially unchanged with I
B_ind/B_coil essentially unchanged with I
J_real/J_imag scale linearly with I
Q scales with I²
Picard still converges
Q spatial distribution still reasonable
```

This is reduced French EM solver scale-consistency acceptance, not full French literature quantitative validation.

---

## 2. Why Strict Quantitative Comparison with French Literature Is Not Possible Now

Current code already has:

```text
complex edge-flux
A_ind one-way
A_ind Picard prototype
JouleHeat complex
one-way heat
cold wall diffusion prototype
```

But French literature has many physics and geometry factors absent in the current reduced case:

```text
Real cold-crucible geometry
Segmented metal crucible walls
Skin-depth surface current approximation or conductor skin layer
Free glass surface
Skull layer / frozen glass layer
Strong temperature-dependent σ(T)
Natural convection / Boussinesq buoyancy
Cooling water / wall convective boundary
Multi-physics iteration under-relaxation
```

Therefore we can only say:

```text
We are doing a French literature-inspired reduced EM/thermal case.
```

Not:

```text
We have already reproduced the full French literature OPHELIE-FLUENT coupled natural-convection case.
```

Hence TEAM7 and analytic benchmarks are needed for EM solver credibility.

---

## 3. Why TEAM7 Should Be Introduced

TEAM7 is a standard eddy-current benchmark: coil excitation on a conductive plate with hole, eddy currents and magnetic response. It is especially valuable because it tests:

```text
coil source absolute amplitude
nonuniform alternating magnetic field
eddy currents in conductor
self-induced A/B
complex phasor real/imag phase
probe-line B_real / B_imag distribution
frequency response, e.g. 50 Hz / 200 Hz
```

This maps to our complex edge-flux mainline:

```text
A_coil_real / A_coil_imag
→ phi_real / phi_imag
→ JEdgeReconReal / JEdgeReconImag
→ A_ind_real / A_ind_imag
→ B_total_real / B_total_imag
```

TEAM7 does not validate cold-crucible thermal processes but is the standard external validation for the current EM solver. It is better than the French reduced case for judging:

```text
complex edge-flux formulation correctness
A_ind phase correctness
Picard self-induction reasonableness
B_real/B_imag trustworthiness
```

---

## 4. Recommended TEAM7 Implementation Roadmap

Do not start with full TEAM7 production. Split into three levels.

### 4.1 TEAM7 Level 1: Coil-Only Source-Field Validation

Goal: validate coil Biot-Savart source only; no plate solve, no A_ind.

Input:

```text
TEAM7 geometry: coil + conductive plate + probe line
Materials: plate σ, μ, etc.
Frequency: benchmark common values, e.g. 50 Hz / 200 Hz
Coil current: benchmark ampere-turns
```

Compute:

```text
A_coil_real
B_coil_real
B probe line
```

Output:

```text
team7_coil_only_probe.csv
x, y, z, Bx_src_real, By_src_real, Bz_src_real, |B_src|
```

Acceptance:

```text
coil geometry, units, ampere-turns, probe line positions correct
B_src peak location reasonable
B_src amplitude order of magnitude reasonable
if reference data available, RMS / peak error acceptable
```

This excludes:

```text
coil geometry error
mm/m unit error
ampere-turns error
probe line coordinate error
Biot-Savart sign/direction error
```

### 4.2 TEAM7 Level 2: Plate Edge-Flux One-Way

Goal: validate edge-flux solve in conductive plate and one-way \(A_{ind}\).

Flow:

```text
1. A_active = A_coil
2. solve complex edge-flux on plate
3. reconstruct J_real/J_imag
4. compute A_ind/B_ind from J
5. output B_total = B_coil + B_ind on probe line
```

No Picard initially, or one-way only. Output:

```text
team7_edge_flux_one_way_probe.csv
team7_edge_flux_plate_joule.csv
team7_edge_flux_current_stats.csv
```

Key metrics:

```text
phi_eq_res_real/imag
edge_res_red_real/imag
q_antisym_real/imag
P_recon_complex
max_J_real/max_J_imag
A_ind/A_coil
B_ind/B_coil
B_total_real/B_total_imag probe RMS error
```

### 4.3 TEAM7 Level 3: Full Complex Edge-Flux + A_ind Picard

Goal: full self-induction feedback validation.

Flow:

```text
A_total^k = A_coil + A_ind^k
solve complex edge-flux
reconstruct J_real/J_imag
A_ind^{k+1} = K[J_real/J_imag]
under-relax
repeat until convergence
sample B_total probe line
```

Output:

```text
team7_picard_convergence.csv
team7_picard_probe_final.csv
team7_picard_plate_fields.vtp
```

Acceptance:

```text
J_rel_raw / J_rel_relaxed converged
P_rel converged
Aind_rel converged
B_real probe RMS error
B_imag probe RMS error
frequency response: 50 Hz / 200 Hz trends correct
```

After Level 3 pass, complex edge-flux + self-induction has credibility on standard eddy-current benchmark.

---

## 5. TEAM7 Geometry Preparation Recommendations

User can prepare TEAM7 geometry bodies. Geometry requirements for Cursor:

### 5.1 Required Geometry Objects

```text
coil / excitation winding
conducting plate
hole / slot geometry, if benchmark requires
probe line / observation points
optional air box for visualization only, not for solving
```

### 5.2 Recommended Geometry Files

Prepare two sets with priority:

```text
1. CAD/mesh geometry: for particle relaxation or reload
2. probe point CSV: fixed-position B_real/B_imag sampling
```

Recommended files:

```text
team7_plate_geometry.step or .stl
team7_coil_centerline.csv
team7_probe_line.csv
team7_material_and_frequency.md
```

Recommended fields for `team7_coil_centerline.csv`:

```text
segment_id, x0, y0, z0, x1, y1, z1, current, turns
```

Recommended fields for `team7_probe_line.csv`:

```text
probe_id, x, y, z
```

`team7_material_and_frequency.md` should record:

```text
units: m / mm must be explicit
plate conductivity σ
relative permeability μ_r
frequency f
coil current / ampere-turns
reference probe data source
```

### 5.3 Unit Risk

TEAM7’s biggest risk is unit confusion. Cursor must print after reading geometry:

```text
geometry_unit=m or mm
plate bounding box
coil bounding box
probe bounding box
min/max coordinate
estimated plate size
estimated coil size
```

If user supplies mm geometry, convert uniformly to SI:

```text
m, kg, s, A, V
```

---

## 6. Additional Validation Cases Beyond TEAM7

### 6.1 Conducting Cylinder Harmonic Skin-Depth Benchmark (Strongly Recommended)

Besides TEAM7, this is the most worthwhile case—closer to cylindrical glass induction heating in a cold crucible.

Goal: validate in cylindrical conductor under harmonic B/A drive:

```text
skin depth δ
Jθ(r) radial distribution
Q(r) radial distribution
frequency scaling
conductivity scaling
```

Typical skin depth:

\[
\delta = \sqrt{\frac{2}{\omega \mu \sigma}} .
\]

Note: analytic solution depends on excitation form (axial B, azimuthal current, external coil, etc.); do not apply formulas blindly. Start with simplest known analytic/semi-analytic setup.

Recommended test levels:

```text
Cylinder Level 1: prescribed analytic A_real/A_imag, no coil Biot
Cylinder Level 2: simple circular coil source
Cylinder Level 3: self-induction Picard
```

Output:

```text
r, J_real_theta, J_imag_theta, Q, Q_normalized
P_total
skin_depth_estimate
```

Acceptance:

```text
J/Q enhanced outside, decay inside—correct trend
higher frequency → smaller skin depth
higher σ → smaller skin depth
P ∝ I²
```

### 6.2 Uniform Complex Edge-Flux Manufactured Tests

Keep existing internal tests and extend to complex:

```text
A_real-only → imag chain active, real chain zero
A_imag-only → real chain active, imag chain zero
A_real + A_imag → total Joule heat = real contribution + imag contribution
I×2 → P×4
```

These are not external benchmarks but protect the complex edge-flux kernel from regressions.

### 6.3 Other TEAM / Eddy-Current Benchmarks

Later options:

```text
TEAM3: simpler coil/plate magnetic field benchmark
TEAM11: conducting hollow sphere eddy-current benchmark, analytic/reference value
TEAM15: eddy-current testing slot/crack benchmark, complex geometry for later
induction heating billet benchmark: magnetothermal coupling for later thermal-EM validation
```

Current priority remains:

```text
TEAM7 + conducting cylinder skin-depth
```

---

## 7. Current OPHELIE Solver Validation Roadmap Recommendations

### 7.1 Do Now: Picard @ 50 kW Reduced French Scale

First priority for current complex edge-flux EM.

Must support:

```text
--target-power=50000
--ophelie-current-form=edge-flux
--ophelie-edge-flux-complex=1
--self-induction-picard=1
```

If current Picard test lacks target-power calibration, Cursor should add support first.

Output:

```text
P_picard_raw
P_picard_scaled
ampere_turns_eff
A_ind/A_coil
B_ind/B_coil
J_rel_raw
J_rel_relaxed
P_rel
Aind_rel
phi_eq_res_real/imag
Qmax/Qmean
Q_outer/center
```

Recommended gates:

```text
Smoke:
  J_rel_relaxed < 0.05
  P finite
  phi_eq_res < 1e-2

Production candidate:
  J_rel_relaxed < 0.01
  P_rel < 0.01
  Aind_rel < 0.01
  phi_eq_res < 1e-2 initially, later tighten
  Q spatial soft gate passed
```

### 7.2 Picard Convergence CSV + VTP Standardization

Each Picard iteration output:

```text
iter
J_rel_raw_before_relax
J_rel_after_relax
Aind_rel
P_complex
P_rel
phi_eq_res_real
phi_eq_res_imag
max_J_real
max_J_imag
A_ind/A_coil
B_ind/B_coil
```

VTP:

```text
coil-only initial
final Picard
optional every N iterations
```

Fields:

```text
PhiReal / PhiImag
EReal / EImag
JReal / JImag
EEdgeReconReal / EEdgeReconImag
JEdgeReconReal / JEdgeReconImag
JouleHeatEdgeReconComplex
AIndReal / AIndImag
ATotalReal / ATotalImag
```

### 7.3 Reload dp Convergence

Prefer reload/relaxed particles; do not use lattice as production accuracy conclusion.

Recommendation:

```text
dp = 0.02 anchor
dp = 0.015 next
dp = 0.01 later if resource allows
```

Output:

```text
n_particles
runtime
P_scaled
phi_eq_res_real/imag
edge_res_red_real/imag
q_antisym
A_ind/A_coil
Picard iterations
Qmax/Qmean
Q_outer/center
```

### 7.4 Thermal One-Way Stays Qualitative Demo

Thermal one-way is already closed; do not do literature quantitative comparison. Next steps only:

```text
DeltaT energy closure
explicit diffusion stability log
cold wall / Robin wall boundary sensitivity
Q/T VTP visualization
```

Not yet:

```text
natural convection
σ(T)-EM fully coupled
thermal-fluid coupling
```

---

## 8. Cursor Immediate Code Changes and Tests to Add

### P0-1: Picard @ 50 kW Calibration Regression

Add or modify test:

```text
test_3d_ophelie_french_self_induction_picard
```

Must support:

```text
--target-power=50000
--ophelie-edge-flux-complex=1
--ophelie-current-form=edge-flux
--reload=1
```

If only fixed current scale is available, first use:

\[
I_{scale}=\sqrt{50000/P_{raw}}
\]

then rerun. Final log must state:

```text
picard_target_power=50000
picard_power_before_scale
picard_current_scale
picard_power_after_scale
```

### P0-2: Strengthen Picard Convergence Gate

Current `J_rel≈0.039` is too loose. Add:

```text
J_rel_raw_before_relax
J_rel_after_relax
P_rel
Aind_rel
```

Split gates:

```text
picard_smoke_passed
picard_production_candidate_passed
```

### P0-3: edge_res_red_real n/a Mechanism

Do not use:

```text
edge_res_red_real = -1
```

Change to:

```text
real_chain_active=0/1
edge_res_real_gate=skipped/pass/fail
edge_res_red_real=n/a or numeric
```

Same for imag chain.

### P1-1: Unify Complex Power Metrics

When `edge_flux_complex=1`, all Q statistics must come from:

```text
JouleHeatEdgeReconComplex
```

Not imag-only:

```text
JouleHeatEdgeReconImag
```

Check:

```text
computeHostEdgeFluxPowerMetrics
hostEdgeFluxReconPower
Q spatial metrics
final log summary
```

Ensure consistent definition / scope.

### P1-2: TEAM7 Geometry and Test Framework

Recommended directory:

```text
tests/extra_source_and_tests/test_3d_ophelie_team7_complex_edge_flux/
```

Or extend existing TEAM7 test, but recommend separating complex edge-flux validation to avoid confusion with old TEAM7 smoke.

Add files:

```text
test_3d_ophelie_team7_complex_edge_flux.cpp
team7_geometry_loader.h
team7_probe_sampler.h
team7_validation_metrics.h
```

CLI:

```text
--team7-geometry-dir=...
--team7-probe-csv=...
--team7-frequency=50|200
--team7-level=coil-only|one-way|picard
--ophelie-edge-flux-complex=1
--ophelie-current-form=edge-flux
```

Output:

```text
team7_probe_coil_only.csv
team7_probe_one_way.csv
team7_probe_picard.csv
team7_picard_convergence.csv
```

Reference data may live under `docs/TEAM7-reference/`.

### P1-3: Conducting Cylinder Skin-Depth Benchmark

Add test:

```text
test_3d_ophelie_conducting_cylinder_skin_depth
```

CLI:

```text
--radius
--height
--sigma
--frequency
--mu-r
--source=analytic-a|coil
--ophelie-edge-flux-complex=1
```

Output:

```text
r_bin, J_real_theta, J_imag_theta, Q_mean, Q_max
skin_depth
P_total
```

### P2: Thermal One-Way Enhancement

Add log:

```text
thermal_diffusivity = k/(rho*cp)
dt_stable_estimate
dt/dt_stable
E_joule
E_thermal
E_lost_boundary
```

Later add Robin wall:

\[
q_{wall}=h(T-T_c).
\]

---

## 9. Code Structure Review and Recommendations

### 9.1 Keep Edge-Flux Production Mainline

Default development direction:

```text
edge-flux complex = production candidate
div-grad = fallback/reference
compatible-div-grad / Neumann / grad correction = archived diagnostic
```

Do not treat div-grad `phi_eq_res≈0.49` as mainline issue anymore.

### 9.2 TEAM7 Must Not Reuse Old Phi-Only Edge-Flux Route

New TEAM7 validation must use:

```text
--ophelie-current-form=edge-flux
--ophelie-edge-flux-complex=1
```

Do not use old:

```text
--phi-projection-operator=edge-flux
```

Unless explicitly as historical diagnostic.

### 9.3 Active Chain and Complex Field Naming

Continue:

```text
phi_real / phi_imag
EEdgeReconReal / EEdgeReconImag
JEdgeReconReal / JEdgeReconImag
JouleHeatEdgeReconComplex
AIndReal / AIndImag
ATotalReal / ATotalImag
```

Do not revert to unsuffixed fields, e.g.:

```text
EEdgeRecon
JEdgeRecon
```

If aliases kept, mark deprecated.

### 9.4 Pair Conductance Symmetry

Later add helper:

```cpp
computeOphelieEdgeFluxSymmetricPairConductance(...)
```

For consistent residual, q antisymmetry, and edge reconstruction pair-weight interpretation. Current q antisymmetry is small (case OK), but complex TEAM7 geometry / reload / multi-resolution needs more stability.

### 9.5 Preserve Thermal Delta Field

Thermal module continues:

```text
OphelieThermalDeltaT
Temperature = T0 + DeltaT
```

Do not accumulate small temperature rise directly on high benchmark baseline temperature.

---

## 10. Direct Execution Summary for Cursor

Copy-paste for Cursor:

```text
Current conclusion: French reduced 50 kW is power-scale acceptance, not French literature pointwise validation.
Next EM standard validation: TEAM7, plus conducting cylinder skin-depth benchmark.

P0:
1. Run complex Picard @ 50 kW regression.
   Support --target-power=50000 + --ophelie-edge-flux-complex=1 + --ophelie-current-form=edge-flux.
   Output P_raw/P_scaled/current_scale/A_ind/A_coil/B_ind/B_coil/J_rel/P_rel/Aind_rel/Q statistics.

2. Strengthen Picard convergence gate.
   Add J_rel_raw_before_relax, J_rel_after_relax, P_rel, Aind_rel.
   Split smoke gate and production candidate gate.

3. Fix edge_res_red_real=-1.
   Change to real_chain_active=0/1 and gate=skipped/pass/fail.

P1:
4. Unify complex power metrics.
   In complex mode all Qmax/Qmean/Q spatial metrics use JouleHeatEdgeReconComplex, no longer imag-only.

5. Prepare TEAM7 Level 1/2/3.
   Level 1 coil-only B probe.
   Level 2 complex edge-flux one-way A_ind.
   Level 3 full complex Picard B_real/B_imag probe.
   Add team7 geometry/probe loader and probe CSV output.

6. Add conducting cylinder harmonic skin-depth benchmark.
   Output Jθ(r), Q(r), P_total, skin_depth, frequency scaling.

P2:
7. Thermal one-way: qualitative demo and energy accounting only.
   Add explicit diffusion stability log.
   Later consider Robin wall; no natural convection or σ(T) fully coupled yet.

Forbidden:
- Do not revisit div-grad residual floor tuning.
- Do not do P_graph κ calibration.
- Do not use old phi-only edge-flux for TEAM7 production.
- No natural convection, full thermo-flow coupling, or σ(T)-EM feedback yet.
```

---

## 11. Final Roadmap

Recommended fixed route:

```text
A. complex edge-flux internal consistency
   - sign tests
   - uniform field tests
   - I² scaling
   - q antisymmetry

B. French reduced power-scale EM
   - Picard @ 50 kW
   - convergence CSV/VTP
   - reload dp sensitivity

C. Standard EM validation
   - TEAM7 Level 1 coil-only
   - TEAM7 Level 2 edge-flux one-way
   - TEAM7 Level 3 complex Picard
   - conducting cylinder skin-depth benchmark

D. one-way thermal demonstration
   - Joule → DeltaT energy closure
   - diffusion + cold wall / Robin wall
   - qualitative VTP only

E. future coupled physics
   - σ(T) diagnostic
   - thermal-EM under-relaxation
   - natural convection / full thermal-fluid coupling
```

This upgrades the solver from “French cold-crucible inspired demo” to:

```text
complex edge-flux OPHELIE-like EM solver validated by analytic benchmark, TEAM benchmark, and reduced cold-crucible scale.
```
