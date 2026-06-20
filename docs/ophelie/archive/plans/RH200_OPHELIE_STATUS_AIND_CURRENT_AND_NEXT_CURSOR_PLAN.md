# RH200 / OPHELIE-inspired EM Joule-Heat Solver Status, `A_ind` Self-Induction Feedback, Current Field Display, and Next-Step Tasks

> Target audience: Cursor / follow-on developers  
> Current goal: first confirm whether RH200 midterm case EM-thermal chain is self-consistent and reportable; do not use "reach 1500 K" as current acceptance target yet.

---

## 0. Current conclusion summary

The current RH200 case should not be interpreted as strict reproduction of French cold-crucible literature conditions, but understood as:

```text
RH200 stirring geometry
+ OPHELIE/French-inspired external induction coil source
+ complex phasor EM solve
+ Joule heat target-power scaling
+ Eulerian fixed-grid Joule heat source
+ SPH stirring/heat transfer coupling
```

From current `rh200_energy_budget.csv` 10 s data:

```text
P_joule_W ≈ 49.5–49.7 kW
T_mean: 373.15 K -> 373.998 K, ΔT_mean ≈ 0.85 K
T_max: ~375.02 K, local ΔT_max ≈ 1.9 K
E_joule_integrated ≈ 4.96e5 J
thermal_energy increase ≈ 4.81e5 J
energy closure ≈ 97%
```

Consistent with energy conservation estimate:

\[
\Delta \bar T \approx \frac{P t}{M c_p}.
\]

Therefore, **~1 K mean temperature rise in 10 s is not unreasonable**. More important work now is not forcing 10 s to 1500 K, but auditing:

```text
coil excitation
  -> A_coil
  -> phi
  -> E
  -> J
  -> Q_raw
  -> target-power scaling
  -> Q_grid(x)
  -> sampled JouleHeat_i(t)
  -> temperature rise
```

---

## 1. Relationship between current RH200 results and French literature

### 1.1 Typical French literature conditions

Core OPHELIE/French cold-crucible literature information:

```text
cold crucible diameter: ~650 mm
external high-frequency generator: ~400 kW, 300 kHz
glass bath temperature level: ~1500 K
total Joule power in glass body: ~50 kW
```

Note:

```text
400 kW is high-frequency generator/supply-side power;
50 kW is Joule power absorbed inside glass body;
they are not the same quantity.
```

Literature 1500 K is closer to high-temperature operating state or coupled-calculation condition temperature, not "from 373 K cold start, reach 1500 K in seconds with 50 kW heating".

### 1.2 Current RH200 geometry volume

From current STL estimate:

| Geometry | Approx. bbox size | Volume |
|---|---:|---:|
| `oil.stl` / initial glass domain | `0.748 x 0.748 x 0.526 m` | `~0.188 m^3` |
| `full_oil.stl` / full glass usable domain | `0.748 x 0.748 x 0.715 m` | `~0.266 m^3` |

If current `rho=2500 kg/m^3`, initial glass mass:

\[
M \approx 2500 \times 0.188 \approx 470~\mathrm{kg}.
\]

If current `cp=1200 J/(kg K)`:

\[
M c_p \approx 5.64\times 10^5~\mathrm{J/K}.
\]

Therefore 50 kW heating 10 s mean temperature rise:

\[
\Delta \bar T
\approx
\frac{50,000\times 10}{5.64\times 10^5}
\approx 0.89~\mathrm{K}.
\]

Very close to simulated `ΔT_mean≈0.85 K`.

### 1.3 Current parameters vs French literature

| Parameter | French literature typical | Current RH200 code/run | Assessment |
|---|---:|---:|---|
| frequency | `300 kHz` | `300 kHz` | match |
| absorbed Joule power in glass | `~50 kW` | `--target-power=50000`, sampled `~49.5–49.7 kW` | same order |
| electrical conductivity `sigma` | `16 S/m @ ~1473 K` | `16 S/m` | value match, but if `T0=373 K` physical definition / scope inconsistent |
| density `rho` | `2750 kg/m^3` | `2500 kg/m^3` | close, current ~9% lower |
| heat capacity `cp` | `1150 J/(kg K)` | `1200 J/(kg K)` | close, current ~4% higher |
| thermal conductivity `k` | `4 W/(m K)` | `0.8 W/(m K)` | current lower, affects temperature distribution smoothing |
| dynamic viscosity `mu` | `4 Pa s @ ~1473 K` | `0.05 Pa s` | very different, current more like numerical demo fluid |
| initial temperature | high-temperature glass bath / ~1500 K condition | `373.15 K` | definition / scope completely different |
| heat loss boundary | water-cooled wall, bottom, free-surface loss | current midterm mainly no external heat loss | inconsistent |

---

## 2. Why current RH200 midterm uses one-way / reduced EM, not immediate `A_ind` self-induction feedback?

### 2.1 What does `A_ind` self-induction feedback mean?

Full EM phasor chain should be:

\[
A_{total}=A_{coil}+A_{ind},
\]

where:

```text
A_coil: source vector potential from external coil prescribed current;
A_ind : vector potential from glass interior induced current J feeding back.
```

Full self-consistent iteration:

```text
A_coil
  -> solve phi
  -> compute E
  -> compute J
  -> compute A_ind = K[J]
  -> A_total = A_coil + A_ind
  -> solve phi again
  -> repeat until E/J/Q convergence
```

Needs Picard / Krylov / nonlinear fixed-point loop:

\[
J^{(m)} \rightarrow A_{ind}^{(m)} \rightarrow A_{total}^{(m)} \rightarrow E^{(m+1)} \rightarrow J^{(m+1)}.
\]

### 2.2 Why midterm uses one-way first?

Current one-way / reduced version:

\[
A_{total}\approx A_{coil},
\]

i.e.:

```text
coil field drives glass;
glass interior induced current no longer feeds back to total magnetic field.
```

Not because `A_ind` can never be implemented, but midterm goal priority differs:

1. **Midterm must first prove EM Joule heat → thermal → stirring coupling route works.**  
   Current 50 kW heat budget already ~97% closed; heat-source coupling chain usable.

2. **`A_ind = K[J]` introduces new integral operator and iteration convergence issues.**  
   Must implement complex current-to-potential pairwise/integral operator and self-consistent iteration with current `phi/E/J/Q` operators.

3. **TEAM7 high-conductivity metal benchmark not fully closed.**  
   Strict accuracy for high-conductivity metal, boundary no-flux, self-induction feedback not fully validated. RH200 glass is low conductivity, target-power scaled engineering demo; different priority.

4. **Current uses target-power scaling.**  
   In target-power mode, total absorbed power calibrated to 50 kW. Then `A_ind` mainly changes heat-source spatial distribution and phase, not total input power. Midterm validating energy chain more important first.

5. **Higher computational cost and implementation risk.**  
   Self-induction feedback needs multiple EM solves, possibly rebuilding `A_ind` each time; much costlier than one-way.

### 2.3 How to implement `A_ind` later?

Can be Stage-2 / Stage-3:

```text
Stage A: one-way A_coil -> phi/E/J/Q     [current midterm]
Stage B: frozen one-step A_ind = K[J] diagnostic
Stage C: Picard feedback loop A_total = A_coil + A_ind
Stage D: sigma(T) + periodic EM update
Stage E: cold-crucible metal segment / wall current coupling
```

Recommend diagnostic version first:

```text
1. compute J from current one-way result;
2. compute once A_ind = K[J];
3. output |A_ind|/|A_coil|, phase, spatial distribution;
4. judge whether self-induction feedback important.
```

If:

\[
\frac{\|A_{ind}\|}{\|A_{coil}\|}\ll 1,
\]

midterm one-way approximation easier to justify; if ratio not small, enter self-consistent feedback development.

---

## 3. Can we display glass interior current?

Yes, and recommended.

Current solver already has:

```text
EReal, EImag
JReal, JImag
JouleHeat
Sigma
```

Glass interior current density is complex phasor:

\[
J = J_r + iJ_i.
\]

Can output/display:

```text
JRealMagnitude
JImagMagnitude
JAmplitude = sqrt(|JReal|^2 + |JImag|^2)
JRMS = JAmplitude / sqrt(2)
JouleHeat = 0.5 * (|JReal|^2 + |JImag|^2) / sigma
```

Also vector output:

```text
JReal_x, JReal_y, JReal_z
JImag_x, JImag_y, JImag_z
```

ParaView display recommendations:

```text
1. JImagMagnitude or JAmplitude contour;
2. Glyph for JImag vector or phase-dominant vector;
3. overlay JouleHeat contour to show current concentration matches heat concentration;
4. if fixed Eulerian grid, also output grid-sampled JouleHeat;
5. threshold-filter current vectors to avoid low-value arrow pollution.
```

### 3.1 Does French literature display current?

In French OPHELIE literature, current density `j` is core EM unknown; Joule heat formula directly uses:

\[
Q_{th}=\frac{|j|^2}{2\sigma}.
\]

Therefore literature certainly computes glass interior current density field. Public papers emphasize temperature, velocity, solidification layer, Joule power density etc.; not every paper shows `j` vector arrows as main figure. From model perspective, displaying current density fully reasonable and better than JouleHeat alone to show what EM solver computes.

Midterm report recommend display:

```text
A_coil magnitude
PhiImag
EImagMagnitude
JImagMagnitude or JAmplitude
JouleHeat
Temperature
Velocity
```

---

## 4. What is currently reliable? What cannot be overstated?

### 4.1 Currently more reliable

| Content | Current judgment |
|---|---|
| 50 kW target-power can be applied | reliable, sampled power ~49.5–49.7 kW |
| JouleHeat -> temperature energy chain | reliable, ~97% energy closure in 10 s |
| 10 s temperature rise 0.85 K | reliable, consistent with `P t/(M cp)` |
| phasor phase relations | basically reliable, real coil source `AReal -> EImag/JImag/Q` reasonable |
| Eulerian grid fixed spatial heat source | physically correct direction, must replace particle-sticking mode |
| RH200 midterm coupling demo | already reportable stage |

### 4.2 Cannot yet overstate

| Content | Current limitation |
|---|---|
| French literature quantitative reproduction | not yet; geometry, coil, cooling, initial temperature, material properties not fully aligned |
| raw power for real coil current | still need output `P_raw` and `power_scale_factor` |
| `A_ind` self-induction feedback | self-consistent feedback not enabled |
| `sigma(T)`, `mu(T)` high-temperature glass feedback | not yet formal chain |
| cold-crucible water-cooled wall/free-surface heat loss | current midterm basically none |
| TEAM7 high-conductivity metal quantitative benchmark | not fully closed, cannot claim universal EM benchmark passed |

---

## 5. Next-step Cursor task manifest

### Task 1: complete EM raw/scaled audit

Add or improve:

```text
RH200_EXCITATION_TO_JOULE_AUDIT.csv
```

After each EM solve output:

```text
frequency
omega
coil_turns
coil_radius
coil_z_min
coil_z_max
coil_current_scale
sigma_min
sigma_max
sigma_mean
ACoilReal_rms
ACoilReal_max
ACoilImag_rms
ACoilImag_max
PhiReal_rms
PhiImag_rms
EReal_rms
EImag_rms
JReal_rms
JImag_rms
Q_raw_min
Q_raw_max
Q_raw_mean
P_joule_raw_particle
P_joule_raw_grid_sample_initial
target_power
power_scale_factor = target_power / P_joule_raw
equivalent_current_scale = sqrt(power_scale_factor)
P_joule_scaled_particle
P_joule_scaled_grid_sample_initial
```

Goal: clarify whether current 50 kW is natural from raw EM or from target-power scaling calibration.

### Task 2: complete Joule formula residual

Check raw and scaled separately:

\[
Q_{raw,ref}=0.5\sigma(|E_r|^2+|E_i|^2)
\]

or equivalent:

\[
Q_{raw,ref}=\frac{|J_r|^2+|J_i|^2}{2\sigma}.
\]

Output:

```text
joule_formula_rel_l2_raw
joule_formula_rel_max_raw
joule_formula_rel_l2_scaled
joule_formula_rel_max_scaled
```

### Task 3: complete heating-rate audit

Add or improve:

```text
RH200_HEATING_RATE_AUDIT.csv
```

Each output:

```text
time
glass_volume
glass_mass
cp
Mcp
target_power
P_grid_sample_current
integrated_joule_energy
thermal_energy
T_mean
T_min
T_max
dT_mean_measured
dT_mean_expected = integrated_joule_energy / Mcp
heating_rate_expected = P_grid_sample_current / Mcp
time_to_1500K_no_loss = Mcp * (1500 - T_mean) / P_grid_sample_current
out_of_grid_particle_count
```

### Task 4: output glass interior current field

VTP must output:

```text
JReal
JImag
JRealMagnitude
JImagMagnitude
JAmplitude
JRMS
JouleHeat
Sigma
```

If VTP scalar only, at least:

```text
JRealMagnitude
JImagMagnitude
JAmplitude
JRMS
```

If vectors possible, output `JReal` and `JImag` for ParaView glyph / streamlines.

### Task 5: align parameters with French literature, short-run check first

Do not pursue 1500 K first; build "French-like parameter preset":

```text
frequency = 300 kHz
sigma = 16 S/m
rho = 2750 kg/m^3
cp = 1150 J/(kg K)
k = 4 W/(m K)
mu = 4 Pa s  [physical high-temperature glass version]
T0 = 1473 K or 1500 K  [high-temperature operating version]
target_power = 50 kW
```

Also keep current numerically stable version:

```text
mu = 0.05 Pa s
T0 = 373.15 K
```

Recommend two presets:

```text
--preset=rh200-demo-current
--preset=rh200-french-like
```

Note: `mu=4 Pa s` may significantly change flow and time step; short-run stability first, not long run directly.

### Task 6: 10–30 s short-run check

Short runs for combinations:

```text
A. current parameters + 50 kW
B. french-like rho/cp/k/sigma + current mu + 50 kW
C. french-like full parameters including mu=4 + 50 kW
```

Acceptance:

```text
P_grid_sample_current ≈ target_power
energy closure > 95%
dT_mean_measured ≈ dT_expected
no NaN/Inf
out_of_grid_particle_count acceptable
```

### Task 7: longer runs after short-run pass

Without pursuing 1500 K, recommend first:

```text
60 s
120 s
300 s
```

Priority goals:

```text
1. temperature field transported by stirring trend;
2. heat source fixed in space not following particles;
3. T_mean/T_max/T_min vs time;
4. Joule energy and thermal energy long-time closure;
5. current stirring effect on hotspot diffusion and mixing.
```

### Task 8: optional target-power sweep

For client magnitude explanation:

```text
50 kW
200 kW
500 kW
1 MW
```

Each group 10–30 s sufficient, output:

```text
heating_rate_expected
heating_rate_measured
T_mean rise
T_max rise
```

Expected:

\[
\frac{d\bar T}{dt}\propto P.
\]

---

## 6. Recommended staged client statement

```text
RH200 case has completed OPHELIE-inspired EM Joule heat and SPH stirring/heat transfer coupling prototype.
Current uses external equivalent induction coil source, computes glass interior E/J/JouleHeat, normalizes total absorbed power to 50 kW.
~1 K mean and ~2 K peak temperature rise in 10 s consistent with P·t/(M cp) energy estimate; heat application and temperature response chain basically closed.

Model is not yet strict quantitative reproduction of French cold-crucible literature.
Next will complete raw EM power, power scaling factor, current density field, Joule formula residual and heating-rate audit tables,
and establish French-like parameter preset, then longer RH200 stirring heating runs.
```

---

## 7. Final judgment

Current OPHELIE-inspired RH200 solver:

```text
as midterm EM-Joule-thermal-stirring coupling prototype: credible, can continue;
as strict French cold-crucible 1500 K condition reproduction: not enough yet;
as complete benchmark-grade EM solver: still needs TEAM7 / A_ind / boundary / sigma(T) follow-on validation.
```

Next steps: do not rush to explain 1500 K or tune material properties to fit results. Cursor should first complete:

```text
raw/scaled power audit
Joule formula residual
heating-rate audit
current-density VTP output
French-like preset
short-run validation
long-run RH200 demo
```
