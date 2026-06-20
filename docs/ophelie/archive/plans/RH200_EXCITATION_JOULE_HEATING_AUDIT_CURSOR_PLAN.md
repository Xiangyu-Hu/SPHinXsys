# RH200 OPHELIE-style EM-Joule-Heating Self-Check and Cursor Task Plan

## 0. Purpose

The immediate goal is not to claim that the current RH200 case is a strict reproduction of the French cold-crucible paper. The goal is to make sure the current midterm-delivery calculation is not physically or numerically unreasonable.

The current chain to verify is:

\[
\text{coil excitation}
\rightarrow A_{\rm coil}
\rightarrow \phi
\rightarrow E
\rightarrow J
\rightarrow Q_{\rm Joule}
\rightarrow \Delta T .
\]

The main questions Cursor must answer are:

1. How much raw Joule power does the present coil excitation produce before target-power scaling?
2. How large is the post-scaling factor if `--target-power=50000` is used?
3. Does the stored Joule heat satisfy the formula
   \[
   Q_{\rm Joule}=\frac12\sigma\left(|E_r|^2+|E_i|^2\right)?
   \]
4. Does the measured temperature rise agree with the energy balance
   \[
   \Delta \bar T \approx \frac{\int_0^t P_{\rm Joule}(\tau)d\tau}{M_{\rm glass}c_p}?
   \]
5. Are current temperature values, units, material properties, total power, and glass mass consistent with the customer's statement that experimental temperatures can reach about 1500 K?

---

## 1. Current case status and terminology

### 1.1 Current case is not a strict French literature reproduction

The current RH200 case should be described as:

> RH200 stirred-vessel geometry + OPHELIE/French-inspired complex edge-flux electromagnetic model + equivalent multiloop coil + target-power scaling.

It should **not** be described as:

> strict reproduction of the French cold-crucible literature geometry and excitation.

Reasons:

- RH200 uses `rotor.stl`, `oil.stl`, and `full_oil.stl`; the glass bath geometry differs from the French cylindrical cold-crucible geometry.
- RH200 includes a mechanical stirrer; the natural-convection French case does not necessarily have the same stirrer geometry.
- The current inductor is an equivalent analytical multiloop coil generated from the RH200 fluid/full-oil geometry, not the exact literature coil geometry.
- The current version neglects the water-cooled wall thermal losses and does not explicitly solve the metal crucible/segment currents.
- The current midterm goal is a coupled prototype: EM Joule heating -> thermal response -> stirring transport.

### 1.2 Meaning of `--target-power=50000`

`--target-power=50000` means:

\[
P_{\rm target}=50,000\, {\rm W}=50\,{\rm kW}.
\]

The solver first computes a raw Joule heat field:

\[
Q_{\rm raw}(x)=\frac12\sigma\left(|E_r|^2+|E_i|^2\right),
\]

then integrates it over the glass domain:

\[
P_{\rm raw}=\int_\Omega Q_{\rm raw}(x)dV \approx \sum_i Q_{{\rm raw},i}V_i.
\]

If target-power scaling is active, compute:

\[
s_P=\frac{P_{\rm target}}{P_{\rm raw}},
\]

and set:

\[
Q_{\rm scaled}(x)=s_P Q_{\rm raw}(x).
\]

This is a **power normalization / post-scaling**, not a physical boundary condition. It preserves the spatial pattern of Joule heat but changes the global intensity.

The equivalent current-amplitude scaling is approximately:

\[
s_I=\sqrt{s_P},
\]

because in the one-way linear regime:

\[
Q\propto E^2\propto A^2\propto I^2.
\]

Therefore, every run must report both raw and scaled quantities.

---

## 2. Literature loading and how it differs from the current RH200 setup

### 2.1 What the French cold-crucible literature says

The French cold-crucible modelling papers describe a direct-induction process where the inductor is supplied by a sinusoidal high-frequency current. For the cold crucible described in the natural-convection paper, the crucible diameter is about 650 mm, and the inductor is connected to a high-frequency generator delivering 400 kW at 300 kHz. The process directly induces currents in the molten glass, and the dissipated Joule power heats/melts the glass.

The thermohydrodynamic model uses the Joule power density from the electromagnetic solver as a thermal source:

\[
Q_{\rm th}=\frac{|j|^2}{2\sigma}.
\]

Since \(j=\sigma E\), this is equivalent to:

\[
Q_{\rm th}=\frac12\sigma |E|^2.
\]

The electromagnetic part is solved by OPHELIE using integral methods. In the cited modelling description:

- The unknowns include scalar electric potential \(V\) and current density \(j\).
- Current conservation is imposed as \(\nabla\cdot j=0\).
- The current density is related to the potential and time-varying vector potential.
- The magnetic vector potential is computed using the Biot-Savart integral.
- Conducting elements such as crucible, coil, and base use a surface treatment with an exponential skin-depth law.
- The molten glass uses a volume mesh because the skin depth is not extremely thin in the glass.
- The electromagnetic and thermal/hydrodynamic meshes can differ, and Joule power is interpolated/coupled between OPHELIE and the thermal/fluid solver.

The literature also mentions results where the glass bath temperature is around 1500 K and the total Joule power in the glass is about 50 kW. This means that the generator output, the actual coil current, and the absorbed Joule power in the glass are not necessarily the same quantity.

### 2.2 Difference from our current RH200 midterm case

Our current RH200 setup is closer to an engineering prototype:

- We generate an equivalent multiloop coil around the RH200 `full_oil`/glass domain.
- The coil is not a meshed physical conductor with its own surface current solution.
- We compute \(A_{\rm coil}\), solve the glass response, compute \(E/J/Q\), then optionally scale \(Q\) to a prescribed target power.
- We currently do not explicitly compute the metal crucible/segment currents or water-cooled wall heat losses.

Therefore, the current result is acceptable for a midterm coupled demonstration, but any claim about quantitative agreement with 1500 K experiments requires further information:

- exact input power or absorbed power,
- real coil geometry and current,
- initial temperature,
- glass volume/mass,
- \(\sigma(T)\), \(c_p(T)\), \(k(T)\), viscosity,
- wall/free-surface losses,
- measurement time and measurement location.

---

## 3. Can we change material properties to alter Joule heat and temperature rise?

Yes, but the effect depends on whether target-power scaling is active.

### 3.1 If target-power scaling is OFF

With no target scaling, physical parameters directly affect raw Joule heat:

| Parameter | Main effect |
|---|---|
| coil current amplitude \(I\) | \(P_{\rm raw}\propto I^2\) approximately |
| frequency \(f\) | in one-way regime, \(P_{\rm raw}\propto f^2\) approximately |
| electrical conductivity \(\sigma\) | affects both current and Joule heat; approximately increases raw power in the weak-coupling regime |
| coil radius/height/location | changes both total raw power and spatial distribution |
| glass volume/mass | changes total absorbed power if the EM field region changes; also changes heating rate |

### 3.2 If target-power scaling is ON

If `--target-power` is used, the final integrated heating power is forced to match the target:

\[
P_{\rm scaled}=P_{\rm target}.
\]

Then:

- Changing \(\sigma\), frequency, or coil geometry mostly changes the **shape** of \(Q(x)\) and the required scaling factor, not the final total power.
- Changing \(c_p\), density, or glass mass changes the mean heating rate:
  \[
  \frac{d\bar T}{dt}=\frac{P_{\rm target}}{M_{\rm glass}c_p}.
  \]
- Changing thermal conductivity \(k\) changes the spatial smoothing of temperature, not the global mean temperature increase if there is no heat loss.
- Changing viscosity affects stirring/mixing, not the direct energy input.

### 3.3 What not to do

Do not arbitrarily tune properties only to force the temperature to reach 1500 K in 10 seconds. Instead, run a transparent sensitivity study:

- target power sweep,
- \(\sigma\) sweep,
- frequency sweep,
- \(c_p\)/mass audit,
- optional \(\sigma(T)\) weak coupling.

---

## 4. Mandatory self-check outputs

Cursor should implement the following diagnostics before using the result in the client presentation.

### 4.1 `RH200_EXCITATION_TO_JOULE_AUDIT.csv`

Output once after each EM solve.

Required fields:

```text
case_name
joule_mode
em_mode
frequency
omega
coil_current_scale
coil_turns
coil_radius
coil_z_min
coil_z_max
coil_center_x
coil_center_y
coil_center_z
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
EReal_max
EImag_rms
EImag_max
JReal_rms
JReal_max
JImag_rms
JImag_max
JouleHeat_raw_min
JouleHeat_raw_max
JouleHeat_raw_mean
P_joule_raw_particle
P_joule_raw_grid_sample_initial
target_power
power_scale_factor
equivalent_current_scale
P_joule_scaled_particle
P_joule_scaled_grid_sample_initial
glass_volume
glass_mass
cp
expected_heating_rate_raw
expected_heating_rate_scaled
```

Definitions:

\[
P_{\rm raw}=\sum_i Q_{{\rm raw},i}V_i.
\]

\[
s_P=P_{\rm target}/P_{\rm raw}.
\]

\[
s_I=\sqrt{s_P}.
\]

\[
\left.\frac{d\bar T}{dt}\right|_{\rm expected}=\frac{P}{M_{\rm glass}c_p}.
\]

Acceptance:

- `P_joule_scaled_*` should match `target_power` if target-power scaling is active.
- `power_scale_factor` should be reported transparently. Large values are not automatically wrong, but they mean the raw coil excitation is far from the target absorbed power.
- `ACoilReal` should be nonzero for real-coil-source runs.
- `ACoilImag`, `EReal`, `JReal` may be near zero in the current one-way real-coil phasor setup.

---

### 4.2 Joule heat formula residual

Add a check:

\[
Q_{\rm ref,i}=\frac12\sigma_i\left(|E_{r,i}|^2+|E_{i,i}|^2\right).
\]

Output:

```text
joule_formula_raw_rel_l2
joule_formula_raw_rel_max
joule_formula_scaled_rel_l2
joule_formula_scaled_rel_max
```

Important:

- If E/J are raw but JouleHeat is scaled, the formula residual will not match unless the comparison uses consistent raw/scaled fields.
- Cursor must explicitly distinguish raw and scaled fields or document exactly which variables are raw and which are scaled.

Acceptance:

- Formula residual should be near machine precision if all compared fields are consistent.
- If target scaling is applied only to JouleHeat, then do not compare scaled JouleHeat against raw E/J without accounting for the scaling factor.

---

### 4.3 `RH200_HEATING_RATE_AUDIT.csv`

Output during the coupled heat-flow run.

Required fields:

```text
time
glass_volume
glass_mass
cp
target_power
P_grid_sample_current
integrated_joule_energy
thermal_energy
T_mean
T_min
T_max
T_mean_initial
dT_mean_measured
dT_mean_expected
heating_rate_expected
out_of_grid_particle_count
out_of_grid_particle_fraction
```

Definitions:

\[
\Delta \bar T_{\rm measured}=\bar T(t)-\bar T(0).
\]

\[
\Delta \bar T_{\rm expected}=\frac{\int_0^t P_{\rm Joule}(\tau)d\tau}{M_{\rm glass}c_p}.
\]

\[
\left.\frac{d\bar T}{dt}\right|_{\rm expected}=\frac{P_{\rm grid\_sample\_current}}{M_{\rm glass}c_p}.
\]

Acceptance:

- If no external heat loss is included, `dT_mean_measured` should be close to `dT_mean_expected`.
- If the two differ strongly, check units, density, volume, cp, duplicate/lost source terms, grid sampling, and temperature units.

---

## 5. Required parameter sweeps

### 5.1 Raw coil current scale sweep

Run with target-power scaling OFF first.

Suggested sweep:

```text
coil_current_scale = 1, 2, 5, 10
```

Expected approximate relation:

\[
P_{\rm raw}\propto I^2.
\]

So:

```text
I/I0 = 1  -> P/P0 ≈ 1
I/I0 = 2  -> P/P0 ≈ 4
I/I0 = 5  -> P/P0 ≈ 25
I/I0 = 10 -> P/P0 ≈ 100
```

If this fails, inspect source scaling, A-field assignment, E/J reconstruction, JouleHeat formula, and volume integration.

### 5.2 Frequency sweep

Suggested sweep:

```text
frequency = 150 kHz, 300 kHz, 600 kHz
```

In the current one-way regime, approximately:

\[
P_{\rm raw}\propto f^2.
\]

This is not a strict law once strong feedback or skin effects are included, but it should be a useful check for the current one-way prototype.

### 5.3 Electrical conductivity sweep

Suggested sweep:

```text
sigma = 4, 8, 16, 32 S/m
```

Expected: raw absorbed power should increase with \(\sigma\) in the present weakly coupled one-way regime.

### 5.4 Target power sweep for client-facing heating rates

Suggested sweep:

```text
target_power = 50 kW, 200 kW, 500 kW, 1 MW
```

For each run, short simulation time is enough, e.g. 10-30 s.

Expected:

\[
\frac{d\bar T}{dt}\propto P_{\rm target}
\]

if there is no external heat loss.

This sweep is useful to explain why 50 kW over 10 s may produce only a few K in a large glass volume, while higher input powers or longer heating times are required to reach 1500 K.

---

## 6. Temperature unit and initial condition checks

Cursor must explicitly report and document:

```text
Temperature_internal_unit = K or C
T_initial_input
T_initial_K
T_initial_C
T_target_customer_K = 1500
T_target_customer_C = 1226.85
```

Recommendation:

- Use Kelvin internally.
- Output a derived `TemperatureC` field only for visualization.

If the case starts from 375 °C:

\[
T_0=648.15\,K.
\]

Then reaching 1500 K requires:

\[
\Delta T=851.85\,K.
\]

The no-loss time estimate is:

\[
t_{1500K}=\frac{M_{\rm glass}c_p(1500-T_0)}{P_{\rm Joule}}.
\]

Output this estimate in the heating audit.

---

## 7. Clarification for current observed result: 10 s gives about 2 K rise

A 2 K mean rise over 10 s can be physically reasonable.

Example if \(P=50\,kW\):

\[
M c_p=\frac{P t}{\Delta T}=\frac{50000\times10}{2}=2.5\times10^5\,J/K.
\]

If \(c_p\approx1200\,J/(kg\,K)\), then:

\[
M\approx208\,kg.
\]

This is plausible for a large RH200 full-oil/glass domain. Therefore, do not treat the 2 K rise as an error until the energy audit proves inconsistency.

---

## 8. How to explain this to the client

Suggested wording for the report:

> The current RH200 calculation is a coupled prototype. The electromagnetic solver provides a spatial Joule heat distribution, which is then normalized to a prescribed absorbed power and coupled into the SPH thermal-flow solver. Since the RH200 geometry, equivalent coil, and thermal boundary conditions are not yet the same as the French cold-crucible experiment, the current short-time temperature rise should be interpreted through the energy balance \(\Delta T=\int Pdt/(Mc_p)\), rather than as a direct reproduction of the reported 1500 K operating temperature.

---

## 9. Immediate Cursor task list

### Task A - Output raw/scaled EM audit

Implement `RH200_EXCITATION_TO_JOULE_AUDIT.csv` with all fields in Section 4.1.

### Task B - Add Joule heat formula residual

Check raw and scaled consistency separately.

### Task C - Output heating-rate audit

Implement `RH200_HEATING_RATE_AUDIT.csv` and compare measured vs expected mean temperature rise.

### Task D - Add no-loss time-to-target estimate

Output:

```text
time_to_1500K_no_loss = glass_mass * cp * (1500 - T_initial_K) / P_grid_sample_current
```

### Task E - Add parameter sweeps

Add CLI/script support for:

```text
coil_current_scale sweep
afrequency sweep
sigma sweep
target_power sweep
```

Correct typo in script names/labels if using `afrequency`; use `frequency`.

### Task F - Preserve transparency

Every output must state whether the displayed JouleHeat is:

```text
raw EM JouleHeat
scaled JouleHeat
sampled Eulerian-grid JouleHeat
particle-tagged debug JouleHeat
```

### Task G - Keep current official mode

For official midterm results, use:

```text
--joule-mode=em-grid
```

not particle-tagged `em-fixed`.

---

## 10. What not to do right now

Do not:

- tune material properties only to force 1500 K;
- claim strict French geometry reproduction using RH200 STL;
- compare 10 s transient directly against a final/operating 1500 K experimental temperature;
- hide the target-power scaling factor;
- use particle-carried JouleHeat as official result;
- re-enter TEAM7 boundary debugging for this midterm deliverable.

---

## 11. Summary decision

The next step is to verify the quantitative chain:

\[
\text{excitation} \rightarrow P_{\rm raw} \rightarrow P_{\rm scaled} \rightarrow \Delta T.
\]

If this chain closes, the current calculation is not numerically absurd even if the 10 s temperature rise is small. The remaining difference from 1500 K is then a matter of input power, glass mass, initial temperature, geometry, and missing experimental boundary conditions.
