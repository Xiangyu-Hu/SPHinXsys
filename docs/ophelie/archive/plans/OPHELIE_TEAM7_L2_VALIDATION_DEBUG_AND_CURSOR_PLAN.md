# OPHELIE / TEAM7 L2 Edge-Flux Validation Issue Review and Next-Steps Execution Plan

**File purpose**  
This document consolidates this round's discussion conclusions on TEAM7 L2 complex edge-flux validation, mainly for Cursor follow-on code changes, re-runs, and validation definition / scope usage.  
Current goal is not to "tune" a CSV to pass with post-hoc scaling, but to establish a repeatable, physically clear TEAM7 standard validation workflow that supports later L3 Picard.

---

## 0. Current conclusion summary

This round's TEAM7 L2 progress is valuable, but **TEAM7 L2 has not truly passed validation yet**.

Currently confirmed:

- TEAM7 native/reload case can run;
- complex edge-flux one-way pipeline can run;
- coil source field phase0 aligns peak under `source_scale=0.754` peak-fit definition / scope;
- at `a_sign=+1`, `J_edge` and `σE` basically self-consistent;
- A_ind / B_ind decomposition, probe CSV, x-band RMS, depth profile, post-scale diagnostic framework already built;
- NaN/overflow and other basic stability issues largely excluded.

But still not closed:

- TEAM7 phase90 not physically passed;
- `j_post_scale=auto` is post-hoc scaling, diagnostic only, not acceptance;
- `a_sign=-1` can flip phase90 sign but breaks edge-drop and `σE` consistency, cannot be production default;
- `Bind/Bcoil≈3.9` amplitude-too-strong issue not yet explained;
- L2 one-way is not TEAM7 standard coupled eddy-current solution; must enter L3 complex edge-flux + A_ind Picard;
- current `passed=1` should be read as smoke passed, not TEAM7 validation passed.

Current correct status:

```text
TEAM7 L1/L2 framework: built
TEAM7 coil-source diagnostic: partial
TEAM7 L2 one-way edge-flux smoke: passed
TEAM7 phase90 validation: failed / not closed
TEAM7 full validation: not passed
```

---

## 1. This round key numbers and interpretation

From Cursor summary and uploaded logs, key numbers:

| Scenario | phase90 RMS | Peak sim/ref | Description |
|---|---:|---:|---|
| baseline `a_sign=+1` | 11.69 | `-17.87 / 1.415 mT` | solver internally more self-consistent, but TEAM7 phase90 sign/amplitude wrong |
| `a_sign=-1` | 9.81 | `+17.87 / 1.415 mT` | phase90 sign flipped, amplitude still too large, constitutive consistency worse |
| `a_sign=-1 + j_post_scale=auto` | 0.41 | `1.17 / 1.415 mT` | post-hoc scaling makes curve close, but not physical solution |
| coil phase0 | 0.32 | `7.81 / 7.81 mT` | aligned under `source_scale=0.754` peak-fit, does not mean physical source fully correct |
| `Bind/Bcoil` | ~3.91 | — | induced field relative to coil field too strong, current main bottleneck |
| `p_graph/p_recon` | ~646 | — | graph energy still not physical Joule power |

These results mean:

1. **phase90 deviation is not simple sign error.**  
   If only sign issue, `a_sign=-1` should be close to reference; actually only flips negative peak to positive, amplitude still ~order of magnitude too large.

2. **`j_post_scale≈1/3.91` success means induced response amplitude overall too strong.**  
   Diagnostic value, but cannot be used as validation correction.

3. **`a_sign=+1` is more like internally physically self-consistent sign.**  
   Because at `a_sign=+1`, `J_edge≈σE`, while `a_sign=-1` significantly increases edge-EMF mismatch.

---

## 2. Correct understanding of TEAM7 phase0 / phase90

TEAM7 reference data usually gives Bz at `ωt=0°` and `ωt=90°`. Distinguish:

- solver internal phasor real/imag components;
- instantaneous times in TEAM7 tables;
- phase mapping chosen when code outputs probe.

If convention:

\[
B(t)=\Re\{(B_r+iB_i)e^{i\omega t}\}
     = B_r\cos(\omega t)-B_i\sin(\omega t),
\]

then:

\[
B(0^\circ)=B_r,
\]

\[
B(90^\circ)=-B_i.
\]

If:

\[
B(t)=\Re\{(B_r+iB_i)e^{-i\omega t}\}
     = B_r\cos(\omega t)+B_i\sin(\omega t),
\]

then:

\[
B(90^\circ)=+B_i.
\]

Therefore, TEAM7 phase90 comparison should provide phase mapping at **probe comparison layer**:

```text
B_phase90_sim = +B_imag
```

or:

```text
B_phase90_sim = -B_imag
```

Not by modifying edge-flux solver internal `a_sign` to fit reference tables.

### Conclusion

**Do not use `a_sign` to fix TEAM7 phase90 sign.**  
`a_sign` is equation definition; phase90 sign belongs to probe post-processing convention.

---

## 3. Decision on imag chain `a_sign`

Should freeze:

```text
imag a_sign = +1
```

Reasons:

- at `a_sign=+1`, `J_edge` and `σE` basically consistent;
- `a_sign=-1` makes phase90 sign positive but `J_em/J_edge` and `e_edge_em_mismatch` clearly worse;
- `a_sign=-1` is more like putting probe phase convention error into solver equation.

### Cursor should modify

Keep underlying parameter:

```text
--ophelie-edge-flux-imag-a-sign
```

But must label:

```text
debug-only, not for production TEAM7 validation
```

Add formal probe-layer parameter:

```text
--team7-phase90-convention=minus-imag|plus-imag
```

or reuse existing but clearly named:

```text
--team7-probe-phase90-neg-ind=1
```

Recommend default:

```text
--team7-phase90-convention=minus-imag
```

Then re-run:

```bash
test_3d_ophelie_team7_complex_edge_flux \
  --team7-level=one-way \
  --ophelie-edge-flux-complex=1 \
  --ophelie-edge-flux-imag-a-sign=1 \
  --team7-phase90-convention=minus-imag \
  --team7-output-tag=asign_p1_phase90_minus
```

Goal: validate TEAM7 phase90 via probe phase convention only, without changing solver equation.

---

## 4. Current main bottleneck: `Bind/Bcoil≈3.9`

Current core bottleneck:

```text
Bind/Bcoil ≈ 3.91
Aind/Acoil ≈ 3.96
```

Induced current/magnetic field amplitude overall too strong.

Cannot formally fix with `j_post_scale`. `j_post_scale=auto` only means:

```text
if J and B_ind post-hoc scaled down ~3.9x, phase90 shape can approach reference.
```

Reveals amplitude scale issue, does not explain physical cause.

### Possible sources by priority

Check in order:

1. **L2 one-way lacks self-induction feedback / shielding**  
   TEAM7 is standard eddy-current coupled problem, not pure one-way post-processing.  
   one-way computes:

   ```text
   A_coil -> phi -> J -> A_ind/B_ind
   ```

   without re-solving:

   ```text
   A_total = A_coil + A_ind
   ```

   for `phi/J`, no Picard iteration. For high-conductivity aluminum plate, may greatly exaggerate induced response.

2. **coil source still not physically closed**  
   Current `source_scale=0.754` is peak-fit, not standard TEAM7 source validation. Source shape, peak location, full-line RMS still need recheck.

3. **TEAM7 phase90 convention not yet separated from solver sign**  
   Must fix `a_sign=+1`, handle `±B_imag` only at probe layer.

4. **edge-flux amplitude normalization / pair conductance / reconstruction scale on high-conductivity plate**  
   If still too large after L1 source and L3 Picard, revisit this layer.

5. **Biot/probe definition or units**  
   Check probe coordinates, mT/T units, plate/coil geometry position, relative coordinate system.

---

## 5. Why L2 one-way is not TEAM7 standard solution

Current L2 physical chain:

```text
A_coil -> edge-flux phi -> J -> Biot A_ind/B_ind -> probe
```

This is one-way induced-field diagnostic.

TEAM7 standard eddy-current problem essentially needs self-consistent induced-field feedback. At least:

```text
A_total^k = A_coil + A_ind^k
edge-flux solve phi/J using A_total^k
A_ind^{k+1} = K[J^k]
Picard iterate
```

So true TEAM7 validation should be:

```text
TEAM7 L3 = complex edge-flux + A_ind Picard
```

L2 one-way only for locating:

- whether coil source reasonable;
- whether edge-flux can generate finite J;
- whether Biot post-processing can output B_ind;
- whether phase convention / output / probe pipeline usable.

But L2 phase90 should not be final hard acceptance.

---

## 6. Coil source issue: `source_scale=0.754` cannot be formal validation

Currently using:

```text
--team7-coil-source-scale=0.754
```

aligns phase0 peak to:

```text
peak 7.81 / 7.81 mT
```

But only peak-fit. Does not prove coil source correct.

Note:

- this scale is not part of TEAM7 physical input;
- TEAM7 standard input should use fixed ampere-turns, e.g. 2742 AT;
- `source_scale=0.754` masks coil geometry, probe position or Biot source implementation issues;
- phase0 peak alignment cannot replace full-line source-only Bz validation;
- if peak location offset, geometry/coordinates/probe line may still have issues.

### Cursor should modify

1. Change default back to:

```text
--team7-coil-source-scale=1.0
```

2. If user explicitly sets non-1:

```text
--team7-coil-source-scale=0.754
```

must print warning:

```text
[team7][warning] coil_source_scale != 1.0: this is a diagnostic amplitude fit and cannot be used as strict TEAM7 validation.
```

3. For TEAM7 validation gate:

```text
source_scale != 1.0 -> validation_passed = 0
```

Can allow:

```text
smoke_passed = 1
diagnostic_only = 1
```

---

## 7. TEAM7 should not use 50 kW calibration

TEAM7 validation definition / scope should be:

```text
fixed ampere-turns, e.g. 2742 AT
fixed frequency 50 / 200 Hz
fixed aluminum conductivity
compare Bz / Jey / phase0 / phase90
```

French reduced cold-crucible 50 kW calibration only for French literature-inspired case. Should not be used for TEAM7.

### Cursor should explicitly forbid

In TEAM7 case:

```text
--target-power=50000
```

should not be validation entry.  
If parameter kept, only debug/sensitivity, not TEAM7 validation.

Recommend log:

```text
[team7] TEAM7 validation uses fixed coil excitation, not target-power calibration.
```

---

## 8. Skin depth and depth profile interpretation

50 Hz aluminum skin depth estimate:

\[
\delta=\sqrt{\frac{2}{\omega\mu_0\sigma}}
\]

With:

```text
σ = 3.526e7 S/m
f = 50 Hz
μ = μ0
```

Gives:

```text
δ ≈ 12 mm
```

TEAM7 plate thickness ~19 mm, so at 50 Hz should not be extremely surface-localized; current through substantial thickness is possible.

But:

- if J at all depth layers too similar;
- and all layers `Bind/Bcoil≈3.5–4.9` too high;

then not just skin effect, more like overall J response amplitude too strong, or missing self-induction feedback shielding.

At 200 Hz:

```text
δ ≈ 6 mm
```

Should see clearer depth decay. TEAM7 must add 200 Hz test later.

---

## 9. About `P_graph / P_recon`

In TEAM7:

```text
P_graph ≈ 26 kW
P_recon ≈ 40 W
P_graph/P_recon ≈ 646
```

This means:

```text
P_graph = Σ 1/4 C_ij edge_drop^2
```

still not physical Joule power.

Consistent with French case conclusion:

```text
P_graph_edge = graph/Laplace energy diagnostic
P_recon = physical Joule power
```

### Cursor should modify log language

Do not continue citing in TEAM7 output:

```text
French soft gate [0.5, 2]
```

Unify as:

```text
P_graph_edge_diagnostic
P_recon_physical
```

And label:

```text
P_graph is diagnostic-only and is not used for TEAM7 validation.
```

### About `P_vol≈80 W` vs `P_recon≈40 W`

Likely harmonic average power 1/2 factor:

\[
P_{\rm avg}=\frac12\int J\cdot E\,dV.
\]

Therefore:

- `∫J·E dV ≈ 80 W`
- `0.5∫σ|E|^2dV ≈ 40 W`

Not necessarily contradictory. Unify log naming.

---

## 10. Code structure review

### 10.1 What is done well

Current TEAM7 validation framework already fairly complete, especially:

```text
electromagnetic_ophelie_team7_validation.h
```

Already has:

- probe decomposition;
- phase90 variant comparison;
- x-band RMS;
- depth profile;
- J_em vs J_edge audit;
- post-scale diagnostic;
- validation history.

`test_3d_ophelie_team7_complex_edge_flux.cpp` level structure also valuable:

```text
L1 coil-only
L2 one-way
L3 picard
```

Recommend keeping this structure.

Core solver still SPHinXsys/SYCL style:

```text
Interaction<Inner<>>
InteractionDynamicsCK
StateDynamics
DiscreteVariable
DelegatedData
SPH neighbor loop
```

No external mesh, CSR or CPU-only graph; correct direction.

---

### 10.2 Structural issues to fix immediately

#### Issue 1: output files get overwritten

Current filenames like:

```text
team7_probe_one-way.csv
team7_probe_phase90_xbands_one-way.csv
team7_plate_depth_profile_one-way.csv
```

baseline, `a_sign=-1`, `jscale_auto` etc. runs overwrite same files. Upload bundle CSV may not be the run user thinks; easily misleads analysis.

##### Cursor modification requirements

Add CLI:

```text
--team7-output-tag=<tag>
```

Output filenames must include tag:

```text
team7_probe_one-way_<tag>.csv
team7_probe_phase90_xbands_one-way_<tag>.csv
team7_plate_depth_profile_one-way_<tag>.csv
team7_bz_one-way_<tag>.csv
team7_validation_history_<tag>.txt
```

If user does not provide tag, auto-generate:

```text
tag = level + frequency + sign + phase90Convention + sourceScaleMode + jPostScaleMode
```

For example:

```text
oneway_f50_asign_p1_phase90_minus_source1_jraw
oneway_f50_asign_m1_phase90_plus_source0754_jauto
picard_f50_asign_p1_phase90_minus_source1
```

---

#### Issue 2: `j_post_scale` runs may still be mislabeled validation

If enabled:

```text
--team7-ind-j-post-scale=auto
```

must:

```text
diagnostic_only = 1
team7_validation_passed = 0
```

Even if RMS very low cannot pass.

##### Cursor modification requirements

In summary/history print both raw and post-scale:

```text
Bind_over_Bcoil_raw
Bind_over_Bcoil_post
P_recon_raw
P_recon_post
J_norm_raw
J_norm_post
phase90_RMS_raw
phase90_RMS_post
diagnostic_only=1
team7_validation_passed=0
```

Do not let post-scale data overwrite raw data.

---

#### Issue 3: `a_sign` parameter easily misused

`--ophelie-edge-flux-imag-a-sign=-1` changes solver equation, not probe comparison. Already seen it breaks constitutive consistency.

##### Cursor modification requirements

1. In CLI help and log label:

```text
--ophelie-edge-flux-imag-a-sign is debug-only.
Do not use it for TEAM7 validation phase convention.
```

2. TEAM7 validation entry should check:

```text
if imag_a_sign != +1:
    diagnostic_only = 1
    team7_validation_passed = 0
```

3. After adding `--team7-phase90-convention`, TEAM7 phase mapping only via that parameter.

---

#### Issue 4: `source_scale=0.754` should not default into validation

##### Cursor modification requirements

1. TEAM7 validation default:

```text
source_scale = 1.0
```

2. If:

```text
source_scale != 1.0
```

then:

```text
diagnostic_only=1
team7_validation_passed=0
```

3. Log warning: amplitude-fit diagnostic, not standard TEAM7 input.

---

#### Issue 5: `passed=1` naming too broad

Current L2 `passed=1` cannot mean TEAM7 validation passed; no phase90 gate, no source default / no post-scale / Picard etc.

##### Cursor modification requirements

Layer passes:

```text
smoke_passed
source_probe_passed
edge_flux_internal_passed
phase90_probe_passed
team7_validation_passed
diagnostic_only
```

Recommend initial definitions:

```text
smoke_passed = finite fields && no NaN && case completed
edge_flux_internal_passed = J_edge/σE audit ok && q antisym ok
source_probe_passed = L1 source-only RMS within tolerance
phase90_probe_passed = phase90 RMS within tolerance
team7_validation_passed =
    smoke_passed
    && source_probe_passed
    && edge_flux_internal_passed
    && phase90_probe_passed
    && !diagnostic_only
```

L2 one-way can have:

```text
one_way_probe_diagnostic_passed
```

But should not equal `team7_validation_passed`.

---

## 11. TEAM7 follow-on validation route

### P0: fix validation definition / scope first, not core solver

#### Step 1: freeze solver sign

```text
imag_a_sign=+1
```

`a_sign=-1` debug only, not validation.

#### Step 2: add phase90 convention

```text
--team7-phase90-convention=minus-imag|plus-imag
```

Re-run:

```bash
test_3d_ophelie_team7_complex_edge_flux \
  --team7-level=one-way \
  --team7-frequency=50 \
  --ophelie-edge-flux-complex=1 \
  --ophelie-edge-flux-imag-a-sign=1 \
  --team7-phase90-convention=minus-imag \
  --team7-coil-source-scale=1.0 \
  --team7-output-tag=oneway_f50_asign_p1_phase90_minus_source1
```

#### Step 3: fix output file tag

All CSV/log/history must carry tag to avoid overwrite.

#### Step 4: mark post-scale / source-scale / a-sign debug runs diagnostic-only

If any condition holds:

```text
j_post_scale != 1
source_scale != 1
imag_a_sign != +1
```

then:

```text
diagnostic_only=1
team7_validation_passed=0
```

---

### P1: fix TEAM7 L1 coil source

#### Step 5: implement or enable filament racetrack source

volume-racetrack + peak-scale cannot be standard TEAM7 coil source.  
Must implement/enable:

```text
TEAM7 filament racetrack source
```

Use real TEAM7 geometry and 2742 AT, no source_scale amplitude tuning.

Output:

```text
Bz_source_only_f0
Bz_phase0_total_f50
Bz_phase90_total_f50
```

If TEAM7 f=0 source-only reference available, prioritize f=0 comparison.

#### Step 6: L1 source validation

L1 gate should include:

```text
source_peak_ratio
source_RMS_abs
source_RMS_rel
source_peak_location_error
```

Do not use peak value only.

---

### P2: run true TEAM7 L3 Picard

After L1 source and phase90 convention clear, enter:

```text
TEAM7 L3 = complex edge-flux + A_ind Picard
```

Recommended command:

```bash
test_3d_ophelie_team7_complex_edge_flux \
  --team7-level=picard \
  --team7-frequency=50 \
  --ophelie-edge-flux-complex=1 \
  --ophelie-edge-flux-imag-a-sign=1 \
  --team7-phase90-convention=minus-imag \
  --team7-coil-source-scale=1.0 \
  --team7-output-tag=picard_f50_source1_phase90_minus
```

Picard output CSV:

```text
iter
J_rel_raw
J_rel_relaxed
P_recon
Aind_over_Acoil
Bind_over_Bcoil
Bz_phase0_RMS
Bz_phase90_RMS
max_J_real
max_J_imag
phi_eq_res_real
phi_eq_res_imag
```

Expect L3 Picard to significantly change `Bind/Bcoil`; must validate.

---

### P3: add 200 Hz and Jey reference

TEAM7 standard validation cannot look only at 50 Hz.  
Must add:

```text
--team7-frequency=200
```

Compare:

```text
Bz_200_phase0
Bz_200_phase90
```

And Jey reference:

```text
Jey_surface
Jey_depth_profile
Jey_probe_line
```

Jey comparison very important because it distinguishes:

```text
J itself too large
```

vs:

```text
Biot/probe coupling too strong
```

---

## 12. Direct execution summary for Cursor

Can send directly to Cursor:

```text
Current TEAM7 L2 one-way is only smoke passed, not TEAM7 validation passed.
phase90 not closed; j_post_scale=auto diagnostic only, not acceptance.

Modify in order:

P0-1. Freeze solver sign:
  imag a_sign default +1.
  --ophelie-edge-flux-imag-a-sign=-1 debug only.
  if imag_a_sign != +1 then diagnostic_only=1, team7_validation_passed=0.

P0-2. Add probe-layer phase90 convention:
  --team7-phase90-convention=minus-imag|plus-imag
  TEAM7 wt=90 sign mapping only at probe comparison layer, no longer via a_sign.

P0-3. Tag all TEAM7 output:
  --team7-output-tag=<tag>
  output:
    team7_probe_one-way_<tag>.csv
    team7_probe_phase90_xbands_one-way_<tag>.csv
    team7_plate_depth_profile_one-way_<tag>.csv
    team7_bz_one-way_<tag>.csv
    team7_validation_history_<tag>.txt
  prevent baseline / a_sign / jscale overwrite.

P0-4. j_post_scale=auto mark diagnostic-only:
  even if phase90 RMS good, cannot team7_validation_passed=1.
  print raw and post-scaled fields:
    Bind_over_Bcoil_raw/post
    P_recon_raw/post
    J_norm_raw/post
    phase90_RMS_raw/post

P0-5. source_scale default back to 1.0:
  source_scale != 1.0 diagnostic only.
  print warning:
    coil_source_scale != 1.0 is an amplitude-fit diagnostic and cannot be used as strict TEAM7 validation.

P0-6. layer passes:
  smoke_passed
  source_probe_passed
  edge_flux_internal_passed
  phase90_probe_passed
  team7_validation_passed
  diagnostic_only
  current L2 one-way at most one_way_probe_diagnostic_passed, not team7_validation_passed.

P1. Fix TEAM7 L1 source:
  implement/enable filament racetrack source.
  use 2742 AT and source_scale=1.0.
  compare source-only/f=0 or coil-only Bz; full-line RMS, peak, peak location.

P2. Run L3 Picard:
  after L1 source and phase90 convention clear, run complex edge-flux + A_ind Picard.
  output each Picard CSV:
    iter, J_rel_raw, J_rel_relaxed, P_recon,
    Aind/Acoil, Bind/Bcoil,
    Bz_phase0_RMS, Bz_phase90_RMS,
    max_J_real/imag, phi_eq_res_real/imag.

P3. Add 200 Hz and Jey reference:
  TEAM7 cannot be 50 Hz only.
  add --team7-frequency=200.
  add Jey reference comparison to judge J too large vs Biot/probe coupling too large.

Forbidden:
  Do not use TEAM7 50 kW power scaling;
  Do not use j_post_scale as validation;
  Do not use a_sign for TEAM7 phase90;
  Do not use source_scale=0.754 as standard input;
  Do not treat P_graph as physical power.
```

---

## 13. Final advancement route

TEAM7 should proceed:

```text
1. fix solver physical sign: imag_a_sign=+1
2. handle phase90 convention at probe layer
3. fix output tag and pass layering to avoid false positives
4. L1 coil source with source_scale=1.0 + filament racetrack re-validation
5. L2 one-way as diagnostic, not final validation
6. L3 complex edge-flux + A_ind Picard as standard TEAM7 validation main line
7. add 200 Hz and Jey reference
8. only then add phase90 RMS to hard gate
```

This round TEAM7 is not failure, but exposed key engineering and physical definition / scope issues required by standard benchmark. Following above route, TEAM7 can become a very strong standard validation case for complex edge-flux solver.
