# OPHELIE-like TEAM7 Native STL Test Review and Next-Steps Correction Plan

## 1. Overall Conclusion

This TEAM7 native STL test shows:

- `PlateBody` and `CoilSourceBody` reload successful, particle counts `n_plate=7700`, `n_coil=7712` consistent with log;
- `--no-phi` path as expected only runs `sigma -> coil source -> Biot-Savart -> Level0 E/J/Q -> power scaling`;
- `passed=1` only means smoke pipeline runs, fields nonzero, power scaling finite;
- cannot yet say "already quantitatively validated against TEAM7/COMSOL";
- current main bottlenecks are not relax, but **reference comparison sample points not on PlateBody**, **coil J0 area definition still coarse / may underestimate ~7–8×**, **power scaling breaks direct amplitude comparison with 1 A/turn reference**, **phi/self-induction still not usable TEAM7 quantitative path**.

Most important: **TEAM7 Bz A1–B1 reference line is in air, z=34 mm, above aluminum plate, below coil; current code only computes B/A/E/J/Q on PlateBody particles, therefore cannot interpolate Bz from current PlateBody VTP to compare A1–B1 reference.** Must add `ProbeBody` or host-side probe evaluator, sample Biot-Savart directly at probe points.

---

## 2. Are Current Case Operation Details Reasonable?

## 2.1 Geometry and reload

Current native STL setup is reasonable:

- STL coordinates in mm;
- `OphelieTeam7NativeMesh::stl_scale_to_meter_ = 1e-3`, imported to EM test as m;
- native small air box: `[-50,350] x [-50,350] x [-50,200] mm`, only system bounds / relax air, not EM PDE air;
- EM only uses `PlateBody` and `CoilSourceBody`.

But note hidden pitfall:

- standalone `particle_generation_em.cpp` prototype STL scale is `1.0`, domain/dp also mm, body names `Coil`, `Plate`, `Air`;
- `test_3d_ophelie_team7.cpp` native STL mode scale is `1e-3`, body names `CoilSourceBody`, `PlateBody`;
- so standalone `particle_generation_em` generated `Reload.xml` **cannot directly be used by** `test_3d_ophelie_team7 --reload=1` unless unit m conversion and body name rename already done.

Recommended fixed workflow:

```bash
./test_3d_ophelie_team7 --native-stl --relax=1 ...
./test_3d_ophelie_team7 --native-stl --reload=1 ...
```

External relax case if integrated must satisfy:

```text
body name: PlateBody / CoilSourceBody
position unit: meter
VolumetricMeasure unit: m^3
geometry exactly matches native STL scale and dp
```

## 2.2 Current `--no-phi` EM computation chain

Current run chain:

```text
Assign sigma on PlateBody
Initialize JSrcReal on CoilSourceBody: J0 * e_theta
Compute coil -> plate Biot-Savart: ACoilReal, BCoilReal
Combine -> ASrcReal, BSrcReal
Level0: EImag = -omega * ASrcReal
JImag = sigma * EImag
Q = 0.5 * JImag · EImag
Power scaling: A/B/E/J/Phi scaled by field_scale, Q scaled by power_scale
VTP/log output
```

Reasonable for smoke.

Two limits for TEAM7 reference comparison:

1. `Bz` reference points in air line, not on PlateBody; current PlateBody fields not directly usable.
2. After power scaling `max_BSrc` etc. no longer correspond to physical amplitude of 2742 turns × 1 A/turn.

---

## 3. Current Main Errors, Bottlenecks and Risks

## 3.1 P0: Bz reference line not actually sampled

`TEAM7_Bz_A1_B1_reference_mT.csv` points:

```text
x = 0..288 mm
y = 72 mm
z = 34 mm
```

Aluminum plate STL bbox:

```text
x: 0..294 mm
y: 0..294 mm
z: 0..19 mm
```

Coil STL bbox:

```text
x: 94..294 mm
y: 0..200 mm
z: 49..149 mm
```

So Bz A1–B1 lies in **air between plate and coil**, not conductor particle positions.

Current `ComputeOphelieCoilToGlassBiotSavartCK` only computes B on plate particles. Even if VTP has `BSrcReal`, it is B on plate interior/surface particles, not reference line B.

Conclusion:

```text
Do not do Bz reference CSV comparison before adding ProbeBody / ProbeEvaluator.
```

Recommend implementation:

- `ProbeBody`: 17 points or line sample points, coordinates from CSV mm->m;
- register only `ACoilReal/BCoilReal/ASrcReal/BSrcReal` on ProbeBody;
- same Biot-Savart source kernel, compute coil -> probe;
- when self-induction available later, compute plate J -> probe induced B;
- output `Bz_raw_mT`, then compare with CSV.

## 3.2 P0: power scaling must be disableable

Current `ScaleOphelieElectromagneticFieldsCK` already scheme B, synchronously scales A/B/E/J/Phi/Q, fields internally self-consistent. Good for thermal coupling.

But TEAM7 reference is 2742 turns, 1 A/turn EM amplitude. If `target_P=50kW` scaling enabled, printed `B/E/J` are under equivalent current `I_eff`, not comparable to COMSOL/muFEM raw reference.

Must add CLI:

```text
--no-power-scaling
```

or:

```text
--target-power=0
```

And:

```cpp
if (target_joule_power_ <= 0 || !enable_power_scaling) {
    power_scale = 1.0;
    field_scale = 1.0;
    effective_current_amplitude = current_amplitude_;
}
```

TEAM7 Bz comparison command must use:

```bash
--no-power-scaling --no-phi
```

## 3.3 P0: native coil J0 area definition still not accurate enough

Current native derived geometry uses:

```cpp
radial_extent = max(coil_bbox_extent_x, coil_bbox_extent_y)
A_cross = radial_extent * coil_extent_z
```

For STL coil, bbox extent x/y ~200 mm, z extent 100 mm, so A_cross ~`0.02 m^2`, J0 ~`2742/0.02 = 1.37e5 A/m^2`.

From STL geometry, coil is annular/ring-like:

```text
bbox center approx: (194,100,99) mm
radial range approx: 90..121 mm
z range: 49..149 mm
STL volume approx: 0.001588 m^3
```

If estimate cross section axisymmetric ring:

```text
A_cross ≈ volume / (2*pi*r_mean)
        ≈ 0.001588 / (2*pi*0.105)
        ≈ 0.0024 m^2
```

then:

```text
J0 ≈ 2742 / 0.0024 ≈ 1.1e6 A/m^2
```

~7–8× current J0. Log scaled `max_BSrc≈9.74`, `field_scale≈7895`, so raw `max_BSrc≈1.23 mT`; after J0 fix 7–8× raw B magnitude ~`9–10 mT`, comparable to A1–B1 CSV mT level. J0 definition likely main cause of current B amplitude too small.

Recommend change native `coil_current_cross_section_m2_` to particle volume integral estimate:

```cpp
A_cross = sum_j Vol_j / (2*pi*r_j)
```

where:

```cpp
r_j = sqrt((x_j - coil_center_x)^2 + (y_j - coil_center_y)^2)
```

Much more reasonable than bbox area for STL/reload particles. For analytic annulus continue:

```cpp
A_cross = (Rout - Rin) * H
```

## 3.4 P1: no-phi Bz can only compare coil field, not full TEAM7 EM

Under `--no-phi`:

```text
B = B_coil only
E/J/Q from A_coil only
no induced current reaction field
```

Means:

- `Bz` does not vary with frequency;
- `Bz_phase90` basically zero;
- cannot capture B phase and frequency difference from eddy current in COMSOL/muFEM.

CSV 50/200 Hz real Bz difference small but exists; phase90 max ~1.4 mT, induced reaction not zero.

So P0 Bz comparison goal only:

```text
validate coil geometry + J0 + Biot-Savart main real Bz magnitude and sign.
```

Cannot claim full TEAM7 EM validation.

## 3.5 P1: phi currently no help for Bz reference, key for Jy reference

`phi` correction affects:

```text
E, induced J, JouleHeat, divJ
```

But without self-induction, `phi` does not change `Bsrc` because `Bsrc` only from coil Biot-Savart. Therefore:

- for Bz line: phi not needed first;
- for Jy line / Joule heat: phi necessary;
- for B phase90 / frequency dependence: need self-induction, i.e. induced current -> B probe.

Current native 7700 particles phi unstable/non-convergent, do not block TEAM7 first stage on phi.

## 3.6 P1: self-induction still not for reference

`ComputeOphelieGlassSelfInducedBiotSavartCK` still writes A from `JImag` integral into `AIndReal`, but in frequency domain should write `AIndImag`. Therefore:

```text
--self-induction still smoke only, not for TEAM7 B phase90/reference validation.
```

After fix also need:

```text
coil -> probe B real
plate JImag -> probe BImag
combine BReal/BImag
```

## 3.7 P1: frequency CLI missing

Native TEAM7 default `frequency=50 Hz`, CSV has 50/200 Hz. Current CLI no `--frequency=`.

Recommend add:

```text
--frequency=50
--frequency=200
```

Note: without self-induction `B_coil` frequency-independent; frequency only affects `E/J/Q`.

---

## 4. Recommended Next-Step Route

## 4.1 Fix/add P0 immediately

### Task 1: add `--no-power-scaling` / `--target-power=`

Goal: reference comparison uses raw fields.

Recommended command:

```bash
./test_3d_ophelie_team7 \
  --native-stl --reload=1 --no-phi --sigma=1e4 \
  --no-power-scaling --state_recording=1
```

### Task 2: fix native J0 cross-section

Implement:

```cpp
inline Real particleAxisymmetricCurrentCrossSectionArea(BaseParticles &coil_particles, const Vecd &coil_center)
{
    sync Position and VolumetricMeasure to host if needed;
    Real area = 0.0;
    for each coil particle j:
        Real r = sqrt((x_j-cx)^2 + (y_j-cy)^2);
        area += Vol_j / (2.0 * Pi * max(r, TinyReal));
    return max(area, TinyReal);
}
```

Then:

```cpp
J0 = Nturns * I_per_turn / A_cross;
```

Log must print:

```text
coil_cross_section_area
J0
coil_volume
mean_radius
```

### Task 3: add probe sampling, do not interpolate Bz from PlateBody

Implement `Team7BzProbeBody` or host evaluator:

```text
17 points from CSV:
(x_mm, 72 mm, 34 mm) * 1e-3
```

First compute coil-only raw B:

```text
BProbeCoilReal
```

Output CSV:

```text
x_mm, Bz_raw_mT, Bz_ref_50Hz_phase0_mT, abs_err_mT, rel_err
```

P0 ignore phase90 first.

### Task 4: add `--compare-team7-bz`

Read:

```text
TEAM7_Bz_A1_B1_reference_mT.csv
```

Do not hardcode build-relative fragile path; support:

```text
--team7-reference-dir=...
```

## 4.2 P0 acceptance criteria

First stage only accept coil-source Bz:

```text
mode: native STL + reload + no-phi + no-power-scaling
frequency: 50 Hz, but B_coil itself frequency-independent
compare: Bz_phase0_mT only
metrics: signed shape + magnitude order
```

Recommended temporary thresholds:

```text
RMS relative error < 50%: can continue optimization
RMS relative error 20–40%: current OPHELIE-light acceptable to enter next stage
RMS relative error < 20%: coil/J0/source geometry very good
```

Do not set 10% initially. Current model has no induced reaction B.

## 4.3 P1: phi/Jy direction

After P0 Bz raw probe quantified, then Jy:

- must add `PhiImag`;
- physical Al `sigma=3.54e7` may be too stiff, continuation from `1e4 -> 1e6 -> 1e7 -> 3.54e7`;
- compare trends first, do not directly require COMSOL Jy values.

Acceptance still:

```text
divJ_phi < divJ_level0
```

But for TEAM7 native first as warning, not passed hard fail.

## 4.4 P2: self-induction/phase90

To compare B phase90 or frequency difference, must:

1. fix `JImag -> AIndImag/BIndImag`;
2. compute plate-induced B at probe points;
3. synthesize `BProbeReal/BProbeImag`;
4. then compare CSV phase0/phase90.

---

## 5. Current Command Recommendations

### 5.1 reload/no-phi smoke check only

```bash
cd /home/yyc/SPHinXsysSYCL/build
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/test_3d_ophelie_team7 \
  --native-stl --reload=1 --no-phi --sigma=1e4 \
  --ophelie-smoke --state_recording=0
```

### 5.2 generate visualization VTP, still no reference

```bash
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/test_3d_ophelie_team7 \
  --native-stl --reload=1 --no-phi --sigma=1e4 \
  --ophelie-smoke --state_recording=1
```

Note: current B/E/J/Q scaled, not comparable to COMSOL.

### 5.3 Bz raw compare after no-power-scaling fix

```bash
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/test_3d_ophelie_team7 \
  --native-stl --reload=1 --no-phi --sigma=1e4 \
  --no-power-scaling --compare-team7-bz \
  --team7-reference-dir=/home/yyc/SPHinXsysSYCL/docs/TEAM7-reference/reference_data/team7 \
  --state_recording=1
```

---

## 6. Direct Task Manifest for Cursor

Can send directly to Cursor:

```text
Current TEAM7 native reload smoke passed, but cannot do reference validation yet. Next steps by priority:

P0-1: Add --no-power-scaling or --target-power=0; when target<=0 power_scale=field_scale=1, I_eff=I0. TEAM7 reference comparison must disable 50 kW scaling.

P0-2: Fix native coil J0 area. Do not use bbox radial_extent*z_extent. For reload/STL particles use A_cross=sum(Vol/(2*pi*r)), r = particle xy radius to coil_center. Print A_cross, coil_volume, mean_radius, J0.

P0-3: Add TEAM7 Bz probe sampling. A1-B1 reference points in air: (x_mm,72,34)*1e-3, cannot interpolate from PlateBody VTP. Add ProbeBody or host evaluator, coil->probe Biot-Savart at these points, output Bz_raw_mT and CSV phase0 comparison.

P0-4: Add --compare-team7-bz and --team7-reference-dir. Compare 50Hz phase0 only first, no-phi, no-power-scaling. RMS rel error smoke threshold 50%, 20-40% iterable, <20% very good.

P1: Phi/Jy do not block Bz P0. Phi only affects J/Q, not coil-only Bz. After Bz raw probe OK, handle PhiImag convergence and Jy line.

P2: self-induction still experimental. For phase90 or 200Hz Bz, first fix JImag -> AIndImag/BIndImag, superpose plate-induced B at probe points.

P3: Add --frequency=, but note no-self B_coil frequency-independent.
```

---

## 7. Final Judgment

Current case operation reasonable as engineering smoke; reload, no-phi, Biot, Level0, scaled heat source pipeline already runs. Real bottlenecks are "comparison definition" not yet implemented: reference Bz line not on plate, power scaling not disabled, J0 area underestimates source strength, self-induction/phi cannot yet support full TEAM7 phase/eddy comparison.

Next steps: do not blindly run more reload. Should first make **raw Bz probe comparison**; that will first tell how far coil source + Biot-Savart is from COMSOL/muFEM main magnetic field.
