# Stage 10.13 — cold crucible demonstration case record

> **Date**: 2026-05-21  
> **Nature**: pipeline smoke test / **not quantitative validation**  
> **Handoff**: [`CURSOR_APHI_STAGE10_13_CHATGPT_HANDOFF.md`](CURSOR_APHI_STAGE10_13_CHATGPT_HANDOFF.md)

---

## 1. New files

| File | Purpose |
|------|---------|
| [`diagnostics/aphi_cold_crucible_demo_helpers.h`](../extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_cold_crucible_demo_helpers.h) | impressed A, E/J/Joule/B/H, VTP, gate |
| [`test_3d_aphi_ck_cold_crucible_demo/`](test_3d_aphi_ck_cold_crucible_demo/) | executable test + README |
| this document | Stage 10.13 P1/P2 record |

---

## 2. Geometry and materials

- **Single body** `AphiLhsTestBody` + cold crucible **region tagging** (reuses `AphiColdCrucibleUnitBoxLayout`)
- MVP-A **does not** use multibody Contact (impressed field + post-processing does not depend on Contact apply)
- Materials (demo defaults):

| Region | sigma | nu | MaterialRegionId |
|------|-------|-----|------------------|
| air | 1e-8 | 1 | 0 |
| crucible wall | 0.1 | 1 | 1 |
| melt | 1.0 | 1 | 2 |
| coil | 1e-8 | 1 | 3 |

`dp=0.1`, `omega=100` (nondimensional demo values).

---

## 3. Impressed A formula

Axisymmetric toroidal vector potential (Cartesian):

```text
r = sqrt(x^2 + y^2)
A_theta(r,z) = A0 * exp(-((r - r_coil)^2 / w_r^2 + (z - z_center)^2 / w_z^2))
A_x = -A_theta * y / r
A_y =  A_theta * x / r
A_z = 0
phi = 0
```

Defaults: `A0=0.05`, `r_coil=0.5*L`, `w_r=0.12*L`, `w_z=0.2*W`, `z_center=0.5*W`.

---

## 4. Post-processing chain (CK)

```text
AssignImpressedCoilVectorPotentialCK
-> B = curl A  (LinearCorrectionMatrix + LinearGradient + AphiVectorGradientCurlCK)
-> H = nu * B
-> grad phi (Inner) -> E = -i omega A - grad phi -> J, Joule
-> magnitude scalars -> BodyStatesRecordingToVtp
```

---

## 5. Output fields (VTP)

`AReal/AImag`, `PhiReal/PhiImag`, `BReal/BImag`, `HReal/HImag`, `ElectricField*`, `CurrentDensity*`, `JouleHeatSource`, `Sigma`, `Nu`, `AMagnitude`, `BMagnitude`, `EMagnitude`, `JMagnitude`, `HMagnitude`, `MaterialRegionId`.

Path (run from `bin/`): `output/AphiLhsTestBody_ite_0000000000.vtp`

---

## 6. Terminal acceptance (user reproducible)

```bash
cd build
ninja test_3d_aphi_ck_cold_crucible_demo
cd tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_cold_crucible_demo/bin
./test_3d_aphi_ck_cold_crucible_demo
```

Example output (2026-05-21):

```text
test_3d_aphi_ck_cold_crucible_demo mode=impressed dp=0.1 particles=1000 omega=100
  max_A=0.0469 max_B=0.256 max_E=4.69 max_J=1.29 max_Joule=1.10
  melt_Joule_integral=0.00249 melt_J_max=1.29
  demonstration_only=1 passed=1
```

Gate: global max finite and >0; melt `J`/Joule integral > threshold; air `J`/Joule below demo tolerance (`sigma_air=1e-8`).

---

## 7. Explicitly not validated

- GMRES coil–melt coupling solve (MVP-B `--mode solve` not implemented)
- Real slit geometry, TEAM7 quantitative accuracy, projection, three-body Contact A-penalty

---

## 8. Recommended next steps

1. **P4** InnerOnly graddiv PC consistency supplementary test  
2. **P3** B dp=0.05 refinement (optional)  
3. MVP-B: add `--mode solve` in same directory (optional, does not block demo)  
4. Multibody Contact + impressed/solve (later)

---

## 9. Status summary for ChatGPT

```text
Stage 10.13 completed a cold-crucible-like A-phi electromagnetic demonstration case. The case uses analytic/impressed vector potential to represent outer coil excitation, outputs A/phi/B/H/E/J/Joule, and can display air-domain magnetic field, conductive melt induced current, and Joule heating distribution in ParaView. This case is a pipeline smoke test, not quantitative validation. Existing MMS gates are not broken; Contact A-penalty remains research-only, production lambda_A off; projection continues defer.
```
