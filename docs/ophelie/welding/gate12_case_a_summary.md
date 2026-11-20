# Gate 1–2 Case A Summary

**Date:** 2026-10-10  
**Status:** **PASS**  
**Nature:** demo / not validated for process design  
**Approximation:** coil Biot–Savart source field only (no φ / J / Joule yet)

## Commands

```bash
cd /home/ubuntu/SPHinXsys/build
cmake --build . --target test_3d_ophelie_induction_demo -j$(nproc)

# or:
../tests/extra_source_and_tests/3d_examples/test_3d_ophelie_induction_demo/scripts/run_case_a_gate12.sh geometry
../tests/extra_source_and_tests/3d_examples/test_3d_ophelie_induction_demo/scripts/run_case_a_gate12.sh biot
```

Binary cwd output:  
`build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_induction_demo/bin/output/`

Archived smoke: `docs/ophelie/welding/gate12_smoke/`

## Results

| Check | Value |
|---|---|
| Particles | 1971 (`dp=0.004 m`) |
| Radial gap coil–work | 0.015 m |
| Axis filament vs analytic `max_rel_err` | ≈4.3×10⁻⁴ |
| `max \|B_src\|` | ≈0.037 T |
| `I→2I` volume-norm ratio | **2.0** |
| Finite fields | yes |
| `passed` | **1** |

## Artifacts

- `WorkBody_ite_0000000000.vtp` — `BSrcReal`, `BSrcImag`, `BSrcMag`, `ASrc*`, `Sigma`
- `CoilFilament_case_a.vtp`
- `case_a_axis_Bz.csv`
- `case_a_gate12_params.csv`
- logs: `case_a_geometry.log`, `case_a_biot.log`

## Known limits

- Skin depth at 10 kHz / Al ≪ `dp` — Gate 3 Joule maps will be qualitative only on this mesh.
- No workpiece reaction field (`A_ind=0`).
- Case B not started.

## Next (Gate 3 only when instructed)

Edge-flux complex φ → recon `J` / `JouleHeatEdgeReconComplex` with residual gates; keep B-only slides if φ fails.
