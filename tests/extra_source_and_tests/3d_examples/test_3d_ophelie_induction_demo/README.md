# OPHELIE Induction Welding Demo (Case A)

**Status:** Gate 1–3 (geometry + Biot `|B_src|` + edge-flux φ / `|J|` / `q_J` → **VTP**)  
**Nature:** demonstration — **not** validated for welding process design.  
**Approximation:** source-field only (`A_src`, `B_src` ← coil Biot–Savart). No workpiece magnetic feedback.

## Build

From the existing Release SYCL build tree:

```bash
cd /path/to/SPHinXsys/build
cmake --build . --target test_3d_ophelie_induction_demo -j$(nproc)
```

## Run

默认写出完整场 VTP（`--stage=em`）：

```bash
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_induction_demo/bin/test_3d_ophelie_induction_demo \
  --case=a --stage=em
```

仅几何 / 仅 Biot：`--stage=geometry` / `--stage=biot`  
脚本：`scripts/run_case_a_em.sh`（会把 VTP 拷到 `docs/ophelie/welding/case_a_vtp/`）

## Case A parameters (SI)

| Quantity | Value |
|---|---:|
| Workpiece | Al cylinder, non-magnetic, constant σ |
| `R_work` | 0.020 m |
| `L_work` | 0.100 m (centered at origin, axis z) |
| `R_coil` | 0.035 m |
| Turns | 5 |
| Pitch | 0.012 m |
| `I_peak` | 500 A (peak phasor) |
| `f` | 10 kHz (stored; unused until Gate 3) |
| `σ` | 3.5×10⁷ S/m |
| `dp` | 0.004 m (smoke mesh; skin depth ≪ dp) |

## Outputs (under binary `output/`；归档 `docs/ophelie/welding/case_a_vtp/`)

- `WorkBody_ite_0000000000.vtp` — `BSrcMag`, `JEdgeReconMag`, `JouleHeatEdgeReconComplex`, `Phi*`, …
- `CoilFilament_case_a.vtp` — 5 匝线圈折线
- `case_a_axis_Bz.csv` / `case_a_params.csv`

## Acceptance

- Gate 2：有限 `|B_src|`；轴上 rel_err&lt;2%；`I→2I` ≈ 2
- Gate 3：有限 `P_recon`、`|J|`、`q_J≥-tol`；φ 方程残差有限（演示级，非工艺验证）

## Allowed / forbidden claims

Allowed: demo of coil source magnetic field on an Al cylinder geometry.  
Forbidden: welding success, melt temperature, validated process prediction.


## Case B (pipe joint)

```bash
./test_3d_ophelie_induction_demo --case=b --stage=em --coil-z-center=0
./test_3d_ophelie_induction_demo --case=b --stage=em --coil-z-center=0.025
# scripts/run_case_b_em.sh → docs/ophelie/welding/case_b_vtp/
```

Continuous hollow Al pipe + joint marker (not a weld). Fields: `JointMarker`, `PipeSegment`, `P'(z)` CSV, `eta_joint`.
