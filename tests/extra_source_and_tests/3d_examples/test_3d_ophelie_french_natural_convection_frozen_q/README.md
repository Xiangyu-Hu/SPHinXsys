# test_3d_ophelie_french_natural_convection_frozen_q

Stage 4.1 frozen `Q_50kW` (default) **or** periodic thermo–EM coupling.

- Default `--thermo-em-coupling=off`: one-shot EM, Q copied onto particles (regression).
- `--thermo-em-coupling=periodic`: T→σ_bar(r,z)→EM→Euler Q, updates on **physical time**.

Natural absorbed-power target remains **50 kW** (not 60 kW). Vertical axis is **+z**.

WCSPH hydro: the melt top is a free surface (no lid). The crucible wall is an open cup with a short **rim above the melt** so startup slosh stays inside. Density is initialized hydrostatically, then gravity-only settle runs before Q / Boussinesq. `flow_ok` now requires `U_max ≤ c0` and escaped fraction ≤ 1%.

Full notes: [`docs/ophelie/10_THERMO_EM_COUPLING_ROUND1.md`](../../../../docs/ophelie/10_THERMO_EM_COUPLING_ROUND1.md).

## Defaults

| Flag | Default |
|------|---------|
| `--thermo-em-coupling` | `off` |
| `--em-update-interval` | 10 s |
| `--sigma-under-relaxation` | 0.3 |
| `--em-control` | `fixed-current` |
| `--aind` | frozen-Q: on; periodic: **off** unless `--aind=on` |
| `--enable-thermal-diffusion` | **periodic: on**; frozen-Q: off |
| `--no-thermal-diffusion` | keep Q + BC only on the periodic path |
| `--q-conservative-remap` | off (log errors; do not silently scale Q) |
| `--c0=` | `10√(gH)` (~13.5 m/s at H=0.185 m) |
| `--u-ref=` | `√(gH)` (so advection starts immediately) |
| `--wall-rim=` | `6 dp` (open rim, not a lid) |
| `--wall-thickness-factor=` | `4` (dummy wall width in dp; `2` leaks through the floor) |
| `--hydro-settle-time=` | 0.15 s (clock reset to 0 afterwards) |
| `--state-record-interval=` | 0 (no VTP) |

## Frozen-Q regression (cwd=`build/`)

```bash
cmake --build . --target test_3d_ophelie_french_natural_convection_frozen_q -j$(nproc)

RELAX_RELOAD=./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_glass_relax/bin/reload

./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/bin/test_3d_ophelie_french_natural_convection_frozen_q \
  --reload-dir="$RELAX_RELOAD" \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-edge-flux-complex=1 \
  --target-power=50000 \
  --end-time=0.02
```

## Periodic coupling smoke

```bash
# from the test bin/ directory; interval < end-time so several EM updates occur
./test_3d_ophelie_french_natural_convection_frozen_q \
  --reload-dir=../../test_3d_ophelie_french_natural_glass_relax/bin/reload \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-edge-flux-complex=1 \
  --thermo-em-coupling=periodic --em-update-interval=0.005 \
  --em-control=fixed-current --aind=off \
  --target-absorbed-power=50000 --end-time=0.03
```

CSV: `output/french_thermo_em_coupling.csv`, `french_sigma_rz.csv`, `french_q_rz.csv`, `french_energy_budget.csv`, `french_natural_convection_monitor.csv`, `french_spatial_stats.csv`, `french_T_rz.csv`. SSH: keep `--state-record-interval=` at 0 (default); do not download VTP.

Hydro check (expect `escaped=0/...`, `U_max` well below `c0`, not `U_max ≈ g t`):

```bash
./test_3d_ophelie_french_natural_convection_frozen_q \
  --reload-dir=../../test_3d_ophelie_french_natural_glass_relax/bin/reload \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-edge-flux-complex=1 \
  --thermo-em-coupling=periodic --em-update-interval=0.05 \
  --em-control=fixed-current --aind=off \
  --target-absorbed-power=50000 --end-time=0.5
```

Periodic coupling now enables SPH conduction (constant k at 1473 K). Add `--no-thermal-diffusion` to match the old Q+BC-only path. Energy residual is written but **not** a pass gate.

Conduction + longer time (still not a literature cell; watch `escaped`, `U_max` vs `U_buoy`, `dT`):

```bash
./test_3d_ophelie_french_natural_convection_frozen_q \
  --reload-dir=../../test_3d_ophelie_french_natural_glass_relax/bin/reload \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-edge-flux-complex=1 \
  --thermo-em-coupling=periodic --em-update-interval=0.5 \
  --em-control=fixed-current --aind=off \
  --target-absorbed-power=50000 --end-time=5
```
