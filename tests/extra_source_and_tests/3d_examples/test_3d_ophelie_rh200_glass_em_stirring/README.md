# test_3d_ophelie_rh200_glass_em_stirring

Midterm coupled prototype: **RH200 STL stirring geometry** + (later) OPHELIE Joule heating.

## Step 0 scope

- Single-phase **glass** fluid (`oil.stl` → glass domain), **wall**, **rotor** (Simbody spin).
- **Thermal**: glass inner diffusion + **glass–rotor** contact heat transfer.
- **Wall**: hydrodynamic boundary only — **no** glass–wall thermal contact (adiabatic wall).

## Step 3 scope (EM fixed JouleHeat + stirring flow)

Single command: OPHELIE EM solve once, then 1 s stirring with frozen spatial JouleHeat:

```bash
# 1) Re-relax at dp=0.008 if reload was from another dp
./test_3d_ophelie_rh200_glass_em_stirring --relax=1 --dp=0.008 --geometry-scale=1.0

# 2) Coupled run (EM + flow + thermal)
./test_3d_ophelie_rh200_glass_em_stirring --reload=1 --run=1 --dp=0.008 --end-time=1 \
  --joule-mode=em-fixed --target-power=50000 --state-recording=1
```

Main loop: flow → rotor → **frozen JouleHeat → T** → thermal diffusion.

### Resolution (`--dp`)

| Body | How `dp` is used |
|------|------------------|
| **glass** | `SPHSystem` global spacing; lattice at relax, reload fixes particle count |
| **wall** | Same global `dp` + extrude shell `4*dp` in `WallBoundary` |
| **rotor** | Same global `dp` (lattice / level-set relax) |

**Important**: `--relax`, `--reload`, `--em-solve`, and `--run` must all use the **same** `--dp`. Default is now **0.008**. If you change `dp`, re-run `--relax=1`.

## Step 2 scope (OPHELIE EM solve once)

After `--relax=1` + `./reload`, run complex edge-flux EM on glass particles:

```bash
./test_3d_ophelie_rh200_glass_em_stirring --reload=1 --em-solve=1 --dp=0.008 --target-power=50000
```

- Glass particle bbox → auto multiloop coil (`--coil-radius-factor=1.15`, `--coil-turns=6`, etc.)
- OPHELIE: Biot-Savart A_coil → complex edge-flux φ → E/J → JouleHeat, coil-current calibration to `target-power`
- Outputs: `rh200_coil_geometry_log.csv`, `rh200_em_solve_summary.csv`, glass VTP (`Sigma`, `JouleHeat`, `EReal`, `EImag`, `JReal`, `JImag`, `PhiReal`, `PhiImag`)

**Note**: EM path uses `SolidBody` + OPHELIE API (same `GlassGeometry` reload name). Flow run (`--run=1`) still uses `FluidBody`; Step 3 will couple EM JouleHeat into stirring.

## Flow Joule modes (production)

- `--joule-mode=off|em-fixed|em-grid`
- `em-grid` (recommended): fixed-space Euler grid sampling each step for `JouleHeat` and `JImag/JReal`.
- `em-fixed` (debug only): frozen particle fields from one-shot EM solve.

**Execution policy (Step 0)**:

- `--relax=1` particle relaxation: **SYCL/GPU** (`par_ck` level-set + `KernelGradientIntegral` / `PositionRelaxationCK`; wall/rotor relaxed; glass lattice only).
- `--reload=1 --run=1` flow + thermal: **SYCL/GPU** (`MainExecutionPolicy` / `par_ck`).

## Geometry / units

STL from RH200 (`data/` → build `input/`). ParaView audit (meters):

| STL | Approx. bounds (m) |
|-----|-------------------|
| `oil.stl` (glass fluid) | x: [-0.39, 0.36], y: [-0.36, 0.39], z: [0.024, 0.55] |

Default: `--geometry-scale=1.0` (meters).

## Build

From repo `build/`:

```bash
ninja test_3d_ophelie_rh200_glass_em_stirring
```

Run from executable directory (CMake sets `bin/` as working directory).

## Run (three-phase plan; Step 0 = relax + flow)

```bash
cd build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_rh200_glass_em_stirring/bin

# Phase 0: particle relaxation (wall + rotor; glass lattice + reload)
./test_3d_ophelie_rh200_glass_em_stirring --relax=1 --dp=0.005 --geometry-scale=1.0

# Phase 2 (Step 0 flow): reload + stirring + thermal
./test_3d_ophelie_rh200_glass_em_stirring --reload=1 --run=1 --dp=0.005 --state-recording=1

# Phase 3 (EM + stirring + thermal, recommended)
./test_3d_ophelie_rh200_glass_em_stirring --reload=1 --run=1 --dp=0.005 \
  --preset=rh200-french-like-full --end-time=60 \
  --joule-mode=em-grid --target-power=50000 \
  --state-recording=1 --state-record-interval=1.0
```

Particle count is large at `dp=0.005`; increase `dp` for faster smoke tests.

## Outputs

- `./reload/` — relaxation particle data
- `./output/rh200_geometry_audit.csv` — STL bounding boxes
- `./output/rh200_energy_budget.csv` — Joule power / thermal energy (Step 1, when `--joule-mode` enabled)
- VTP: glass `Temperature`, `Pressure`, `JouleHeat`, `JReal`, `JImag`; rotor proxy surface (if `--state-recording=1`)

## Later steps

- Periodic EM update (`--em-update-interval=...`) when σ(T) feedback is needed

See `OPHELIE_RH200_GLASS_EM_STIRRING_CURSOR_TASKS.md` at repo root.
