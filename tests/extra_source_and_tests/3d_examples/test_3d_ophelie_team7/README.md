# test_3d_ophelie_team7

OPHELIE TEAM7 case in two geometry modes:

| Mode | Flag | Geometry |
|------|------|----------|
| Analytic scaffold | (default) | 1 m cube + cylinder plate + annular coil |
| **Native TEAM7 STL** | `--native-stl` | `coil.stl` / `plate.stl` from [particle_generation_em](../particle_generation_em) (mm mesh, scaled to **m**) |

Native mode matches [COMSOL TEAM Problem 7](https://www.comsol.com/blogs/solving-team-problem-7-acdc-module/) (2742 turns, 1 A/turn, 50/200 Hz benchmark). Probe lines and reference CSVs live under `../reference_data/team7/` (validation not wired yet).

## Build

```bash
cd build
ninja test_3d_ophelie_team7
```

## Run

**Run from `~/SPHinXsysSYCL/build`** (do not run from `.../bin/` or `output/` will land in `bin/output/`).

Executable path (copy-paste ready):

```text
/home/yyc/SPHinXsysSYCL/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/test_3d_ophelie_team7
```

### Native STL (real geometry, `particle_generation_em`)

STLs are copied at build time to `.../test_3d_ophelie_team7/bin/input/`; run from **`~/SPHinXsysSYCL/build`** (auto-finds `bin/input`).

```bash
cd /home/yyc/SPHinXsysSYCL/build

# 1) Lattice + EM smoke only (dp_mm=6 → dp=0.006 m, ~7.7k plate particles)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/test_3d_ophelie_team7 \
  --native-stl --skip-relax --no-phi --sigma=1e4 --ophelie-smoke --state_recording=0

# 2) SYCL body-fitted relax (default includes small air; Reload writes CoilSourceBody + PlateBody only)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/test_3d_ophelie_team7 \
  --native-stl --relax=1 --state_recording=0 --relax-log-every=100 --sigma=1e4

# 3) Run EM on relaxed particles (same as mr cylinder case: --relax=1 first, then --reload=1)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/test_3d_ophelie_team7 \
  --native-stl --reload=1 --sigma=1e4 --ophelie-smoke --state_recording=1
```

Options:

| Option | Meaning |
|------|------|
| `--native-dp-mm=6` | Reference particle spacing (mm); default 6 |
| `--native-standard-air` | Use large air box (default small air, same as `particle_generation_em`) |
| `--native-no-air-relax` | Relax coil+plate only (no air particles) |
| `--reload-dir=...` | External `Reload.xml` (body names must be **CoilSourceBody**, **PlateBody**) |
| `--no-power-scaling` | Disable 50 kW field scaling (required for TEAM7 Bz benchmark) |
| `--compare-team7-bz` | A1–B1 air-probe coil-only Biot vs CSV phase0 |
| `--compare-team7-bz-loop` | Same + circular equivalent-coil diagnostic (R=mean_radius) |
| `--compare-team7-bz-rect-loop` | Same + rectangular STL footprint coil diagnostic (z=bottom/mid/top) |
| `--team7-coil-turns=N` | Override default 2742 turns (sensitivity) |
| `--coil-j-outer-shell[=0.85]` | Assign J on outer xy shell only (experimental; usually worse) |
| `--team7-reference-dir=...` | TEAM7 CSV directory |

**TEAM7 Bz benchmark (native reload + no-phi + no scaling):**

```bash
cd /home/yyc/SPHinXsysSYCL/build
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/test_3d_ophelie_team7 \
  --native-stl --reload=1 --no-phi --sigma=1e4 --no-power-scaling \
  --compare-team7-bz --compare-team7-bz-loop --compare-team7-bz-rect-loop \
  --team7-reference-dir=/path/to/tests/extra_source_and_tests/3d_examples/reference_data/team7 \
  --state_recording=0
```

Typical results today (dp=6 mm reload): full-line RMS ≈66%; point at x=126 mm ~15% error; sim peak at x≈198 (coil centroid side) vs ref peak at x=126. STL coil mesh x∈[94,294] mm; probes with x<94 are far-field. CSV may be **total B** (includes plate induction); `--no-phi` is **coil-only**.

Output: `build/output/team7_bz_A1B1_compare.csv`

If you use `Reload.xml` from `particle_generation_em` (body names `Coil`/`Plate`), rename bodies or re-run `--relax=1` in this program to write reload.

### Analytic TEAM7 smoke (1 m scaffold)

```bash
cd /home/yyc/SPHinXsysSYCL/build

# 1) No φ, fastest smoke (expect passed=1)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/test_3d_ophelie_team7 \
  --state_recording=0 --skip-relax --no-phi --sigma=1e4

# 2) With φ PCG (coarse-grid smoke: add --ophelie-smoke so phi_rel_res does not fail the run)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/test_3d_ophelie_team7 \
  --state_recording=0 --skip-relax --sigma=1e4 --ophelie-smoke

# 3) View fields: add state_recording=1; VTP in build/output/
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/test_3d_ophelie_team7 \
  --state_recording=1 --skip-relax --no-phi --sigma=1e4
```

**How to read your log:**

| Run | passed | Notes |
|------|--------|------|
| `--no-phi` | 1 | Main path OK: Biot, E/J, power scaled to 50 kW |
| default φ | 0 | `phi_rel_res` uses **L∞** (`max|rhs−lhs|/max|rhs|`); PCG log `rel_res_l2_vol` is reference only |
| φ after `--reload` | 0 | On coarse grid + PCG at 6000 steps, `rel_res_linf` may still ≫ tolerance; use `--ophelie-smoke` or `--phi-solver=GMRES` |
| `divJ_red=0.72` | — | On skip-relax coarse grid, φ **does not** improve divJ; TEAM7 smoke **warning only** |

`P_raw≈0.037`, `power_scale≈1.3e6`: on coarse grid + `sigma=1e4`, raw Joule heat is small; scaling to `target_P=50000` is expected (scheme B active: `field_scale≈1155`).

### Operator tests (run before formal TEAM7 acceptance)

```bash
cd /home/yyc/SPHinXsysSYCL/build
ninja test_3d_ophelie_biot_savart_direction test_3d_ophelie_coil_source_direction \
  test_3d_ophelie_level0_uniform_A test_3d_ophelie_scaling

./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_biot_savart_direction/bin/test_3d_ophelie_biot_savart_direction
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_coil_source_direction/bin/test_3d_ophelie_coil_source_direction
```

### Box regression (divJ / Level0 vs φ)

```bash
cd /home/yyc/SPHinXsysSYCL/build
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie/bin/test_3d_ophelie \
  --state_recording=0 --ophelie-compare-level0
```

Default full run:

```bash
cd /home/yyc/SPHinXsysSYCL/build
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/test_3d_ophelie_team7 --state_recording=0
```

Particle positions during relax are written to `build/output/` as `PlateBody_ite_0000000000.vtp`, `PlateBody_ite_0000000100.vtp`, … (default every **100** steps, same as `test_3d_particle_relaxation_single_resolution_sycl`). This works even with `--state_recording=0`. Use `--relax-vtp-every=25` for denser snapshots or `--relax-vtp-every=0` to turn off.

Final electromagnetic fields need `--state_recording=1` (one snapshot `*_ite_0000000000.vtp` at end of run).

Progress lines look like `[ophelie] relax:PlateBody | step 50/400 | elapsed_s=...`.

### Useful flags

| Flag | Effect |
|------|--------|
| `--relax-steps=200` | Fewer relaxation iterations (default 400) |
| `--relax-log-every=25` | Print relax progress every N steps (default 50) |
| `--relax-vtp-every=100` | Write relax particle VTP every N steps (default 100; `0` = off) |
| `--skip-relax` | Skip relax (fast debug, worse particle quality) |
| `--sigma=1e4` | Plate conductivity (default `1e4`; physical Al `3.54e7` needs finer mesh) |

SYCL builds use the same CK relax pipeline as `test_3d_particle_relaxation_single_resolution_sycl` (cell list → relation → residual → scaling → position → levelset bounding each step). Legacy `RelaxationStepInner` is not used on device.

| `--self-induction` | Enable self-induction loop |
| `--reload=1` | Load relaxed particles from `reload/Reload.xml` (skip in-code relax) |
| `--reload-dir=PATH` | Folder containing `Reload.xml` (default `./reload` under cwd) |

### Relax alignment with official SYCL example

Same core steps as `tests/tests_sycl/3d_examples/test_3d_particle_relaxation_single_resolution_sycl`:

| Step | Official SYCL | OPHELIE `relaxSolidBodyParticles` |
|------|-----------|-----------------------------------|
| Level set | `defineBodyLevelSetShape(par_ck).correctLevelSetSign()` | same (team7) |
| Randomize | `RandomizeParticlePositionCK` | same |
| Each step | cell list → relation → residual → scaling → position → bounding | same |
| Progress VTP | `BodyStatesRecordingToVtpCK` every 100 steps | `BodyStatesRecordingToVtpCK` (`--relax-vtp-every`) |
| Reload | one `writeToFile` for single body | **one** `ReloadParticleIO(sph_system)` writes Plate+Coil |

TEAM7 geometry is still a scaffold; you can replace `Reload.xml` for your own case—the relax pipeline is unchanged.

### Pre-relaxed particles (recommended for your geometry case)

1. Generate reload under **build** (cwd = test `bin/`):

```bash
cd ~/SPHinXsysSYCL/build
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/test_3d_ophelie_team7 \
  --relax=1 --state_recording=0 --relax-log-every=50
# → build/.../test_3d_ophelie_team7/bin/reload/Reload.xml (PlateBody AND CoilSourceBody)
```

2. Run EM on those particles (from the same `bin/` directory, or pass `--reload-dir=`):

```bash
cd build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin
./test_3d_ophelie_team7 --reload=1 --state_recording=0 --sigma=1e4
```

Do **not** put `Reload.xml` under the source case folder (`tests/.../test_3d_ophelie_team7/reload/`).

Body names in XML must be `PlateBody` and `CoilSourceBody`. Geometry parameters in `electromagnetic_ophelie_team7_geometry.h` must match your relax case.

## Defaults

| Item | Value |
|------|--------|
| `dp` | 0.04 m |
| `frequency` | 50 Hz (TEAM7-like) |
| `sigma` (plate) | 3.54e7 S/m (aluminum) |
| Plate | cylinder, r=0.36 m, half-thickness 0.02 m |
| Coil | annulus r_in=0.12, r_out=0.18, half-height 0.12 m |

Same OPHELIE CLI flags as `test_3d_ophelie` (`--phi-solver=`, `--self-induction`, etc.).

If `sigma=3.54e7` makes φ solve stiff on your mesh, try lowering `sigma_glass_` in code or ask for a CLI override.
