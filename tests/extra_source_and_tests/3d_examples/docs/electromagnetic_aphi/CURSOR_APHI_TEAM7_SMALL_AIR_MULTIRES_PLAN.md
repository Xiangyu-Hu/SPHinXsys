# TEAM7 small air — uniform vs multiresolution validation

## Goal

On **small** air preset, verify **air multiresolution** (fine near coil/plate, coarse far field) gives **similar Bz** to **uniform dp=6 mm**, then use the same MR recipe on **big** air to save particles.

## Air presets (config header)

`tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_team7_native_geometry_config.h`

| Preset | Box [m] | Notes |
|--------|---------|--------|
| small | (−0.05…0.35) × 0.25 h | dev / MR check (~173k air @ uniform 6 mm) |
| legacy | (−0.20…0.50) × 0.50 h | medium (~1.12M air @ uniform 6 mm) |
| big | (−1.353…1.647) × 0.749 h | COMSOL far field |

## MR discretization (actual code behaviour)

**There is no `--team7-dp-air-finest` CLI.** Spacing is controlled only by:

| Symbol | Formula | Default (dp_reference=0.006 m) |
|--------|---------|--------------------------------|
| `dp_reference` | finest spacing (coil, plate, air near conductors) | 6 mm |
| `local_refinement_levels` (air only) | dyadic coarsening levels | 0 = uniform |
| coarsest air spacing | `dp_reference × 2^levels` | L1 → 12 mm; L2 → 24 mm |

Coil/plate always stay uniform at `dp_reference` (`solid_refinement_levels=0`).

`SmoothingLengthRatio` is written at lattice generation **and** updated each relax step via `UpdateSmoothingLengthRatio` in `particle_generation_em.cpp`.

## Canonical reload cases (`./reload_cases/`)

| Case id | CLI preset | Air | levels | coarsest air |
|---------|------------|-----|--------|--------------|
| `uniform_legacy_dp6` | `--team7-case=uniform-legacy-dp6` | legacy | 0 | 6 mm |
| `multires_small_dp6_l2` | `--team7-case=multires-small-dp6-l2` | small | 2 | **24 mm** |

Probe/solve after staging:

```bash
--team7-reload-case=uniform_legacy_dp6 --team7-write-vtp
```

## Step A — uniform baseline (small air + dp=6 mm)

```bash
cd /home/yongchuan/sphinxsys/build/tests/extra_source_and_tests/3d_examples/particle_generation_em/bin
./particle_generation_em --team7-dp=0.006 --team7-air=small --team7-discretization=uniform
```

(`--team7-case=uniform-legacy-dp6` uses **legacy** air, not small — use explicit flags above for small-air uniform.)

## Step B — multires (small air, levels=1, coarsest 12 mm)

**Validated recipe** (not the `multires-small-dp6-l2` preset, which uses levels=2 / 24 mm coarsest):

```bash
/home/yongchuan/sphinxsys/build/tests/extra_source_and_tests/3d_examples/particle_generation_em/bin/particle_generation_em \
  --team7-dp=0.006 --team7-air=small \
  --team7-local-refinement-level=1
```

Check metadata in `./reload/team7_native_geometry.txt`:

```text
air_local_refinement_levels,1
dp_air_coarsest_m,0.012        # = 0.006 × 2^1
dp_0_m / dp_reference_m,0.006  # finest (coil/plate + air near conductors)
air_multiresolution,1
```

## Step C — P3 Bz compare

```bash
cd /home/yongchuan/sphinxsys/build
cmake ..

cd /home/yongchuan/sphinxsys/build/tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_team7_native_geometry_bz_reference_probe/bin
./test_3d_aphi_ck_team7_native_geometry_bz_reference_probe \
  --team7-reload-case=<case_id> --team7-write-vtp
```

**Pass criterion (engineering):**

- `A1_B1_profile_rel_err_real` within ~5% absolute of uniform run on the **same** air preset
- Peak `A1_B1_max_abs_Bz_mT` within ~10% of uniform and ref (~7.8 mT)

## Step D — big air + MR (after small MR matches)

```bash
./particle_generation_em \
  --team7-dp=0.006 --team7-air=big \
  --team7-local-refinement-level=3
```

`levels=3` → coarsest 0.006×2³ = **48 mm** far from coil (fewer particles than uniform 6 mm filling entire big box).

## Bz comparison results

### Small air (2026-06-04)

| Case | Air particles | A1-B1 L2 real | Peak Bz | P3 |
|------|---------------|---------------|---------|-----|
| uniform 6 mm | 173133 | 20.4% | 7.89 mT | pass |
| MR levels=1 (coarsest 12 mm) | 233004 | 24.1% | 7.90 mT | pass |

### Legacy air uniform (2026-06-05, default P3 reload)

| Case | Air particles | A1-B1 L2 real | Peak Bz | P3 |
|------|---------------|---------------|---------|-----|
| uniform 6 mm (`uniform_legacy_dp6`) | 1122763 | **13.1%** | 7.11 mT | pass |

### Known failure (open)

| Case | Result |
|------|--------|
| legacy air + `multires-small-dp6-l2` (levels=2, coarsest 24 mm) | GMRES `nan`, Bz=0, `solver_ok=0` — **not validated** |

## Implementation notes

Reload/solve: `AdaptiveBody` + `AphiTeam7NativeAdaptiveNearInnerSurface` + framework `CellLinkedList<AdaptiveSmoothingLength>` and CK `UpdateRelation` (same pattern as `test_2d_free_stream_around_cylinder_mr_sycl`).

Snapshots under `particle_generation_em/bin/reload_cases/` and P3 `bin/reload_cases/`.

## Naming debt (fix after GPT review)

- Preset `multires-small-dp6-l2` implies “finest 6 mm / level 2” but does **not** mean 3–6 mm band; coarsest is 24 mm.
- Consider renaming preset or adding `multires-small-dp6-l1` as the validated small-air MR case.
