# test_3d_ophelie_team7_complex_edge_flux

TEAM7 validation with **volume-racetrack** coil source (Laplace-aligned) and **complex edge-flux** on the plate.

**Development log**: see [TEAM7_VALIDATION_LOG.md](TEAM7_VALIDATION_LOG.md) (milestones + open items); each run appends to `output/team7_validation_history.txt`.

## Prerequisites

1. Generate reload (coil + plate only, no air), e.g. 3 mm:

```bash
cd build/.../particle_generation_TEAM7/bin
./particle_generation_TEAM7 --team7-case=uniform-small-dp3 --team7-no-air
# → reload_cases/uniform_small_dp3/  (CoilSourceBody + PlateBody in Reload.xml or *_rld.xml)
```

2. Copy or symlink `Reload.xml` into this test's **build** `bin/reload/`, or pass `--reload-dir=` pointing at another build reload folder.

Body names: `CoilSourceBody`, `PlateBody`.

## Levels

| `--team7-level=` | Description |
|------------------|-------------|
| `coil-only` (L1) | Coil Biot Bz vs muFEM reference; no phi / no A_ind |
| `one-way` (L2) | Complex edge-flux + A_ind one-way; probes decomposed as B_coil + B_ind |
| `picard` (L3) | Self-induction Picard + B_total probes |

### Pass semantics (2026-06-02)

| Flag | Meaning |
|------|---------|
| `smoke_passed` | Case completes; coil Bz gate (L2); EM finite; exit code |
| `team7_phase0_source_passed` | Coil phase0 RMS within threshold |
| `team7_phase90_probe_passed` | **Diagnostic only** (default threshold 0.5 mT abs RMS); not hard gate |
| `diagnostic_only` | `source_scale≠1`, `j_post_scale`, `imag_a_sign≠1`, or L2 one-way |
| `team7_validation_passed` | Strict Picard + physical source only; **always 0 for L2** today |

### L2 gates vs logged metrics

| Metric | Gated? | Notes |
|--------|--------|-------|
| `bz_rms_coil` | **yes** (smoke, < 0.70) | Coil-only Biot on A1–B1 line |
| `coil_path_audit` | **yes** | NI vs integrated J_src |
| EM finite (`max_J`, φ res) | **yes** | one-way only |
| `bz_rms_total` (phase0) | log only | Equals coil when J_real≈0 |
| `bz_rms_total_phase90` | log only | Induction in imag channel |
| `bz_rms_total_mag` | log only | √(phase0²+phase90²) |
| `bz_rms_ind_imag_skin` | log only | Biot(J) from plate-top skin (z ≥ z_top−2dp) |
| imag-chain audit | log only | `\|J\|`, `\|E\|`, `P_vol`, `cos(J,E)` before probes |
| `--team7-ind-j-post-scale=auto` | off | Diagnostic J scale → recompute A_ind; improves phase90 RMS ~12→~1.7 |
| x-band phase90 CSV | log only | Per-x-span RMS (`team7_probe_phase90_xbands_*.csv`) |
| plate depth profile CSV | log only | J/B vs z (`team7_plate_depth_profile_*.csv`) |

## Example commands

`--reload-dir=` must be a **real path** (do not paste `...` placeholders).

From `build/` (reload auto-found if you omit `--reload-dir` after particle generation):

```bash
# L2 — one-way edge-flux (default level, coil source_scale=0.754)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7_complex_edge_flux/bin/test_3d_ophelie_team7_complex_edge_flux \
  --reload=1 \
  --reload-dir=$HOME/SPHinXsysSYCL/build/tests/extra_source_and_tests/3d_examples/particle_generation_TEAM7/bin/reload \
  --team7-level=one-way --native-dp-mm=3 --ophelie-edge-flux-complex=1 \
  --team7-reference-dir=$HOME/SPHinXsysSYCL/tests/extra_source_and_tests/3d_examples/reference_data/team7
```

From test `bin/` (default search: `./reload`, then `../reload`, then `build/reload`):

```bash
cd build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7_complex_edge_flux/bin
./test_3d_ophelie_team7_complex_edge_flux --reload=1 \
  --team7-level=one-way --native-dp-mm=3 --ophelie-edge-flux-complex=1
```

```bash
# L1 — coil-only Bz
./test_3d_ophelie_team7_complex_edge_flux --reload=1 --team7-level=coil-only --native-dp-mm=3

# L3 — Picard
./test_3d_ophelie_team7_complex_edge_flux --reload=1 --team7-level=picard --native-dp-mm=3 --self-induction-max-iter=8

# Physical NI (no coil Bz calibration): source_scale=1
./test_3d_ophelie_team7_complex_edge_flux --reload=1 --team7-coil-source-scale=1 ...

# Optional: feedback re-solve with A_total (experimental)
./test_3d_ophelie_team7_complex_edge_flux --reload=1 --ophelie-aind-one-way-feedback=1 ...
```

## CLI reference

| Option | Default | Description |
|--------|---------|-------------|
| `--team7-coil-source-scale=` | `0.754` | Calibrate volume-racetrack peak Bz vs TEAM7 filament reference |
| `--team7-reference-dir=` | auto | Directory with `TEAM7_Bz_A1_B1_reference_mT.csv` |
| `--ophelie-aind-one-way-feedback=` | `0` | L2: skip A_total second edge-flux solve by default |
| `--team7-ind-j-post-scale=` | off | `auto` or factor: scale J/E and refresh A_ind (diagnostic) |
| `--team7-phase90-convention=` | `minus-imag` | Probe mapping: `minus-imag` → `B_phase90=−B_imag`; `plus-imag` → `+B_imag` |
| `--team7-probe-phase90-neg-ind=1` | (alias) | Same as `--team7-phase90-convention=minus-imag` |
| `--team7-validation-mode=` | `smoke` | `smoke` allows legacy `source_scale=0.754`; `strict` requires `source_scale=1` |
| `--team7-output-tag=` | auto | Suffix for all CSV/log (`team7_probe_one-way_<tag>.csv`, …) |
| `--ophelie-edge-flux-imag-a-sign=` | `1` | **Debug-only** solver sign; use `phase90-convention` for TEAM7 probes |
| `--team7-bz-rms-threshold=` | level-specific | Override Bz RMS smoke thresholds |

## Outputs under `./output/`

All artifacts include `<tag>` (auto-generated or `--team7-output-tag=`):

- `team7_probe_<level>_<tag>.csv` — probe decomposition
- `team7_bz_<level>_<tag>_coil_only.csv` / `_total.csv` — phase0 compare
- `team7_summary_<level>_<tag>.txt` — pass report snapshot
- `team7_validation_history_<level>_<tag>.txt` — per-run record
- `team7_validation_history.txt` — **also appended** (aggregate)
- `team7_probe_phase90_xbands_<level>.csv` — phase90 RMS by x-band
- `team7_plate_depth_profile_<level>.csv` — vol-weighted J/B per z-bin
- `PlateBody_ite_0000000000.vtp` — plate fields for ParaView (one-way / picard)

Reference CSV copied to `bin/reference_data/team7/` at configure time.
