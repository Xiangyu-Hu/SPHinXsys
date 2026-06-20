# TEAM7 Source Smoke Archive

TEAM7 is **not** the current main validation target. Results below are archived for source-field smoke and future total-field benchmark.

## Geometry

- coil.stl: x=94–294 mm, y=0–200 mm, z=49–149 mm
- plate.stl: x=0–294 mm, y=0–294 mm, z=0–19 mm
- A1–B1 probe: (x, 72, 34) mm, x=0..288

## Key results (dp=6 mm reload, no-phi, no-power-scaling)

| Model | RMS rel err | peak_x_sim | peak_x_ref | Notes |
|-------|-------------|------------|------------|-------|
| volume-e_theta | 66.5% | 198 mm | 126 mm | x=126: sim/ref ≈ 8.95/7.81 mT (~15%) |
| circular loop (R=mean) | 56.0% | 198 mm | 126 mm | diagnostic |
| filament-racetrack z=149 | 44.1% | 198 mm | 126 mm | lowest RMS, amplitude low |
| filament-racetrack z=99 | 61.7% | 198 mm | 126 mm | |
| inset=20, z=49 | 275% | 126 mm | 126 mm | peak aligns but RMS unusable |

## Conclusions

1. Simple coil-only source path fitting (e_theta, racetrack inset/z sweep) does **not** stably reproduce TEAM7 Bz reference.
2. Reference likely includes **total B** (coil + plate-induced), not coil-only.
3. dp=3 vs dp=6 does not materially change RMS (~66%).
4. **Do not continue** racetrack/shell path fitting as mainline work.

## Future TEAM7 work (P2)

- probe total B = coil_B + plate_induced_B
- phase0 / phase90, 50 Hz / 200 Hz
- after self-induction complex channel is fixed

## Current mainline

See `OPHELIE_FRENCH_REDUCED_PRE_GEOMETRY_TASKS.md` and `test_3d_ophelie_french_reduced`.
