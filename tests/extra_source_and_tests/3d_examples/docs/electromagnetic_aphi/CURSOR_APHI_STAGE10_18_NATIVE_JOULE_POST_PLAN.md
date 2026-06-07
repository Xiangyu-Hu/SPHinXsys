# Stage 10.18 — Native TEAM7 Joule post (entry)

> Parent: Stage 10.17 native SI TEAM7 Bz reference **CLOSED**  
> Record: `CURSOR_APHI_STAGE10_17_TEAM7_NATIVE_RECORD.md`

## Goal

On the **same SI native reload** as 10.17, verify eddy-current **Joule** post-processing:

- plate (conductor) Joule integral > 0  
- air Joule / plate Joule below diagnostic ratio (not physical heating target)  
- coil Joule reported but not used as heat target  

## Test

`test_3d_aphi_ck_team7_native_joule_post_smoke` — 50 Hz, TEAM7 coil-path RHS, writes `team7_native_joule_post_report/summary.csv`.

## Not in scope yet

- One-way thermal coupling / temperature solve  
- Jy probe reference comparison  
- Multiresolution / A-penalty production tuning  

## Smoke

Included in `run_team7_native_si_smoke.sh` (with P0, P1, source_driven, P3).
