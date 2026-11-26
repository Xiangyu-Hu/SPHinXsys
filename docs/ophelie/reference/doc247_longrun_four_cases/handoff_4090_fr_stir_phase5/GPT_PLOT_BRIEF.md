# Brief for plotting assistant

Compare two 10 s French-stirring SPH cases (same IC/Q/EM; only rheology differs).

## Files
- `01_MU_const_rpm10_.../french_T_rz.csv` vs `02_MU_varT_rpm10_.../french_T_rz.csv`
- `derived_midplane_velocity_rz/*_U_rz_midplane.csv` for streamlines
- time series in each case's `french_energy_budget.csv` (has mu_min/mean/max) and `french_stirring_em_monitor.csv`
- Fig.8a digitization overlay: `00_prototype_literature_seeds/STIR_Fig8a_contour_points.csv`

## Suggested figures
1. Side-by-side midplane T contours (const vs varT) at t≈10 s, optional Fig.8a overlay.
2. Side-by-side midplane Ur–Uz streamlines.
3. Time history: mu_mean, U_max, U_rms, T_mean for both cases.

## Caveats to put on figures
ASSUMED / engineering demo; varT had U_max≫c0 for most of the run; not literature PASS.
