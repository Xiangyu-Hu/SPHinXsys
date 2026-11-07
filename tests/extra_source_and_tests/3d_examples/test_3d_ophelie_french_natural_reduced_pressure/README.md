# test_3d_ophelie_french_natural_reduced_pressure (Doc126 / Doc131)

French cold-crucible melt cylinder (R = 0.25 m, H = 0.185 m, bottom at z = 0) in a **closed flat-lid**
domain, WCSPH with the **reduced-pressure (p')** formulation and the production
`BoussinesqBuoyancyForceCK`. Doc131 adds the thermal completion (power contract, boundary-area audit,
conduction + Robin/radiation boundary, T0-T5 tests) **in this same binary**; R0-R3 keep the Doc126 flow path.

## Flow contract (unchanged from Doc126)

| Item | Setting |
| --- | --- |
| Momentum pressure | `PerturbationMomentum = 1` (acoustic uses p') |
| Background | `BackgroundPressure = rho0 g (H - z)`, z from bottom, rewritten every advection step |
| Wall density | `AcousticDensMode = 1`; default `ReferenceBodyAcceleration = 0`; optional `--cref-buoyancy` sets `a_ref := F_boussinesq/m` each advection step (cancels buoyancy wall dens jump; Doc134) |
| Gravity | **No `GravityForceCK`** (Gravity object only feeds `BoussinesqBuoyancyForceCK`) |
| Buoyancy | `BoussinesqBuoyancyForceCK` kept (beta, Tref) |
| Density | continuity-only (default) or `--density-summation` + `DensityRegularization<Internal>`; no FreeSurface, no dens-sum cap, no analytic density rewrite per step |
| Acoustic walls | `OphelieFrenchClosedCrucibleWallShape` (side + bottom + lid) |
| Viscous walls | `OphelieFrenchNaturalCrucibleWallShape` (open cup, rim = 0): side + bottom no-slip |
| Lid | free-slip: no viscous contact on the lid (`top_bc = flat_free_slip_viscous_omit`) |
| Wall thickness | `4 dp` (`--wall-thickness-factor`), mesh resolution 20 |

## Modes (`--mode=`)

Flow mechanism (T frozen, no Q, no conduction; short-time mechanism smoke, not accuracy validation):

- **R0** - `T = Tref`, quiet hold. Pass: `U_rms < 5e-4` at end.
- **R1** - `T = Tref + a z` (`--a-strat`, default 50 K/m). Density initialised once from
  `p'_strat = 1/2 rho0 g beta a (z^2 - H^2)` (Mass synced once); no per-step rewrite. Pass: `U_rms < 1e-3`.
- **R2** (Doc131 D2: gate is now sign-correct, see "Doc131 D additions"; the pass text on the next lines is the OLD gate) - frozen radial driver `T = Tref + dT (r/R - <r/R>)` (`--dT-radial`, default 10 K; hot wall); mean removed
  by default (`--r2-keep-mean` keeps raw `dT r/R`). Pass: `U_rms > 1e-4` and (circulation proxy nonzero or buoyancy
  present). Proxies: `circulation_uz_proxy` (>0 upflow at wall), `circulation_ur_proxy` (<0 inflow aloft).
- **R3** - quiet R0 hold, then a divergence-free impermeable axisymmetric meridional roll of amplitude
  `--perturb-u` (default 1e-3 m/s). Pass: perturbed `U_rms > 0`, `U_rms(end)/U_rms(perturb) < 5`.

Thermal (Doc131):

- **baseline** - **`thermal_smoke_baseline = 1`** (was `literature_baseline = 1`; `literature_baseline` is now always `0`,
  `calibrated_match = 0`). Fixed discrete Q (`--q-power`, default 50 kW) + production conduction + Robin/radiation BC
  (side/bottom `h` with `t_cool`, top `h_top` + radiation with `t_ambient`) + thawed T (`--thaw-T`, default on) +
  Boussinesq buoyancy. Pass: finite, no runaway, `P_applied = P_target`, energy ledger `error_input <= 1%`.
  `--no-thaw-T` freezes T (flow diagnostics only). `--legacy-rad-clamp` restores the OLD host path (top-layer-only radiation
  with `max(T, T_amb)` clamp, no conduction, no side/bottom loss) for comparison only.
- **T0** no Q, adiabatic, uniform T, frozen flow: T stays constant (`|dT|max <= 1e-6 K`).
- **T1** 50 kW discrete uniform Q, no loss, frozen flow: `DeltaH = P_target * t` and `DeltaT = P t/(M cp)` (rel <= 1e-4).
- **T2** no Q, adiabatic, non-uniform start (`--a-strat` z-gradient + `--dT-radial` radial part), frozen flow: total energy
  conserved (`|dH| <= 1e-5 sum m cp |dT0|`) and temperature variance / spread decrease. Default end time 600 s.
- **T3** isothermal start, **one face enabled** (`--t3-face=top|side|bottom`), Q = 0, frozen flow: field loss `-DeltaH`
  vs the independent analytic `flux(T0) * A_geom * t` (top uses convection + radiation). Tolerance 2 % with
  `--area-normalize`, 15 % without (the raw SPH shell areas are 12-19 % below the geometric areas at dp = 0.025, so the
  un-normalised run is expected to fail and shows the area gap). Default end time 0.2 s, thermal dt 0.01 s.
- **T4** 1D cosine decay `dT = A cos(pi z/H)` (`--t4-amp`, default 10 K), all faces adiabatic, frozen flow: LS-fitted
  amplitude vs `exp(-alpha (pi/H)^2 t)`, decay-rate error <= 10 % and energy conserved (engineering gate). Default
  end time 1000 s, thermal dt 5 s.
- **T5** full BC + fixed discrete Q, **buoyancy off** (`beta = 0`), flow runs: energy ledger `error_input <= 1%` and
  `P_applied = P_target`.

`--freeze-flow` is forced for T0-T4 (no viscous/acoustic/buoyancy update, positions/velocity/density untouched, settle
phase skipped); it may be added to other modes. For frozen runs `--thermal-dt` is checked against the explicit stability
limit (exit 1 if exceeded).

## Properties and CLI (Doc131 A3)

| Option | Default | Note |
| --- | --- | --- |
| `--cp` | **1150** J/(kg K) | CEP contract value (was 1300; 1300 is historical only) |
| `--k` | 4 W/(m K) | registered reference |
| `--h-side` / `--h-bottom` / `--h-top` | 300 / 35 / 20 W/(m^2 K) | **CODE_PROVISIONAL** (not re-certified) |
| `--t-cool` | 300 K | side/bottom sink, **CODE_PROVISIONAL**; independent of `--t-ambient` |
| `--emissivity`, `--t-ambient` | 0.8, 300 K | top radiation + top convection sink |
| `--rho0 --mu --beta --t-ref --c0 --g` | 2750, 4, 3e-5, 1473, 20, 9.81 | |
| `--power-mode=` | `fixed_discrete_absorbed_power` | or `fixed_volumetric_density` (diagnostic only) |
| `--area-normalize` | off | normalise B3 face areas to pi R^2 / 2 pi R H |
| `--freeze-flow`, `--t3-face=`, `--thermal-dt=`, `--t4-amp=` | | thermal tests |
| `--output-dir=` | none | creates the dir; relative `--summary-json` / `--out-csv` are prefixed (defaults inside it: `summary.json`, `doc131_<mode>.csv`) |
| `--legacy-rad-clamp`, `--legacy-thermal-sort` | off | legacy/diagnostic switches (see below) |

Other options: `--mode= --dp=0.025 --settle-time= --end-time= --u-ref=0.5 --a-strat= --dT-radial= --perturb-u=
--kernel-correction --density-summation --r2-keep-mean --glass-radius=0.25 --glass-height=0.185 --q-power=50000
--thaw-T/--no-thaw-T --no-top-loss --wall-thickness-factor --mesh-resolution --summary-json= --out-csv=` (`--help`).

Parsing is strict: empty / non-numeric / trailing-garbage / NaN / Inf values are rejected; `dp, c0, cp, k, rho0, R, H > 0`,
times `>= 0`, `emissivity in [0,1]`; unknown arguments exit with 1. The effective configuration is printed at start
(`[doc131][config]`) and the main values are in the JSON.

## Power contract (Doc131 B)

`fixed_discrete_absorbed_power` (default): `V_i^Q = m_i/rho0`, `N_Q = sum(s_i V_i^Q)`, `q_i = P_target s_i / N_Q`,
`DeltaT_i = q_i dt/(rho0 cp)` (so `P_applied = sum q_i V_i^Q = P_target`). `fixed_volumetric_density` uses
`q = P/V_geom` (diagnostic; reproduces the old 47.57 kW gap: `sum V_Q / V_geom = 0.9515`). `N_Q <= 0` is refused.
Reported: `P_target_W, P_applied_W, P_gap_rel, V_geom_m3, sum_V_Q_m3, mass_sum_kg, sum_V_Q_over_V_geom, power_norm_coeff,
q_min/q_max`.

## Boundary area audit (Doc131 B3)

Shell thickness `dp`; `A_i = m_i/(rho0 dp)`. Top `z > z_top - dp` (target pi R^2), bottom `z < z_bottom + dp` (pi R^2),
side `r > R - dp` (2 pi R H). Corner particles are assigned to **one** primary face (closest surface), no double count.
With `--area-normalize` each face weight is scaled so `sum A` equals the target. JSON: `A_{top,side,bottom}_{raw,target,
rel_err_raw,after,rel_err,scale,ai_min,ai_max}`, `n_{top,side,bottom}`, `area_n_corner_primary_assigned`. The weights are
written to the thermal boundary weight fields consumed by the production Robin/radiation kernel (so T3/T5/baseline use
exactly the audited areas).

## Conduction / boundary path (Doc131 C)

`conduction_path = production_ck:OpheliePairwiseLaplaceCK + ApplyLaplaceDiffusionDtCK`, boundary
`ApplyOphelieThermalFrenchNaturalBcCK` with `OphelieFrenchBoundaryLossPowerReduceCK` as the loss ledger (same kernels as the
natural/stirring cases, run on the fluid body with the fluid `Inner<>` relation). Q uses a small local per-particle kernel
(dt at exec). `OphelieThermalDeltaT` is the authoritative thermal variable (`T = t_ref + DeltaT`); the energy ledger uses it.
Order per thermal step: Laplace -> Q -> diffusion -> boundary loss ledger -> BC. Heat flows both ways (no `max(T, Tamb)`
clamp); radiation only on the top face; top sink `t_ambient`, side/bottom sink `t_cool`; shell thickness `dp`.
Energy ledger (JSON): `E_source_J, E_boundary_J, DeltaH_J, energy_residual_J, energy_error_input / _exchange / _net`.

**Particle sort**: only `evolving` variables are permuted by `ParticleSortCK`. The thermal fields (Temperature, DeltaT,
Laplace, conductivity, boundary mask/face/weights, Q) are now registered as evolving so they follow the particles.
Before Doc131 Temperature was not permuted, which scrambled T against positions after each sort (R1/R2 numbers change
accordingly, e.g. R1 `U_rms` 5.77e-5 -> 1.75e-5, R2 1.87e-4 -> 2.27e-4; `--legacy-thermal-sort` reproduces the old numbers).

## Output (JSON `--summary-json`, CSV `--out-csv`)

- Legacy fields: `case, mode, dp, R, H, perturbation_momentum, acoustic_dens_mode, gravity_force_ck=0, boussinesq=1,
  U_rms, Uz_rms, Ur_rms, Utheta_rms, K, U_mean_*, T_*, passed, note, top_bc`, mode diagnostics.
- **A2 EOS semantics**: `eos_initialization_ok` (init residual), `aux_density_range_ok` (|rho - rho0| < 10 % at end),
  `eos_field_consistency_ok` (Pressure vs EOS, sampled read-only right after the acoustic 1st half of the last production
  substep, where Pressure is current; `null` + reason for frozen runs), `pressure_storage_semantics`
  (**actual**: `Pressure = EOS(rho_half) + BackgroundPressure`, i.e. total pressure; `p' = Pressure - BackgroundPressure`),
  `pressure_sample_phase`. Legacy `eos_identity_ok = eos_initialization_ok AND aux_density_range_ok` is kept for old scripts and
  is **not** field consistency. `pi_mean, pi_rms_demean, pi_min, pi_max` come from the Density EOS reconstruction at production
  end; `pressure_field_pi_*` come from `Pressure - BackgroundPressure` at the sample phase (plus `eos_field_residual_*`).
- **A4 reproducibility**: `doc131 = 1, real_bytes, precision_note, argv, effective_cp, effective_k, power_mode,
  diagnostic_consistent`. Production end state (U, K, T, t, pi, dH, mass) is read before the diagnostic force evaluation.
- **Naming**: `thermal_smoke_baseline` (baseline mode only), `literature_baseline = 0`, `calibrated_match = 0`.
- **Doc131 gates** (`value` is `1`, `0` or `null` + `<name>_reason`; never defaulted to true): `source_power_verified`
  (Q modes), `boundary_area_verified` (2 % after normalisation), `conduction_checked` (T2/T4), `thermal_boundary_checked`
  (T3 per face); OPEN (null): `weak_buoyancy_response_checked, precision_sensitivity_checked, slip_boundary_checked`;
  `full_thermal_flow_baseline = 0` and `literature_targets_available = 0` (not implemented / no data).
- CSV: the flow columns plus `dH_J, E_source_J, E_boundary_J`.

## Doc131 D additions (weak buoyancy / R2 sign matrix)

- `--r2-reverse` (negate R2 dT; same as `--dT-radial=-10`), `--r2-fixed-spatial-T` (diagnostic: T(r) re-applied from the
  current position every advection step; default = Lagrangian carry, stated in JSON `r2_T_driver_mode`),
  `--acoustic-cfl=0.6`, `--advection-cfl=0.25` (library defaults).
- R2 pass gate (replaces the old "circulation nonzero OR buoyancy present"): driven (`beta>0`, `dT!=0`) requires
  `U_rms>1e-4`, `sign(dT)*circulation_uz_proxy>0` and `sign(dT)*circulation_ur_proxy<0`; `beta=0`/`dT=0` requires
  `U_rms<1e-5`. JSON: `r2_gate_kind, r2_sign_uz_ok, r2_sign_ur_ok, r2_dT_sign`. Cross-run reversal / beta=0 comparison:
  `docs/ophelie/reference/doc131_runs/analyze_doc131_D.py`.
- CSV gains `pi_mean, pi_rms_demean, pi_min, pi_max, rho_min, rho_max, circ_uz_proxy, circ_ur_proxy, U_max, top_uz_rms,
  top_uz_maxabs, n_outside_strict, n_outside_half_dp, z_max` (appended after the old columns).
- JSON `d4_*`: top/side/bottom band velocity and escape counts (behaviour report only). Gates:
  `precision_sensitivity_checked = null` (`build_is_float_only_no_double_binary`), `slip_boundary_checked = null`
  (`no_independent_plate_slip_mms_yet`).
- Results: `docs/ophelie/132_DOC131_A_TO_D_PROGRESS_2026-10-14.md`.

## Doc133/134/136 additions

- `--cref-buoyancy` (candidate): `ReferenceBodyAcceleration := F_boussinesq/m` each advection step (buoyancy-aware wall
  dens). Default remains off. `--probe-cref-buoyancy` aliases it and enables `--drift-probe`.
- `--density-increment-audit` (Doc136 P0): half-step Density write budget vs store; OriginalID samples;
  `--density-audit-stride=` / `--density-audit-max-rows=`.
- `--q-shape=uniform|radial_skin` and `--q-from-rz-csv=path` (`em_rz_csv`: prior EM Q(r,z) CSV with `r_m,z_m,Q_Wpm3`,
  bilinear sample, discrete-renormalise to `P_target`; `em_q_live_solve=0`).
- `--state-record-interval=` VTP of fluid fields to `./output`.
- Progress: `docs/ophelie/133_*.md` … `137_*.md`; runs under `reference/doc136_runs/`.

## OPEN (not finished for Doc131 A–F)

D1 contamination budget still narrowly FAIL with cref (U_R1/U_R2 ≈ 10.4% > 10%); D3 float/double; D4 plate slip MMS;
live EM solve on this lattice; CEP digitisation / F calibration. `h_*` / `t_cool` remain CODE_PROVISIONAL.

## Run

```
./test_3d_ophelie_french_natural_reduced_pressure --mode=R0 --dp=0.025 --summary-json=r0.json
./test_3d_ophelie_french_natural_reduced_pressure --mode=T1 --output-dir=out_T1
./test_3d_ophelie_french_natural_reduced_pressure --mode=T3 --t3-face=side --area-normalize --output-dir=out_T3_side
./test_3d_ophelie_french_natural_reduced_pressure --mode=baseline --area-normalize --output-dir=out_baseline
```
