# Thermo–EM coupling after the ChatGPT task (2026-09-14)

> Branch: `feature/EM-Biot-Savart`  
> Period: 2026-09-14 (same calendar day as the ChatGPT implementation brief).  
> Companion: [`10_THERMO_EM_COUPLING_ROUND1.md`](10_THERMO_EM_COUPLING_ROUND1.md) (code audit + CLI), [`09_FRENCH_LITERATURE_REPRODUCTION_STATUS_2026-09.md`](09_FRENCH_LITERATURE_REPRODUCTION_STATUS_2026-09.md) (pre-task snapshot).  
> Midterm claim remains an **auditable SPH–EM–thermal–flow prototype**, not Jacoutot paper-grade reproduction.

This note is the close-out of everything done **after** the user pasted the ChatGPT brief (“periodic thermo–EM, natural first, then 10 rpm stirring”). It is meant to be pasted back to ChatGPT / an advisor without the chat log.

---

## 0. One-paragraph verdict

The ChatGPT architecture was implemented and **run**, not only sketched. A shared periodic controller now does `T → σ(T) → azimuthal σ(r,z) → re-solve EM → Euler Q` on physical time. Natural convection (CEP 50 kW / 282 kHz, `dp = 15 mm`) is contained and couples for 50 s (`U_rms` stays ~2 cm/s from t ≈ 0; **no literature buoyancy cell**). Mechanical stirring (CES 60 kW / 300 kHz / 10 rpm, `dp = 5 mm`) completed **5 revolutions** with the same controller: absorbed power fell 60.00 → 54.81 kW at locked coil current because σ dropped with cooling. Heat loss still exceeds Joule input (~2×) because there is **no skull**. SSH VTP download is impractical; spatial scalars (`U_rms`, `T_out`/`T_ctr`, …) were added so ParaView is not required for this prototype gate.

**Gate B-light (Fig. 4-like dataflow in SPH) is demonstrated. Gate C (figures, skull, Np/Nu, 1000 s thermal steady state) is not.**

---

## 1. What the ChatGPT brief asked, and what we accepted

### 1.1 Accepted north star

Starting point of the brief: the code already had Biot–Savart `A_coil`, complex edge-flux φ, E/J/Q reconstruction, `P_recon` current calibration, and **one-shot** EM → Euler Q → SPH (plus 10 rpm kinematics). That is frozen-Q / one-shot, not Jacoutot Fig. 4.

The brief’s next milestone (quoted in spirit):

> On the French geometry and properties, drive **periodic** thermo–EM coupling from temperature-dependent conductivity; **natural convection first**, then the **same** controller on 10 rpm stirring.

Ordered work list from the brief, with status at the end of this day:

| # | Brief item | Status |
|---|------------|--------|
| 1 | Unify French case parameters | **Done** (`FrenchLiteratureParameters`) |
| 2 | Periodic thermo–EM coupling | **Done** (`shouldUpdateEM` on physical time) |
| 3 | `T → σ(T) → EM → Q` | **Done** (azimuthal average + under-relax) |
| 4 | Q map power book (no silent rescale) | **Done**; `--q-conservative-remap` is opt-in and logged |
| 5 | SPH internal conduction | **Done on periodic path** (constant k, cp, μ @ 1473 K) |
| 6 | Validate natural convection | **Path + hydro OK; not a literature cell** |
| 7 | Attach 10 rpm stirring | **Done**; 1 s / 1 rev / 5 rev periodic runs |
| 8 | Later: skull, `A_metal`, segmented crucible | **Still deferred** (as requested) |

### 1.2 Explicit non-goals (kept)

- No new pile of standalone tests  
- No segmented cold crucible, no `A_metal`  
- Do not mix RH200 client CAD into the French baseline (Euler CIC **numerics** may be reused)  
- Do **not** advertise `--balance-heat-loss` as paper physics  
- No large SPH rewrite  
- TEAM7 remains a side-track, **not** a French delivery gate  
- Full `k(T)/cp(T)/μ(T)` and a blocking 2–5% energy residual were **cut from round 1** (brief §8 vs §13 conflict)

### 1.3 Corrections applied to the brief

| Brief / ChatGPT item | Decision |
|----------------------|----------|
| Stirring “maybe y-up” | **Wrong for current meshes.** `glass-z.stl` / paddle `_z` are **+z**. Averager uses `FrenchVerticalAxis` (default Z). |
| Silent Q renormalization on stirring | **Removed** as default. |
| Periodic A_ind | Default **off** unless `--aind=on`. Frozen-Q / one-shot stirring still default A_ind on. |
| Natural cup as a lid | **Rejected.** Open free surface; short **rim** only (splash guard). Stirring walls stay `rim_height = 0`. |
| Generator 400 kW | Label only. Glass absorbed targets remain **50 kW** (natural) and **60 kW** (stirring). |

---

## 2. Code that landed

### 2.1 Shared

| File | Role |
|------|------|
| `tests/.../electromagnetic_ophelie_french_literature_parameters.h` | Two regimes (natural 50 kW / 282 kHz vs stirring 60 kW / 300 kHz); `EMControlMode`; `ThermoEMCouplingMode`; OPEN_PARAMETER notes |
| `tests/.../electromagnetic_ophelie_thermo_em_coupling.h` | r–z σ average, empty-bin fill, σ under-relax, `shouldUpdateEM`, power CSV, melt **spatial stats** |

Periodic path: SPH conduction **on** unless `--no-thermal-diffusion`. Energy CSV writes a residual; it is **not** a pass gate.

### 2.2 Cases

| Executable | Change |
|------------|--------|
| `test_3d_ophelie_french_natural_convection_frozen_q` | `--thermo-em-coupling=periodic`; default remains frozen-Q for regression; hydro rim + hydrostatic settle; spatial stats on the existing monitor cadence |
| `test_3d_ophelie_french_stirring_em` | Same controller; default still one-shot EM + 10 rpm + A_ind on; periodic: A_ind default off; `em_ok` = finite φ residual and ≥1 update (not a 1% power gate) |

### 2.3 Reloads (do not mix)

| Case | Reload | `dp` | Typical `n_glass` |
|------|--------|------|-------------------|
| Natural | `test_3d_ophelie_french_natural_glass_relax/bin/reload` | **15 mm** (this case’s production spacing) | ~10 474 |
| Stirring | `test_3d_ophelie_french_stirring_glass_relax/bin/reload` | **5 mm** (paddle blades ~8 mm; do not coarsen) | ~322 281 |

Natural `dp = 15 mm` is **not** a stirring-style production mesh. It is the reload that exists for the CEP cup. Stirring cannot use that XML.

### 2.4 Spatial diagnostics (SSH / no ParaView)

Glass VTP frames are tens of MB; remote copy is too slow. Added host-side melt stats instead of asking for ParaView:

| Quantity | Definition |
|----------|------------|
| `U_rms`, `U_th`, `U_z` | Volume RMS of speed / azimuthal / vertical |
| `T_std` | Temperature std. dev. (see §5.3: must accumulate in `double`) |
| `T_out` / `T_ctr` | Mean T for `r > 0.75 R` vs `r < 0.25 R` |
| `T_bot` / `T_top` | Mean T in lower / upper 25% of melt height |

CSV: `output/french_spatial_stats.csv` (overwrite `french_T_rz.csv`). Screen fields were **appended** to existing lines; cadence unchanged (natural and stirring: every 20 advection steps). This **lengthens lines**, it does **not** add more `[stirring-em]` / `[stage4.1]` rows. Terminal flood is still the pre-existing φ PCG / `edge_flux_*` dump on each EM update.

---

## 3. Natural-convection hydro (blocking bug, then fix)

A previous natural run lost many particles. The 2026-09-14 periodic smoke reproduced the same class of failure until hydro was fixed.

**Cause (not “wrong Boussinesq sign” as the first guess):** WCSPH melt in an open cup, dummy wall too thin / no rim, density not hydrostatic, so startup slosh escaped and `U_max` could look like free fall (`~ g t`).

**Fix (kept an open top; no lid):**

- Dummy wall thickness factor **4** (factor 2 leaked through the floor)  
- Short **rim** `6 dp` **above** the free surface (splash guard, not a lid)  
- Hydrostatic density, then **gravity-only settle** (~0.15 s) before Q / Boussinesq; clock reset to 0  
- `c0 = 10 √(gH)` (~13.5 m/s at H = 0.185 m); `u_ref = √(gH)` so advection starts immediately  
- `flow_ok`: `U_max ≤ c0` and escaped fraction ≤ 1%

After the fix: `escaped = 0 / 10474` on 0.5 s, 5 s, and 50 s periodic runs.

---

## 4. Runs actually completed (user machine)

All commands were run by the user; the agent does not compile or execute these tests. Periodic settings unless noted: `--em-control=fixed-current`, `--aind=off`, SPH conduction on, no `--balance-heat-loss`.

### 4.1 Natural (`frozen_q`, natural reload, `dp = 15 mm`)

| Run | End time | EM interval | Result |
|-----|----------|-------------|--------|
| Smoke B too short | 0.03 s | 0.005 s | `passed` path, but **no advection** → `coupling_updates = 1` (not a coupling demo) |
| Hydro check | 0.5 s | 0.05 s | `escaped = 0`, `U_max ≪ c0` |
| Conduction short | 5 s | 0.5 s | 10 EM updates; P 50.00 → 49.68 kW; `T_mean` −3 K; `dT` → 22 K; loss ~119 kW vs 50 kW Joule |
| Longer periodic (spatial stats + `double` `T_std`) | **50 s** | **2 s** | `passed = 1`; 25 EM updates; I **168.46 A**; P 50.00 → **46.46 kW**; `T_std` 0.09 → **35.2 K** (fix works) |

Natural 50 s endpoint (`french_spatial_stats.csv` + final screen line):

```text
escaped=0/10474  coupling_updates=25  phi_eq_res_vol=2.0e-4
P=46460 W  power_rel_err=0.071  σ_mean 16→14.67  σ_min→11.46
U_max=0.089  U_rms=0.022  U_th=0.006  U_z=0.017  U_buoy=0.053
T_min=1340  T_mean=1443  T_max=1496  dT=156  T_std=35.2
T_out=1447  T_ctr=1416  T_bot=1439  T_top=1441
loss=111 kW (side 70 / bottom 6.8 / free conv 3.7 / rad 31) vs Joule 46 kW
```

- `U_rms` is **0.020 already at t ≈ 0.07 s** and stays **0.021–0.022** for 50 s. It does **not** grow with `dT` or `U_buoy`.  
- `U_max` remains ~9 cm/s WCSPH residual. `U_z_rms ≈ 1.7 cm/s` is the same hydrostatic/acoustic leftover, not a roll cell.  
- `U_buoy` **did** grow (0 → 5.3 cm/s) because `dT` grew; it is a scale only.  
- `T_out > T_ctr` here: the outer bin (`r > 0.75 R`) is mostly induction skin on a 15 mm mesh; the cold wall is ≲ one `dp` and does not dominate that average. Stirring at 5 mm had `T_out < T_ctr` because the wall BL is better resolved and mixed.  
- Q_outer/center stayed ~17–18 (skin in **Q**).  
- **Do not lengthen this run hoping a Jacoutot convective cell appears.**

### 4.2 Stirring (`french_stirring_em`, stirring reload, `dp = 5 mm`)

`c0 = 8 m/s` is a **stability cap** (reload overlap → pressure spikes ∝ `c0²`; `c0 = 10` blows up). `U_tip = ω r_tip ≈ 0.209 m/s` at 10 rpm. `U_max > U_tip` is contact-spike behaviour, not bulk fluid faster than the blade. Use `U_rms` / `U_th` for stirring strength.

| Run | End time | EM interval | Result |
|-----|----------|-------------|--------|
| Periodic smoke | 1 s | 0.5 s | `passed = 1`, `update = 2` at t ≈ 0.50, P ≈ 59.96 kW, I = 132.68 A |
| One revolution | 6 s | 1 s | `passed = 1`; P 60.00 → 59.20 kW (−1.3%); σ_mean 24.75 → 24.53; `T_mean` −4 K; `T_min` −48 K; loss ~146 kW vs 59 kW; wall ~54 s |
| **Five revolutions** | **30 s** | **2 s** | See table below |

Five-revolution endpoint (screen final line):

```text
passed=1  physical_time=30.000  revolutions=5.000  coupling_updates=15
I=132.68 A locked  P_joule=54806 W  power_rel_err=0.0866  phi_eq_res_vol=4.98e-4
U_max=0.576  U_rms=0.123  U_th=0.074  U_z=0.069  U_tip=0.209
T_min=1339  T_mean=1453  T_max=1477  T_out=1441  T_ctr=1461  T_bot=1457  T_top=1435
loss≈138 kW (side 93 / bottom 6.8 / free conv 4.2 / rad 34) vs Joule 55 kW
wall_clock_h=0.067  acoustic_steps=79968  diverged=0  out_of_grid=0
```

Power vs time (I fixed; from coupling CSV):

| t (s) | P (kW) | σ_mean | σ_min |
|-------|--------|--------|-------|
| 0 | 60.00 | 24.75 | 24.75 |
| 6 | 59.24 | 24.54 | 23.45 |
| 16 | 57.24 | 23.95 | 20.83 |
| 28 | 54.81 | 23.18 | 18.74 |

Q_outer / Q_center stayed ~15–16 (skin in **Q**). **`T_out < T_ctr`** because the outer bin includes the **cold wall boundary layer**, not the induction-skin temperature peak. **`T_top < T_bot`** from free-surface radiation. `U_rms ~ 0.08–0.13` m/s ≈ 40–60% of tip → bulk motion exists; `Re_imp ≈ 19` (creeping). Mechanical power `Tz·ω` is ~0.4 W vs ~55 kW EM.

**Do not extend to 1200 s without skull:** bulk cooling continues at ~0.7 K/s.

---

## 5. Diagnostics lessons

### 5.1 What `U_max` is allowed to mean

| Case | Do not use `U_max` as | Use instead |
|------|----------------------|-------------|
| Natural | Convective cell speed | `U_rms` vs `U_buoy`; `escaped`; `dT`; `T_out`/`T_ctr` |
| Stirring | Impeller / mixing speed | `U_rms`, `U_th`, `U_z` vs `U_tip` |

### 5.2 Screen volume

Spatial fields do **not** print extra rows. Quieter SSH: `--screen-every=100` on stirring (natural cadence is hard-coded every 20 advection steps). φ PCG blocks once per EM update are unchanged.

### 5.3 `T_std = 0` bug (fixed; confirmed on this 50 s run)

One-pass `Var(T) = E[T²] − E[T]²` in `Real` (float32) at T ≈ 1470 K cancelled. Accumulators are now **`double`**. The 50 s natural rerun shows a smooth `T_std` 0.09 → 35.2 K, not a string of zeros.

---

## 6. Claim ladder (honest)

| Gate | Meaning | After this work |
|------|---------|-----------------|
| **A** — EM operator trust | Local J/E vs closed-form / TEAM7 | **Unchanged.** Skin-depth `J_rel ≈ 0.87` still open. Not blocking this coupling prototype. |
| **B-light** — Fig. 4-like loop | Periodic σ(T) → EM → Q in SPH | **Closed as a prototype:** both cases re-solve EM, book `P_recon`, lock I. |
| **B-full** | Paper 3D FLUENT ↔ 2D axisymmetric OPHELIE | **Open.** We average σ to (r,z) then still solve **3D** SPH EM. |
| **C** — quantitative paper | Figs, skull, Np/Nu, ~1000 s | **Open.** Over-cooling without skull; no convective cell; creeping mix. |

`--balance-heat-loss` still exists as an **engineering** heat-budget knob. It was **not** used on the runs in §4 and must not be cited as Jacoutot physics.

---

## 7. Still deferred (same list as the brief)

- Skull / free-surface crust (required for any long-time 1473 K claim)  
- `A_metal`, segmented cold-crucible carriers  
- TEAM7 as a French blocker  
- Blocking 2–5% energy residual  
- `k(T)`, `μ(T)`, paper Np / Nu  
- Natural convection at stirring `dp = 5 mm` (would need a **new** natural relax; not started)  
- Periodic stirring with `--aind=on`  
- ParaView campaigns on SSH  

---

## 8. Status after the spatial 50 s natural run

The command below **has been run** (`passed=1`, 25 updates). Item 6 of the ChatGPT brief is closed as: **coupling works; there is no literature convective cell at `dp = 15 mm` / current WCSPH `c0`.** Do not add physical time. Do not re-relax natural to 5 mm unless a new campaign is explicitly started.

The remaining physics bottleneck is **thermal closure** (skull or an explicitly disclaimed effective loss), not sound speed vs tip speed, and not another 50 s.

Historical command (already executed):

### Commands (user runs)

```bash
cd /home/yongchuan/sphinxsys/build
ninja test_3d_ophelie_french_natural_convection_frozen_q
```

```bash
cd /home/yongchuan/sphinxsys/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/bin

./test_3d_ophelie_french_natural_convection_frozen_q \
  --reload-dir=../../test_3d_ophelie_french_natural_glass_relax/bin/reload \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-edge-flux-complex=1 \
  --thermo-em-coupling=periodic \
  --em-update-interval=2 \
  --em-control=fixed-current \
  --aind=off \
  --target-absorbed-power=50000 \
  --end-time=50 \
  --state-record-interval=0
```

Copy (KB, not VTP): `output/french_spatial_stats.csv`, `french_thermo_em_coupling.csv`, `french_energy_budget.csv`.

---

## 9. File index

| Path | Role |
|------|------|
| `docs/ophelie/09_FRENCH_LITERATURE_REPRODUCTION_STATUS_2026-09.md` | Pre-task status vs CES 2008 |
| `docs/ophelie/10_THERMO_EM_COUPLING_ROUND1.md` | Round-1 audit, CLI, early commands |
| `docs/ophelie/11_THERMO_EM_COUPLING_AFTER_CHATGPT_TASK.md` | This close-out |
| `docs/ophelie/12_THERMO_EM_COUPLING_CHATGPT_DISCUSSION_2026-09.md` | Shorter paste pack for ChatGPT (questions at end) |
| `tests/.../test_3d_ophelie_french_natural_convection_frozen_q/README.md` | Natural flags |
| `tests/.../test_3d_ophelie_french_stirring_em/README.md` | Stirring flags |

Papers: Jacoutot *Chem. Eng. Process.* 47 (2008) natural ~50 kW / 282 kHz; Jacoutot *Chem. Eng. Sci.* 63 (2008) stirring ~60 kW / 10 rpm / 300 kHz.
