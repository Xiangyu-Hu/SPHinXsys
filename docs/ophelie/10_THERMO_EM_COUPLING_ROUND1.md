# Thermo–EM coupling round 1 (2026-09-14)

> Branch: `feature/EM-Biot-Savart`  
> Implements the **minimum closed loop** from the ChatGPT roadmap: unified French parameters, periodic thermo–EM controller, σ(T) azimuthal average, two EM control modes, Q-map power books, natural-convection wiring. Frozen-Q remains the default regression.

---

## 0. Review of the ChatGPT plan (accepted, with scope cuts)

The direction is correct (Gate B-light: periodic SPH thermo–EM, natural first, then stirring). Two adjustments were applied so round 1 stays implementable:

| ChatGPT item | Decision |
|--------------|----------|
| §8 full k(T)/cp(T)/μ(T) + 2–5% energy residual as blocking | **Deferred.** Round 1 wires optional SPH conduction (`--enable-thermal-diffusion`) with **constant** k, cp, μ at 1473 K. Energy CSV is written; residual gate is not enforced yet. |
| Stirring “maybe y-up” | **Incorrect for current meshes.** CAD was −y; `glass-z.stl` / paddle `_z` are **+z**. `FrenchVerticalAxis::Z` is the default, and the axis is stored in the unified parameter struct (not hard-coded in the averager). |
| Silent Q renormalization | **Removed as default** on stirring; opt-in `--q-conservative-remap` with explicit P_before / scale / P_after. |
| Round 1 = 3–5 natural coupling updates | Code path is in. **You must run** it (commands below). Default `--thermo-em-coupling=off` does not change the frozen-Q regression. |
| Do not mix RH200 geometry | Euler CIC grid is **reused** from `rh200_joule_heat_grid.h` (numerics only), not RH200 STL/physics. |

`--balance-heat-loss` is unchanged and is **not** used to claim paper reproduction.

---

## 1. Audit (before / as implemented)

### 1.1 What already existed (reused)

| Component | Location |
|-----------|----------|
| Complex edge-flux + `P_recon` coil-current calibration | `electromagnetic_ophelie_french_literature.h` (`calibrateFrenchCoilCurrentToTargetPower`) |
| EM or A_ind Picard handoff | `runFrenchReducedEmOrSelfInductionForThermalHandoff` |
| σ(T) laws (CEP-anchored 16 S/m, thesis III.12 / IV.10) | `electromagnetic_ophelie_french_material_laws.h` |
| Stage 3.2 σ under-relaxation (reduced one-way, not Fig. 4) | `electromagnetic_ophelie_sigma_t_coupling.h` |
| Robin / radiation BC | `electromagnetic_ophelie_thermal_diffusion_one_way.h` |
| Pairwise thermal Laplace | `OpheliePairwiseLaplaceCK` |
| Euler CIC Q grid + trilinear sample | `rh200_joule_heat_grid.h` |
| Natural frozen-Q + WCSPH + Boussinesq | `test_3d_ophelie_french_natural_convection_frozen_q` |
| Stirring one-shot EM → Euler Q → SPH | `test_3d_ophelie_french_stirring_em` |

### 1.2 Data flow (before)

**Natural frozen_q**

```text
reload glass (SolidBody EM)
  → uniform σ → EM (+ A_ind default ON) → calibrate to 50 kW
  → copy Q by particle index onto FluidBody
  → WCSPH + Boussinesq + Robin/radiation
  → Q never updated; no SPH conduction (enable_diffusion=false)
```

**Stirring**

```text
reload → EM once (A_ind default ON) → calibrate 60 kW
  → CIC Euler Q, then silent scale so sampled power = target
  → each advection: sample Q(x) from grid
  → no SPH conduction
```

### 1.3 Vertical axis

Both production French cases are **+z cylinders** after CAD rotation. Averager takes `FrenchVerticalAxis` (default Z). Do not assume x/y without reading that field.

### 1.4 Material laws actually used

| Case | σ(T) | k, cp, μ |
|------|------|----------|
| natural frozen_q | CEP journal-anchored law, but **EM used constant σ=16** | constants at 1473 K |
| stirring | IV.10 at T0 for σ_glass assign; **not updated in time** | constants; ρ=2800 vs Table-1 2750 |
| reduced EM | typically constant 16 | n/a |

### 1.5 Thermal conduction

**Confirmed missing in the long stir and in frozen_q production:** `thermal_bc.enable_diffusion = false`; only Q source + boundary losses. Round 1 adds `--enable-thermal-diffusion` on the natural case only.

### 1.6 Q grid integral

CIC deposit + trilinear sample. `P_euler` in the book is **particle-sampled** grid power, not a naive ∑Q h³. Axisymmetric power is π(r₁²−r₀²)Δz ∑ Q.

### 1.7 Calibration quantity

**`P_recon` / `probe.joule_power_w`**, then `I ← I √(P_target/P_raw)`. Not graph Laplace energy.

### 1.8 A_ind default

| Mode | Default |
|------|---------|
| frozen-Q regression | **ON** (unchanged) |
| `--thermo-em-coupling=periodic` | **OFF** unless `--aind=on` |
| stirring | still ON unless `--no-self-induction` |

A_ind on does **not** mean self-induction is physically validated.

### 1.9 OPEN_PARAMETER (not guessed)

- Natural f = **282 kHz** vs paper 300 kHz  
- Natural H = **0.185 m** vs stirred **0.21 m**  
- ρ **2750** vs **2800** kg/m³  
- σ(1473 K) **16** vs IV.10 **~24.75** S/m  
- Generator **400 kW** ≠ glass absorbed power  
- Np diameter / paper Nu formula  

---

## 2. New / modified files

| File | Role |
|------|------|
| `extra_src/.../electromagnetic_ophelie_french_literature_parameters.h` | Unified two-regime parameters + OPEN notes |
| `extra_src/.../electromagnetic_ophelie_thermo_em_coupling.h` | r–z averager, under-relax, power book, CSV, `shouldUpdateEM` |
| `test_3d_ophelie_french_natural_convection_frozen_q.cpp` | Periodic coupling CLI + loop |
| `test_3d_ophelie_french_stirring_em.cpp` | No silent Q scale; `--q-conservative-remap` |
| `docs/ophelie/10_THERMO_EM_COUPLING_ROUND1.md` | this file |
| natural / stirring READMEs | CLI defaults |

Stirring does **not** yet call the periodic controller (same header is ready to attach next).

---

## 3. New data flow (periodic natural)

```text
t = 0: EM at uniform σ → calibrate I (FixedCoilCurrent freezes I)
     → σ_bar(r,z), Q(r,z), Euler CIC grid, power book (no silent rescale)

each advection:
  sample Q from Euler grid → particle JouleHeat   (lab frame)
  Q heating (+ optional Laplace diffusion)
  Robin/radiation, Boussinesq

if physical_time >= last_EM + interval:
  T → σ(T) → azimuthal σ_bar(r,z) → under-relax onto EM particles
  re-solve EM (I fixed OR re-calibrate if FixedAbsorbedPower)
  refresh Q(r,z) + Euler grid + CSV
```

Power labels in logs/CSV:

```text
generator_power                 (400 kW, not a gate)
target_glass_absorbed_power     (50 kW natural / 60 kW stirred)
reconstructed_glass_power       (P_recon)
coil_current                    (A/loop)
```

---

## 4. CLI (natural convection)

Defaults: `--thermo-em-coupling=off` (old frozen-Q), A_ind on, 50 kW, no diffusion.
`--thermo-em-coupling=periodic` turns SPH conduction **on** unless `--no-thermal-diffusion`.

```text
--thermo-em-coupling=off|periodic     default off
--em-update-interval=<s>              default 10
--sigma-under-relaxation=<0..1>       default 0.3
--em-control=fixed-current|fixed-power   default fixed-current
--target-absorbed-power=<W>           alias of --target-power=
--aind=off|on                         periodic default off
--enable-thermal-diffusion            periodic: on by default
--no-thermal-diffusion                keep Q + BC only
--constant-sigma
--zero-q
--no-boussinesq
--q-conservative-remap
--state-record-interval=<s>           0 = no VTP (default)
```

---

## 5. Commands for you to run (not run by the agent)

From `build/`:

```bash
# compile
cmake --build . --target test_3d_ophelie_french_natural_convection_frozen_q \
  test_3d_ophelie_french_stirring_em -j$(nproc)
```

**A. Frozen-Q regression (must still pass)**

```bash
RELAX=./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_glass_relax/bin/reload

./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/bin/test_3d_ophelie_french_natural_convection_frozen_q \
  --reload-dir="$RELAX" \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-edge-flux-complex=1 \
  --target-power=50000 --end-time=0.02
```

**B. Periodic coupling smoke (need ≥3 EM updates: interval 0.005, end-time ~0.03)**

```bash
cd tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/bin

./test_3d_ophelie_french_natural_convection_frozen_q \
  --reload-dir=../../test_3d_ophelie_french_natural_glass_relax/bin/reload \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-edge-flux-complex=1 \
  --thermo-em-coupling=periodic \
  --em-update-interval=0.005 \
  --sigma-under-relaxation=0.3 \
  --em-control=fixed-current \
  --aind=off \
  --target-absorbed-power=50000 \
  --end-time=0.03
```

Look for `output/french_thermo_em_coupling.csv`, `french_sigma_rz.csv`, `french_q_rz.csv`, `french_energy_budget.csv`.

**C. Constant-σ check:** add `--constant-sigma` to B; reconstructed power should stay near the first update.

If relax reload is missing, generate it with `test_3d_ophelie_french_natural_glass_relax` first (same dp/geometry as the README).

---

## 6. Round 1 acceptance

| Item | Status |
|------|--------|
| Unified parameters + OPEN list | done |
| Controller + physical-time `shouldUpdateEM` | done |
| σ(T) + under-relax + empty-bin fill | done |
| FixedCoilCurrent / FixedAbsorbedPower | done |
| Q map power book, no silent rescale | done (natural + stirring) |
| Natural periodic path | hydro contained (`escaped=0`, `U_max≪c0` at 0.5 s) |
| Frozen-Q default preserved | yes |
| 3–5 coupling updates | smoke B / 0.5 s hydro check |
| Energy / σ / Q / EM CSV | written in periodic mode |
| SPH conduction on periodic | **round 2** (default on; `--no-thermal-diffusion` to disable) |
| Energy residual columns | **round 2** (reported, not a pass gate) |
| Stirring periodic attach | **round 2** (default still one-shot) |
| k(T), μ(T), full energy 2–5% gate | still deferred |
| Skull / A_metal / TEAM7 / RH200 | still deferred |

---

## 7. Remaining risks

- EM SolidBody and WCSPH FluidBody are two SPHSystems; GPU must keep both alive (same pattern as stirring).  
- Particle-index Q freeze is still used when coupling=off.  
- Axisymmetric σ is applied to a **3D** EM solve (Fig. 4-like dataflow, not 2D OPHELIE).  
- Skin-depth J_rel ≈ 0.87 unchanged.  
- Heat-loss vs 50 kW without skull still over-cools on long runs. Do **not** treat `--balance-heat-loss` as paper physics.  
- Natural WCSPH: open-top cup now has a rim + hydrostatic settle + `U_max ≤ c0` / escaped-particle gate. Do **not** add a lid. Stirring still uses a flush rim (`rim_height=0`).  
- `energy_residual` on 20-advection windows is SPH acoustic noise; look at a long-time mean, not a 2–5% gate.  
- `U_buoy = √(g β ΔT H)` is a scale, not a literature velocity. At 0.5 s, `U_max` is still residual WCSPH, not a convective cell.

## 8. Round 2 (this drop)

Hydro containment is in. Round 2 does **not** add a lid or thicker walls.

1. Periodic natural convection: SPH conduction **on by default**; `--no-thermal-diffusion` restores the old Q+BC-only path. Frozen-Q regression stays off.  
2. Energy CSV now has `energy_residual`, `energy_residual_rel`, `T_min/mean/max`, `U_max`, `U_buoy`. Screen prints the same. Not a pass/fail gate.  
3. Monitor CSV: `output/french_natural_convection_monitor.csv`. Optional VTP: `--state-record-interval=1`.  
4. Same controller on `test_3d_ophelie_french_stirring_em` (`--thermo-em-coupling=periodic`). One-shot 10 rpm default is unchanged, including A_ind on.

**D. Natural convection with conduction (you run this)** — from `frozen_q/bin`, same geometry/reload as C. Conduction is implied by `periodic`. Five physical seconds is still short vs a convective cell; watch `escaped`, `U_max` vs `U_buoy`, `dT`.

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
  --target-absorbed-power=50000 --end-time=5 \
  --state-record-interval=1
```

**E. Stirring periodic smoke** — uses the stirring relax reload, not the natural one. Keep `--end-time` small; EM re-solves are expensive at dp=5 mm.

```bash
cd tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_em/bin
./test_3d_ophelie_french_stirring_em \
  --reload-dir=../../test_3d_ophelie_french_stirring_glass_relax/bin/reload \
  --thermo-em-coupling=periodic --em-update-interval=0.5 \
  --em-control=fixed-current --aind=off \
  --end-time=1 --max-wall-hours=1 --no-state-recording
```

## 9. After this round

Remote / SSH: **do not download VTP**. Screen lines now include `U_rms`, `U_th`, `U_z`, `T_std`, `T_out`, `T_ctr`. Copy `output/french_spatial_stats.csv` (KB) and optionally `french_T_rz.csv` / `french_thermo_em_coupling.csv`. Use `--no-state-recording`.

1. Natural 50 s periodic already contained; `U_max` is residual WCSPH, not a cell.  
2. Stirring 1-rev periodic smoke passed; σ(T) dropped P by 1.3%.  
3. Longer stirring only with `--no-state-recording`; do not treat `--balance-heat-loss` as paper physics.  
4. Gate A (skin / probes) stays on a parallel track.  
5. Still deferred: skull, A_metal, TEAM7 blocking, 2–5% energy gate, k(T)/μ(T).

Runs, hydro fix, 5-rev stirring, and spatial diagnostics after this file: [`11_THERMO_EM_COUPLING_AFTER_CHATGPT_TASK.md`](11_THERMO_EM_COUPLING_AFTER_CHATGPT_TASK.md).
