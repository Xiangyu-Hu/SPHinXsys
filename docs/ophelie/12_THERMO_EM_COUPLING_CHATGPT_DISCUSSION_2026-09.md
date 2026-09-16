# ChatGPT discussion pack — thermo–EM work since your brief (2026-09-16)

> **Paste this file to ChatGPT** to continue the Jacoutot / OPHELIE roadmap discussion.  
> Branch: `feature/EM-Biot-Savart` (pushed: `d0096015e`).  
> Full close-out with tables: [`11_THERMO_EM_COUPLING_AFTER_CHATGPT_TASK.md`](11_THERMO_EM_COUPLING_AFTER_CHATGPT_TASK.md).  
> Pre-brief snapshot: [`09_FRENCH_LITERATURE_REPRODUCTION_STATUS_2026-09.md`](09_FRENCH_LITERATURE_REPRODUCTION_STATUS_2026-09.md).  
> Midterm claim: **auditable SPH–EM–thermal–flow prototype**, not paper-grade Jacoutot reproduction.

---

## 0. Ask of ChatGPT

Please treat this as a **status update after implementing your brief**, then recommend **one primary next step** (4–8 weeks), with explicit deferrals. Prefer honest gates over “make it look like the paper.”

Questions at the end (§7).

---

## 1. What your brief asked (reminder)

Starting point (already true before the brief):

- Biot–Savart `A_coil`, complex edge-flux φ, E/J/Q, `P_recon` current calibration  
- One-shot EM → Euler Q → SPH; 10 rpm paddle kinematics  

Your next milestone:

> Periodic thermo–EM driven by σ(T); **natural convection first**, then the **same** controller on 10 rpm stirring.

Ordered list you gave:

1. Unify French parameters  
2. Periodic thermo–EM coupling  
3. `T → σ(T) → EM → Q`  
4. Q-map power book (no silent rescale)  
5. SPH conduction  
6. Validate natural convection  
7. Attach stirring  
8. Later: skull / `A_metal` / segmented crucible  

Explicit non-goals (kept): no RH200 mixed into French baseline; no `--balance-heat-loss` as paper physics; no TEAM7 as French blocker; no large SPH rewrite.

---

## 2. What we implemented

### 2.1 Shared controller

| Component | Role |
|-----------|------|
| `electromagnetic_ophelie_french_literature_parameters.h` | Natural 50 kW / 282 kHz vs stirring 60 kW / 300 kHz; `EMControlMode`; `ThermoEMCouplingMode` |
| `electromagnetic_ophelie_thermo_em_coupling.h` | Azimuthal σ(r,z), empty-bin fill, under-relax, `shouldUpdateEM` on **physical time**, power CSV, melt spatial stats |

Dataflow (periodic):

```text
t=0: EM @ uniform σ → calibrate I → freeze I (fixed-current) → Euler CIC Q
each advection: sample Q → heat (+ optional Laplace) → Robin/radiation → Boussinesq
every Δt_EM: T→σ(T)→σ_bar(r,z)→under-relax→re-solve EM→refresh Q grid
```

Corrections vs the brief:

| Item | Decision |
|------|----------|
| Vertical axis | **+z** (not y) for current meshes |
| Silent Q rescale on stirring | **Removed**; `--q-conservative-remap` opt-in |
| Periodic A_ind | Default **off** |
| Natural “lid” | **Rejected** — open free surface + short rim only |
| Full k(T)/μ(T) + 2–5% energy gate | **Deferred** (round-1 cut) |

### 2.2 Cases

| Case | Reload | dp | Periodic status |
|------|--------|-----|-----------------|
| `french_natural_convection_frozen_q` | natural relax only | **15 mm** | Wired; frozen-Q remains default regression |
| `french_stirring_em` | stirring relax only | **5 mm** | Wired; default still one-shot + A_ind on |

**Do not mix** the two `Reload.xml` files.

### 2.3 Hydro fix (natural — was blocking)

Particles previously escaped. Fix: thicker dummy wall (factor 4), short rim `6 dp` above melt, hydrostatic density + gravity settle before Q/Boussinesq. Still an **open cup**, not a lid. After fix: `escaped=0/10474` through 50 s.

### 2.4 Diagnostics without ParaView

SSH VTP download is too slow. Added host spatial stats on existing screen cadence (longer lines, **not** more lines):

- `U_rms`, `U_th`, `U_z`  
- `T_std` (fixed: accumulate in `double`; float32 one-pass variance was falsely 0)  
- `T_out` / `T_ctr` (`r > 0.75R` vs `r < 0.25R`)  
- `T_bot` / `T_top`  

CSV: `french_spatial_stats.csv`, `french_thermo_em_coupling.csv`, `french_energy_budget.csv`.

Git: pushed to `origin/feature/EM-Biot-Savart` as `d0096015e`. Jacoutot PDF was **not** committed (Elsevier copyright on a public repo).

---

## 3. Runs completed (honest numbers)

Settings unless noted: `--thermo-em-coupling=periodic`, `--em-control=fixed-current`, `--aind=off`, SPH conduction on, **no** `--balance-heat-loss`.

### 3.1 Natural — 50 s, `dp=15 mm`, EM every 2 s

```text
passed=1  escaped=0/10474  coupling_updates=25  I=168.46 A locked
P: 50.00 → 46.46 kW   σ_mean: 16 → 14.67   φ_res ~ 2e-4
T_mean: 1473 → 1443 K   dT → 156 K   T_std → 35 K
T_out=1447  T_ctr=1416   (outer hotter — skin on coarse mesh)
loss ≈ 111 kW vs Joule 46 kW
U_max ≈ 0.09 m/s (flat)   U_rms ≈ 0.022 m/s (flat from t≈0)
U_buoy = √(g β ΔT H) → 0.053 m/s  (scale only)
```

**Verdict on brief item 6:** coupling path works; **there is no Jacoutot-like buoyancy cell**.  
Evidence: `U_rms` does not grow with `dT` / `U_buoy`; residual WCSPH noise dominates (~2 cm/s from the first advection). Viscous buoyancy for μ=4 Pa·s is likely mm/s and buried. **Do not lengthen physical time hoping a cell appears.**

### 3.2 Stirring — 30 s = 5 rev, `dp=5 mm`, EM every 2 s

```text
passed=1  revolutions=5  coupling_updates=15  I=132.68 A locked
P: 60.00 → 54.81 kW   (σ drop with cooling)
U_tip=0.209   U_rms≈0.08–0.13   U_max≈0.3–0.6 (contact spikes — ignore for mix)
T_mean: 1473 → 1453   T_out < T_ctr   T_top < T_bot
loss ≈ 138 kW vs Joule 55 kW
Re_imp ≈ 19 (creeping)   wall_clock ≈ 4 min
```

**Verdict on brief item 7:** Gate B-light on stirring is demonstrated (σ→P feedback at fixed I; bulk motion exists). Not paper mixing / not 1000 s steady state. **Do not run 1200 s without skull** (~0.7 K/s bulk cool).

### 3.3 Speed bookkeeping (common confusion)

| Quantity | Role |
|----------|------|
| `c0` | Artificial sound speed (natural ~13.5 m/s; stirring **8** = stability cap) |
| `U_tip` | Paddle tip (~0.21 m/s) — kinematics |
| `U_max` | Peak particle speed — spikes / residual; **not** cell or mix strength |
| `U_rms` / `U_th` | Bulk motion diagnostics |
| `U_buoy` | Buoyancy **scale** from ΔT, not measured velocity |

Sound speed vs tip is **not** the open issue. Over-cooling without skull is.

---

## 4. Claim ladder (do not oversell)

| Gate | Meaning | Status |
|------|---------|--------|
| **A** EM field trust | Local J/E vs analytic / TEAM7 | Unchanged; skin `J_rel≈0.87` still open |
| **B-light** Fig.4-like loop | Periodic σ(T)→EM→Q in SPH | **Closed as prototype** (natural + stirring) |
| **B-full** Paper 3D↔2D OPHELIE | Exact architecture | Open (we keep 3D SPH EM after (r,z) average) |
| **C** Quantitative paper | Figs, skull, Np/Nu, ~1000 s | Open |

`--balance-heat-loss` remains an engineering knob only; unused on the runs above.

---

## 5. Still deferred (same as brief)

- Skull / crust (needed for any long-time 1473 K claim)  
- `A_metal`, segmented cold crucible  
- TEAM7 as French blocker  
- Blocking 2–5% energy residual  
- Full `k(T)` / `μ(T)`, paper Np / Nu  
- Natural at `dp=5 mm` (needs **new** natural relax — not started)  
- Periodic stirring with `--aind=on`  

---

## 6. Our current read (for debate)

1. Brief items 1–5 and 7 are done; item 6 is **architecture OK, physics cell missing**.  
2. Next science bottleneck is **thermal closure** (skull or explicitly disclaimed effective loss), not another 50 s natural run and not tip/sound-speed tuning.  
3. Gate A (skin / probes) can stay parallel; it does not block B-light prototype claims.  
4. Client / RH200 demo polish is orthogonal to French reproduction claims.

---

## 7. Questions for ChatGPT

1. Given no buoyancy cell at `dp=15 mm` WCSPH, should we **declare natural cell out of midterm scope**, try **finer natural dp**, or change hydro (c0 / viscosity / settle) before skull?  
2. Is **fixed skull / effective wall resistance** the correct next primary, or should we push **Gate A** (skin-depth J) first?  
3. For stirring, is **5 rev + σ→P feedback** enough Gate B-light evidence, or do you want Np-like / longer time with `--balance-heat-loss` under a **clear disclaimer**?  
4. Keep midterm sentence as “SPH Fig.4-like prototype” and leave cold-crucible `A_metal` **out of scope**?  
5. Any red flags in `T_out`/`T_ctr` opposite signs between natural (15 mm) and stirring (5 mm)?

---

## 8. Key paths

```text
docs/ophelie/12_THERMO_EM_COUPLING_CHATGPT_DISCUSSION_2026-09.md  ← this pack
docs/ophelie/11_THERMO_EM_COUPLING_AFTER_CHATGPT_TASK.md          ← full close-out
docs/ophelie/10_THERMO_EM_COUPLING_ROUND1.md                      ← CLI / audit
docs/ophelie/09_FRENCH_LITERATURE_REPRODUCTION_STATUS_2026-09.md  ← pre-brief
tests/.../test_3d_ophelie_french_natural_convection_frozen_q/
tests/.../test_3d_ophelie_french_stirring_em/
```

Papers: Jacoutot CEP 47 (2008) natural ~50 kW / 282 kHz; Jacoutot CES 63 (2008) stirring ~60 kW / 10 rpm / 300 kHz. Generator 400 kW ≠ glass absorbed power.
