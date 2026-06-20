# RH200 Glass EM Stirring (Midterm Case)

> Consolidated from: `OPHELIE_RH200_GLASS_EM_STIRRING_CURSOR_TASKS.md`, `RH200_OPHELIE_STATUS_AIND_CURRENT_AND_NEXT_CURSOR_PLAN.md`, `RH200_JOULE_HEAT_EULERIAN_GRID_CURSOR_PLAN.md`, `RH200_EXCITATION_JOULE_HEATING_AUDIT_CURSOR_PLAN.md`, client LaTeX report.  
> Archive: [`archive/plans/RH200_*`](archive/plans/), [`archive/plans/OPHELIE_RH200_*`](archive/plans/)

---

## 1. Case definition

```text
RH200 STL geometry (oil.stl → glass fluid, rotor, wall)
+ OPHELIE/French-style analytic equivalent coil
+ complex edge-flux EM (one-way A_coil)
+ target-power Joule calibration (50 kW glass domain)
+ Eulerian fixed-grid Q(x) and J(x) sampling (em-grid)
+ SPH weakly compressible flow + Simbody rotor + thermal diffusion
```

**Not:** strict French cold-crucible reproduction; **not** 1500 K in 10 s acceptance.

---

## 2. Executable and modes

**Path:** `tests/extra_source_and_tests/3d_examples/test_3d_ophelie_rh200_glass_em_stirring/`

| CLI | Purpose |
|-----|---------|
| `--relax=1 --dp=...` | Particle relaxation → `./reload` |
| `--reload=1 --em-solve=1` | EM-only diagnostic VTP |
| `--reload=1 --run=1` | Flow + thermal + Joule |

| `--joule-mode` | Description |
|----------------|-------------|
| `off` | No Joule heating |
| `em-fixed` | Debug: frozen particle Q (wrong Lagrangian tag) |
| **`em-grid`** | **Production:** fixed Euler grid sample each step |

| `--preset` | T₀ | μ [Pa·s] | Notes |
|------------|-----|----------|-------|
| `rh200-demo-current` | 373 K | 0.05 | Default smoke |
| `rh200-french-like-stable` | 1473 K | 0.05 | Jacoutot props |
| `rh200-french-like-full` | 1473 K | 4 | High viscosity demo |

---

## 3. em-grid coupling (Step 3b)

**Problem:** em-fixed attaches Q to particles → heat source moves with stirred glass (incorrect).

**Solution:**

1. EM solve once on glass reload particles
2. CIC deposit scaled Q, φ, |E|, J components to Euler grid on `full_oil.stl` bbox
3. Each acoustic step: trilinear sample `Q(x_i(t))`, `JReal/JImag(x_i(t))`
4. Explicit: `T_i += Q_i Δt / (ρ₀ cp)`

Grid spacing default: `1.0 × dp`. Domain margin: `--joule-grid-bbox-margin=2`.

---

## 4. Excitation–field–Joule consistency (audit conclusion)

From EM logs and audit CSV:

```text
P_raw ≈ 4.736 W
I_scale ≈ 102.746  →  em_power_scale ≈ 10557 ≈ I_scale²
P_scaled = 50000 W (matches --target-power)
grid_rescale_factor ≈ 1 (em-grid)
phi_eq_res_vol ~ 4e-4
```

**Conclusion:** coil excitation, E/J fields, and Joule heat are **consistent at power-book level**. No obvious unit or scaling order error. Pointwise J·E residual audit (Task B) still optional.

---

## 5. Representative results

### Demo preset (T₀=373 K, 10 s, em-grid)

```text
P_grid ≈ 50 kW
ΔT_mean ≈ 0.85 K
energy closure ≈ 97%
```

### French-like-full (T₀=1473 K, μ=4, 60 s)

```text
P_grid ≈ 49.5–49.7 kW
ΔT_mean ≈ 0.75 K (1473 → 1473.75 K)
T_max ≈ 1476 K
heating_rate_expected ≈ 0.083 K/s (no-loss estimate)
```

---

## 6. Audit outputs (`./output/`)

| File | Content |
|------|---------|
| `RH200_EXCITATION_TO_JOULE_AUDIT.csv` | raw/scaled EM stats, scale factors |
| `RH200_HEATING_RATE_AUDIT.csv` | dT measured vs expected |
| `rh200_energy_budget.csv` | P_joule, thermal energy, integrated Joule |
| `rh200_joule_heat_grid_summary.csv` | Grid build report |
| `JouleHeatGrid.vti` | Q, φ, |E|, J component grids |
| `GlassGeometry_ite_*.vtp` | T, P, JouleHeat, **JReal, JImag** (needs `--state-recording=1`) |

---

## 7. Recommended server run (high resolution)

```bash
# Re-relax when changing dp
./test_3d_ophelie_rh200_glass_em_stirring --relax=1 --dp=0.005

./test_3d_ophelie_rh200_glass_em_stirring \
  --preset=rh200-french-like-full \
  --reload=1 --run=1 --dp=0.005 --end-time=60 \
  --joule-mode=em-grid --target-power=50000 \
  --state-recording=1 --state-record-interval=1.0 \
  --joule-grid-bbox-margin=2
```

---

## 8. Task checklist (from cursor plans)

| Task | Status |
|------|--------|
| A — EXCITATION_TO_JOULE audit CSV | ✅ |
| C — HEATING_RATE audit CSV | ✅ |
| em-grid Eulerian Q | ✅ |
| French-like presets | ✅ |
| JReal/JImag VTP + live em-grid mapping | ✅ |
| B — J·E formula residual audit | ⏸ Optional |
| A_ind Picard diagnostic | ⏸ Future |

Detail: [`archive/plans/RH200_OPHELIE_STATUS_AIND_CURRENT_AND_NEXT_CURSOR_PLAN.md`](archive/plans/RH200_OPHELIE_STATUS_AIND_CURRENT_AND_NEXT_CURSOR_PLAN.md).

---

## 9. Client report

Full technical report (Chinese LaTeX, equations, validation matrix):

`RH200_OPHELIE_CLIENT_PROGRESS_REPORT.tex` / `.pdf` at repo root.

English consolidated summary in [01_MASTER_DEVELOPMENT_TIMELINE.md](01_MASTER_DEVELOPMENT_TIMELINE.md).
