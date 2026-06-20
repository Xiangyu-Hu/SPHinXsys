# French Reduced Case and Thermal Coupling

> Consolidated from: `FRENCH_OPHELIE_PARTICLE_DEVELOPMENT_PLAN.md`, `OPHELIE_FRENCH_REDUCED_*`, `OPHELIE_FRENCH_GLASS_EM_STIRRING_MIDTERM_CURSOR_PLAN.md`, `docs/ophelie/FRENCH_*`.  
> Archive: [`archive/plans/OPHELIE_FRENCH_*`](archive/plans/), [`archive/plans/FRENCH_OPHELIE_PARTICLE_DEVELOPMENT_PLAN.md`](archive/plans/)

---

## 1. Literature reference (Jacoutot et al.)

| Quantity | Typical value | Notes |
|----------|---------------|-------|
| Crucible diameter | ~650 mm | Cold crucible |
| Generator power | ~400 kW | **Supply-side**, not glass absorption |
| Glass Joule power | ~50 kW | **Target for `--target-power=50000`** |
| Frequency | 300 kHz | |
| σ @ ~1500 K | ~16 S/m | |
| Bath temperature | ~1500 K | Operating state, not cold-start in seconds |

**400 kW ≠ 50 kW** — midterm RH200 uses 50 kW as **absorbed glass-domain power**.

---

## 2. French reduced EM case

**Executable:** `test_3d_ophelie_french_reduced`

**Geometry:** simplified cylinder glass + analytic multiloop coil (not full French CAD).

**Production acceptance (edge-flux + literature calibrate @ 50 kW):**

```text
edge_res_red > 1e3
P_recon ≈ 50 kW
phi_eq_res_vol ~ 1e-4
q_antisym_rel_l2 ~ 1e-7
```

**Docs:**

- [`FRENCH_REDUCED_CASE_README.md`](FRENCH_REDUCED_CASE_README.md)
- [`FRENCH_LITERATURE_MODE.md`](FRENCH_LITERATURE_MODE.md)
- [`FRENCH_REDUCED_CASE_ASSUMPTIONS.md`](FRENCH_REDUCED_CASE_ASSUMPTIONS.md)

---

## 3. Geometry implementation tasks

| Task doc (archive) | Content |
|--------------------|---------|
| `OPHELIE_FRENCH_REDUCED_PRE_GEOMETRY_TASKS.md` | Pre-geometry checklist |
| `OPHELIE_FRENCH_REDUCED_GEOMETRY_IMPLEMENTATION.md` | STL/relax/reload pipeline |
| `OPHELIE_FRENCH_REDUCED_04_06_REVIEW.md` | 2026-04-06 review notes |

Particle relaxation helper: `electromagnetic_ophelie_french_geometry_sphinxsys_style_relax.cpp`.

---

## 4. Thermal one-way coupling

**Executable:** `test_3d_ophelie_french_complex_joule_to_heat_one_way`

```text
complex EM → JouleHeatEdgeReconComplex → explicit dT/dt → optional isotropic diffusion
```

**Material presets (Jacoutot Table 1 @ 1473 K):**

| Preset | ρ [kg/m³] | cp [J/kg/K] | k [W/m/K] | T₀ [K] |
|--------|-----------|-------------|-----------|--------|
| reduced (regression default) | 2500 | 1200 | 1 | 300 |
| literature | 2750 | 1150 | 4 | 1473 |

RH200 reuses French-like presets via `--preset=rh200-french-like-stable|full` (see [06_RH200_GLASS_EM_STIRRING.md](06_RH200_GLASS_EM_STIRRING.md)).

---

## 5. Midterm pivot: RH200 geometry

French reduced remains the **EM regression** case. Midterm **delivery geometry** switched to RH200 STL stirring (see consolidated RH200 doc). TEAM7 is side-track only.

Plan: [`archive/plans/OPHELIE_FRENCH_GLASS_EM_STIRRING_MIDTERM_CURSOR_PLAN.md`](archive/plans/OPHELIE_FRENCH_GLASS_EM_STIRRING_MIDTERM_CURSOR_PLAN.md).

---

## 6. Related tests

```text
test_3d_ophelie_french_reduced
test_3d_ophelie_french_complex_joule_to_heat_one_way
test_3d_ophelie_french_glass_relax
test_3d_ophelie_french_em_dp_scan
test_3d_ophelie_french_aind_diagnostic
test_3d_ophelie_french_self_induction_picard
```
