# Edge-Flux Production Path

> Consolidated from: `OPHELIE_EDGE_FLUX_SOLVER_DEVELOPMENT_PLAN.md`, `OPHELIE_EDGE_FLUX_POWER_FIX_AND_NEXT_DEVELOPMENT_PLAN.md`, `OPHELIE_EDGE_FLUX_STAGE2_CLEANUP_AND_CURSOR_PLAN.md`, `OPHELIE_COMPLEX_EDGE_FLUX_SOLVER_PLAN.md`, `discussion_bundle/OPHELIE_COMPLEX_EDGE_FLUX_RUN_RECORD.md`, `discussion_bundle/OPHELIE_STAGE2_CURSOR_WORK_RECORD.md`, `discussion_bundle/OPHELIE_GPT_ROUND_SUMMARY_2026-06-08.md`.  
> Archive: [`archive/plans/OPHELIE_EDGE_FLUX_*`](archive/plans/), [`archive/plans/OPHELIE_COMPLEX_EDGE_FLUX_*`](archive/plans/)

---

## 1. Why edge-flux replaced div-grad φ

| Route | French real Biot field | Issue |
|-------|------------------------|-------|
| div-grad φ projection | phi_eq_res ≈ 0.49 | Better divJ metric but weak edge physics |
| **edge-flux production** | phi_eq_res ~ 1e-4, P_recon @ 50 kW | **Production choice** (2026-06-08) |

Legacy div-grad code: `electromagnetic_ophelie/legacy/div_grad/` and `docs/ophelie/archive_phi_divgrad_residual_floor/`.

---

## 2. Discrete edge-flux equations

Edge electromotive drop:

```text
e_ij = (φ_j - φ_i) + ω · Ā_ij · (x_j - x_i)
```

Edge conductance: `C_ij = σ̄_ij · w_ij^Lap`

Edge flux: `q_ij = -C_ij · e_ij`

Particle residual: `r_i = Σ_j q_ij → 0`

Primary fields after solve: **EEdgeRecon**, **JEdgeRecon**, **JouleHeatEdgeReconComplex**.

---

## 3. Power calibration (critical fix)

**Wrong (pre-2026-06-08):** calibrate coil current so graph edge energy `P_graph = Σ ¼ C_ij e_ij²` hits 50 kW.

**Correct:** calibrate on reconstructed Joule body integral:

```text
P_recon = ∫ Q_recon dV
I_new = I_old · sqrt(P_target / P_raw)
```

Observed before fix: `P_graph/P_recon ≈ 1.8×10⁶`, `production_literature_passed = 0`.  
After fix: French H **`production_literature_passed = 1`**.

---

## 4. Complex edge-flux mode

Stores separate real/imag fields: `PhiReal/PhiImag`, `EReal/EImag`, `JReal/JImag`.

Acceptance (French H @ 50 kW):

```text
edge_res_red > 1e3
P_recon ≈ 50 kW
phi_eq_res_vol ~ 1e-4
q_antisym_rel_l2 ~ 1e-7
```

---

## 5. French three-way comparison (H / A / G)

Summaries (English): [`discussion_bundle/FRENCH_THREE_WAY_SUMMARY_EN.md`](../../discussion_bundle/FRENCH_THREE_WAY_SUMMARY_EN.md)

Original text files (Chinese headers preserved):

| Version | File | Notes |
|---------|------|-------|
| v3 (latest Stage 2 + user reload) | `discussion_bundle/FRENCH_THREE_WAY_SUMMARY_v3.txt` | H passed @ 50 kW |
| v2 (post power-fix) | `discussion_bundle/FRENCH_THREE_WAY_SUMMARY_v2.txt` | |
| v1 (pre power-fix) | `discussion_bundle/FRENCH_THREE_WAY_SUMMARY.txt` | Wrong P_edge calibration |

| Label | Mode | Expected |
|-------|------|----------|
| **H** | edge-flux + literature calibrate | **PASS** |
| **A** | div-grad baseline | Reference only |
| **G** | edge-flux phi-only | **FAIL** (by design) |

CSV diagnostics: `discussion_bundle/ophelie_edge_flux_*.csv`, `discussion_bundle/ophelie_phi_p0_*.csv`.

---

## 6. Source layout (post Stage 2 cleanup)

| Location | Content |
|----------|---------|
| `electromagnetic_ophelie/` (root) | edge-flux + φ + French + shared infra |
| `legacy/div_grad/` | div-grad fallback |
| `diagnostics/` | MMS / operator audit / A_ind |
| `team7/` | TEAM7-specific helpers |
| `french_extras/` | glass relaxation |
| `stage2/` | self-induction Picard experiments |

---

## 7. Core tests

```text
test_3d_ophelie_edge_flux_sign
test_3d_ophelie_edge_flux_power_uniform_field
test_3d_ophelie_edge_flux_scaling
test_3d_ophelie_french_reduced
test_3d_ophelie_french_complex_joule_to_heat_one_way
```

Full run narrative: [`archive/discussion/OPHELIE_COMPLEX_EDGE_FLUX_RUN_RECORD.md`](archive/discussion/OPHELIE_COMPLEX_EDGE_FLUX_RUN_RECORD.md).
