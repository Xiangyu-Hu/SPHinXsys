# French Reduced Three-Way Comparison (H / A / G) — English Summary

> Consolidated English summary of `FRENCH_THREE_WAY_SUMMARY.txt` (v1), `_v2.txt`, `_v3.txt`.  
> All three `.txt` files are now English; this file remains the single cross-version index.

---

## v3 (latest — Stage 2 + user reload, 2026-06-08)

Same geometry: glass R=0.325 m, H=0.50 m, dp=0.02, n_glass≈20500, f=300 kHz, σ=16 S/m.  
50 kW calibration runs: `--literature-mode --reload=1`.

| Case | Primary gate | phi_eq_res | edge_res_red | P @ 50 kW | divJ_L2_red* | production_lit | passed |
|------|--------------|------------|--------------|-----------|--------------|----------------|--------|
| **H** full edge-flux | edge_res + P_recon + q_antisym + Q_spatial | **9.9e-5** | **~9340** | P_recon | 0.13† | **1** | **1** |
| **A** div-grad baseline | divJ_L2_red + P_particle | **0.427** | ~0.02‡ | P_particle | **2.34** | **1** | **1** |
| **G** phi-only edge-flux | divJ_L2_red (legacy gate) | **9.4e-5** | ~477 | P_particle | **0.93** | 0 | **0** |

\* H divJ is diagnostic only (`--ophelie-particle-gradient-diagnostics=1`); **not** a production gate.  
† H production path does not use particle divJ; with particle diag, divJ_L2_red≈0.13.  
‡ A edge_res is legacy diagnostic; not directly comparable to H.

### H additional acceptance (Stage 2)

| Item | Value |
|------|-------|
| q_antisym_rel_l2 | ~7.1e-8 |
| Q_outer/center | ~14.1 |
| P_graph/P_recon | ~1.15e5 (diagnostic only) |
| French I scaling | P(0.5I)/P(I)=0.25, P(2I)/P(I)=4.0 ✓ |

### Log files (local tee, may not be in repo)

- H 50 kW v4: `french_H_edge_flux_production_v4.log`
- H scaling: `french_H_edge_flux_scaling_reload.log`
- A: `french_A_baseline.log`
- G: `french_G_edge_flux_phi_only.log`

---

## v2 (post power-fix, build default reload)

Same table structure as v3; H production **passed=1** after switching calibration from graph energy to **P_recon**.

See original: `FRENCH_THREE_WAY_SUMMARY_v2.txt`

---

## v1 (pre power-fix)

H run used **wrong P_graph calibration** → `production_literature_passed=0`, `P_graph/P_recon ≈ 1.8e6`.

See original: `FRENCH_THREE_WAY_SUMMARY.txt`

---

## Interpretation

- **H** = production edge-flux + P_recon calibration → **use this route**
- **A** = div-grad reference; good divJ metric but not edge-physics production
- **G** = phi-only edge-flux ablation → **expected fail**

Full context: [`../docs/ophelie/03_EDGE_FLUX_PRODUCTION.md`](../docs/ophelie/03_EDGE_FLUX_PRODUCTION.md)
