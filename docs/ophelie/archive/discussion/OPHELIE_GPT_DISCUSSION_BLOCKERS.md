# OPHELIE Edge-Flux — Pending GPT discussion Summary

> Updated: 2026-06-08. Completed items: see `OPHELIE_PHI_DEVELOPMENT_LOG.md` §4.14. 
> **Intent**: The items below **do not block** current Stage 1 production / Stage 2 A_ind diagnostic work; collect them for a single GPT discussion later.

---

## Closed (this round)

| Topic | Conclusion |
|------|------|
| Can `P_total_edge` be used for 50 kW calibration? | **No**; switched to `P_recon`; H production `passed=1` |
| Did edge φ / edge residual fail? | **No**; `edge_res_red≈10⁴`, `phi_eq_res≈8e-5` |
| Is `P_recon` trustworthy on uniform field? | **Yes**; `P_recon/P_exact≈1` (potential + pure induction cases) |

---

## Deferred discussion (non-blocking)

### B1. `P_graph/P_recon ≈ 1.16×10⁵` (diagnostic only)

- **Observation**: After H production calibration, `P_graph_edge≈5.8×10⁹ W`, `P_recon≈50 kW`.
- **Done**: No longer a production gate; log `P_graph_over_recon_warning=0`.
- **Open**: Study `P_graph_norm = κ·P_graph` (calibrate κ on uniform-field MMS)? Or permanently deprecate graph Joule heat?

### B2. H particle `divJ_L2_red ≈ 0.116` (baseline A ≈ 2.06)

- **Observation**: Edge-flux production does not improve particle-gradient `divJ`; Jn_post_phi_rel≈0.0084.
- **Cause (known)**: E/J post-processing still uses particle grad, not edge `JEdgeRecon` continuity.
- **Open**: Should Stage 2 switch E/J post-processing to `JEdgeRecon`? Or accept divJ as baseline comparison only?

### B3. `I×2 → P_recon×4` scaling — **closed (2026-06-08)**

- French reload: `P(0.5I)/P(I)=0.25`, `P(2I)/P(I)=4.0` (`--no-literature-calibrate`).
- Uniform field: `test_3d_ophelie_edge_flux_scaling` passed.
- Log: `discussion_bundle/french_H_edge_flux_scaling_reload.log`; §4.16.

### B4. §4.11 historical decision points (div-grad era)

- Is div-grad `phi_eq_res≈0.49` a formulation floor?
- Should edge-flux abandon particle divJ as primary metric?
- See dev log §4.9–§4.11; edge-flux Stage 1 partially closed, div-grad fallback still preserved.

### B5. Code structure (low priority)

- `edge_flux_diagnostics.h` should not include CLI → split to `parameters.h` (Plan §9.2.1)
- Legacy edge diagnostic reuses `DivJImag` → dedicated fields (Plan §9.2.2)

---

## Blockers (none currently)

None. Stage 1 production passed; Step 5 H/A/G v2 and Step 6 A_ind baseline diagnostic completed.

---

## Step 6 record (2026-06-08)

| Quantity | Value |
|----|-----|
| `A_ind/A_coil` | 0.436 |
| `B_ind/B_coil` | 0.366 |
| `phi_eq_res_vol` | 0.410 (div-grad path) |
| `P_joule_W` | 7.34 (no literature calibration) |
| `passed` | 1 |

Log: `discussion_bundle/french_aind_diagnostic_v1.log`

**Open**: Should A_ind use `JEdgeRecon` as K[J] source (related to B2)?

---

## Run artifact paths

| Case | Log |
|------|------|
| H edge-flux production (v2) | `discussion_bundle/french_H_edge_flux_production_v2.log` |
| H (v1, wrong calibration) | `discussion_bundle/french_H_edge_flux_production.log` |
| A baseline | `discussion_bundle/french_A_baseline.log` |
| G phi-only | `discussion_bundle/french_G_edge_flux_phi_only.log` |
| Three-way summary (old) | `discussion_bundle/FRENCH_THREE_WAY_SUMMARY.txt` |
| Three-way summary (v2) | `discussion_bundle/FRENCH_THREE_WAY_SUMMARY_v2.txt` |
| A_ind diagnostic | `discussion_bundle/french_aind_diagnostic_v1.log` |
