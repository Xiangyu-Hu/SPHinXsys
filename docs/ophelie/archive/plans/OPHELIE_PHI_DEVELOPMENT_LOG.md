# OPHELIE φ Boundary and Accuracy: Development Log

> Companion to [`OPHELIE_PHI_BOUNDARY_AND_NEXT_DEVELOPMENT_PLAN.md`](OPHELIE_PHI_BOUNDARY_AND_NEXT_DEVELOPMENT_PLAN.md).  
> Records implemented code, run instructions, acceptance numbers, and **deferred/optional** items (retry later).

Last updated: 2026-06-02 (TEAM7 L2 complex edge-flux validation + record)

---

## 1. Current stage

| Sprint | Status | Content |
|--------|------|------|
| **Remediation Sprint 1** | ✅ Done | CLI prefix, `solvePhiImagWithCurrentRhs`, RHS fingerprint, operator linear/discrete MMS tests |
| **Remediation Sprint 2** | ✅ Discrete MMS passed | `discrete_self_consistency_rel≈3.3e-7` |
| **Remediation Sprint 3** | ✅ Implemented / ❌ French failed | paired G_c/D_c; **MMS self-consistent passed, Biot field failed** → needs ChatGPT |
| **A — P0 Diagnostic** | ✅ Done | RHS integral, Jn, CSV; zero-mean conclusion pending re-check after RHS fix |
| **B — P1 one-sided Neumann** | ✅ First version | RHS flux correction + optional grad postprocess + MMS tests |
| **B+ — LHS grad Neumann** | ❌ Deferred | in-operator projection ineffective |
| **C — P2 corrected gradient** | ❌ Deferred | standalone enable worsens `phi_eq_res`; wait for Sprint 3 full set |
| **D/E** | ⏸ Deferred | virtual shell, A_ind, thermal coupling |

---

## 2. Sprint A (P0 Diagnostic) — Completed

See git history / 2026-06-02 entry. Core conclusions:

- RHS zero-mean projection **does not change** `phi_eq_res` / `Jn_post` → **2026-06-02 remediation**: GMRES entry previously called `setupOpheliePhiImagRhsProblem` overwriting RHS; Neumann/finalize **did not enter solve**; conclusion must be re-verified after RHS lifecycle fix.
- `Jn_post_phi_rel≈0.0038` → boundary normal current leakage is the main observable metric.

---

## 3. Sprint B (P1 one-sided Neumann) — First version completed

### 3.1 Added files

| File | Role |
|------|------|
| `electromagnetic_ophelie_phi_boundary.h` | Boundary mode enum, analytic box/cylinder normals, RHS Neumann flux correction, optional grad Neumann projection |
| `test_3d_ophelie_phi_neumann_slab/` | Quadratic φ MMS (E=0), acceptance `Jn_post` |
| `test_3d_ophelie_phi_neumann_cylinder/` | Cylinder quadratic φ MMS |

### 3.2 Added CLI

| Option | Default | Description |
|------|------|------|
| `--phi-boundary-mode=none\|one-sided-neumann\|virtual-shell-diagnostic` | `none` | literature **default unchanged** |
| `--phi-boundary-normal-source=analytic-box\|analytic-cylinder\|level-set` | `analytic-cylinder` | French uses analytic cylinder normal |
| `--phi-compatible-correction=0\|1` | `0` | **full set** G_c/D_c (French acceptance failed, default off) |
| `--phi-gradient-correction=0\|1` | `0` | **diagnostic only** (grad only, worsens eq_res) |
| `--phi-boundary-grad-neumann=0\|1` | `0` | **experimental** grad normal projection (MMS) |

Existing P0 options still valid: `--phi-rhs-project-zero-mean`, `--phi-boundary-diagnostics`, `--phi-p0-csv`.

### 3.3 Implementation notes

**RHS flux correction** (`applyOpheliePhiOneSidedNeumannRhsCorrection`):

```text
g_n = −ω n·A
rhs_i += (V_i / d_b) · σ_i · g_n     (boundary-shell particles)
d_b = max(dist_to_surface, 0.25·dp)
```

**Grad Neumann projection** (`applyOpheliePhiBoundaryGradNeumannProjection`, requires `--phi-boundary-grad-neumann=1`):

```text
grad_phi_i ← grad_phi_i − n · (n·grad_phi_i − g_n)
```

**Pipeline order**: `setupOpheliePhiImagRhsFromASrc` → `finalizeOpheliePhiImagRhsHost` (Neumann + optional projection) → GMRES → grad_phi → optional grad projection → E/J.

### 3.4 P0 case labels (CSV)

| Label | Switches |
|------|------|
| `A_baseline` | Default |
| `B_rhs_proj` | `--phi-rhs-project-zero-mean=1` |
| `C_neumann` | `--phi-boundary-mode=one-sided-neumann` |
| `D_rhs_proj_neumann` | B + C |

### 3.5 Run examples

```bash
cd ~/SPHinXsysSYCL/build

# MMS (auto-enable grad projection)
./tests/.../test_3d_ophelie_phi_neumann_slab/bin/test_3d_ophelie_phi_neumann_slab
./tests/.../test_3d_ophelie_phi_neumann_cylinder/bin/test_3d_ophelie_phi_neumann_cylinder

# French C: RHS Neumann only (does not change literature default)
./tests/.../test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced \
  --reload=1 --literature-mode --phi-boundary-mode=one-sided-neumann --state_recording=1

# French: experimental grad projection (destroys divJ; do not use for literature_passed)
./tests/.../test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced \
  --reload=1 --literature-mode --phi-boundary-mode=one-sided-neumann \
  --phi-boundary-grad-neumann=1 --state_recording=1
```

---

## 4. Acceptance record

| Date | case | phi_eq_res | Jn_post_rel | divJ_L2_red | P_raw | literature | Notes |
|------|------|------------|-------------|-------------|-------|------------|------|
| 2026-06-02 | A_baseline | 0.486 | 0.00378 | 2.06 | 50000 | ✅ | P0 |
| 2026-06-02 | B_rhs_proj | 0.486 | 0.00378 | 2.06 | 50000 | ✅ | RHS projection ineffective |
| 2026-06-07 | C_neumann (RHS only) | 0.485 | 0.00385 | 2.06 | 50000 | ✅ | **almost no improvement** |
| 2026-06-07 | C + grad_proj | 0.486 | **2e-8** | **0.40** | 50000 | ❌ | Jn↓ but divJ collapse |
| 2026-06-07 | neumann_slab MMS | — | **2e-8** | — | — | passed=1 | grad_proj=1 |
| 2026-06-07 | neumann_cylinder MMS | — | **6e-8** | — | — | passed=1 | grad_proj=1 |
| 2026-06-08 | LHS grad Neumann (slab) | — | 0.36 | — | — | passed=0 | in-operator projection ineffective |
| 2026-06-08 | MMS restored host grad projection | — | **2e-8** | — | — | passed=1 | device CK postprocess bug |
| 2026-06-08 | A_baseline | 0.486 | 0.00379 | 2.06 | 50000 | ✅ | regression |
| 2026-06-08 | C_neumann RHS | 0.486 | 0.00380 | 2.06 | 50000 | ✅ | almost no improvement |
| 2026-06-08 | C + grad_proj | 0.486 | **2e-8** | **0.40** | 50000 | ❌ | Jn↓ divJ collapse |
| 2026-06-08 | E_gradcorr | 0.685 | 0.00240 | 1.46 | 50000 | ❌ | phi_eq_res↑ |
| 2026-06-08 | C_neumann_gradcorr | 0.685 | 0.00240 | 1.46 | 50000 | ❌ | combination no benefit |
| 2026-06-02 | **Sprint1** A_baseline (after RHS lifecycle fix) | 0.486 | 0.00378 | 2.06 | 50000 | ✅ | `finalize`/`gmres_entry`/`post_solve` fingerprint **consistent** |
| 2026-06-02 | **Sprint3** F_compatible | **0.745** | **7.9e-5** | **1.34** | 50000 | ❌ | `phi_eq_res`↑; `divJ_L2_red`↓; **default off** |

### 4.7 Remediation Sprint 3 (2026-06-02, paired G_c/D_c)

**Implementation (Steps 5–6)**

1. **B matrix inversion**: `det≤0` / `|det|<1e-6` / nonfinite → `Identity` (remove negative-det weighted explosion).
2. **`ComputeOphelieVecdDivergenceCorrectedCK`**: `D_c(F)=Σ(F_i-F_j)·C_i∇W_ijV_j`, shares `C_i=B_i^{-1}` with `G_c`.
3. **Full-set switch** (`--phi-compatible-correction=1`):
   - LHS: `G_c(φ)` + `D_c(σ G_c φ)`
   - RHS: `D_c(σ A)`
   - E/J postprocess: `G_c(φ)`
   - divJ diagnostic: `D_c(J)`
4. **Diagnostics**: `phi_grad_B_stats` (det_min/max, negative/small/fallback counts).
5. **P0 label**: `F_compatible` (and Neumann/proj variants).

**Difference from `--phi-gradient-correction=1`**

| Switch | G_c | D_c | Purpose |
|------|-----|-----|------|
| `gradient-correction` | ✅ | ❌ | Diagnostic only (known to worsen eq_res) |
| `compatible-correction` | ✅ | ✅ | Full-set acceptance |

**Regression (dp=0.02, reload, literature)**

| case | phi_eq_res | Jn_post | divJ_L2_red | B det_min | literature |
|------|------------|---------|-------------|-----------|------------|
| A_baseline | 0.486 | 0.00378 | 2.06 | — | ✅ |
| F_compatible | **0.745** | 7.9e-5 | **1.34** | 0.075 | ❌ |

```text
MMS operator_linear_consistency:     passed=1  (uncorrected, default)
MMS discrete_self_consistency:       passed=1  (rel≈3.4e-7)
F_compatible pre_solve_eq_res_vol:   1.0       (LHS/RHS severely incompatible)
F_compatible rhs_div_vs_legacy_vol:  ~10.8      (D_c RHS vs legacy one order of magnitude)
phi_grad_B_stats: det_negative=0 fallback=0   (inversion healthy; issue in operator pairing/real Biot field)
```

**Conclusion**: full D_c/G_c on **MMS box** (with corrected `A=-G_c(φ)/ω`) passes linearity + discrete self-consistency; fails on **French real Biot**. **Keep default off**.

### 4.8 Sprint 3 follow-up (2026-06-02, MMS dedicated + diagnostic fix)

**Fixes**

1. **French regression bug**: `evaluateOpheliePhiCompatibleVsUncorrectedOperators` wrongly placed between `finalize` and `solve`, **overwriting RHS** → moved before `setupOpheliePhiImagRhsFromASrc` (diagnostic only).
2. **MMS `A_src` construction**: `assignManufacturedASrcFromPhiExact(DiscreteGrad)` changed to `execOphelieScalarPhiGradient` (uses `G_c` when compatible).

**Added**

| File | Role |
|------|------|
| `test_3d_ophelie_phi_compatible_operator_mms/` | uncorrected + compatible dual-mode MMS acceptance |
| `evaluateOpheliePhiCompatibleVsUncorrectedOperators` | `D` vs `D_c` RHS comparison + phi=0 eq_res |

**MMS acceptance (`test_3d_ophelie_phi_compatible_operator_mms`)**

| Mode | linearity | discrete_self_consistency_rel | passed |
|------|-----------|-------------------------------|--------|
| uncorrected | ✅ | ~3.0×10⁻⁷ | ✅ |
| compatible | ✅ | ~3.3×10⁻⁷ | ✅ |

On MMS box `rhs_unc_vs_compatible_vol≈0.22` (same field `A=-G_c(φ)/ω`, D≠D_c still 22% apart).

**French real Biot (`phi_compatible_ops`, phi=0, pre-solve)**

```text
rhs_unc_vs_compatible_vol ≈ 0.97     (D and D_c almost completely different)
rhs_cosine_unc_compatible   ≈ 0.38     (low correlation, not simple scaling)
eq_res_vol_uncorrected        ≈ 1.0      (Biot A does not satisfy uncorrected discrete equation)
eq_res_vol_compatible       ≈ 1.0      (D_c/G_c also cannot make Biot A satisfy Lφ=b)
```

| case | phi_eq_res | divJ_L2_red | literature |
|------|------------|-------------|------------|
| A_baseline (after fix) | 0.486 | 2.06 | ✅ |
| F_compatible | 0.745 | 1.34 | ❌ |

**→ Recommend pause coding and ChatGPT discussion** (see §9).

### 4.10 Continuous vector-field divergence MMS + dp scan + projection comparison (2026-06-02)

**Added**

| Component | Description |
|------|------|
| `electromagnetic_ophelie_vector_divergence_diagnostics.h` | Analytic field MMS, Biot `D(σA)` diagnostic |
| `test_3d_ophelie_vector_divergence_mms/` | 5 cases × dp scan + CSV |
| CLI | `--phi-biot-divergence-diagnostics`, `--vector-divergence-mms-scan`, `--phi-projection-operator=` |
| Bug fix | Zero-denominator explosion in scattered-field L2 metric; literature no longer overrides explicit user `--phi-*-operator` |

**dp scan (interior L2, excerpt)**

| case | dp=0.04 | dp=0.02 | dp=0.015 | Trend |
|------|---------|---------|----------|------|
| constant_x | 0 | 0 | 0 | Correct |
| linear_xyz (D_unc) | ~1.99 | ~1.99 | ~2.00 | **no convergence** |
| rotational_xy (D_unc) | ~1.4e-9 | ~1.3e-9 | ~2.0e-9 | already excellent, dp-independent |
| rotational_xy (D_c) | ~7.3e-9 | ~7.2e-9 | ~9.2e-9 | **D better than D_c** |
| biot_like (D_unc) | 0.137 | 0.122 | 0.398 | no monotonic convergence |

**Biot diagnostic (French real A, dp=0.02)**

```text
A_theta_fraction=1.0
D_unc_vs_D_c_l2≈0.97, cosine≈0.38
normalized_div_unc≈0.008, normalized_div_comp≈0.00076
```

**French projection three-way comparison (literature-mode, reload)**

| Case | LHS/RHS | phi_eq_res_vol | divJ_L2_red | literature |
|------|---------|----------------|-------------|------------|
| A_baseline | DivSigmaGrad / DivSigmaA | **0.486** | **2.06** | ✅ |
| G_edge_flux | LegacyPairwise / LegacyFlux | **9.7e-5** | **0.93** | ❌ |
| F_compatible | G_c/D_c | 0.745 | 1.34 | ❌ |

**Conclusion**

1. **Implementation bugs ruled out**: constant field, rotational D_unc, RHS fingerprint, baseline regression all normal.
2. **Method bottleneck**: `D(σA)` does not converge on linear/quadratic analytic fields; on Biot, D and D_c directions nearly orthogonal.
3. **edge-flux (existing legacy pairing)**: φ equation residual very low but **divJ physical metric worsens** → cannot simply replace baseline.
4. **production remains A_baseline**; `compatible-div-grad` stays frozen.

**Pending discussion (see §4.11)**: does edge-flux low `phi_eq_res` come from E/J still using DivSigmaGrad grad path without self-consistency?

### 4.12 Phase 0 + Sprint A/B implementation (2026-06-02)

**Phase 0 — route gating**

- `logOpheliePhiProjectionRouteWarnings()`: startup warnings for `compatible` / `edge-flux-phi-only`
- `production_literature_passed`: only `div-grad` production route = 1
- literature default `applyOpheliePhiProjectionOperatorKind(DivGrad)`

**Sprint A — vector MMS v2**

- Added `sign_alpha`, `flipped_l2`, `best_sign_l2`; CSV `ophelie_vector_divergence_mms_scan_v2.csv`
- `singular_toroidal_legacy` (/r_xy, diagnostic only) + `smooth_toroidal` (smooth toroidal field)
- **Sign conclusion**: `linear_xyz` `l2_interior≈2`, `flipped≈0.0057`, `sign_alpha≈-0.91` → prior "non-convergence" was sign diagnostic issue
- **smooth_toroidal**: `D_unc interior≈1e-7` (dp=0.04), `singular_legacy≈0.14` → non-smooth field cannot be main conclusion
- **rotational_xy**: `D_unc≈1e-9` still better than `D_c≈7e-9`; `validation_passed=1`

**Sprint B — edge-current residual**

- New file `electromagnetic_ophelie_edge_flux_diagnostics.h`; CLI `--phi-edge-flux-diagnostics=1`
- French three-way (scaled P=50kW):

| Case | edge_res_red_l2 | phi_eq_res | divJ_L2_red | production_lit |
|------|-----------------|------------|-------------|----------------|
| A_baseline | **0.024** | 0.486 | **2.06** | ✅ |
| G_edge_flux | **341** | 9.6e-5 | **0.93** | ❌ |
| F_compatible | 1.04 | 0.745 | 1.34 | ❌ |

**Decision (settled 2026-06-02)**: production stays `div-grad`; `edge-flux` self-consistent in edge-current sense (`edge_res_red≈341`) but particle divJ postprocess not companion → **freeze old phi-only edge-flux** (see §4.13 supersede).

### 4.13 Route pivot — SPH edge-flux production (2026-06-02, post-GPT confirmation)

**Background**: div-grad `phi_eq_res≈0.49` treated as formulation floor; further tuning limited. GPT discussion confirmed: edge means **SPH neighbor pairs i–j**, must use SYCL CK pairwise kernel-sum; no external mesh/CSR.

**New mainline** (`OPHELIE_EDGE_FLUX_SOLVER_DEVELOPMENT_PLAN.md`):

```text
--ophelie-current-form=edge-flux   (main switch)
  → legacy-pairwise φ LHS
  → pair electromotive drop RHS: e_ij = (φ_j-φ_i) + ω Ā_ij·(x_j-x_i)
  → JouleHeatEdge + JEdgeRecon
  → P_total_edge for 50 kW calibration
```

**Relation to §4.12**: §4.12 froze **phi-only** (`--phi-projection-operator=edge-flux` without current-form); **not** full edge-flux solver freeze.

**Sign**: production fixed legacy-compatible convention; MMS `best_sign_l2` diagnostic only. First pass `test_3d_ophelie_edge_flux_sign` (constant-A + φ=-ωA·x → edge_drop≈0).

**Stage 1 acceptance** (do not use particle `divJ_L2_red` as primary):

```text
edge_res_red > 100
P_total_edge finite
I×2 → P_total_edge×4
JouleHeatEdge outer strong, inner weak
P_edge/P_recon ∈ [0.5, 2] soft gate
```

**Stage 1 forbidden**: A_ind, thermal coupling, full J–φ monolithic, external mesh.

**Code**: `electromagnetic_ophelie_edge_flux.h`; archive index `docs/ophelie/archive_phi_divgrad_residual_floor/README.md`.

> **§4.14 supersedes Stage 1 calibration/acceptance on `P_total_edge` and `P_edge/P_recon` gate**; edge residual targets remain valid.

### 4.14 Edge-flux power fix + Stage 1 production passed (2026-06-08)

**Basis**: `OPHELIE_EDGE_FLUX_POWER_FIX_AND_NEXT_DEVELOPMENT_PLAN.md` (ChatGPT analysis + Cursor implementation).

**Root cause**: edge φ / edge residual succeeded; failure was using graph edge energy `P_graph_edge = Σ 0.25·C_ij·edge_drop²` as physical Joule power for 50 kW calibration. `C_ij` from pairwise Laplace weight can self-consistently serve φ equation but **cannot** directly be physical conductance.

**Code changes (Steps 1–4)**

| Change | File |
|------|------|
| Calibration/acceptance use `P_total_recon`; `P_graph_edge` diagnostic only | `electromagnetic_ophelie_french_literature.h` |
| `p_graph_edge` / `p_graph_over_recon` naming; `JouleHeatEdgeRecon` CK | `electromagnetic_ophelie_edge_flux.h` |
| Fields `JouleHeatEdgeRecon`, `EdgeDropAbsMax/SqMean` | `field_names.h`, `register_fields.h` |
| Pair-level edge_drop stats | `test_3d_ophelie_edge_flux_sign` |
| Uniform field power MMS | `test_3d_ophelie_edge_flux_power_uniform_field/` (new) |
| Source directory classification | `electromagnetic_ophelie/README.md` + subdirs `legacy/`, `diagnostics/`, `team7/`, … |

**New Stage 1 acceptance (edge-flux production)**

```text
edge_res_red > 100
P_recon ≈ 50 kW (calibration target)
P_recon finite; phi_res ok; fields ok
P_graph/P_recon warning only, does not fail production
```

**Unit / MMS tests**

| Test | Result |
|------|------|
| `test_3d_ophelie_edge_flux_sign` | passed=1 (+ edge_drop_linf_rel/l2_rel) |
| `test_3d_ophelie_edge_flux_power_uniform_field` case=potential_field | P_recon/P_exact=1.0; P_graph/P_exact≈12232 |
| same case=induction_field | P_recon/P_exact=1.0; P_graph/P_exact≈12232 |

**French H production (after calibration @50 kW, v2 log)**

Log: `discussion_bundle/french_H_edge_flux_production_v2.log`

| Metric | before power-fix | after power-fix |
|------|--------------|--------------|
| `P_recon` | 0.43 W | **50000 W** |
| `ampere_turns_eff` | ~1.9 | **~660** |
| `phi_eq_res_vol` | 8e-5 | 8e-5 |
| `edge_res_red` | 10920 | **10955** |
| `production_literature_passed` | 0 | **1** |
| `P_graph/P_recon` | 116402 | 116402 (warning only) |
| `divJ_L2_red` | — | **0.116** (particle grad, not edge gate) |

**Directory archive (same day)**: `electromagnetic_ophelie/` main dir keeps edge-flux + shared; legacy routes moved to `legacy/div_grad/`, `diagnostics/`, `team7/`, `french_extras/`, `stage2/`.

**Deferred / pending GPT discussion** (does not block Stage 2 prep diagnostics): see `discussion_bundle/OPHELIE_GPT_DISCUSSION_BLOCKERS.md`.

**Step 5 — H/A/G v2 comparison (2026-06-08 completed)**

Summary: `discussion_bundle/FRENCH_THREE_WAY_SUMMARY_v2.txt`

| Case | production_lit | passed | Notes |
|------|----------------|--------|------|
| H edge-flux | **1** | **1** | P_recon calibration |
| A div-grad | 1 | 1 | fallback stable |
| G phi-only | 0 | 0 | **expected** (diagnostic_only) |

**Step 6 — A_ind one-way diagnostic (2026-06-08 baseline run)**

- Test: `test_3d_ophelie_french_aind_diagnostic --reload=1` → `passed=1`
- Log: `discussion_bundle/french_aind_diagnostic_v1.log`
- Result: `A_ind/A_coil≈0.44`, `B_ind/B_coil≈0.37`; still uses **div-grad φ + particle J** (not yet `JEdgeRecon`)

**Subsequent deferred**: divJ comparison after full E/J switch to JEdgeRecon, thermal coupling / TEAM7 — see blockers B2.

**GPT discussion pack (2026-06-08)**: [`discussion_bundle/OPHELIE_GPT_ROUND_SUMMARY_2026-06-08.md`](discussion_bundle/OPHELIE_GPT_ROUND_SUMMARY_2026-06-08.md) — post-Stage 1 GPT discussion on Stage 2 direction.

### 4.15 Edge-flux Stage 2 definition / scope unification — primary E/J/Q (2026-06-08)

**Basis**: `OPHELIE_EDGE_FLUX_STAGE2_CLEANUP_AND_CURSOR_PLAN.md` Steps 1–2.

- Under edge-flux, `EImag/JImag/JouleHeat` synced to edge reconstruction (`CopyOphelieEdgeReconToPrimaryEJQCK`)
- particle-gradient default off; `--ophelie-particle-gradient-diagnostics=1` writes `*ParticleDiag`
- VTP: `JouleHeatEdgeGraph`, `PowerEdgeGraphParticle`
- H v3: `production_literature_passed=1` — `discussion_bundle/french_H_edge_flux_production_v3.log`

**Stage 2 Steps 2.5–6 (2026-06-08)**:
- `q_antisymmetry` diagnostic + gate (`q_antisym_rel_l2 < 1e-5`)
- Q spatial soft gate (outer/center, Qmax/Qmean)
- scaling regression: `test_3d_ophelie_edge_flux_scaling`
- active A: `a_coil` / `a_src` via `--ophelie-use-a-total-for-edge-flux=1`
- A_ind one-way uses `JEdgeRecon`
- umbrella header slimmed (team7/self-induction explicit include)

**Pending**: archive/ large migration, French target-power scaling sweep (12.5/50/200 kW).

### 4.16 French reload full-case acceptance + I² scaling (2026-06-08)

**Particles**: user SYCL relax → `build/reload/Reload.xml` (`n_glass=20500`), slightly different from v2 `n_glass=20625`; conclusions consistent.

**H edge-flux @50 kW** (`french_H_edge_flux_production_v4.log`):

| Metric | Value |
|------|-----|
| production_literature_passed | **1** |
| P_recon_edge | 49999.9 W |
| edge_res_red | 9348 |
| q_antisym_rel_l2 | 7.13e-8 |
| Q_outer/center | 14.14 |
| max_EImag / max_JImag | 336 / 5379 |

**French fixed-I scaling** (`--no-literature-calibrate`, `discussion_bundle/french_H_edge_flux_scaling_reload.log`):

| coil-current-scale | P_recon (W) | P/P_ref |
|--------------------|-------------|---------|
| 1.0 | 7.27889 | 1.00 |
| 0.5 | 1.81972 | **0.250** |
| 2.0 | 29.1156 | **4.00** |

`edge_flux_scaling_passed=1`; B3 (blockers) **closed**.

**Same-batch three-way comparison** (user reload):

| Case | production_lit | Key difference |
|------|----------------|----------|
| H edge-flux | **1** | edge mainline + P_recon |
| A div-grad | **1** | phi_eq_res≈0.427, divJ_L2_red≈2.34 |
| G phi-only | **0** | divJ_L2_red≈0.93 (expected fail) |

Summary: `discussion_bundle/FRENCH_THREE_WAY_SUMMARY_v3.txt`

### 4.17 Complex edge-flux Stage 3.0a / 3.1 (2026-06-02)

**Basis**: `OPHELIE_COMPLEX_EDGE_FLUX_SOLVER_PLAN.md`; user chose **Option A** explicit field names + `--ophelie-edge-flux-complex=0|1` (default 0).

**Stage 3.0a (completed)**

- `test_3d_ophelie_french_aind_diagnostic` CLI wired to French/Ophelie filters
- `runFrenchReducedAIndOneWayDiagnostic` edge-flux branch adds `setupOpheliePhiImagRhsFromASrc` + `finalize`
- `test_3d_ophelie_edge_flux_sign`: imag/real dual-chain edge-drop sign tests both `passed=1`
- A_ind phase audit: `A_ind_real=0`, `A_ind_imag/A_src_real≈0.43` → needs complex edge-flux

**Stage 3.1 (completed)**

- Fields: `e/j_edge_recon_{real,imag}`, `phi_{real,imag}`, `edge_flux_residual_{real,imag}`, `joule_heat_edge_recon_{real,imag,complex}`
- `OphelieEdgeFluxComponent` + `makeOphelieEdgeFlux{Imag,Real}Component`
- CK parameterization: RHS / residual / recon / complex Joule heat (`Q=0.5σ(|E_r|²+|E_i|²)`)
- `electromagnetic_ophelie_phi_component.h`: componentized phi RHS/solve/LHS
- French pipeline: imag GMRES → real PCG → complex Joule heat when `edge_flux_complex_=1`

**Regression (build/reload, GPU)**

| Test | complex=0 | complex=1 |
|------|-----------|-----------|
| `test_3d_ophelie_edge_flux_sign` | passed=1 | passed=1 |
| `test_3d_ophelie_french_reduced` @50kW | P_recon=49999.9 W, prod_lit=1 | P_recon=49999.9 W, prod_lit=1 |
| `test_3d_ophelie_french_aind_diagnostic` | — | A_ind_imag/A_src≈0.43, passed=1 |

**real chain under complex=1**: `A_imag≈0` → real RHS≈0, PCG residual≈0 (expected); power still dominated by imag chain.

**Stage 3.2 / 3.5 (2026-06-02, continued)**

- `evaluateOphelieEdgeFluxQAntisymmetryForComponent` + `evaluateOphelieEdgeFluxResidualForComponent`
- Q antisymmetry CK parameterized (`phi` / `active_a` / `a_sign`); `ji` direction passes sign explicitly for `a_sign<0`
- `solvePhiComponentGMRESWithCurrentRhs`: real chain GMRES (imag still `solvePhiImagGMRESWithCurrentRhs`)
- French complex=1: imag/real dual-chain q_antisym logs; real chain `phi_gmres_real`
- A_ind one-way: complex mode solves real chain + `JEdgeReconReal/Imag` into Biot–Savart

**Regression**

| Test | Result |
|------|------|
| `test_3d_ophelie_edge_flux_sign` | imag+real q_antisym≈9.5e-8, `passed=1` |
| `test_3d_ophelie_french_reduced` complex=1 | P_recon=49999.9 W, real q_antisym=0, `prod_lit=1` |
| `test_3d_ophelie_french_aind_diagnostic` complex=1 | `passed=1`, dual-chain J source |

**Stage 3.4 / 3.5 / 3.6 infrastructure (2026-06-02, code ready, cases pending)**

- `solveOphelieComplexEdgeFluxWithCurrentA` / `execOphelieComplexEdgeFluxSolveReconAndPower`
- A_ind one-way feedback: `--ophelie-aind-one-way-feedback=0|1` (default 1); after `K[J]`, `A_total` re-solves complex edge-flux
- `OphelieAIndOneWayDiagnostic` extended: `P_complex_coil_only`, `P_complex_total_A`, `edge_res_red_{imag,real}`
- `runOphelieComplexEdgeFluxSelfInductionWithPhiSolve`: complex edge-flux Picard (auto when `edge_flux_complex=1`)
- `test_3d_ophelie_edge_flux_scaling`: added `complex` I² scaling suite (manufactured φ + `P_complex`)

**First run data (reload, 2026-06-02)** — see [`discussion_bundle/OPHELIE_COMPLEX_EDGE_FLUX_RUN_RECORD.md`](discussion_bundle/OPHELIE_COMPLEX_EDGE_FLUX_RUN_RECORD.md)

| Test | Result | Description |
|------|------|------|
| aind_diagnostic complex=1 | **passed=1** | P_coil=7.28W → P_total=8.64W; max_J_real=21.9 |
| edge_flux_scaling | I² OK / exit=1 | P(0.5I)/P(1I)=0.25, P(2I)/P(1I)=4; target_P gate N/A for MMS |
| self_induction_picard | false passed=1 | CLI not wired; actually imag-only Picard; **fixed** |

### 4.18 Test fixes (2026-06-02)

1. **`test_3d_ophelie_french_self_induction_picard`**: wired French/Ophelie CLI filters; removed unconditional `applyFrenchLiteratureMode`; outputs `complex_picard_path`, `max_J_{real,imag}`; complex path gates `phi_eq_res_vol<0.01`.
2. **`test_3d_ophelie_edge_flux_scaling`**: `passed` gates I² ratio only; `target_P` changed to diagnostic (French absolute calibration via `french_reduced --no-literature-calibrate`).

**Rerun (after fix, 2026-06-02)** — see [`discussion_bundle/OPHELIE_COMPLEX_EDGE_FLUX_RUN_RECORD.md`](discussion_bundle/OPHELIE_COMPLEX_EDGE_FLUX_RUN_RECORD.md) §3

| Test | Result |
|------|------|
| `test_3d_ophelie_edge_flux_scaling` | `edge_flux_scaling_passed=1` |
| `test_3d_ophelie_french_self_induction_picard` complex=1 | `passed=1`, `self_induction_complex` + `phi_pcg_real`, 7 iter, J_rel=0.039 |
| `test_3d_ophelie_french_reduced` complex=1 @50kW | P_scaled=50000 W, `demo_passed=1` |
| `test_3d_ophelie_french_aind_diagnostic` complex=1 | regression `passed=1` |

**Minor fix (§4.18 cont.)**: `ophelieEdgeResidualReductionRatio` — when level0≈0, `edge_res_red=-1` (n/a), avoids real chain false zero report.

### 4.19 Picard joint convergence + Stage 3.7 power test (2026-06-02)

**Picard joint gating (Stage 3.6 wrap-up)**

- Parameter: `self_induction_phi_eq_res_tolerance_` (default 0.01); CLI `--self-induction-phi-tol=`
- Helper: `ophelieSelfInductionPicardConverged(j_rel, phi_eq_res_vol, params)` — complex uses `phi_tol=0.01`, legacy uses `phi_eq_res_vol_gate_=0.65`
- Outer loop: `picard_converged=1` and break only when **J_rel < j_tol and phi_eq_res < phi_tol**
- Complex path records `phi_eq_res_vol=max(imag, real)` each step; legacy applies `applyOpheliePhiImagLhsOperator` + `hostPhiEqResVolFromCurrentLhsRhs` each step
- `OphelieFrenchSelfInductionPicardResult::picard_converged`; test `passed` requires `converged=1` (joint)

**Rerun (build/reload)**

| Test | Result |
|------|------|
| `test_3d_ophelie_french_self_induction_picard` complex=1 | `picard_converged=1`, 7 iter, J_rel=0.039, phi_eq_res=3.4e-4, `passed=1` |
| `test_3d_ophelie_edge_flux_power_uniform_field` | 5 cases (H v4×2 + complex A/B/C), `summary passed=1` |

**Stage 3.7 / Plan §8.2 (completed)**: `test_3d_ophelie_edge_flux_power_uniform_field` extended complex uniform-field power closure:
- `complex_potential_real`: PhiReal=-E0·x → P=0.5σ|E0|²V
- `complex_induction_imag_chain`: AReal=A0 → P=0.5σω²|A0|²V
- `complex_induction_real_chain`: AImag=A0 → P=0.5σω²|A0|²V

### 4.20 Picard parameter sweep + thermal one-way (2026-06-02)

**Picard convergence sweep** — `test_3d_ophelie_french_self_induction_picard_sweep`
- Grid: relax ∈ {0.15, 0.3, 0.45} × max_iter ∈ {6, 8, 12} (reload, complex=1)
- Gate: reference (0.3, 8) converges + all "sufficient budget" grid points converge (relax≥0.3 and iter≥8; or relax<0.3 and iter≥12)
- Diagnostic: low relax + few iter expected non-convergence (e.g. 0.15/6 → J_rel≈0.14)

**Joule → heat one-way (Stage 4 start, no feedback)**
- Header: `diagnostics/electromagnetic_ophelie_joule_to_heat_one_way.h`
- Explicit Euler: `OphelieThermalDeltaT += Q·dt/(ρ·cp)`, `Temperature = T₀ + OphelieThermalDeltaT` (see §4.22)
- `syncOphelieJouleHeatPrimaryForThermalOneWay`: complex Q → `JouleHeat` (VTP / primary field)
- Heat source field: `JouleHeatEdgeReconComplex` (device-authoritative; **no upload** host Q before thermal step)
- `test_3d_ophelie_complex_joule_to_heat_one_way`: uniform induction field → Q_complex → one temperature step, `passed=1`
- Original `test_3d_ophelie_joule_to_heat_one_way` uses shared helper

### 4.21 Stage 4 — French reload thermal one-way end-to-end (2026-06-02, late fix closed)

**Pipeline** `runFrenchReducedEmThenJouleHeatOneWay` (`electromagnetic_ophelie_joule_to_heat_one_way.h`)
- EM: `runFrenchReducedEmPipeline` (complex edge-flux) → `syncOphelieJouleHeatPrimaryForThermalOneWay`
- Thermal: `ApplyOphelieJouleHeatOneWayTemperatureStepCK` (**all device**: freeze Q, `ReduceDynamicsCK` closure diagnostics)
- Multi-step: `--thermal-steps` × `--thermal-dt` (default 3×1 s)

**New test** `test_3d_ophelie_french_complex_joule_to_heat_one_way` (`--reload=1 --ophelie-edge-flux-complex=1`)

| Metric | Before fix (reload) | After fix (reload) |
|------|------------------|------------------|
| `P_joule_W` | 7.28 | 7.28 |
| `phi_eq_res_vol` | ~3.2e-4 | ~3.2e-4 |
| `vol_weighted_rel_err` | **~0.12** | **0** |
| `energy_balance_rel_err` | **~0.12** | **0** |
| `energy_vs_power_rel_err` | **~0.12–0.18** | **~1.7e-7** |
| `E_joule_J` / `E_thermal_J` | 21.84 / **19.12** | **21.8365 / 21.8365** |
| `closure_mismatch_vol_frac` | ~0.66 | **0** |
| `passed` | 0 (5% gate) | **1** |

**Root cause and conclusion (§4.22)**: ~12% gap **not** EM/Q distribution error or SYCL double-buffer main cause, but **IEEE-754 precision loss from accumulating ΔT~10⁻⁵ directly on T₀=300 K** (~11799/20500 particle temperature rises swallowed by rounding). Fixed with independent `OphelieThermalDeltaT` accumulating from zero; **issue closed**.

### 4.22 Thermal one-way energy closure: floating-point root cause + all-device path (2026-06-02)

**Observed (French reload, n=20500, 3×1 s)**
- EM normal: `P_joule_W≈7.28`, `phi_eq_res_vol≈3.2e-4`; ParaView Q spatial distribution reasonable (outer high, inner low).
- Thermal closure failed: `E_joule_J≈21.84`, `E_thermal_J≈19.12` → **relative gap ~12.4%**; `thermal_max_rel_err≈2`; `closure_mismatch_vol_frac≈0.66`.
- MMS small case (uniform Q≈5e4, `ΔT≈2e-3`) always `passed=1`, masking French large-case issue.

**Investigation (ruled out)**
1. ~~Particle/Vol/reload~~: confirmed 20500 particles, uniform Vol.
2. ~~SYCL host direct T write~~: all-device CK + device reduce still had gap → not main cause.
3. ~~`syncVariableToDevice(Q)` overwriting device Q~~: removed Q upload, gap remained → not main cause.
4. **Host control experiment**:
   - Pull Q/T from device, host `T=300` then `T+=Q·dt/(ρ·cp)` → reproduced ~12–18% gap; ~**11799 particles Q>0 but ΔT≈0**.
   - Same Q, **accumulate from T=0** → `energy_balance_rel_err=0`, `E_joule=E_thermal=21.8365 J` **exact closure**.

**Root cause**
- Under French non-uniform Q, many particles have `ΔT ~ 10⁻⁵–10⁻⁴`, too small vs `T₀=300`; `double` gives `300 + 1e-5 → 300` (increment lost).
- Joule integral `∫Q·Vol·dt` still counts full energy; thermal storage `∫ρ·cp·ΔT·Vol` under-counts → apparent **~12% energy non-conservation**.
- Ratio **1050/1200=0.875** matches observed ~12.5% relative error (not material parameter mismatch; ρ·cp cancels in closure formula).

**Fix (`electromagnetic_ophelie_joule_to_heat_one_way.h`)**
1. New field **`OphelieThermalDeltaT`**: host/device accumulate from **0**: `ΔT += Q·dt/(ρ·cp)`.
2. **`Temperature = T₀ + OphelieThermalDeltaT`** (output/ParaView still shows physical temperature).
3. **Closure diagnostics** (`ReduceDynamicsCK`) read `OphelieThermalDeltaT`, not `T−T₀`.
4. **Q stays device-authoritative**: after EM only `ensureVariableDelegatedOnDevice(Q)`; **forbid** `syncVariableToDevice(Q)`.
5. Register sync `OphelieThermalDeltaT` (0) and `Temperature` (T₀) to device.

**Acceptance (2026-06-02 build run)**

```text
test_3d_ophelie_joule_to_heat_one_way:          max_rel_err=0           passed=1
test_3d_ophelie_complex_joule_to_heat_one_way:  max_rel_err≈1.2e-7      passed=1
test_3d_ophelie_french_complex_joule_to_heat_one_way:
  vol_weighted_rel_err=0  energy_balance_rel_err=0  closure_mismatch_vol_frac=0
  E_joule_J=E_thermal_J=21.8365  energy_vs_power_rel_err≈1.7e-7  passed=1
```

**Next steps (Stage 4.1)**: add real thermal diffusion + crucible BC on closed foundation; still forbid Picard thermal feedback. ChatGPT discussion **not needed now** (Q distribution and power already reasonable).

### 4.23 Stage 4.1 — Isotropic thermal diffusion + cold crucible Dirichlet (2026-06-02)

**Implementation** `diagnostics/electromagnetic_ophelie_thermal_diffusion_one_way.h`

- Equation (explicit, one-way): `OphelieThermalDeltaT += Q·dt/(ρ·cp) − (k/(ρ·cp))·L_pair(T)·dt`, `T = T₀ + OphelieThermalDeltaT`
- `L_pair(T)`: reuses `OpheliePairwiseLaplaceCK` (same family as φ operator discretization)
- Cold crucible BC: French cylindrical shell layer `OphelieThermalBoundaryMask` → Dirichlet `T = T₀` (width `phi_boundary_distance_factor_ · dp`)
- Material defaults: `ρ=2500`, `cp=1200`, `k=1 W/(m·K)` (Jacoutot order of magnitude)

**New test** `test_3d_ophelie_thermal_diffusion_mms`: uniform Q + box cold walls, `boundary_compliance=1`, `passed=1`

**French e2e extension** `--thermal-diffusion=1` on `test_3d_ophelie_french_complex_joule_to_heat_one_way`

| Mode | Typical result (reload) | passed |
|------|-------------------|--------|
| No diffusion | `E_joule=E_thermal=21.84 J`, closure 0 | 1 |
| With diffusion | `E_thermal≈11.85 J < E_joule` (cold-wall loss), `boundary_compliance=1` | 1 |

**Note**: diffusion mode **does not require** pointwise Joule closure (heat redistributed/lost at boundary); gate changed to `E_thermal ≤ E_joule` + boundary compliance.

**Pending (Stage 4.2+)**: natural convection coupling, σ(T) feedback, literature full-field transient quantitative comparison.

### 4.24 Stage 4.2 — Literature thermal material properties + Temperature VTP + sub-step recording (2026-06-02)

**New files**

| File | Role |
|------|------|
| `diagnostics/electromagnetic_ophelie_french_thermal_material.h` | Jacoutot 2008 Table 1 @ 1473 K and reduced prototype presets |
| `diagnostics/electromagnetic_ophelie_thermal_vtp.h` | device→host sync then write VTP; optional record every N thermal steps |

**Jacoutot Table 1 (@ 1473 K) vs reduced default**

| Quantity | literature preset | reduced (regression default) |
|----|-------------------|---------------------|
| ρ | 2750 kg/m³ | 2500 |
| cp | 1150 J/(kg·K) | 1200 |
| k | 4 W/(m·K) | 1 |
| T₀ | 1473 K | 300 K |

**CLI** (`test_3d_ophelie_french_complex_joule_to_heat_one_way`)

- `--thermal-material=reduced|literature` or `--use-literature-thermal=1`
- Overrides: `--thermal-t0=`, `--rho=`, `--cp=`, `--k=`
- Longer integration: `--thermal-steps=N`, `--thermal-dt=`
- VTP: `--thermal-state-recording=1`, `--thermal-record-interval=N` (0 = final step only)

**VTP fields**: `Temperature`, `OphelieThermalDeltaT`, `JouleHeat`, `JouleHeatEdgeReconComplex`; diffusion mode also `OphelieThermalLaplaceT`, `OphelieThermalConductivity`, `OphelieThermalBoundaryMask`.

**Note**: regression gates still use **reduced + T₀=300**; literature preset for ParaView long runs and future literature comparison. At T₀=1473, `OphelieThermalDeltaT` accumulation path still valid with larger relative rise.

**Typical command**

```bash
./test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload=1 --ophelie-edge-flux-complex=1 --thermal-diffusion=1 \
  --use-literature-thermal=1 --thermal-steps=30 --thermal-state-recording=1 --thermal-record-interval=5
```

**Pending (Stage 4.2+)**: natural convection, σ(T), Jacoutot full-field transient quantitative comparison.

### 4.25 EM side: A_ind Picard independent gate + dp convergence scan (2026-06-02)

**Existing (Stage 3.6–3.7, §4.19–§4.20)**

| Test | Role |
|------|------|
| `test_3d_ophelie_french_aind_diagnostic` | one-way `A_ind=K[J]`, no feedback |
| `test_3d_ophelie_french_self_induction_picard` | Picard joint gate: `J_rel` + `phi_eq_res_vol` (**not** in `literature_passed`) |
| `test_3d_ophelie_french_self_induction_picard_sweep` | relax × max_iter sweep |

**reload + complex edge-flux reference**

```bash
test_3d_ophelie_french_self_induction_picard --reload=1 \
  --ophelie-current-form=edge-flux --ophelie-edge-flux-complex=1
# → self_induction_complex, ~7 iter, picard_converged=1, passed=1
```

**Added** `test_3d_ophelie_french_em_dp_scan` (lattice, decoupled from thermal one-way)

- `--em-dp-scan`: dp ∈ {0.08, 0.06, 0.04}, complex edge-flux EM; gate **finest** `phi_eq_res_vol < 0.01` (lattice may not be monotonic; quality benchmark uses reload @ dp=0.02)
- `--em-dp-scan-picard`: same + Picard; finest dp must `picard_converged`
- CSV: `output/ophelie_french_em_dp_scan.csv`

**CLI minor fix**: when `--ophelie-edge-flux-complex=1` and `--ophelie-current-form=` not explicit, auto-set `edge-flux` (avoid accidental particle-gradient Picard).

**Relation to thermal / natural convection**: this line is EM + optional A_ind Picard only; **no** fluid, σ(T) feedback, or Jacoutot convection comparison.

### 4.26 TEAM7 L2 complex edge-flux (native STL, dp=3 mm, 2026-06-02)

**Test**: `test_3d_ophelie_team7_complex_edge_flux`  
**Dedicated log**: `tests/.../test_3d_ophelie_team7_complex_edge_flux/TEAM7_VALIDATION_LOG.md`  
**Run record**: `output/team7_validation_history.txt` (appended each run)

| Milestone | Conclusion |
|--------|------|
| L2 NaN @ dp=3 mm | edge-flux double accumulation + RHS normalization `max(L2,max)` + fine grid `pair_weight` |
| feedback secondary solve | default **off**; skip real φ when imag dominates; else `bz_rms_total` collapse |
| Coil Bz | volume-racetrack ~33% high; **`--team7-coil-source-scale=0.754`** default; peak 7.81/7.81 mT |
| phase0 total field | ≈ coil field (`J_real≈0`); phase90 induction in `B_ind_imag`, **not gated** |
| after restore | recompute `A_ind/B_ind`; probe Biot excludes fallback; skin-layer z diagnostic |

**Typical L2 numbers** (scale=0.754): `bz_rms_coil≈0.32`, `bz_rms_total_phase90≈11.7`, `bz_rms_ind_imag_skin≈5.75`, `P_recon≈40 W`, `Bind/Bcoil≈3.9`, `passed=1`.

**x-segment phase90 (2026-06-02 evening)**: `x_peak_ref` RMS_p90≈6.3, `x_high_bias`(144–216)≈12.5; depth profile shows uniform `J_imag` ~2×10⁴ A/m² through thickness, `Bind/Bcoil` ~4 at each layer.

**imag-chain audit + J post-calibration (2026-06-02)**: `cos(J,E)≈1`, `P_vol≈80 W`; `--team7-ind-j-post-scale=auto` (s≈1/`Bind/Bcoil`) brings `phase90_RMS` 11.7→**1.68**, skin layer 5.75→**1.30** — confirms **J/B_ind amplitude ~4× high** is main phase90 cause; fix needed in restore/edge-recon chain, not post-hoc scaling.

**`imag_a_sign` + power audit (2026-06-02 cont.)**: `--ophelie-edge-flux-imag-a-sign=-1` flips induction B **sign** (phase90 peak −17.9→+17.9 mT), `phase90_RMS` 11.7→**9.8**, skin layer 5.75→**3.9**; `Bind/Bcoil` unchanged. `p_graph/p_recon`: pre-restore ~10⁷ (misleading) → post-restore **~646** (graph `Σ¼c_ijΔ²` vs vol `J·E` recon ~40 W).

**Pending**: edge-flux power definition closure (~650×); phasor default sign; racetrack peak x offset. **Recommend GPT discussion** (see TEAM7_VALIDATION_LOG §4).

### 4.11 Decision points needing discussion (2026-06-02, partially closed by §4.13–§4.16)

1. **Primary acceptance metric**: keep `divJ_L2_red + P_raw` (current plan §10), or also require `phi_eq_res` drop?
2. **edge-flux route**: `phi_eq_res→0` but `divJ_L2_red<1` — continue investing in unified edge-flux + E/J grad operator set, or abandon edge-flux?
3. **phi_eq_res saturation ~0.49**: with vector MMS showing D_unc excellent on rotational fields, accept as **range projection error of Biot A under DivSigmaGrad/DivSigmaA pairing** (not GMRES bug)?

### 4.9 Topics for recommendation and ChatGPT discussion (2026-06-02)

1. **MMS self-consistent but Biot failed**: `D_c/G_c` closes on manufactured `A=-G(φ)/ω`, but on **multiloop Biot `A_src`** (non-gradient, strong curl) `D_c(σA)` and `D(σA)` cosine only ~0.38. Does this mean **compatible correction applies only to near-potential fields**, not induction vector potential?
2. **baseline `phi_eq_res≈0.49` with `eq_res(phi=0)=1` on Biot**: shift acceptance from "equation residual" to "divJ improvement + power closure", no longer pursue strong φ equation consistency on real fields?
3. **Next discretization direction** (pick one):
   - keep uncorrected `D/G`, change **RHS/LHS sign/legacy-flux pairing** or **pairwise Laplace** route;
   - **Helmholtz/Hodge projection** of Biot `A` before `D_c`;
   - **kernel support/boundary ghost** instead of B-matrix;
   - accept `phi_eq_res` as diagnostic, optimize **divJ_L2_red** and **Jn** tradeoff.

### 4.6 Remediation Sprint 1 (2026-06-02, ChatGPT plan Steps 1–4)

**Engineering bug fixes**

1. **CLI prefix offset**: `--phi-gmres-*`, `--phi-eq-res-gate` use `strlen` helper; added `logOphelieFinalParams()`.
2. **RHS lifecycle (Option B)**: `solvePhiImagGMRES/PCG/Jacobi` split to `*WithCurrentRhs` (zero φ + solve only); French / Neumann MMS use `solvePhiImagWithCurrentRhs`; RHS after `finalize` **no longer overwritten inside solve**.

**Added diagnostics**

| File | Role |
|------|------|
| `electromagnetic_ophelie_phi_rhs_diagnostics.h` | RHS fingerprint (sum/l2/xor) |
| `electromagnetic_ophelie_phi_operator_diagnostics.h` | LHS linearity + discrete self-consistency MMS |
| `test_3d_ophelie_phi_operator_linear_consistency/` | repeat/add/scale (float tolerance 1e-6) |
| `test_3d_ophelie_phi_discrete_self_consistency/` | `A=-G_h(phi)/ω` → `‖L_h(phi)-b_h‖` |

**Acceptance**

```text
operator_linear_consistency:  passed=1  (repeat~3e-8, add/scale~2e-7)
discrete_self_consistency:    passed=1  (rel≈3.3e-7)
French A_baseline:            phi_eq_res_vol=0.486, literature_passed=1
RHS fingerprint:              finalize == gmres_entry == post_solve
```

**Conclusion**: RHS enters GMRES, but `phi_eq_res≈0.49` **unchanged** → bottleneck in discrete operator LHS/RHS inconsistency (Sprint 3 paired D/G), not RHS overwrite bug.

### 4.2 Device boundary projection root cause (2026-06-08, fixed)

1. **SYCL `finalizeLoadIn` no-op before variable delegate** → `setupOpheliePhiBoundaryParticleFields` end forces `DelegatedData(par_device)` + re-upload.
2. **`prepareForOutput` on host-written boundary fields overwrites host with device zeros** → `FromFields` projection only `prepareForOutput(grad_phi)`; mask/normal/gn read host directly.
3. **Device CK in-place `grad_phi` path still unreliable** → unified **host geom projection** (pull `grad` from device → project → push to device). MMS `Jn~10⁻⁸` restored.

### 4.3 LHS grad Neumann conclusion (2026-06-08)

- `--phi-boundary-lhs-grad-neumann=1` on MMS slab **worsens residual** (`res_neumann > res_no_neumann`), `Jn_post` still ~0.36.
- Double-counting with RHS flux correction avoided (skip RHS when LHS on), no improvement.
- **Device CK postprocess projection** (`applyOpheliePhiBoundaryGradNeumannProjectionDynamics`) **ineffective**; restored **host projection** → MMS `Jn~10⁻⁸`.
- LHS path kept in CLI/code for diagnostic, **default off**, not main route.

### 4.4 P2 corrected gradient (2026-06-08)

- New file `electromagnetic_ophelie_phi_gradient.h`: B matrix + linear corrected gradient.
- CLI: `--phi-gradient-correction=0|1` (default off).
- Wired into `DivSigmaGrad` LHS and French E/J postprocess grad path.
- French `E_gradcorr` / `C_neumann_gradcorr`: `phi_eq_res` **0.685** (baseline 0.486), `literature_passed=0`; `Jn` slightly better (0.0024) but **not worth default on**.
- May need **RHS DivSigmaA also corrected grad** for LHS/RHS consistency — pending theory discussion.

### 4.5 Issues for recommendation and ChatGPT discussion

1. **φ equation residual saturation ~0.49**: with `phi_gauge_penalty=0`, constant σ, both RHS/LHS Div form — mainly **discrete Neumann compatibility / null space** issue, not insufficient GMRES iterations?
2. **Grad postprocess Neumann**: MMS `Jn→10⁻⁸`, but French real Biot **divJ_L2_red 2.0→0.4** — does postprocess fix **normal only** while breaking **divergence and Laplace consistency**? How to embed ghost in **grad operator itself** not E/J postprocess?
3. **P2 corrected gradient**: standalone B-matrix worsens `phi_eq_res` because **LHS DivSigmaGrad and RHS DivSigmaA gradient discretizations still inconsistent**? Which operators need paired correction?
4. **French boundary `|n·σA|~10⁻¹⁰`**: RHS Neumann flux ~10⁻⁸ — does this mean **induction field naturally satisfies approximate Neumann on boundary**, and real bottleneck is **kernel-truncated grad/div consistency**?

### 4.1 Key findings (2026-06-07)

1. **French real Biot boundary `|n·σA|≈10⁻¹⁰`** → RHS Neumann correction ~10⁻⁸, **almost no effect on phi_eq_res / Jn**.
2. **Grad postprocess projection** drops MMS `Jn_post` to ~10⁻⁸, but French **divJ_L2_red 2.0→0.4** (literature failed) → default **off**, MMS/diagnostic only.
3. Next priority: **embed Neumann ghost in DivSigmaGrad LHS** (in-operator grad, not postprocess), or **P2 corrected gradient**.

---

## 5. Deferred / Optional (retry later)

### 5.0 Archived (div-grad residual-floor route, see archive README)

- LHS ghost DivSigmaGrad, grad correction default on, French D group — **not Stage 1 investment**
- Focus **edge-flux Stage 1** (§4.13)

### 5.1 High priority (Sprint C+, archived, diagnostic preserved only)

- [ ] **LHS ghost in DivSigmaGrad** (active during solve, not grad postprocess) — archived
- [ ] `--phi-gradient-correction=0|1` (P2) — archived
- [ ] LevelSet normal (replace analytic cylinder)
- [ ] French **D group**: B + C + corrected gradient — archived

### 5.2 Medium/low priority

- [ ] `--phi-boundary-mode=virtual-shell-diagnostic` (P3, cylinder only)
- [ ] Particle fields `PhiBoundaryFlag` / `PhiBoundaryNormal` to VTP
- [x] Real thermal diffusion + crucible Dirichlet BC (Stage 4.1 §4.23)
- [x] Literature thermal material preset + Temperature VTP + sub-step thermal recording (Stage 4.2 §4.24)
- [ ] Literature full-field transient comparison / natural convection (Stage 4.2+)
- [x] A_ind Picard independent gate (§4.19–§4.20; `test_3d_ophelie_french_self_induction_picard` + `_sweep`)
- [x] dp convergence scan (§4.25; `test_3d_ophelie_french_em_dp_scan` lattice EM / Picard)
- [ ] device_krylov performance
- [ ] literature default Neumann (needs LHS ghost effective first)
- [ ] complex dummy shell (**not doing**)

### 5.3 Tried but not main route for now

| Method | French Jn | French divJ | French eq_res | Conclusion |
|------|-----------|-------------|---------------|------|
| RHS zero-mean | — | — | no change | not main cause |
| RHS Neumann | no change | no change | no change | ineffective when n·A≈0 |
| Grad postprocess Neumann | **large↓** | **collapse** | no change | MMS acceptance only |

---

## 6. Changed file manifest (Sprint A + B + Remediation Sprint 1)

| File | Sprint |
|------|--------|
| `electromagnetic_ophelie_phi_rhs_diagnostics.h` | Remediation 1 new |
| `electromagnetic_ophelie_phi_operator_diagnostics.h` | Remediation 1 new |
| `electromagnetic_ophelie_phi_gmres.h` | Remediation 1 `*WithCurrentRhs` |
| `electromagnetic_ophelie_phi.h` | Remediation 1 Jacobi `*WithCurrentRhs` |
| `electromagnetic_ophelie_cli.h` | A+B+Remediation 1 |
| `electromagnetic_ophelie_french_literature.h` | A+B+Remediation 1 |
| `test_3d_ophelie_phi_operator_linear_consistency/` | Remediation 1 new |
| `test_3d_ophelie_phi_discrete_self_consistency/` | Remediation 1 new |
| `electromagnetic_ophelie_phi_gradient.h` | C+Remediation3 G_c/D_c, inversion fix |
| `electromagnetic_ophelie_diagnostics.h` | Remediation3 divJ passes params |
| `electromagnetic_ophelie_parameters.h` | Remediation3 `phi_compatible_correction_` |
| `test_3d_ophelie_phi_compatible_operator_mms/` | Sprint3+ dual-mode MMS |
| `electromagnetic_ophelie_phi_mms_helpers.h` | Sprint3+ `A_src` via `execOphelieScalarPhiGradient` |
| `electromagnetic_ophelie_phi_boundary.h` | B new, B+ fix |
| `electromagnetic_ophelie_phi_boundary_diagnostics.h` | A new, B extended |
| `electromagnetic_ophelie_parameters.h` | A+B |
| `electromagnetic_ophelie_cli.h` | A+B |
| `electromagnetic_ophelie_french_literature.h` | A+B |
| `electromagnetic_ophelie.h` | A+B |
| `electromagnetic_ophelie_progress.h` | A fix |
| `test_3d_ophelie_french_reduced.cpp` | A |
| `test_3d_ophelie_phi_neumann_slab/` | B new |
| `test_3d_ophelie_phi_neumann_cylinder/` | B new |
| `OPHELIE_PHI_DEVELOPMENT_LOG.md` | this file |
| `electromagnetic_ophelie/README.md` | source directory classification (2026-06-08) |
| `electromagnetic_ophelie_edge_flux.h` | §4.14 power fix + JouleHeatEdgeRecon |
| `electromagnetic_ophelie_thermal_diffusion_one_way.h` | §4.23 thermal diffusion + cold-wall BC |
| `electromagnetic_ophelie_french_thermal_material.h` | §4.24 Jacoutot Table 1 preset |
| `electromagnetic_ophelie_thermal_vtp.h` | §4.24 Temperature VTP + sub-step recording |
| `test_3d_ophelie_french_em_dp_scan/` | §4.25 lattice dp scan + optional Picard |
| `test_3d_ophelie_thermal_diffusion_mms/` | §4.23 MMS diffusion + Dirichlet |
| `test_3d_ophelie_french_complex_joule_to_heat_one_way/` | §4.21–§4.24 French thermal one-way e2e |
| `test_3d_ophelie_team7_complex_edge_flux/` | §4.26 TEAM7 L2 edge-flux + validation log |
| `.../TEAM7_VALIDATION_LOG.md` | §4.26 dedicated milestone record |
| `test_3d_ophelie_complex_joule_to_heat_one_way/` | §4.20 complex MMS thermal one-way |
| `test_3d_ophelie_joule_to_heat_one_way/` | §4.20 uniform Q MMS |
| `test_3d_ophelie_edge_flux_power_uniform_field/` | §4.14 new |
| `discussion_bundle/OPHELIE_GPT_DISCUSSION_BLOCKERS.md` | pending GPT discussion summary |

---

## 7. Decision record

1. **literature default**: `phi_boundary_mode=none`, `phi_boundary_grad_neumann_projection=false`.
2. **French C (RHS only)**: runnable, does not break acceptance, but **no clear accuracy gain yet**; keep CLI for regression comparison.
3. **Grad projection**: only in MMS tests with `phi_boundary_grad_neumann_projection_=true`; French requires explicit `--phi-boundary-grad-neumann=1`.
4. **production default**: `--ophelie-current-form=particle-gradient` + `--phi-projection-operator=div-grad` (A_baseline fallback).
5. **edge-flux phi-only (2026-06-02)**: `phi_eq_res≈1e-4` but particle `divJ_L2_red≈0.93` — E/J postprocess not companion; **diagnostic only** (§4.12).
6. **edge-flux production (2026-06-02)**: `--ophelie-current-form=edge-flux` starts Stage 1 — see §4.13 and `OPHELIE_EDGE_FLUX_SOLVER_DEVELOPMENT_PLAN.md`.
7. **div-grad baseline**: preserved fallback/regression, no longer primary push to lower `phi_eq_res`.
8. **edge-flux calibration (2026-06-08)**: 50 kW calibration and production gate use **`P_total_recon`** (edge-reconstructed E/J); `P_graph_edge` diagnostic only — see §4.14.
9. **P_graph/P_recon≈1e5**: known; `C_ij` not physical conductance; does not block production — pending GPT discussion on future graph normalization MMS.
10. **Thermal one-way energy closure (2026-06-02)**: French ~12% gap from **small ΔT floating-point loss on T₀=300**, not EM physics error; fixed with **`OphelieThermalDeltaT` accumulate from zero** + all-device path; French/MMS tests `passed=1`.

---

## 8. Related documents

- [`OPHELIE_COMPLEX_EDGE_FLUX_SOLVER_PLAN.md`](OPHELIE_COMPLEX_EDGE_FLUX_SOLVER_PLAN.md) — GPT complex edge-flux baseline plan
- [`discussion_bundle/OPHELIE_GPT_DISCUSSION_SINCE_COMPLEX_EDGE_FLUX_PLAN.md`](discussion_bundle/OPHELIE_GPT_DISCUSSION_SINCE_COMPLEX_EDGE_FLUX_PLAN.md) — **GPT discussion pack** (full summary + issue manifest since plan)
- [`OPHELIE_PHI_DISCRETIZATION_NEXT_SPRINT_PLAN.md`](OPHELIE_PHI_DISCRETIZATION_NEXT_SPRINT_PLAN.md)
- [`OPHELIE_PHI_BOUNDARY_AND_NEXT_DEVELOPMENT_PLAN.md`](OPHELIE_PHI_BOUNDARY_AND_NEXT_DEVELOPMENT_PLAN.md)
- [`OPHELIE_PHI_OPERATOR_AUDIT.md`](OPHELIE_PHI_OPERATOR_AUDIT.md)
