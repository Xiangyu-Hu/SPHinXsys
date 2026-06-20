# TEAM7 complex edge-flux validation log

> Dedicated log: `test_3d_ophelie_team7_complex_edge_flux`
> Automatically added after each run: `./output/team7_validation_history.txt`
> Overview is also written to [`OPHELIE_PHI_DEVELOPMENT_LOG.md`](../../../../../../OPHELIE_PHI_DEVELOPMENT_LOG.md) §4.26

Last updated: **2026-06-02** (P5.6-lite finishing effect evaluation; P6a/P6b; P5-fix)

---

## 1. Calculation examples and levels

| item | value |
|------|-----|
| Geometry | TEAM7 native STL reload (`CoilSourceBody` + `PlateBody`) |
| coil source | volume-racetrack (default) or `--coil-source-model=filament-racetrack` (P1) |
| Plate solution | complex edge-flux + A_ind one-way |
| Reference | `TEAM7_Bz_A1_B1_reference_mT.csv` (50 Hz phase0 / phase90) |
| typical dp | 3 mm (`--native-dp-mm=3`) |

| Level | CLI | Access Control |
|-------|-----|------|
| L1 `coil-only` | `--team7-level=coil-only` | `bz_rms_coil` < 0.70 |
| L2 `one-way` | `--team7-level=one-way` (default) | L1 coil Bz + EM finiteness (`max_J`, φ residual) |
| L3 `picard` | `--team7-level=picard` | L3 total field Bz threshold |

**L2 Description**: `bz_rms_total` (phase0) and `phase90` **Log only**, no hard access control (the induction field is still being calibrated).

---

## 2. Milestone Timeline

### 2026-06-02 — L2 NaN fix (dp=3 mm)

- **Phenomena**: NaN appears in `particle_generation_TEAM7` / L2 one-way.
- **Root cause**: The edge-flux pairwise conductance `c_ij` accumulates overflow in the float neighborhood under fine meshes (it has nothing to do with French `dp=20 mm`).
- **FIX**:
- `ReconstructOphelieEdgeFluxElectricCurrentCK`: double accumulation + Tikhonov ε ∝ `condition_proxy²`
  - `applyOphelieEdgeFluxInputNormalization`：`scale = max(scale_l2, scale_max)`（p99.5 \|rhs\|）
- TEAM7 `pair_weight_regularization_` Zoom in by `(dp_ratio)²` relative to French dp

### 2026-06-02 — feedback Quadratic solution crash

- **Phenomena**: `bz_rms_total≈10.3` (probe `B_ind_real` explodes) when `aind_one_way_feedback_resolve_=1`.
- **Root cause**: `A_total` lower real part φ produces false `j_real` on imaginary part dominant J.
- **FIX**:
- L2 **off by default** feedback (`--ophelie-aind-one-way-feedback=1` optional)
- feedback path: **skip real phi** when imaginary part dominates; optional `A_total` RHS additional normalization
- **Conclusion**: True one-way = φ(A_coil) → J → A_ind, **should not** do A_total quadratic solution by default.

### 2026-06-02 — Coil Bz calibration (volume-racetrack vs filament reference)

- **Phenomenon**: peak Bz ≈ 10.36 mT for `source_scale=1` vs reference 7.81 mT (+33%).
- **Handling**: `--team7-coil-source-scale=`, L2 **Default 0.754** (≈ 7.81/10.36).
- **Audit**: `ophelieTeam7CoilPathAmpereTurnsAuditPassed` checked against `NI × source_scale`.
- **Results** (scale=0.754): `bz_rms_coil≈0.318`, peak 7.81/7.81 mT, `passed=1`.

### 2026-06-02 — phase0 / phase90 Decomposition of cognition

- First round one-way: `max_J_real=0`, induced current in **J_imag** → probe `B_ind_real=0`.
- Therefore **phase0 total field ≈ coil field** (`bz_rms_total ≈ bz_rms_coil`) is expected, non-inductive missing.
- **phase90** corresponds to `B_ind_imag` (coil `j_src` is real only, coil imag B≈0).
- x=288 mm: sim `+1.31 mT` vs ref `+1.24 mT` (close); x=126–216 mm still has a large deviation.

### 2026-06-02 — x segment phase90 + plate depth profile + VTP

- **x segment** (`team7_probe_phase90_xbands_<level>.csv`):

| band | x [mm] | n | RMS coil phase0 | RMS phase90 full | RMS phase90 skin |
|------|--------|---|-----------------|------------------|------------------|
| x_left | 0–95 | 6 | 0.73 | 8.78 | 4.43 |
| x_under_coil | 95–293 | 11 | 0.29 | 11.81 | 5.81 |
| x_peak_ref | 108–144 | 3 | **0.12** | **6.28** | **3.26** |
| x_high_bias | 144–216 | 5 | 0.35 | **12.54** | 6.15 |
| x_right | 216–288 | 5 | 0.31 | 12.07 | 5.90 |

- **Conclusion**: phase90 near the reference peak (108–144 mm) is the best but still too large; **144–216 mm** (coil B overshoot area) phase90 is the worst → consistent with the coil peak offset (sim peak ~216 mm).
- **Depth Profile** (`team7_plate_depth_profile_<level>.csv`): `J_imag_L2_vol` **Approximately uniform across layers** (~2.0–2.5×10⁴) at z≈2.7–17.7 mm, non-apparent surface skin; `Bind/Bcoil` by layer **3.5–4.9**. This shows that the phase90 deviation of the probe is not only a "deep J leak", but the **J amplitude is overall high**.
- **VTP**: `PlateBody_ite_0000000000.vtp` includes `B_ind_*`, `J_edge_recon_imag`, `Team7DepthBinZ_mm` (ParaView shaded by depth).

### 2026-06-02 — Code side root cause troubleshooting (restore / J source / scan parameters)

**restore chain**: `j_linearity=1.0`, `j_edge/j_imag=1`; `safe_rhs_l2` from 1e4→2e6 (barely normalized) **J, Bind/Bcoil, phase90 completely unchanged** → **Exclude input normalize/restore amplification J**.

**J source comparison @ x=162 mm**:

| quantity | value |
|----|-----|
| `\|J_edge\|_vol` vs `\|J_em\|_vol` (σ(−∇φ−ωA)) | 53408 vs 53569, **ratio 1.003** |
| `Bz_ind` Biot(J_edge) | **−12.86 mT** |
| `Bz_ind` Biot(J_em) | **−13.60 mT** |
| Reference phase90 | **+1.415 mT** |
| After negation `Bz_ind` | +13.6 mT (**~9.6× reference amplitude**) |

→ edge-recon and EM constitutive are **self-consistent**; the problem lies in the magnitude/sign of **the entire φ→J→Biot induction chain relative to the TEAM7 reference**, not the restore or J_edge calculation alone.

**Other signals**:
- `p_graph/p_recon ≈ 2.3×10⁷` (graph theory edge-drop power vs `J·E` reconstruction power is severely disconnected)
- `e_edge_em_mismatch ≈ 0.14`（14%）
- `ind_imag_particle_RMS=1` is an illusion that the **air probe SPH interpolation is always 0** and cannot be accepted.

**Scan parameters**: `pair_weight` 0.01→0.444 only slightly moves phase90 (11.69→11.69), not the main reason.

**Biot Geometric/Symbolic Variants** (phase90 RMS vs reference):

| Variants | RMS | Description |
|------|-----|------|
| `volume_full` | 11.69 | Default full board volume score |
| `volume_neg` | **9.81** | **B_ind_imag negation** (sign) |
| `skin` | 5.75 | Top of board z≥z_top−2dp |
| `sheet_top` | 13.59 | Source point raised to top of plate z (worse) |
| `zmid_band` | 4.56 | z_mid±dp layer only |

**Combined (diagnostic)**: `--team7-ind-j-post-scale=auto` + negated `volume_neg` RMS≈**0.41**, peak sim≈**+1.17 mT** vs ref≈+1.42 mT (magnitude ~82%).
→ **Sign + ~1/Bind/Bcoil Amplitude** Two lines can explain most of the phase90 deviation; it is still physically necessary to compress the amplitude from the φ/J chain, and cannot rely on post-scale for a long time.

### 2026-06-02 — `imag_a_sign` Physical symbol experiment + power audit after restore

**CLI**: `--ophelie-edge-flux-imag-a-sign=` (default `1`; real part chain fixed `a_sign=-1`).

| Indicators | `a_sign=+1` | `a_sign=-1` |
|------|-------------|-------------|
| `bz_rms_total_phase90` | 11.69 | **9.81** |
| `ind_imag_skin` RMS | 5.75 | **3.88** |
| `zmid_band` RMS | 4.56 | **2.68** |
| phase90 peak sim/ref | **−17.87** / 1.415 | **+17.87** / 1.415 |
| `Bind/Bcoil` / `J_imag` | 3.91 / 5.3e4 | **unchanged** (restore linear normalization) |
| `e_edge_em_mismatch` | 0.14 | **2.37** (variation) |

- **Conclusion**: `a_sign=-1` is equivalent to the **inversion** of `B_ind_imag` in the probe (`volume_full(+1)` ↔ `volume_neg(+1)` swap), **correct the phase90 sign**, slightly reduce the RMS, but **do not suppress the amplitude** (still ~9× reference); and the edge–EMF constitutive consistency becomes worse → **do not change the default for now**, and will be determined after the phasor agreement is discussed with GPT.
- **Combination still valid**: `a_sign=-1` + `--team7-ind-j-post-scale=auto` can further compress the amplitude (for diagnostic purposes).

**`p_graph/p_recon` Root cause (partial clarification)**:

| stage | `p_graph` | `p_recon` | ratio |
|------|-----------|-----------|------|
| restore **before** (old audit) | ~26 kW | ~1 mW | **~2.3×10⁷** (misleading) |
| restore **after** (new audit) | ~26 kW | ~40 W | **~650** |

→ 10⁷ Mainly from **audit before restore**, recon is still in coil-only scale; after restore, recon≈40 W is consistent with `P_recon_W`. graph power (`Σ ¼ c_ij Δ²`) and vol integral `J·E` recon **still differ by ~650×**, which belongs to **two power definitions/normalization paths** and is not closed, not a simple restore bug. Power auditing has been moved to printing after **restore**.

### 2026-06-02 — imag-chain power audit + J post-calibration experiment

**imag-chain audit** (after restore, before probe):

| Quantity | Typical value |
|----|--------|
| `\|J_imag\|_vol` | ~5.34×10⁴ |
| `\|E_imag\|_vol` | ~1.51×10⁻³ |
| `Re(J·E)_vol` / `P_vol_W` | ~80.6 W |
| `cos(J,E)` | ~1 (σ·E=J holds in the sense of vol) |
| `P_recon_W` (edge-recon) | ~40 W (about half P_vol, complex power agreement difference) |

**Diagnostic CLI**: `--team7-ind-j-post-scale=auto` scales `J/E` to `Bind/Bcoil≈1` and recalculates `A_ind` (**does not change L2 gate**, only probe/profile comparison).

| Metrics | Baseline | `auto` (s≈0.256) |
|------|------|-------------------|
| `bz_rms_total_phase90` | 11.69 | **1.68** |
| `bz_rms_ind_imag_skin` | 5.75 | **1.30** |
| `x_peak_ref` phase90 | 6.28 | **1.34** |
| `x_high_bias` phase90 | 12.54 | **1.74** |
| `\|B\|_RMS` | 1.46 | **0.30** |

- **Conclusion**: The main reason for the phase90 deviation is that the **induced current amplitude (→ B_ind) is generally higher ~4×**, not just the skin layer or the shape of the x segment; uniform scaling J can greatly bring the reference closer, but the **peak position/sign** (sim peak≈−1.17 mT vs ref≈+1.42 mT) still needs to be clarified by the phasor agreement.
- **Next step (physics)**: Find root causes on edge-flux normalization/restore chains or Tikhonov scaling to avoid post hoc `post-scale`.

### 2026-06-02 — Probe and logging enhancements

- **Recalculate A_ind/B_ind** after restore (to avoid scaling chain drift)
- Probe Biot(J) **Exclude** `edge_recon_fallback` particles
- Probe **Plate Top Skin Layer** Biot (z ≥ z_plate_top − 2dp) Diagnostic Depth Leak
- Logs: `coil_x_span_RMS`, `|B|_RMS`, `edge_flux_scale`, `J_imag_L2_vol`, `P_recon_W`
- Output: `team7_validation_history.txt` appends running records

### 2026-06-02 — P0 verification caliber frozen (GPT decision)

- **phase90**: only `team7_phase90_probe_passed` (threshold 0.5 mT abs RMS), **no hard access control**; default `B_phase90 = −B_ind_imag` (`--team7-phase90-convention=minus-imag`).
- **L2 one-way**: always `diagnostic_only=1`, `team7_validation_passed=0`; `source_scale≠1` also diagnostic.
- **exit code**: based on `smoke_passed`, not `team7_validation_passed`.
- **Output tag**: `f50_vol_asign_p1_p90_minus_src0.754_jraw` and other automatic naming; CSV/history with tag.

### 2026-06-02 — P1 L1 source model comparison (`source_scale=1`, dp=3 mm)

Command (from `build/`):

```bash
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7_complex_edge_flux/bin/test_3d_ophelie_team7_complex_edge_flux \
  --reload=1 --team7-level=coil-only --native-dp-mm=3 \
  --team7-coil-source-scale=1.0 \
  --team7-reference-dir=$HOME/SPHinXsysSYCL/tests/extra_source_and_tests/3d_examples/reference_data/team7 \
--coil-source-model=volume-racetrack # or filament-racetrack
```

| source model | RMS_rel | RMS_abs [mT] | coil_x_span | peak sim/ref [mT] | peak_x sim/ref [mm] | best_fit_scale_peak |
|--------|---------|--------------|-------------|-------------------|---------------------|---------------------|
| volume-racetrack | 0.630 | 2.946 | 0.629 | 10.36 / 7.81 | 216 / 126 | **0.754** |
| filament-racetrack | **0.527** | **2.464** | **0.499** | 9.85 / 7.81 | 198 / 126 | 0.793 |

- **Conclusion**: filament is closer to the TEAM7 filament reference** than volume ** (RMS −16%, peak position 198 vs 216 mm), but still **~26% overshoot** with `source_scale=1` (best_fit ≈ 0.75–0.79). Volume's 0.754 is **consistent** with volume best_fit 0.754, indicating that the L2 default scale mainly comes from the geometric difference between volume vs filament reference, rather than NI audit errors.
- **Peak Position**: Both model sim peaks are skewed to the right (198–216 mm vs ref 126 mm) → Coil field spatial distribution remains to be audited by P3 (not just amplitude).
- **Next step**: ~~L2 one-way use `filament-racetrack` + `source_scale=1` to find out `Bind/Bcoil`~~ (already run, see below); then open P2 Picard.

### 2026-06-12 — P3 L1 source strict caliber (document frozen)

**Important**: The current reference CSV for L1 `coil-only` comparison is the **phase0 total field probe** (`TEAM7_Bz_A1_B1_reference_mT.csv`), **not** the f=0 pure source field reference.

Logs/access must specify:

```text
phase0 total reference is not a pure source-only reference.
```

meaning:

| Comparison item | Whether the current state can be strictly closed | Description |
|--------|------------------|------|
| coil-only vs ref phase0 RMS | **Part** | RMS mixes source field + environment when ref includes plate/medium effects |
| peak fit / best_fit_scale | diagnostic use | estimable geometric scale, not source-only acceptance |
| L2 hard access control | **Not used** phase0 total Impersonates L1 source pass | See §P0 Freeze |

**Follow-up preparations required** (P3 document item, non-code blocking):

- f=0 source-only B reference, or
- no-plate / σ=0 trusted FEM source field reference

Before the L1 source is fully closed, phase90 / Bind/B should not be treated as a hard validation gate (consistent with the Round 3 decision).

**L2 one-way + filament + `source_scale=1`**（2026-06-02）：

| Indicators | filament s=1 | volume s=0.754 (historical) |
|------|--------------|------------------------|
| `bz_rms_coil` | 0.527 | ~0.32 |
| `Bind/Bcoil` | **3.99** | ~3.91 |
| `phase90_RMS` | 11.81 | ~11.7 |
| `phase90_rms_abs_mT` | 10.45 | ~8.7 |
| `J_imag_L2_vol` | 6.4e4 | ~5.3e4 |

→ Replace filament and remove source_scale **not improved** `Bind/Bcoil≈4`; induction chain amplitude problem and coil source geometric scale **basic decoupling**, root cause still in edge-flux / φ→J→Biot chain (observed in P3 or Picard iteration).

### 2026-06-02 — P2 L3 Picard bottom line (filament + `source_scale=1`)

Added `runTeam7ComplexEdgeFluxPicardWithLog`: record `J_rel_raw/relaxed`, `Bind/Bcoil`, `Bz_phase90_RMS`, `P_recon`, etc. in each round, and output `team7_picard_<tag>.csv`.

```bash
./test_3d_ophelie_team7_complex_edge_flux \
  --reload=1 --team7-level=picard --native-dp-mm=3 --ophelie-edge-flux-complex=1 \
  --team7-coil-source-scale=1.0 --coil-source-model=filament-racetrack \
  --self-induction-max-iter=8 \
  --team7-reference-dir=$HOME/SPHinXsysSYCL/tests/extra_source_and_tests/3d_examples/reference_data/team7
```

| iter | J_rel_relaxed | Bind/Bcoil | Bz_phase90_RMS | P_recon [W] |
|------|---------------|------------|----------------|-------------|
| 1（≈L2 one-way） | 1.00 | **3.99** | 11.81 | 1240 |
| 2 | 0.86 | 4.89 | 7.01 | 1493 |
| 3 | 1.69 | 5.66 | **1.94** | 1776 |
| 4 | 0.66 | 6.19 | 9.90 | 2134 |
| 5 | 1.84 | 6.52 | 14.66 | 2638 |
| 6 | 0.81 | 6.80 | 13.18 | 3371 |
| 7 (final) | **2.12** | **7.26** | 5.58 | 4402 |

- **Not converged** (`picard_converged=0`, `J_rel_final=2.12`); `Bind/Bcoil` diverges from 3.99 ** to 7.3**, not the "Picard pressure amplitude" expected by GPT.
- phase90 RMS briefly dropped to 1.94 in iter 3, but worsened as `Bind/Bcoil` increased → **Unstable physical solution**, cannot be regarded as acceptance.
- Final state: `bz_rms_total_phase90=5.58`, `phase90_rms_abs_mT=4.93`, `smoke_passed=0` (Picard not converged).
- **Conclusion**: The current Picard (`A_total` + relaxation=default) is **oscillating/diverging** on the TEAM7 high guide board; you need to adjust `self_induction_relaxation_factor`, convergence criteria or change the coupling strategy before running. P3 edge-flux scale audit priority increased.

### 2026-06-02 — P3 multi-frequency Bz + Jey diagnosis (GPT solution)

**New CLI**:
- `--team7-frequency=50|200` → update `params.frequency_` with reference CSV column (50 Hz: col 2–3; 200 Hz: col 4–5)
- `--team7-probe-line=a1-b1|a2-b2` → A1–B1 (y=72,z=34) or A2–B2 (y=144,z=34)
- Jey: `TEAM7_Jey_A1_surface_reference_Am2.csv` (**not yet in stock**) → automatically degrades to **sim-only** probe line output

**L2 one-way + filament + `source_scale=1` frequency comparison** (A1–B1):

| f [Hz] | `bz_rms_coil` | `bz_rms_total` | `phase90_RMS` | `phase90_rms_abs_mT` | `Bind/Bcoil` | `J_imag_L2_vol` |
|--------|---------------|----------------|---------------|----------------------|--------------|-----------------|
| 50 | 0.527 | 0.527 | 11.81 | 10.45 | **3.99** | 6.4×10⁴ |
| 200 | 0.658 | 0.658 | **147.5** | **44.8** | **15.9** | 2.6×10⁵ |

- `Bind/Bcoil≈16` at 200 Hz (non-physical ~4× linear expectation at 50 Hz) → **ω scale/edge-flux chain goes out of control at high frequencies**, qualitatively compared with skin depth (δ_200≈6 mm) and needs to be compared with the depth profile.
- Depth profile: 50 Hz `J_imag` is relatively uniform between layers; 200 Hz absolute value ~4× but **still approximately uniform** between layers** (not obvious surface skin) → the problem is still **overall J amplitude scaling**, not pure deep leakage.
- **Jey** (z≈19 mm probe line, sim-only): `Jy_surface/mid_phase90≈0.76` (50/200 Hz is the same, because one-way only has `J_imag` chain); **J itself vs Biot coupling** can be determined after COMSOL table2 CSV is put into the library.
- Output: `team7_jey_probe_<tag>.csv` (sim Jy edge x); fixed excitation log: `TEAM7 validation uses fixed coil excitation, not target-power calibration`.

### 2026-06-02 — P3c edge-flux ω scale audit (GPT scheme, **complete**)

**New capabilities**:
- `OphelieAIndOneWayDiagnostic` / `Team7OneWayEdgeFluxSummary` extensions: `rhs_l2_pre_norm`, `grad_phi_imag_vol`, `omega_a_coil_vol`, `e_imag_vol`, `e_edge_em_mismatch`, `j_restore_linearity`, `p_graph_over_recon`
- CLI: `--team7-omega-sweep=1` → Sweep **25 / 50 / 100 / 200 Hz** (filament + `source_scale=1`), write `team7_omega_scaling_<tag>.csv`
- Console printing `Bind/B(f)/Bind/B(50)`, `J(f)/J(50)` vs `f/50`, `(f/50)²`

**Command** (from `build/`):

```bash
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7_complex_edge_flux/bin/test_3d_ophelie_team7_complex_edge_flux \
  --reload=1 --team7-level=one-way --native-dp-mm=3 --ophelie-edge-flux-complex=1 \
  --coil-source-model=filament-racetrack --team7-coil-source-scale=1.0 \
  --team7-frequency=50 --team7-omega-sweep=1 \
  --team7-reference-dir=$HOME/SPHinXsysSYCL/tests/extra_source_and_tests/3d_examples/reference_data/team7
```

**Sweep results** (L2 one-way, filament, `source_scale=1`, σ_glass fixed):

| f [Hz] | δ_skin [mm] | `rhs_l2_pre` | `input_scale` | `ωA/∇φ` | `J_imag_vol` | `Bind/B` | `P_recon` [W] | `e_edge_em_mis` |
|--------|-------------|--------------|---------------|---------|--------------|----------|---------------|-----------------|
| 25 | 16.9 | **0**† | 0.00825 | 1.958 | 3.20×10⁴ | **1.99** | 14.5 | 0.144 |
| 50 | 12.0 | 1.0×10⁴ | 0.00413 | 1.958 | 6.40×10⁴ | **3.99** | 57.8 | 0.144 |
| 100 | 8.46 | 1.0×10⁴ | 0.00206 | 1.958 | 1.28×10⁵ | **7.97** | 231 | 0.144 |
| 200 | 5.98 | 1.0×10⁴ | 0.00103 | 1.958 | 2.56×10⁵ | **15.9** | 925 | 0.144 |

†At 25 Hz `applyOphelieEdgeFluxInputNormalization` reports `rhs_l2=0` (RHS volume L2 value 0 before normalization), low frequency boundary/structure pending; 50–200 Hz behavior consistent.

**Scaling Law (relative to 50 Hz)**:

| Quantity Ratio | 100 Hz / 50 Hz | 200 Hz / 50 Hz | Theory `f` Linear | Theory `f²` |
|------|----------------|----------------|---------------|-----------|
| `Bind/B` | **2.0** | **4.0** | 2 / 4 | 4 / 16 |
| `J_imag_vol` | **2.0** | **4.0** | 2 / 4 | — |
| `P_recon` | **4.0** | **16.0** | — | 4 / 16 |
| `J_imag/ω` | ~204 (constant) | ~204 | — | — |

**Physical Conclusion**:
1. **200 Hz `Bind/B≈16` is not "high frequency out of control"**, but the entire imag chain is approximately linear amplification of **ω** (`Bind/B ∝ f`, `P_recon ∝ f²`), which is consistent with the two points of P3.
2. **`ωA_coil/∇φ_imag ≈ 1.96` Constant at full frequency** → The ratio of the induced electromotive force term and the gradient term in edge-flux RHS is stable, and the non-frequency dependent phase/sign is disordered.
3. **`e_edge_em_mismatch ≈ 14.4%` Constant over frequency** → edge EMF has a systematic deviation from volume `E_imag`, but is independent of ω.
4. **`j_restore_linearity ≈ 1`** → restore `J` is a linear scaling, non-linear amplification root cause.
5. **The root cause is still at the 50 Hz baseline**: The absolute value of `Bind/B≈4` is too large (relative to the TEAM7 reference induction field), and the ratio is amplified in proportion to ω; changing the frequency cannot "wash out" the ~4× deviation.

**Output**: `output/team7_omega_scaling_one-way_f50_fil_asign_p1_p90_minus_src1.000_jraw.csv`

### 2026-06-02 — P3b Jey Reference + P4 Baseline J/B Separation (GPT Scheme)

**Jey Reference Inventory** (50 Hz):
- File: `reference_data/team7/TEAM7_Jey_A1_surface_reference_Am2.csv`
- Source: Fujiwara & Nakata 1990 / NGSolve [TEAM-7/team7.ipynb](https://github.com/NGSolve/TEAM-problems) Measured phasor `Jy_A4B4_50Hz_ref` (top of plate z=19 mm, y=72 mm), ×10⁶ → A/m²
- 200 Hz column is temporarily 0; can be completed with `reference_data/team7/tools/import_team7_jey_comsol_table2.py` from COMSOL `multiturn_coil_asymmetric_conductor_table2.txt`

**New Diagnosis**:
- `team7_j_vs_b_probe_split_<tag>.csv`: Compares `Jy_sim/Jy_ref` to `B_ind_sim/B_ref` along x (phase90/minus-imag convention)
- `team7JeyReferenceHasFrequencyData`: 200 Hz skip RMS when reference is missing

**L2 one-way filament `source_scale=1`，50 Hz**：

| Indicator | Value | Description |
|------|-----|------|
| `Jy_phase90_RMS` | **4.88** | sim only imag chain; phase0 RMS≈1 (J_real≈0) |
| `median \|J_sim/J_ref\|` | **~9.8** | Probe line Jy is overall too large ~10× |
| `median \|B_ind_sim/B_ref\|` | **~11.1** | phase90 induction B (minus-imag) |
| `@x=162mm` `Jy_sim/ref` | 1.94e6 / −1.49e5 → **~13×** | |
| `@x=162mm` `B_ind_sim/ref` | 15.33 / 1.415 mT → **~10.8×** | |
| `(B/J)_sim` vs `(B/J)_ref` | **7.91 vs 9.50** | Biot coupling in J is slightly lower than reference |

**in conclusion**:
1. **J and B are of the same magnitude (~10×)**, and `J_sim/J_ref` is slightly larger than `B_ind_sim/B_ref` → The main reason is that the J amplitude recovered by the **edge-flux chain is too large** and is not independently amplified by the Biot probe.
2. `(B/J)_sim < (B/J)_ref` → If J is reduced to the reference, the **unit J coupling from Biot to the probe is slightly weaker** (it is not inconsistent with the volume ratio of `Bind/B_coil`: the latter is the intra-board field ratio, and the former is the air probe line).
3. Consistent with P3c: **50 Hz baseline ~4× `Bind/B` + ~10× Jy** Need to return to edge-flux RHS/restore/one-way missing mask instead of adjusting the Biot core alone.

**Output**: `team7_jey_probe_*`, `team7_j_vs_b_probe_split_*`

### 2026-06-02 — P5 A_total feedback comparison (GPT: one-way missing mask)

**Settings**: filament, `source_scale=1`, 50 Hz; compare `--ophelie-aind-one-way-feedback=0` (default) vs `=1`.

| Indicators | feedback=0 | feedback=1 |
|------|------------|------------|
| `Bind/Bcoil` (after restore) | 3.987 | **3.987** (same) |
| `phase90_RMS` | 11.81 | **11.81** |
| `J_imag_L2_vol` | 63993 | **63993** |
| `Jy_phase90_RMS` | 4.88 | **4.88** |
| `max_J_imag` (after feedback, before restore) | 7.32×10⁶ | **3.02×10⁴** |
| `P_total` / `P_coil_only` | 0 / 0.000985 | **0.000985 / 0.000985** |

- **skip real phi** (imag dominant) when feedback=1; `edge_res_red_imag≈0.88`.
- **restore linear descaling** (`j_restore_linearity≈1`) pushes down the feedback `J` and pulls it back → **Probe/Bind/B indicators are completely consistent with pure one-way**.
- **Conclusion**: A single quadratic solution of **A_total in the current pipeline cannot be used as a shielding repair**; if you want to do shielding diagnosis, you should record **before restore** `Bind/B` or change the restore strategy (only diagnose branch).

### 2026-06-02 — P2b Picard Relaxation Sweep (filament, `source_scale=1`, 50 Hz, 6 rounds)

| `self-induction-relax` | Final state `Bind/B` | Final state `phase90_RMS` | `J_imag_L2_vol` | Remarks |
|------------------------|---------------|-------------------|-----------------|------|
| — (L2 one-way) | 3.99 | 11.81 | 63993 | Baseline |
| **0.05** | 4.32 | **7.11** | 36638 | phase90 **best**; `J_rel` still 0.30; not converged |
| 0.15 | 6.52 | 14.66 | 94327 | Oscillation |
| 0.35 | **42.6** | **112** | 1.09×10⁶ | Divergence |

- phase90 goes from ~12→~7 (iter 5–6) with **relax=0.05**, but `Bind/B` is still ~4.3, **not closed**.
- Default relax=0.15 (historical bottom line) and larger relaxation both **worse**.
- CSV: `team7_picard_picard_f50_fil_pic_r0.05.csv` etc.

**Stage conclusion (recommended to align with GPT)**:
1. **The main reason is still the J amplitude recovered by edge-flux** (P3b J/B separation), not Biot independent amplification.
2. Neither **A_total feedback** nor **Mild Picard** can pull `Bind/B` to ~1; feedback is also offset by restore.
3. Next steps should be to discuss: **edge-flux RHS/conductance scale**, whether **restore should be involved in acceptance**, or **TEAM7-specific A-φ coupling** (non-French scale direct migration).

### 2026-06-11 — P0 normalization/restore audit (after GPT decision, **first version completed**)

**New CLI**:
- `--ophelie-edge-flux-normalization-mode=off|field-scale-restore|solver-local`
- `--team7-normalization-sweep=1` (scan `safe_rhs_l2` = 1e4/1e5/1e6/1e12, force `field-scale-restore`)

**Acceptance Logic Fix**:
- `team7_validation_passed` now includes `team7_phase90_probe_passed` (strict L3)
- `field-scale-restore` + Picard/feedback → `diagnostic_only`
- `solver-local` (P0 partial implementation: no scale/restore physics) → `diagnostic_only`

**Command** (filament, `source_scale=1`, 50 Hz):

```bash
./test_3d_ophelie_team7_complex_edge_flux \
  --reload=1 --team7-level=one-way --native-dp-mm=3 --ophelie-edge-flux-complex=1 \
  --coil-source-model=filament-racetrack --team7-coil-source-scale=1.0 \
  --team7-normalization-sweep=1 \
  --team7-reference-dir=.../reference_data/team7
```

**Sweep results (post-restore sensitivity to `safe_rhs_l2`)**:

| `safe_rhs_l2` | `input_scale` | `J_post` | `Bind/B_post` | `P_post` [W] | `phase90_RMS_post` | `inv_err_J` |
|---------------|---------------|----------|---------------|--------------|-------------------|-------------|
| 1e4 | 0.00413 | 63993 | 3.987 | 57.8 | 11.81 | ~0 |
| 1e5 | 0.0413 | 63993 | 3.987 | 57.8 | 11.81 | ~0 |
| 1e6 | 0.413 | 63993 | 3.987 | 57.8 | 11.81 | ~0 |
| 1e12 | 1.0 | 63993 | 3.987 | 57.8 | 11.81 | 0 |

**in conclusion**:
1. **post-restore acceptance quantity is not sensitive to `safe_rhs_l2`** → one-way linear restore is at least self-consistent; `Bind/B≈4` / `phase90≈12` are not accidental products of normalization parameters.
2. **restore still covers the changes in J before feedback** (the conclusion of P5 remains unchanged); Picard/feedback continues to be marked `diagnostic_only` under `field-scale-restore`.
3. **`solver-local` completely implements P2**; the next step is **P1a** high-σ uniform/harmonic-A σ sweep.

**Output**: `team7_normalization_sweep_<tag>.csv`

### 2026-06-11 — P1a high-σ edge-flux benchmark (**Complete v3**)

> **⚠️ v3 conclusion superseded**: see below **P1a v4 + P1 rotational-A** (2026-06-11).
> `harmonic_A_imag_phi_solve failed` is **misjudgment**; constant A is pure gauge, and J≈0 after full φ-solve is gauge cancellation.

**New test**: `test_3d_ophelie_high_sigma_edge_flux_scaling` (box dp=40 mm, n=512)

**Three groups of use cases** (σ ∈ {16, 1e3, 1e5, 1e7, 3.526e7} S/m; f ∈ {50, 200} Hz):

| case | path | result (v3, deprecated caliber) |
|------|------|------|
| `uniform_E_real` | Default φ=−E₀·x, `execOphelieEdgeFluxPostPhiPipeline` (do not resolute φ) | **20/20 passes** (E/J rel ~1e-7, P≈1) |
| `harmonic_A_imag_edge_only` | Default A_coil_real=A₀, φ=0, post-pipeline only | **20/20 passed** |
| `harmonic_A_imag_phi_solve` | normalize → full φ solve → restore | ~~0/10 failed~~ → **v4 renamed gauge_cancellation, 20/20 pass** |

### 2026-06-11 — P1a v4 caliber correction + P1 rotational-A (**Complete, user rerun confirmation**)

**CSV**：`discussion_bundle/team7_l2_outputs/high_sigma_edge_flux_scaling_P1a_v4.csv`  
**Summary**: `rows=60 pass_rows=60 passed=1`

| case | acceptance criteria | result |
|------|----------|------|
| `constant_A_gauge_cancellation` | `gauge_pass=1`, `J_phi/J_edge≈3.6e−5` | **20/20 pass** (J≈0 is **expected**) |
| `constant_A_edge_only` | edge reconstruction diagnosis, non-production solution | **20/20 pass** |
| `rotational_A_uniform_B_phi_solve` | non-conservative A，`J_phi/J_edge≈0.91` | **10/10 pass** |
| `uniform_E_real` | E/J/P closed | **10/10 pass** |

**P1 rotational-A key number** (full σ/f stable):

- `E_rel ≈ 0.304`，`J_phi/J_edge ≈ 0.910`，`P_ratio ≈ 0.828`
- simple box **No 4× J deviation** → TEAM7 Bind/B≈4 more likely to come from hole wall/source field/feedback

**Command** (recommended to refresh the screen less):

```bash
cd build
./tests/.../test_3d_ophelie_high_sigma_edge_flux_scaling/bin/... \
  --output-csv=../discussion_bundle/team7_l2_outputs/high_sigma_edge_flux_scaling_P1a_v4.csv \
  2>&1 | tee /tmp/p1a_v4.log | grep '^\[p1a\]\|rows=\|passed='
```

Detailed record: `HIGH_SIGMA_EDGE_FLUX_P1A_LOG.md` §v4

**Old v3 conclusion withdrawn**: constant A + φ-solve is not "path damage"; v3 uses E/J rel≈1 to judge failure and is not applicable to pure gauge case.

### 2026-06-11 — P2β solver-local complete implementation (**closed, equivalent to field-scale-restore**)

**Implementation** (physical `A_coil` / `A_ind` without scale/restore):
1. **φ Solve**: `phi_rhs` multiplies `scale` → PCG → `phi` on particles and `phi_rhs` divides `scale` to restore
2. **edge reconstruction**: `b_acc` is multiplied by `scale` in `ReconstructOphelieEdgeFluxElectricCurrentCK`, and `E` is divided by `scale` after solving.
3. Share `params.edge_flux_solver_local_rhs_scale_` (consistent with P0 measurement formula)

**TEAM7 one-way comparison** (filament, `source_scale=1`, 50 Hz):

| mode | Bind/B | J_imag_L2_vol | phase90 RMS | phi_eq_res_imag | smoke |
|------|--------|---------------|-------------|-----------------|-------|
| field-scale-restore | 3.987 | 63993 | 11.81 | ~4.3e-4 | 1 |
| **solver-local** | **3.987** | **63993** | **11.81** | **~4.3e-4** | **1** |

→ **Values ​​are consistent with post-restore field-scale**; provide a correct iteration basis for physical-scale Picard.

**CLI**：`--ophelie-edge-flux-normalization-mode=solver-local`

**Still open**: `Bind/B≈4` / `phase90≈12` are not improved (the two normalizations are equivalent, the bottleneck is not in the restore mechanism itself).

### 2026-06-11 — P2 Picard physical-scale + A_ind relaxation (**First version**)

**accomplish**:
- `runTeam7ComplexEdgeFluxPicardWithLog`: `combine → solve(A_total) → Biot → relax A_ind` per round (physical scale, no field restore)
- default `self_induction_relax_aind_=true`; legacy J slack: `--team7-picard-relax-aind=0`
- Picard CSV extensions: `normalization_mode/scale`, `aind_rel_*`, `J_imag_vol_norm`, etc.
- `solver-local` + L3 Picard is no longer marked `diagnostic_only` (L2 one-way is still diagnostic)

**Recommended Command** (filament, `source_scale=1`, 50 Hz):

```bash
./test_3d_ophelie_team7_complex_edge_flux \
  --reload=1 --team7-level=picard --native-dp-mm=3 --ophelie-edge-flux-complex=1 \
  --coil-source-model=filament-racetrack --team7-coil-source-scale=1.0 \
  --ophelie-edge-flux-normalization-mode=solver-local \
  --self-induction-relax=0.05 --self-induction-max-iter=8 \
  --team7-reference-dir=.../reference_data/team7
```

**First Run** (`relax=0.05`, 6 rounds):

| iter | Aind_rel | Bind/B | Bz_phase90_RMS |
|------|----------|--------|----------------|
| 0 | — | 3.987 | 11.81 |
| 2 | 0.655 | 4.195 | 11.28 |
| 5 | 0.288 | 4.340 | **7.11** |

`diagnostic_only=0`; not converged; phase90 improved, Bind/B increased slightly. Log: `discussion_bundle/team7_l2_outputs/logs/p2_picard_solver_local_r005.log`

### 2026-06-11 — P4 edge-flux operator scale audit (**first run**)

**accomplish**:
- New diagnosis: `electromagnetic_ophelie_edge_flux_operator_audit.h`
- CLI: `--team7-operator-audit=1` (triggered after L2 one-way / L3 Picard is solved)
- **P4.1**: Particle-by-particle `C_ij` statistics → `sum_j C_ij |r_ij|² / (σ_i Vol_i)`
- **P4.4**: Plate particle partitioning (`skin_h=dp`): interior/top/bottom/exterior_lateral/hole_lateral_proxy
- P4.2/P4.3 has been covered by P1a box benchmark (uniform E, harmonic A edge-only are both closed)

**Command** (filament, `source_scale=1`, solver-local, 50 Hz):

```bash
./test_3d_ophelie_team7_complex_edge_flux \
  --reload=1 --team7-level=one-way --native-dp-mm=3 --ophelie-edge-flux-complex=1 \
  --coil-source-model=filament-racetrack --team7-coil-source-scale=1.0 \
  --ophelie-edge-flux-normalization-mode=solver-local \
  --team7-operator-audit=1 --team7-output-tag=p4_audit \
  --team7-reference-dir=.../reference_data/team7
```

**P4.1 Global** (n_plate=49580):

| quantity | value |
|----|-----|
| `C_ij` [min, max, mean] | 685, 1.26e12, 2.33e11 |
| `sum(C·r²)` [min, max, mean] | 4.5e7, 1.57e8, 1.30e8 |
| `σ·Vol` [min, max, mean] | 0.956 (even) |
| **ratio** [min, max, mean, median] | 4.7e7, 1.64e8, **1.36e8**, 1.50e8 |

**P4.4 Partition** (Bind/B and `e_edge_em_mismatch`):

| Partition | n | Bind/B | \|J\|_vol | e_edge_em_mis | \|edge_res\|_vol |
|------|---|--------|-----------|---------------|------------------|
| interior | 31637 | **3.79** | 49814 | **0.068** | 854 |
| top_surface | 8143 | 3.87 | 28752 | 0.198 | 358 |
| bottom_surface | 7993 | **4.99** | 23516 | **0.239** | 389 |
| exterior_lateral | 1807 | 4.19 | 15300 | 0.202 | 235 |
| all | 49580 | 3.99 | 63993 | 0.144 | 1032 |

**First run conclusion**:
1. **Bind/B≈4 non-uniform**: The bottom/side boundary Bind/B of the board is higher (~4.2–5.0), and the interior ~3.8 → has both **boundary edge-flux effect** and global J amplitude.
2. **`e_edge_em_mismatch` is enlarged by 3×** at the surface layer (0.07 interior → 0.20–0.24 surface) → edge reconstruction is less closed than `−∇φ−ωA` at the boundary.
3. **`hole_lateral_proxy` n=0** (`skin_h=dp` heuristic misses the hole edge); requires hole bbox or level-set distance to subdivide.
4. **P4.1 ratio ~1.4e8 stable** (σ·Vol uniform) → Next step is to compare **continuous admittance scale** and `pairwiseNegativeLaplaceWeight` dimension (**It is recommended to discuss with GPT whether to change the C_ij scale**).

**Output**: `team7_edge_conductance_audit_*.csv`, `team7_edge_partition_audit_*.csv`; log `discussion_bundle/team7_l2_outputs/logs/p4_operator_audit_oneway.log`

### 2026-06-11 — P2 hole wall partition + Jn/Jt audit (**BFS hole grid, user rerun confirmation**)

**Changes relative to first run P4.4**:
- `hole_lateral_proxy` → **`hole_lateral`** (occupancy grid + BFS hole void detection)
- New partition records: `j_normal_l2`, `j_tangential_l2`, `Jn/Jt`, `e_normal_l2`
- Still **not patch C_ij**; P4.1 `moment_trace_mean` diagnostic only

**Command** (consistent with P1a v4 contemporaneous configuration):

```bash
cd build
./test_3d_ophelie_team7_complex_edge_flux \
  --reload=1 \
  --reload-dir=tests/extra_source_and_tests/3d_examples/particle_generation_TEAM7/bin/reload \
  --team7-level=one-way --native-dp-mm=3 --ophelie-edge-flux-complex=1 \
  --coil-source-model=filament-racetrack --team7-coil-source-scale=1.0 \
  --ophelie-edge-flux-normalization-mode=solver-local \
  --team7-operator-audit=1 \
  --team7-reference-dir=tests/extra_source_and_tests/3d_examples/reference_data/team7
```

**Global one-way** (tag=`f50_fil_asign_p1_p90_minus_src1.000_jraw`):

| quantity | value |
|----|-----|
| Bind/Bcoil | **3.986** |
| J_imag_L2_vol | 63941 |
| e_edge_em_mismatch | 0.143 |
| omegaA/gradPhi | 1.957 |
| P_recon | 57.7 W |

**hole grid**：`nx=98 ny=98 hole_void_cells=1295`

**P4.4 partition + Jn/Jt** (skin_h=3 mm):

| Partition | n | Bind/B | \|J\|_vol | e_edge_em_mis | \|edge_res\|_vol | Jn/Jt | \|Jn\| | \|Jt\| |
|------|---|--------|-----------|---------------|------------------|-------|--------|--------|
| interior | 32474 | 3.83 | 50155 | **0.050** | 981 | 0.84 | 32233 | 38426 |
| top_surface | 8138 | 3.87 | 28759 | 0.197 | 414 | **0.010** | 285 | 28758 |
| bottom_surface | 6896 | **5.04** | 21720 | **0.254** | 416 | **0.011** | 245 | 21718 |
| exterior_lateral | 1874 | 4.31 | 15110 | 0.206 | 267 | 0.038 | 579 | 15098 |
| **hole_lateral** | **199** | **2.76** | 6765 | **0.352** | 79 | **0.51** | 3069 | 6029 |
| all | 49581 | 3.986 | 63941 | 0.143 | 1176 | 0.59 | 32386 | 55132 |

**P4.1 conductance** (invariant level, diagnostics only):

- ratio median ≈ **1.50e8**，moment_trace_mean ≈ **1.30e8**
- **Prohibited** Directly change the C_ij experience coefficient accordingly

**P2 Conclusion**:
1. **`hole_lateral n=199`**: The hole wall partition gap has been closed (first run proxy n=0 has been repaired).
2. **Bind/B is still ~4×**, but **uneven partitioning**: bottom ~5.0, hole_lateral ~2.76, interior ~3.83 → Problem **non-uniform scaling**, it is not suitable to use the global J× constant to fix it.
3. **Hole wall Jn/Jt≈0.51**: The normal/tangential currents at the hole edge are of the same magnitude, `e_edge_em_mismatch` is the highest at the hole wall (0.35) → the hole boundary edge reconstruction is audited first.
4. **Top/bottom surface Jn/Jt≈0.01**: The upper and lower surface currents are nearly tangential (in line with thin plate induction), but the bottom Bind/B is still high.
5. Compare with **P1 rotational-A** (simple box J preserved ~91%, no 4×) → TEAM7 4× is more likely to come from **geometric boundaries/source fields/feedback** rather than systematic amplitude errors on the box by the generic operator.
6. probe CSV is still skipped (reference path has no Bz CSV) → L1/L2 probe comparison is not closed; operator audit CSV is valid.

**Output**:
- `discussion_bundle/team7_l2_outputs/team7_edge_partition_audit_one-way_f50_fil_asign_p1_p90_minus_src1.000_jraw.csv`
- `build/output/team7_edge_conductance_audit_one-way_f50_fil_asign_p1_p90_minus_src1.000_jraw.csv`
- Full log: `/tmp/team7_audit.log`

---

### 2026-06-12 — P4 A_ind / Lenz-law sign audit (**First version completed**)

**Code**: `diagnostics/electromagnetic_ophelie_aind_lenz_audit.h`

- `computeOphelieAIndLenzAuditFromFields`: Plate volume weighted B/A related
- `runTeam7AIndLenzAuditAfterOneWay`: TEAM7 one-way post-auto Biot + audit
- CLI: `--team7-aind-lenz-audit=1`; **Automatically enabled** when `--team7-operator-audit=1`

**Criterion Description** (to avoid misjudgment):

| Quantity | Meaning | TEAM7 Expectation |
|----|------|------------|
| `corr(B_ind, B_coil)` | In-phase complex-vector correlation | **≈0** (quadrature, non-sign error) |
| `corr_phasor(B_ind_i, B_coil_r)` | imag chain vs coil real cross correlation | **<0** ⇒ Lenz direction correct |
| `corr(B_ind_z, B_drive_z)` | box rotational-A only | **<0** ⇒ against uniform B₀ |
| `lenz_opposes_coil_b` | phasor or signed-dot criterion | 1 = pass |
| `lenz_audit_passed` | General | 1 = **Sign chain not reversed** (≠ Bind/B magnitude closed) |

**Command** (TEAM7, same configuration as P2):

```bash
cd build
./test_3d_ophelie_team7_complex_edge_flux \
  --reload=1 --reload-dir=tests/.../particle_generation_TEAM7/bin/reload \
  --team7-level=one-way --native-dp-mm=3 --ophelie-edge-flux-complex=1 \
  --coil-source-model=filament-racetrack --team7-coil-source-scale=1.0 \
  --ophelie-edge-flux-normalization-mode=solver-local \
  --team7-operator-audit=1 \
  --team7-reference-dir=tests/.../reference_data/team7
```

**TEAM7 one-way** (tag=`f50_fil_asign_p1_p90_minus_src1.000_jraw`, rerun on 2026-06-12):

| Quantity | Value | Interpretation |
|----|-----|------|
| Aind/Acoil | 3.92 | The amplitude of the induction chain is still ~4× |
| Bind/Bcoil | **3.986** | Consistent with P2 |
| Btotal/Bcoil | 4.11 | The total field is still too large |
| corr(Bind, Bcoil) | **0** | Orthogonal component, **cannot** be used to determine Lenz failure |
| corr(Btotal, Bcoil) | +0.24 | The total field is still partially in phase coil |
| corr_phasor(Bind_i, Bcoil_r) | **−0.588** | imag chain **against** coil real |
| mean(Bind_z_i) | −7.75e−3 | The induced z component is dominated by imag |
| lenz_opposes_coil | **1** | |
| lenz_pass | **1** | |

**P1 rotational-A box** (high_sigma, `rotational_A_uniform_B_phi_solve`, full σ/f):

| quantity | value |
|----|-----|
| corr(Bind_z, B_drive_z) | **−0.802** (stable) |
| lenz_vs_drive | **1** |
| lenz_pass | **1** |
| J_phi/J_edge | ~0.91 (no 4×) |

**P4 Conclusion**:

1. **Lenz symbol chain is correct on both box and TEAM7** (against driving field/coil quadrature), **exclude** "A_ind feedback enhancement rather than shielding" as the main cause of Bind/B≈4.
2. TEAM7 `corr(B_ind,B_coil)=0` **Not a bug**: The induced field and the coil field are nearly 90° out of phase, and the in-phase correlation fails; **phasor cross-correlation** must be used.
3. Bind/B≈4 is still a **amplitude/boundary/source field scaling** problem (consistent with the P2 partition conclusion), **not** Picard or Lenz sign alone can be repaired.
4. **Picard relax sweep is still frozen** (P5); in the next step, the hole wall edge reconstruction and L1 source reference are given priority instead of relaxation.

**Output**:

- `discussion_bundle/team7_l2_outputs/team7_aind_lenz_audit_one-way_f50_fil_asign_p1_p90_minus_src1.000_jraw.csv`
- `discussion_bundle/team7_l2_outputs/high_sigma_edge_flux_scaling_P1a_v4_lenz.csv` (including `corr_b_ind_z_b_drive_z`, `lenz_vs_drive` columns)

---

### 2026-06-12 — P4.1 conductance moment audit (**Completed**)

**Code**: `electromagnetic_ophelie_edge_flux_operator_audit.h`

- Particle-by-particle **M_i = Σ_j C_ij r_ij ⊗ r_ij** (SYCL kernel accumulation Matd)
- Host characteristic value: `getPrincipalValuesFromMatrix` → λ₁≥λ₂≥λ₃
- Diagnosis volume: `trace(M_i)`, `anisotropy=(λ₁−λ₃)/trace`, `trace/(σ·Vol)`, `ratio_median/inv(h³)`
- **trace consistency**: `|trace(M)−Σ C|r|²|/Σ C|r|²` max ≈ **6.6×10⁻⁷** (implemented correctly)
- Partition CSV adds moment column; **does not patch C_ij**

**Command** (same as P2/P4):

```bash
cd build
./test_3d_ophelie_team7_complex_edge_flux ... --team7-operator-audit=1 \
  --coil-source-model=filament-racetrack --team7-coil-source-scale=1.0 \
  --ophelie-edge-flux-normalization-mode=solver-local
```

**Global** (n=49581, dp≈3 mm, h_typ≈0.003 m):

| Quantity | Value | Interpretation |
|----|-----|------|
| trace(M) mean / median | 1.30e8 / 1.43e8 | ≈ Σ C|r|² (same magnitude as old ratio) |
| λ₁, λ₂, λ₃ mean | 4.61e7, 4.43e7, 3.98e7 | Nearly isotropic (λ₁/λ₃≈1.16) |
| anisotropy mean / median | **0.057 / 0.032** | Overall low anisotropy |
| ratio_med / inv(h³) | **4.05** | vs 1/h³≈3.7×10⁷ → ~4× dimensionless ratio |
| trace/(σ·Vol) median | 1.50e8 | Consistent with old conductance_ratio |

**Partition moment + border comparison** (skin_h=3 mm):

| Partition | trace/(σ·Vol) | anisotropy | λ₁ / λ₃ | Remarks |
|------|---------------|------------|---------|------|
| interior | 1.53e8 | **0.026** | ~1.08 | nearly isotropic |
| top_surface | 1.04e8 | **0.108** | ~1.40 | Boundary moment lower, more anisotropic |
| bottom_surface | 1.02e8 | **0.129** | ~1.50 | Same as above |
| exterior_lateral | 1.21e8 | 0.093 | ~1.30 | |
| **hole_lateral** | 1.03e8 | **0.129** | ~1.50 | Similar to bottom; e_edge_em is still the highest |
| all | 1.36e8 | 0.057 | ~1.16 | |

**P4.1 Conclusion**:

1. **M_i eigenvalue audit closure**: trace is consistent with Σ C|r|² value; old `moment_trace_mean` placeholder has been replaced by real tensor trace.
2. **Overall nearly isotropic** (aniso~0.03–0.06), **boundary/pore wall layer more anisotropic** (0.10–0.13) and trace/(σ·Vol) **lower than interior ~30%** → boundary pair geometry/neighbor distribution is different, **not** uniform C_ij scaling error.
3. `ratio_med/inv(h³)≈4` shows that the conductance moment is of the same order as kernel 1/h³, 1.5e8 ** alone cannot prove that C_ij requires empirical scaling** (consistent with the Round 3 decision).
4. **Prohibited** Based on this, patch C_ij; Bind/B≈4 still needs to be reconstructed from **boundary edge + source field/probe** direction.

**Output**:

- `discussion_bundle/team7_l2_outputs/team7_edge_conductance_audit_one-way_f50_fil_*.csv` (contains eigenvalue/anisotropy line)
- `discussion_bundle/team7_l2_outputs/team7_edge_partition_audit_one-way_f50_fil_*.csv` (with moment column)

---

### 2026-06-12 — P4.5 hole wall edge-reconstruction audit (**complete**)

**Code**: `electromagnetic_ophelie_edge_flux_operator_audit.h` — partition CSV extension + `printTeam7HoleLateralEdgeReconSummary`; operator audit path added `evaluateOphelieEdgeFluxImagEdgeDropMetrics`.

**New partition field**:

| Field | Meaning |
|------|------|
| `omega_a_over_grad_phi` | \|ωA\|/\|∇φ\| Volume norm ratio |
| `edge_recon_fallback_fraction` | edge LS fallback volume ratio |
| `edge_drop_l2` | paired edge_drop RMS |
| `j_sigma_e_mismatch` | \|J−σE\|/\|J\| |
| `j_edge_recon_mismatch` | \|J_edge−J\|/\|J\| |

**hole_lateral vs interior**（2026-06-12，filament, solver-local）：

| Quantity | hole_lateral | interior | Interpretation |
|----|--------------|----------|------|
| e_edge_em_mismatch | **0.352** | **0.050** | Hole wall **7×** — Main local notch |
| Bind/Bcoil | **2.76** | **3.83** | Hole wall amplitude **lower**, non-4× source |
| omegaA/gradPhi | 1.39 | 1.80 | Pore wall EMF balance slightly different |
| fallback_frac | **0** | **0** | LS pathological **not** main cause |
| j_sigmaE_mis | **0** | **0** | σE=J exactly on the storage field |
| j_edge_mis | **0** | **0** | J has been synchronized to edge solution |
| edge_drop_l2 | 2.69×10⁻⁴ | 1.51×10⁻⁴ | Slightly higher, still small in magnitude |
| Jn/Jt | 0.51 | 0.84 | The hole wall normal current is significant |

**P4.5 Conclusion**:

1. The hole wall problem is concentrated in the fact that **E_edge is inconsistent with E_em=−∇φ−ωA** (`e_edge_em`), **not** fallback, **not** σE≠J.
2. The global Bind/B≈4 **cannot** be attributed to the amplification of the pore wall alone (the pore wall Bind/B is lower); bottom ~5.0, interior ~3.8 → **Partitioning mechanism split**.
3. Round 3 **Diagnostic chain closed**; next step requires GPT decision: generic boundary BC vs RHS/restore vs L1 reference (see `discussion_bundle/OPHELIE_GPT_DISCUSSION_TEAM7_ROUND3_CLOSURE.md`).

**GPT Discussion Package**: `discussion_bundle/OPHELIE_GPT_DISCUSSION_TEAM7_ROUND3_CLOSURE.md` (§4 Six Final Questions)

---

### 2026-06-14 — P5 no-flux boundary (SPHinXsys normal + diagnostic chain)

> Plan: `OPHELIE_TEAM7_NO_FLUX_BOUNDARY_NORMAL_CURSOR_PLAN.md`
> Principle: **It is forbidden** to use analytical normal as the main solution; `NormalFromBodyShapeCK` → `NormalDirection` / `SignedDistance` must be used.

#### P5.0 — SPHinXsys normal + relax/reload integration ✅ **Success**

**Code**: `electromagnetic_ophelie_team7_boundary_normal.h`; `particle_generation_team7.cpp` / `test_3d_ophelie_team7_complex_edge_flux.cpp` access.

**Order**:

```bash
cd build
# Only relax + write reload (do not run EM)
./test_3d_ophelie_team7_complex_edge_flux --relax=1 --reload=0 --reload-dir=./reload_team7 --native-dp-mm=3
# normal audit
./test_3d_ophelie_team7_complex_edge_flux --reload=1 --reload-dir=./reload_team7 --native-dp-mm=3 \
  --team7-level=coil-only --team7-boundary-normal-audit=1
```

**acceptance**:

| items | results |
|----|------|
| `./reload_team7/Reload.xml` | Contains `NormalDirection` + `SignedDistance` |
| Particle number | coil **57947** / plate **49583** |
| `\|n\|` | ≈ **1**（audit CSV） |
| CLI | `--team7-boundary-normal-audit=1` has been added to `isOphelieTestCommandLineOption` filtering |

**Output**: `output/team7_boundary_normal_audit_coil-only_default.csv`

---

#### P5.1 — project-normal post-projection diagnosis ❌ **no-flux OK, EMF closure failed**

**CLI**: `--ophelie-edge-recon-boundary-mode=project-normal` (**diagnostic-only**, does not cover production E/J)

Boundary layer: `E ← E - n(n·E)`, `J ← σE`.

| partition | raw `e_edge_em` | projected | raw `\|Jn\|` → projected |
|------|-----------------|-----------|--------------------------|
| all | 0.137 | **0.142 ↑** | 17729 → 17678 |
| interior | 0.047 | **0.056 ↑** | — |
| hole_lateral | 0.338 | **0.380 ↑** | 484 → **≈0** |
| top/bottom | 0.19–0.23 | slightly worse | → **≈0** |

**Conclusion**: Normal projection **can force n·J≈0**, but **cannot** reduce `e_edge_em_mismatch`, but worsens slightly → **cannot be upgraded to production**.

**Output**: `output/team7_edge_recon_boundary_one-way_*.csv`

---

#### P5.2 — tangent-ls tangent plane LS diagnostic ❌ **Same as above mode**

**CLI**：`--ophelie-edge-recon-boundary-mode=tangent-ls`

| partition | raw | tangent-ls | `\|Jn\|` raw → tangent-ls |
|------|-----|------------|---------------------------|
| all | 0.137 | **0.234 ↑** | 17729 → 17678 |
| interior | 0.047 | **0.174 ↑** | — |
| hole_lateral | 0.338 | **0.443 ↑** | 484 → **≈0** |
| P_recon | 40.3 W | **30.8 W** | fallback_frac=**0** |

**Conclusion**: The tangent plane LS is also **no-flux effective, but the EMF is worse**; fallback=0 excludes LS pathology as the main cause.

---

#### P5.3 — Boundary operator consistency audit ✅ **Run through** (31 pass / 35 fail)

**CLI**: `--team7-boundary-consistency-audit=1` (can be run in conjunction with P5.1/P5.2)

**2026-06-14 running number** (L2 one-way, dp=3 mm, scale=0.754, `tangent-ls` joint running):

| partition | e_edge_em raw | EMF ωA/∇φ | σE=J | no-flux raw `\|Jn\|/\|Jt\|` | project/tangent no-flux |
|------|---------------|-----------|------|---------------------------|-------------------------|
| interior | **0.047 PASS** | PASS | PASS | **0.47 FAIL** | interior global Jn is still large; surface partition → ≈0 |
| top_surface | 0.186 FAIL | PASS | PASS | **0.016 PASS** | project/tangent **≈0 PASS** |
| bottom_surface | 0.235 FAIL | PASS | PASS | **0.021 PASS** | project/tangent **≈0 PASS** |
| hole_lateral | 0.338 FAIL | PASS | PASS | **0.087 PASS** | project/tangent **≈0 PASS** |
| all | 0.137 FAIL | PASS | PASS | **0.352 FAIL** | interior shell dominates the global Jn |

**Normal Mass** (SPHinXsys n vs Partition Resolution n):

| Partition | mean angle | max angle |
|------|------------|-----------|
| top/bottom/exterior | **4–6° PASS** | 55–86° FAIL (corner/edge) |
| hole_lateral | 22.7° FAIL | 38.7° PASS |
| interior (shell mislabeling) | 85.6° FAIL | 90° FAIL |

**Interpretation**: P5.3 to P4.5 Conclusion **Structured** - σE=J is always closed; **Shell partition + interior mislabeling** makes raw no-flux fail in all/interior; **Exterior surface/hole wall** raw no-flux has <0.1. Afterwards project/tangent passed no-flux on **labeled boundary shell**, but EMF mismatch still FAIL.

**Output**: `output/team7_boundary_consistency_one-way_*.csv`

---

#### P5.4 — Operator-level no-flux ghost edge ⏳ **Partially successful (no-flux↑, EMF is not closed)**

**CLI**: `--ophelie-edge-recon-boundary-mode=no-flux-ghost-edge` (**production**, modified `ReconstructOphelieEdgeFluxElectricCurrentCK`)

**Mechanism**: Boundary shell particles (`|SignedDistance| ≤ width`) append **e_ig=0** ghost edge constraints in LS reconstruction (`w_ig·n n^T` added to M, b unchanged), weight = neighbor maximum conductance; normals from SPHinXsys `NormalDirection`.

**2026-06-14 running numbers** (L2 one-way, dp=3 mm, scale=0.754, same as baseline reload):

| Indicators | baseline (none) | no-flux-ghost-edge | Judgment |
|------|------------------|---------------------|------|
| `e_edge_em` all | 0.137 | **0.138** | ❌ slightly worse |
| `e_edge_em` hole_lateral | 0.338 | **0.351** | ❌ Worse |
| `e_edge_em` interior | 0.047 | 0.048 | ≈ unchanged |
| `\|Jn\|/\|Jt\|` top | 0.016 | **0.010** | ✅ Improvement |
| `\|Jn\|/\|Jt\|` bottom | 0.021 | **0.014** | ✅ Improvement |
| `\|Jn\|/\|Jt\|` hole_lateral | 0.087 | **0.057** | ✅ Improvement |
| `\|Jn\|/\|Jt\|` all | 0.352 | 0.351 | ≈ unchanged (interior shell mislabeling dominates) |
| `P_recon` | 40.28 W | 40.27 W | ≈ unchanged |
| `Bind/Bcoil` | 3.91 | 3.91 | ≈ unchanged |
| P5.3 pass/fail (no post-hoc item) | — | **22 / 20** | — |

**in conclusion**:

1. **Operator-level ghost edge effectively reduces the outer surface/hole wall `\|Jn\|/\|Jt\|`** (production field, non-post-hoc).
2. **`e_edge_em_mismatch` The global situation is not improved**, but the hole wall is slightly worse → ghost constraint alone **is not enough to close the EMF**.
3. Consistent with P5.1/P5.2: **no-flux and EMF balancing are decoupling issues**; the next step requires **φ-RHS / missing moment** (level-set static confinement style) instead of just E rebuilding the LS patch.
4. **Disable** using ghost edge as final production default; left in opt-in diagnostic/experimental mode.

**Order**:

```bash
cd build
./test_3d_ophelie_team7_complex_edge_flux --reload=1 --reload-dir=./reload_team7 --native-dp-mm=3 \
  --team7-level=one-way --ophelie-edge-recon-boundary-mode=no-flux-ghost-edge \
  --team7-boundary-consistency-audit=1
```

**Output**: `output/team7_boundary_consistency_one-way_*.csv` (direct audit of production farm)

---

#### P5.5 — φ-RHS ghost + missing moment prototype ❌ **Invalid alone; combination≈P5.4**

**New CLI mode** (production, requires `NormalDirection`/`SignedDistance`):

| Mode | Function |
|------|------|
| `no-flux-missing-moment` | E reconstruct LS plus `(I - n n^T)` tangential moment compensation |
| `no-flux-phi-rhs-ghost` | φ-RHS assembly plus outward ghost pair: `rhs += c_ig ω a·(h n)` |
| `no-flux-full` | ghost edge + missing moment + phi-RHS ghost (all three are on) |

**Code**: `electromagnetic_ophelie_edge_flux_boundary_closure.h`; change `ComputeOphelieEdgeFluxPhiRhsFromASrcCK` + `ReconstructOphelieEdgeFluxElectricCurrentCK`.

**2026-06-14 Frequency sweep** (L2 one-way, dp=3 mm, scale=0.754, same as reload):

| mode | `e_edge_em` all | hole_lateral `e_edge_em` | hole `\|Jn\|/\|Jt\|` | top `\|Jn\|/\|Jt\|` | P5.3 pass/fail |
|------|-----------------|--------------------------|---------------------|---------------------|----------------|
| none（baseline） | 0.137 | 0.338 | 0.087 | 0.016 | 22 / 20 |
| `missing-moment` | **0.137** | **0.338** | 0.087 | 0.016 | 22 / 20 |
| `phi-rhs-ghost` | **0.137** | **0.338** | 0.087 | 0.016 | 22 / 20 |
| `no-flux-full` | **0.138** | **0.351** | **0.057** ✅ | **0.010** ✅ | 22 / 20 |

**in conclusion**:

1. **missing moment alone**: **No measurable changes** to `e_edge_em`, Jn/Jt, Bind/B** → ❌ Failure (weight/geometry model is insufficient or not the main cause).
2. **φ-RHS ghost alone**: Same as above, φ residuals and EMF indicators **unchanged** → ❌ failed (RHS first-order ghost is not enough to change the φ solution).
3. **`no-flux-full`**: The effect **completely comes from P5.4 ghost edge** (surface Jn/Jt↓, EMF is not closed); the other two items are superimposed **no additional benefits**.
4. **TEAM7 hole wall `e_edge_em≈0.35` is still a blocking item**; the next step requires **LHS/φ solver level** boundary closure (level-set static confinement v2) or GPT coupling instead of continuing to patch E/RHS.

**Order**:

```bash
cd build
#sub-test
./test_3d_ophelie_team7_complex_edge_flux ... --ophelie-edge-recon-boundary-mode=no-flux-missing-moment --team7-boundary-consistency-audit=1
./test_3d_ophelie_team7_complex_edge_flux ... --ophelie-edge-recon-boundary-mode=no-flux-phi-rhs-ghost --team7-boundary-consistency-audit=1
./test_3d_ophelie_team7_complex_edge_flux ... --ophelie-edge-recon-boundary-mode=no-flux-full --team7-boundary-consistency-audit=1
```

```

---

#### P5-fix — partition fix + tangent-LS distance_t A/B + EMF component ⏳ **Diagnostic infrastructure ready; tangent normalization part improves EMF**

**Motivation** (aligned `OPHELIE_TEAM7_P5_NO_FLUX_NEXT_CURSOR_PLAN.md`): legacy `interior` partitioning mislabels shell particles as in vivo, resulting in false positives for `|Jn|/|Jt|≈0.47`; tangent-LS uses 3D `distance` instead of tangential `distance_t` which may worsen shell EMF; fair compare (tangential EMF vs. tangential EMF) is required.

**Code changes**:

| subkeys | files | content |
|------|------|------|
| P5-fix-2 | `electromagnetic_ophelie_edge_flux_operator_audit.h` | `true_interior` / `boundary_shell_all` / `corner_or_edge` partition + `classifyTeam7ParticlePartition` |
| P5-fix-2 | `electromagnetic_ophelie_team7_boundary_consistency.h` | audit access the new partition; `true_interior`/`all` skip no-flux Jn/Jt evaluation |
| P5-fix-3 | `electromagnetic_ophelie_team7_edge_recon_boundary.h` | `em_tangential_mismatch`, `e_edge_em_mismatch_tangent_pair`, `projected_pair`, etc. |
| P5-fix-1 | Same as above + `electromagnetic_ophelie_parameters.h` | `--ophelie-tangent-ls-distance-norm=3d\|tangent` |

**2026-06-02 running numbers** (L2 one-way, dp=3 mm, scale=0.754, reload includes `NormalDirection`+`SignedDistance`, `--ophelie-edge-recon-boundary-mode=tangent-ls --team7-boundary-consistency-audit=1`):

**Partition Fix (P5-fix-2)** — legacy `interior` vs new partition:

| Partition | `e_edge_em` raw | `\|Jn\|/\|Jt\|` raw | Description |
|------|-----------------|---------------------|------|
| `true_interior` | **0.035** ✅ | 0.82 (in vivo, **not involved** no-flux rating) | Away from borders, good EMF closure |
| `boundary_shell_all` | 0.158 | **0.030** ✅ | Shell no-flux physics established |
| `interior`(legacy) | 0.047 | **0.467** ❌ | Still mixed with 16248 shell particles → False positive |

**Conclusion**: P5.3's previous `Jn/Jt≈0.47` of `interior` was a **partition mislabeling**, not a no-flux failure; the boundary should be evaluated as `boundary_shell_all`.

**tangent-LS distance normalized A/B (P5-fix-1)**:

| `distance_norm` | `all` tangent-ls `e_edge_em` | `boundary_shell` tangent-ls | legacy `interior` tangent-ls | hole_lateral tangent-ls `em_tan` |
|-----------------|------------------------------|----------------------------|------------------------------|----------------------------------|
| `3d` (default) | 0.234 ❌ | 0.290 ❌ | 0.174 ❌ | 0.140 |
| `tangent` | **0.145** ❌ | **0.168** ❌ | **0.062** ✅ | 0.197 |

- `tangent` normalization: **Significant improvement** in EMF of shell/body tangent-LS (`boundary_shell` 0.290→0.168, `interior` 0.174→0.062).
- **no-flux still holds**: `boundary_shell` / `Jn/Jt` for each surface tangent-ls ≈ 10⁻⁹ under both normalizations.
- **hole wall still poor**: hole_lateral `e_edge_em` tangent-ls 0.443→0.393; fair `em_tan` 0.140→0.197 (slightly worse).
- **Prohibit promotion to production**; `distance_norm=tangent` is only a default candidate for diagnosis.

**EMF component (P5-fix-3)** — hole_lateral example (`distance_norm=3d`):

| indicators | raw | tangent-ls |
|------|-----|------------|
| `e_edge_em` (unfair) | 0.338 | 0.443 |
| `em_tangential_mismatch`（fair） | 0.180 | **0.140** |
| `\|Jn\|/\|Jt\|` | 0.087 ✅ | ~10⁻⁸ ✅ |

→ tangent-LS **worse the overall e_edge_em but improve the tangential EMF fair compare**; the hole wall contradiction remains in the tangential component (0.14 still >0.1).

**P5.3 audit count changes**: New partition + more post-hoc items → pass/fail increased from 31/35 to **49/62** (3d) / **52/59** (tangent); `true_interior` In vivo EMF full PASS.

**Order**:

```bash
cd build
#Default 3d distance
./tests/.../test_3d_ophelie_team7_complex_edge_flux \
  --reload=1 --reload-dir=./reload_team7 --native-dp-mm=3 \
  --team7-level=one-way --ophelie-edge-recon-boundary-mode=tangent-ls \
  --team7-boundary-consistency-audit=1

# A/B: Tangential distance_t normalization
... --ophelie-tangent-ls-distance-norm=tangent
```

**decision making**:

1. ✅ Zoning infrastructure **available**; subsequent boundary evaluation is based on `boundary_shell_all` / each surface, no legacy `interior` is used to evaluate no-flux.
2. ⏳ tangent-LS + `distance_norm=tangent` **Improves shell EMF but does not close**; hole wall is still a blocking term.
3. ❌ Still **forbidden** P5.1/P5.2/tangent-LS is upgraded to production; production defaults to `boundary-mode=none`.
4. ~~**Next step** (GPT plan): P6a L1 source/probe audit; P6b box/slab φ-LHS/Neumann verification~~ → **P6a/P6b first run completed** (see below); TEAM7 large φ-LHS is still behind.

---

#### P6a — L1 source / reference / probe audit ✅ **Infrastructure + one-way first run**

**Implementation** (`team7/electromagnetic_ophelie_team7_l1_source_audit.h`):

- Pearson correlation, `best_fit_scale_l2/peak`, geometry/probe coordinate audit
- volume-racetrack vs filament-racetrack cross-check
- L2: `ref_phase0` vs `coil+ind_skin` comparison
- CSV：`team7_l1_source_audit_*.csv` + `team7_l1_source_audit_summary_*.csv`
- `coil-only` / `one-way` automatically enabled; CLI: `--team7-l1-source-audit=1`

**2026-06-02 running numbers** (L2 one-way, dp=3 mm, `source_scale=0.754`, f=50 Hz, volume-racetrack, a1-b1):

| Metrics | Values ​​| Access | Results |
|------|-----|------|------|
| `Bz_source RMS_rel` | **0.318** | <0.70 | ✅ |
| Pearson `corr` | **0.951** | >0.90 | ✅ |
| `best_fit_scale_peak` | **1.000** | ≈1 | ✅ |
| `best_fit_scale_l2` | 0.830 | ≈1 | Log |
| `volume_vs_filament_rms_rel` | **0.140** | <0.30 | ✅ |
| `ref_vs_coil+ind_skin_rms_rel` | **0.627** | <0.318 | ❌ (expected: phase0 ref non-pure source field) |

**Summary**: `pass=20` `fail=2` (diagnostic-only; phase0 reference contains sensing component, non-source-only fails).

**Order**:

```bash
cd build
./tests/.../test_3d_ophelie_team7_complex_edge_flux \
  --reload=1 --reload-dir=./reload_team7 --native-dp-mm=3 \
  --team7-level=one-way --team7-coil-source-scale=0.754 \
  --team7-reference-dir=.../reference_data/team7
```

**decision making**:

1. ✅ L1 source field and volume/filament cross-check **Available**; `scale=0.754` peak alignment confirmation (`best_fit_peak≈1`).
2. ❌ `ref_phase0≈coil+skin` Assumption **not true** (rms_rel=0.63); the phase0 reference still contains extra components and is not suitable for pure source field access control.
3. Coil peak x≈216 mm vs ref x≈126 mm still visible in P6a audit (geometry/probe segment).

---

#### P6b — box/slab border mode sweep ✅ **8/8 pass; SignedDistance fixed**

**Motivation**: Verify P5 boundary closure (none / ghost-edge / full / φ-Neumann) in isolation on a high-σ analytic box, without mixing in TEAM7 hole wall geometry.

**accomplish**:

| Documentation | Content |
|------|------|
| `diagnostics/electromagnetic_ophelie_p6b_box_boundary.h` | `initBoxEdgeFluxBoundaryNormalsFromShape`（`NormalFromBodyShapeCK` host + device sync） |
| `test_3d_ophelie_high_sigma_edge_flux_scaling.cpp` | `--p6b-boundary-sweep=1`：σ=1e7, f=50, dp=0.04; 4 boundaries × 2 case |

**Root cause fix**: Manual `registerStateVariableData(SignedDistance)` device delegate on SYCL is null → use the same TEAM7 model `SPHSolver` + `par_host` + `NormalFromBodyShapeCK` instead.

**2026-06-02 running number** (σ=1e7, f=50 Hz, dp=0.04, n=512):

| border mode | case | `gauge_pass` | `J_phi/J_edge` | `E_rel` | `boundary_Jn_rel` | `row_pass` |
|----------|------|--------------|----------------|---------|-------------------|------------|
| `none` | constant_A gauge | 1 | 3.57e-05 | 1.00 | 0.702 | ✅ |
| `none` | rotational_A | 0 | 0.910 | 0.304 | 0.274 | ✅ |
| `no-flux-ghost-edge` | constant_A gauge | 1 | 2.97e-05 | 1.00 | 0.589 | ✅ |
| `no-flux-ghost-edge` | rotational_A | 0 | 0.904 | **0.335** | **0.251** | ✅ |
| `no-flux-full` | constant_A gauge | 1 | 2.97e-05 | 1.00 | 0.589 | ✅ |
| `no-flux-full` | rotational_A | 0 | 0.903 | **0.335** | **0.251** | ✅ |
| `phi-neumann` | constant_A gauge | 1 | 3.57e-05 | 1.00 | 0.702 | ✅ |
| `phi-neumann` | rotational_A | 0 | 0.910 | 0.304 | 0.274 | ✅ |

**in conclusion**:

1. ✅ **gauge cancellation** (constant A + φ-solve) **remains closed** (`gauge_pass=1`, `J_phi/J_edge<10⁻⁴`) under all 4 boundaries.
2. ⏳ **rotational A** (non-conservative driver): ghost-edge/full is slightly worse than `E_rel` (0.304→0.335) and `boundary_Jn_rel` (0.274→0.251); `phi-neumann` is **same** as `none` (this case does not benefit from Neumann).
3. ✅ P5 ghost-edge/full closure **works on box**; does not damage gauge case; has slight side effects on rotational case.
4. The default value of production is still `boundary-mode=none`; the box result **cannot** extrapolate the TEAM7 hole wall.

**Order**:

```bash
cd build
./tests/.../test_3d_ophelie_high_sigma_edge_flux_scaling/bin/... \
  --p6b-boundary-sweep=1 --dp=0.04
# CSV: ./output/p6b_box_boundary_sweep.csv  (8 rows, pass=8)
```

---

#### P5.6-lite Finale - Are this round of modifications effective? **（2026-06-02）

**Short answer**: It is effective for **diagnosis and misjudgment removal**; it is basically invalid for **TEAM7 L2 core indicator closure**. In line with the expectations of the GPT program in the P5.6-lite phase, it is not sufficient to launch TEAM7 full φ-LHS static confinement.

| Work package | What works | What doesn't work | Judgment |
|--------|-----------|-----------|------|
| **P5-fix-2 partition** | Correct legacy `interior` treats shell as body → `Jn/Jt` False positive 0.47→True value ~0.03 | No improvement `e_edge_em`, Bind/B, phase90 | ✅ Required for diagnosis |
| **P5-fix-1 tangent distance_t** | tangent-ls shell EMF improvement (0.29→0.17) | hole wall `e_edge_em` still ~0.35–0.40 | ⏳ local |
| **P5-fix-3 EMF components** | Distinguish raw vs fair tangential mismatch | Unclosed hole wall | ✅ Diagnosis |
| **P6a L1 Audit** | Confirm scale≈0.754 peak alignment; volume/filament cross-check passed | phase0 ref impure source field (fail=2); peak position x 216 vs 126 mm still | ✅ Audit |
| **P6b box sweep** | gauge case 4 boundary full pass; prove closure can run | rotational slightly worse; φ-Neumann is the same as none | ✅ isolation verification |
| **P5.4 ghost-edge / full (TEAM7)** | Hole wall `Jn/Jt` 0.087→0.057 | ** Hole wall e_edge_em 0.338→0.351 (worse)**; Bind/B, phase90, J are almost unchanged | ❌ Do not advance production |

**TEAM7 one-way A/B**（dp=3 mm，scale=0.754，f=50 Hz，`--team7-boundary-consistency-audit=1`）：

| edge mode | `Bind/B` | `e_edge_em` all | `hole_lateral` e_edge_em | `hole` Jn/Jt | phase90 RMS | `J_imag_L2` |
|----------|----------|-----------------|--------------------------|--------------|-------------|-------------|
| `none` | 3.913 | 0.137 | **0.338** | 0.087 | 9.80 | 53405 |
| `no-flux-ghost-edge` | 3.911 | 0.138 | **0.351** ↑ | 0.057 ↓ | 9.80 | 53395 |
| `no-flux-full` | 3.911 | 0.138 | **0.351** ↑ | 0.057 ↓ | 9.80 | 53395 |

→ ghost-edge only suppresses the surface normal current component, **EMF consistency (hole wall) is slightly deteriorated**; L2 hard blocking term (Bind/B≈4, phase90≈10) **zero improvement**.

**Decision (aligned `OPHELIE_TEAM7_P5_NO_FLUX_NEXT_CURSOR_PLAN.md` §3.2)**:

1. ✅ **P5.6-lite List 1–5 has been completed**; the simple benchmark **does not provide evidence of the benefits of going back to TEAM7 to do large-scale P5.6**.
2. ❌ **Do not advance the full φ-LHS static confinement** on the TEAM7 hole (large investment, variable aliasing, tight 6-month window).
3. ✅ **Reserved** P5 no-flux modes are opt-in diagnostic; production defaults to `boundary-mode=none`.
4. **Choose one of the two next steps** (parallel lightweight is recommended):
- **A**: P6a deepening — `source_scale=1` + filament strict L1; prepare f=0 / source-only reference (phase0 ref cannot be used as pure source field confirmed);
- **B**: It is planned to evaluate **French reduced glass delivery** as the main line before July (TEAM7 is reduced to the benchmark branch line).

---

#### P5.4 — Decision (blocking production boundary patch)

1. **Disable** to upgrade P5.1/P5.2 post-projection to production; P5.4 ghost edge **opt-in test**, non-default.
2. **Still prohibited**: C_ij empirical scaling, J×constant, cancel φ-solve, edge-only production, a_sign=−1, Picard relax sweep.
3. ~~**Next step**: φ-RHS boundary closure + level-set missing moment~~ → P5.6-lite has been certified **No TEAM7 core income**; the hole wall `e_edge_em≈0.35` is still a contradiction, but it is not appropriate to pile up E/J patches.
4. SPHinXsys normal direction **ready** (P5.0); top/bottom/exterior normal direction quality is acceptable, and the interior shell label needs to be refined.

---

## 3. Current acceptance number (L2 one-way, dp=3 mm, scale=0.754)

See [README.md](README.md) for the command.

| Indicators | Typical values ​​| Access control |
|------|--------|------|
| `bz_rms_coil` | ~0.32 | ✓ < 0.70 |
| `bz_rms_coil_x_span` | ~0.29 | Log |
| `bz_rms_total` (phase0) | ~0.32 | Log |
| `bz_rms_total_phase90` | ~11.7 | Log |
| `bz_rms_total_mag` | ~1.46 | Log |
| `bz_rms_ind_imag_skin` | ~5.75 | Log(board top z≥z_top−2dp skin layer Biot) |
| `Bind/Bcoil` | ~3.91 | Log |
| `edge_flux_input_scale` | ~0.0053 | Log |
| `J_imag_L2_vol` | ~5.3e4 A/m²·√m³ | Log |
| `P_recon_W` | ~40 W | log (edge-recon power, not 50 kW calibrated) |
| `max_J_imag` | ~6.2e6 | EM Finiteness |
| `passed` | 1 | |

**Skin layer experiment (2026-06-02)**: Full board Biot phase90 RMS≈11.7 → only skin layer on top of board **≈5.75** (about half, still much higher than the reference). It shows that the phase90 deviation of the probe partly comes from **Biot superposition of the air probe by deep J in the board**, but it is not the only root cause.

---

## 4. Suggest discussion with GPT (2026-06-02)

0. **P3c–P5 update**: J driver; **feedback and restore indicators are the same**; Picard r=0.05 only phase90 7.1. GPT should make a decision: **Whether to change the edge-flux conductance/RHS scale**, or **TEAM7 dedicated coupling** instead of continuing to adjust the Biot/feedback/Picard super parameters.
1. **J amplitude vs `Bind/Bcoil≈4`**: Uniform `J×(1/3.9)` makes phase90 RMS from ~12→~1.7, but `σ·E=J` is self-consistent with `P_vol≈80 W` - is the higher one **Biot coupling to probe/reference** or restore's non-linear amplification of `J`? (P3c: `j_restore_linearity≈1`, exclude restore non-linearity)
2. **Reference phasor convention**: `auto` rear peak sim≈−1.17 mT vs ref≈+1.42 mT (phase90); is the reference CSV `B_ind_imag` or the total field orthogonal component?
3. **`P_recon` vs `P_vol`**: 40 W vs 80 W. Is factor 2 a complex power/active power definition difference? Should the TEAM7 50 kW target be connected to L2?
4. **racetrack + `source_scale=0.754`**: The coil phase0 is aligned, but the induction is still ~4× - is it necessary to filament the source or change `A_eff` to close the L2 coil and induction at the same time?

---

## 5. Known issues/To-do

| ID | Issue | Priority |
|----|------|--------|
| ~~P1~~ | ~~filament source access L1~~ | **Complete** (see table above) |
| ~~P2~~ | ~~L2 one-way filament bottoming out~~ | **Completed** |
| ~~P2~~ | ~~L3 Picard CSV per round~~ | **Complete** (not converged, see table above) |
| ~~P3~~ | ~~200 Hz Bz + Jey sim diagnosis~~ | **Complete** (Jey ref to be warehoused) |
| ~~P3c~~ | ~~edge-flux ω Scaling Audit~~ | **Completed** (`Bind/B∝f`, see table above) |
| ~~P3b~~ | ~~Jey reference storage + J/B separation diagnosis~~ | **Completed** (50 Hz; 200 Hz to be COMSOL table2) |
| P3b+ | COMSOL `table2.txt` → 200 Hz Jy column (`tools/import_team7_jey_comsol_table2.py`) | Low |
| ~~P2b~~ | ~~Picard relaxation sweep~~ | **Basing completed** (r=0.05 slightly improved phase90; Bind/B still ~4) |
| ~~P1a~~ | ~~high-σ uniform E / harmonic A benchmark~~ | **Complete** (edge-only closed; phi_solve failed, see above) |
| ~~P2β~~ | ~~solver-local φ+edge internal scaling~~ | **Done** (equivalent to field-scale post-restore, see above) |
| ~~P2 Picard~~ | ~~physical-scale + A_ind relax + extended CSV~~ | **First version completed** (see above; not converged) |
| ~~P4 Lenz~~ | ~~A_ind / Lenz-law sign audit~~ | **Complete** (2026-06-12; lenz_pass=1, Bind/B still ~4) |
| P4.1 moment | ~~conductance moment eigenvalue/anisotropy~~ | **Complete** (2026-06-12; aniso~0.06, boundary~0.13) |
| ~~P4.5~~ | ~~hole wall edge-recon audit~~ | **Complete** (e_edge_em hole wall 0.35; Bind/B hole wall 2.76) |
| P4+ | **GPT Round 3 Finale** → Boundary BC/RHS/L1 ref | **Blocking Implementation** |
| P5+ | edge-flux RHS/conductivity scaling or restore strategy (relative to TEAM7) | **High** |
| ~~P5.0~~ | ~~SPHinXsys normal direction + reload~~ | **Complete** (NormalDirection write reload) |
| ~~P5.1~~ | ~~project-normal diagnostics~~ | **Complete** (no-flux OK, EMF worse) |
| ~~P5.2~~ | ~~tangent-ls diagnosis~~ | **Complete** (same as above; fallback=0) |
| ~~P5.3~~ | ~~Boundary consistency audit~~ | **Complete** (31 pass / 35 fail) |
| P5.4 ghost | **Operator level no-flux ghost edge** | **Partial success** (Jn/Jt↓ surface; e_edge_em is not closed) |
| ~~P5.5a~~ | ~~missing moment recon~~ | **Failed** (Same indicator as baseline) |
| ~~P5.5b~~ | ~~φ-RHS ghost~~ | **Failed** (Same indicator as baseline) |
| ~~P5.5 full~~ | ~~combination of the three~~ | **≈P5.4** (no additional EMF income) |
| P5-fix | Partition correction + tangent distance_t A/B + EMF component | **Partially successful** (Partition clarification; tangent norm improves shell EMF; hole wall still poor) |
| ~~P6a~~ | L1 source/reference/probe audit | **Complete** (pass=20/fail=2; rms_rel≈0.32; corr≈0.95) |
| ~~P6b~~ | box/slab boundary sweep (none/ghost/full/φ-Neumann) | **Complete** (8/8 pass; gauge fully closed; rotational boundary slightly deteriorated E_rel) |
| P5.6 | φ-LHS level-set static confinement / GPT coupling | **Block L2 closure** |
| P1 | phase90: `x_high_bias` RMS≈12.5; `J_imag` full thickness uniformly larger | High |
| P1 | `Bind/Bcoil≈3.9` Each depth layer is ~3.5–4.9 | High |
| P2 | Coil peak position x≈216 mm vs reference x≈126 mm; `x_peak_ref` segment coil RMS is already better | Medium |
| P2 | Reduce full-plate `J_imag` amplitude (not skin Biot filter alone) | Medium |
| P3 | phase90 / \|B\| Incorporated into L2 hard access control | Low (reopen after repair) |
| P3 | L3 Picard Total Game Benchmark | Low |

---

## 6. Key CLI

| Options | Default | Description |
|------|------|------|
| `--team7-level=` | `one-way` | L1/L2/L3 |
| `--coil-source-model=` | `volume-racetrack` | `filament-racetrack` with TEAM7 centerline Biot |
| `--team7-coil-source-scale=` | **0.754** | volume-racetrack peak calibration to filament reference; `≠1` → diagnostic |
| `--team7-phase90-convention=` | `minus-imag` | probe layer phase90 mapping |
| `--team7-frequency=` | **50** | **50** or **200** (refer to CSV with multiple columns in the same file) |
| `--team7-omega-sweep=1` | off | Sweep frequency 25/50/100/200 Hz, write `team7_omega_scaling_*.csv` |
| `--team7-normalization-sweep=1` | off | Scan `safe_rhs_l2`, write `team7_normalization_sweep_*.csv` |
| `--ophelie-edge-flux-normalization-mode=` | `field-scale-restore` | `solver-local` (Picard recommended) |
| `--team7-picard-relax-aind=` | **1** | `0` = legacy J slack |
| `--self-induction-relax=` | 0.15 | Picard A_ind/J relaxation factor |
| `--team7-probe-line=` | `a1-b1` | `a2-b2` is the second muFEM probe line |
| `--team7-output-tag=` | automatic | like `f50_fil_asign_p1_p90_minus_src1.000_jraw` |
| `--team7-reference-dir=` | Automatic search | Reference CSV directory |
| `--ophelie-aind-one-way-feedback=` | **0** | L2 is turned off by default A_total quadratic solution |
| `--team7-ind-j-post-scale=` | off | `auto` or value: J amplitude diagnostic calibration (do not rerun edge-flux) |
| `--native-dp-mm=3` | 3 | Same as reload |
| `--ophelie-edge-recon-boundary-mode=` | `none` | Diagnosis: `project-normal`/`tangent-ls`; production: `no-flux-ghost-edge`/`no-flux-missing-moment`/`no-flux-phi-rhs-ghost`/`no-flux-full` |
| `--team7-boundary-normal-audit=1` | off | P5.0: SPHinXsys normal `\|n\|`, SignedDistance audit |
| `--team7-boundary-consistency-audit=1` | off | P5.3: Partitioned EMF/no-flux/normally consistent CSV |
| `--ophelie-edge-flux-complex=1` | on | complex edge-flux |
| `--ophelie-edge-flux-imag-a-sign=` | **1** | `ωA` item symbol in imag chain edge-drop (`-1` flip sense B phase) |
| `--team7-probe-phase90-neg-ind=1` | off | Probe total field uses `−B_ind_imag` (equivalent combination with `a_sign=-1`) |

---

## 7. Output file

| path | content |
|------|------|
| `output/team7_probe_<level>.csv` | Probe decomposition (including skin / particle columns) |
| `output/team7_bz_<level>_coil_only.csv` | Coil phase0 comparison |
| `output/team7_bz_<level>_total.csv` | Total field phase0 comparison |
| `output/team7_validation_history.txt` | **Additional** running record (including `phase90_x_bands` segment) |
| `output/team7_probe_phase90_xbands_<level>.csv` | x segment phase90 RMS |
| `output/team7_plate_depth_profile_<level>.csv` | Board z Stratified J/B Statistics |
| `output/team7_picard_<tag>.csv` | Picard per round `Bind/Bcoil`, phase90 RMS, J_rel |
| `output/team7_omega_scaling_<tag>.csv` | ω Sweep: `Bind/B`, `J_imag`, `P_recon`, `ωA/∇φ`, etc. |
| `output/team7_jey_probe_<tag>.csv` | Jy probe line sim vs ref |
| `output/team7_j_vs_b_probe_split_<tag>.csv` | J/B separation: judge whether J is too large vs Biot is too large |
| `output/high_sigma_edge_flux_scaling_P1a.csv` | P1a high-σ box benchmark (30 lines) |
| `output/team7_boundary_normal_audit_*.csv` | P5.0 SPHinXsys normal audit |
| `output/team7_edge_recon_boundary_*.csv` | P5.1/P5.2 raw vs project vs tangent-ls comparison |
| `output/team7_boundary_consistency_*.csv` | P5.3 Partition operator consistency (pass/fail matrix) |
| `output/PlateBody_ite_0000000000.vtp` | Plate EM field + depth bin (L2/L3) |
