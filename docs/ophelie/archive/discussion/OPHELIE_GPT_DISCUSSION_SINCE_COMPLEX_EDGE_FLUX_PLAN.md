# OPHELIE — full work record since `OPHELIE_COMPLEX_EDGE_FLUX_SOLVER_PLAN.md` and GPT discussion bundle

> **Purpose**：upload for ChatGPT discussion under one iteration direction。
> **baseline document**：repository root [`OPHELIE_COMPLEX_EDGE_FLUX_SOLVER_PLAN.md`](../OPHELIE_COMPLEX_EDGE_FLUX_SOLVER_PLAN.md)（GPT issued complex edge-flux primary path decision and Stage 3.0–3.6 plan）。
> **detailed development log**：[`OPHELIE_PHI_DEVELOPMENT_LOG.md`](../OPHELIE_PHI_DEVELOPMENT_LOG.md) §4.14–§4.25。
> **updated**：2026-06-02（includes Stage 4 thermal one-way / diffusion / VTP + EM dp sweep + Picard gate）。

---

## 0. three-sentence summary for GPT

1. **Complex edge-flux primary pathalready landed per plan**（Stage 3.0a–3.7）：dual imag/real chains、`Q=0.5σ(|E_r|²+|E_i|²)`、A_ind one-way feedback、complex Picard joint convergence（J_rel + phi_eq_res），French reload `passed=1`。
2. **thermal coupling one-way closed loop**（Stage 4.0–4.2）：frozen Q → explicit thermal step → optional diffusion + cold-wall BC；~12% energy gap as **T₀=300 floating-point rounding**，not EM physics error；already use `OphelieThermalDeltaT` fix。
3. **explicit not doing / deferred**：natural convection（needs fluid-body–EM–thermal full coupling）、σ(T) feedback、Jacoutot full-field transient quantity comparison；next discussion should focus on **Picard vs literature 50 kW calibration relationship、acceptance definition/scope、role of real chain under coil-only、dp/reload convergence policy**。

---

## 1. GPT plan decision review (baseline)

[`OPHELIE_COMPLEX_EDGE_FLUX_SOLVER_PLAN.md`](../OPHELIE_COMPLEX_EDGE_FLUX_SOLVER_PLAN.md) coredecision：

| decision | content |
|------|------|
| abandon real-only as final OPHELIE route | self-induction feedback `A_ind=K[J]` need **A_real A_imag** enter edge-drop simultaneously |
| A_ind phase audit | Stage 3.0: validate `K[J_imag]→A_ind_imag` dominance; **not** replace complex implementation |
| complex equations | `E_r=-∇φ_r+ωA_i`，`E_i=-∇φ_i-ωA_r`；`Q=0.5σ(|E_r|²+|E_i|²)` |
| Stage order | 3.0 audit → 3.1 Fields → 3.2 componentized kernel → 3.3 dual-scalar solve → 3.4 complex power → 3.5 A_ind one-way → 3.6 Picard |
| P0 fix | A_ind diagnostic CLI + edge-flux RHS initialization order |

**prerequisites completed before plan**（real-only edge-flux Stage 1–2，see development log §4.13–§4.16）：

- Production Calibrationswitch to **`P_recon`**（not graph edge energy `P_graph`）
- French H @50 kW：`production_literature_passed=1`，`edge_res_red≈10⁴`
- I² scaling：`P(0.5I)/P(I)=0.25`，`P(2I)/P(I)=4.0`（B3 closure）

---

## 2. stage overview (after plan → today)

| Stage | Status | one sentence |
|------|------|--------|
| **3.0a** A_ind phase audit + CLI/RHS fix | ✅ | `A_ind_imag/A_src≈0.43`，imag dominant |
| **3.1** Complex Fieldsregister | ✅ | `phi/e/j_edge_recon_{real,imag}`、`joule_heat_edge_recon_complex` |
| **3.2** Edge-flux componentized | ✅ | `OphelieEdgeFluxComponent`，`a_sign=±1` |
| **3.3** dual-scalar solve | ✅ | imag GMRES + real PCG/GMRS |
| **3.4** Complex Joule heat + I² | ✅ | MMS + scaling `passed=1` |
| **3.5** A_ind one-way + feedback | ✅ | `P: 7.28→8.64 W`（+18.8%），`max_J_real≈22` |
| **3.6** Picard self-induction | ✅ | joint gate，reload ~7 iter，`picard_converged=1` |
| **3.7** complex uniform-field power MMS | ✅ | Plan §8.2 five cases full `passed=1` |
| **4.0** thermal one-way | ✅ | French end-to-end，`E_joule=E_thermal` |
| **4.1** thermal diffusion + cold wall Dirichlet | ✅ | `E_thermal<E_joule`（heat loss expected） |
| **4.2** literature thermal material properties + temperature VTP | ✅ | Jacoutot Table 1 preset + stepped VTP |
| **4.25** EM dp sweep + Picard independent gate consolidation | ✅ | lattice dp sweep； and `literature_passed` decoupled |
| **4.2+** natural convection / σ(T) / Literaturefullfieldcomparison | ❌ not yetdo | needs fluid solver; large engineering |

---

## 3. detailed work record

### 3.1 Stage 3.0a — A_ind phaseauditP0 fix

**plan requirement**：fix `test_3d_ophelie_french_aind_diagnostic` CLI；edge-flux branch first `setup/finalize RHS` then `solvePhiImagWithCurrentRhs`。

**Completed**：

- hook into `filterFrenchReducedCommandLine` + `filterOphelieTestCommandLine`
- `runFrenchReducedAIndOneWayDiagnostic` complete RHS pipeline
- `test_3d_ophelie_edge_flux_sign` add imag/real dual-chain edge-drop zero test

**key observations（reload, complex=1）**：

| Quantity | value | meaning |
|----|-----|------|
| `A_ind_real` | **0** | coil-only source: real chain not yet active |
| `A_ind_imag/A_src_real` | **≈0.43** | consistent with plan expected O(0.3–0.4) |
| `B_ind/B_coil` | **≈0.36** | same order of magnitude |

**conclusion**：phase audit **supports** complex edge-flux necessity；real-only active-A cannot feedback `A_ind_imag`。

---

### 3.2 Stage 3.1–3.3 — Complex Fields and dual-chain solve

**added/extended fields**（`electromagnetic_ophelie_field_names.h` / `register_fields`）：

```text
PhiReal, PhiImag
EEdgeReconReal, EEdgeReconImag
JEdgeReconReal, JEdgeReconImag
EdgeFluxResidualReal, EdgeFluxResidualImag
JouleHeatEdgeReconComplex (+ real/imag component fields)
ATotalReal/Imag, AIndReal/Imag,...
```

**componentized**（`electromagnetic_ophelie_edge_flux.h`）：

```text
Imag: active_a=ATotalReal, a_sign=+1 → PhiImag → E_i, J_i
Real: active_a=ATotalImag, a_sign=-1 → PhiReal → E_r, J_r
edge_drop = (φ_j-φ_i) + a_sign·ω·Ā_ij·Δx_ij
```

**French pipeline** (`electromagnetic_ophelie_french_literature.h`): with `--ophelie-edge-flux-complex=1`, imag GMRES → real PCG → complex Joule heat.

**regression (reload @50 kW)**：

| Test | complex=0 | complex=1 |
|------|-----------|-----------|
| `test_3d_ophelie_french_reduced` | P_recon=50 kW, prod_lit=1 | same as left |
| `test_3d_ophelie_edge_flux_sign` | passed=1 | passed=1（dual-chain q_antisym≈1e-7） |

**expected behavior**：coil-only `A_imag≈0` → real RHS≈0 → `PhiReal≈0`（not a bug）。

---

### 3.3 Stage 3.4 — Complex PowerI² law

**Implementation**：`joule_heat_complex = 0.5·σ·(|E_r|²+|E_i|²)`；scalar `P_complex = Σ Q·Vol`。

**Test**：

| Test | result |
|------|------|
| `test_3d_ophelie_edge_flux_scaling` | imag + complex suite：`P(0.5I)/P(I)=0.25`，`P(2I)/P(I)=4.0`，`passed=1` |
| `test_3d_ophelie_edge_flux_power_uniform_field` | 5 cases（see §3.6），`P_recon/P_exact=1.0` |

**fix**：absolute calibration on MMS box with `target_P` not appropriate → scaling Testonly gate I² ratio。

---

### 3.4 Stage 3.5 — A_ind one-way + feedback resolution

**CLI**：`--ophelie-aind-one-way-feedback=0|1`（default 1）

**flow**：

```text
A_coil → complex edge-flux → JEdgeRecon{Real,Imag}
→ Biot K[J] → A_ind{Real,Imag}
→ A_total = A_coil + A_ind → re-solve complex edge-flux（one-time feedback）
```

**reload typical values（complex edge-flux）**：

| metric | coil-only | after feedback |
|------|-----------|----------------|
| `P_complex` (W) | 7.28 | **8.64** (+18.8%) |
| `max_J_real` | 0 | **21.9** |
| `max_J_imag` | ~64 | ~64 |
| `A_ind_imag/A_src` | — | 0.43 |
| `edge_res_red_real` | — | **-1 (n/a)**（denominator guard when level0≈0） |

**Test**：`test_3d_ophelie_french_aind_diagnostic --reload=1 --ophelie-edge-flux-complex=1` → `passed=1`

---

### 3.5 Stage 3.6 — Picard self-induction prototype

**Implementation**：`stage2/electromagnetic_ophelie_self_induction.h` + `runFrenchReducedSelfInductionPicard`（`aind_diagnostic.h`）

**outer loop**：

```text
A_total^k = A_coil + A_ind^k
→ complex edge-flux solve + J reconstruction
→ A_ind^{k+1} = K[J]（under-relaxation）
→ until J_rel and phi_eq_res both met
```

**joint gate (Stage 3.6 closure, §4.19)**：

- `ophelieSelfInductionPicardConverged(j_rel, phi_eq_res_vol, params)`
- complex：`phi_tol=0.01`（`--self-induction-phi-tol=`）；legacy：`phi_tol=0.65`
- **independent of `literature_passed` (experimental path)

**reload + complex edge-flux reference run data**：

```bash
test_3d_ophelie_french_self_induction_picard --reload=1 \
 --ophelie-edge-flux-complex=1 --self-induction-max-iter=8
```

| metric | value |
|------|-----|
| `complex_picard_path` | 1 |
| `self_induction_iters` | 7 |
| `final_J_rel` | 0.039 |
| `phi_eq_res_vol` | 3.4×10⁻⁴ |
| `picard_converged` | 1 |
| `P_joule_W` | 6.39（**not yet**do 50 kW literature Calibration） |
| `A_ind/A_coil` | 0.422 |
| `max_J_real` / `max_J_imag` | 18.2 / 59.3 |
| `passed` | 1 |

**Picard parameter sweep**（`test_3d_ophelie_french_self_induction_picard_sweep`）：

- relax ∈ {0.15, 0.3, 0.45} × max_iter ∈ {6, 8, 12}
- reference (0.3, 8) convergence + enough budget grid points fully converge → `sweep_passed=1`

**engineering mistake (already fixed)**：Picard test first iteration CLI not hooked → ran imag-only path but `passed=1`（§4.18）。

**CLI minor fix (2026-06-02)**: when `--ophelie-edge-flux-complex=1` and user did not explicitly set `--ophelie-current-form=`, auto-select `edge-flux`.

---

### 3.6 Stage 3.7 — complex uniform-field power（Plan §8.2）

**Test**：`test_3d_ophelie_edge_flux_power_uniform_field`

| Case | physical setup | P_recon/P_exact |
|------|----------|-----------------|
| potential_field | H v4 potential field | 1.0 |
| induction_field | H v4 pure induction | 1.0 |
| complex_potential_real | PhiReal=-E0·x | 1.0 |
| complex_induction_imag_chain | AReal=A0 | 1.0 |
| complex_induction_real_chain | AImag=A0 | 1.0 |

**summary passed=1**

---

### 3.7 Stage 4.0 — Joule → thermal one-way (no EM feedback)

**header**：`diagnostics/electromagnetic_ophelie_joule_to_heat_one_way.h`

**pipeline**：

```text
runFrenchReducedEmPipeline (complex edge-flux)
→ syncOphelieJouleHeatPrimaryFor ThermalOneWay
→ ApplyOphelieJouleHeatOneWayTemperatureStepCK（full device）
→ ReduceDynamicsCK closure diagnostic
```

**French reload acceptance (before fix → after fix)**：

| metric | fixbefor e | fixafter |
|------|--------|--------|
| `E_joule_J` / `E_thermal_J` | 21.84 / **19.12** | **21.8365 / 21.8365** |
| `vol_weighted_rel_err` | ~0.12 | **0** |
| `closure_mismatch_vol_frac` | ~0.66 | **0** |
| `passed` | 0 | **1** |

**Test**：

- `test_3d_ophelie_joule_to_heat_one_way`（uniform Q MMS）
- `test_3d_ophelie_complex_joule_to_heat_one_way`（complex Q MMS）
- `test_3d_ophelie_french_complex_joule_to_heat_one_way`（French end-to-end）

---

### 3.8 Stage 4.22 — thermal closure ~12% gap root cause（Important）

**observation**：~11799/20500 particles Q>0 but apparent ΔT≈0；ParaView Q distribution reasonable。

**root cause**: direct add on `T₀=300 K` via `T += ΔT` (ΔT~10⁻⁵); IEEE-754 rounding swallows increment. **Not** EM/Q error; also **not** SYCL double-buffer main cause.

**fix**：

1. new field **`OphelieThermalDeltaT`**：accumulate from 0
2. **`Temperature = T₀ + OphelieThermalDeltaT`**
3. closure diagnostic reads `OphelieThermalDeltaT`，not read `T−T₀`
4. **forbidden** EM after `syncVariableToDevice(Q)` (Q stays device-authoritative)

---

### 3.9 Stage 4.1 — thermal diffusion + cold crucible Dirichlet

**header**：`diagnostics/electromagnetic_ophelie_thermal_diffusion_one_way.h`

**equation (explicit one-way)：

```text
OphelieThermalDeltaT += Q·dt/(ρ·cp) − (k/(ρ·cp))·L_pair(T)·dt
T = T₀ + OphelieThermalDeltaT
```

**BC**：French cylindrical shell `OphelieThermalBoundaryMask` → Dirichlet `T=T₀`

**French reload（`--thermal-diffusion=1`）**：

| mode | E_joule | E_thermal | Note |
|------|---------|-----------|------|
| no diffusion | 21.84 J | 21.84 J | closure 0 |
| with diffusion | 21.84 J | **11.85 J** | cold-wall heat loss, **expected** E_thermal < E_joule |

**gate**：diffusion mode does not require pointwise Joule closure；require `E_thermal ≤ E_joule` + `boundary_compliance > 90%`

**MMS**：`test_3d_ophelie_thermal_diffusion_mms` → `passed=1`

---

### 3.10 Stage 4.2 — literature thermal material properties + temperature VTP

**new files**：

- `electromagnetic_ophelie_french_thermal_material.h` — Jacoutot 2008 Table 1 @ 1473 K preset
- `electromagnetic_ophelie_thermal_vtp.h` — device→host sync afterwrite VTP

**preset comparison**：

| | reduced (regression default) | literature |
|--|---------------------|------------|
| ρ | 2500 | 2750 |
| cp | 1200 | 1150 |
| k | 1 | 4 |
| T₀ | 300 K | 1473 K |

**VTP Fields**：`Temperature`, `OphelieThermalDeltaT`, `JouleHeat`, `JouleHeatEdgeReconComplex`；diffusion mode includes LaplaceT / Conductivity / BoundaryMask

**CLI**：`--thermal-state-recording=1`, `--thermal-record-interval=N`, `--use-literature-thermal=1`

---

### 3.11 Stage 4.25 — EM dp sweep + Picard gate consolidation

**new test**：`test_3d_ophelie_french_em_dp_scan`

| mode | content |
|------|------|
| smoke | single dp=0.08 lattice |
| `--em-dp-scan` | dp∈{0.08,0.06,0.04} EM；finest `phi_eq_res<0.01` |
| `--em-dp-scan-picard` | same + Picard；finest must `picard_converged` |

**CSV**：`output/ophelie_french_em_dp_scan.csv`

**observation**：lattice φ **not necessarily** monotone with dp（0.06→0.04 rebound）；**reload @ dp=0.02** as quality anchor。

---

## 4. current production acceptance definition/scope (H vs A vs G)

from [`FRENCH_THREE_WAY_SUMMARY_v3.txt`](FRENCH_THREE_WAY_SUMMARY_v3.txt)（real-only full edge-flux; pre-complex baseline; complex φ/P same order）：

| Case | main gate | phi_eq_res | edge_res_red | P@50kW | divJ_L2_red* | prod_lit |
|------|---------|------------|--------------|--------|--------------|----------|
| **H** full edge-flux | edge_res + P_recon + q_antisym | ~1e-4 | ~9340 | P_recon | 0.13† | **1** |
| **A** div-grad | divJ + P_particle | 0.427 | ~0.02 | P_particle | 2.34 | 1 |
| **G** phi-only | divJ（old） | ~1e-4 | ~477 | P_particle | 0.93 | **0** |

\* H particle divJ **diagnostic only, not production gate。

**Complex=1 under French @50 kW**：`demo_passed=1`，`P_scaled=50000 W`，`phi_eq_res_vol≈6e-4`。

---

## 5. known open issues (for GPT discussion)

### 5.1 still-effective blockers（from [`OPHELIE_GPT_DISCUSSION_BLOCKERS.md`](OPHELIE_GPT_DISCUSSION_BLOCKERS.md)）

| ID | observation | current handling |
|----|------|----------|
| **B1** | `P_graph/P_recon ≈ 1.15×10⁵` | diagnostic only, not a gate |
| **B2** | particle `divJ_L2_red ≈ 0.13` on H | E/J post-processing still uses particle grad；edge Main pathuse `JEdgeRecon` |
| **B4** | div-grad `phi_eq_res≈0.49` floor | Preserved A as fallback |

### 5.2 issues appearing after plan

| Issue | details |
|------|------|
| Picard `P_joule≈6.4 W` vs EM one-way `7.28 W` | Picard path not yet literature-calibrated；comparable to 50 kW case？ |
| Picard convergence `J_rel≈0.039` borderline | `j_tol=0.05` loose？tighten or switch relative/absolute mixed criterion？ |
| A_ind one-way `edge_res_red_real=-1` | real chain level0≈0；feedback meaningful only after real chain active |
| lattice dp sweep φ not monotone | 0.04 worse than 0.06；use reload for dp study？ |
| thermal diffusion `E_thermal≈54% E_joule`（3×1s） | cold wall + short time constant；**not** bug；literature transient comparison still early |
| natural convection | **explicit not doing (needs NS + Boussinesq + coupling) |

---

## 6. explicit not doing / next-stage boundary

| item | cause |
|----|------|
| natural convection | melt as fluid body，buoyancy + momentum equation；incompatible with current solid-body + one-way thermal |
| σ(T) ↔ EM Picard | temperature changes conductivity then Joule；needs coupled architecture |
| Jacoutot full-field transient quantity comparison | needs flow + correct BC + possibly finer grid |
| Picard ↔ thermal feedback | plan and stage 4 both forbid |

**one-way milestone already complete**：EM → frozen Q → thermal（+diffusion）→ VTP ParaView sanity check。

---

## 7. file manifest (recommended GPT upload)

### 7.1 Plan and Record（Priority）

| File | Note |
|------|------|
| [`OPHELIE_COMPLEX_EDGE_FLUX_SOLVER_PLAN.md`](../OPHELIE_COMPLEX_EDGE_FLUX_SOLVER_PLAN.md) | **baseline** GPT direction |
| [`OPHELIE_PHI_DEVELOPMENT_LOG.md`](../OPHELIE_PHI_DEVELOPMENT_LOG.md) | §4.14–§4.25 detailed record |
| [`OPHELIE_EDGE_FLUX_SOLVER_DEVELOPMENT_PLAN.md`](../OPHELIE_EDGE_FLUX_SOLVER_DEVELOPMENT_PLAN.md) | real-only Stage 1 Plan（prerequisite） |
| [`OPHELIE_EDGE_FLUX_STAGE2_CLEANUP_AND_CURSOR_PLAN.md`](../OPHELIE_EDGE_FLUX_STAGE2_CLEANUP_AND_CURSOR_PLAN.md) | Stage 2 definition / scope |
|this document | GPT discussion bundle |
| [`OPHELIE_COMPLEX_EDGE_FLUX_RUN_RECORD.md`](OPHELIE_COMPLEX_EDGE_FLUX_RUN_RECORD.md) | run data and fix record |
| [`OPHELIE_GPT_DISCUSSION_BLOCKERS.md`](OPHELIE_GPT_DISCUSSION_BLOCKERS.md) | deferred topics |
| [`FRENCH_THREE_WAY_SUMMARY_v3.txt`](FRENCH_THREE_WAY_SUMMARY_v3.txt) | H/A/G comparison |
| [`docs/ophelie/FRENCH_REDUCED_CASE_ASSUMPTIONS.md`](../docs/ophelie/FRENCH_REDUCED_CASE_ASSUMPTIONS.md) | geometry/material assumptions |

### 7.2 core source code

| File | Stage |
|------|-------|
| `tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_edge_flux.h` | 3.2–3.4 complex kernels |
| `.../electromagnetic_ophelie_phi_component.h` | 3.3 componentized phi solve |
| `.../electromagnetic_ophelie_field_names.h` | 3.1 Fields |
| `.../electromagnetic_ophelie_french_literature.h` | French EM pipeline |
| `.../electromagnetic_ophelie_cli.h` | CLI（includes complex auto edge-flux） |
| `.../diagnostics/electromagnetic_ophelie_aind_diagnostic.h` | 3.5 + Picard entry |
| `.../stage2/electromagnetic_ophelie_self_induction.h` | 3.6 Picard |
| `.../diagnostics/electromagnetic_ophelie_joule_to_heat_one_way.h` | 4.0 thermal one-way |
| `.../diagnostics/electromagnetic_ophelie_thermal_diffusion_one_way.h` | 4.1 diffusion |
| `.../diagnostics/electromagnetic_ophelie_french_thermal_material.h` | 4.2 literature material properties |
| `.../diagnostics/electromagnetic_ophelie_thermal_vtp.h` | 4.2 VTP |
| `.../electromagnetic_ophelie_biot_savart.h` | K[J] Biot–Savart |

### 7.3 tests / cases (per stage)

| test directory | role |
|----------|------|
| `test_3d_ophelie_edge_flux_sign` | 3.0 edge-drop sign |
| `test_3d_ophelie_edge_flux_power_uniform_field` | 3.7 Power MMS |
| `test_3d_ophelie_edge_flux_scaling` | 3.4 I² |
| `test_3d_ophelie_french_aind_diagnostic` | 3.5 A_ind one-way |
| `test_3d_ophelie_french_self_induction_picard` | 3.6 Picard gate |
| `test_3d_ophelie_french_self_induction_picard_sweep` | 3.6 parameter sweep |
| `test_3d_ophelie_french_reduced` | French EM main case + VTP |
| `test_3d_ophelie_french_em_dp_scan` | 4.25 dp sweep |
| `test_3d_ophelie_joule_to_heat_one_way` | 4.0 MMS |
| `test_3d_ophelie_complex_joule_to_heat_one_way` | 4.0 complex MMS |
| `test_3d_ophelie_french_complex_joule_to_heat_one_way` | 4.0–4.2 French thermal e2e |
| `test_3d_ophelie_thermal_diffusion_mms` | 4.1 diffusion MMS |

**reload particles**：`build/reload/Reload.xml`（n≈20500，dp=0.02）

### 7.4 run artifacts / CSV (optional upload)

| path | Note |
|------|------|
| `discussion_bundle/french_H_edge_flux_production_v4.log` | H @50 kW |
| `discussion_bundle/french_H_edge_flux_scaling_reload.log` | I² scaling |
| `discussion_bundle/french_aind_diagnostic_v1.log` | early A_ind (div-grad path) |
| `build/output/ophelie_french_em_dp_scan.csv` | dp sweep |
| `build/output/GlassBody_ite_*.vtp` | temperature/Q fields (stage 4.2) |
| `discussion_bundle/ophelie_phi_p0_H_v2.csv` | φ P0 diagnostic |

---

## 8. recommendedreproducecommand （build/ under）

```bash
# --- Complex EM @ reload ---
./tests/.../test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced \
 --reload=1 --ophelie-edge-flux-complex=1 --target-power=50000

# --- A_ind one-way + feedback ---
./tests/.../test_3d_ophelie_french_aind_diagnostic/bin/test_3d_ophelie_french_aind_diagnostic \
 --reload=1 --ophelie-edge-flux-complex=1

# --- Picard (independent gate)---
./tests/.../test_3d_ophelie_french_self_induction_picard/bin/test_3d_ophelie_french_self_induction_picard \
 --reload=1 --ophelie-edge-flux-complex=1 --self-induction-max-iter=8

# --- thermal one-way closure ---
./tests/.../test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
 --reload=1 --ophelie-edge-flux-complex=1

# --- thermal diffusion + cold wall ---
./tests/.../test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
 --reload=1 --ophelie-edge-flux-complex=1 --thermal-diffusion=1

# --- literature material properties + temperature VTP ---
./tests/.../test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
 --reload=1 --ophelie-edge-flux-complex=1 --thermal-diffusion=1 \
 --use-literature-thermal=1 --thermal-steps=30 --thermal-state-recording=1 --thermal-record-interval=5

# --- EM dp sweep ---
./tests/.../test_3d_ophelie_french_em_dp_scan/bin/test_3d_ophelie_french_em_dp_scan --em-dp-scan-picard

# --- Stage 3.7 Power MMS ---
./tests/.../test_3d_ophelie_edge_flux_power_uniform_field/bin/test_3d_ophelie_edge_flux_power_uniform_field
```

---

## 9. issue manifest for ChatGPT

### 9.1 complex edge-flux and Picard (EM core)

1. **Picard power definition/scope**：after Picard convergence on reload `P_joule≈6.4 W`（no 50 kW calibration），while one-way EM `≈7.28 W`、feedback one-way `≈8.64 W`。how to interpret the three values？should Picard re-accept under `literature-mode` + 50 kW？
2. **J_rel gate**：current `J_rel≈0.039` passes under `j_tol=0.05`；for induction should switch to **absolute residual**、**relative L2(J)** or **P_complex rate of change**as convergence criterion？
3. **real chain role**：coil-only `A_imag≈0` → real chain dormant；after feedback/Picard `max_J_real≈18` becomes significant。does this match phasor physics expectation？PCG vs GMRES preference for real chain？
4. **A_ind amplitude**：`A_ind/A_coil≈0.42` compare magnitude to literature or TEAM7？next: segmented crucible / finer dp A_ind validation, or defer until σ(T)？
5. **edge_res_red_real=-1**：ratio undefined when level0≈0。should acceptance change to「check edge_res on real chain only when `‖J_r‖>ε`」？

### 9.2 acceptance definition/scope and divJ (continuing B2/B4)

6. **production main metrics**：Stage 1–2 already decided `edge_res_red + P_recon + q_antisym`；maintain in complex + Picard era？include `phi_eq_res_vol<0.01` in production (currently Picard gate only)？
7. **particle divJ**：H on `divJ_L2_red≈0.13` permanently diagnostic-only？if divJ must return ~1, fix full edge divergence operator not grad post-processing？
8. **P_graph**：`P_graph/P_recon~10⁵` worth conductance normalization MMS，or or permanently deprecate graph Joule heat？

### 9.3 thermal one-way and literature (stage 4)

9. **thermal closure fix**：`OphelieThermalDeltaT` accumulation insufficient for long steps at T₀=1473 K literature preset，or wait for flow coupling; use double-precision accumulator / incremental energy bookkeeping？
10. **diffusion + cold wall**：3×1 s internal `E_thermal≈54% E_joule` reasonable？for Jacoutot sanity: extend `--thermal-steps` + VTP first, or require flow？
11. **σ(T) EM feedbackPriority**：before doing NS，worth fixed σ(T) lookup one-way (still no EM feedback)，or wait for fluid module？

### 9.4 next-step route choices (GPT please prioritize)

12. under no-natural-convection constraint, rank the four items below？
 - (a) Picard @ 50 kW literature Calibration + Jacoutot power comparison
 - (b) reload dp convergence（relax Particle @ 0.015/0.01）
 - (c) A_ind / Picard convergence curve CSV + standardized VTP field output
 - (d) thermal long-time steps + literature T₀=1473 K ParaView sanity report
 - (e) pause EM, wait for fluid module

13. **when to enable σ(T)–EM Picard**？minimal MMS needed before hooking feedback？

14. **div-grad fallback (case A)**：with complex edge-flux production, freeze case A or keep investing？

---

## 10. Recommended GPT discussion order

1. first confirm **§9.1 Picard Power and convergencedefinition / scope**（most affects EM next steps）
2. reconfirm **§9.2 production gate** need tuning
3. revisit **§9.3 thermal** only maintain one-way sanity
4. finally **§9.4 priority ranking**

---

## 11. related document links

- development log：[`OPHELIE_PHI_DEVELOPMENT_LOG.md`](../OPHELIE_PHI_DEVELOPMENT_LOG.md)
- Complex run data:[`OPHELIE_COMPLEX_EDGE_FLUX_RUN_RECORD.md`](OPHELIE_COMPLEX_EDGE_FLUX_RUN_RECORD.md)
- Blockers：[`OPHELIE_GPT_DISCUSSION_BLOCKERS.md`](OPHELIE_GPT_DISCUSSION_BLOCKERS.md)
- French thermal e2e README：`tests/.../test_3d_ophelie_french_complex_joule_to_heat_one_way/README.md`
- Picard README：`tests/.../test_3d_ophelie_french_self_induction_picard/README.md`
