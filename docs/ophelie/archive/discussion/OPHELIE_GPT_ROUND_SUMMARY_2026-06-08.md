# OPHELIE Edge-Flux Round work summary & GPT discussion bundle

> **round**：2026-06-02 ~ 2026-06-08
> **Status**：Stage 1 edge-flux production **alreadypassed**；**need GPT discussion befor e fixing Stage 2 direction**
> **Purpose**：Upload this file + the attachment manifest belowto ChatGPT；for each later round copy this template to create `OPHELIE_GPT_ROUND_SUMMARY_YYYY-MM-DD.md`

---

## 1. Round goal and conclusion (one sentence)

**Goal**: from div-grad baseline（`phi_eq_res≈0.49`）pivot to SPH pairwise **edge-flux production**；fixPowerCalibrationError 。
**conclusion**: edge φ / edge residual **succeeded**; failure was treating graph edge energy as physical Joule power; calibration changed to **`P_recon`**; H production **`production_literature_passed=1`**.

---

## 2. Timeline and completed work

### Phase A — Edge-flux Stage 1 Implementation（2026-06-02）

| item | content |
|----|------|
| coreoperator | `ComputeOphelieEdgeFluxPhiRhsFromASrcCK`、`ComputeOphelieEdgeFluxResidualCK`、`ComputeOphelieEdgeFluxJouleHeatCK`、`ReconstructOphelieEdgeFluxElectricCurrentCK` |
| Main switch | `--ophelie-current-form=edge-flux` |
| edge-drop | `e_ij = (φ_j-φ_i) + ω·Ā_ij·(x_j-x_i)`（not legacy `g_ij·(A_i-A_j)`） |
| document | `OPHELIE_EDGE_FLUX_SOLVER_DEVELOPMENT_PLAN.md`；Archive `docs/ophelie/archive_phi_divgrad_residual_floor//README.md` |
| Test | `test_3d_ophelie_edge_flux_sign`（sign MMS，passed=1） |

**first H run data issue**（power-fix befor e）：

```text
P_graph_edge = 50000 W （Error Calibrationtarget）
P_total_recon = 0.43 W
P_graph/P_recon ≈ 116402
ampere_turns_eff ≈ 1.9 A （baseline need ~660 A）
production_literature_passed = 0
```

### Phase B — source code directory consolidation（2026-06-08）

[`electromagnetic_ophelie/README.md`](file:///home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/README.md) Notes on classification：

| directory | content |
|------|------|
| main directory | edge-flux + φ + French + shared infrastructure |
| `legacy/div_grad/` | div-grad fallback |
| `diagnostics/` | MMS / operator audit / A_ind diagnostic |
| `team7/` | Team7 Case |
| `french_extras/` | glass relax |
| `stage2/` | self-induction Picard |

### Phase C — Power fix（2026-06-08，per ChatGPT Plan Step 1–4）

**Code change points**：

1. **Calibration / Acceptance**：`joule_power_raw` switch to `P_total_recon`；`P_graph_edge` diagnostic only
2. **Fields**：`JouleHeatEdgeRecon = 0.5·σ·|EEdgeRecon|²`；VTP Output
3. **Acceptance**：`P_graph/P_recon` only warning， not fail production
4. **Test**：
 - extend `test_3d_ophelie_edge_flux_sign`（pair-level `edge_drop_*`）
 - new `test_3d_ophelie_edge_flux_power_uniform_field` (both cases `P_recon/P_exact≈1`)

**Modified source files**（see §5）：

- [`electromagnetic_ophelie_edge_flux.h`](file:///home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_edge_flux.h)
- [`electromagnetic_ophelie_french_literature.h`](file:///home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_french_literature.h)
- [`electromagnetic_ophelie_field_names.h`](file:///home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_field_names.h)
- [`electromagnetic_ophelie_register_fields.h`](file:///home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_register_fields.h)

### Phase D — Validationrun data Step 5–6（2026-06-08）

#### French reduced three-way v2（@50 kW）

| Case | CLI | phi_eq_res | edge_res_red | P_raw | ampere_turns | divJ_L2_red | production_lit | passed |
|------|-----|------------|--------------|-------|--------------|-------------|----------------|--------|
| **H** | `--ophelie-current-form=edge-flux` | 8.0e-5 | **10955** | 50000 (**P_recon**) | **660** | 0.116 | **1** | **1** |
| **A** | `--ophelie-current-form=particle-gradient --phi-projection-operator=div-grad` | 0.486 | 0.024* | 50000 (particle) | 660 | **2.06** | **1** | **1** |
| **G** | `--ophelie-current-form=particle-gradient --phi-projection-operator=edge-flux` | 9.6e-5 | **341** | 50000 (particle) | 660 | 0.93 | 0 | 0** |

\* A edge_res as legacy diagnosticin div-grad solution。
\** G failed as**expected**（phi-diagnostic only）。

#### H power-fix befor eaftercomparison

| metric | v1（Error Calibration） | v2（P_recon Calibration） |
|------|----------------|-------------------|
| `P_recon` | 0.43 W | **50000 W** |
| `ampere_turns_eff` | ~1.9 | **~660** |
| `edge_res_red` | 10920 | 10955 |
| `production_literature_passed` | 0 | **1** |
| `P_graph/P_recon` | 116402 | 116402（warning only） |

#### Unit / MMS Test

| Test | result |
|------|------|
| `test_3d_ophelie_edge_flux_sign` | passed=1 |
| `test_3d_ophelie_edge_flux_power_uniform_field` potential | P_recon/P_exact=1.0；P_graph/P_exact≈12232 |
| same induction | P_recon/P_exact=1.0；P_graph/P_exact≈12232 |

#### A_ind one-way baseline（Step 6）

| Quantity | value |
|----|-----|
| `A_ind/A_coil` | 0.436 |
| `B_ind/B_coil` | 0.366 |
| `phi_eq_res_vol` | 0.410（div-grad φ，not yetCalibration 50 kW） |
| `passed` | 1 |

stilluse **div-grad φ + particle J**，not yethook up `JEdgeRecon`。

---

## 3. currentarchitecture（for GPT context）

```text
--ophelie-current-form=edge-flux （production）
 → legacy-pairwise φ LHS
 → electromotive-drop RHS: e_ij = (φ_j-φ_i) + ω Ā_ij·(x_j-x_i)
 → edge residual CK（SYCL Inner neighbor pair）
 → JouleHeatEdgeGraph（diagnostic）+ JouleHeatEdgeRecon（physical）
 → P_recon = Σ 0.5·σ·|E_edge|²·V_i → 50 kW Calibration

--ophelie-current-form=particle-gradient --phi-projection-operator=div-grad （A fallback）
--ophelie-current-form=particle-gradient --phi-projection-operator=edge-flux （G phi-diagnostic only）
```

**forbidden**（Stage 1 already observed）：external mesh/CSR、host-only pair loop、use `P_graph_edge` as production Calibration。

---

## 4. Issues for GPT to answer（P0 → P2）

### P0 — Decide Stage 2 first (please answer first)

**Q1. E/J post-processingswitch to `JEdgeRecon`？**

- observation：H on `divJ_L2_red≈0.116`，A baseline `≈2.06`；E/J stilluse particle gradient。
- Option：
 - (a) Stage 2 first `EImag/JImag` switch to edge reconstruction output, and use edge divJ as acceptance；
 - (b) accept particle divJ only as baseline comparison; edge production no longer chase divJ；
 - (c) otheroption？

**Q2. A_ind = K[J] J sourceusewhich？**

- current A_ind diagnostic use particle `JImag` + div-grad φ。
- edge-flux production already has `JEdgeRecon`。
- should should ：**first Q1 then A_ind Picard**？or still A_ind can continue use particle J？

**Q3. `P_graph_edge` handle going forward？**

- `P_graph = Σ 0.25·C_ij·edge_drop²`；`C_ij` = pairwise Laplace weight。
- uniform MMS on `P_graph/P_exact≈10⁴`，`P_recon/P_exact≈1`。
- Option：permanently diagnostic / do κ Calibration / fromcodeDelete？

### P1 — Acceptance and regression

**Q4.** French onshouldneed to **`I×2 → P_recon×4`** as as literature regression？or still uniform-field MMS enough？

**Q5.** edge-flux production literature gate should**permanently remove** `divJ_L2_red`？still add back edge-version divJ gate in Stage 2？

**Q6.** `JouleHeatEdgeRecon` spatial distribution (stronger outer ring, weaker interior): need a magnitude gate, or only look at total `P_recon`？

### P2 — historicalroute（div-grad era，optional）

**Q7.** div-grad `phi_eq_res≈0.49` treat as formulation floor, preserve only as fallback？

**Q8.** compatible-div-grad / Neumann / grad correction etc. archived routeshouldstopinvest？

---

## 5. Recommended GPT output for mat

Ask GPT to reply in the following structure so Cursor can hook up execution directly:

```text
1. decisiontable：Q1–Q6 each item → Option + rationale
2. Stage 2 implementationorder（3–5 step，includes CLI/Fields/test names）
3. explicit「not doing」manifest（avoid scope creep）
4. If Q1=(a): give minimal diff scope for E/J switch (which CK / fields to change)
```

---

## 6. Attachment manifest — files to upload to GPT

Paths are relative to repository root **`/home/yyc/SPHinXsysSYCL/`**. Links are clickable in Cursor.

### 6.1 Must -read this round (highest priority)

| File | Note |
|------|------|
| [ this File](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/OPHELIE_GPT_ROUND_SUMMARY_2026-06-08.md) | work summary + issue manifest |
| [FRENCH_THREE_WAY_SUMMARY_v2.txt](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/FRENCH_THREE_WAY_SUMMARY_v2.txt) | H/A/G three-waySummary |
| [OPHELIE_EDGE_FLUX_POWER_FIX_AND_NEXT_DEVELOPMENT_PLAN.md](file:///home/yyc/SPHinXsysSYCL/OPHELIE_EDGE_FLUX_POWER_FIX_AND_NEXT_DEVELOPMENT_PLAN.md) | ChatGPT oniterationPlan（alreadyexecution Step 1–6） |
| [OPHELIE_PHI_DEVELOPMENT_LOG.md §4.14](file:///home/yyc/SPHinXsysSYCL/OPHELIE_PHI_DEVELOPMENT_LOG.md) | development log (search `4.14`） |

### 6.2 Run logs (v2 primary; v1 for comparison)

| File | Note |
|------|------|
| [french_H_edge_flux_production_v2.log](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/french_H_edge_flux_production_v2.log) | **H production passed** |
| [french_H_edge_flux_production.log](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/french_H_edge_flux_production.log) | H v1（Error P_edge Calibration） |
| [french_A_baseline_v2.log](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/french_A_baseline_v2.log) | A div-grad baseline |
| [french_G_edge_flux_phi_only_v2.log](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/french_G_edge_flux_phi_only_v2.log) | G phi-only（expected fail） |
| [french_aind_diagnostic_v1.log](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/french_aind_diagnostic_v1.log) | A_ind one-way baseline |
| [FRENCH_THREE_WAY_SUMMARY.txt](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/FRENCH_THREE_WAY_SUMMARY.txt) | v1 three-way（includes F_compatible） |

### 6.3 Run-data CSV (v2)

| File | Case |
|------|------|
| [ophelie_edge_flux_H_v2.csv](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/ophelie_edge_flux_H_v2.csv) | H |
| [ophelie_phi_p0_H_v2.csv](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/ophelie_phi_p0_H_v2.csv) | H |
| [ophelie_edge_flux_A_v2.csv](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/ophelie_edge_flux_A_v2.csv) | A |
| [ophelie_phi_p0_A_v2.csv](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/ophelie_phi_p0_A_v2.csv) | A |
| [ophelie_edge_flux_G_v2.csv](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/ophelie_edge_flux_G_v2.csv) | G |
| [ophelie_phi_p0_G_v2.csv](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/ophelie_phi_p0_G_v2.csv) | G |

### 6.4 Core source (power fix + edge-flux primary path)

| File | Note |
|------|------|
| [electromagnetic_ophelie_edge_flux.h](file:///home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_edge_flux.h) | edge CK + P_recon + JouleHeatEdgeRecon |
| [electromagnetic_ophelie_french_literature.h](file:///home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_french_literature.h) | Calibration + acceptance |
| [electromagnetic_ophelie_edge_flux_diagnostics.h](file:///home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_edge_flux_diagnostics.h) | edge diagnostic |
| [electromagnetic_ophelie_phi.h](file:///home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_phi.h) | φ solve + edge RHS |
| [electromagnetic_ophelie_cli.h](file:///home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_cli.h) | CLI |
| [electromagnetic_ophelie_parameters.h](file:///home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_parameters.h) | parameter / threshold |
| [electromagnetic_ophelie_field_names.h](file:///home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_field_names.h) | field names |
| [electromagnetic_ophelie_register_fields.h](file:///home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_register_fields.h) | Fieldsregister / VTP |
| [electromagnetic_ophelie_postprocess.h](file:///home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_postprocess.h) | E/J particle post-processing（**Stage 2 may change**） |
| [electromagnetic_ophelie/README.md](file:///home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/README.md) | directory classification |

### 6.5 Test cases

| path | Note |
|------|------|
| [test_3d_ophelie_french_reduced.cpp](file:///home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/test_3d_ophelie_french_reduced.cpp) | French literature Primary acceptance |
| [test_3d_ophelie_edge_flux_sign.cpp](file:///home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_edge_flux_sign/test_3d_ophelie_edge_flux_sign.cpp) | edge-drop sign MMS |
| [test_3d_ophelie_edge_flux_power_uniform_field.cpp](file:///home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_edge_flux_power_uniform_field/test_3d_ophelie_edge_flux_power_uniform_field.cpp) | uniform field Power MMS |
| [test_3d_ophelie_french_aind_diagnostic.cpp](file:///home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_aind_diagnostic/test_3d_ophelie_french_aind_diagnostic.cpp) | A_ind one-way |

### 6.6 Plan / archive documents (background)

| File | Note |
|------|------|
| [OPHELIE_EDGE_FLUX_SOLVER_DEVELOPMENT_PLAN.md](file:///home/yyc/SPHinXsysSYCL/OPHELIE_EDGE_FLUX_SOLVER_DEVELOPMENT_PLAN.md) | Stage 1–3 master plan |
| [docs/ophelie/archive_phi_divgrad_residual_floor//README.md](file:///home/yyc/SPHinXsysSYCL/docs/ophelie/archive_phi_divgrad_residual_floor//README.md) | div-grad ArchiveIndex |
| [OPHELIE_GPT_DISCUSSION_BLOCKERS.md](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/OPHELIE_GPT_DISCUSSION_BLOCKERS.md) | deferred items quick reference |

### 6.7 optional reference (historical v1 CSV / old logs)

<details>
<summary>Expand — not required for upload</summary>

| File |
|------|
| [ophelie_edge_flux_H.csv](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/ophelie_edge_flux_H.csv) |
| [ophelie_phi_p0_H.csv](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/ophelie_phi_p0_H.csv) |
| [french_A_baseline.log](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/french_A_baseline.log) |
| [french_G_edge_flux_phi_only.log](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/french_G_edge_flux_phi_only.log) |
| [french_compatible-div-grad.log](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/french_compatible-div-grad.log) |
| [ophelie_vector _divergence_mms_scan_v2.csv](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/ophelie_vector _divergence_mms_scan_v2.csv) |

</details>

### 6.8 VTP output (ParaView, optional)

H v2 run data Glass VTP (if in `build/output/`）：

```text
/home/yyc/SPHinXsysSYCL/build/output/GlassBody_ite_0000000000.vtp
```

includes `JouleHeatEdgeRecon`, `EEdgeRecon`, `JEdgeRecon`, etc.

---

## 7. Reproduce commands (GPT or manual validation)

W or king directory: `/home/yyc/SPHinXsysSYCL/build`

```bash
# Build
cmake --build. --target test_3d_ophelie_french_reduced \
 test_3d_ophelie_edge_flux_sign \
 test_3d_ophelie_edge_flux_power_uniform_field \
 test_3d_ophelie_french_aind_diagnostic -j$(nproc)

# MMS
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_edge_flux_sign/bin/test_3d_ophelie_edge_flux_sign
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_edge_flux_power_uniform_field/bin/test_3d_ophelie_edge_flux_power_uniform_field

# H production (needs build/reload/Reload.xml)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced \
 --reload=1 --literature-mode --state_recording=0 \
 --ophelie-current-form=edge-flux --phi-edge-flux-diagnostics=1

# A baseline
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced \
 --reload=1 --literature-mode --state_recording=0 \
 --ophelie-current-form=particle-gradient --phi-projection-operator=div-grad

# G phi-only
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced \
 --reload=1 --literature-mode --state_recording=0 \
 --ophelie-current-form=particle-gradient --phi-projection-operator=edge-flux

# A_ind
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_aind_diagnostic/bin/test_3d_ophelie_french_aind_diagnostic --reload=1
```

---

## 8. Paste template for GPT (opening message)

```text
I am working on OPHELIE-like French reduced edge-flux on SPHinXsysSYCL.
Stage 1 complete: edge-flux production literature_passed=1 (P_recon calibration @50 kW).
Please read OPHELIE_GPT_ROUND_SUMMARY_2026-06-08.md and FRENCH_THREE_WAY_SUMMARY_v2.txt.
Focus on answering §4 Q1–Q6 and give Stage 2 implementation order in §5 for mat.
Background: SPH/SYCL inner-neighbor pairs; no mesh/CSR; P_graph_edge demoted to diagnostic.
```

---

## 9. Next-round discussion document template (Cursor internal)

For each later round create:

```text
discussion_bundle/OPHELIE_GPT_ROUND_SUMMARY_YYYY-MM-DD.md
```

Structure: §1 conclusion → §2 completed work → §3 architecture → §4 issues → §5 expected reply for mat → §6 attachment manifest → §7 reproduce commands

Also update:[`OPHELIE_PHI_DEVELOPMENT_LOG.md`](file:///home/yyc/SPHinXsysSYCL/OPHELIE_PHI_DEVELOPMENT_LOG.md) new section; [`OPHELIE_GPT_DISCUSSION_BLOCKERS.md`](file:///home/yyc/SPHinXsysSYCL/discussion_bundle/OPHELIE_GPT_DISCUSSION_BLOCKERS.md) blocker items.

---

*Generated: 2026-06-08 | Cursor agent | dev log §4.14*
