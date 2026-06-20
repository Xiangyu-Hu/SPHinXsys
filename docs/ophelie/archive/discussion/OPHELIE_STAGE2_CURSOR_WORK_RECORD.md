# OPHELIE Stage 2 — Cursor execution full record (since `OPHELIE_EDGE_FLUX_STAGE2_CLEANUP_AND_CURSOR_PLAN.md`)

> **time span**：2026-06-08（ChatGPT issued Stage 2 plan → Cursor implementation → user French reload run data → A_ind supplementary run data）
> **basis plan**：[`/home/yyc/SPHinXsysSYCL/OPHELIE_EDGE_FLUX_STAGE2_CLEANUP_AND_CURSOR_PLAN.md`](../OPHELIE_EDGE_FLUX_STAGE2_CLEANUP_AND_CURSOR_PLAN.md)
> **development log index**：`OPHELIE_PHI_DEVELOPMENT_LOG.md` §4.15–§4.16
> **purpose ofthis document**：before GPT discussion — single-file full handoff——includes code changes, run data, conclusions, open items, questions, upload manifest。

---

## 0. pre-execution baseline（Stage 1 Completed）

before Stage 2 plan，edge-flux Stage 1 already achieved（after power fix）：

| metric | typical value |
|------|--------|
| `production_literature_passed` | 1 |
| `phi_eq_res_vol` | ~8e-5 |
| `edge_res_red` | ~1.1e4 |
| `P_recon_edge` @50 kW | ~50000 W |
| Calibration | `P_recon`（not `P_graph`） |

three-way comparison v2：`discussion_bundle/FRENCH_THREE_WAY_SUMMARY_v2.txt`（H passed=1，A passed=1，G expected failed=0）。

---

## 1. ChatGPT confirmation on Cursor issues (before Stage 2 start)

user raised two points in Cursor，ChatGPT replied as follows（**already used as implementation basis**）：

### 1.1 ParticleDiag disabled by default — **option A confirmed**

- edge-flux under **`EImag/JImag/JouleHeat` = `EEdgeRecon/JEdgeRecon/JouleHeatEdgeRecon`**
- old particle-gradient **only**in `--ophelie-particle-gradient-diagnostics=1` when computing，write to `EImagParticleDiag` etc.
- **keep as-is**

### 1.2 `q_antisymmetry_error ` — **required; insert Step 2.5**

- in scaling regression **before**completed
- neighbor loop compute simultaneously `q_ij`、`q_ji`
- Output：`q_antisym_l1/l2/linf`、`q_antisym_rel_l2`、`q_nonfinite_count`
- initial gate：`q_nonfinite_count=0`，`q_antisym_rel_l2 < 1e-5`（~1e-6）

### 1.3 execution order(ChatGPT confirmed)

1. Step 2.5 — q_antisymmetry diagnostic
2. Step 3 — scaling regression（fixed I + target P）
3. Step 4 — Q spatial soft gate
4. Step 5 — active A field + A_ind use `JEdgeRecon`
5. Step 6 — umbrella/header slimming

### 1.4 explicitly deferred

- A_ind Picard
- thermal coupling
- div-grad residual fix
- compatible/Neumann productionization
- `P_graph` κ Calibration

---

## 2. Cursor Completed code work （per Stage 2 plan items）

### 2.1 Step 1 — Production field unification（Plan §3.1）

**target**：edge-flux underprimary physical quantitiesfrom edge reconstruction，no longer mixed with particle-gradient。

| item | Implementation |
|----|------|
| main E/J/Q | `CopyOphelieEdgeReconToPrimaryEJQCK`：`EImag/JImag/JouleHeat` ← edge recon |
| call site | `runFrenchReducedEmPipeline` in edge post-φ after `syncOphelieEdgeReconToPrimaryEJQ` |
| particle path | default** not run**；CLI `--ophelie-particle-gradient-diagnostics=1` runs grad E/J/Q and copy to `*ParticleDiag` |
| Log | `logOphelieEdgeFluxProductionFieldPolicy()` explicit particle divJ diagnostic only |

**main files**：

- `electromagnetic_ophelie_edge_flux.h` — `CopyOphelieEdgeReconToPrimaryEJQCK`、`CopyOpheliePrimaryEJQToParticlediagnosticCK`
- `electromagnetic_ophelie_french_literature.h` — pipeline branch
- `electromagnetic_ophelie_parameters.h` — `output_particle_gradient_diagnostics_`
- `electromagnetic_ophelie_cli.h` — `--ophelie-particle-gradient-diagnostics=`

**Acceptance**：H production v3/v4 `production_literature_passed=1`；`max_EImag/max_JImag` as edge recon order of magnitude。

---

### 2.2 Step 2 — naming cleanup（Plan §3.3）

| old/ambiguous name | new name | meaning |
|-----------|------|------|
| `JouleHeatEdge`（easily confused） | `JouleHeatEdgeGraph` | graph/Laplace energy density，**not**physicalJoule heat |
| `PowerEdgeParticle` | `PowerEdgeGraphParticle` | particle-graph power sum diagnostic |
| CalibrationPower | `P_recon` / `P_total_recon` | physicalPower，50 kW gate |

**File**：`electromagnetic_ophelie_field_names.h`、`electromagnetic_ophelie_register_fields.h`

---

### 2.3 Step 2.5 — `q_antisymmetry` diagnostic（ChatGPT inserted item）

**physical**：conservation form requires `q_ij = -q_ji`；Implementation on-the-fly must validate pair antisymmetry。

**Implementation**：

- new CK：`ComputeOphelieEdgeFluxQAntisymmetryCK`（neighbor loop internally computes `q_ij`、`q_ji`）
- new per-particle Fields：`EdgeQAntisymMax`、`EdgeQAntisymSqSum`、`EdgeQScaleSqSum`、`EdgeQNeighbor Count`、`EdgeQNonfiniteCount`
- aggregate：`OphelieEdgeFluxQAntisymMetrics` + `evaluateOphelieEdgeFluxQAntisymmetry`
- hook into French pipeline：`runFrenchReducedEmPipeline` post-φ evaluation and `logOphelieEdgeFluxQAntisymMetrics`
- **acceptance gate**：`evaluateEdgeFluxLiteratureAcceptance` add `q_antisym_ok`（`q_nonfinite_count==0` and `q_antisym_rel_l2 < q_antisym_rel_l2_max_`，default `1e-5`）

**File**：`electromagnetic_ophelie_edge_flux.h`、`electromagnetic_ophelie_field_names.h`、`electromagnetic_ophelie_register_fields.h`、`electromagnetic_ophelie_parameters.h`、`electromagnetic_ophelie_french_literature.h`

**Test**：`test_3d_ophelie_edge_flux_sign` add q_antisym gate。

**French H v4 measured**：`q_antisym_rel_l2≈7.1e-8`，`q_nonfinite_count=0`，`q_antisym_gate=1`。

---

### 2.4 Step 3 — Scaling regression（Plan §3.6）

#### 2.4.1 new unit test（uniform-field MMS）

- **path**：`tests/extra_source_and_tests/3d_examples/test_3d_ophelie_edge_flux_scaling/`
- **content**：A×0.5/1/2 → P_recon×0.25/1/4；target 12.5/50/200 kW → P_recon≈target（uniform box，σ=16）
- **build validation**：`build` underbuildRun `edge_flux_scaling_passed=1`

#### 2.4.2 French CLI

- `--coil-current-scale=` → `params.coil_current_scale_`
- `applyOphelieCoilCurrentScale()` in `test_3d_ophelie_french_reduced` should be applied at start
- and existing `--no-literature-calibrate` combined with fixed-I sweep

#### 2.4.3 user French reload run data（2026-06-08）

Particle：`build/reload/Reload.xml`，`n_glass=20500`，`dp=0.02`（user SYCL relax generated）。

| scale | I/loop | P_recon (W) | P/P_ref |
|-------|--------|-------------|---------|
| 1.0 | 1 A | 7.27889 | 1.00 |
| 0.5 | 0.5 A | 1.81972 | **0.250** |
| 2.0 | 2 A | 29.1156 | **4.00** |

A/E/J field peaks scale linearly with I。**B3 blocker closure**。

**Record**：`discussion_bundle/french_H_edge_flux_scaling_reload.log`（Summary）；completeterminalOutputinuserlocal `tee`。

**not yetdo**：French on `--target-power=12500/50000/200000` + calibrate systematic sweep（uniform MMS already covers target-power logic）。

---

### 2.5 Step 4 — Q spatial soft gate（Plan §3.7）

**Implementation**：`computeHostEdgeFluxQSpatialMetrics` + `logOphelieEdgeFluxQSpatialMetrics`

statistics：`Q_min/max/mean`、`Qmax/Qmean`、`Q_nonfinite/negative`、`Q_outer_mean`、`Q_center_mean`、`outer/center`

Gate (parameterized)：

- `q_nonfinite_count=0`，`q_negative_count=0`
- `Qmax > Qmean > 0`
- `1 < Qmax/Qmean < q_spatial_max_over_mean_max_`（default 1e4）
- `outer/center > q_spatial_outer_over_center_min_`（default 1）

hook into `evaluateEdgeFluxLiteratureAcceptance` → `q_spatial_ok`

**French H v4**：`Q_outer/center≈14.14`，`soft_gate=1`

---

### 2.6 Step 5 — Active A field + A_ind source（Plan §3.4–§3.5）

#### 2.6.1 Active A for edge kernels

- `getOphelieActiveARealFieldName(names, params)`：
 - default：`a_coil_real`（Coil Biot， not includes A_ind）
 - `--ophelie-use-a-total-for-edge-flux=1`：`a_src_real`（A_coil + A_ind）
- edge CK constructors add `a_real_field` parameter；all read active A：
 - `ComputeOphelieEdgeFluxPhiRhsFromASrcCK`
 - `ComputeOphelieEdgeFluxResidualCK`
 - `ComputeOphelieEdgeFluxJouleHeatCK`
 - `ReconstructOphelieEdgeFluxElectricCurrentCK`
 - `ComputeOphelieEdgeFluxEdgeDropPairStatsCK`
 - `ComputeOphelieEdgeFluxQAntisymmetryCK`
- `electromagnetic_ophelie_phi.h`：edge RHS pass active A into

#### 2.6.2 A_ind Biot source

- `getOphelieAIndJImagFieldName`：edge-flux → `j_edge_recon`，nothen `j_imag`
- `ComputeOphelieGlassSelfInducedBiotSavartCK` supports optional `j_real_field` / `j_imag_field`
- `runFrenchReducedAIndOneWayDiagnostic`：edge-flux branch uses `solvePhiImagWithCurrentRhs` + `execOphelieEdgeFluxPostPhiPipeline` + `syncOphelieEdgeReconToPrimaryEJQ` + Biot(`JEdgeRecon`)

**⚠️ Importantgap**：`test_3d_ophelie_french_aind_diagnostic.cpp` **not yethook into** `filterOphelieTestCommandLine`，userupload's `--ophelie-current-form=edge-flux` ** ignored**（see §5.2）。

---

### 2.7 Step 6 — Umbrella header slimming (Plan §3.8 partial)

`electromagnetic_ophelie.h` removed from default includes：

- team7 Geometry/Probe
- self-induction
- large phi boundary / vector divergence diagnostics

Changed to comment: include explicitly as needed。already changed tests：

- `test_3d_ophelie_team7.cpp` — explicit include team7 + self_induction
- `test_3d_ophelie.cpp`、`test_3d_ophelie_french_reduced.cpp`、`test_3d_ophelie_french_self_induction_picard.cpp` — explicit include self_induction / aind

**not yetdo**：`archive/` large directory migration; decouple diagnostics from CLI（B5）。

---

### 2.8 document and discussion bundleupdated

| File | content |
|------|------|
| `OPHELIE_PHI_DEVELOPMENT_LOG.md` | §4.15 Stage 2 steps；§4.16 French reload + scaling |
| `discussion_bundle/FRENCH_THREE_WAY_SUMMARY_v3.txt` | H/A/G + scaling Summary |
| `discussion_bundle/french_H_edge_flux_scaling_reload.log` | I² scaling Summary |
| `discussion_bundle/OPHELIE_GPT_DISCUSSION_BLOCKERS.md` | B3 closure |
| `discussion_bundle/README.md` | v3 Index |

**not yetupdated**：`discussion_bundle/OPHELIE_GPT_ROUND_SUMMARY_2026-06-08.md` still defined as "Stage 1 awaiting GPT Stage 2 direction"，**update before GPT upload or attachthis document**。

---

## 3. user run data full record (Stage 2 period)

### 3.1 Particlegenerated（user first time did not see relax in H log because）

- command ：`test_3d_ophelie_french_reduced --relax=1...`（SYCL mesh relax）
- artifact：`build/reload/Reload.xml`，`n_glass≈20500`
- backup：`cp -a reload../tests/.../test_3d_ophelie_french_reduced/reload`
- **subsequent `--reload=1` runs skip relax**，hence log only has `loaded relaxed GlassBody from./reload/Reload.xml`

### 3.2 H — edge-flux production @50 kW（v4）

```bash
$BIN --glass-radius=0.325 --glass-height=0.50 --dp=0.02 \
 --coil-radius=0.40 --coil-num-loops=8 --coil-segments-per-loop=256 \
 --frequency=300000 --sigma=16 --target-power=50000 \
 --reload=1 --literature-mode --state_recording=0 \
 --ophelie-current-form=edge-flux \
 2>&1 | tee../discussion_bundle/french_H_edge_flux_production_v4.log
```

| metric | value |
|------|-----|
| `production_literature_passed` | **1** |
| `P_recon_edge` | 49999.9 W |
| `edge_res_red` | 9347.63 |
| `phi_eq_res_vol` | 9.90e-5 |
| `q_antisym_rel_l2` | 7.15e-8 |
| `Q_outer/center` | 14.14 |
| `P_graph/P_recon` | ~1.15e5 |
| `max_EImag` / `max_JImag` | 336 / 5379 |
| `ampere_turns_eff` | ~663 |

### 3.3 French fixed-I scaling（userSecondtimeCorrectcommand ）

| scale | P_recon | ratio |
|-------|---------|------|
| 1.0 | 7.27889 W | 1 |
| 0.5 | 1.81972 W | 0.25 |
| 2.0 | 29.1156 W | 4.0 |

**mis-run lesson**：once used `$BIN...` placeholder caused no reload、not yet literature-mode、re-ran relax — **ineffective**。

### 3.4 three-way comparison（user reload）

| Case | CLI points | production_lit | passed | notes |
|------|----------|----------------|--------|------|
| **H** | `--ophelie-current-form=edge-flux` | 1 | 1 | primary path |
| **A** | `--ophelie-current-form=particle-gradient --phi-projection-operator=div-grad` | 1 | 1 | `phi_eq_res≈0.427`，divJ≈2.34 |
| **G** | `--phi-projection-operator=edge-flux`（phi-only） | 0 | 0 | divJ≈0.93，**expected** |

### 3.5 H + particle-gradient diagnostic

`--ophelie-particle-gradient-diagnostics=1`：

- `production_literature_passed=1`
- `divJ_L2_red≈0.13`（**diagnostic only**，confirms not suitable as H production gate）
- Log prints `Particle-gradient fields saved to EImagParticleDiag/...`

### 3.6 A_ind one-way（user supplementary run data）

```bash
$BIN --reload=1 --literature-mode --ophelie-current-form=edge-flux \
 2>&1 | tee../discussion_bundle/french_aind_edge_flux_v1.log
```

| Quantity | value |
|----|-----|
| `phi_eq_res_vol` | **0.290** (GMRES stuck at ~0.29, not yet converged) |
| `P_joule_W` | 7.28 |
| `A_ind/A_coil` | 0.433 |
| `B_ind/B_coil` | 0.364 |
| `max_J_imag` | 62.9 |
| `passed` | 1（only checks finiteness） |

**Cursor analysis**: because aind test does not parse `--ophelie-current-form`, actually still **div-grad φ + particle J**, same class as v1 (`A_ind/A_coil≈0.44`). **Cannot** count as edge-flux + `JEdgeRecon` acceptance.

---

## 4. conclusion：edge-flux vs old route

| dimension | H full edge-flux | A div-grad | G phi-only |
|------|------------------|------------|------------|
| production primary path | **is** | fallback | diagnostic only |
| `phi_eq_res` | ~1e-4 | ~0.43 | ~1e-4 |
| `edge_res_red` | ~10³–10⁴ | ~0.02 | ~500 |
| 50 kW Calibration | **P_recon** | P_particle | P_particle |
| legacy divJ gate | not applicable on main path | **passed** | **Failed** |
| q_antisym / Q spatial | **passed** | no | no |
| I² scaling | **passed** | not yet systematically tested | — |

**one sentence**： as OPHELIE production solver，**full edge-flux (H) beats G (phi-only) and A on φ residual, edge conservation, physical power, and scaling G（phi-only） and A φ residual performance**；A can still pass under legacy divJ gate, but is not the main direction。

---

## 5. open items, questions, and risks (for GPT discussion)

### 5.1 already closed（closed this round）

| ID | topic | conclusion |
|----|------|------|
| — | ParticleDiag off by default、option A | ChatGPT confirmed, already implemented |
| — | `P_recon` vs `P_graph` Calibration | only `P_recon` |
| — | `divJ_L2_red` as H production gate | no；diagnostic only |
| B3 | French `I×2 → P_recon×4` | **already accepted**（reload） |

### 5.2 codegap（Cursor questions — **recommended fix before or immediately after GPT**）

| ID | Issue | impact | Recommended |
|----|------|------|------|
| **C1** | `test_3d_ophelie_french_aind_diagnostic` does not parse `--ophelie-current-form` / most ophelie CLI | user thought they ran edge-flux A_ind, actually still div-grad | hook into `filterOphelieTestCommandLine` |
| **C2** | `OPHELIE_GPT_ROUND_SUMMARY_2026-06-08.md` does not yet reflect Stage 2 completed state | GPT may discuss outdated Q1–Q6 | update v2 or attachthis document |
| **C3** | partial `french_*.log` may not be added to repo（`.gitignore`） | discussion bundlemissing attachments | user confirms locally `tee` file exists |

### 5.3 Deferred GPT discussion（blockers still open）

| ID | topic | Status |
|----|------|------|
| **B1** | `P_graph/P_recon ≈ 1.15×10⁵` | permanently diagnostic? study κ? |
| **B2** | H on particle `divJ_L2_red≈0.13`（when enabling diag） | main path already `JEdgeRecon`；divJ always only as A comparison? |
| **B4** | div-grad `phi_eq_res≈0.43–0.49` floor | accept as formulation limit? keep A long-term? |
| **B5** | `archive/` migration；decouple diagnostics from CLI | low-priority engineering debt |

### 5.4 Planinternalnot yetdo（ not GPT blockers, but need scheduling）

| item | Note |
|----|------|
| French target-power sweep | 12.5/50/200 kW + calibrate on reload |
| `use_a_total_for_edge_flux=1` run data | A_ind feedback into edge φ one-way acceptance |
| A_ind Picard | plan explicitly deferred; revisit after E/J/Q + one-way stable |
| thermal coupling、div-grad fix、compatible production | forbidden |

### 5.5 technical questions (for GPT discussion manifest)

1. **A_ind acceptance criteria**：one-way only look at `A_ind/A_coil` ratio enough？shouldneed to `P_recon(A_coil)` vs `P_recon(A_coil+A_ind one-way)`？
2. **`phi_eq_res` in aind test GMRES stuck at 0.29**：is test not using edge RHS, or aind pipeline missing `finalizeOphelieCurrentFor mConfiguration`？
3. **`P_graph` huge**：preserve in any paper/engineering output?？or hide in VTP `JouleHeatEdgeGraph`？
4. **Stage 3 direction**：before Picard, must first run `use_a_total_for_edge_flux` + recompute edge φ one-way report？
5. **TEAM7 / thermal coupling**：when to switch from French reduced to more native geometry？

---

## 6. Modified source file manifest (complete)

### 6.1 core production（edge-flux）

| File | change summary |
|------|----------|
| `tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_edge_flux.h` | primary E/J/Q copy；q_antisym CK/metrics；Q spatial metrics；active A in all edge CKs；helpers |
| `tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_french_literature.h` | pipeline q/Q diagnostic；acceptance gate；`applyOphelieCoilCurrentScale`；`OphelieFrenchEmSolveResult` extended |
| `tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_phi.h` | edge RHS uses `getOphelieActiveARealFieldName` |
| `tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_field_names.h` | renamed + q antisym Fields + particle diag Fields |
| `tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_register_fields.h` | register new fields |
| `tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_parameters.h` | q gate、scaling、active A、particle diag parameter |
| `tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_cli.h` | `--coil-current-scale`、`--ophelie-use-a-total-for-edge-flux`、`--ophelie-particle-gradient-diagnostics` |
| `tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie.h` | umbrella slimming |

### 6.2 Biot / A_ind

| File | change summary |
|------|----------|
| `electromagnetic_ophelie_biot_savart.h` | SelfInduced Biot optional J field names |
| `electromagnetic_ophelie_biot_savart.hpp` | Implementation |
| `diagnostics/electromagnetic_ophelie_aind_diagnostic.h` | edge-flux branch + `JEdgeRecon` source（**test CLI not yet hooked up**） |

### 6.3 Tests and examples

| File | change summary |
|------|----------|
| `test_3d_ophelie_french_reduced/test_3d_ophelie_french_reduced.cpp` | `applyOphelieCoilCurrentScale`；acceptance new parameters；include self_induction |
| `test_3d_ophelie_edge_flux_sign/test_3d_ophelie_edge_flux_sign.cpp` | `a_coil_real` initialization；q_antisym gate |
| `test_3d_ophelie_edge_flux_power_uniform_field/...cpp` | `a_coil_real` |
| `test_3d_ophelie_edge_flux_scaling/` | **new** scaling regression |
| `test_3d_ophelie_team7.cpp` | explicit includes |
| `test_3d_ophelie.cpp` | explicit include self_induction |
| `test_3d_ophelie_french_self_induction_picard.cpp` | explicit includes |

### 6.4 document

| File |
|------|
| `OPHELIE_PHI_DEVELOPMENT_LOG.md` |
| `discussion_bundle/FRENCH_THREE_WAY_SUMMARY_v3.txt` |
| `discussion_bundle/french_H_edge_flux_scaling_reload.log` |
| `discussion_bundle/OPHELIE_GPT_DISCUSSION_BLOCKERS.md` |
| `discussion_bundle/README.md` |
| **this document** `discussion_bundle/OPHELIE_STAGE2_CURSOR_WORK_RECORD.md` |

---

## 7. GPT upload file manifest

### 7.1 required — main discussion documents (choose 1–2)

| Priority | path | Note |
|--------|------|------|
| **P0** | `discussion_bundle/OPHELIE_STAGE2_CURSOR_WORK_RECORD.md` | **this document, full handoff** |
| P0 | `OPHELIE_EDGE_FLUX_STAGE2_CLEANUP_AND_CURSOR_PLAN.md` | ChatGPT original Stage 2 plan |
| P1 | `discussion_bundle/OPHELIE_GPT_DISCUSSION_BLOCKERS.md` | open blocker |
| P1 | `discussion_bundle/FRENCH_THREE_WAY_SUMMARY_v3.txt` | H/A/G + scaling table |
| P2 | `discussion_bundle/OPHELIE_GPT_ROUND_SUMMARY_2026-06-08.md` | Stage 1 round（**slightly outdated**） |
| P2 | `OPHELIE_PHI_DEVELOPMENT_LOG.md` | §4.13–§4.16 long-term log |

### 7.2 required — run logs (user local `tee`; confirm saved before upload)

| Priority | path | Note |
|--------|------|------|
| **P0** | `discussion_bundle/french_H_edge_flux_production_v4.log` | H @50 kW reload **passed=1** |
| **P0** | `discussion_bundle/french_H_edge_flux_scaling_reload.log` | I² scaling Summary |
| P0 | `discussion_bundle/french_A_baseline.log` | A div-grad（user reload） |
| P0 | `discussion_bundle/french_G_edge_flux_phi_only.log` | G expected fail |
| P1 | `discussion_bundle/french_aind_edge_flux_v1.log` | A_ind（**note: CLI not yet effective**） |
| P1 | `discussion_bundle/french_aind_diagnostic_v1.log` | A_ind div-grad v1 baseline |
| P2 | `discussion_bundle/french_H_edge_flux_production_v2.log` | early H after power fix |
| P2 | `discussion_bundle/french_H_edge_flux_production_v3.log` | Cursor aut or un v3 |

### 7.3 recommended upload — CSV diagnostics (optional)

| path | Note |
|------|------|
| `discussion_bundle/ophelie_phi_p0_H_v2.csv` | H phi P0 |
| `discussion_bundle/ophelie_edge_flux_H_v2.csv` | H edge flux diag |
| `discussion_bundle/ophelie_phi_p0_A_v2.csv` | A |
| `discussion_bundle/ophelie_phi_p0_G_v2.csv` | G |

### 7.4 recommended upload — key source files (if GPT needs to change code)

**minimal set**：

```
tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/
 electromagnetic_ophelie_edge_flux.h
 electromagnetic_ophelie_french_literature.h
 electromagnetic_ophelie_parameters.h
 electromagnetic_ophelie_cli.h
 electromagnetic_ophelie_field_names.h
 diagnostics/electromagnetic_ophelie_aind_diagnostic.h
 electromagnetic_ophelie.h
```

**Test**：

```
tests/extra_source_and_tests/3d_examples/
 test_3d_ophelie_french_reduced/test_3d_ophelie_french_reduced.cpp
 test_3d_ophelie_french_aind_diagnostic/test_3d_ophelie_french_aind_diagnostic.cpp
 test_3d_ophelie_edge_flux_scaling/test_3d_ophelie_edge_flux_scaling.cpp
 test_3d_ophelie_edge_flux_sign/test_3d_ophelie_edge_flux_sign.cpp
```

### 7.5 not required (unless GPT asks historical)

- `legacy/div_grad/`、`compatible` related archive
- TEAM7 STL / thermal couplingdocument
- `build/output/*.vtp`（large volume；unless need Q spatial distribution）

### 7.6 reproducecommand quick reference（paste GPT）

```bash
cd ~/SPHinXsysSYCL/build
BIN=./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced
GEO="--glass-radius=0.325 --glass-height=0.50 --dp=0.02 --coil-radius=0.40 --coil-num-loops=8 --coil-segments-per-loop=256 --frequency=300000 --sigma=16"

# relax（one-time）
$BIN $GEO --relax=1 --state_recording=0 --relax-log-every=50

# H production
$BIN $GEO --target-power=50000 --reload=1 --literature-mode --state_recording=0 --ophelie-current-form=edge-flux

# I scaling
$BIN $GEO --reload=1 --literature-mode --no-literature-calibrate --ophelie-current-form=edge-flux --coil-current-scale={0.5,1,2}

# UnitTest
./tests/.../test_3d_ophelie_edge_flux_scaling/bin/test_3d_ophelie_edge_flux_scaling
./tests/.../test_3d_ophelie_edge_flux_sign/bin/test_3d_ophelie_edge_flux_sign
```

---

## 8. recommended GPT discussion opening (can copy)

```text
OPHELIE edge-flux Stage 2 alreadyin Cursor l and ed（see OPHELIE_STAGE2_CURSOR_WORK_RECORD.md）。
French reduced cylindrical glass reload：H production passed=1，I² scaling passed，q_antisym/Q_spatial passed。
please confirm Stage 3 priority：A_ind one-way（fix aind CLI + use_a_total）、B1 P_graph、abandon div-grad fix?。
attachments:this document + v3 Summary + french_H v4 log + scaling log + blockers。
```

---

## 9. version

| Fields | value |
|------|-----|
| Recordauthor | Cursor Agent |
| record date | 2026-06-08 |
| corresponding plan | `OPHELIE_EDGE_FLUX_STAGE2_CLEANUP_AND_CURSOR_PLAN.md` |
| Stage 2 Steps | 1, 2, 2.5, 3, 4, 5（code）, 6（partial） |
| userrun data | H v4, scaling, A/G, particle diag, aind（CLI not yet effective） |
