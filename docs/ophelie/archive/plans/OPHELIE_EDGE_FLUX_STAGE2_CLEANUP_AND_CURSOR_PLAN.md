# OPHELIE Edge-Flux Mainline Cleanup and Stage 2 Correction Development Plan

> For Cursor execution. This document unifies the current OPHELIE edge-flux solver development roadmap, clears ambiguity from legacy div-grad/compatible/Neumann paths, and specifies the code changes that must be made in the next round.

---

## 0. Current Decision Summary

The following is already settled:

1. **edge-flux is the new OPHELIE mainline**. It is not the old `phi-only edge-flux`, but the full SPH neighbor-pair current-flux form:

   ```text
   phi edge equation
   -> edge residual / pairwise current conservation
   -> EEdgeRecon / JEdgeRecon
   -> JouleHeatEdgeRecon
   -> P_recon 50 kW calibration
   -> A_ind = K[JEdgeRecon] one-way / Picard
   ```

2. **div-grad baseline is preserved only as fallback / reference / historical comparison**. Do not keep patching around `phi_eq_res≈0.49`. It is not the primary direction for the current true OPHELIE-like solver.

3. **compatible-div-grad, Neumann, grad projection, D_c/G_c, and related paths are all downgraded to archived diagnostics**. Do not let them appear in production docs, the main umbrella header, the main CLI default flow, or primary acceptance tables.

4. **Legacy particle-gradient E/J/JouleHeat/divJ must not be mixed with edge-flux production**. Under `--ophelie-current-form=edge-flux`, the production primary physical fields must be:

   ```text
   EEdgeRecon
   JEdgeRecon
   JouleHeatEdgeRecon
   P_recon_edge
   EdgeFluxResidual
   ```

5. **`P_graph_edge` is permanently diagnostic**. It must not be called `P_total_edge`, must not be used for 50 kW calibration, and must not participate in production pass/fail. It is only a graph/Laplace energy diagnostic.

6. **Next steps do not add thermal coupling, do not build a full J-phi monolithic solve, and do not continue div-grad residual patching**. Next steps only unify edge-flux production fields, clean up diagnostics, add scaling gates, then run A_ind = K[JEdgeRecon] one-way diagnostic.

---

## 1. Why a Clean Break Is Needed Now

To explore the OPHELIE reduced approach, the codebase previously carried multiple parallel routes:

- `div-grad baseline`: `DivSigmaGrad + DivSigmaA + particle-gradient E/J`
- `compatible-div-grad`: `G_c/D_c`
- `Neumann / boundary correction`
- `grad postprocess projection`
- `legacy phi-only edge-flux`
- new `edge-flux current-form`

These routes were valuable during R&D, but they now cause three serious problems:

### 1.1 Terminology confusion

The same name `edge-flux` can mean:

1. Old `--phi-projection-operator=edge-flux`: replaces only the φ equation, not E/J/Q;
2. New `--ophelie-current-form=edge-flux`: full edge-flux current-flux solver.

This leads to logs like:

```text
phi_eq_res very low, but divJ_L2_red very poor
```

and misinterpretation as edge-flux failure. In reality, the failure is from mixing old particle-gradient postprocessing with the edge-flux equation.

### 1.2 Acceptance metric confusion

`divJ_L2_red` is a particle-gradient route diagnostic, not the edge-flux production current-continuity gate. Edge-flux continuity should use:

```text
edge_res_red
post_phi_edge_res_l2
edge_global_conservation
q_antisymmetry_error
```

### 1.3 Power definition confusion

Two edge-flux power quantities exist:

```text
P_graph_edge = graph/Laplace energy diagnostic
P_recon_edge = physical Joule power from EEdgeRecon/JEdgeRecon
```

Only `P_recon_edge` is physical power and may be used for 50 kW calibration.

---

## 2. Current Edge-Flux Status

From this round's H production v2 log, edge-flux Stage 1 status is:

```text
phi_eq_res_vol              ≈ 8e-5
edge_res_red                ≈ 1.1e4
P_recon_edge                ≈ 50000 W
ampere_turns_eff            ≈ 660 A-turn scale
max_BSrc                    ≈ 1.27e-3 T
production_literature_passed = 1
```

This means:

- the edge-flux φ equation already solves to high accuracy;
- pairwise edge residual is significantly reduced;
- after `P_recon_edge` calibration, coil current, magnetic field, and electric field return to reasonable physical magnitudes;
- Stage 1 is essentially successful.

But key issues remain:

```text
legacy particle-gradient E/J/JouleHeat/divJ still mixed in the pipeline;
A_ind does not yet use JEdgeRecon;
edge-flux field naming and primary output definition / scope are not fully unified.
```

---

## 3. Changes Required Immediately

The following are changes Cursor must execute in the next round. Do not start A_ind Picard first, and do not add thermal coupling.

---

## 3.1 Production field unification: primary E/J/Q in edge-flux mode must come from edge reconstruction

### Current issue

Under `--ophelie-current-form=edge-flux`, the pipeline already computes:

```text
EEdgeRecon
JEdgeRecon
JouleHeatEdgeRecon
P_recon_edge
```

but still runs the legacy path:

```text
execOphelieScalarPhiGradient(...)
EImag = -GradPhi - omega A
JImag = sigma EImag
JouleHeat = 0.5 sigma |EImag|^2
DivJImag = div(JImag)
```

This produces `divJ_L2_red≈0.1` in logs and is mistaken for edge-flux production failure.

### Target

In edge-flux production, primary physical fields must be:

```text
primary E = EEdgeRecon
primary J = JEdgeRecon
primary Q = JouleHeatEdgeRecon
primary P = P_recon_edge
primary continuity = edge_res_red
```

### Recommended implementation

#### Option A: keep dual fields (safest)

Add/confirm fields:

```text
EImagParticle
JImagParticle
JouleHeatParticle
DivJParticle

EEdgeRecon
JEdgeRecon
JouleHeatEdgeRecon
EdgeFluxResidual
```

In edge-flux mode:

```cpp
if (params.ophelie_current_form_ == OphelieCurrentForm::EdgeFlux)
{
    // Production metrics read edge fields.
    result.max_e = maxNorm(EEdgeRecon);
    result.max_j = maxNorm(JEdgeRecon);
    result.joule_power_raw = sum(JouleHeatEdgeRecon * Vol);
    result.continuity_metric = edge_res_red;

    // Particle-gradient fields are only diagnostic.
    if (params.output_particle_gradient_diagnostics_)
    {
        computeParticleGradientEJQDiagnostic();
    }
}
else
{
    // Old fallback mode.
    computeParticleGradientEJQProduction();
}
```

#### Option B: alias legacy fields to edge fields in edge-flux mode

If downstream VTP/thermal modules still only recognize `EImag/JImag/JouleHeat`, copy in edge-flux mode:

```cpp
EImag      = EEdgeRecon;
JImag      = JEdgeRecon;
JouleHeat  = JouleHeatEdgeRecon;
```

Store legacy particle-gradient results separately as:

```text
EImagParticleGradientDiagnostic
JImagParticleGradientDiagnostic
JouleHeatParticleGradientDiagnostic
```

**Recommend Option A plus alias when needed.** Do not leave `EImag/JImag/JouleHeat` ambiguous in edge-flux mode.

---

## 3.2 `divJ_L2_red` is no longer the edge-flux production gate

### Current issue

Legacy `divJ_L2_red` computes divergence of particle-gradient reconstruction, not edge-flux pairwise current conservation.

### Target

Edge-flux production gate becomes:

```text
edge_res_red_l2 > threshold
edge_res_post_phi_l2 finite
edge_global_conservation finite
q_antisymmetry_error small
P_recon_edge target match
fields finite
Q distribution soft gate pass
```

Legacy `divJ_L2_red` is preserved as:

```text
particle_divJ_diagnostic_only
```

### Log must print

```text
[ophelie] edge-flux mode: particle divJ is diagnostic-only and is not a production continuity gate.
```

### Acceptance change recommendation

Old:

```cpp
production_literature_passed = passed && divJ_L2_red >= threshold;
```

New:

```cpp
if (current_form == EdgeFlux)
{
    production_literature_passed =
        edge_res_ok &&
        p_recon_ok &&
        fields_finite &&
        q_distribution_ok;
}
else
{
    production_literature_passed =
        divj_red_ok &&
        p_particle_ok &&
        fields_finite;
}
```

---

## 3.3 `P_graph_edge` is permanently diagnostic; must rename or strongly label

### Current issue

The prior error came from:

```text
P_graph_edge = sum(0.25 * C_ij * edge_drop^2)
```

being mistaken for physical Joule power. Actual results show:

```text
P_graph_edge / P_recon_edge ≈ 1e5
```

This means `C_ij` is a Laplace residual weight, not physical conductance.

### Target

1. Do not use ambiguous names in code and logs:

   ```text
   P_total_edge
   JouleHeatEdge
   ```

2. Use instead:

   ```text
   P_graph_edge
   JouleHeatEdgeGraph
   P_graph_over_recon_warning
   ```

3. Production power uses:

   ```text
   P_recon_edge
   JouleHeatEdgeRecon
   ```

### Not recommended now

Do not calibrate with a single coefficient `kappa`:

```text
P_physical ≈ kappa * P_graph_edge
```

Reason: `P_graph/P_recon` ratio differs greatly between uniform field and French case; it is not a reliable constant.

---

## 3.4 A field in edge-flux must not keep hardcoding `a_src_real`

### Current issue

Edge-flux RHS/residual/reconstruction reads:

```cpp
names.a_src_real
```

directly in many places. Stage 1 with coil source only is fine, but Stage 2 needs:

```text
A_total = A_src + A_ind
```

If edge-drop still hardcodes `A_src`, computed A_ind cannot enter the solver.

### Target

Add active A field helper:

```cpp
inline const std::string& getOphelieActiveARealFieldName(
    const OphelieGlassFieldNames& names,
    const OphelieParameters& params)
{
    if (params.use_a_total_for_edge_flux_)
        return names.a_total_real;
    return names.a_src_real;
}
```

Or more explicitly:

```text
AFieldMode = SourceOnly | SourcePlusInduced
```

All edge-flux related classes must read the active A field:

```text
ComputeOphelieEdgeFluxPhiRhsFromASrcCK       -> rename FromAField
ComputeOphelieEdgeFluxResidualCK             -> active A
ReconstructOphelieEdgeFluxElectricCurrentCK  -> active A
ComputeOphelieEdgeFluxJouleHeatCK            -> active A if needed
```

Recommend renaming the class:

```text
ComputeOphelieEdgeFluxPhiRhsFromAFieldCK
```

Do not keep the name `FromASrc`.

---

## 3.5 A_ind one-way diagnostic must use `JEdgeRecon`

### Current issue

Existing A_ind diagnostic may still default to legacy particle-gradient J.

### New rule

```text
if current_form == edge-flux:
    A_ind source = JEdgeRecon
else:
    A_ind source = JImagParticle / fallback J
```

### Stage 2 one-way diagnostic output

Do not Picard yet; one-way only:

```text
A_ind = K[JEdgeRecon]
A_total = A_src + A_ind
```

Output:

```text
A_ind/A_src
B_ind/B_src
P_recon(A_src only)
P_recon(A_src + A_ind one-way)
edge_res_red(A_src only)
edge_res_red(A_src + A_ind one-way)
Q distribution change
```

### Forbidden

Do not run A_ind Picard until E/J/Q definition / scope is unified.

---

## 3.6 Edge-flux scaling regression must be added

### Why it is needed

Uniform field MMS already validates the local power formula direction, but the French case still needs validation:

```text
multiloop Biot source
reload cylinder particles
GMRES edge-flux φ solve
edge reconstruction
literature current calibration
```

### Required tests

#### Test 1: fixed-current scaling

Example commands:

```bash
test_3d_ophelie_french_reduced \
  --reload=1 \
  --ophelie-current-form=edge-flux \
  --no-literature-current-calibration \
  --coil-current-scale=1.0

test_3d_ophelie_french_reduced \
  --reload=1 \
  --ophelie-current-form=edge-flux \
  --no-literature-current-calibration \
  --coil-current-scale=2.0
```

Check:

```text
A_src, B_src, EEdgeRecon, JEdgeRecon ∝ I
P_recon_edge ∝ I^2
phi ∝ I
edge residual level0 ∝ I
```

#### Test 2: target-power scaling

```bash
--target-power=12500
--target-power=50000
--target-power=200000
```

Check:

```text
ampere_turns_eff ∝ sqrt(P_target)
P_recon_edge ≈ P_target
Qmax ∝ P_target
```

### Acceptance

Add:

```text
edge_flux_scaling_passed
```

Do not force the scaling gate into every individual French run; a standalone regression test is clearer.

---

## 3.7 JouleHeatEdgeRecon spatial distribution soft gate

### Current issue

Total power at 50 kW alone is insufficient; spatial heat distribution must be checked for physical behavior.

### Recommended statistics

Add diagnostics:

```text
Q_edge_recon_min
Q_edge_recon_max
Q_edge_recon_mean
Q_edge_recon_max_over_mean
Q_edge_recon_nonfinite_count
Q_edge_recon_negative_count
Q_edge_recon_outer_mean
Q_edge_recon_center_mean
Q_edge_recon_outer_over_center
```

### Initial soft gate

```text
nonfinite_count = 0
negative_count = 0
Qmax > Qmean > 0
1 < Qmax/Qmean < 1e4
outer_mean / center_mean > 1
```

This gate is for early-stage sanity check only; do not set it too strict initially.

---

## 3.8 Clean up includes and directory structure

### Current issue

Directories are partially organized but need further cleanup.

### Recommended structure

Under:

```text
tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/
```

Keep edge-flux mainline and necessary shared infrastructure; archive the rest.

Recommended layout:

```text
electromagnetic_ophelie/
├── README.md
├── electromagnetic_ophelie.h                    # edge-flux production umbrella
├── electromagnetic_ophelie_edge_flux.h/.hpp     # main production
├── electromagnetic_ophelie_french_literature.h  # French edge-flux pipeline
├── electromagnetic_ophelie_phi.h/.hpp           # scalar solve infrastructure
├── electromagnetic_ophelie_phi_gmres.h
├── electromagnetic_ophelie_phi_krylov_ck.h
├── electromagnetic_ophelie_biot_savart.*
├── electromagnetic_ophelie_multiloop_source.h
├── electromagnetic_ophelie_parameters.h
├── electromagnetic_ophelie_cli.h
├── electromagnetic_ophelie_field_names.h
├── electromagnetic_ophelie_register_fields.h
├── electromagnetic_ophelie_observables.h
├── electromagnetic_ophelie_diagnostics.h
│
├── archive/
│   ├── div_grad_baseline/
│   ├── compatible_divgrad/
│   ├── neumann_boundary/
│   ├── old_phi_only_edge_flux/
│   └── old_docs_index.md
│
├── diagnostics/
│   ├── vector_divergence_mms/
│   ├── rhs_fingerprint/
│   └── operator_linearity/
│
├── benchmarks/
│   └── team7/
│
└── stage2/
    └── self_induction/
```

### Umbrella header rule

`electromagnetic_ophelie.h` must not include all historical routes. It should include only:

```text
production edge-flux
shared source/solver/field infrastructure
French pipeline
explicitly needed diagnostics
```

Historical routes must not be auto-included by the production umbrella.

### Diagnostic headers should not include CLI

If:

```cpp
#include "electromagnetic_ophelie_cli.h"
```

still appears in diagnostic headers, remove it.

Put enum-name helpers in:

```text
electromagnetic_ophelie_parameters.h
```

or create:

```text
electromagnetic_ophelie_enums.h
```

and have both CLI and diagnostics include it.

---

## 4. Recommended documentation updates

### 4.1 `OPHELIE_PHI_DEVELOPMENT_LOG.md`

Add:

```text
§4.15 Edge-flux Stage 2 decision: primary E/J/Q switched to edge reconstruction
```

Content must state:

```text
H v2 confirmed edge-flux Stage 1 passed.
P_recon_edge is the physical calibration power.
P_graph_edge is diagnostic only.
Particle divJ is not production continuity gate in edge-flux mode.
Next step: primary fields switch to EEdgeRecon/JEdgeRecon/JouleHeatEdgeRecon.
```

### 4.2 `README.md` in `electromagnetic_ophelie/`

Must clearly state:

```text
Current production route: edge-flux.
Fallback route: div-grad baseline.
Archived diagnostic routes: compatible-div-grad, Neumann, grad correction, phi-only edge-flux.
```

### 4.3 `discussion_bundle/README.md`

Need to list each round's log:

```text
H production v1: graph power calibration bug
H production v2: P_recon power fix, Stage 1 pass
A baseline: fallback/reference
G phi-only: archived diagnostic
```

---

## 5. Next development order

### Step 1: unify field and acceptance definition / scope

Highest priority.

```text
edge-flux mode uses EEdgeRecon/JEdgeRecon/JouleHeatEdgeRecon as primary.
particle-gradient E/J/Q becomes diagnostic-only.
```

### Step 2: update logs and field naming

```text
P_total_edge -> P_graph_edge
JouleHeatEdge -> JouleHeatEdgeGraph
JouleHeatEdgeRecon -> production Q
```

### Step 3: add scaling regression

```text
fixed current scaling
power target scaling
```

### Step 4: add Q spatial soft gate

```text
outer/center heat ratio
Qmax/Qmean
nonfinite/negative count
```

### Step 5: A_ind one-way using JEdgeRecon

```text
A_ind source = JEdgeRecon
A_total support in edge-flux active A field
no Picard yet
```

---

## 6. Forbidden items

The next round must not do:

```text
1. Do not go back to fix div-grad phi_eq_res≈0.49.
2. Do not continue compatible-div-grad production.
3. Do not continue Neumann / grad projection production.
4. Do not use P_graph_edge for 50 kW calibration.
5. Do not reject edge-flux production using particle divJ_L2_red.
6. Do not run A_ind Picard before primary E/J/Q is unified.
7. Do not add thermal coupling.
8. Do not build a full J-phi monolithic solve.
9. Do not keep ambiguous field names that mix edge-flux and particle-gradient.
```

---

## 7. Direct execution summary for Cursor

```text
Current H edge-flux production v2 passed Stage 1:
- phi_eq_res≈8e-5
- edge_res_red≈1.1e4
- P_recon≈50 kW
- production_literature_passed=1

But before Stage 2, edge-flux production definition / scope must be unified.

Immediate changes:
1. In edge-flux mode, primary E/J/Q use EEdgeRecon/JEdgeRecon/JouleHeatEdgeRecon.
2. Legacy particle-gradient E/J/Q/divJ all downgraded to diagnostic.
3. particle divJ_L2_red no longer the edge-flux production gate.
4. P_graph_edge permanently diagnostic; not used for calibration or pass/fail.
5. 50 kW calibration and production pass use P_recon_edge.
6. All edge-flux RHS/residual/reconstruction stop hardcoding A_src; switch to active A field for A_ind one-way prep.
7. A_ind one-way must use JEdgeRecon as K[J] input.
8. Add French edge-flux scaling regression and Q spatial soft gate.
9. Clean directories and umbrella header; archive legacy routes to avoid continued confusion.
```

---

## 8. Final assessment

Edge-flux has moved from early exploration to the eve of Stage 2:

```text
Stage 1 phi/residual: passed
Stage 1 P_recon power calibration: passed
Stage 1 production acceptance: passed
Stage 2 primary-field consistency: not yet complete
Stage 2 A_ind source consistency: not yet complete
```

The next key step is not chasing lower φ residual, but unifying all edge-flux production physical quantities:

```text
edge phi
edge residual
edge E/J reconstruction
edge JouleHeat
edge A_ind source
```

Only after this unification does subsequent `A_ind = K[JEdgeRecon]` and Picard self-induction have clear physical meaning.
