# OPHELIE — TEAM7 P5 no-flux boundary closure and GPT decision request

> **Purpose**: P5.0–P5.5 systematically run; ask GPT to decide **continue φ-LHS boundary closure** or **pivot to L1 reference / French primary path**.  
> **Date**: 2026-06-14  
> **Repository**: `/home/yyc/SPHinXsysSYCL`  
> **Original plan**: [`OPHELIE_TEAM7_NO_FLUX_BOUNDARY_NORMAL_CURSOR_PLAN.md`](../OPHELIE_TEAM7_NO_FLUX_BOUNDARY_NORMAL_CURSOR_PLAN.md)  
> **Milestone log**: [`TEAM7_VALIDATION_LOG.md`](../tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7_complex_edge_flux/TEAM7_VALIDATION_LOG.md) §P5  
> **Round 3 prerequisite**: [`OPHELIE_GPT_DISCUSSION_TEAM7_ROUND3_CLOSURE.md`](OPHELIE_GPT_DISCUSSION_TEAM7_ROUND3_CLOSURE.md)

---

## 0. Five-sentence summary (read first)

1. **SPHinXsys normals hooked up** (P5.0): `NormalFromBodyShapeCK` → reload with `NormalDirection` + `SignedDistance`; coil 57947 / plate 49583 particles, `|n|≈1`.
2. **Post-hoc no-flux projection ineffective for EMF closure** (P5.1/P5.2): `project-normal` / `tangent-ls` both drive **`|Jn|→0`**, but **`e_edge_em_mismatch` rises 0.137 → 0.142 / 0.234** → forbid as production.
3. **Operator-level ghost edge partially succeeds** (P5.4): production `no-flux-ghost-edge` lowers top/bottom/hole **`|Jn|/|Jt|`**, but **global `e_edge_em` unchanged or slightly worse** (0.137→0.138; hole 0.338→0.351).
4. **φ-RHS ghost and missing moment alone ineffective** (P5.5): identical to baseline; `no-flux-full` **≈ P5.4 ghost edge alone**.
5. **Ask GPT to decide**: by end of June continue **P5.6 φ-LHS level-set static confinement**, or **downgrade TEAM7 to benchmark branch, pivot primary path to French reduced**? (see §4)

---

## 1. P5 stage table (successes and failures recorded)

| ID | Content | CLI / mechanism | Result | Verdict |
|----|------|------------|------|------|
| **P5.0** | SPHinXsys normals + relax/reload | `--relax=1` + `NormalFromBodyShapeCK` | reload has normals; audit `\|n\|≈1` | ✅ **success** |
| **P5.1** | post-hoc project-normal | `--ophelie-edge-recon-boundary-mode=project-normal` | `\|Jn\|→0`; `e_edge_em` 0.137→**0.142** | ❌ **no-flux OK, EMF fail** |
| **P5.2** | post-hoc tangent-ls | `--ophelie-edge-recon-boundary-mode=tangent-ls` | `\|Jn\|→0`; `e_edge_em` 0.137→**0.234**; `P_recon` 40→31 W | ❌ **same** |
| **P5.3** | boundary consistency audit | `--team7-boundary-consistency-audit=1` | 31 pass / 35 fail (incl. post-hoc); σE=J all PASS | ✅ **working** |
| **P5.4** | operator ghost edge (production) | `no-flux-ghost-edge` | surface `\|Jn\|/|Jt|`↓; `e_edge_em` not improved | ⏳ **partial** |
| **P5.5a** | missing moment (production) | `no-flux-missing-moment` | identical to baseline | ❌ **fail** |
| **P5.5b** | φ-RHS ghost (production) | `no-flux-phi-rhs-ghost` | identical to baseline | ❌ **fail** |
| **P5.5 full** | all three combined | `no-flux-full` | ≈ P5.4; no extra EMF gain | ❌ **EMF not closed** |

**Still forbidden** (not violated in P5): empirical C_ij scaling, J×constant, cancel φ-solve, edge-only production, a_sign=−1, Picard relax sweep, conductive dummy params and Biot/Joule.

---

## 2. Core numbers (L2 one-way, dp=3 mm, scale=0.754, volume-racetrack, field-scale-restore)

### 2.1 Global (production fields, no post-hoc)

| Mode | `e_edge_em` all | hole_lateral `e_edge_em` | `Bind/Bcoil` | `P_recon` |
|------|-----------------|--------------------------|--------------|-----------|
| baseline (none) | **0.137** | **0.338** | 3.91 | 40.28 W |
| ghost edge (P5.4) | 0.138 | 0.351 | 3.91 | 40.27 W |
| missing moment (P5.5a) | **0.137** | **0.338** | 3.91 | 40.28 W |
| phi-RHS ghost (P5.5b) | **0.137** | **0.338** | 3.91 | 40.28 W |
| no-flux-full (P5.5) | 0.138 | 0.351 | 3.91 | 40.27 W |

### 2.2 Normal no-flux diagnostic `|Jn|/|Jt|` (SPHinXsys n, raw production field)

| Partition | baseline | ghost edge / full |
|------|----------|-------------------|
| top_surface | 0.016 | **0.010** |
| bottom_surface | 0.021 | **0.014** |
| hole_lateral | 0.087 | **0.057** |
| all | 0.352 | 0.351 |
| interior (shell mislabel) | 0.467 FAIL | 0.466 FAIL |

### 2.3 P5.1/P5.2 post-hoc diagnostics (production unchanged)

| Partition | raw `e_edge_em` | project-normal | tangent-ls |
|------|-----------------|----------------|------------|
| all | 0.137 | **0.142 ↑** | **0.234 ↑** |
| hole_lateral | 0.338 | **0.380 ↑** | **0.443 ↑** |
| hole `\|Jn\|` | 484 | **≈0** | **≈0** |

**Interpretation**:

- **no-flux (n·J=0) and EMF balance (E_edge ≈ −∇φ−ωA) are decoupled in this discretization**; forcing no-flux alone cannot fix `e_edge_em`.
- Hole-wall **`e_edge_em≈0.35`** vs interior **≈0.05** remains Round 3 top contradiction; P5 did not shrink that gap.
- **Bind/B≈4** unchanged across P5 modes → boundary patch did not touch probe amplitude root cause.

---

## 3. Implemented code (Cursor side)

| File | Role |
|------|------|
| `team7/electromagnetic_ophelie_team7_boundary_normal.h` | P5.0 normal audit / reload |
| `team7/electromagnetic_ophelie_team7_edge_recon_boundary.h` | P5.1/P5.2 post-hoc diagnostic CK |
| `team7/electromagnetic_ophelie_team7_boundary_consistency.h` | P5.3 partition pass/fail audit |
| `electromagnetic_ophelie_edge_flux_boundary_closure.h` | P5.4/P5.5 closure flags |
| `electromagnetic_ophelie_edge_flux.h` | ghost edge + missing moment + φ-RHS ghost kernels |
| `test_3d_ophelie_team7_complex_edge_flux.cpp` | CLI hookup |
| `particle_generation_team7.cpp` | relax writes normals |

**Case**: `test_3d_ophelie_team7_complex_edge_flux`; reload `./reload_team7/Reload.xml` (dp=3 mm).

---

## 4. Eight issues for GPT to decide

### Q1: End-of-June strategy

Does P5 conclusion support **downgrading TEAM7** (high-σ benchmark only), **primary path to French reduced glass**? Or **1–2 weeks** worth P5.6 (φ-LHS static confinement)?

### Q2: Is no-flux BC the correct physics target?

At conductor–air interface, should edge-flux use **`n·J=0`**? Or is mismatch from **gauge / 2D hole geometry / missing displacement**, so changing BC is wrong direction?

### Q3: Continue P5.6 (φ-LHS level)?

If yes, which **non-empirical** path to prioritize?

- **A**: pairwise Laplacian LHS mirror ghost (`φ_g=φ_i` adiabatic analogy)
- **B**: level-set missing support adds **LHS matrix moment** (not just E-recon M)
- **C**: hole **Neumann on φ** explicit: `n·∇φ = −ωA·n`
- **D**: stop boundary operator changes; pivot to **L1 source-only reference** audit

### Q4: Keep P5.4 ghost edge?

Production `no-flux-ghost-edge` improves surface `\|Jn\|/|Jt|` but worsens hole `e_edge_em`. **Keep opt-in experiment** or **rollback**?

### Q5: Interior shell partition mislabel

P5.3 shows interior shell **normal angle 85°**, raw `\|Jn\|/|Jt\|≈0.47`. Refine partition labels before evaluating BC? (else all/interior metrics dominated by mislabeled particles)

### Q6: Hole 0.35 vs interior 0.05

Must hole gap close before L2 smoke? Or acceptable as **boundary-layer artifact**, close **whole-plate Bind/B 4×** (amplitude/source/reference) first?

### Q7: Relation to Round 3 Q2

Round 3 asked generic operator vs RHS/restore vs L1. P5 tried **E/RHS-level generic BC**. Agree **next generic code change worth trying is φ-LHS**? Else pivot to reference/coupling?

### Q8: Acceptance gates

If continuing TEAM7, add **hard gates** (e.g. hole_lateral `e_edge_em<0.1` or `\|Jn\|/|Jt\|<0.1`), or stay **smoke-only** until French primary path ready?

---

## 5. Cursor-side recommendations (not conclusions; for GPT to challenge)

1. **P5.6 low ROI**: E-recon + φ-RHS first-order patches exhausted; φ-LHS large scope; **French primary path** may be better value by end of June.
2. **If GPT chooses continue TEAM7**: prioritize **C (φ Neumann)** or **A (Laplacian mirror)**, validate on **box benchmark** before TEAM7 hole.
3. **If GPT chooses downgrade TEAM7**: keep P5 audit CLI as regression, **freeze** production boundary mode to `none`.

---

## 6. Attachment manifest

See same directory [`TEAM7_P5_UPLOAD_MANIFEST.md`](TEAM7_P5_UPLOAD_MANIFEST.md) and zip `p5_outputs/`, `p5_sources/`.

**Recommended GPT first prompt**: [`TEAM7_P5_GPT_FIRST_PROMPT.txt`](TEAM7_P5_GPT_FIRST_PROMPT.txt)
