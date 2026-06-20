### 2026-06-14 — P5 no-flux boundary (SPHinXsys normal + diagnostic chain)

> Plan: `OPHELIE_TEAM7_NO_FLUX_BOUNDARY_NORMAL_CURSOR_PLAN.md`  
> Intent: **forbid** analytic normals as primary; must use `NormalFromBodyShapeCK` → `NormalDirection` / `SignedDistance`.

#### P5.0 — SPHinXsys normal + relax/reload unified ✅ **success**

**Code**: `electromagnetic_ophelie_team7_boundary_normal.h`; hooked up in `particle_generation_team7.cpp` / `test_3d_ophelie_team7_complex_edge_flux.cpp`.

**Commands**:

```bash
cd build
# relax only + write reload (no EM run)
./test_3d_ophelie_team7_complex_edge_flux --relax=1 --reload=0 --reload-dir=./reload_team7 --native-dp-mm=3
# normal audit
./test_3d_ophelie_team7_complex_edge_flux --reload=1 --reload-dir=./reload_team7 --native-dp-mm=3 \
  --team7-level=coil-only --team7-boundary-normal-audit=1
```

**Acceptance**:

| Item | Result |
|----|------|
| `./reload_team7/Reload.xml` | includes `NormalDirection` + `SignedDistance` |
| Particle count | coil **57947** / plate **49583** |
| `\|n\|` | ≈ **1** (audit CSV) |
| CLI | `--team7-boundary-normal-audit=1` added to `isOphelieTestCommandLineOption` filter |

**Output**: `output/team7_boundary_normal_audit_coil-only_default.csv`

---

#### P5.1 — project-normal post-hoc projection diagnostic ❌ **no-flux OK, EMF closure fail**

**CLI**: `--ophelie-edge-recon-boundary-mode=project-normal` (**diagnostic-only**, does not overwrite production E/J)

Boundary layer: `E ← E - n(n·E)`, `J ← σE`.

| Partition | raw `e_edge_em` | projected | raw `\|Jn\|` → projected |
|------|-----------------|-----------|--------------------------|
| all | 0.137 | **0.142 ↑** | 17729 → 17678 |
| interior | 0.047 | **0.056 ↑** | — |
| hole_lateral | 0.338 | **0.380 ↑** | 484 → **≈0** |
| top/bottom | 0.19–0.23 | slightly worse | → **≈0** |

**Conclusion**: Normal projection **can force n·J≈0**, but **cannot** lower `e_edge_em_mismatch`; slightly worsens → **do not promote to production**.

**Output**: `output/team7_edge_recon_boundary_one-way_*.csv`

---

#### P5.2 — tangent-ls tangent-plane LS diagnostic ❌ **same pattern**

**CLI**: `--ophelie-edge-recon-boundary-mode=tangent-ls`

| Partition | raw | tangent-ls | `\|Jn\|` raw → tangent-ls |
|------|-----|------------|---------------------------|
| all | 0.137 | **0.234 ↑** | 17729 → 17678 |
| interior | 0.047 | **0.174 ↑** | — |
| hole_lateral | 0.338 | **0.443 ↑** | 484 → **≈0** |
| P_recon | 40.3 W | **30.8 W** | fallback_frac=**0** |

**Conclusion**: Tangent-plane LS same: **no-flux effective, EMF worse**; fallback=0 rules out LS ill-conditioning as main cause.

---

#### P5.3 — Boundary operator consistency audit ✅ **end-to-end** (31 pass / 35 fail)

**CLI**: `--team7-boundary-consistency-audit=1` (can run with P5.1/P5.2)

**2026-06-14 run** (L2 one-way, dp=3 mm, scale=0.754, with `tangent-ls`):

| Partition | e_edge_em raw | EMF ωA/∇φ | σE=J | no-flux raw `\|Jn\|/\|Jt\|` | project/tangent no-flux |
|------|---------------|-----------|------|---------------------------|-------------------------|
| interior | **0.047 PASS** | PASS | PASS | **0.47 FAIL** | interior global Jn still large; surface partition → ≈0 |
| top_surface | 0.186 FAIL | PASS | PASS | **0.016 PASS** | project/tangent **≈0 PASS** |
| bottom_surface | 0.235 FAIL | PASS | PASS | **0.021 PASS** | project/tangent **≈0 PASS** |
| hole_lateral | 0.338 FAIL | PASS | PASS | **0.087 PASS** | project/tangent **≈0 PASS** |
| all | 0.137 FAIL | PASS | PASS | **0.352 FAIL** | interior shell dominates global Jn |

**Normal quality** (SPHinXsys n vs partition analytic n):

| Partition | mean angle | max angle |
|------|------------|-----------|
| top/bottom/exterior | **4–6° PASS** | 55–86° FAIL (corners/edges) |
| hole_lateral | 22.7° FAIL | 38.7° PASS |
| interior (shell mislabel) | 85.6° FAIL | 90° FAIL |

**Interpretation**: P5.3 structures P4.5 conclusions — σE=J always closed; **shell partition + interior mislabel** makes raw no-flux fail on all/interior; **outer surface/hole wall** raw no-flux already <0.1. Post-hoc project/tangent on labeled boundary shell passes no-flux but EMF mismatch still FAIL.

**Output**: `output/team7_boundary_consistency_one-way_*.csv`

---

#### P5.4 — Operator-level no-flux ghost edge ⏳ **partial success (no-flux↑, EMF not closed)**

**CLI**: `--ophelie-edge-recon-boundary-mode=no-flux-ghost-edge` (**production**, modifies `ReconstructOphelieEdgeFluxElectricCurrentCK`)

**Mechanism**: On boundary shell particles (`|SignedDistance| ≤ width`), LS reconstruction adds **e_ig=0** ghost edge constraint (`w_ig·n n^T` to M, b unchanged); weight = max neighbor conductance; normal from SPHinXsys `NormalDirection`.

**2026-06-14 run** (L2 one-way, dp=3 mm, scale=0.754, same reload as baseline):

| Metric | baseline (none) | no-flux-ghost-edge | Verdict |
|------|------------------|---------------------|------|
| `e_edge_em` all | 0.137 | **0.138** | ❌ slightly worse |
| `e_edge_em` hole_lateral | 0.338 | **0.351** | ❌ worse |
| `e_edge_em` interior | 0.047 | 0.048 | ≈ unchanged |
| `\|Jn\|/\|Jt\|` top | 0.016 | **0.010** | ✅ improved |
| `\|Jn\|/\|Jt\|` bottom | 0.021 | **0.014** | ✅ improved |
| `\|Jn\|/\|Jt\|` hole_lateral | 0.087 | **0.057** | ✅ improved |
| `\|Jn\|/\|Jt\|` all | 0.352 | 0.351 | ≈ unchanged (interior shell mislabel dominates) |
| `P_recon` | 40.28 W | 40.27 W | ≈ unchanged |
| `Bind/Bcoil` | 3.91 | 3.91 | ≈ unchanged |
| P5.3 pass/fail (no post-hoc items) | — | **22 / 20** | — |

**Conclusion**:

1. **Operator ghost edge lowers outer surface/hole wall `\|Jn\|/\|Jt\|`** on production field (not post-hoc).
2. **Global `e_edge_em_mismatch` not improved**; hole wall slightly worse → ghost constraint alone **insufficient for EMF closure**.
3. Consistent with P5.1/P5.2: **no-flux vs EMF balance is decoupled**; next steps need **φ-RHS / missing moment** (level-set static confinement style), not only E-reconstruction LS patch.
4. **Do not** make ghost edge final production default; keep as opt-in diagnostic/experiment mode.

**Command**:

```bash
cd build
./test_3d_ophelie_team7_complex_edge_flux --reload=1 --reload-dir=./reload_team7 --native-dp-mm=3 \
  --team7-level=one-way --ophelie-edge-recon-boundary-mode=no-flux-ghost-edge \
  --team7-boundary-consistency-audit=1
```

**Output**: `output/team7_boundary_consistency_one-way_*.csv` (production field direct audit)

---

#### P5.5 — φ-RHS ghost + missing moment prototype ❌ **alone ineffective; combined ≈ P5.4**

**New CLI modes** (production; requires `NormalDirection`/`SignedDistance`):

| Mode | Role |
|------|------|
| `no-flux-missing-moment` | E-recon LS adds tangential moment `(I - n n^T)` |
| `no-flux-phi-rhs-ghost` | φ-RHS assembly adds outward ghost pair: `rhs += c_ig ω a·(h n)` |
| `no-flux-full` | ghost edge + missing moment + phi-RHS ghost (all three) |

**Code**: `electromagnetic_ophelie_edge_flux_boundary_closure.h`; changes `ComputeOphelieEdgeFluxPhiRhsFromASrcCK` + `ReconstructOphelieEdgeFluxElectricCurrentCK`.

**2026-06-14 sweep** (L2 one-way, dp=3 mm, scale=0.754, same reload):

| Mode | `e_edge_em` all | hole_lateral `e_edge_em` | hole `\|Jn\|/\|Jt\|` | top `\|Jn\|/\|Jt\|` | P5.3 pass/fail |
|------|-----------------|--------------------------|---------------------|---------------------|----------------|
| none (baseline) | 0.137 | 0.338 | 0.087 | 0.016 | 22 / 20 |
| `missing-moment` | **0.137** | **0.338** | 0.087 | 0.016 | 22 / 20 |
| `phi-rhs-ghost` | **0.137** | **0.338** | 0.087 | 0.016 | 22 / 20 |
| `no-flux-full` | **0.138** | **0.351** | **0.057** ✅ | **0.010** ✅ | 22 / 20 |

**Conclusion**:

1. **Missing moment alone**: no measurable change to `e_edge_em`, Jn/Jt, Bind/B → ❌ fail.
2. **φ-RHS ghost alone**: same; φ residual EMF metric unchanged → ❌ fail.
3. **`no-flux-full`**: effect **entirely from P5.4 ghost edge** (surface Jn/Jt↓, EMF not closed); other two add **no extra benefit**.
4. **TEAM7 hole wall `e_edge_em≈0.35` remains blocking**; next steps need **LHS/φ solver-level** boundary closure or GPT decision to change coupling, not more E/RHS patches.

**Commands**:

```bash
cd build
./test_3d_ophelie_team7_complex_edge_flux ... --ophelie-edge-recon-boundary-mode=no-flux-missing-moment --team7-boundary-consistency-audit=1
./test_3d_ophelie_team7_complex_edge_flux ... --ophelie-edge-recon-boundary-mode=no-flux-phi-rhs-ghost --team7-boundary-consistency-audit=1
./test_3d_ophelie_team7_complex_edge_flux ... --ophelie-edge-recon-boundary-mode=no-flux-full --team7-boundary-consistency-audit=1
```

---

#### P5.4 — Decision (blocks production boundary patch)

1. **Forbid** promoting P5.1/P5.2 post-hoc projection to production; P5.4 ghost edge **opt-in experiment**, not default.
2. **Still forbidden**: empirical C_ij scaling, J×constant, cancel φ-solve, edge-only production, a_sign=−1, Picard relax sweep.
3. **Next steps**: φ-RHS boundary closure + level-set **missing moment** (static confinement category); hole wall `e_edge_em≈0.35` remains top contradiction.
4. SPHinXsys normals **ready** (P5.0); top/bottom/exterior normal quality acceptable; interior shell labels need refinement.

---

## 3. Current acceptance numbers (L2 one-way, dp=3 mm, scale=0.754)

(See full validation log for probe tables.)
