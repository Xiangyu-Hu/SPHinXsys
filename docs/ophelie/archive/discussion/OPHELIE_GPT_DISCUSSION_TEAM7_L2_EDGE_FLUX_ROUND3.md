# OPHELIE — TEAM7 L2 edge-flux push Round 3 (P0→P4 full record + GPT discussion pack)

> **Purpose**: Upload to ChatGPT; on Round 2 basis discuss **normalization audit, P1a high-σ benchmark, solver-local, physical-scale Picard, P4 operator audit** next steps. 
> **Date**: 2026-06-11 
> **Repository**: `/home/yyc/SPHinXsysSYCL` 
> **Round 2 main document**: [`OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX.md`](OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX.md) 
> **Milestone log**: [`tests/.../TEAM7_VALIDATION_LOG.md`](../tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7_complex_edge_flux/TEAM7_VALIDATION_LOG.md) 
> **Push master plan**: [`OPHELIE_TEAM7_HIGH_SIGMA_EDGE_FLUX_BOTTLENECK_AND_NEXT_PLAN.md`](../OPHELIE_TEAM7_HIGH_SIGMA_EDGE_FLUX_BOTTLENECK_AND_NEXT_PLAN.md) 
> **Upload manifest**: [`TEAM7_L2_UPLOAD_MANIFEST_ROUND3.md`](TEAM7_L2_UPLOAD_MANIFEST_ROUND3.md)

---

## 0. Five-sentence summary for GPT (Round 3)

1. **P0 normalization sweep**: `safe_rhs_l2` from 1e4→1e12 changes `input_scale`, but **post-restore `J`, `Bind/B`, `phase90` unchanged** → `Bind/B≈4` **not** an accidental normalization artifact.
2. **P1a analytic box**: `uniform_E` and `harmonic_A edge_only` **20/20 passed**; `harmonic_A + phi_solve + restore` **0/10 all failed** → **pure edge reconstruction correct**, **φ re-solve + field-scale-restore destroys harmonic A** (same path as TEAM7 L2).
3. **P2β solver-local** (φ PCG internal scaling + edge recon internal scaling, **no A_coil scaling**): **bitwise identical** to `field-scale-restore` (Bind/B=3.987, J_vol=63993, phase90=11.81) → bottleneck **not in restore mechanism**.
4. **P2 Picard physical-scale** (`solver-local`, `relax A_ind=0.05`): phase90 11.8→**7.1**, `Bind/B` still **~4.3**, not converged → Picard **cannot alone close** J amplitude.
5. **P4 operator audit**: `Bind/B` **partition non-uniform** (interior 3.79, bottom 4.99); surface `e_edge_em_mismatch` **3× interior**; `C_ij` conductance ratio median **~1.5×10⁸** → ask GPT: **change C_ij/RHS scale** vs **boundary Neumann/edge-flux** vs **accept discrete ~4× as TEAM7-specific**?

---

## 1. Frozen consensus (Round 2 continues, not overturned)

| Decision | Content |
|------|------|
| phase90 | **not L2 hard gate**; `B_phase90 = −B_ind_imag` |
| L2 one-way | always `diagnostic_only=1` |
| `imag_a_sign` | **frozen +1**; forbid `a_sign=-1` to fix phase90 |
| Main cause (Round 2) | edge-flux **J amplitude** (J-driven), not independent Biot amplification |
| `Bind/B ∝ f` | P3c proved, not ω² runaway |
| feedback | post-restore metrics **identical** to one-way |
| Picard (field-scale) | relax=0.05 only improves phase90, Bind/B ~4 |
| Forbidden | `J×0.25`, TEAM7-specific equations, 50 kW scaling, `a_sign=-1` |

---

## 2. Round 3 added work (GPT plan P0→P4)

### 2.1 P0 — normalization / restore audit (2026-06-11)

**CLI**:
- `--ophelie-edge-flux-normalization-mode=off|field-scale-restore|solver-local`
- `--team7-normalization-sweep=1`

**Results** (filament, 50 Hz, `source_scale=1`):

| `safe_rhs_l2` | `input_scale` | `J_post` | `Bind/B_post` | `phase90_RMS` |
|---------------|---------------|----------|---------------|---------------|
| 1e4 | 0.00413 | 63993 | 3.987 | 11.81 |
| 1e5 | 0.0413 | 63993 | 3.987 | 11.81 |
| 1e6 | 0.413 | 63993 | 3.987 | 11.81 |
| 1e12 | 1.0 | 63993 | 3.987 | 11.81 |

**Conclusion**: Post-restore acceptance quantities **insensitive** to normalization.

**Output**: `team7_normalization_sweep_*.csv`, `logs/team7_run_normalization_sweep.log`

---

### 2.2 P1a — high-σ edge-flux benchmark (2026-06-11)

**New test**: `test_3d_ophelie_high_sigma_edge_flux_scaling` (box dp=40 mm, n=512)

| Case | Result |
|------|------|
| `uniform_E_real` + post-pipeline | **20/20 passed** (E/J rel ~1e-7) |
| `harmonic_A_imag_edge_only` | **20/20 passed** |
| `harmonic_A_imag_phi_solve` + field-scale-restore | **0/10 failed** (E/J rel≈1, P≈0) |

**Conclusion**:
- Edge operator **no systematic amplitude bug** on analytic fields
- **φ re-solve + restore** destroys preset harmonic A (TEAM7 L2 uses same path)

**Output**: `high_sigma_edge_flux_scaling_P1a_v3.csv`, `HIGH_SIGMA_EDGE_FLUX_P1A_LOG.md`, `logs/p1a_high_sigma_sweep_v6.log`

---

### 2.3 P2β — solver-local full implementation (2026-06-11)

**Implementation** (physical `A_coil` not scaled/restored):
1. `phi_rhs` × scale → PCG → `phi`/`phi_rhs` ÷ scale
2. Edge reconstruction: `b_acc` × scale, then `E` ÷ scale
3. Shared `params.edge_flux_solver_local_rhs_scale_`

**TEAM7 one-way comparison** (50 Hz, filament, `source_scale=1`):

| mode | Bind/B | J_imag_L2_vol | phase90 RMS |
|------|--------|---------------|-------------|
| field-scale-restore | 3.987 | 63993 | 11.81 |
| **solver-local** | **3.987** | **63993** | **11.81** |

→ Both normalizations **numerically equivalent**; physical-scale Picard base ready.

---

### 2.4 P2 — physical-scale Picard + A_ind relax (2026-06-11)

**Implementation**: `runTeam7ComplexEdgeFluxPicardWithLog` each iter `combine → solve → Biot → relax A_ind` (no field restore)

**First run** (`solver-local`, `relax=0.05`, 6 iters):

| iter | Aind_rel | Bind/B | Bz_phase90_RMS |
|------|----------|--------|----------------|
| 1 | 0.00 | 3.987 | 11.81 |
| 3 | 0.66 | 4.195 | 11.28 |
| 6 | 0.29 | 4.340 | **7.11** |

- `diagnostic_only=0` (L3 Picard + solver-local)
- phase90 improved; **Bind/B still ~4.3**; not converged

**Output**: `logs/p2_picard_solver_local_r005.log`

---

### 2.5 P4 — edge-flux operator scale audit (2026-06-11)

**New diagnostic**: `electromagnetic_ophelie_edge_flux_operator_audit.h` 
**CLI**: `--team7-operator-audit=1`

**P4.1 global** (n_plate=49580, solver-local, 50 Hz):

| Quantity | Value |
|----|-----|
| `C_ij` [min, max, mean] | 685, 1.26e12, 2.33e11 |
| `sum(C·r²)/(σ·Vol)` [mean, median] | **1.36e8**, **1.50e8** |

**P4.4 partitions** (`skin_h=dp=3 mm`):

| Partition | n | Bind/B | e_edge_em_mis |
|------|---|--------|---------------|
| interior | 31637 | **3.79** | **0.068** |
| top_surface | 8143 | 3.87 | 0.198 |
| bottom_surface | 7993 | **4.99** | **0.239** |
| exterior_lateral | 1807 | 4.19 | 0.202 |
| hole_lateral_proxy | 0 | — | — (heuristic missed) |

**Conclusions**:
1. Bind/B **non-uniform** → boundary effects + global J too large coexist
2. Surface **e_edge_em_mismatch 3×** interior
3. conductance ratio ~1.5e8 needs comparison to continuum conductance

**Output**: `team7_edge_conductance_audit_*.csv`, `team7_edge_partition_audit_*.csv`, `logs/p4_operator_audit_oneway.log`

---

## 3. Core numbers (filament L2, dp=3 mm, 50 Hz, `source_scale=1`)

| Metric | Typical value |
|------|--------|
| Bind/B | ~3.99 |
| phase90 RMS | ~11.8 mT (Picard r=0.05 → ~7.1) |
| J_imag_L2_vol | ~6.4e4 |
| omegaA/gradPhi | ~1.96 |
| e_edge_em_mismatch | ~0.14 global; interior ~0.07; surface ~0.20–0.24 |
| P_recon_W | ~58 W |
| phi_eq_res_imag | ~4.3e-4 |

---

## 4. Issues for GPT (Round 3, 12 items)

### A. Scaling and operator (P4)

1. **P4.1 `C_ij` ratio ~1.5×10⁸**: vs continuum conductance `σ·Vol/h²` or Laplace stencil, is magnitude reasonable? Change `pairwiseNegativeLaplaceWeight` or `C_ij` definition?
2. P1a proves **edge-only correct, phi_solve+restore wrong** — should TEAM7 **skip φ re-solve** (edge recon only) or is **solver-local enough**? Third path?
3. **solver-local ≡ field-scale-restore** — deprecate field-scale-restore as default, keep only solver-local?

### B. Boundary vs global (P4.4)

4. **bottom_surface Bind/B≈5 vs interior≈3.8**: prioritize **boundary edge-flux / Neumann** or global **C_ij scale**?
5. **e_edge_em_mismatch surface 3×**: edge recon **E** vs **−∇φ−ωA** not closed at boundary? Fix priority?
6. **`hole_lateral_proxy` n=0**: provide TEAM7 hole bbox or level-set distance partition?

### C. Picard / acceptance

7. **physical-scale Picard** phase90 7.1 but Bind/B ~4.3 — continue relax sweep (0.02/0.1) or abandon Picard for J amplitude?
8. L3 Picard + solver-local has `diagnostic_only=0` — **when `team7_validation_passed=1`**? Which hard metrics?
9. Round 3 proves restore insensitive — formally require **TEAM7 acceptance on physical-scale (solver-local)**?

### D. Round 2 open issues

10. Round 2: change edge-flux conductance/RHS vs TEAM7-specific coupling — is Round 3 enough to **lock C_ij/RHS change**?
11. **Jey ~10× reference** and **Bind/B≈4** same root cause (J-driven) under P4 partitions?
12. Code priority: **(a) C_ij scale patch (b) boundary Neumann (c) hole partition+boundary (d) dp sweep P1b** — rank them.

---

## 5. Reproduce commands

### L2 one-way + P4 audit (recommended baseline)

```bash
cd /home/yyc/SPHinXsysSYCL/build

./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7_complex_edge_flux/bin/test_3d_ophelie_team7_complex_edge_flux \
 --reload=1 \
 --reload-dir=tests/extra_source_and_tests/3d_examples/particle_generation_TEAM7/bin/reload \
 --team7-level=one-way --native-dp-mm=3 --ophelie-edge-flux-complex=1 \
 --coil-source-model=filament-racetrack --team7-coil-source-scale=1.0 \
 --ophelie-edge-flux-normalization-mode=solver-local \
 --team7-operator-audit=1 --team7-output-tag=p4_audit \
 --team7-reference-dir=tests/extra_source_and_tests/3d_examples/reference_data/team7
```

### P0 normalization sweep

```bash
 --team7-normalization-sweep=1
```

### P1a high-σ benchmark

```bash
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_high_sigma_edge_flux_scaling/bin/test_3d_ophelie_high_sigma_edge_flux_scaling \
 --output-csv=./output/high_sigma_edge_flux_scaling_P1a.csv
```

### L3 Picard physical-scale

```bash
 --team7-level=picard \
 --ophelie-edge-flux-normalization-mode=solver-local \
 --self-induction-relax=0.05 --self-induction-max-iter=8
```

---

## 6. Attachment index (in zip)

| Category | Path |
|------|------|
| Discussion | this doc, `OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX.md` (Round 2) |
| Plan | `OPHELIE_TEAM7_HIGH_SIGMA_EDGE_FLUX_BOTTLENECK_AND_NEXT_PLAN.md` |
| Milestones | `TEAM7_VALIDATION_LOG.md`, `HIGH_SIGMA_EDGE_FLUX_P1A_LOG.md` |
| P0 CSV | `team7_normalization_sweep_*.csv` |
| P1a CSV | `high_sigma_edge_flux_scaling_P1a_v3.csv` |
| P4 CSV | `team7_edge_conductance_audit_*.csv`, `team7_edge_partition_audit_*.csv` |
| Round 2 baseline | `team7_probe_*f50_fil*`, `team7_omega_scaling_*`, `team7_j_vs_b_probe_split_*` |
| Logs | `team7_l2_outputs/logs/` all |
| Source | see `TEAM7_L2_UPLOAD_MANIFEST_ROUND3.md` §C |

---

*Round 3 pack: 2026-06-11*
