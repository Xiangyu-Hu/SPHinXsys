# OPHELIE — TEAM7 Round 3 Diagnostic Closure and GPT Decision Request

> **Purpose**: Round 3 basic physics diagnostics are complete; ask GPT to decide **next operator direction** (generic boundary vs RHS/restore vs diagnostics only).  
> **Date**: 2026-06-12  
> **Repository**: `/home/yyc/SPHinXsysSYCL`  
> **Execution plan**: [`OPHELIE_TEAM7_ROUND3_GAUGE_BOUNDARY_NEXT_CURSOR_PLAN.md`](../OPHELIE_TEAM7_ROUND3_GAUGE_BOUNDARY_NEXT_CURSOR_PLAN.md)  
> **Milestone log**: [`TEAM7_VALIDATION_LOG.md`](../tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7_complex_edge_flux/TEAM7_VALIDATION_LOG.md)  
> **Round 2 discussion pack** (historical): [`OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX.md`](OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX.md)

---

## 0. Five-sentence summary (read first)

1. **Round 3 diagnostic chain closed**: P0 normalization, P1a gauge, P1 rotational-A, P2 hole-wall partition+Jn/Jt, P3 L1 scope, P4 Lenz, P4.1 moment, **P4.5 hole-wall edge-recon** all done; **Picard relax sweep frozen**.
2. **`Bind/Bcoil≈4` is not a Lenz sign error**: `lenz_pass=1`, `corr_phasor(B_ind_i,B_coil_r)≈−0.59`; simple box rotational-A **no 4×** (J preserved ~91%).
3. **Global `C_ij` scaling cannot explain 4×**: moment audit shows `ratio_med/inv(h³)≈4` (same order as kernel 1/h³); **forbid** patching `C_ij` directly with `conductance_ratio≈1.5e8`.
4. **Hole-wall `e_edge_em_mismatch≈0.35` (interior≈0.05) is a boundary EMF closure gap**, but hole-wall `Bind/B≈2.76` **below** interior `≈3.83` — 4× is **not** hole-wall alone; it is **whole-plate J amplitude + non-uniform partition**.
5. **Ask GPT to decide**: next step is (A) generic boundary/hole edge-flux BC, (B) RHS/restore/conductance scaling (not empirical J×constant), (C) close L1 source-only reference first, or (D) reopen physical-scale Picard (A_ind relax) only after the above?

---

## 1. Round 3 completed items (frozen decisions unchanged)

| ID | Content | Conclusion |
|----|------|------|
| P0 | safe_rhs_l2 sweep + solver-local | post-restore insensitive to threshold; solver-local can be default |
| P1a | constant_A_gauge_cancellation | full φ-solve J≈0 is **gauge cancellation**, not failure |
| P1 | rotational_A_uniform_B_phi_solve | J_phi/J_edge≈0.91, **no 4×** |
| P2 | hole_lateral BFS + Jn/Jt | n=199; hole-wall Jn/Jt≈0.51 |
| P3 | L1 documentation | phase0 total ref **≠** pure source-only |
| P4 | A_ind / Lenz audit | lenz_pass=1; in-phase corr≈0 is orthogonal, not bug |
| P4.1 | M_i eigenvalues | aniso~0.06; boundary~0.13; trace and ΣC|r|² consistent |
| P4.5 | hole-wall edge-recon | see §3; **fallback=0**, σE=J exact |
| P5 | Picard relax sweep | **frozen** (cannot fix Bind/B) |

**Still forbidden**: J×0.25, empirical C_ij scaling, cancel φ-solve, edge-only production, a_sign=−1, L2 hard acceptance.

---

## 2. Core numbers (filament, solver-local, source_scale=1, 50 Hz, dp=3 mm)

| Metric | Value | Notes |
|------|-----|------|
| Bind/Bcoil (all) | **3.986** | whole plate |
| phase90 RMS | ~11.8 | log only, not hard gate |
| corr_phasor(B_ind_i,B_coil_r) | **−0.588** | Lenz direction correct |
| conductance ratio median | 1.50×10⁸ | ratio/inv(h³)≈**4.05** |
| moment anisotropy (all) | 0.057 | nearly isotropic |
| rotational-A box J_phi/J_edge | 0.91 | control: operator has no generic 4× |

### 2.1 Partition Bind/B and EMF (P2 + P4.5)

| Partition | n | Bind/B | e_edge_em_mis | ωA/∇φ | fallback | edge_drop_l2 |
|------|---|--------|---------------|-------|----------|--------------|
| interior | 32474 | 3.83 | **0.050** | 1.80 | 0 | 1.51×10⁻⁴ |
| bottom | 6896 | **5.04** | 0.254 | — | 0 | — |
| hole_lateral | 199 | **2.76** | **0.352** | **1.39** | **0** | 2.69×10⁻⁴ |
| all | 49581 | 3.986 | 0.143 | — | 0 | — |

**Interpretation**:
- Global 4× **not** dominated by hole-wall amplitude (hole Bind/B lower).
- Hole-wall **E_edge vs E_em=−∇φ−ωA mismatch** is local maximum (0.35 vs 0.05), but **not** LS fallback, **not** σE≠J (j_sigmaE_mis=0).
- bottom Bind/B highest (5.0) → **partition non-uniformity**; avoid global J scaling.

---

## 3. Ruled out / narrowed issue scope

```text
✗ normalization accidental threshold      → P0 ruled out
✗ φ-solve path broken                     → P1a gauge cancellation ruled out
✗ generic edge-flux 4× on box             → P1 rotational-A ruled out
✗ Lenz / A_ind sign flip                  → P4 ruled out
✗ global C_ij scale clearly wrong         → P4.1 moment ruled out (same order 1/h³)
✗ hole-wall LS fallback dominant          → P4.5 fallback_frac=0 ruled out
✗ σE=J broken on stored fields            → P4.5 j_sigmaE_mis=0 ruled out
✗ Picard relaxation alone                 → P2b/P5 ruled out
```

**Still open**:

```text
? E_edge vs E_em at boundary/hole (especially hole_lateral, bottom)
? whole-plate J_imag amplitude ~4× vs reference (interior also ~3.8× Bind/B)
? L1 coil field vs phase0 total reference mixing
? probe/reference phasor convention (phase90 peak sign)
? need generic boundary conditions (normal current / hole Neumann) vs changing C_ij
```

---

## 4. Six issues for GPT to decide

### Q1: Root-cause ranking for 4×

After ruling out Lenz, global C_ij scale, box generic 4×, does GPT agree the main-cause order is:

1. **Boundary EMF / edge reconstruction** (e_edge_em hole 7× interior)
2. **Whole-plate induction amplitude scale** (interior still Bind/B≈3.8)
3. **Source / probe / reference alignment** (L1 not strict; peak x offset)
4. Picard not converged (**secondary**, cannot fix alone now)

### Q2: Allow changing generic edge-flux operator?

If yes, which **non-empirical** path?

- **A**: hole/outer-surface **BoundaryParticle** edge-drop / LS weight fix (geometry-aware)
- **B**: φ RHS / conductance **continuum-limit scaling** then change `pairwiseNegativeLaplaceWeight` (needs analytic proof)
- **C**: close **L1 source-only reference** first, then judge if 4× is reference mixing
- **D**: diagnostics only; no code change before thermal-flow coupling

**Explicitly forbidden**: `C_ij×k` or `J×0.25` for Bind/B→1.

### Q3: bottom Bind/B≈5.0 vs hole Bind/B≈2.76 split

Support **partition physics** hypothesis (top/bottom near-tangential J, hole-wall normal component significant), hence oppose any **global** amplitude fix?

### Q4: Picard reopen conditions

Re-run physical-scale Picard (A_ind relax) only after **all** of:

```text
boundary edge-recon audit has clear fix direction
rotational-A benchmark still passes
L1 source convention documented
Lenz audit pass (already satisfied)
```

### Q5: L1 reference strategy

Prioritize **f=0 source-only B** or **no-plate σ=0 FEM**, and **forbid** phase0 total as L1 pass until then?

### Q6: TEAM7 acceptance target

Is L2 one-way **permanently** `diagnostic_only=1` until L3 Picard + reference closed? Keep **50 kW power calibration out of L2**?

---

## 5. Recommended upload files (Round 3 closure)

| # | Path |
|---|------|
| 1 | This document |
| 2 | `OPHELIE_TEAM7_ROUND3_GAUGE_BOUNDARY_NEXT_CURSOR_PLAN.md` |
| 3 | `tests/.../TEAM7_VALIDATION_LOG.md` |
| 4 | `discussion_bundle/team7_l2_outputs/team7_edge_partition_audit_one-way_f50_fil_*.csv` |
| 5 | `discussion_bundle/team7_l2_outputs/team7_edge_conductance_audit_one-way_f50_fil_*.csv` |
| 6 | `discussion_bundle/team7_l2_outputs/team7_aind_lenz_audit_one-way_f50_fil_*.csv` |
| 7 | `discussion_bundle/team7_l2_outputs/high_sigma_edge_flux_scaling_P1a_v4_lenz.csv` |
| 8 | `tests/.../HIGH_SIGMA_EDGE_FLUX_P1A_LOG.md` |

**First prompt (copy-paste)**:

```text
Read OPHELIE_GPT_DISCUSSION_TEAM7_ROUND3_CLOSURE.md.
Round 3 diagnostics closed: Lenz/C_ij/box/Picard cannot alone explain Bind/B≈4.
Decide §4 six issues and give Cursor execution order (forbid J×constant, empirical C_ij scaling).
```

---

## 6. Recommended Cursor next steps (after GPT reply)

| Priority | Task | Depends on GPT |
|--------|------|----------|
| — | **Upload §5 bundle, wait for GPT Q1–Q6** | **now** |
| 1 | If C: prepare f=0 / no-plate reference pipeline | Q5 |
| 2 | If A: hole boundary edge LS / normal projection prototype | Q2 |
| 3 | If B: continuum conductance vs discrete C_ij derivation + unit test | Q2 |
| 4 | Re-run Picard A_ind relax | Q4 |

---

*Round 3 closure pack v1 — 2026-06-12*
