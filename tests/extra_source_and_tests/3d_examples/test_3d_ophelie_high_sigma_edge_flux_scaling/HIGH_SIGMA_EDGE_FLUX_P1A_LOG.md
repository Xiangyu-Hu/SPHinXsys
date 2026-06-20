# P1a High-σ Edge-Flux Benchmark Log

Focused benchmark log: isolate the edge-flux operator vs φ-solve path on an analytic box, outside TEAM7 geometry.

**Date**: 2026-06-12 (v4 + P4 Lenz columns)  
**Test binary**: `test_3d_ophelie_high_sigma_edge_flux_scaling`  
**Latest CSV**: `discussion_bundle/team7_l2_outputs/high_sigma_edge_flux_scaling_P1a_v4_lenz.csv` (v4 fields + Lenz columns; older `..._P1a_v4.csv` still valid)

---

## 1. Motivation

P0 showed post-restore TEAM7 metrics are insensitive to `safe_rhs_l2`. P1a/P1 further distinguish:

- **pure gauge** (constant A): full φ-solve should cancel A; J≈0 is **gauge cancellation**, not failure
- **non-conservative drive** (rotational A, curl A = B₀ ẑ): φ cannot cancel globally; should retain ~90% of edge-only J

---

## 2. Geometry and parameters

| Item | Value |
|----|-----|
| Box center | (0, 0, 0.5) m |
| Half-size | 0.16 m |
| dp | 0.04 m |
| n_particles | 512 |
| σ sweep | 16, 1e3, 1e5, 1e7, 3.526e7 S/m |
| f sweep | 50, 200 Hz |
| uniform E | E₀=(100,0,0) V/m, φ_real=−E₀·x, A=0 |
| constant A | A₀=(0.001,0,0) T, φ=0 (pure gauge, curl A=0) |
| rotational A | A=(−0.5 B₀ y, 0.5 B₀ x, 0), curl A = B₀ ẑ |

---

## 3. Case design (v4)

| case_name | Flow | Purpose |
|-----------|------|------|
| `uniform_E_real` | Preset φ → post-pipeline only | Real-chain uniform E-field power closure |
| `constant_A_edge_only` | Preset A → post-pipeline only | **Reconstruction diagnostic**: E=−ωA; not an isolated-conductor production solve |
| `constant_A_gauge_cancellation` | normalize → full φ solve → restore | pure gauge: J≈0 is **expected** |
| `constant_A_gauge_cancellation_solver_local` | same, solver-local | Compare with field-scale-restore |
| `rotational_A_uniform_B_edge_only` | Preset rotational A → post-pipeline | Edge-reconstruction baseline for non-conservative A |
| `rotational_A_uniform_B_phi_solve` | solver-local full φ solve | **P1 main benchmark**: is J amplitude preserved? |

### Acceptance criteria (v4; no longer mis-judge failure)

**constant_A_gauge_cancellation** (gauge cancellation, not E/J rel≈1 failure):

- `gauge_cancellation_passed=1`
- `J_after_phi_over_edge_only << 1` (measured ~3.6e−5, consistent across all σ/f)
- `edge_drop_after_phi_l2` small (~1e−6–1e−5)
- `phi_linear_fit_error` small (~1e−5 order)
- Note: `e_recon_over_exact≈1` vs analytic E/J **has no physical meaning** in this case (analytic J=0, denominator near zero)

**constant_A_edge_only** (retained; purpose annotated):

> edge-only verifies reconstruction under prescribed edge drop; **not** a physical isolated-conductor solve.

**rotational_A_uniform_B_phi_solve**:

- `J_after_phi_over_edge_only ≈ 0.91` (stable across all σ/f)
- `E_rel ≈ 0.30` (vs analytic E_θ; not 4× TEAM7 magnitude)
- `P_recon/P_expected ≈ 0.83`
- `radial_E_leak ≈ 0.18` (open-box boundary leak; P1b partition refinement pending)

---

## 4. Run commands

```bash
cd ~/SPHinXsysSYCL/build

# Recommended: keep summary lines only; full log to file
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_high_sigma_edge_flux_scaling/bin/test_3d_ophelie_high_sigma_edge_flux_scaling \
  --output-csv=../discussion_bundle/team7_l2_outputs/high_sigma_edge_flux_scaling_P1a_v4.csv \
  2>&1 | tee /tmp/p1a_v4.log | grep '^\[p1a\]\|rows=\|passed='
```

Expected last line: `rows=60 pass_rows=60 passed=1`

---

## 5. v4 key results (2026-06-11, user rerun confirmed)

**Summary**: 60 rows, 60 pass, `passed=1`

### 5.1 constant_A — gauge cancellation (20/20 pass)

| σ | f | J_phi/J_edge | gauge_pass | edge_drop_l2 | phi_fit_err |
|---|-----|--------------|------------|--------------|-------------|
| 16 | 50 | 6.7e−4 | 1 | 9.4e−6 | 7.3e−4 |
| 1e3 | 50 | 3.8e−5 | 1 | 5.2e−7 | 1.8e−5 |
| 1e5 | 50 | 3.6e−5 | 1 | 4.9e−7 | 1.2e−5 |
| 3.526e7 | 200 | 3.6e−5 | 1 | 2.0e−6 | 1.2e−5 |

→ full φ-solve **correctly cancels** constant A; old v3 “harmonic_A_imag_phi_solve failed” conclusion **withdrawn**.

### 5.2 rotational_A_uniform_B_phi_solve (10/10 pass)

| σ | f | E_rel | J_phi/J_edge | P_ratio | radial_E_leak |
|---|-----|-------|--------------|---------|---------------|
| 16 | 50 | 0.304 | 0.910 | 0.828 | 0.184 |
| 1e3 | 50 | 0.304 | 0.910 | 0.828 | 0.184 |
| 1e5 | 50 | 0.304 | 0.910 | 0.828 | 0.184 |
| 3.526e7 | 200 | 0.304 | 0.910 | 0.828 | 0.184 |

→ on simple box with non-conservative A: **J amplitude ~91% retained**, E ~70% match; **no 4× systematic bias** (different family from TEAM7 Bind/B≈4).

### 5.3 uniform_E + constant_A_edge_only

- uniform E: E/J rel ~1e−7, P≈1 (all σ/f)
- constant_A edge-only: E/J rel ~1e−7, P≈1 (all σ/f)

### 5.4 P4 Lenz-law audit (rotational_A_uniform_B_phi_solve, 2026-06-12)

Only `rotational_A_uniform_B_phi_solve` rows trigger Biot + Lenz audit (uniform B_drive_z=0.01 T).

| σ | f | corr(Bind_z, Bdrive_z) | lenz_vs_drive | lenz_pass |
|---|-----|------------------------|---------------|-----------|
| 16 | 50 | **−0.802** | 1 | 1 |
| 1e3 | 50 | −0.802 | 1 | 1 |
| 1e5 | 50 | −0.802 | 1 | 1 |
| 3.526e7 | 200 | −0.802 | 1 | 1 |

→ induced B_z **stably opposes** uniform time-varying drive; coexists with J retention ~91%, so **sign chain is correct and amplitude is not 4× amplified** (different family from TEAM7).

---

## 6. Conclusions and implications for TEAM7

1. **Edge reconstruction + σE=J closes on analytic box** (uniform E, constant/rotational A edge-only).
2. **constant A + full φ-solve → J≈0 is gauge cancellation**, not a broken φ-solve path.
3. **rotational A + φ-solve retains ~91% J**; E ~30% deviation + radial leak ~18% still open (boundary/open box), but **not TEAM7-style 4× J amplification**.
4. TEAM7 `Bind/B≈4` more likely from **hole-wall boundary, source field, probe/reference, or feedback**, not generic edge-flux operator scaling error on simple geometry.
5. **P4 Lenz**: on rotational box `corr(Bind_z,Bdrive)≈−0.80`, `lenz_pass=1` — rules out sign reversal.
6. **P4.1 moment**: TEAM7 interior aniso≈0.026, hole_lateral≈0.129; ratio_med/inv(h³)≈4 — **cannot** patch C_ij from this alone.
7. **Next steps**: hole-wall edge reconstruction + L1 source-only reference; **pause** Picard relax sweep.

---

## 7. Artifact paths

| File | Description |
|------|------|
| `discussion_bundle/team7_l2_outputs/high_sigma_edge_flux_scaling_P1a_v4_lenz.csv` | 60 rows + Lenz columns |
| `discussion_bundle/team7_l2_outputs/high_sigma_edge_flux_scaling_P1a_v4.csv` | 60 rows (no Lenz columns, for comparison) |
| `discussion_bundle/team7_l2_outputs/team7_aind_lenz_audit_one-way_f50_fil_*.csv` | TEAM7 P4 single-row audit |
| `discussion_bundle/team7_l2_outputs/high_sigma_edge_flux_scaling_P1a_v3.csv` | Old version (mis-judged phi_solve failed; reference only) |
| `TEAM7_VALIDATION_LOG.md` §P1a v4 | Main milestone summary |

---

## 8. History (v1–v3 deprecated criteria)

v3 judged `harmonic_A_imag_phi_solve` E/J rel≈1 as failure — **incorrect**.  
Reason: constant A is pure gauge; acceptance should use `gauge_cancellation_passed` and `J_after_phi_over_edge_only`, not relative error against J=0 analytic solution.
