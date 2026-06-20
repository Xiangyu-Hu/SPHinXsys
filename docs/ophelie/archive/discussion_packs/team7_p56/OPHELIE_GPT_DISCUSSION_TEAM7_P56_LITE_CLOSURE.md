# TEAM7 P5.6-lite closure — GPT discussion main document

**Date**: 2026-06-02  
**Prerequisite**: `OPHELIE_TEAM7_P5_NO_FLUX_NEXT_CURSOR_PLAN.md` (GPT Cursor execution plan)  
**This round scope**: P5-fix (1–3) + P6a + P6b + TEAM7 boundary A/B acceptance  
**Conclusion upfront**: P5.6-lite **diagnostics effective, L2 closure ineffective**; do not recommend full φ-LHS static confinement on TEAM7 hole geometry.

---

## 1. What we did (timeline)

### 1.1 P5-fix-2: partition fix ✅

**Issue**: Legacy `interior` partition mixed shell particles with `|SignedDistance|≤2dp` (~16248) into “bulk”, causing false-positive `no_flux_jn_over_jt` **~0.47** on `interior`.

**Implementation**:
- `electromagnetic_ophelie_edge_flux_operator_audit.h`: `true_interior` / `boundary_shell_all` / `corner_or_edge`
- `electromagnetic_ophelie_team7_boundary_consistency.h`: new partitions in audit; skip no-flux Jn/Jt on `true_interior` / `all`

**Effect**:
| Partition | `e_edge_em` | `\|Jn\|/\|Jt\|` |
|------|-------------|-----------------|
| `true_interior` | 0.035 ✅ | 0.82 (bulk; no no-flux eval) |
| `boundary_shell_all` | 0.158 | **0.030** ✅ |
| `interior` (legacy) | 0.047 | **0.467** ❌ mislabel |

### 1.2 P5-fix-1: tangent-LS distance_t A/B ✅

**CLI**: `--ophelie-tangent-ls-distance-norm=3d|tangent`

**Effect** (tangent-ls mode, `boundary_shell`):
- `3d`: tangent-ls `e_edge_em` ≈ 0.290
- `tangent`: ≈ **0.168** (improved)
- Hole wall `hole_lateral`: 0.443→0.393 (still FAIL); fair `em_tan` 0.140→0.197
- **Do not promote to production**

### 1.3 P5-fix-3: EMF component diagnostics ✅

`electromagnetic_ophelie_team7_edge_recon_boundary.h`: `em_tangential_mismatch`, `e_edge_em_mismatch_tangent_pair`, `projected_pair`, etc.

Hole example: raw `e_edge_em` 0.338 vs fair `em_tan` 0.140 — tangent-LS worsens global metric but improves tangential fair compare.

### 1.4 P6a: L1 source/reference/probe audit ✅

**New file**: `team7/electromagnetic_ophelie_team7_l1_source_audit.h`  
**Hookup**: `test_3d_ophelie_team7_complex_edge_flux.cpp` (auto on coil-only / one-way)

| Run | rms_rel | corr | best_fit_peak | volume_vs_filament | pass/fail |
|------|---------|------|---------------|-------------------|-----------|
| one-way, volume, scale=0.754 | 0.318 | 0.951 | 1.000 | 0.140 | 20/2 |
| coil-only, filament, scale=1.0 | 0.528 | 0.922 | 0.793 | — | 18/2 |

**fail=2 cause**: `ref_phase0_vs_coil+ind_skin_rms_rel` over threshold (phase0 muFEM reference **not** pure source-only).

**Still open**: peak x_sim≈198–216 mm vs x_ref≈126 mm (geometry/probe segment).

### 1.5 P6b: box/slab boundary sweep ✅

**New file**: `diagnostics/electromagnetic_ophelie_p6b_box_boundary.h`  
**Case**: `test_3d_ophelie_high_sigma_edge_flux_scaling --p6b-boundary-sweep=1`  
**Fix**: `SignedDistance` must use `NormalFromBodyShapeCK` (host) + device sync; cannot manually `registerStateVariableData`.

**Result**: 8/8 pass (4 boundaries × 2 cases)
- constant_A gauge: all 4 boundaries `gauge_pass=1`
- rotational_A: ghost/full slightly worsen `E_rel` (0.304→0.335); `phi-neumann` same as `none`

### 1.6 TEAM7 boundary mode A/B (key acceptance this round)

L2 one-way, dp=3 mm, scale=0.754, f=50 Hz:

| Boundary | Bind/B | e_edge_em all | hole_lateral e_edge_em | hole Jn/Jt | phase90 RMS |
|------|--------|---------------|------------------------|------------|-------------|
| none | 3.913 | 0.137 | **0.338** | 0.087 | 9.80 |
| no-flux-ghost-edge | 3.911 | 0.138 | **0.351** ↑ | 0.057 ↓ | 9.80 |
| no-flux-full | 3.911 | 0.138 | **0.351** ↑ | 0.057 ↓ | 9.80 |

→ Ghost-edge only suppresses surface normal J; **hole-wall EMF slightly worse**; zero L2 hard-metric improvement.

---

## 2. Pass/fail verdict (aligned with GPT plan §3.2)

| GPT plan item | Status | L2 improved? |
|-----------|------|-------------|
| P5-fix-1 tangent distance_t | ✅ done | ❌ hole wall still bad |
| P5-fix-2 partition | ✅ done | ❌ diagnostic only |
| P5-fix-3 EMF components | ✅ done | ❌ diagnostic only |
| P6a L1 audit | ✅ done | ❌ reveals ref not pure source field |
| P6b box benchmark | ✅ done | ❌ no TEAM7 benefit evidence |
| TEAM7 full P5.6 φ-LHS | ❌ **not started** | — |

**P5.6-lite items 1–5 completed**; simple benchmark **did not show** return worth large φ-LHS on TEAM7.

---

## 3. Current production decision (Cursor recommendation; pending GPT)

1. **Default** `--ophelie-edge-recon-boundary-mode=none` (keep)
2. **Do not pursue full φ-LHS static confinement on TEAM7 hole**
3. Keep P5 no-flux modes as **opt-in diagnostics**
4. **Fork**:
   - **A**: Deepen P6a — f=0 / source-only reference; probe/peak geometry audit
   - **B**: July primary path → **French reduced glass delivery**; TEAM7 as benchmark branch

---

## 4. Reproduce commands

```bash
cd build

# P6a one-way
./tests/.../test_3d_ophelie_team7_complex_edge_flux \
  --reload=1 --reload-dir=./reload_team7 --native-dp-mm=3 \
  --team7-level=one-way --team7-coil-source-scale=0.754 \
  --team7-reference-dir=.../reference_data/team7

# P6a filament strict L1
./tests/.../test_3d_ophelie_team7_complex_edge_flux \
  --reload=1 --reload-dir=./reload_team7 --native-dp-mm=3 \
  --team7-level=coil-only --coil-source-model=filament-racetrack \
  --team7-coil-source-scale=1.0 \
  --team7-reference-dir=.../reference_data/team7

# P6b box sweep
./tests/.../test_3d_ophelie_high_sigma_edge_flux_scaling/bin/... \
  --p6b-boundary-sweep=1 --dp=0.04

# TEAM7 boundary A/B
for m in none no-flux-ghost-edge no-flux-full; do
  ./tests/.../test_3d_ophelie_team7_complex_edge_flux \
    --reload=1 --reload-dir=./reload_team7 --native-dp-mm=3 \
    --team7-level=one-way --team7-coil-source-scale=0.754 \
    --ophelie-edge-recon-boundary-mode=$m \
    --team7-boundary-consistency-audit=1 \
    --team7-reference-dir=.../reference_data/team7
done
```

---

## 5. Ten issues for GPT decision

See `TEAM7_P56_LITE_GPT_FIRST_PROMPT.txt` § open issues.

---

## 6. Bundle file index

See `TEAM7_P56_LITE_UPLOAD_MANIFEST.md`.
