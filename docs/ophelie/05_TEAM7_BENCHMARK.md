# TEAM7 Benchmark (Side-Track)

> Consolidated from: all `OPHELIE_TEAM7_*` root plans, `OPHELIE_COMPLEX_EDGE_FLUX_VALIDATION_AND_TEAM7_PLAN.md`, `tests/.../TEAM7_VALIDATION_LOG.md`, `docs/ophelie/TEAM7_GPT_DISCUSSION_PACKAGE.md`.  
> Archive: [`archive/plans/OPHELIE_TEAM7_*`](archive/plans/), [`archive/discussion/`](archive/discussion/), [`archive/discussion_packs/`](archive/discussion_packs/)

**Status:** Side-track benchmark — **not** midterm delivery gate. RH200 is the production demo path.

---

## 1. Purpose

Validate edge-flux complex solver on **high-conductivity metal** (TEAM7 workshop geometry) with reference B/J probes. Exposed scaling and boundary issues not visible in French glass (σ ≈ 16 S/m).

---

## 2. Workstream summary

| Workstream | Plan file (archive) | Focus |
|------------|---------------------|-------|
| L2 probe validation | `OPHELIE_TEAM7_L2_VALIDATION_DEBUG_AND_CURSOR_PLAN.md` | Bz, Jey vs reference CSV |
| High-σ scaling | `OPHELIE_TEAM7_HIGH_SIGMA_EDGE_FLUX_BOTTLENECK_AND_NEXT_PLAN.md` | P1a scaling sweeps |
| Round 3 gauge/boundary | `OPHELIE_TEAM7_ROUND3_GAUGE_BOUNDARY_NEXT_CURSOR_PLAN.md` | Gauge + boundary closure |
| P5 no-flux normal | `OPHELIE_TEAM7_P5_NO_FLUX_NEXT_CURSOR_PLAN.md` | Normal/tangent boundary |
| P5 normal cursor detail | `OPHELIE_TEAM7_NO_FLUX_BOUNDARY_NORMAL_CURSOR_PLAN.md` | tangent-LS distance, partitions |
| Coil racetrack source | `OPHELIE_TEAM7_COIL_SOURCE_RACETRACK_PLAN.md` | Racetrack coil geometry |
| Review / action | `OPHELIE_TEAM7_REVIEW_AND_ACTION_PLAN.md` | Cross-cutting review |

---

## 3. GPT discussion rounds (merged index)

| Round | Date | Main archive doc |
|-------|------|------------------|
| L2 Round 2 | 2026-06 | `archive/discussion/OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX.md` |
| L2 Round 3 | 2026-06-11 | `archive/discussion/OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX_ROUND3.md` |
| Round 3 closure | 2026-06 | `archive/discussion/OPHELIE_GPT_DISCUSSION_TEAM7_ROUND3_CLOSURE.md` |
| P5 no-flux | 2026-06-14 | `archive/discussion/OPHELIE_GPT_DISCUSSION_TEAM7_P5_NO_FLUX_BOUNDARY.md` |
| P56 lite closure | 2026-06 | `archive/discussion_packs/OPHELIE_GPT_DISCUSSION_TEAM7_P56_LITE_CLOSURE.md` |

Upload manifests (paths + zip commands): `archive/discussion/TEAM7_L2_UPLOAD_MANIFEST*.md`, `archive/discussion/TEAM7_P5_UPLOAD_MANIFEST.md`.

---

## 4. Key outputs (CSV)

Location: `discussion_bundle/team7_l2_outputs/`

| File pattern | Content |
|--------------|---------|
| `team7_probe_one-way*.csv` | Probe B/T/J one-way runs |
| `team7_bz_one-way*.csv` | Bz vs reference |
| `team7_jey_probe*.csv` | Surface Jey |
| `high_sigma_edge_flux_scaling_P1a*.csv` | High-σ scaling |
| `team7_edge_partition_audit*.csv` | Edge partition diagnostics |
| `TEAM7_*_reference_*.csv` | Reference data |

Reference data: `tests/extra_source_and_tests/3d_examples/reference_data/team7/` (see [`docs/TEAM7-reference/README.md`](../TEAM7-reference/README.md)).

---

## 5. P5 fixes (completed in code)

- **P5-fix-1:** tangent-LS distance A/B CLI + kernel
- **P5-fix-2:** true_interior / boundary_shell / corner partitioning
- **P5-fix-3:** normal/tangent EMF mismatch diagnostics

Validation log: `tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7_complex_edge_flux/TEAM7_VALIDATION_LOG.md` (1297 lines — kept at case path, not duplicated).

---

## 6. Executable

```text
test_3d_ophelie_team7_complex_edge_flux   # primary TEAM7 case
test_3d_ophelie_team7                     # legacy entry
test_3d_ophelie_high_sigma_edge_flux_scaling
```

Smoke archive: [`TEAM7_SOURCE_SMOKE_ARCHIVE.md`](TEAM7_SOURCE_SMOKE_ARCHIVE.md).

---

## 7. Known limitations

- Strict L2 metrics not at delivery quality for client sign-off
- No-flux boundary still under active diagnostic (P5/P56 packs)
- Not coupled to RH200 stirring midterm

Do **not** block RH200 midterm on TEAM7 closure.
