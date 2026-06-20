# TEAM7 P5 no-flux boundary — GPT discussion bundle manifest

> Packaging date: 2026-06-14  
> Zip: `discussion_bundle/TEAM7_P5_NO_FLUX_BOUNDARY_BUNDLE.zip`

---

## 1. Recommended reading order

| Order | Path | Description |
|------|------|------|
| 1 | `OPHELIE_GPT_DISCUSSION_TEAM7_P5_NO_FLUX_BOUNDARY.md` | **Main document**: summary + pass/fail table + 8 decision issues |
| 2 | `TEAM7_P5_GPT_FIRST_PROMPT.txt` | First prompt to copy-paste to GPT |
| 3 | `p5_logs/TEAM7_VALIDATION_LOG_P5_EXCERPT.md` | Validation log P5 excerpt |
| 4 | `p5_logs/P5_SWEEP_SUMMARY.txt` | P5.5 four-mode sweep grep summary |
| 5 | `p5_outputs/*.csv` | Boundary consistency / audit CSV |
| 6 | `p5_sources/` | Key source copies |
| 7 | `OPHELIE_TEAM7_NO_FLUX_BOUNDARY_NORMAL_CURSOR_PLAN.md` | Original Cursor plan |
| 8 | `OPHELIE_GPT_DISCUSSION_TEAM7_ROUND3_CLOSURE.md` | Round 3 prerequisite (background) |

---

## 2. Case and commands

| Item | Value |
|------|-----|
| Case | `test_3d_ophelie_team7_complex_edge_flux` |
| Reload | `reload_team7/Reload.xml` (dp=3 mm, includes NormalDirection) |
| Typical command | See main document §2 or `p5_logs/P5_SWEEP_SUMMARY.txt` |

```bash
cd build
./test_3d_ophelie_team7_complex_edge_flux \
  --reload=1 --reload-dir=./reload_team7 --native-dp-mm=3 \
  --team7-level=one-way \
  --ophelie-edge-recon-boundary-mode=no-flux-full \
  --team7-boundary-consistency-audit=1
```

---

## 3. `p5_sources/` file list

| File | P5 stage |
|------|---------|
| `electromagnetic_ophelie_team7_boundary_normal.h` | P5.0 |
| `electromagnetic_ophelie_team7_edge_recon_boundary.h` | P5.1/P5.2 |
| `electromagnetic_ophelie_team7_boundary_consistency.h` | P5.3 |
| `electromagnetic_ophelie_edge_flux_boundary_closure.h` | P5.4/P5.5 |
| `electromagnetic_ophelie_edge_flux.h` | Production kernel (large file) |
| `electromagnetic_ophelie_parameters.h` | boundary mode enum |
| `electromagnetic_ophelie_phi_component.h` | φ-RHS hookup |
| `test_3d_ophelie_team7_complex_edge_flux.cpp` | Case main executable |
| `particle_generation_team7.cpp` | relax normals |
| `test_3d_ophelie_team7_complex_edge_flux/README.md` | Case description |

---

## 4. `p5_outputs/`

| File | Description |
|------|------|
| `team7_boundary_consistency_one-way_*.csv` | P5.3/P5.4/P5.5 pass-fail matrix |
| `team7_boundary_normal_audit_*.csv` | P5.0 normal audit (if coil-only run) |

---

## 5. Not included in zip (size / not P5-specific)

- `reload_team7/Reload.xml` (~18 MB) — reproduce locally with `--reload-dir=./reload_team7`
- Full `TEAM7_VALIDATION_LOG.md` — excerpt in `p5_logs/`
- VTP output — ParaView optional
