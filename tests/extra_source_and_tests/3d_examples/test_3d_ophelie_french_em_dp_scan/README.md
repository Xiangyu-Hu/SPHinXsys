# test_3d_ophelie_french_em_dp_scan

French reduced **lattice** cases: complex edge-flux EM vs `dp`, optional **Picard A_ind** (independent of `literature_passed` / thermal).

## Run (cwd = `build/`)

```bash
cmake --build . --target test_3d_ophelie_french_em_dp_scan -j$(nproc)

# Smoke — single coarse dp (fast CI)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_em_dp_scan/bin/test_3d_ophelie_french_em_dp_scan

# EM refinement scan: dp ∈ {0.08, 0.06, 0.04}
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_em_dp_scan/bin/test_3d_ophelie_french_em_dp_scan --em-dp-scan

# Picard scan (same dp grid, joint J_rel + phi_eq_res gate)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_em_dp_scan/bin/test_3d_ophelie_french_em_dp_scan --em-dp-scan-picard

# CSV (default ./output/ophelie_french_em_dp_scan.csv)
# --em-dp-scan-csv=path/to.csv
```

Reload-quality check at `dp=0.02` remains on `test_3d_ophelie_french_self_induction_picard --reload=1 --ophelie-current-form=edge-flux --ophelie-edge-flux-complex=1`.

## Gates

| Mode | passed when |
|------|-------------|
| smoke | finite fields; `P_joule > 0`; `phi_eq_res_vol < 0.05` |
| `--em-dp-scan` | smoke + **finest dp** `phi_eq_res_vol < 0.01` (on lattice, φ may not decrease monotonically with dp; reload @0.02 see Picard test) |
| `--em-dp-scan-picard` | smoke + finest dp: `picard_converged` (J_rel + phi_eq_res joint) + `A_ind/A_coil > 0` |

## Related

- One-way A_ind: `test_3d_ophelie_french_aind_diagnostic`
- Picard gate @ reload: `test_3d_ophelie_french_self_induction_picard`
- Picard relax×iter sweep: `test_3d_ophelie_french_self_induction_picard_sweep`
