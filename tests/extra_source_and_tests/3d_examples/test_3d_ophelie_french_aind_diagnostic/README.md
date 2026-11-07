# test_3d_ophelie_french_aind_diagnostic

One-way **A_ind = K[J_glass]** diagnostic on French reduced geometry.

Pipeline:

```text
multiloop filament Biot -> A_coil
phi solve with A_src = A_coil only
E / J from conductive glass
glass-to-glass Biot -> A_ind, B_ind
optional one-way feedback: A_src = A_coil + A_ind -> full complex recon (real + imag)
report ||A_ind||/||A_coil||, ||B_ind||/||B_coil||, P_coil vs P_total
```

Batch B (P0-03): complex edge-flux feedback **always** runs the real phi/recon chain even if pre-feedback `J_real~0` (nonzero `A_imag` still drives `E_real`). Older Round A logs with `max_J_real_feedback=0` and `P_total=P_coil` are superseded — see `docs/ophelie/25_BATCH_B_...`.

## Build & run (cwd = `build/`)

```bash
cmake --build . --target test_3d_ophelie_french_aind_diagnostic -j$(nproc)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_aind_diagnostic/bin/test_3d_ophelie_french_aind_diagnostic --reload=1
```

Requires `build/reload/Reload.xml` from `test_3d_ophelie_french_reduced --relax=1`.
