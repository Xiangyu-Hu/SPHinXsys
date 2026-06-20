# test_3d_ophelie_french_aind_diagnostic

One-way **A_ind = K[J_glass]** diagnostic on French reduced geometry (no A_ind feedback).

Pipeline:

```text
multiloop filament Biot -> A_coil
phi solve with A_src = A_coil only
E / J from conductive glass
glass-to-glass Biot -> A_ind, B_ind (JReal/JImag -> matching A/B channels)
report ||A_ind||/||A_coil||, ||B_ind||/||B_coil||
```

## Build & run (cwd = `build/`)

```bash
cmake --build . --target test_3d_ophelie_french_aind_diagnostic -j$(nproc)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_aind_diagnostic/bin/test_3d_ophelie_french_aind_diagnostic --reload=1
```

Requires `build/reload/Reload.xml` from `test_3d_ophelie_french_reduced --relax=1`.
