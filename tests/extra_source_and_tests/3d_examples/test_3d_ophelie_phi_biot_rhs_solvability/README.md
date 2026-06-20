# test_3d_ophelie_phi_biot_rhs_solvability

Compare `div-sigma-a` vs `legacy-flux` RHS under `DivSigmaGrad` LHS (French reduced + multiloop Biot).

## Build

```bash
cd ~/SPHinXsysSYCL/build
cmake --build . --target test_3d_ophelie_phi_biot_rhs_solvability -j$(nproc)
```

## Run (cwd = `build/`)

```bash
BIN=tests/extra_source_and_tests/3d_examples/test_3d_ophelie_phi_biot_rhs_solvability/bin/test_3d_ophelie_phi_biot_rhs_solvability

$BIN
$BIN --reload=1
$BIN --reload=1 --reload-dir=$HOME/SPHinXsysSYCL/build/reload

# Host vs device GMRES dot/norm on full reload (2× GMRES; long)
$BIN --reload=1 --gmres-device-ops-parity
```

`--reload=1` uses `dp=0.02`, GMRES 120 outer. Requires `build/reload/Reload.xml` from:

```bash
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced --relax=1
```

Related: `test_3d_ophelie_phi_rhs_flux_sign_audit` (linear A sign check).
