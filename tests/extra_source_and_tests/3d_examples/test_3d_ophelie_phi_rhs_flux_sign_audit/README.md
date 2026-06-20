# test_3d_ophelie_phi_rhs_flux_sign_audit

Linear `A(x)` field: verify `legacy-flux RHS ≈ −div(σA) RHS`.

```bash
cd ~/SPHinXsysSYCL/build
cmake --build . --target test_3d_ophelie_phi_rhs_flux_sign_audit -j$(nproc)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_phi_rhs_flux_sign_audit/bin/test_3d_ophelie_phi_rhs_flux_sign_audit
```
