# test_3d_aphi_ck_source_driven_em_solve

Stage 10.14 **P1**: single-body TEAM7-like layout, impressed **current RHS** + GMRES (`lambda_A=off`).

## Run

Build via parent `3d_examples`, then:

```bash
./bin/test_3d_aphi_ck_source_driven_em_solve
```

VTP under `output/`.

## Gate

`passed=1`: GMRES converged, finite fields, conductor J/Joule > thresholds, air Joule small.
