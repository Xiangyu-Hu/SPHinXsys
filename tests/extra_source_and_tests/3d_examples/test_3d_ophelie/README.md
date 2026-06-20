# test_3d_ophelie

OPHELIE-like particle induction heating on a box glass + coil shell (Stages 0–5).

## Build

```bash
cmake .. -DSPHINXSYS_BUILD_EXTRA_SOURCE_AND_TESTS=ON
ninja test_3d_ophelie
```

## Run

```bash
# Use full path (do not copy literal "...")
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie/bin/test_3d_ophelie --state_recording=0
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie/bin/test_3d_ophelie --state_recording=0 --ophelie-compare-level0
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie/bin/test_3d_ophelie --state_recording=0 --self-induction
```

TEAM7-like annular geometry: see `test_3d_ophelie_team7/`.

## Options (filtered before SPH `handleCommandlineOptions`)

| Flag | Effect |
|------|--------|
| `--no-phi` | Level 0 only (no φ solve) |
| `--ophelie-compare-level0` | Run Level 0 then φ; print delta line |
| `--self-induction` | OPHELIE-like `A_ind` iteration + φ |
| `--phi-solver=PCG\|GMRES\|Jacobi` | Scalar φ solver |
| `--self-induction-relax=<0..1>` | J under-relaxation per outer iter (default 0.15) |
| `--self-induction-max-iter=<N>` | Max outer iterations (default 8) |

## Expected (dp=0.05, default params)

- Baseline φ: `passed=1`, `max_AInd=0`, `P_raw≈41`
- Self-induction: `passed=1`, `self_ind_J_rel<0.05`, `max_AInd>0`, `P_raw` lower before power scaling

VTP under `output/` when `--state_recording=0`.
