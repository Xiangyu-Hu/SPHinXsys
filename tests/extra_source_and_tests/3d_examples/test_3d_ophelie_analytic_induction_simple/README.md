# test_3d_ophelie_analytic_induction_simple

Batch E weak-feedback induction check (doc 23 §14.1).

## Physics

Prescribed rotational gauge for uniform axial `B₀`:

```text
A = ½ B₀ (−y, x, 0)   ⇒   B = B₀ ê_z,  div A = 0
E = −i ω A            ⇒   |E| = ω B₀ r / 2
Q = ½ σ |E|²          ⇒   Q ∝ r²
```

Volume-mean ratio `⟨Q⟩(r>0.75R) / ⟨Q⟩(r<0.25R) = 25` for a right circular cylinder.

**Does not** validate self-consistent `A_ind` / Picard. That needs an independent FEM cylinder (GetDP scaffold under `docs/ophelie/reference/em_rz_getdp/`).

## Build & run

```bash
cd /home/ubuntu/SPHinXsys/build
cmake ..   # once, to pick up the new subdirectory
ninja test_3d_ophelie_analytic_induction_simple

./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_analytic_induction_simple/bin/test_3d_ophelie_analytic_induction_simple \
  --dp=0.05
```

Optional: `--radius=`, `--b0=`, `--sigma=`, `--frequency=`.

Expect `passed=1`, `Q_outer_over_center` near 25 (coarse lattice tolerance 35%), `corr_Q_r2 > 0.85`.
