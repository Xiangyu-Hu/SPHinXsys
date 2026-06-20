# test_3d_ophelie_french_self_induction_picard

Experimental **Picard** loop: `A_total = A_coil + K[J]`, re-solve φ/J each outer step.  
**Not** part of `literature_passed`.

## Run (cwd = `build/`)

```bash
cmake --build . --target test_3d_ophelie_french_self_induction_picard -j$(nproc)

# Lattice smoke (~few minutes)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_self_induction_picard/bin/test_3d_ophelie_french_self_induction_picard --dp=0.08 --self-induction-max-iter=3

# Reload full case (~long: O(N^2) Biot per outer iter × GMRES)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_self_induction_picard/bin/test_3d_ophelie_french_self_induction_picard --reload=1 --self-induction-max-iter=5
```

Prerequisite one-way check: `test_3d_ophelie_french_aind_diagnostic --reload=1`.

**Reload reference (imag-only, legacy):** 2 outer iters, `final_J_rel≈4e-6`, `A_ind/A_coil≈0.371`, `phi_eq_res_vol≈0.574`.

**Complex edge-flux (2026-06-02):**
```bash
./tests/.../test_3d_ophelie_french_self_induction_picard --reload=1 \
  --ophelie-current-form=edge-flux --ophelie-edge-flux-complex=1 --self-induction
```
Expect `self_induction_complex`, `picard_converged=1` (joint **J_rel** + **phi_eq_res**), ~7 iters @ reload.

Tolerances: `--self-induction-tol=0.05` (J), `--self-induction-phi-tol=0.01` (phi_eq_res_vol).

Same Picard path on main case: `test_3d_ophelie_french_reduced --reload=1 --literature-mode --self-induction` (does not set `literature_passed=1`).
