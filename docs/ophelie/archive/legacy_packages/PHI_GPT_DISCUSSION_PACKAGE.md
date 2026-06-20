# OPHELIE φ Operator / RHS Solvability — ChatGPT Discussion Material Pack

> Repository: `/home/yyc/SPHinXsysSYCL`  
> Build: `/home/yyc/SPHinXsysSYCL/build`  
> Main record: [`OPHELIE_PHI_TEST_UPDATE_AFTER_GRADIENT_RESULT.md`](../../OPHELIE_PHI_TEST_UPDATE_AFTER_GRADIENT_RESULT.md) (§6–§9)

---

## 1. Three sentences for ChatGPT to read first (paste at conversation start)

1. **GradPhi sign already ruled out with `test_3d_ophelie_phi_gradient_linear`; do not recommend flipping `(phi_i - phi_j)`.**
2. **Production code unified to `DivSigmaGrad` LHS + `DivSigmaA` RHS; MMS passed; on reload `reconstructed_divJ_red≈1.74`, `literature_passed=1`.**
3. **LegacyFlux RHS and DivSigmaA on Biot fields satisfy `b_legacy ≈ −b_div` (cosine=−1, machine ε); both RHS `eq_res` stall at ~0.57 (reload); `cross_eq_res≈1.79` is the inevitable result of opposite-sign RHS, not independent physics error.**

---

## 2. User measurement summary (2026-06-02, run under `build/`)

### `test_3d_ophelie_phi_rhs_flux_sign_audit`

```
n=125  rhs_cosine_div_legacy=-1  rhs_cosine_div_neg_legacy=1
rhs_div_vs_neg_legacy_vol=1.3e-07  passed=1
```

### `test_3d_ophelie_phi_biot_rhs_solvability` — lattice

```
n=572 dp=0.06  eq_res(div/legacy)=0.331/0.331  cross~1.92
divJ_red(div/legacy)=3.01/0.52  sign_flip_ok=1  passed=1
```

### `test_3d_ophelie_phi_biot_rhs_solvability` — reload

```
n=16066 dp=0.02  eq_res(div/legacy)=0.574/0.574  cross~1.79
divJ_red(div/legacy)=1.74/0.56  sign_flip_ok=1 (neg_legacy~3e-07)  passed=1
```

### `test_3d_ophelie_french_reduced --reload=1 --literature-mode`

```
phi_eq_res_vol≈0.574  reconstructed_divJ_red≈1.74  literature_passed=1
```

---

## 3. Recommended upload: documents / records (priority)

| File | Role |
|------|------|
| [`OPHELIE_PHI_TEST_UPDATE_AFTER_GRADIENT_RESULT.md`](../../OPHELIE_PHI_TEST_UPDATE_AFTER_GRADIENT_RESULT.md) | **Main timeline**: gradient ruled out → Option A → RHS solvability → sign audit |
| [`OPHELIE_PHI_OPERATOR_AUDIT.md`](../../OPHELIE_PHI_OPERATOR_AUDIT.md) | Operator inconsistency root cause and DivSigmaGrad fix |
| [`docs/ophelie/FRENCH_LITERATURE_MODE.md`](../FRENCH_LITERATURE_MODE.md) | Literature acceptance scope, CLI |
| [`docs/ophelie/FRENCH_REDUCED_CASE_ASSUMPTIONS.md`](../FRENCH_REDUCED_CASE_ASSUMPTIONS.md) | Reduced geometry assumptions |
| This file `PHI_GPT_DISCUSSION_PACKAGE.md` | Index + discussion topics |
| (Optional) [`OPHELIE_CODE_REVIEW_NEXT_STEPS.md`](../../OPHELIE_CODE_REVIEW_NEXT_STEPS.md) | Broader task background; φ section may be stale — prefer PHI_TEST_UPDATE |

**Terminal log**: Upload lines **358–634** of `terminals/8.txt` or copy §2 summary.

**CSV (optional)**: `build/output/ophelie_phi_rhs_solvability.csv` (if solvability test was run).

---

## 4. Recommended upload: core source files

Path prefix: `tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/`

| File | Content |
|------|------|
| `electromagnetic_ophelie_parameters.h` | `OpheliePhiLhsOperatorKind`, `OpheliePhiRhsOperatorKind`, GMRES/eq_res params |
| `electromagnetic_ophelie_phi.h` + `electromagnetic_ophelie_phi.hpp` | DivSigmaGrad LHS, DivSigmaA/LegacyFlux RHS, legacy flux kernel |
| `electromagnetic_ophelie_phi_gmres.h` | GMRES, eq_res stopping, preconditioner |
| `electromagnetic_ophelie_phi_solvability.h` | RHS alignment, cross residual, sign-flip metrics |
| `electromagnetic_ophelie_phi_mms_helpers.h` | MMS, `hostPhiEqResVol`, shell stats |
| `electromagnetic_ophelie_laplace.h` | PairwiseLaplace / DivSigmaGrad diagonal |
| `electromagnetic_ophelie_diagnostics.h` | `ComputeOphelieVecdDivergenceCK`, divJ |
| `electromagnetic_ophelie_french_literature.h` | Literature pipeline, acceptance, RHS log |
| `electromagnetic_ophelie_french_reduced_geometry.h` | Cylinder + multiloop Biot geometry |
| `electromagnetic_ophelie_multiloop_source.h` | Biot A_src |
| `electromagnetic_ophelie_cli.h` | `--phi-lhs-operator`, `--phi-rhs-operator`, literature |
| `electromagnetic_ophelie_postprocess.h` | E/J/Q, GradPhi post-processing |

---

## 5. Recommended upload: test cases (diagnostic chain order)

Path prefix: `tests/extra_source_and_tests/3d_examples/`

| Test | Validates |
|------|----------|
| `test_3d_ophelie_phi_gradient_linear` | **Do not change GradPhi sign** |
| `test_3d_ophelie_phi_laplace_rhs_consistency` | MMS: old PairwiseLaplace ≠ RHS |
| `test_3d_ophelie_phi_solve_manufactured_A` | DivSigmaGrad + MMS φ solve |
| `test_3d_ophelie_phi_grad_div_laplace_consistency` | GradPhi+DivJ and Laplace consistency |
| `test_3d_ophelie_phi_rhs_constant_A` | Uniform A → div RHS ≈ 0 |
| `test_3d_ophelie_phi_rhs_flux_sign_audit` | **legacy ≈ −div** (linear A) |
| `test_3d_ophelie_phi_biot_rhs_solvability` | Biot real field div vs legacy + reload |
| `test_3d_ophelie_phi_reduces_divJ` | divJ decreases with φ (simplified field) |
| `test_3d_ophelie_french_reduced` | End-to-end literature / demo |

Each test directory includes `CMakeLists.txt`; solvability / sign_audit include `README.md`.

**Run (cwd = `build/`)**:

```bash
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_phi_rhs_flux_sign_audit/bin/test_3d_ophelie_phi_rhs_flux_sign_audit
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_phi_biot_rhs_solvability/bin/test_3d_ophelie_phi_biot_rhs_solvability
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_phi_biot_rhs_solvability/bin/test_3d_ophelie_phi_biot_rhs_solvability --reload=1
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced --reload=1 --literature-mode --state_recording=0
```

---

## 6. Issues for ChatGPT discussion (copy as prompt)

### A. Confirm our conclusions (do not repeat closed recommendations)

- Should LegacyFlux and DivSigmaA under SPH discretization **satisfy** `b_legacy ≈ −b_div`? Our derivation: div uses `−ω·div(σA)`, legacy uses `−ω σ g·(A_i−A_j)` ≈ `+ω σ div(A)`.
- With `cross_eq_res≈1.79` at `eq_res≈0.57` and `b_legacy≈−b_div`, is it **approximately** `||Lφ−b_div||/||b_legacy|| + ||b_div+b_legacy||/...` scale? Cleaner PDE explanation?

### B. Meaning of eq_res ≈ 0.57 (reload, Biot A)

- With **unified operators** and MMS passed, is GMRES stall at `||Lφ−b||_vol/||b||_vol≈0.57` explained by:
  - (i) Biot discrete A RHS not in div-grad operator range?
  - (ii) Missing weak-form boundary flux terms?
  - (iii) Pure grid/SPH consistency error floor?
- If projecting RHS onto Ran(L) or Galerkin weak-form fix, minimal change path?

### C. Acceptance and delivery (client)

- Current scope: `literature_passed` requires `divJ_L2_red≥1.25`, `phi_eq_res_vol<0.65`, no field scaling; demo uses `--no-phi`.
- Without improving eq_res, is `divJ_red≈1.74` + power calibration enough for French reduced delivery? Missing independent validation?

### D. Do **not** recommend again (unless new evidence)

- Flip `(phi_i - phi_j)` / GradPhi sign;
- Use LegacyFlux as literature RHS;
- Only increase GMRES outer iterations to push eq_res from 0.57 to <0.1 (120+ outer stalled);
- Self-induction, expanded geometry, STL (unrelated to this φ diagnostic chain).

---

## 7. Do not upload (unless ChatGPT asks TEAM7)

- `biot_reference_selection_for_cursor/` (old Biot reference snapshot)
- `build/` VTP, full `Reload.xml` (>1MB)
- TEAM7 STL / `test_3d_ophelie_team7` (separate line)

---

## 8. One-click packaging command (optional)

```bash
cd /home/yyc/SPHinXsysSYCL
zip -r /tmp/ophelie_phi_gpt_package.zip \
  OPHELIE_PHI_TEST_UPDATE_AFTER_GRADIENT_RESULT.md \
  OPHELIE_PHI_OPERATOR_AUDIT.md \
  docs/ophelie/archive/legacy_packages/PHI_GPT_DISCUSSION_PACKAGE.md \
  docs/ophelie/FRENCH_LITERATURE_MODE.md \
  docs/ophelie/FRENCH_REDUCED_CASE_ASSUMPTIONS.md \
  tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_parameters.h \
  tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_phi.h \
  tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_phi.hpp \
  tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_phi_gmres.h \
  tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_phi_solvability.h \
  tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_phi_mms_helpers.h \
  tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_laplace.h \
  tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_diagnostics.h \
  tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_french_literature.h \
  tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_french_reduced_geometry.h \
  tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_multiloop_source.h \
  tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_cli.h \
  tests/extra_source_and_tests/3d_examples/test_3d_ophelie_phi_gradient_linear \
  tests/extra_source_and_tests/3d_examples/test_3d_ophelie_phi_laplace_rhs_consistency \
  tests/extra_source_and_tests/3d_examples/test_3d_ophelie_phi_solve_manufactured_A \
  tests/extra_source_and_tests/3d_examples/test_3d_ophelie_phi_grad_div_laplace_consistency \
  tests/extra_source_and_tests/3d_examples/test_3d_ophelie_phi_rhs_constant_A \
  tests/extra_source_and_tests/3d_examples/test_3d_ophelie_phi_rhs_flux_sign_audit \
  tests/extra_source_and_tests/3d_examples/test_3d_ophelie_phi_biot_rhs_solvability \
  tests/extra_source_and_tests/3d_examples/test_3d_ophelie_phi_reduces_divJ \
  tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/test_3d_ophelie_french_reduced.cpp \
  tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/README.md
```
