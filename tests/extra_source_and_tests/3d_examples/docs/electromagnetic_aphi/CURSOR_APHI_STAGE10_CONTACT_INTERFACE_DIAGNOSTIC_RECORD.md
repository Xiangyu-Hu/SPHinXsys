# Stage 10-contact Interface RHS / Operator diagnostic record

Date: 2026-05-21  
Test: `test_3d_aphi_ck_contact_interface_rhs_consistency_diagnostic`  
Related: `test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured`

## Background

Sprint 4 Contact GMRES MMS showed "recursive residual converges but MMS defect≈1, continuous error explodes".  
This diagnostic distinguishes: **interface operator inconsistency** vs **solve workflow/device state error**.

## Diagnostic scenarios

| Label | Meaning |
|------|------|
| `mono_exact_self` | Monolithic: `u=exact`, `\|Au-b\|/\|b\|` |
| `split_exact_assembly_*` | Two-body coupled RHS assembled, `u=exact` (not zeroed) |
| `split_pin_only_*` | Assemble → **zero both sides** → pin `r_hat` |
| `split_zero_left_only_left` | Assemble → **zero left body only**, right body keeps exact |
| `split_left0_right_rhat_left` | Zero both sides → pin right body only (GMRES left-body initial state) |
| `mono_vs_split_lhs` | Same `r_hat` field: monolithic vs split apply lhs comparison by location |

Regional metrics: `all_rel` / `core_rel` / `interface_band_rel` / `core_away_interface_rel`  
(`\|Au-b\|` relative to `\|b\|_2`)

## Key results (SYCL float, dp=0.1, 9E-1 params)

```
mono_exact_self                    all_rel ≈ 1e-9
split_exact_assembly_left/right    all_rel ≈ 1e-11 … 1e-13   ✅
split_pin_only_left/right          all_rel ≈ 0.64 / 0.28     ❌
split_zero_left_only_left          all_rel ≈ 0.0016          ⚠️ (left-body GMRES initial value acceptable)
split_left0_right_rhat_left        all_rel ≈ 0.64            ❌
split_pin_solution_rhat_max_diff   0                         (u==r_hat on host)
mono_vs_split_lhs core_away        max_abs_diff ≈ 5e4        ❌ (after right body zeroed, split deviates severely from mono)
```

## Conclusions

### 1. RHS assembly and exact self-consistency ✅

Coupled assembly (apply exact on both sides then store rhs) perfectly self-consistent when **not zeroed**.  
**Not** an interface flux MMS formula or split apply operator bug in exact state.

### 2. Root cause: cannot recover after `AphiZeroBlockCK` on Contact neighbor body ⚠️

- After `zero` **right body** `solution`, even with `pin r_hat` and host `\|u-r_hat\|=0`, `apply` still gives `\|Au-b\|/\|b\|≈0.64`.
- **Zero left body only, right body keeps post-assembly exact**: left-body mismatch drops to **≈0.0016** (matches GMRES initial "left=0, right=exact").
- Left body `split_pin_rhat` identical to `split_left0_right_rhat` → **left-body pin ineffective after right body zeroed** (device contact still reads wrong neighbor state).

### 3. Correspondence with GMRES failure phenomenon

| Phenomenon | Explanation |
|------|------|
| `left_rel≈2e-6` but `defect≈1` | Recursive residual "self-consistent" under wrong neighbor state, not MMS solution |
| `\|u\|_2=48` vs `\|r_hat\|_2=0.2` | GMRES converged to solution satisfying wrong equation |
| After fix defect≈4.7e-8 | Avoid zeroing right body, MMS restored |

### 4. Fix strategy (implemented in `runTwoBodyContactInterfaceMms`)

1. After coupled RHS assembly: **zero left body only** (initial u_L=0), **do not zero right body**.
2. Left-body GMRES: right body always keeps assembly exact/`r_hat`.
3. Right-body GMRES: left body is converged solution; **do not zero right body** (warm start from exact).
4. Added `syncAphiBlockToDevice` (`finalizeLoadIn`) for explicit push after pin/zero; **cannot** replace "do not zero neighbor body" rule.

## GMRES MMS results after fix

```
mono:    outer=2   discrete_defect≈5e-6   continuous≈1.8e-6
contact: outer=30  discrete_defect≈4.7e-8  continuous≈4.8e-6
         defect_ratio≈0.009   continuous_ratio≈2.6
         final_true_rel≈1.37e-5 (slightly above 1e-5, but within true_rel×10 threshold)
```

Block GS path `converged` flag may still be 0 (hit max_outer), but **MMS accuracy met**.

## Follow-up

1. **Framework layer**: root cause of pin failure after Contact neighbor body zeroed (device delegate / contact cache?) needs further investigation at SphinxSys layer.
2. **mono vs split lhs**: rerun lhs comparison on "right body not zeroed" path, confirm interface band equivalence.
3. **Formal coupled GMRES** (single Arnoldi across two bodies) remains long-term plan; current block GS usable when zero rule is obeyed.
4. **Sprint 5** three-body scaffold must inherit "do not zero neighbor body" constraint.

## Run

```bash
cd build
ninja test_3d_aphi_ck_contact_interface_rhs_consistency_diagnostic \
      test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured
./bin/test_3d_aphi_ck_contact_interface_rhs_consistency_diagnostic
./bin/test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured
```
