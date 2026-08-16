# OPHELIE Stage 3.1：A_glass Picard → Q → T（2026-08-11）

## 状态
**2026-08-11 用户普通终端：编译成功 + `passed=1`。**

## 变更
1. `runFrenchReducedEmOrSelfInductionForThermalHandoff`：
   - 默认：coil-only `runFrenchReducedEmPipeline`
   - `--self-induction`：`runFrenchReducedSelfInductionPicard`，再做 Q handoff
2. frozen-Q one-way 与 diffusion one-way 两条热路径共用上述 EM 入口
3. CLI 增加 `--self-induction-tol=`（全局）
4. Stage3 测试输出 Picard 指标：`self_induction`、`final_J_rel`、`A_ind_over_A_coil`、`em_phi_gate`

## 不包含
- 尚无 `σ(T)` 反馈
- 不改变 Stage2 Picard 求解器本身
- 不在 Cursor 沙箱内跑 GPU

## 用户编译 / 运行

```bash
cd ~/SPHinXsysSYCL/build
cmake --build . --target test_3d_ophelie_french_complex_joule_to_heat_one_way -j"$(nproc)"

RELAX_RELOAD=~/SPHinXsysSYCL/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_glass_relax/bin/reload

./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload-dir="$RELAX_RELOAD" \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-current-form=edge-flux --ophelie-edge-flux-complex=1 \
  --self-induction --self-induction-max-iter=25 --self-induction-relax=0.3 \
  --self-induction-tol=1e-3 --self-induction-phi-tol=2e-4 \
  --thermal-dt=0.1 --thermal-steps=1 \
  2>&1 | tee stage3_1_aglass_to_thermal.log
```

## 期望关注行
```text
self_induction=1
final_J_rel=...
j_ok=1
phi_eq_res_vol≈1.78e-4
phi_ok=1          # with phi_tol=2e-4
A_ind_over_A_coil≈0.16
P_joule_W≈2.24    # should move toward Picard raw power, not coil-only 2.30
E_joule_J ≈ E_thermal_J ≈ E_power_expected_J
em_ok=1 thermal_ok=1 passed=1
```

## 实测结果（2026-08-11）
编译：仅既有 warning（`CL/sycl.hpp` deprecated、`sscanf` float/double、C++20 extension），无 error，链接成功。

运行最终行：
```text
self_induction=1 self_induction_iters=17
final_J_rel=0.000698801 j_ok=1 picard_converged=1
A_ind_over_A_coil=0.160286 B_ind_over_B_coil=0.13538
P_joule_W=2.24307
phi_eq_res_vol=0.000178259 em_phi_gate=0.0002 phi_ok=1
E_joule_J=E_thermal_J=E_power_expected_J=0.224307
em_ok=1 thermal_ok=1 passed=1
```

对照：
- Picard-only Stage2：`P≈2.2428 W`，`A_ind/A_coil≈0.160`
- coil-only Stage3.0：`P≈2.30193 W`
- 本 run：`P≈2.24307 W`，与 Stage2 一致，说明热源确实来自 `A_coil+A_glass`

Picard 在 outer iter 17 因 `J_rel=6.99e-4 < 1e-3` 且 `phi_eq_res_vol≈1.78e-4 < 2e-4` 收敛退出。
