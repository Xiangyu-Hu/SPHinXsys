# OPHELIE Stage 3.2：σ(T) ↔ EM ↔ Q ↔ T 最小闭环（2026-08-11）

## 状态
**2026-08-11 用户普通终端：`passed=1`（coil-only + σ(T) 外环）。**

## 做了什么
1. `makeProvisionalJacoutotLikeSigmaTemperatureLaw()`：Arrhenius 近似，锚定 `σ(1473K)=16`，`Ea/R=15000 K`
   - **不是**论文 Fig.2 数字化；`paper_digitized_=false`
2. `UpdateOphelieGlassSigmaFromTemperatureCK` + radial T seed（machinery）
3. `runFrenchReducedSigmaTEmThermalCoupling`：外层
   `T → σ(T) → EM(+optional A_glass) → Q → T`
4. Stage3 测试 CLI：`--sigma-t`、`--sigma-t-max-iter=`、`--sigma-t-relax=`、`--sigma-t-tol=`、`--sigma-t-seed-delta-k=`

## 实测结果（无 self-induction）
```text
sigma_t_iters=6 sigma_t_converged=1
sigma_mean=17.489  sigma_min=16.019  sigma_max=18.355
sigma_rel=0.0394   power_rel=0.0452
P_joule: 3.915 → 3.431 → 3.092 → 2.855 → 2.689 → 2.573 W
P_graph/P_recon≈0.87–0.876  P_undirected/P_graph≈1
E_joule=E_thermal=E_power=0.257306 J
em_ok=1 thermal_ok=1 sigma_spatial_ok=1 sigma_coupling_ok=1 passed=1
paper_digitized_sigma=0
```

解读：
- 径向 T seed（+100 K）使首步 `σ` 空间非均匀，外环欠松弛后向较低、更均匀的 `σ` 收敛
- `P` 随 `σ` 外环单调下降并稳定，说明 `T→σ→Q` 反馈链路工作正常
- 尚未接 `A_glass` Picard；也尚未替换为论文 Fig.2 数字化 `σ(T)`

## 用户编译 / 运行

```bash
cd ~/SPHinXsysSYCL/build
cmake --build . --target test_3d_ophelie_french_complex_joule_to_heat_one_way -j"$(nproc)"

RELAX_RELOAD=~/SPHinXsysSYCL/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_glass_relax/bin/reload

# 建议先不加 self-induction，验证 σ(T) 外环本身（更快）
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload-dir="$RELAX_RELOAD" \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-current-form=edge-flux --ophelie-edge-flux-complex=1 \
  --sigma-t --sigma-t-max-iter=8 --sigma-t-relax=0.3 --sigma-t-tol=0.05 \
  --sigma-t-seed-delta-k=100 \
  --thermal-dt=0.1 --thermal-steps=1 \
  2>&1 | tee stage3_2_sigma_t_coupling.log
```

通过后可选加 `--self-induction ... --self-induction-phi-tol=2e-4`。

## Stage 3.2b：σ(T) + A_glass Picard（待跑）

```bash
cd ~/SPHinXsysSYCL/build
# 若未改代码可跳过 rebuild；改过 header 则重建
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
  --sigma-t --sigma-t-max-iter=8 --sigma-t-relax=0.3 --sigma-t-tol=0.05 \
  --sigma-t-seed-delta-k=100 \
  --thermal-dt=0.1 --thermal-steps=1 \
  2>&1 | tee stage3_2b_sigma_t_aglass.log
```

期望：
```text
self_induction=1  j_ok=1  phi_ok=1  A_ind_over_A_coil≈0.1–0.2
sigma_t=1  sigma_t_converged=1  sigma_max>sigma_min*1.05
E_joule≈E_thermal  passed=1
```

注意：每一步 σ 外环内会再跑完整 Picard，总时间大约是 Stage3.1×σ外环次数。

## Stage 3.2b 实测（2026-08-11，σ(T)+A_glass）
```text
self_induction=1 self_induction_iters=16
final_J_rel=9.58e-4 j_ok=1 picard_converged=1
A_ind_over_A_coil=0.177311 B_ind_over_B_coil=0.149803
phi_eq_res_vol=1.996e-4 em_phi_gate=2e-4 phi_ok=1
sigma_t_iters=6 sigma_t_converged=1
sigma_mean=17.489 sigma_min=16.019 sigma_max=18.355
sigma_rel=0.0394 power_rel=0.0423
P_joule_W=2.49376
E_joule=E_thermal=E_power=0.249376 J
em_ok=1 thermal_ok=1 sigma_spatial_ok=1 sigma_coupling_ok=1 passed=1
total_s≈4.83
```

对照：
| 路径 | P_joule [W] | A_ind/A_coil |
|------|-------------|--------------|
| coil-only + σ(T) | 2.573 | — |
| A_glass only (Stage3.1) | 2.243 | 0.160 |
| σ(T)+A_glass (本 run) | 2.494 | 0.177 |

自感略增强（相对 Stage3.1），功率介于 coil-only+σ(T) 与纯 A_glass 之间，符合双反馈叠加预期。

## 下一步（验证优先）
见 [`OPHELIE_STAGE3_3_THERMAL_DIFFUSION_COLD_WALL_2026-08-11.md`](OPHELIE_STAGE3_3_THERMAL_DIFFUSION_COLD_WALL_2026-08-11.md)：先 3.3a（A_glass+扩散/冷壁），再 3.3b（σ 收敛后冻结 Q 扩散）。**暂不做 Fig.2 数字化。**

## 期望关注
```text
sigma_t=1
sigma_max/sigma_min > 1.05
sigma_t_converged=1
E_joule≈E_thermal
em_ok=1 thermal_ok=1 sigma_spatial_ok=1 sigma_coupling_ok=1 passed=1
paper_digitized_sigma=0
```
