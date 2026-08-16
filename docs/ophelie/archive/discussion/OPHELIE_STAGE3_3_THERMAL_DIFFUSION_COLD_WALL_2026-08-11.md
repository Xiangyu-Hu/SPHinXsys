# OPHELIE Stage 3.3：热扩散 + 冷壁 Dirichlet（可验证优先）

## 状态
**2026-08-11 用户普通终端：Stage 3.3a / 3.3b 均 `passed=1`。Stage 3.3 关闭。**

## Fig.2 是什么（说明，本阶段不做）
Jacoutot et al., *Chemical Engineering and Processing* 47 (2008) 的 **Figure 2** 是玻璃电导率 **σ(T)** 实验/文献曲线。  
当前 `--sigma-t` 用的是锚定 `σ(1473K)=16` 的 **临时 Arrhenius**（`paper_digitized_sigma=0`），**不是** Fig.2 数字化。  
按稳妥顺序：先把 **传导 / 冷壁 / 能量不虚增** 验绿，再谈 Fig.2 替换。

## 本阶段做了什么
1. **3.3a**（已有管线，补跑法）：`--self-induction --thermal-diffusion`  
   `A_glass Picard → Q → 各向同性扩散 + 冷埚壳 Dirichlet`
2. **3.3b**（新）：允许 `--sigma-t --thermal-diffusion`  
   - σ 外环仍用集总加热（保持 σ/P 收敛可审计）  
   - 收敛后：**T 重置到 T₀**，对 **冻结最终 Q** 跑扩散+冷壁  
   - 验收用扩散门：`max_T>T0`、`boundary_compliance>0.9`、`E_thermal ≤ E_joule`
3. Helper：`applyOphelieFrozenQDiffusionFromUniformT0`（`electromagnetic_ophelie_thermal_diffusion_one_way.h`）

## 用户编译 / 运行（cwd = `build/`，普通终端）

```bash
cd ~/SPHinXsysSYCL/build
cmake --build . --target test_3d_ophelie_french_complex_joule_to_heat_one_way -j"$(nproc)"

RELAX_RELOAD=~/SPHinXsysSYCL/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_glass_relax/bin/reload

# --- 3.3a：先验证 A_glass + 扩散/冷壁（无 σ(T)，更快定位热侧问题）---
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload-dir="$RELAX_RELOAD" \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-edge-flux-complex=1 \
  --self-induction --self-induction-max-iter=25 --self-induction-relax=0.3 \
  --self-induction-tol=1e-3 --self-induction-phi-tol=2e-4 \
  --use-literature-thermal=1 --thermal-diffusion=1 \
  --thermal-dt=0.1 --thermal-steps=30 \
  2>&1 | tee stage3_3a_aglass_diffusion.log

# --- 3.3b：σ(T)+A_glass 收敛后，冻结 Q 做扩散/冷壁 ---
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload-dir="$RELAX_RELOAD" \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-edge-flux-complex=1 \
  --self-induction --self-induction-max-iter=25 --self-induction-relax=0.3 \
  --self-induction-tol=1e-3 --self-induction-phi-tol=2e-4 \
  --sigma-t --sigma-t-max-iter=8 --sigma-t-relax=0.3 --sigma-t-tol=0.05 \
  --sigma-t-seed-delta-k=100 \
  --use-literature-thermal=1 --thermal-diffusion=1 \
  --thermal-dt=0.1 --thermal-steps=30 \
  2>&1 | tee stage3_3b_sigma_t_frozen_q_diffusion.log
```

## 期望摘要行
```text
# 3.3a
self_induction=1 thermal_diffusion=1
boundary_compliance>0.9 max_T>T0 energy_cap_rel_err<=0
em_ok=1 thermal_ok=1 passed=1

# 3.3b
sigma_t=1 sigma_t_converged=1 thermal_diffusion=1
boundary_compliance>0.9 max_T>T0 energy_cap_rel_err<=0
em_ok=1 thermal_ok=1 sigma_spatial_ok=1 sigma_coupling_ok=1 passed=1
```

## 刻意不做（本阶段）
- 扩散进 σ 外环每一步（T 与能量账更难审计）
- Jacoutot Fig.2 数字化 σ(T)
- 自然对流 / 50 kW 全耦合

## 实测

### Stage 3.3a（2026-08-11，A_glass + 扩散/冷壁）
```text
self_induction=1 iters=17 final_J_rel=6.99e-4 j_ok=1 picard_converged=1
A_ind_over_A_coil=0.160286 B_ind_over_B_coil=0.13538
P_joule_W=2.24307 phi_eq_res_vol=1.78e-4 phi_ok=1
thermal_diffusion=1 thermal_steps=30 T0=1473 k=4
boundary_compliance=1
max_delta_T=1.02e-4 mean_delta_T=2.43e-5
E_joule_J=6.72922 E_thermal_J=2.71488 energy_cap_rel_err=-0.597
em_ok=1 thermal_ok=1 passed=1
```

解读：
- Picard / φ 门与 Stage3.1 一致量级（`P≈2.24 W`，`A_ind/A_coil≈0.16`）
- 冷壁 Dirichlet 全合规（`boundary_compliance=1`）
- `E_thermal < E_joule`：热量经冷壁离开，符合扩散+Dirichlet；点式焦耳闭合误差大是预期，不作为门
- 日志里 `max_T=1473` 是打印精度；`max_delta_T>0` 说明内部有升温

### Stage 3.3b（2026-08-11，σ(T)+A_glass → 冻结 Q 扩散/冷壁）
```text
sigma_t_iters=6 sigma_t_converged=1
sigma_mean=17.489 sigma_min=16.019 sigma_max=18.355
sigma_rel=0.0394 power_rel=0.0423
self_induction_iters=16 final_J_rel=9.58e-4 j_ok=1
A_ind_over_A_coil=0.177311 B_ind_over_B_coil=0.149803
P_joule_W=2.49376 phi_eq_res_vol=1.996e-4 phi_ok=1
thermal_diffusion=1 thermal_steps=30
boundary_compliance=1
max_delta_T=1.15e-4 mean_delta_T=2.65e-5
E_joule_J=7.48127 E_thermal_J=2.96508 energy_cap_rel_err=-0.604
em_ok=1 thermal_ok=1 sigma_spatial_ok=1 sigma_coupling_ok=1 passed=1
```

对照 3.2b / 3.3a：
| 路径 | P [W] | A_ind/A_coil | boundary | E_th/E_j |
|------|-------|--------------|----------|----------|
| 3.2b σ+A_glass（无扩散） | 2.494 | 0.177 | — | ≈1 |
| 3.3a A_glass+扩散 | 2.243 | 0.160 | 1 | ≈0.40 |
| 3.3b σ+A_glass→冻结Q扩散 | 2.494 | 0.177 | 1 | ≈0.40 |

σ/EM 与 3.2b 一致；扩散/冷壁与 3.3a 同类（约 60% 焦耳能经壁离开）。Stage 3.3 可关闭。
