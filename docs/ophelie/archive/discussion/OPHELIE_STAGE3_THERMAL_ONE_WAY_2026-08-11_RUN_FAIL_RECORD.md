# OPHELIE Stage 3（热链 one-way）运行记录（2026-08-11）

## 1. 目的
验证 French 端到端热链：
`complex edge-flux EM → JouleHeat → one-way thermal（no feedback）`

目标测试：
`test_3d_ophelie_french_complex_joule_to_heat_one_way`

## 2. Cursor 沙箱误报（已关闭）
在 Cursor 受限 Shell 中运行时出现：

```text
CUDA_ERROR_OPERATING_SYSTEM = 304
sycl::exception → Segmentation fault
```

同一命令在用户普通终端直接运行后 **通过**。结论：CUDA 304 来自 Cursor 沙箱无法正常访问 GPU，不是 Stage3 代码或硬件故障。

后续 GPU 测试请在普通终端运行，不要依赖 Cursor Shell。

## 3. coil-only one-way 成功结果（重建二进制后）
几何：`R=0.25`，`H=0.185`，`dp=0.015`，`n=10474`，coil `R=0.285`，`z∈[-0.0225,0.2075]`，`f=282 kHz`，`σ=16`。

```text
P_graph_over_recon=0.877279
P_undirected_over_graph=1
P_legacy_unweighted_over_recon=259935   # diagnostic only
weak_current_antisym_rel_l2≈5.31e-8
phi_eq_res_vol≈1.770e-4
P_joule_W=2.30193
E_joule_J=E_thermal_J=E_power_expected_J=0.230193
em_ok=1 thermal_ok=1 passed=1
```

注意：旧二进制（2026-06-10）曾打印 `P_graph_over_recon≈2.6e5`；重建后才进入 Stage1 功率闭合口径。

## 4. 当前语义
本 case 默认仍是：

```text
coil-only complex edge-flux → Q → T
```

未默认打开 `A_glass` Picard。下一步（Stage 3.1）通过 CLI `--self-induction` 接入：

```text
A_coil + A_glass Picard → Q → T
```

## 5. 下一步
1. 重建 `test_3d_ophelie_french_complex_joule_to_heat_one_way`
2. 用 `--self-induction` + engineering `phi_tol=2e-4` 跑 A_glass → Q → T
3. 通过后再进入 `σ(T)` 双向反馈
