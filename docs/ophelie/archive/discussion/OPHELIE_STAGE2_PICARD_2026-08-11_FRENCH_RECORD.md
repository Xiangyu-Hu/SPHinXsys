# OPHELIE Stage 2 Picard 记录（2026-08-11）

## 1. 背景
本次仅针对 French glass 的 `A_glass` self-induction Picard（complex edge-flux 路径）做一次收敛性审计，重点检查：

- `J_rel` 是否收敛（`j_ok`）
- `phi_eq_res_vol` 是否过你设置的 `phi_tol` 收敛门槛（`phi_ok`）
- 若仅 `phi` 未过门槛，是否仍可作为“可放行”的工程前提继续往热链推进

## 2. 运行命令（原始）
在 `build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_self_induction_picard/bin/` 下执行：

```bash
./test_3d_ophelie_french_self_induction_picard \
  --reload-dir="$RELAX_RELOAD" \
  --dp=0.015 \
  --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-current-form=edge-flux --ophelie-edge-flux-complex=1 \
  --self-induction-max-iter=25 --self-induction-relax=0.3 \
  --self-induction-tol=1e-3 --self-induction-phi-tol=1e-4
```

## 3. 关键日志指标（最终）
```text
[ophelie] self_induction_complex | outer iter 25
  J_rel=0.000014
  phi_eq_res_vol=0.000178
  phi_imag_rel_res=0.000172
  phi_real_rel_res=0.000139
  picard_converged=0
...
test_3d_ophelie_french_self_induction_picard
  final_J_rel=1.44154e-05
  j_ok=1
  phi_eq_res_vol=0.000178261
  phi_ok=0
  picard_converged=0
  converged=0
  P_joule_W=2.2428
  A_ind_over_A_coil=0.160174
  B_ind_over_B_coil=0.135285
  passed=0
```

## 4. 结论与当前决策
1. **`J_rel` 已满足收敛门槛**：`j_ok=1`，`final_J_rel=1.44e-5 < self-induction-tol=1e-3`。
2. **`phi_eq_res_vol` 稳定卡在约 `1.78e-4` 的平台**：未过你设置的 `self-induction-phi-tol=1e-4`，因此 `passed=0`。
3. 这次没有出现发散、非有限值或明显的迭代失控迹象；`phi_eq_res_vol` 在外迭代推进时几乎不继续下降。
4. 按你的指示：本次**不再为该 phi gate 重跑调参**，交由后续与 GPT 讨论是否能放行。

## 5. 下一步
在确认“除 phi gate 外无其它红旗”的前提下，继续往下推进热链/热耦合相关的 Stage（一旦你要求再次收敛性讨论，我们再停下来让你做抉择）。

