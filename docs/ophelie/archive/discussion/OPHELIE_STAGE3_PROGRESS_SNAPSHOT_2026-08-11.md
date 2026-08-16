# OPHELIE Stage 3 进度快照（2026-08-11）

## 总览
French natural 几何（`R=0.25`, `H=0.185`, `dp=0.015`, `n=10474`）上，**电磁 → 焦耳 → σ(T) → 热扩散/冷壁** 可验证链已打通。当前无阻塞。

| 阶段 | 内容 | 结果 |
|------|------|------|
| 3.0 | coil-only → Q → T（能量闭合） | `passed=1` |
| 3.1 | + A_glass Picard | `passed=1`，`P≈2.24 W`，`A_ind/A≈0.16` |
| 3.2 / 3.2b | provisional Arrhenius σ(T) ± A_glass | `passed=1`，`σ~17.5`，`P≈2.49–2.57` |
| 3.3a | A_glass + 扩散/冷壁 | `passed=1`，`boundary=1`，`E_th/E_j≈0.40` |
| 3.3b | σ 外环后冻结 Q 扩散 | `passed=1` |
| **3.4a** | **thesis III.12 σ(T)** coil-only | **`passed=1`，`σ~46.5`，`P≈6.80 W`** |
| 3.4b | thesis III.12 + A_glass | **首跑 `passed=0`：φ=2.10e-4 > 门 2.0e-4；J/σ/热 OK。建议 `phi_tol=2.5e-4` 工程复跑** |
| 3.4c | thesis + A_glass + 冻结 Q 扩散 | 可选 |
| Stage 4 | 自然对流 / 50 kW 外层 | 未开始 |

## 关键更正
- 期刊 **Fig.2 ≠ σ(T)**（耦合流程图）。
- 真 σ(T)：论文 **III.12** `log10 σ = 3.7921 − 3179.8/T`（`σ(1473K)≈43`）。
- Table 1 的 `16@1473K` 是代表点，**不是** III.12 在该温度的值。

## 记录文档
- `archive/discussion/OPHELIE_STAGE3_*_2026-08-11.md`
- `reference/jacoutot_sigma_t/`
- 测试 README：`test_3d_ophelie_french_complex_joule_to_heat_one_way`

## 下一步（稳妥）
1. **现在**：3.4b = `--sigma-law=thesis-iii12 --self-induction ...`
2. 通过后可选 3.4c（扩散）或进入 Stage 4 准备：`target-power=50 kW` + thesis σ + Picard，再谈 Boussinesq / 流动。
