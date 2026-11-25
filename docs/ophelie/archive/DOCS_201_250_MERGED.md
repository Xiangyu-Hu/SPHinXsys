# OPHELIE Docs 201–250 (merged archive)

Generated for GitHub hygiene. Original per-doc filenames listed in TOC.
Coverage: Docs 201–247 (incl. 239A; Doc241–245 plans; Doc246 calibration plan; Doc247 execution plan + Doc241–243 work summaries). Standalone `201_`–`247_` plan files removed after merge.

## Table of contents

- [201_OPHELIE_UNIFIED_COMPLEX_EM_AND_LIVE_BUOYANCY_CURSOR_PLAN_2026-10-30.md](#201-ophelie-unified-complex-em-and-live-buoyancy-cursor-plan-2026-10-30-md)
- [202_DOC201_ROUND_SUMMARY_AND_CLOSE_2026-10-31.md](#202-doc201-round-summary-and-close-2026-10-31-md)
- [202_OPHELIE_B2_LID_SIDE_SPLIT_ACOUSTIC_CURSOR_PLAN_2026-10-31.md](#202-ophelie-b2-lid-side-split-acoustic-cursor-plan-2026-10-31-md)
- [203_DOC201_DOC202_ROUND_PACK_2026-10-31.md](#203-doc201-doc202-round-pack-2026-10-31-md)
- [203_OPHELIE_WCSPH_2024_ISOLATED_IMPLEMENTATION_AND_VALIDATION_CURSOR_PLAN_2026-10-31.md](#203-ophelie-wcsph-2024-isolated-implementation-and-validation-cursor-plan-2026-10-31-md)
- [204_OPHELIE_VALIDATION_CORRECTION_FIG7_FIG8_CURSOR_PLAN_2026-11-01.md](#204-ophelie-validation-correction-fig7-fig8-cursor-plan-2026-11-01-md)
- [205_OPHELIE_MINIMAL_CYLINDER_CLOSED_CIRCULATION_CURSOR_PLAN_2026-11-01.md](#205-ophelie-minimal-cylinder-closed-circulation-cursor-plan-2026-11-01-md)
- [206_OPHELIE_BUOYANCY_CLOSURE_AND_THERMAL_HANDOFF_CURSOR_PLAN_2026-11-02.md](#206-ophelie-buoyancy-closure-and-thermal-handoff-cursor-plan-2026-11-02-md)
- [207_OPHELIE_ROBIN_VALIDATION_CALIBRATION_AND_BUOYANCY_CLOSE_CURSOR_PLAN_2026-11-02.md](#207-ophelie-robin-validation-calibration-and-buoyancy-close-cursor-plan-2026-11-02-md)
- [208_OPHELIE_THERMAL_FIELD_AND_BUOYANCY_CLOSE_CURSOR_PLAN_2026-11-03.md](#208-ophelie-thermal-field-and-buoyancy-close-cursor-plan-2026-11-03-md)
- [209_OPHELIE_CLEAN_RIEMANN_CAVITY_ROOT_CAUSE_AND_CLOSURE_CURSOR_PLAN_2026-11-03.md](#209-ophelie-clean-riemann-cavity-root-cause-and-closure-cursor-plan-2026-11-03-md)
- [210_OPHELIE_PARTICLE_MOTION_NATIVE_WORKFLOW_AND_CAVITY_TEST_PLAN_2026-11-04.md](#210-ophelie-particle-motion-native-workflow-and-cavity-test-plan-2026-11-04-md)
- [211_OPHELIE_CAVITY_CAUSAL_CLOSURE_AND_EM_VALIDATION_CURSOR_PLAN_2026-11-06.md](#211-ophelie-cavity-causal-closure-and-em-validation-cursor-plan-2026-11-06-md)
- [212_OPHELIE_TARGETED_HYDRO_EM_REPAIR_CURSOR_PLAN_2026-11-06.md](#212-ophelie-targeted-hydro-em-repair-cursor-plan-2026-11-06-md)
- [213_OPHELIE_LEGACYFLUX_EMF_SKIN_CLOSURE_CURSOR_PLAN_2026-11-07.md](#213-ophelie-legacyflux-emf-skin-closure-cursor-plan-2026-11-07-md)
- [214_OPHELIE_CONSISTENT_EM_REPAIR_AND_CONVECTION_VALIDATION_CURSOR_PLAN_2026-11-08.md](#214-ophelie-consistent-em-repair-and-convection-validation-cursor-plan-2026-11-08-md)
- [215_OPHELIE_CAVITY_NEIGHBOR_AND_EM_ADJOINT_REPAIR_CURSOR_PLAN_2026-11-08.md](#215-ophelie-cavity-neighbor-and-em-adjoint-repair-cursor-plan-2026-11-08-md)
- [216_OPHELIE_3D_EM_REFERENCE_AND_CAVITY_TIMESTEP_CURSOR_PLAN_2026-11-09.md](#216-ophelie-3d-em-reference-and-cavity-timestep-cursor-plan-2026-11-09-md)
- [217_OPHELIE_GENERAL_EM_REPAIR_AND_VALIDATION_CURSOR_PLAN_2026-11-09.md](#217-ophelie-general-em-repair-and-validation-cursor-plan-2026-11-09-md)
- [218_OPHELIE_D2_BOUNDARY_COMPATIBILITY_AND_FRENCH_THERMAL_CURSOR_PLAN_2026-11-10.md](#218-ophelie-d2-boundary-compatibility-and-french-thermal-cursor-plan-2026-11-10-md)
- [219_OPHELIE_VALIDATION_CLOSURE_AND_FRENCH_CONVECTION_RESTART_CURSOR_PLAN_2026-11-10.md](#219-ophelie-validation-closure-and-french-convection-restart-cursor-plan-2026-11-10-md)
- [220_OPHELIE_SORTING_THERMAL_IDENTITY_AND_HYDROSTATIC_REPAIR_CURSOR_PLAN_2026-11-12.md](#220-ophelie-sorting-thermal-identity-and-hydrostatic-repair-cursor-plan-2026-11-12-md)
- [221_OPHELIE_FORCEPRIOR_LIFECYCLE_REPAIR_AND_BUOYANCY_RECOVERY_CURSOR_PLAN_2026-11-12.md](#221-ophelie-forceprior-lifecycle-repair-and-buoyancy-recovery-cursor-plan-2026-11-12-md)
- [222_OPHELIE_FREE_SURFACE_PRESSURE_AND_FORCE_LIFECYCLE_CURSOR_PLAN_2026-11-15.md](#222-ophelie-free-surface-pressure-and-force-lifecycle-cursor-plan-2026-11-15-md)
- [223_OPHELIE_FORCE_RESET_REPAIR_AND_PAPER_TOP_BOUNDARY_CURSOR_PLAN_2026-11-17.md](#223-ophelie-force-reset-repair-and-paper-top-boundary-cursor-plan-2026-11-17-md)
- [224_OPHELIE_PAPER_TREND_THERMAL_ALIGNMENT_AND_ACCEPTANCE_CURSOR_PLAN_2026-11-17.md](#224-ophelie-paper-trend-thermal-alignment-and-acceptance-cursor-plan-2026-11-17-md)
- [225_OPHELIE_ABSOLUTE_TEMPERATURE_SLIP_AUDIT_AND_WARM_RESTART_CURSOR_PLAN_2026-11-18.md](#225-ophelie-absolute-temperature-slip-audit-and-warm-restart-cursor-plan-2026-11-18-md)
- [226_OPHELIE_TOP_SLIP_FIG8_AND_LOCAL_VISCOSITY_CURSOR_PLAN_2026-11-19.md](#226-ophelie-top-slip-fig8-and-local-viscosity-cursor-plan-2026-11-19-md)
- [227_OPHELIE_ACTUAL_TIMESTEP_PERIODIC_NEIGHBORS_AND_VISCOSITY_REPAIR_CURSOR_PLAN_2026-11-19.md](#227-ophelie-actual-timestep-periodic-neighbors-and-viscosity-repair-cursor-plan-2026-11-19-md)
- [228_OPHELIE_LINUX_ELASTIC_PADDLE_SOFTWARE_ADDENDUM_2026-11-20.md](#228-ophelie-linux-elastic-paddle-software-addendum-2026-11-20-md)
- [229_CLIENT_MATERIAL_AUDIT_AND_CURSOR_INPUT_RULES_2026-11-20.md](#229-client-material-audit-and-cursor-input-rules-2026-11-20-md)
- [230_CLIENT_CASE_PHYSICS_LANDING_2026-11-20.md](#230-client-case-physics-landing-2026-11-20-md)
- [231_OPHELIE_PHYSICS_NC_THERMAL_BALANCE_EM_AND_ELASTIC_PADDLE_CURSOR_PLAN_2026-11-21.md](#231-ophelie-physics-nc-thermal-balance-em-and-elastic-paddle-cursor-plan-2026-11-21-md)
- [232_OPHELIE_PRODUCTION_THERMAL_PROBES_BOUNDARY_AUDIT_AND_EM_REPAIR_PLAN_2026-11-22.md](#232-ophelie-production-thermal-probes-boundary-audit-and-em-repair-plan-2026-11-22-md)
- [233_OPHELIE_BUOYANCY_CAUSALITY_BOUNDARY_SCOPE_TIMESTEP_AND_EM_CONVERGENCE_CURSOR_PLAN_2026-11-22.md](#233-ophelie-buoyancy-causality-boundary-scope-timestep-and-em-convergence-cursor-plan-2026-11-22-md)
- [234_OPHELIE_FORCE_PRECISION_FIXED_T_RESPONSE_AND_EM_SPATIAL_REPAIR_CURSOR_PLAN_2026-11-22.md](#234-ophelie-force-precision-fixed-t-response-and-em-spatial-repair-cursor-plan-2026-11-22-md)
- [235_OPHELIE_OFF_FLOW_EM_GEOMETRY_AND_DIAGNOSTIC_CLOSURE_CURSOR_PLAN_2026-11-22.md](#235-ophelie-off-flow-em-geometry-and-diagnostic-closure-cursor-plan-2026-11-22-md)
- [236_OPHELIE_HYDROSTATIC_SUPPORT_AND_EM_GEOMETRIC_QUADRATURE_CURSOR_PLAN_2026-11-23.md](#236-ophelie-hydrostatic-support-and-em-geometric-quadrature-cursor-plan-2026-11-23-md)
- [237_OPHELIE_LOCAL_BACKGROUND_AND_EM_QUADRATURE_CLOSURE_CURSOR_PLAN_2026-11-23.md](#237-ophelie-local-background-and-em-quadrature-closure-cursor-plan-2026-11-23-md)
- [238_OPHELIE_OFF_DECAY_DISCRETE_WORK_AND_VALIDATION_CLOSEOUT_CURSOR_PLAN_2026-11-23.md](#238-ophelie-off-decay-discrete-work-and-validation-closeout-cursor-plan-2026-11-23-md)
- [239_FRENCH_THESES_EXTRACTION_AND_CURSOR_PLAN.md](#239-french-theses-extraction-and-cursor-plan-md)
- [239A_CLIENT_AND_FRENCH_CASE_SCOPE_ADDENDUM_2026-11-24.md](#239a-client-and-french-case-scope-addendum-2026-11-24-md)
- [240_OPHELIE_FLOW_THERMAL_CLIENT_HANDOFF_AND_LINUX_BRANCH_CURSOR_PLAN_2026-11-24.md](#240-ophelie-flow-thermal-client-handoff-and-linux-branch-cursor-plan-2026-11-24-md)
- [241_OPHELIE_FRENCH_TEMPERATURE_FLOW_AND_CLIENT_ELASTIC_RELEASE_CURSOR_PLAN_2026-11-24.md](#241-ophelie-french-temperature-flow-and-client-elastic-release-cursor-plan-2026-11-24.md)
- [242_OPHELIE_FINAL_PHYSICS_PREPARATION_AND_LINUX_ENTRY_CURSOR_PLAN_2026-11-24.md](#242-ophelie-final-physics-preparation-and-linux-entry-cursor-plan-2026-11-24.md)
- [243_OPHELIE_72H_PHYSICS_CLOSEOUT_AND_FRENCH_REPRODUCTION_CURSOR_PLAN_2026-11-24.md](#243-ophelie-72h-physics-closeout-and-french-reproduction-cursor-plan-2026-11-24.md)
- [244_OPHELIE_TARGETED_PHYSICS_DEVELOPMENT_AND_TEST_PLAN_2026-11-24.md](#244-ophelie-targeted-physics-development-and-test-plan-2026-11-24.md)
- [245_OPHELIE_THESIS_PROPERTIES_AND_PAPER_COMPARISON_PLAN.md](#245-ophelie-thesis-properties-and-paper-comparison-plan.md)
- [246_OPHELIE_ORDERED_CALIBRATION_AND_STIRRING_METRICS_PLAN.md](#246-ophelie-ordered-calibration-and-stirring-metrics-plan.md)
- [247_CURSOR_EXECUTION_PLAN.md](#247-cursor-execution-plan.md)
- [OPHELIE_Doc241_round_handoff_20261124_WORK_SUMMARY_ZH.md](#ophelie-doc241-round-handoff-20261124-work-summary-zh.md)
- [OPHELIE_Doc242_round_handoff_lean_20261124_WORK_SUMMARY_ZH.md](#ophelie-doc242-round-handoff-lean-20261124-work-summary-zh.md)
- [OPHELIE_Doc243_round_handoff_lean_20261124_WORK_SUMMARY_ZH.md](#ophelie-doc243-round-handoff-lean-20261124-work-summary-zh.md)

---

<!-- ===== BEGIN 201_OPHELIE_UNIFIED_COMPLEX_EM_AND_LIVE_BUOYANCY_CURSOR_PLAN_2026-10-30.md ===== -->

<a id="201-ophelie-unified-complex-em-and-live-buoyancy-cursor-plan-2026-10-30-md"></a>

# 201_OPHELIE_UNIFIED_COMPLEX_EM_AND_LIVE_BUOYANCY_CURSOR_PLAN_2026-10-30.md

# Doc201：统一复数电磁求解与真实推进下浮力/顶部验证

日期：2026-10-30  
**状态：已关闭（2026-10-31）** — 总结 [`202_DOC201_ROUND_SUMMARY_AND_CLOSE_2026-10-31.md`](202_DOC201_ROUND_SUMMARY_AND_CLOSE_2026-10-31.md)；证据 [`reference/doc201_runs/`](reference/doc201_runs/) · [`DELIVERABLE_CHECKLIST.md`](reference/doc201_runs/DELIVERABLE_CHECKLIST.md)  
执行对象：Cursor，OPHELIE分支  
依据：Doc200交付包及本轮讨论。该计划细化并纠正Doc200中验证范围偏大的结论。
## 1. 目标与本轮交付原则

主线A：使用同一生产候选的算子、边界和电流定义，直接求解完整复数A，不依赖已知梯度增量分拆，完成场精度、相位、守恒的联合验证。

主线B：使用法国实际压力/连续性/壁面推进路径，定位稳定分层残余；在位置不冻结的条件下验证顶部法向和切向响应。

优先复用已有case，新增最小诊断与必要局部修复。不要求本轮把所有项目标绿；明确失败应交付最小例及修复设计，不能靠切换测试专用算子获得整体PASS。

本轮不运行法国20/100s完整长耦合，不展开TEAM7 200Hz和频率/松弛/粒距矩阵。禁止修改独立A–phi Laplace分支。CPU隐式黏性保留，准确优先，不开展GPU隐式迁移或新增物理。

## 2. 先修正验证台账

| Doc200内容 | 正确范围 |
|---|---|
| T1/T3非零φ通过 | legacy-pairwise+legacy-flux的梯度场子问题通过，不等于生产DivSigmaA路径通过 |
| T2规范通过 | 已知gradχ分拆求解的诊断通过，不等于完整A输入求解通过 |
| B2误差约0.78% | 冻结布局的黏性滑移参考通过，不覆盖法国实际声学顶部 |
| 新压力初态 | 5s残余仅降约7%，不再将初值作为主因定论 |
| 粒子散度小 | 保留，但不能关闭边通量连续性失败 |
| 规定场功率3.49% | 已解释为粒子r²矩求积，非有限长自感误差 |

状态至少分：mechanism、field_accuracy、complex_phase、flux_conservation、boundary_accuracy、live_flow_accuracy、literature_alignment。

生产默认暂不静默切换。旧legacy、Doc198混合版及本轮统一候选分别命名、明确配置和hash。

## 3. P0：冻结与公共入口合同

结果根目录：`docs/ophelie/reference/doc201_runs/`。每次运行唯一目录，禁止覆盖；禁止新包递归包含旧PACKAGES。

登记commit、未提交diff、编译器/Real精度、二进制hash、输入/曲线hash、CPU/GPU路径。

建立一份`UNIFIED_EM_CONTRACT.json`，所有电磁测试共用：

- 复相量约定、E/J/Q公式；
- LHS/RHS、边界、规范、求解器、容差；
- 物理电流/边通量定义和重建方式；
- 实虚部相同配置的机器可读echo。

测试只允许改变A、sigma、几何和已声明物理边界，不按case标签自动切换算法。解析φ及χ仅供参考生成器使用，不能传给生产求解器作隐藏分解或答案注入。

## 4. A1：先完成相位旋转回归，隔离实虚部问题

对已有规定旋转A及同一σ/几何，分别输入A和iA。必须经同一完整复数入口求解φ/J/Q。

参考关系：J(iA)=iJ(A)、φ(iA)=iφ(A)（去掉规范常数）、Q不变。先关自感，后续完整自感阶段重复一次。

Doc200实部仍强制历史RHS、虚部使用DivSigmaA；这不算统一候选。先审计并让实虚部共享算子实现/边界规则，符号按相量方程体现，不用两条独立经验实现。

输出实虚φ、J、Q和两套残差；建议相位旋转前后J差≤1%、Q差≤2%。失败先修接线或算子不一致，不展开更多测试。

不得把已知相位解复制/交换实虚部充当第二次求解。

## 5. A2：选择一个相容离散候选

先交付`DISCRETE_FORMULATION.md`，从连续约束写到实际离散：

`E=-iωA-gradφ`，`J=σE`，`divJ=0`，绝缘`n·J=0`。

说明：源A如何进入电动势；φ梯度与散度如何配对；边界是否隐含于通量结构，是否还需显式贡献；常数φ零空间和可解性；J如何进入自感和Q。

可以修改OPHELIE内部φ算子，不等于改用户独立A–phi方法分支。不要因“不能动独立分支”而拒绝必要的OPHELIE局部修复。

主候选只选一套。可基于相容边通量修复，也可基于一致梯度/散度与重建设计；不预设必须保留某个旧模块。禁止针对旋转场用DivSigmaA、针对gradχ用legacy-flux的选择性路线。

若边量被定义为物理通量，φ残差必须约束它；若它只是辅助几何量，应更名、解释所有用途，并提供实际物理电流的守恒结构。不能只改报告标签关闭旧edge失败。

检查矩阵性质；PCG仅在适当内积下满足对称正定/半正定条件且零空间已处理时使用。实现尚未参数化或实虚部不同的模式标INVALID_CONFIGURATION，不用其崩溃否定数学方案。

### 几何事实修正

对同轴圆柱、恒σ，A=B0/2*(-y,x,0)在侧面及上下端面均切向，A·n=0，常数φ合理。不要沿用Doc200“该圆柱旋转场一般需要非零Neumann”的错误解释。

一般几何绝缘条件仍需相容处理；不能把齐次电势法向导数当普遍绝缘条件。显式边界项与隐含通量不能重复施加。

## 6. A3：同一候选的完整输入测试阶梯

所有测试调用同一生产候选入口，固定同一参数文件。参考解析生成器与求解器分离。

| 编号 | 输入 | 参考及目的 |
|---|---|---|
| U1 | 圆柱旋转A、恒σ、无自感 | 常数φ；验证不再制造虚假φ |
| U2 | 纯gradχ、恒σ、无自感 | φ=-iωχ+常数，J=Q=0；验证非零φ和绝缘 |
| U3 | 完整A0+gradχ | 直接输入总A，J/Q与U1相同；禁止分拆求解 |
| U4 | 平滑正σ(x)、纯gradχ | J仍为0；验证变σ一致性 |
| U5 | U3整体相位旋转90° | 检查非零φ下复数同构性 |
| U6 | 完整圆柱自感delta/R=1 | 与对应趋肤解析解比较复J、Q、功率 |

χ先线性，后低阶二次。U2/U3中求解器不能接收χ、解析φ或base/gradient标签；只能接收总A和物理配置。

为防测试特化，增加一次盲化输入检查：参考脚本先将总A写成普通CSV，求解器由通用导入接口读取。输出后由独立脚本比参考；可复用U3，不必增加很多工况。

先主档dp=0.02；出现问题只加一项有辨别力对照，候选通过后最多一档加密。不做无目标角度/粒距矩阵。

### 输出和门槛

phi去一个复常数后比较；J用复矢量体积加权L2；非弱场给相位；Q既给绝对/加权L2，也给峰和空间分区。

项目目标：U2零J相对规定非零尺度≤1%；一般复J L2≤5%；U3相对U1的J变化≤1%、Q≤2%；强场相位≤3°。门槛运行前冻结，不能删除边界区获得全域PASS。

自感检查最终J与最终A_total是否一致、固定点方程残差，不只看松弛增量；规定场使用prescribed参考，自感使用full_skin参考。

## 7. A4：守恒与功率不能选最好的一项

沿用正确无量纲化：epsilon_div=Lstar*RMS(divJ_complex)/Jstar。长度、J尺度和floor固定，零J问题采用外加场规定尺度。

同时交付：

1. 线性系统残差；
2. 物理边电流逐粒子净和（若采用边通量）；
3. 重建复J的独立散度；
4. 外表面和预定义内部控制体通量；
5. 边界法向J/通量与采样误差。

全域反对称抵消不能证明局部守恒。重建散度近零也不能自动抹掉同一物理通量的edge失败。建议无量纲连续性和边界相对通量目标≤1%，绝对值并列。

功率保留P_num、P_ref_discrete、P_ref_continuous、真实边能量表达及几何矩。P_J²/2σ与P_num如果代数恒等，只算公式一致，不算独立验证。

边功率与粒子功率的关系须从离散结构解释；不强行要求任意经验定义相等，但不能把未解释差异当物理守恒通过。

## 8. B1：真实推进下稳定分层残余定位

保留Doc200新初始化，不再预设初值是主因。首先重现mu=4、关闭EM/热演化、稳定分层的5s基线及均匀T零流控制。

不冻结位置，不重写rho/p，不清零速度，不去掉均值运动制造静止。

### 8.1 记录真实阶段，而不是事后混合快照

t=0推进前、最初几个声学步、0.01/0.1/0.7/2.3/5s按预算记录：

- 位置、u、rho、p、背景压、扰动压、T、质量和体积；
- 浮力、压力、壁面、黏性以及其它实际更新项的速度增量；
- 密度/压力更新前后变化及每步dt；
- 法向速度、每高度截面净通量、总动量变化。

CPU隐式黏性无显式力时，以执行前后动量增量除dt表示，明确算子分裂时序。累计实际增量核对Δu，禁止将不同时刻的“力”相加声称平衡。

分内部、顶部、侧壁、底部和角区统计；分类按实际位置及dp定义，保留样本数和原始CSV。

### 8.2 确定第一个异常环节

先回答：初始各粒子是否平衡？第一个明显不平衡发生在压力、密度更新、壁面还是黏性步骤？残余由局部传播还是整体均值驱动？

只允许基于证据增加一个消融：例如短时诊断冻结某项更新以定位，不将该诊断当生产修复。均匀T控制、机制M±和质量/EOS必须同步复核。

若初始局部误差小、推进后增长，优先查更新/壁面一致性；若初始边界已不平衡，查边界压力外推及支持域。

旧1%×M+目标保留。另报告残余相对同几何、同mu模式驱动流的比例；不同工况只作风险提示，不直接相减估计真流。

时间步减半历史无改善，不重复压dt；定位后最多一档加密确定空间趋势。

## 9. B2：真实声学顶部测试（本轮必须实现或给出具体阻碍）

复用Doc200已通过的切向黏性衰减，但将流体压力、密度、位置、真实French wall_closed及黏性omit一起推进。

若继续使用底部无滑移、顶部自由滑移参考：

`ux=U0*sin(pi*z/(2H))*exp[-nu*(pi/(2H))²*t]`，uz=0。

其它方向必须平移不变或有量化的端部隔离。周期SYCL未验证时不可静默使用旧host周期。必要时采用宽域中央观测并给端部污染证据，不把整个有限域当无限通道。

不允许冻结粒子位置、每步投影uz=0或直接重施参考速度。若法向约束本身属于正式边界算法，应如实记录其作用和动量交换，不把测试级清零混入。

记录切向剖面误差、衰减率、法向响应、壁面通量、压力/密度漂移、位移及能量变化。项目建议剖面L2≤5%、法向残余相对U0≤1%。outside=0仅辅助。

该参考覆盖切向剪切下的顶部响应，不自动覆盖所有浮力/自由面条件；仍需B1联合证据。

## 10. 法国与搅拌的恢复条件

本轮以修基础为主，暂不新增法国完整参数矩阵。

统一候选U1–U5、连续性通过后，可回归一组固定法国T/sigma的EM场，保留原始P及50kW归一后的Q；这属于条件项，不替代文献验证。

U6及浮力/顶部证据齐后，再提出法国短耦合配置与时长，不本轮自动长算。原Fig.8重建问题继续保留，不能假定温度形状已对齐。

搅拌可继续独立机械验证，但本任务不声称叶轮移动边界、力矩或热混合已通过。若固定mu/T且无Lorentz，EM不应直接改变纯机械速度；热耦合开启后影响需另核验。

TEAM7须使用统一候选和完整自感回归，不使用针对已知gradχ的测试分解。先50Hz，源场/实验总场口径分开。

## 11. 工程执行与停止规则

- 两线独立推进，不要求并行智能体，不同时占满设备跑大矩阵。
- 每项先写输入合同、阈值、粒子数、迭代/步数上限、内存和墙钟预算。
- 最多选择一个电磁主候选加一个有明确数学目的的对照，不反复换路线直到图好看。
- 不以更换名字、放宽阈值、只统计内部或换算法测试路径解除失败。
- 新候选只解决梯度不解决旋转场，或反之：标PARTIAL，继续修同一离散结构。
- 两次有依据的局部修改仍无改善：提交最小失败例、方程、排除项和下一设计，不扩展成无边界重写。
- 修改后只跑受影响基准和必要控制，不无条件回归全部历史测试。

## 12. 交付清单

1. `DOC201_PROGRESS.md`、`VALIDATION_LEDGER.csv`、合同JSON、环境/版本/hash。
2. `UNIFIED_COMPLEX_PATH_AUDIT.md`：实虚部入口及相位旋转结果。
3. `DISCRETE_FORMULATION.md`：算子/边界/电流/能量完整关系。
4. `UNSPLIT_INPUT_TEST_REVIEW.md`：U1–U6及盲化CSV，未运行项说明。
5. `CONSERVATION_AND_POWER_REVIEW.md`：所有指标分列，不互相替代。
6. `LIVE_STRATIFIED_RESIDUAL_LOCALIZATION.md`：真实阶段增量及最早异常。
7. `LIVE_ACOUSTIC_TOP_REVIEW.md`：位置真实推进的顶部验证。
8. 条件项`FRENCH_FIXED_T_UNIFIED_Q_REVIEW.md`。
9. `NEXT_STAGE_DECISION.md`：明确当前生产候选是否统一、哪些失败、下一唯一主优先。

附实际算子/边界/更新源码与diff、参考生成脚本、唯一边和粒子CSV、复数场、命令、日志、分析及绘图脚本。报告引用使用真实包内路径，打包前检查；禁止多层PACKAGES或同报告重复复制。

## 13. 最终必须回答

- 相同生产入口是否能处理A、iA及完整A0+gradχ？
- 是否仍把已知梯度单独求解，或按case切换RHS？若是，不能宣布统一通过。
- 场、局部连续性、边界和能量是否共同相容？
- 稳定分层残余最先从哪个更新项/区域产生？
- 顶部是否真正运行压力/密度/位置/声学壁面，而非冻结替代？
- 剩余误差对法国弱自然对流和热耦合搅拌的限制是什么？


<!-- ===== END 201_OPHELIE_UNIFIED_COMPLEX_EM_AND_LIVE_BUOYANCY_CURSOR_PLAN_2026-10-30.md ===== -->

---

<!-- ===== BEGIN 202_DOC201_ROUND_SUMMARY_AND_CLOSE_2026-10-31.md ===== -->

<a id="202-doc201-round-summary-and-close-2026-10-31-md"></a>

# 202_DOC201_ROUND_SUMMARY_AND_CLOSE_2026-10-31.md

# Doc201 轮次总结与关闭

日期：2026-10-31。  
计划：[`201_OPHELIE_UNIFIED_COMPLEX_EM_AND_LIVE_BUOYANCY_CURSOR_PLAN_2026-10-30.md`](201_OPHELIE_UNIFIED_COMPLEX_EM_AND_LIVE_BUOYANCY_CURSOR_PLAN_2026-10-30.md)  
证据树：[`reference/doc201_runs/`](reference/doc201_runs/)（单层；无嵌套 PACKAGES）  
核对：[`reference/doc201_runs/DELIVERABLE_CHECKLIST.md`](reference/doc201_runs/DELIVERABLE_CHECKLIST.md)  
Brief：[`reference/doc201_runs/GPT_DISCUSSION_BRIEF_2026-10-30_DOC201.md`](reference/doc201_runs/GPT_DISCUSSION_BRIEF_2026-10-30_DOC201.md)

## 一句话

冻结命名统一电磁候选 `doc201_unified_legacy_flux_norm_off`（A1–U6 PASS）；A4 分列 PARTIAL；B1 分层 hold 改善并冻结 thaw→ForcePrior on；B2 真活声学顶 FAIL 且 dens 局部修停码。**本轮关闭；生产 DivΣA 默认未静默回切。**

## §13 六问（关闭答）

| 问题 | 答 |
| --- | --- |
| 相同生产入口能否处理 A、iA、A0+∇χ？ | **命名统一候选能**；生产混合入口 **不能** 宣称统一（A1 INVALID） |
| 是否仍按 case 切换 RHS / 拆梯度？ | 统一候选：**否**。生产默认仍 imag DivΣA ≠ real LegacyFlux |
| 场、局部连续、边界、能量是否共同相容？ | **否（统一候选）**：场阶梯 PASS；U1 ε_div≈14% FAIL；P1 field_rel FAIL |
| 分层残余最先从哪来？ | **声学推进首漂**；ForcePrior 双重计数放大 hold 残量 |
| 顶部是否真跑压力/密度/位置/声学壁？ | **是**（B2 live）；Doc200 冻结顶不算本项 PASS |
| 对法国弱对流 / 热耦合限制？ | EM 相位可用但 U1 连续与功率场分需解释；活顶切向未收口 → 长耦合仍缓 |

## 门禁摘要

| 线 | 状态 |
| --- | --- |
| A1–U6 统一候选 | **PASS** |
| A4 | **PARTIAL**（ε_div↔边权衡 + DivΣA/dual 已入账） |
| B1 hold | **改善**（5 s U_rms ~1.2e-5→9.8e-7） |
| B1 thaw 政策 | **默认 ForcePrior on**（thaw 关 owns） |
| B2 活声学顶 | **FAIL + dens 停码**（壁声学定位） |
| F1 法国固定 T 统一 Q | **未跑**（条件项） |

## 产物入口

- 进度 / 台账 / 合同：`reference/doc201_runs/{DOC201_PROGRESS,VALIDATION_LEDGER,UNIFIED_EM_CONTRACT}.*`
- 命令：`reference/doc201_runs/COMMANDS.md`、`scripts/`
- 下一轮选项：`reference/doc201_runs/NEXT_STAGE_DECISION.md`（**Doc201 已关闭**）

## 明确不做（关闭后仍成立）

- 不以 DivΣA 对照冒充统一候选 A4 PASS  
- 不以统一候选静默替换生产默认  
- 不把 0.2 s thaw U 分叉当自由对流校准  
- 不再叠 B2 dens 局部修  


<!-- ===== END 202_DOC201_ROUND_SUMMARY_AND_CLOSE_2026-10-31.md ===== -->

---

<!-- ===== BEGIN 202_OPHELIE_B2_LID_SIDE_SPLIT_ACOUSTIC_CURSOR_PLAN_2026-10-31.md ===== -->

<a id="202-ophelie-b2-lid-side-split-acoustic-cursor-plan-2026-10-31-md"></a>

# 202_OPHELIE_B2_LID_SIDE_SPLIT_ACOUSTIC_CURSOR_PLAN_2026-10-31.md

# Doc202：B2 盖/侧分体声学（非 dens）

日期：2026-10-31  
前置：Doc201 关闭（B2 dens 局部修停码；`ablate=wall` 定位壁声学）  
证据树：`docs/ophelie/reference/doc202_runs/`（单层；禁止嵌套 PACKAGES）

## 1. 目标

在**不再叠 dens 开关**的前提下，用可测假说分离「盖」与「侧/底」声学接触对活顶剪切的责任。

## 2. 已否决（不得重做）

| 项 | 证据 |
| --- | --- |
| 盖 free-slip dens 鬼元 + 一步力∥n | L2≈1.05 |
| H-dens-fs | L2≈1.05 |
| 全壁 free-slip dens/力∥n | L2≈0.99 |
| `ablate=acoustic` 作壁对照 | disp=0（无效） |

## 3. 本轮唯一主假说 H1（反臂）

**侧壁+底：整对跳过声学壁接触；盖：保留现有 free-slip 法向不可穿透（二步法向反射 + DissipativePJump·n；一步力∥n）。**

可测：

| 结果 | 判读 |
| --- | --- |
| core L2 → ~0.12（近 `ablate=wall`）且 uz/U0≤1%、disp>0 | 杀剪切的是**侧/底 Riemann**，不是盖 |
| L2 仍 ~1 | 盖接触仍是主因 → 下一刀 H1b（盖整对 continue，侧保留） |
| uz/U0 崩 | 盖法向不可穿透不够 → 记录失败，不叠 dens |

## 4. 明确不做

- dens / H-dens-fs / AcousticDensMode 再扫  
- 静默改生产 DivΣA / 法国长耦合  
- 冻结位置 / 每步 uz=0 投影  
- 无假说连改第二刀（两次无改善才停码）

## 5. 合同（沿用 Doc201 B2）

`french_natural_reduced_pressure --doc201-b2-live-top-shear`，`implicit_split`，`end=0.1 s`，`U0=0.1`，`H=0.1`，`R=0.2`，`dp=0.01`，`c0=20`，`μ=10`。  
门禁：core `L2_rel≤5%`，`amp_rel≤5%`，`uz_rms/U0≤1%`，`disp_rms>0`。  
新 CLI：`--doc202-b2-ablate-side-acoustic`。

## 6. 交付

1. `DOC202_PROGRESS.md`、`VALIDATION_LEDGER.csv`、`METRIC_CONTRACT.md`  
2. `H1_SIDE_ABLATE_REVIEW.md` + 运行目录  
3. `NEXT_STAGE_DECISION.md`  
4. 源码：壁接触 `AcousticWallContactOmit` + 法国标记/CLI


<!-- ===== END 202_OPHELIE_B2_LID_SIDE_SPLIT_ACOUSTIC_CURSOR_PLAN_2026-10-31.md ===== -->

---

<!-- ===== BEGIN 203_DOC201_DOC202_ROUND_PACK_2026-10-31.md ===== -->

<a id="203-doc201-doc202-round-pack-2026-10-31-md"></a>

# 203_DOC201_DOC202_ROUND_PACK_2026-10-31.md

# Doc201+Doc202 轮次打包

日期：2026-10-31。

## 包位置

- 目录：`docs/ophelie/packs/Doc201_Doc202_UNIFIED_EM_LIVE_B2_2026-10-31/`
- tar：`docs/ophelie/packs/Doc201_Doc202_UNIFIED_EM_LIVE_B2_2026-10-31.tar.gz`（≈30 MB）
- zip：同名 `.zip`
- 校验：`docs/ophelie/packs/Doc201_Doc202_UNIFIED_EM_LIVE_B2_2026-10-31.SHA256`  
  - tar `824b669f4b0b098e256316d423493b3ab12d6b005fab2048763b19aeb10ac478`  
  - zip `02f270f02183d942521cb55b7f2abf867f5b14661229962153993753f0ff705f`
- 脚本：`docs/ophelie/packs/build_doc201_doc202_round_pack.sh`
## 含什么

计划与关闭总结；Doc201/Doc202 审查与台账；`french_natural` + `acoustic_ck` + EM 共享源码；GPU 脚本；`results/doc201_runs` + `results/doc202_runs` 全树证据（KEY_LINES / summary / CSV / 评审）。

## 先读（包内）

`README.md` → `docs/ROUND_DISCUSSION_BRIEF_2026-10-31.md` → Doc201 brief / Doc202 `H2_PAIRED_LID_FS_REVIEW.md`。


<!-- ===== END 203_DOC201_DOC202_ROUND_PACK_2026-10-31.md ===== -->

---

<!-- ===== BEGIN 203_OPHELIE_WCSPH_2024_ISOLATED_IMPLEMENTATION_AND_VALIDATION_CURSOR_PLAN_2026-10-31.md ===== -->

<a id="203-ophelie-wcsph-2024-isolated-implementation-and-validation-cursor-plan-2026-10-31-md"></a>

# 203_OPHELIE_WCSPH_2024_ISOLATED_IMPLEMENTATION_AND_VALIDATION_CURSOR_PLAN_2026-10-31.md

# Doc203：OPHELIE 自然对流的 WCSPH 2024 对照实现、验证与法国算例接入任务书

日期：2026-10-31  
执行对象：Cursor；基于当前 OPHELIE 分支及 Doc201/202 结果继续推进。  
本轮主目标：得到一条可验证、可归因、可回退的自然对流计算路径，再判断法国几何内能否产生可信的环流。GPU 优先，准确性优先于速度。

## 0．先读本节：本轮究竟要做什么

1. 保留当前 Riemann WCSPH 基线，新增 **OPHELIE 局部的中心 WCSPH 候选路径**，参考 Reece 等 2024 年论文。
2. 新建小型验证入口；所有入口共用同一套局部算子。法国大 case 保留原入口，通过显式选项接入候选，禁止复制出第二份数千行法国主程序。
3. 逐项审计压力、连续性、密度稳定化、粒子整理、变黏度和壁面。不要把“未通过验证”直接写成“已证明有 bug”。
4. 先做标准方腔的热—流验证，再做法国几何固定温度场，随后做固定非均匀热源的热—流计算。电磁不作为前两步的阻塞条件。
5. 本轮不保证得到法国双环流；必须交付可复现的结果和失败归因。出现任意涡旋不等于复现法国论文。
6. 不改变 SPHinXsys main、其他流体 case、A-phi/Laplace 分支，不默认切换全框架的压力求解方法。

### 0.1 对上一轮的必要纠偏

- Doc201/202 的 B2 初始化横向速度依赖高度，而在封闭圆柱侧壁通常有非零法向速度。它与不穿透边界不相容，因此不能将其原解析衰减误差直接归咎于声学/壁面算法。
- 不再继续围绕该 B2 叠加 omit、半步删除、盖单项消融来追求 PASS。保留历史结果，标记 `REFERENCE_BC_INCOMPATIBLE`。
- 旧 B2 的冻结几何黏性算子测试可以保留其局部用途，但不能升级成完整水动力验证。
- 不能因为关闭 ForcePrior 后静止残差变小就删除生产浮力；压力与体力相互平衡本来是正确解的一部分。
- transport velocity 与 shifting 都承担粒子分布整理功能。不得仅凭名称差异另造一套重复模块，也不得在同一步叠加两套位置修正。
- 黎曼方法的中心压力部分可能已经与论文一致；必须先展开实际公式，再决定新增哪些计算。

## 1．来源、证据与本轮边界

### 1.1 主参考

Reece, Rogers, Fourtakas, Lind (2024), *Buoyancy-driven circulation and multi-component mixing using SPH with a new adiabatic boundary condition*, IJHMT 233, 125904。

- DOI：https://doi.org/10.1016/j.ijheatmasstransfer.2024.125904
- 作者机构全文：https://orca.cardiff.ac.uk/id/eprint/171280/1/1-s2.0-S001793102400735X-main.pdf
- 实施前本地保存 PDF、SHA256；公式以原 PDF 第 2 节为准，不能只复制本文中的摘要。
- 原文基础公式索引：连续性 Eq.(7)、压力 Eq.(8)、导热 Eq.(9)、EOS Eq.(11)、黏性 Eq.(12)/(14)、浮力 Eq.(15)、壁面 Eq.(18)/(19)。
- 本文涉及的阈值、测试拆分、软件结构和执行预算均为本项目建议，不能写成原文要求。

论文事实的最小摘要：该工作采用 WCSPH、Boussinesq、Tait EOS 和 mDBC；方腔采用算术平均的 Monaghan–Gingold 黏性、粒子整理，不采用密度重初始化/扩散；另有圆柱流动应用。其结果对邻域支持范围敏感。不要推导出“所有自然对流必须 δ-SPH”或“mDBC 必然比当前边界准确”。

### 1.2 独立基准

de Vahl Davis (1983), *Natural convection of air in a square cavity: A bench mark numerical solution*。

- DOI：https://doi.org/10.1002/fld.1650030305
- 原始表格 Table V 用于检查方腔整体换热和速度极值；不要混用不同文章的无量纲速度。
- 可读取的原文镜像：https://elearning.univ-bejaia.dz/pluginfile.php/310888/mod_resource/content/0/Ref%20%5B4%5D.pdf

### 1.3 当前工程证据

以 `Doc201_Doc202_UNIFIED_EM_LIVE_B2_2026-10-31.tar.gz`、当前工作树和实际运行命令为起点。已有包不是完整仓库；Cursor 必须在服务器工作树再次确认调用关系，不得照搬旧包行号。

当前可见代码存在 Riemann 声学两半步、背景压力分离、可选密度正则化、粒子对调和黏度，以及显式/隐式黏性路径。某对象被构造不等于它在当前生产命令中被执行。

## 2．代码隔离与 case 组织：明确回答“新建还是原地修改”

### 2.1 推荐：两个新验证入口＋一个原法国入口

| 入口 | 是否新建 | 职责 |
|---|---|---|
| `test_ophelie_wcsph_operator_validation` | 是，小入口 | 按 `--case` 运行压力/连续性单步、声波、静止平衡、黏性和导热验证 |
| `test_2d_ophelie_natural_convection_validation` | 是，小入口 | 2024 方腔及高 Pr 补充测试；调用同一套局部算子 |
| 现有 `test_3d_ophelie_french_natural_reduced_pressure` | 保留并小幅接入 | 基线/候选切换、法国固定温度场、固定热源以及后续电磁接入 |

如果仓库已有功能等价入口，扩展并复用，不再按文档编号新建目标。二维算子必须与三维共用维度模板；禁止“二维验证一套，法国三维手写另一套”。

### 2.2 建议目录（先适配真实工程结构）

```text
tests/extra_source_and_tests/shared/electromagnetic_ophelie/hydrodynamics/
  ophelie_hydro_options.h
  ophelie_central_wcsph_ck.h/.hpp
  ophelie_mdbc_ck.h/.hpp
  ophelie_viscosity_policy.h
  ophelie_hydro_diagnostics.h
tests/extra_source_and_tests/.../test_ophelie_wcsph_operator_validation/
tests/extra_source_and_tests/.../test_2d_ophelie_natural_convection_validation/
validation/ophelie_wcsph_2024/
  configs/ references/ scripts/ reports/
```

以上是建议职责，不强制生成空文件。新代码处于 OPHELIE 命名空间，避免与框架公共类型同名。若现有共享目录不同，沿用现有目录并记录映射。

### 2.3 不允许的改动

- 修改全局 `AcousticStep*` 默认实现以便仅让 OPHELIE 通过。
- 修改全局 EOS 默认值、核函数默认比例、壁面默认规则。
- 为每个验证复制一份力学求解器或法国主程序。
- 使用进程全局环境变量悄悄切换算法。
- 修改或重新验证无关电磁路径，借机合并生产电磁候选。

确需公共层扩展时，先尝试利用现有模板策略；仍不可行，新增默认关闭的扩展点，提供旧配置前后对比。不得改变默认数值行为。

### 2.4 建议配置接口

```text
--hydro-scheme legacy-riemann | central-wcsph
--eos linear | tait
--wall-scheme legacy | mdbc
--density-control none | legacy-regularization | hydro-compatible-diffusion
--particle-regularization legacy-transport | paper-shifting | none
--viscous-operator legacy | paper-mg
--viscosity-mean harmonic | arithmetic
--viscosity-integrator explicit | existing-implicit
--thermal-mode frozen-spatial | evolving
--heat-source off | prescribed-field | electromagnetic
```

名称可以适配现有 CLI，语义必须保持。配置需要校验非法组合，并输出完整 `effective_config.json`。`legacy-riemann` 未指定任何新选项时必须保持原行为。

至少定义三种命名配置：

1. `legacy_frozen`：原基线，参数全部显式冻结。
2. `central_controlled`：仅更换压力—连续性路径，尽量保持其余项一致；是归因实验，不称为完整论文复现。
3. `paper_aligned`：按原文核对的 EOS、黏性、边界、粒子整理和密度设置；不同于原文的项目必须逐项标注。

`central_controlled` 若与旧壁面无法形成自洽组合，标记 `UNSUPPORTED`，不要强行拼接后把失败用于否定中心离散。

## 3．P0：先审计，再写代码

### 3.1 冻结当前环境

记录 Git HEAD、dirty diff、编译器、SYCL 后端、GPU、精度、构建参数、命令、输入哈希。不要替用户清理未提交文件。

为每项已有验证建立复用表：入口、实现版本、配置、参考解、几何兼容性、是否仍有效、需要补哪项。原有测试只有在本轮改变相关路径时才重跑必要部分。

### 3.2 交付 `FORMULA_AND_CALL_GRAPH_AUDIT.md`

每个变量从写入到读取列出：

| 字段 | 必须回答的问题 |
|---|---|
| `Density` | 连续性更新、求和重置、正则化、边界外推分别何时执行？ |
| `Pressure` | 是 EOS 压力、总压力，还是减去背景后的压力？ |
| `BackgroundPressure` | 代表哪部分可平衡体力？是否参与每一项梯度/耗散？ |
| `ForcePrior` | 浮力、黏性、其他外力分别写入几次，何时清除？ |
| `Displacement/Position` | 物理速度与整理位移分别贡献多少？ |
| `Temperature` | 粒子随流携带、重映射、固定空间场覆盖分别在哪发生？ |
| `ParticleViscosity` | 使用动态还是运动黏度？平均、截断、隐式更新顺序是什么？ |

每项给出“论文公式→现有代码公式→候选公式→差异→测试编号”。核对质量、体积、密度因子、力/加速度单位与符号，不能只比类名。

## 4．压力与连续性候选：实现要求

### 4.1 内部粒子对

参考中心离散的数学目标为：

\[
\dot\rho_i=\rho_i\sum_jV_j(\boldsymbol u_i-\boldsymbol u_j)\cdot\nabla_iW_{ij},
\qquad
\boldsymbol a_i^p=-\sum_jm_j\frac{p_i+p_j}{\rho_i\rho_j}\nabla_iW_{ij}.
\]

这些式子不含额外 δ 型密度扩散，不含黎曼跳跃耗散。核修正作为显式选项，首次中心公式核对时关闭。

先展开当前 `AverageP/AverageU`、`DissipativeUJump/DissipativePJump`。若中心项已相同，复用中心项代码并新增 OPHELIE 策略，不复制重复循环。只有完全核对两半步、内接触及壁面后，才能将零耗散策略视为中心候选。

### 4.2 压力语义必须唯一

本轮标准封闭方腔候选使用去除常量参考重力后的 Boussinesq 体力，初始背景压力为零，不叠加旧的均值吸收/扰动动量试验开关。温差形成的平衡压力由求解器生成；这并不意味着删去浮力。

法国接入时允许继续采用减压形式，但必须写清 `p_eos`、`p_background`、实际进入梯度的压力以及相应体力。不同方案比较必须代表相同物理方程。

### 4.3 时间积分

- 第一版可以保留已经核对的预测—校正/两半步骨架；必须通过单步 RHS 和时间步收敛验证。
- “空间公式参考论文”与“逐项复现论文求解器”分开标注；不同积分不能冒充完全一致。
- 每步输出或采样记录压力项、体力项、黏性项和数值耗散项，避免重复累加。
- 新声学路径仍须满足声学 CFL；隐式黏性只解除相应黏性限制，不能解除声学及热扩散限制。

### 4.4 EOS 与声速

优先在相同线性 EOS、相同声速下比较空间离散，随后单独验证 Tait EOS。Tait 的参数按相同参考声速匹配，不同时更改可压缩性与耗散。

禁止用 `c=12*初始最大速度` 导致静止初态声速为零。采用预估速度与正声速下限，记录依据和实际最大 Mach；若动态更新声速，需要与 EOS 参数更新一致，第一版建议固定声速。

目标 `Ma_max<=0.1`、密度相对变化通常不超过 1%；达标不等于低速压力一定准确。对通过的候选额外做一次声速敏感性，不盲目将声速提高几个数量级。

## 5．密度稳定化与粒子整理

### 5.1 密度

- 基线保留原有行为；论文方腔候选使用 `density-control=none`。
- `none` 指无额外重置/扩散，连续性仍然求解。
- 审计并记录密度改变的三个来源：连续性、黎曼耗散、额外正则化。候选不允许留下隐藏的重置调用。
- δ 型扩散只在候选有实测噪声且基础算例已定位后启用，不作为首版默认。
- 不得用均值强制校正掩盖质量/体积误差；若采取校正，要单独报告校正量。

密度改变量对压力的影响按 EOS 量化；线性 EOS 下为 `delta_p=c0^2*delta_rho`。输出 RMS/max 及近壁/核心分区，不能只输出最终密度均值。

### 5.2 transport velocity 与 shifting

先审计现有 transport velocity：位置更新频率、位移限幅、壁面法向投影、是否修改物理速度、是否有对应连续性/动量/温度输运修正。

保留可复用的实现；只有原文与现有方法存在实际不同，才增加 `paper-shifting`。一次只启用一个整理方法。

若粒子沿输运速度移动，但场变量按物理速度的材料导数更新，需要核对所选整理形式要求的修正；不能凭“都是修正粒子位置”跳过此项。固定空间温度场诊断则按新位置重新评价 T，明确这是外部规定，不计为热方程验证。

输出 `regularization_displacement/dp` 的 RMS/max、邻居数、近壁粒子体积分数、整理前后温度/动能改变。对静止常量场和线性温度场做短检验；不能凭粒子图更整齐判定准确。

## 6．黏性：分开检验平均方式、算子和时间积分

### 6.1 禁止一次替换所有黏性细节

当前局部粒子黏度路径的调和平均，与候选算术平均分别定义为：

\[
\mu_H=\frac{2\mu_i\mu_j}{\mu_i+\mu_j},\qquad
\mu_A=\frac{\mu_i+\mu_j}{2}.
\]

以 `mu_i=4, mu_j=400 Pa s` 为单元检查时，分别约为 7.92 与 202 Pa s。这个数值例子是项目分析，不是文献数据。

但原文使用的运动黏度与当前动态黏度不能不经密度换算直接交换。相同密度下再比较平均方式；密度不相同时按原式完整实现。

`paper-mg` 的速度投影、维数系数、核导数、距离正则项以原 PDF 公式及算例说明核对。不要仅将调和平均改为算术平均，就声称复现了原文黏性算子。特别注意常黏度、近不可压缩条件下的等效关系，不代表离散粒子上逐点相等。

### 6.2 显式与隐式

- 首先在常黏度小测试用受 CFL 约束的显式参考；不拿大时间步隐式结果作精确基准。
- 现有隐式路径可以复用，但必须证明其小时间步极限对应选定的半离散 RHS。
- 当前适配中存在遍历次数/内接触/壁面系数补偿。切换黏性算子后禁止沿用未经复核的经验因子。
- 同一时间段不能既通过 ForcePrior 施加显式黏性，又隐式施加一次。
- 输出隐式真实残差、迭代次数、失败状态及动能耗散；分裂算法没有矩阵残差时，改用与同一离散算子的步长收敛对照，不伪造线性求解残差。
- GPU 能复用则复用；短验证允许 CPU 显式/隐式参考。报告数据传输和运行时间，不能混淆准确性与性能。

## 7．壁面候选 mDBC

### 7.1 实施顺序

在独立 OPHELIE 模块新增 mDBC，不替换全局接触壁面。先平直方腔，再圆柱。实现物理壁面位置、法向、镜像点、邻域插值、梯度修正、插值退化及角点规则。

每种规则输出覆盖率与退化计数。角点不能随意平均互相冲突的 Dirichlet 温度；按参考解的边界约定并在误差中单独处理角点。

### 7.2 三类边界分开

1. 流动力学：密度/压力外推与速度条件。
2. 热学绝热：只在绝热验证里使用相应镜像规则。
3. 法国换热：仍是 Robin 通量和辐射，不替换成绝热或固定冷水温度。

标准方腔四壁无滑移。法国侧底无滑移、顶部不穿透且切向滑移。不能认为“顶面不计算黏性”就自动证明整个顶部边界实现正确，仍要检查法向速度、法向压力平衡及近壁切向梯度。

### 7.3 几何一致性

压力、连续性、黏性、温度必须采用同一个物理壁面位置。记录第一层流体与边界粒子距离，不能让不同项对应相差半个 dp 的有效几何。

## 8．测试总表与复用原则

| 编号 | 内容 | 优先级 | 入口 | 通过后支持什么结论 |
|---|---|---|---|---|
| O0 | 单步压力/密度 RHS、粒子对系数核对 | 必做 | operator | 公式翻译正确 |
| O1 | 周期微幅声波 | 必做，小规模 | operator | 压力—连续性耦合、符号和步进正确 |
| O2 | 无驱动静止场与垂直稳定分层平衡 | 必做 | operator | 不应流动时无明显伪流 |
| V1 | 常黏度正弦衰减 | 复用有效旧结果，变更后补跑 | operator | 黏性系数及时间积分正确 |
| V2 | 分层黏度 Couette | 必做 | operator | 强变黏度动量传递可信 |
| T1 | 方形导热 | 必做，可复用现有热核 | operator | 导热与温度边界 |
| T2 | 绝热常量/余弦温度场 | 必做 | operator | mDBC 热边界和能量守恒 |
| C1 | 2024 差温方腔 | 核心必做 | 2D convection | 热—流自然对流整体验证 |
| C2 | 高 Pr 桥接 | C1 后必做一次 | 2D convection | 从空气向高黏玻璃过渡 |
| F1 | 法国固定空间温度场 | 核心必做 | 原法国入口 | 给定温差的真实环流响应 |
| F2 | 法国固定非均匀热源 | F1 后 | 原法国入口 | 自洽温度—流动形成能力 |
| F3 | 恢复电磁反馈 | F2 后续 | 原法国入口 | 完整应用，不是本轮首要门槛 |

2024 圆柱多组分混合不在首轮必做范围。它不能直接替代法国单组分玻璃验证，FTLE 等混合指标也不是本轮需求。C1通过后，如确需 mDBC 曲面验证，可复用法国圆柱几何的短测试；只有决定复现原文圆柱算例时才另建配置，不新造重复求解器。

## 9．O0/O1/O2：低成本排除基础错误

### O0：RHS 对照

- 使用小型规则粒子和轻微扰动粒子；分别给定常压、线性压力、平滑速度、零速度。
- CPU 双精度直接求和参考与 GPU kernel 比较；参考循环独立书写，禁止直接调用被测函数。
- 区分内部完整支持区与近壁区。常压在离散无序/截断邻域不必逐点零梯度，不能假定其机器精度为零；应报告其误差随分辨率/修正变化。
- 比较 force/mass 与 acceleration，不允许混用。
- 正确条件下内部成对压力力总和检查相对误差：双精度目标 1e-10，单精度 1e-5，按非零力绝对和归一并加零场保护。非对称修正若不守恒必须显式报告，不能强行套门槛。
- 相同公式 CPU/GPU RHS 目标相对误差：双精度 1e-9，单精度 1e-4；大消去情形同时报告绝对误差。

### O1：周期微幅声波

- 真正一维周期或真实二维周期方向，禁止用封闭端壁替代周期。
- 建议 L=1、rho0=1000、c0=10，密度相对振幅 1e-4；速度按传播波关系配置，关闭浮力、黏性、热源、整理及额外密度控制。
- 运行两个声学周期；检查波速、振幅、相位、总质量与发散。
- 首轮相位/振幅误差目标各 5% 内，并展示一次空间加密/时间减半改善。中心法可能较少耗散也可能更噪，不能只以振幅较大判优。
- 如果 GPU 周期暂不支持，允许 CPU 复用框架原生周期完成参考；不得为此迁移整个周期边界项目。GPU 用 O0 和封闭平衡验证覆盖本轮；报告该项 GPU 未覆盖。

### O2：平衡场

- O2a：均匀 T=Tref、u=0、p'=0；关闭整理先测，随后开启选定整理。
- O2b：稳定垂直温度分层，u=0。使用同一浮力公式，求得满足 `grad(p)=rho*a_b` 的初始压力，再按 EOS 反求密度。对 Boussinesq 常参考密度解析初始化，同时报告离散残差；有限弱可压缩误差需要区分。
- 在矩形域执行，先不使用圆柱；不靠冻结 Pressure 掩盖连续性漂移。
- 观察 10 个声学穿越时间，另做较长黏性时间尺度短检验。
- 输出伪速度、动能、密度、总质量和法向泄漏；速度用后续 C1 的 Uref 归一，建议平衡伪速度不超过 Uref 的 1%，且末段无持续增长。
- 不设逐点绝对严格零的错误标准；若加密后不改善，定位边界/整理/密度更新。

## 10．V1/V2：黏性测试具体要求

### V1：正弦衰减

使用与边界相容的周期剪切，或平板通道流向周期、跨壁速度为正弦且壁面满足零速度。解析幅值为 `A(t)=A0*exp(-nu*k^2*t)`，必须说明 k 与边界的对应关系。

比较幅值拟合衰减率、速度 L2、黏性耗散以及步长敏感性。目标：中细分辨率衰减率误差 <=3%，加密/减小时间步后改善。关闭压力路径的冻结算子检查与完整流动检查分开报告。

### V2：分层黏度 Couette

建议无量纲 L=1、rho=1，下壁速度 0、上壁速度 1，流向周期；两层等厚黏度分别 1 和 100，关闭热与浮力，界面由材料标记携带。

解析剪应力：

\[
\tau=\frac{U}{H_1/\mu_1+H_2/\mu_2},\quad
u(y)=\begin{cases}
\tau y/\mu_1,&0\le y\le H_1,\\
\tau H_1/\mu_1+\tau(y-H_1)/\mu_2,&H_1<y\le H.
\end{cases}
\]

两种平均方式都运行同一配置；比较分段速度、界面连续性和两侧剪应力，不能靠哪条曲线视觉更平滑选方案。目标速度归一 L2<=5%，远离界面与壁面处剪应力偏差<=5%；全域、界面、近壁误差均单列。

这项是人为分层材料测试，不等于玻璃连续变黏度。若不连续界面测试揭示两种平均方式不同，再补一项光滑 mu(y)，参考 `u(y)=U*integral_0^y(1/mu)dy/integral_0^H(1/mu)dy`；不因论文采用算术平均就忽略界面应力结果。

## 11．T1/T2：热与热边界

T1 使用 2024 的方形导热类型，也可复用几何等价旧 case。建议采用无量纲等价参数 L=1、alpha=1、右壁1其余0、初值0；这是加速后的验证配置，不称为原文物理时间复现。

稳态参考为独立求解的 Laplace Fourier 级数：

\[
\theta(x,y)=\frac4\pi\sum_{n=1,3,5,\ldots}
\frac{\sinh(n\pi x/L)}{n\sinh(n\pi)}\sin(n\pi y/L).
\]

采用指数比值避免 sinh 溢出；级数截断做参考收敛。角点 Dirichlet 不连续，报告全域以及距角点2h外误差，不用裁剪后的误差替代全域。

T2a 为全绝热常量温度保持；T2b 为全绝热余弦场 `theta=theta0+A*cos(pi*x/L)`，解析幅值 `A*exp(-alpha*pi^2*t/L^2)`，检验瞬态导热和边界。

目标：T1 中细网格 L2/温差<=2%；T2 归一衰减误差<=3%、热能漂移<=0.5%。存在变 cp 时按焓积分，不用 `m*cp(T)*T` 冒充严格能量。

## 12．C1：必须做的论文自然对流验证

### 12.1 核对后的物理输入

二维方腔；左热右冷、上下绝热，四壁无滑移。原文参数集中记录如下，启动前通过公式重新计算 Ra/Pr，允许舍入误差，不允许物性抄错后硬写目标数字：

| 参数 | 数值 |
|---|---|
| L | 0.0386 m |
| Th / Tc / T初始 / Tref | 303 / 283 / 293 / 293 K |
| rho | 1.204 kg/m3 |
| nu | 1.506e-5 m2/s |
| cp / k | 1006 J/(kg K) / 0.02587 W/(m K) |
| beta / 重力大小 | 3.43e-3 1/K / 9.81 m/s2 |
| 目标 Pr / Ra | 约0.71 / 约1e5 |

原文 h/dp 设置存在不同实验系列；本项目主系列先采用 h/dp=2，必要时比较2.5。这里 h 指核定义的平滑长度，不是截断半径；首先确认与 SPHinXsys 参数语义一致。

### 12.2 分辨率与运行

- 粗级 L/dp≈40：只做稳定性和调试，不用于最终通过判定。
- 中级≈80、细级≈160：用于定量结论；实际粒子中心与墙的位置按同一几何约定生成。
- 不进行 `分辨率×所有开关` 全排列。先中级筛选命名配置，选出合理候选后只对基线/候选加密。
- 用 `alpha=k/(rho*cp)`，`t*=alpha*t/L^2`，`u*=u*L/alpha`，`theta*=(T-Tc)/(Th-Tc)`；比较时刻 t*=0.5。
- 使用时间步估计与实际 dt 统计，不能把不同配置运行到不同物理时间直接比较。
- 严格二维优先；三维薄层若用来模拟二维，厚度方向必须周期并有充分支持。封闭薄夹层不是二维等价算例。本轮不为薄层绕路。

### 12.3 独立数值参照

de Vahl Davis 原文 Table V、Pr=0.71、Ra=1e5 的参考数：

| 量 | 参考值 | 定义 |
|---|---|---|
| 平均 Nu | 4.519 | 注意与具体壁面/截面平均定义对应 |
| 水平中心线的竖直速度最大值 | 68.59 | uL/alpha 量纲体系，位置 x/L≈0.066 |
| 竖直中心线的水平速度最大值 | 34.73 | 同一无量纲体系，位置 y/L≈0.855 |

先核对轴向、冷热方向、符号和极值位置再比较。2024 图中的温度/速度曲线可数字化补充，但本任务书没有附带已提取曲线，不得编造 CSV。数字化数据必须附图号、坐标变换、提取方法、不确定性及原图哈希；缺失时状态为 `REFERENCE_MISSING`，不阻止使用上述独立表格验证。

### 12.4 验收

核心量：冷热壁热流、平均 Nu、两条中心线、流函数、动能、密度、末段稳定性。Nu 从真实壁面法向温度梯度与独立热流累计两种方式比较，若不能一致，需要解释离散定义。

本项目目标：细级平均 Nu 误差<=5%，中心线速度极值误差<=10%，极值位置偏差<=0.05L；若完整参考曲线可用，速度 L2<=10%、温度以温差归一 L2<=5%。

末段建议取不少于0.05无量纲时间，均值漂移<1%；上下绝热通量接近零，冷热功率不平衡<=2%。若 t*=0.5 未达到该标准，报告同一时刻误差并继续到可判断稳定，不把“没稳态”直接写为通过。

全局能量平衡与热壁功率不平衡不同：瞬态允许能量积累，应满足 dE/dt=Qin-Qout。方腔有一个主涡属于合理目标，不要求它产生法国双环流。

### 12.5 如果失败，最多两轮有依据的修改

依次检查 Ra/Pr/轴向、真实边界位置与 h 支持、压力/密度耗散、黏性系数、温度输运与整理。每轮一次主要变化。

如果全中心配置压力噪声发散，先排除 O0/O1/O2 错误，再考虑必要的稳定化；此时标记“稳定化后的中心 WCSPH 候选”，不冒称原文无扩散配置。两轮仍无定量改善，交付定位报告，停止无限开关搜索。

## 13．C2：高 Pr 桥接，不直接拿空气验证替代玻璃验证

采用与 C1 相同方腔及边界，构造 Pr≈1e3、Ra≈1e4 的常物性工况。先用无量纲参数生成器计算 nu、alpha、beta、温差及尺度，保证 beta*deltaT 足够小；禁止把物性单位混用。

高 Pr 没有本任务书提供的定量参考曲线，不能把它写为已验证。可使用现有成熟有限体积/有限元求解器得到独立参考，并做至少两级网格收敛；若暂时没有，C2 标记 `CONSISTENCY_ONLY`，只支持稳定性/收敛观察。

无参考时先做一次低成本桥接，不为新写完整欧拉求解器停滞。F1 固定温度场诊断可以继续，但最终法国定量结论须保留高 Pr 验证限制。

## 14．F1：法国几何固定空间温度场

### 14.1 几何与参数冻结

- 玻璃直径500 mm已由用户确认；高度185 mm作为博士论文补充几何输入，执行时核对当前 case 实际高度并记录差异，不能静默改 STL。
- 重力方向由当前几何和代码核对，不照搬旧记忆：已有自然对流代码可采用 z 竖直，而其他几何采用 y。统一写入 `coordinate_map.json`。
- 侧底无滑移、顶部不穿透滑移；关闭电磁力、Marangoni、热源与温度方程。
- 先常黏度区分水动力，再恢复期刊物性对应的 mu(T)，不先将黏度调小到容易出涡。
- 建议启动 dp=10 mm 诊断、再5 mm；这是本项目预算建议，不是分辨率已足够的结论。mm级冷边界层若未解析，必须报告。

### 14.2 温度输入

优先复用已有文献温度剖面/二维初猜及其来源。若只有几条剖面，二维重建不是唯一解，至少保存一种合理重建及一个不确定性扰动版本。

输入格式至少含 r、z、T、单位、坐标零点、温度标度、来源、插值方法。只有归一化图片时，禁止直接当作绝对 Kelvin；缺少温差幅值时设为明确校准参数并报告扫描，不称为文献原值。

每次粒子位置变化后按当前位置重评 T(r,z)，再更新浮力和指定物性。这是外部固定温度场，不是冻结粒子携带温度；二者会产生不同问题。输出施加前后温度覆盖量，但不将其作为热方程能量收支。

### 14.3 两个比较，不做全组合

同一初态、同一 T(r,z)、同一物性：基线 vs 通过 C1 的候选。然后只在选定候选上比较常黏度与变黏度，必要时单独比较黏度平均方式。

物理时间以稳定性、动量扩散和周转时间判定，不照搬此前2.3s，也不预设几秒出现双环流。固定温度问题可从零速启动，也可使用有记录的热启动；比较必须控制初态。

### 14.4 必须输出

- 原始三维速度、温度、压力、密度、mu 及粒子位置。
- Python 的子午面温度图、速度矢量、流函数图，标出数据覆盖和插值掩膜。
- 在多个方位子午面检查一致性；先评估非轴对称分量，再决定是否进行方位平均。
- 主回流方向、涡心位置、主/次涡强度、轴线与侧壁上下行速度；小于数值噪声的闭合曲线不算次涡。
- 近壁法向速度、总动能、总质量、密度波动、中心/壁面力平衡。

验收分层：`CIRCULATION_OBSERVED` 只说明有可信环流；`TOPOLOGY_COMPARED` 需方向、位置与文献图逐项对照；`FRENCH_VALIDATED` 还需物性、温度、热源、分辨率及参考不确定性充分覆盖。本轮不因图相似直接给最后标签。

## 15．F2：固定热源的自洽热—流计算

F1 有合理流动后再开启温度方程。输入固定非均匀 Q(r,z)，尽量使用已有可追溯电磁快照或文献重建；保留来源不确定性。

- 积分功率按玻璃体积进行 `sum(Q_i*V_i)` 检查；Q为W/m3，不是K/s。
- 法国期刊目标功率为50 kW，底部h=35、顶部对流h=20，侧壁h有100/200/300/400工况。首轮先选一个有清晰图对应关系的工况，未确认对应图前不得声明对齐某条流线。
- 温度相关mu/cp/k按已核对曲线更新；本阶段sigma不参与反馈，EM不重新求解。
- 顶面按对流+四次方辐射，侧底按等效换热，不能固定边缘流体粒子为冷却水温度。
- 环境温度、发射率和完整物性曲线有博士补充/工程假设时记录；尤其不能混用博士底h=300与期刊底h=35。
- 50kW不是线圈总输入功率的同义词，目标应是玻璃吸收热功率。
- 能量使用 `H(T)=integral(cp(T)dT)`；用累计Qin、Qout与DeltaE检查，初始瞬态允许输入输出不相等。

验收：温度稳定性、能量闭合<=2%、关键温度剖面及流型随加密的变化。热源归一化不能改变“非均匀分布是否正确”的证据等级；均匀热源只能作对照，不保证产生文献双涡。

期刊采用稳态耦合求解，约20次外迭代不是20秒；约1500K不是统一初温。本阶段不能以未在论文给出的加热时长作为验收。

## 16．F3：后续电磁反馈接回原则

只在F2可解释后接回 T→sigma(T)→EM→Q→T,u。保留已冻结的生产EM基线，本轮不自动采用Doc201统一候选。

更新间隔和sigma欠松弛是数值策略，需要温度变化、功率、热源形状及更新间隔敏感性测试，不能从稳态外迭代次数推导物理更新秒数。法国期刊自然对流先不加入洛伦兹力。

## 17．Python 后处理：防止画出不存在的涡

### 17.1 坐标和插值

明确局部重力轴、径向方向及速度投影。径向速度为 `(u_x*x+u_y*y)/r` 的表达只适用于对应轴映射；不能把任意笛卡尔横向分量当成径向速度。

r=0特殊处理；圆柱外、粒子稀疏区、壁外均掩膜。不用无约束散点插值跨越空洞，不把NaN填零后绘制闭合流线。原始矢量与插值矢量抽样对照。

### 17.2 流函数与连续性

二维不可压缩方腔可用平面流函数；轴对称法国子午面使用含r权重的流函数。采用约定：

\[
u_r=-\frac1r\frac{\partial\psi}{\partial z},\qquad
u_z=\frac1r\frac{\partial\psi}{\partial r}.
\]

重建psi后反算速度，报告重建误差；轴对称散度为 `1/r*d(r*u_r)/dr+du_z/dz`，不使用平面散度替代。对于三维非轴对称明显的流，轴对称psi只反映平均流，不能冒充完整三维流线。

原始粒子算子散度与插值网格散度分别报告。实测环流速度需显著高于O2伪速度背景，并在两级分辨率及至少两个后期时刻保持方向/位置。

### 17.3 后处理自检

复用已有解析涡测试；补充无旋平行流/静止场不应识别成闭合涡。检查图上的速度归一化、色标和流线种子，不用每帧单独归一化掩盖衰减。

ParaView为可选复核，不阻塞Python定量报告。不能因为暂停ParaView就省略环流验证。

## 18．执行顺序、预算和失败分支

### 阶段A：审计与最小算子验证

完成P0、O0、V1复用检查；先中等精度CPU参考，再GPU对照。构建问题单独记录，不立即重构整个框架。

### 阶段B：闭合的标准自然对流路径

实现中心WCSPH及必要局部壁面，完成O1/O2、T1/T2、C1。黏性先常量；粒子整理先沿用经过审计的路径。初次跑通后才切换原文黏性/边界形成paper_aligned配置。

### 阶段C：变黏度与法国固定T

完成V2和C2短检验，运行F1基线/候选。C1通过而F1失败，优先查温度输入、物性和曲面/滑移边界，不返回去无目的重写整个压力模块。

### 阶段D：固定Q热流与交付

运行F2，说明距离法国目标尚缺哪项。F3可列为下一轮，不以未完成电磁反馈否定本轮水动力进展。

### 控制工作量

- 一个主验证只先跑一个中级配置；明显失败不先上细网格长时。
- 单GPU首次筛选建议预算每配置30分钟，允许检查点续跑；达到预算只表示 `INCOMPLETE`，不能凭短烟雾测试PASS稳态指标。
- 任一问题连续两轮没有主要指标改善，冻结失败配置并给出可检验的原因与下一决策，不继续堆开关。
- 参考缺失、实现失败、运行未完、精度不足四种状态分开。
- 不允许通过改初始场强迫出现双涡、改绘图或放宽阈值来关闭问题。

## 19．验收状态、输出及复现

每个运行目录至少有：

```text
effective_config.json
provenance.json
command.txt
stdout.log
monitor.csv
summary.json
reference_metadata.json
figures/
```

`summary.json` 至少包含 scheme、EOS、壁面、核支持、邻居数、整理方式、密度控制、黏性公式/平均/积分、精度、实际物理时间、dt统计、收敛状态、误差、门槛来源、reference_id、GPU/CPU、终止原因。

状态统一为 `PASS / FAIL / INCOMPLETE / REFERENCE_MISSING / UNSUPPORTED / CONSISTENCY_ONLY`；不使用一个总PASS掩盖子项缺失。

需要提供：

1. `FORMULA_AND_CALL_GRAPH_AUDIT.md`：每项公式与真实调用。
2. `IMPLEMENTATION_MAP.md`：新增/复用代码，隔离边界及默认行为。
3. `VALIDATION_LEDGER.csv`：已有证据复用、运行结果和缺口。
4. `CASE_REFERENCE_TABLE.md`：原文参数、项目参数、借用参数、工程假设分列。
5. `HYDRO_COMPARISON_REPORT.md`：基线与候选的定量差别，包含失败配置。
6. `FRENCH_CIRCULATION_REPORT.md`：F1/F2温度、矢量、流函数、方向与文献对照。
7. `NEXT_STAGE_DECISION.md`：保留哪个候选、理由、还缺什么、是否可以接EM。
8. 一个串行一键脚本，支持分阶段运行/续跑和跳过有效缓存，不将所有case全部重跑。

缓存键至少包括源码/二进制哈希、有效配置、输入和参考哈希；日志缺失或二进制改变不得复用旧PASS。脚本失败返回非零，不能用 `|| true`吞掉测试失败。

### 隔离回归

相同硬件/编译设置下，原法国入口未指定新选项的短基线结果应与改动前一致；如受并行非确定性影响，用重复运行噪声包络判断。编译一个原有无关流体目标，确认默认模板/接口未破坏；不要为了局部改动重跑整个仓库所有测试。

## 20．最终决策规则

- 如果旧Riemann和中心候选均通过C1：不宣布Riemann错误，按F1/F2准确性、稳定性和成本选用。
- 如果中心候选通过、旧方案不通过：用单项对照定位耗散/边界/密度处理，保留有证据的局部修改。
- 如果两者均失败：先查共同热模块、几何、黏性和后处理，不声称更换方法已解决。
- 如果C1通过、变黏度失败：优先修正黏性传递，不靠降低黏度凑环流。
- 如果F1有环流、F2没有：优先查形成的温度场和散热，不先否定浮力。
- 如果F2有合理环流、接EM后改变：查Q空间分布、sigma曲线与更新策略，不以总功率闭合代替空间正确性。

## 21．Cursor 执行前检查清单

- [ ] 我知道本轮只在OPHELIE新增候选，不替换所有流体算法。
- [ ] 我先检查已有入口和算子，没有复制法国大case。
- [ ] 我没有继续把不满足圆柱边界的旧B2当精确参考。
- [ ] 我核对了transport velocity，没有重复叠加shifting。
- [ ] 我明确压力和密度每个更新来源，没有保留隐藏耗散/重置。
- [ ] 我区分了算术/调和平均、黏性算子及时间积分。
- [ ] 我知道mDBC绝热处理不能替代法国散热边界。
- [ ] 我使用实际物理时间、统一无量纲和有来源的参考。
- [ ] 我知道空气方腔成功不代表高Pr变黏玻璃已验证。
- [ ] 我会提交Python原始矢量、温度剖面及流函数重建误差。
- [ ] 我没有通过扩大阈值、只改图或冻结压力伪造通过。
- [ ] 我会报告已完成、未完成和不确定性，并给出可复现命令。

**本轮完成标志：得到一条经过标准自然对流检验、隔离于OPHELIE内部的候选，并对法国固定温度场给出可信流动结论。完整电磁反馈和法国双环流定量复现作为随后明确的验收层级，不能提前宣称完成。**


<!-- ===== END 203_OPHELIE_WCSPH_2024_ISOLATED_IMPLEMENTATION_AND_VALIDATION_CURSOR_PLAN_2026-10-31.md ===== -->

---

<!-- ===== BEGIN 204_OPHELIE_VALIDATION_CORRECTION_FIG7_FIG8_CURSOR_PLAN_2026-11-01.md ===== -->

<a id="204-ophelie-validation-correction-fig7-fig8-cursor-plan-2026-11-01-md"></a>

# 204_OPHELIE_VALIDATION_CORRECTION_FIG7_FIG8_CURSOR_PLAN_2026-11-01.md

# Doc204：修正验证依据、确认法国环流、接入文献数字化的执行任务书

日期：2026-11-01。接续Doc203，禁止重做整套项目。  
优先范围：OPHELIE局部自然对流；保留现有Riemann候选，准确性优先，GPU优先。  
配套数据：本目录 `reference/`、`qa/`、`scripts/`，先读 `README_FIRST.md`。

## 0．本轮结论和验收目标

当前已证明：中心压力/连续性基础项与现有框架可以对齐；方腔能够产生自然对流；法国固定工程温度场下存在明显的浮力响应。

当前没有证明：C1标准方腔定量全面通过、法国场形成可信闭合环流、法国双环流拓扑与文献一致、完整2024论文配方已经实现。

本轮交付目标按先后为：

1. 修正C1的Ra/Pr与参考对应关系，解释两种Nu的分歧，输出真正的定量验收。
2. 利用已有F1结果/检查点导出空间数据，确认闭合环流与顶部尖峰的关系。
3. 接入本包Fig.7/8数字化数据，先比较形态，再在有明确温标与补全假设的条件下构造F1输入。
4. 如果水动力空间结果可信，再开始F2固定非均匀热源；不等待电磁路线全部验证。

### 禁止的结论跳跃

- `SMOKE_PASS`不等于标准验证PASS。
- β>0比β=0速度大，不等于已形成闭合环流。
- RMS下降不等于物理环流减弱；必须分离相干流、周向流、近壁尖峰。
- 加密后速度下降不等于排除了分辨率影响。
- 自由落体速度是标尺，不是高黏玻璃必须达到的速度目标。
- 常数背景压力/减压处理与壁面是否一致，需要验证，不能只凭类名判断。
- 无耗散候选在当前组合中失败，不等于2024论文方法无效；当前尚非完整原文边界/黏性/整理配置。

## 1．先冻结，不重建项目

### 1.1 复用

保留Doc203新增两个小入口及原法国入口；不新增Doc204版复制case，不复制法国主程序。复用O0结果作为公式证据，仅在相关模板/公式改变时重跑。

C1、F1各保留一个旧二进制/配置/日志作为对照。不要清理用户dirty工作树。记录HEAD、diff、编译参数、精度、GPU、二进制哈希和命令。

### 1.2 旧状态修订，只追加历史

建立 `DOC203_STATUS_CORRECTION.md`：

| 旧项 | 修订标签 | 原因 |
|---|---|---|
| C1整体PASS | PARTIAL_VALIDATION | Ra不匹配、一个速度偏差超门槛、两种Nu未解释 |
| F1 CIRCULATION_OBSERVED | BUOYANCY_RESPONSE_OBSERVED | 代理量不能单独证明闭合拓扑 |
| dp12说明“不是分辨率问题” | NOT_CONVERGED | 缺匹配β0、末段统计和空间比较 |
| t80“弱饱和” | STATIONARITY_NOT_ESTABLISHED | 单个终点与监测仍有波动 |
| 文献T缺失阻塞一切流型检查 | 只阻塞同文献条件比较 | 自身流场可立即后处理 |

旧文件保留，并用新入口状态表明确哪份是当前结论，避免多个文档同时写FAIL/PASS而无时间顺序。

## 2．P0：本轮先补的实现审计

### 2.1 验证配置和法国配置是否真正一致

分别列出：EOS、声速、背景压力、PerturbationMomentum、DissipativeU/PJump、limiter、核支持、密度重置、核修正、壁面、整理、黏性算子和积分、温度更新顺序。

本轮日志显示C1默认 `dens_sum=1, dens_reg=1, dens_norm=1`；源码会按初始密度求和作归一化。F1则执行 `density_summation + Internal regularization`。这两条路径不是同一个密度控制，不得写成“C1通过，所以F1完全相同配置已验证”。

C1源码实际每个外步更新 `BackgroundPressure=rho0*g*(L-y)`，与早期审计文档所写C1背景为零不一致。核对各算子如何使用总压/扰动压，更新公式审计；先不因文档错误直接删除背景压。

### 2.2 显式黏性步长，常物性也必须覆盖

F1代码中黏性时间步估计位于 `if(variable_properties)` 内。常物性分支的旧估计可能由其他时间步算子保证，也可能遗漏；必须追踪实际调用，不能猜。

无论常/变物性，都应输出：实际外步、声学步、黏性上限、黏性上限对应mu/h/系数，以及生效限制。若当前上限已正确，由审计证明即可，不重复开发；若遗漏，将常物性也纳入同一限制逻辑。短对照dt减半验证速度/尖峰是否改变。

不能把GPU隐式作为修复这一检查的前置条件。F1现有mu=4常物性规模可先显式，验证准确性后再复用隐式。

### 2.3 壁面法向挂起

dp10mm曾卡在CPU `NormalFromBodyShapeCK`，不要重复盲跑。加入阶段耗时与超时输出，定位具体shape/粒子/法向计算。

对理想圆柱允许在OPHELIE局部使用解析法向作为诊断候选：侧壁径向、顶/底竖直，角部按明确几何规则处理；检查方向和单位长度。这是局部替代方案，不是迁移整个level-set项目。与原法向在可运行网格比较后再决定采用。

## 3．C1：先修正基准参数

### 3.1 已发现的实际问题

现有输入 L=0.0386、rho=1.204、cp=1006、k=0.02587、nu=1.506e-5、beta=3.43e-3、g=9.81、deltaT=20，计算得到：

```text
alpha = 2.135855960582023e-5 m2/s
Pr    = 0.7051037278701198
Ra    = 120325.61889860145
```

这不是严格Ra=1e5基准，禁止继续引用Nu=4.519等Ra=1e5数字计算正式误差后标PASS。

### 3.2 两个命名配置

1. `paper_table_input`：保留已使用表格值，日志记录实际Ra/Pr，只报告该工况；不删除已有结果。
2. `dvd_exact_ra1e5_pr071`：保持L/rho/cp/beta/g/温差，按目标Ra、Pr同时生成nu和alpha，再令k=alpha*rho*cp。

生成式：

\[
B=g\beta\Delta T L^3,\quad
\nu=\sqrt{Pr\,B/Ra},\quad \alpha=\nu/Pr,\quad k=\alpha\rho c_p.
\]

本配置应约为：

```text
nu    = 1.657702883504525e-5 m2/s
alpha = 2.334792793668345e-5 m2/s
k     = 0.028279570667181474 W/(m K)
```

由程序重新计算Ra/Pr，双精度相对误差目标1e-8、单精度1e-5；不将文本常数覆盖实测值。若现有CLI不支持这些输入，只在验证入口新增目标Ra/Pr或nu/k参数，不修改公共物性默认。

### 3.3 最少运行集

先n≈80、后n≈160，相同命名配置与相同无量纲比较时间。n≈40只做短检查；不要再把全部scheme×网格排列重跑。保留当前Riemann+limiter，纯中心路线本轮不继续堆开关。

比较时间按新alpha重新换算，不能沿用旧34.88s作为所有配置时间。保存t*=0.5快照；另保存不少于0.05无量纲时间的末段统计，以判定稳态。

### 3.4 正式指标

Pr=0.71、Ra=1e5对应参考：平均Nu≈4.519；水平中线竖直速度最大值68.59，x/L≈0.066；竖直中线水平速度最大值34.73，y/L≈0.855。使用u*=uL/alpha。

中细级目标：Nu误差<=5%，两个速度极值各<=10%，位置偏差<=0.05L；均值漂移<1%，相应波动单列。所有量参与定量判定，不能只选接近参考的一项。

参考来源：de Vahl Davis (1983), DOI 10.1002/fld.1650030305，Table V。数据表和参考定义本地保存。原始参考区分壁面/截面平均Nu；指标必须对应。

## 4．C1：解释Nu约2倍分歧

本轮细级日志：Nu_grad_avg≈2.1431，Nu_contact_avg≈4.5034。先查定义，不预设哪条正确，更不能先把某个值乘2来过关。

### 4.1 独立的一维导热检查

复用小算子入口新增/启用已有平板导热子配置，不新建大case。L=1、k=1、左右温度1/0、上下绝热、流动关闭，稳态解析T=1-x、单位面积热流=1。

以与C1完全相同的热接触kernel和墙位置运行：

- 记录实际流体粒子与物理壁面距离、边界粒子距离；不能混用dp/2与dp。
- 展开接触k平均、当前wall conductivity factor、镜像和距离系数，确认单位。
- 对热核直接累积Q；独立边界拟合梯度计算q；检查热端输入=冷端输出，内部能量变化接近0。
- 原始首层差商可保留，但不能未经近壁一致性分析就视为准确参考。
- n=40/80两级足够定位几何因子。目标稳态温度归一L2<=2%、热流<=3%；这是项目门槛。

如果两种热流不一致，用一维解析结果定位到具体项。只修有证据的问题，再补跑C1受影响配置。

### 4.2 恢复C1换热统计

输出Q_hot、Q_cold、面积、Nu定义、温度梯度拟合窗口、瞬态能量积累。角点温度不连续引起的误差分开统计，不能任意剔除低值点。

报告接触Nu与独立梯度Nu的网格趋势；相互差别目标<=5%。若方法误差在粗级较大，以解释和加密改善为证据，不静默删除另一列。

程序新增 `smoke_pass`、`quantitative_pass`、`reference_match`、`stationary`、`energy_balance_pass`；进程返回值和报告状态按本次运行目的明确区分。`reference_match=false`时量化通过必须为null或false，不能给PASS。

## 5．F1：先导出已有空间数据，别先跑更久

### 5.1 选择现有结果

- dp15mm、beta>0，20s及80s附近。
- dp15mm、beta=0，已有20s对照。
- dp12mm、beta>0，20s附近；补相同设置beta=0短对照。

优先从服务器已有VTP或检查点导出；本包没有原始空间场不意味着服务器没有。只有缺失时才重跑最少必要快照；输出前设置专属目录，禁止自动backup覆盖其他长算例结果。

### 5.2 统一采样，修正终点歧义

旧summary与monitor末条不一定是同一个物理时间，不能直接比较。例如15mm案例的monitor末条与summary的U值已有差别。

每个snapshot带 physical_time、phase_time、step、configuration_hash；summary标记采用哪一帧或哪一个时间窗。对照应取相同的激活温度后时间，而不是文件最后一行。

### 5.3 必须导出的字段

```text
particle_id, x,y,z, ux,uy,uz, rho,mass,Vol,
T,mu,p_eos,p_background,p_momentum,
buoyancy_acceleration,viscous_acceleration,
distance_to_side,distance_to_bottom,distance_to_top
```

不要求首次导出所有三维力数组；如果成本高，先输出位置/速度/温度/密度/体积及壁距，足够后处理定位，再针对问题补力账本。

## 6．F1空间检查：把代理量升级成真实流型证据

### 6.1 统计分区

沿用坐标映射，现自然对流为z竖直；不得从其他STL默认y竖直。计算ur、uz、utheta。

输出全域、核心、近侧壁、近底、近顶、顶侧交线六类区域：体积/质量加权RMS、中位数、p95、p99、最大值、质量/体积、越界数和最大越界距离。

原先纯粒子算术平均继续保留以追溯，但不作为粒子体积不均匀后的唯一指标。不要直接删除高速点后称环流正确；先标记位置，再同时展示包含/排除该区的差异。

### 6.2 顶部与周向运动

已见15mm案例顶层uz RMS可大于全域U RMS，且有显著局部尖峰。目标不是让“某层所有粒子uz=0”，而是检查物理平面z=H处不穿透和近壁梯度。

- 区分顶层取样厚度、真实边界面插值速度与越界。
- 单列顶侧角部；比较法向方向、壁面距离、邻居覆盖与Pressure/continuity接触。
- 输出体积流量穿过顶部的正/负积分和净值；净值接近零不能掩盖正负相抵的局部泄漏。
- 轴对称驱动下检查utheta与子午流之比、不同方位切片差异；有限粒子噪声允许存在，但不能在明显非轴对称时只画平均图掩盖。

项目首轮筛查建议：时间平均顶部法向速度不超过相干子午面速度的5%；若只有近壁粒子量无法获得面速度，标为未完成，不假造通过。具体精度需结合分辨率与解析测试，不能靠裁剪调门槛。

### 6.3 绘图与流函数

每级至少画：原始子午面粒子矢量、经掩膜的网格矢量、T/Tmax或工程温度、流函数、utheta、近壁尖峰位置。图上有时间、dp、beta、mu、T来源。

方位平均之前先用至少4个方位检查；不做光滑到只剩两个涡的后处理。插值只能覆盖有邻居支持区域，空洞和壁外mask；不得填零制造闭合线。

轴对称约定：

\[
u_r=-\frac1r\partial_z\psi,\qquad
u_z=\frac1r\partial_r\psi.
\]

用物理r,z而非未换尺度的r/R,z/H计算导数。反算速度并报告重建误差和轴对称连续性残差；明显不满足连续性时，流函数只是投影结果，不能当证据。原始矢量与重建流必须方向一致。

### 6.4 判定标签

- `BUOYANCY_RESPONSE_OBSERVED`：已有beta对照支持。
- `CLOSED_CIRCULATION_CONFIRMED`：矢量、流函数、回流通道及连续性一致，且相干运动显著高于beta0背景。
- `TWO_CELL_TOPOLOGY_COMPARED`：主/次涡、方向、上升带与Fig.7约束逐项比较，注明参数和温度场不同。
- `FRENCH_VALIDATED`：本轮不能仅凭上述三项自动授予，还需文献工况、温度/热源与数值收敛等证据。

## 7．F1质量、密度控制与时间步

### 7.1 quiet初始化的账本

15mm结果mass约97.230→96.379kg，主要对应一次quiet-reset-mass，不能称运行持续漏质量，也不能用总体passed=1掩盖。

记录 `mass_before_reset`、`mass_after_reset`、`mass_end`；分别输出IC修改量、后续质量漂移。beta与beta0必须使用相同重置后位置、mass/rho/Vol，并保存初态哈希。

生产候选优先采用无质量重置的几何一致初态；若继续诊断重置，单独标签，禁止直接当正式守恒验证。一次改变IC后不要再逐步重写质量。

### 7.2 密度控制对照有限开展

先做F1短窗 `当前Internal reset` 与 `C1已验证的归一化处理` 的一致性审计。若要试另一路径，先静止/平衡+同一小网格，不直接拿80s长算做筛选。

密度重置前后记录delta_rho、通过EOS产生的delta_p和动能改变；不能只看rho均值。确定主误差与密度更新相关后才考虑delta-SPH，不再同时叠三种平滑。

### 7.3 dt敏感性

按§2确认常黏度CFL，候选同一小/中网格做dt减半。报告核心子午速度、p99、顶部法向通量及环流代理量；如果仅尖峰显著变化而相干流不变，不能再用全域RMS判物理幅值。

## 8．文献数字化：本包已经完成什么

### 8.1 可以立即使用的数据

Fig.7：20个可读箭头的近似方向、两个涡心搜索区域及上下行区。没有速度大小。搜索框不是精确涡心。

Fig.8：0.5–0.9可读等温线4928点；中央切线10个交点；部分线性重建。所有温度是 `T/Tmax`。

先运行 `scripts/validate_reference.py`，然后检查 `qa/fig8_contour_overlay.png`、`qa/fig7_direction_overlay.png`。Cursor不需要再次从PDF抠同一批曲线；修改标注必须保留原版本与差异说明。

### 8.2 立刻能做的对比

1. 对F1结果的主回流方向、上升带、次回流是否存在作Fig.7形态比较，说明当前工程温度不同，不能量化宣称复现。
2. 对未来F2自洽温度按其真实Tmax归一化，与Fig.8等温线位置和中央交点比较。
3. 对当前纯径向工程T作Fig.8叠图，清楚显示它缺少文献高度结构；这一步不需重新运行流动。

### 8.3 不允许直接加载部分CSV为绝对场

`fig8_partial_linear_reconstruction.csv` 有凸包外空白，热核高于0.9的具体温度未知，不是全域Kelvin输入文件。输入器必须拒绝将其误当 `r_m,z_m,T_K`。

元数据中 `absolute_Tmax_K=null`、`case_side_h=null` 是真实缺口，不允许自动填默认值。

## 9．从归一化约束构造可计算的F1温度：明确的不确定性分支

### 9.1 温标来源优先级

先检查期刊正文/其他图/已有博士提取材料能否提供同工况Tmax或独立绝对温度锚点；记录精确页码与工况。原文“约1500K”不能当精确最大温度。

如果有可靠Tmax：使用 `T_K=Tstar*Tmax_K`。如果只有一个绝对温度锚点且对应Tstar确定，可据此估计Tmax并传播误差。

如果仍缺：允许建立明确工程假设的温标家族用于诊断，文件名和provenance写 `RECONSTRUCTED_ASSUMED_INPUT`。给出中心假设及±10%敏感性只是一种项目探索范围，不是文献置信区间；不得通过扫温标只保留最像图的结果。

### 9.2 完整场重建

- 观测点作为等值约束，文字断口可以插值，但标为重建。
- 域边界处和0.9内部需要独立假设，绝不把线性插值平台0.9当真实最高温。
- 热核的最大值1约束、位置与边界温度仅可作为显式补全假设；不能凭插值算法自动知道。
- 使用平滑、无伪振荡的约束重建，至少两种合理补全进行敏感性。优先现有插值库，不能为此开发完整新求解器。
- 输出每个网格点的来源类型：OBSERVED_CONSTRAINT / INTERPOLATED / ASSUMED_CORE / ASSUMED_BOUNDARY。
- 将重建等温线叠回原图，报告归一化坐标误差；目标中位距离<=0.01、p95<=0.02，属于图形一致性目标，不是物理精度。

### 9.3 加载与力学一致性

明确R、H、原点与重力映射；转换成 `r_m,z_m,T_K` 时写出各参数。用固定空间温度，粒子移位后重新评价T。不要把按粒子初始位置固定的温度当空间固定场。

若选常黏度F1，均匀参考温度变化在理想Boussinesq方程中可由压力吸收；但有变物性/数值边界时不能随意demean绝对温度后又用错误mu(T)。记录浮力Tref和物性温标各自用途。

加载前检查Tmin/max、beta*deltaT、mu/cp/k适用范围。Fig.8温差可能远大于旧40K；不要不检查就沿用beta=3e-5，期刊beta主线为1e-5，博士参数单独标签。

### 9.4 最少F1新运行

在当前最可信网格先跑一个重建温度候选及对应beta0；随后只对一个补全/温标扰动做敏感性。先短窗稳定性，再按相干速度和黏性时间判断是否续跑。不能一开始开启热、电磁、变黏度和新壁面全部耦合。

如果已具备可信现有工程T环流，可先做上述输入对照，再逐项恢复mu(T)；不要求工程T必须长出文献双涡。

## 10．分辨率与运行时长：改正旧判据

- dp15→12使RMS下降，可能是噪声降低，也可能是离散耗散/几何改变。比较相干子午速度、涡心、psi强度、近壁尖峰，而非只比较RMS。
- 同两级在相同物理时刻和相同温度输入比较；12mm补beta0，不必立刻上5mm。
- 当两级相干流变化<=10%且关键方向/位置稳定，可作为初步收敛证据；否则报告未收敛，不能冻结默认15mm并宣布排除网格影响。
- 窄壁面层仍未解析时单列限制；不能以粒子总数一万就认为分辨率足够。
- `sqrt(g*beta*deltaT*H)`只作量纲标尺，不是速度应达到的目标。
- 以R或H计算的L²/nu不同，日志注明L的选择；旧43s主要对应R，而非所有方向的建立时间。
- 每个运行报告末段时间平均、标准差、趋势和实际完成时间，不用一个80s终点值宣称稳定饱和。

## 11．其余验证如何补，避免重走老路

| 项目 | 本轮动作 | 不需要做什么 |
|---|---|---|
| O0 | 复用有效证据；修改声学项才补跑 | 不再重写中心压力循环 |
| O2 | 对最终生产候选跑短平衡，检查与旧NoRiemann配置差别 | 不把旧静止PASS扩大到所有配置 |
| O1声波 | 已有结果就挂入账本；没有则列明确缺口，必要时小测试 | 不为周期边界重构另一项目 |
| V1常黏度 | 复用已有一致算子结果；确认F1 dt限制 | 不先解决所有GPU隐式性能 |
| V2变黏度 | F1恢复mu(T)前检查已有分层Couette结果；缺失才补 | 不把算术平均默认视为比调和正确 |
| T1/T2 | 本轮Nu问题用一维导热+绝热常量/余弦最小验证 | 不复制多个导热主程序 |
| C2高Pr | 保留CONSISTENCY_ONLY；如有成熟外部参考补定量 | 不新写大型欧拉求解器卡住主线 |
| mDBC | 只有定位到现有壁面误差才局部候选对照 | 不立即全框架替换所有边界 |

## 12．F2启动条件与范围

满足：C1严格基准基本通过、F1有空间可信环流、顶部异常已解释或受控、温度输入元数据完整。并不要求先完美匹配法国双涡才能探索F2，否则可能陷入固定T自身不完整的循环。

F2：相同法国几何、物性、散热边界，输入固定非均匀Q；温度与流动自行演化。不再每步覆盖T，不同时给温度演化指令。

Q优先已有可追溯电磁快照；必要时再数字化Fig.6。Fig.6不是本包已交付数字化数据，禁止把Fig.8温度当热源。玻璃吸收功率50kW与期刊保持一致；若采用其他功率，作为另工况。

热边界：期刊侧h按对应100/200/300/400工况、底35、顶对流20加辐射。环境温度/发射率借用博士时注明；不要混用博士底h=300。

能量按焓积分计算，包含瞬态储能；报告Q输入、侧底顶输出和DeltaE。先能量账本闭合，再比较T/Tmax等温线与流型。

电磁反馈F3仍为后续：T→sigma(T)→EM→Q→T,u；本轮不改A-phi分支，不静默切换统一EM候选。

## 13．执行顺序与停止条件

### 第一阶段：不重跑长算也能完成

修订状态；核对C1实测Ra/Pr；修正参数生成；审计Nu定义；导出F1空间数据；加载本包Fig.7/8；画工程T与文献等温线对照。

### 第二阶段：低成本计算定位

一维导热、C1 n80严格配置、F1 dt短对照/顶部诊断、12mm beta0。单项改变，不全排列。

### 第三阶段：形成可信比较

C1 n160；F1空间流型与末段统计；重建文献约束温度家族并做最少对照。

### 第四阶段：根据证据接F2

若已有可信闭合流，推进自洽热—流；若仍是边界噪声，先局部修正，不继续改β/μ把信号强行放大。

每类问题最多两轮有依据的修改后总结主误差；两轮无改善，冻结失败配置，给出下一可检验方案。运行预算不足标INCOMPLETE，不能放宽阈值。可以使用CPU短参考，不以GPU工程细节阻塞准确性验证。

## 14．验收与交付清单

交付一个包，包含：

1. `DOC203_STATUS_CORRECTION.md`。
2. `EFFECTIVE_SCHEME_COMPARISON.md`：C1/F1真实配置差异及默认不变证据。
3. `C1_STRICT_BENCHMARK_REPORT.md`：实测Ra/Pr、两种Nu、两速度极值、稳态和热收支。
4. `THERMAL_FLUX_FACTOR_AUDIT.md`：一维解析与热核系数/壁面位置。
5. `F1_SPATIAL_CIRCULATION_REPORT.md`：原始矢量、psi、周向流、顶部尖峰、连续性与统计。
6. `F1_IC_AND_DENSITY_LEDGER.md`：重置前后质量、实际密度控制、dt敏感性。
7. `LITERATURE_INPUT_PROVENANCE.md`：观测、重建、假设和未知分开，引用本包数据版本。
8. `VALIDATION_LEDGER.csv`：V1/V2/T1/T2等已有证据和未完成项，不只列本轮成功项。
9. 精简空间快照或抽样CSV，使外部能复核流型，不再只交全局summary。
10. 可复现命令、源码diff、配置、二进制哈希与一键分阶段脚本。

所有结果目录独立命名；原始运行数据只读保留。源码、配置、初态、参考数据任一哈希改变时，不能沿用旧的量化PASS。

### 最终判定

- 方腔标准量正确、法国有相干环流但与Fig.7不同：优先查温度结构与变黏度/边界，不否定浮力公式。
- 改温度输入后主次涡形态接近：支持热驱动结构的重要性；仍须F2验证其能自行形成。
- 只有去掉顶部异常后RMS下降，而大尺度流更清楚：属于进步，不要求速度恢复到旧噪声幅值。
- 只匹配流向、未验证速度大小：如实交付拓扑对照，不宣称定量速度复现。

**本轮重点是让证据准确、让真实流场可见，再使用文献数据约束输入。不得把数据数字化本身当成求解器修复，也不得用无限延长/加密替代明确定位。**


<!-- ===== END 204_OPHELIE_VALIDATION_CORRECTION_FIG7_FIG8_CURSOR_PLAN_2026-11-01.md ===== -->

---

<!-- ===== BEGIN 205_OPHELIE_MINIMAL_CYLINDER_CLOSED_CIRCULATION_CURSOR_PLAN_2026-11-01.md ===== -->

<a id="205-ophelie-minimal-cylinder-closed-circulation-cursor-plan-2026-11-01-md"></a>

# 205_OPHELIE_MINIMAL_CYLINDER_CLOSED_CIRCULATION_CURSOR_PLAN_2026-11-01.md

# Doc205：最小圆柱闭合环流定位实验——Cursor详细实施任务书

日期：2026-11-01  
承接Doc203/204及当前ParaView截图。本轮优先级高于Doc204的文献温度重建、F2与电磁反馈。不是取消这些任务，而是先解决“基本回流能否闭合”。

## 1．任务边界与完成定义

本轮只回答三个问题：

1. 已能运行方腔自然对流的力学装配，放入法国尺寸圆柱后，能否形成可信闭合回流？
2. 若可以，只把顶部从无滑移改成不穿透滑移后，是否出现右上角异常？
3. 若不可以，主要差异来自圆柱边界，还是密度—压力更新？

不以复现法国双环流为本轮要求。固定径向温度应产生流动，但不保证产生文献双涡。不得用“文献CSV还缺”“电磁未完善”“还没算到足够长”无限延期判断。

**完成结果必须是：闭合环流的空间与流量证据，或被受控实验缩小到明确模块的失败证据。不得只提交几个PASS标签。**

本任务是诊断与一致性检验，没有独立完整速度参考解；成功不等于法国模型定量验证通过。

## 2．代码隔离：新建一个小入口，不复制法国主程序

建议目标：`test_3d_ophelie_minimal_cylinder_convection`。

- 如果已有等价小入口，复用并改名配置，不重复新建。
- 新入口只组装圆柱、边界、基础水动力、固定温度、诊断输出。
- 从C1复用实际运行的声学/黏性类型与必要初始化，通过维度模板适配3D；不能用“相同scheme名称”代替真正同源实现。
- 从旧法国case只复用经过核对的圆柱几何/字段注册/输出工具，禁止复制数千行主程序和全部历史开关。
- OPHELIE局部目录维护新辅助代码，不修改全局声学默认、不影响其他case、main或A-phi分支。
- C1定量仍有Ra/Nu问题，但它已经能形成自然对流，可作本轮装配参照。不得把它写成“完全验证通过”。

### 2.1 第一份交付：装配差异表

写出真实类型、调用顺序和有效配置，对照：C1、旧法国F1、新最小圆柱。

必须包含：EOS、背景压、扰动动量、limiter、两种跳跃耗散、密度求和/正则化/归一化、核支持、核修正、粒子整理、黏性系数与积分、壁面压力/速度构造、体力写入、位置更新。

只复用“有出处”的组件。不要因为此前审计文档说背景压为零就相信；源代码和有效日志才是依据。

## 3．固定物理问题：所有主要比较只改变顶部边界

### 3.1 几何、坐标和物性

| 项目 | 规定 |
|---|---|
| 半径R | 0.25 m |
| 高度H | 0.185 m，作为博士补充几何，不声称从期刊图读出 |
| 坐标 | 本入口z向上，底z=0，轴x=y=0；其他Y轴STL不混入 |
| rho0 | 2750 kg/m³ |
| 动态黏度mu | 4 Pa·s，空间与时间均为常数 |
| beta | 主实验1e-5 K⁻¹；与期刊一致，旧3e-5仅作后续明确诊断 |
| 重力向量 | (0,0,-9.81) m/s² |
| Tref | 1473 K，仅作本诊断的参考，不声称文献初温 |
| EOS | 首轮复用C1线性EOS；c0固定20 m/s，后续必要时单项敏感性 |
| 初速 | 全部严格为0；质量从创建开始固定 |
| dp | 主诊断15 mm；12 mm作为候选空间敏感性，不直接上5 mm |

明确常物性：不加载mu(T)、cp(T)、sigma(T)，不做黏度截断。所有实际参数输出日志，禁止旧默认静默覆盖。

### 3.2 固定空间温度场

采用轴线上光滑的径向场：

\[
T(r)=T_{ref}+40\left[(r/R)^2-1/2\right]\quad\mathrm{K}.
\]

轴心1453 K、侧壁1493 K，截面体积平均为1473 K，径向温差40 K。它与旧F1线性r场不同，属于新的明确工程诊断输入；不是文献温度数字化场。

选择r²是为避免r在轴心不可微的问题，不是为了凑涡。beta0控制保持完全相同T，只将beta设0。

每次粒子位置更新后重新评价T；每个需要的声学步刷新浮力。不能让温度仅跟随初始粒子标签，变成另一物理问题。为避免执行顺序歧义，本轮直接由当前位置计算浮力使用的T值，并核对输出温度字段一致。

### 3.3 唯一体力

\[
\boldsymbol a_b=-\beta(T-T_{ref})\boldsymbol g.
\]

不另外加完整重力加速度。使用去除常量参考重力后的方程；默认背景压为0、EOS压力作为实际待求压力。不叠旧均值吸收或静水背景实验开关。

**实施前核对壁面构造是否依赖参考体加速度或总压。若需要总压/背景压分离，必须完整成套配置并在A/B中一致，不得只删除其中一项。**

压力不冻结，不每步按解析式覆盖。浮力通过同一ForcePrior机制进入动量且只计一次。位置、质量、密度、压力不能在运行中为降噪反复清零。

### 3.4 本轮关闭项

电磁、热源、热方程、辐射、Robin散热、Marangoni、自由表面变形、固化、变物性、叶轮全部关闭。固定T诊断不做热能守恒验收；温度由外部维持，其能量变化不能当热模块错误。

## 4．边界A/B的严格定义

### A：完全封闭无滑移

侧壁、底面、顶面：物理边界不穿透，切向速度为0。

用于回答圆柱中能否闭合回流，不是法国最终工况。所有面进入相同压力/连续性壁面构造，并进入黏性无滑移处理。

### B：只改变顶部为不穿透滑移

侧底与A完全相同；顶面保留不穿透，切向无应力。

允许采用与当前框架相容的速度镜像/ghost构造，但要说明法向与切向各自如何实现。仅从黏性接触列表移除盖板不能独自证明滑移正确。

**A/B必须共用同一个顶面几何和压力/连续性接触。B不能删除压力盖，也不能改变顶面位置、厚度、粒子间距或同时更换密度算法。**

若框架现有边界无法表达该条件，直接标`UNSUPPORTED`并实现OPHELIE局部策略；禁止以omit整个盖替代滑移。

### 4.1 物理壁面位置与角部

- 壁面位于r=R、z=0、z=H；边界粒子层及首层流体距离都输出。
- 压力、连续性、黏性必须使用一致的物理面。
- 顶侧/底侧角区可以同时接触两个面，但边界粒子不能无意重叠并被重复计数；明确ownership或面贡献分解。
- 壁面法向满足单位长度、方向一致。理想圆柱可局部计算解析法向以避开已出现的NormalFromBodyShape挂起；角部用明确分面规则，不能随便平均后用于所有接触。
- 不得运行后钳制/投影流体粒子回容器内作为正式修复。若做此消融，只能诊断并报告其额外动量修正。

## 5．初态与密度策略：不要再先制造运动再修改质量

### 5.1 初始化

规则粒子生成后检查边界覆盖、邻域数、Vol/m/rho关系和总质量。初始u=0、rho=rho0、EOS压力为0；记录体积误差而不是通过运行后质量重置隐藏它。

先做beta0静止控制。若立刻出现与旧F1相近的0.1m/s运动，不继续加热或长跑，直接定位初始化和压力壁面。

禁止：settle过程中出现大运动后再quiet-reset-mass；重写mass；运行中重置位置或强行压密度到rho0。

### 5.2 首轮只保留一套密度策略

将C1实际可运行的密度求和/归一化组合封装成明确策略，保持m固定、解释Vol更新，并记录每次重置引起的delta_rho和EOS delta_p。若这套方法不适用于当前减压/壁面组合，应在短beta0测试中暴露，不隐式混入旧法国策略。

不得同时叠加新delta-SPH。连续性更新、Riemann密度耗散和额外密度重置分开记账。无额外重置作为后续受控候选，不预设它必然更稳定。

### 5.3 精度

驱动压力标尺仅约rho0*g*beta*40*H≈2Pa（本轮beta=1e-5），不能只用rho相对误差<1%判断压力准确。

记录EOS压力计算精度、rho分辨能力和近壁压力梯度。双精度优先用于小参考；若GPU双精度不可用，可用CPU双精度短对照，但不因此停掉GPU主路径。压力RMS大本身不是错误证明，必须看其梯度/高频变化与体力的关系。

## 6．最少运行矩阵与顺序

| 编号 | 顶面 | beta | 用途 |
|---|---|---:|---|
| A0 | 无滑移 | 0 | 完全封闭静止基线 |
| A1 | 无滑移 | 1e-5 | 基本回流 |
| B0 | 不穿透滑移 | 0 | 滑移边界静止基线 |
| B1 | 不穿透滑移 | 1e-5 | 只改顶部后的回流 |

先15mm；四项使用同一个初始粒子文件/哈希。除了beta和顶部切向条件，其余全部一致。

### 6.1 时间安排

1. A0/B0先运行0.2s短检查（约20个H/c声学穿越时间），健康后至5s，检查静止与边界稳定。
2. A1先5s看方向/尖峰，再到20s看回流是否形成。A通过后再做B1。
3. 已形成回流但仍演化时，可续到80s，保存检查点；没有回流且存在显著噪声/泄漏时禁止直接续跑。
4. nu=mu/rho0；H²/nu约23.5s、R²/nu约43s，只作尺度。稳态依赖数据，不按固定80s自动PASS。

声学子步必须截断到外步剩余时间和目标时刻；输出实际时间，beta对照按施加温度后的同一物理时间配对。

### 6.2 时间步

显式黏性上限必须在常物性分支也生效。输出dt_adv、dt_ac、dt_visc、有效h、mu及生效限制。实际时间步取所有适用上限，不将隐式黏性与显式同时施加。

首轮用现有稳定CFL；只在出现异常时对同一case做dt减半。不能同时改变c0、dp和dt后归因某一项。

## 7．一开始就设置空间输出，禁止只交summary

保存t=0、0.2、1、5、10、20s；续跑保存40、60、80s。日志监测比空间快照更频繁。可以分块输出压缩VTP/CSV，控制文件量。

至少字段：ID、位置、速度、rho、mass、Vol、T、mu、实际动量压力、背景压、浮力加速度、壁距。力项可以先输出统计；定位异常后再补局部数组。

保存初态、有效配置、源码和二进制哈希、命令、硬件及精度；单独输出目录，禁止覆盖旧case。

## 8．核心验收：截面质量流量，而不是速度有多大

### 8.1 三个水平截面

取z/H=0.25、0.5、0.75，分别计算：

\[
\dot m_+(z)=\int_{u_z>0}\rho u_z\,dA,\quad
\dot m_-(z)=\int_{u_z<0}\rho u_z\,dA,\quad
\dot m_{net}=\dot m_++\dot m_-.
\]

从三维粒子插值到完整圆盘并积分；或先验证方位近似轴对称后采用2pi*r权重积分。不可只对一个子午面的箭头求和。每个截面记录插值覆盖率；空洞不能补0混入积分。

准稳态封闭容器应有上升、下降两区，其质量流量基本相抵。定义：

\[
E_Q=|\dot m_++\dot m_-|/\max(\dot m_+,|\dot m_-|).
\]

**分母接近beta0数值背景时，不计算虚假的小误差，而标记NO_RESOLVED_FLOW。**

瞬态净流量允许对应截面上下区域的质量储存变化；此时使用控制体账本验证，不把非零净流量直接当泄漏。

### 8.2 边界通量

在实际顶部/侧壁/底面分别积分局部正、负法向质量流量及净量。正负相抵的净0不等于不穿透。

对同一个控制体检查：

\[
dM_{cv}/dt+\int_{\partial cv}\rho\boldsymbol u\cdot\boldsymbol n\,dA\approx0.
\]

移动粒子计数、面插值、边界表示各有离散误差，需报告估计和覆盖率。不能只凭粒子总mass固定宣布局部连续性正确。

### 8.3 项目筛查门槛（不是文献误差界）

- 激活流在三个截面出现可辨认的上升/下降区，并且空间上连接成返回通道。
- 准稳态三个截面的EQ目标<=5%；若瞬态则用储存项闭合，并明确未稳态。
- 各面正负绝对通量之和，相对典型内部循环通量目标<=5%；无可靠面通量测量则标未覆盖。
- 后续总质量相对变化目标<=1e-6（固定mass求和允许舍入），无质量重置。
- 相干子午运动显著高于beta0背景，建议至少10倍；这个比值只是必要条件，不单独授予环流PASS。
- 末段平均循环通量趋势<5%、相对波动<10%可作为初步稳态指标；未满足时可报告瞬态可信环流，不能称稳态通过。

在已有空间覆盖不足导致误差时，先改采样/诊断；不得通过丢掉近壁数据来降低通量误差。

## 9．空间流型、尖峰和显示检查

1. 至少4个方位原始子午面切片，评估非轴对称与utheta。
2. 同时画等长箭头方向图、按速度着色的矢量图、保留真实比例的箭头图。不得仅用强尖峰缩放的Glyph判断内部静止。
3. 输出全域、核心、近侧壁、近顶、角部的RMS、p95、p99、max；报告高速度粒子位置。
4. Glyph伸出边界不是穿墙证据；用原始位置、面法向流量和越界距离判断。
5. 流函数作为辅助，使用轴对称r权重和物理坐标；反算速度与原始速度比较。若散度过大，不得用光滑psi强行画出闭合涡。
6. 可以画流线，但不能靠种子选择和插值跨空洞造出漂亮闭合线。主证据仍是速度/流量。
7. 不要求速度达到sqrt(g*beta*deltaT*H)；此标尺不是黏性流的目标值。

## 10．失败决策树：最多两次有依据的修改

### 情况1：A0或B0产生明显运动

停止其对应受驱动长算。查初态邻域、密度正则化造成的delta_p、常量压力梯度误差、壁面参考体力/法向。beta0失败时不能以beta>0速度更大就通过。

### 情况2：A1有可信闭合回流，B1没有或顶部尖峰明显

优先定位顶部速度镜像、法向、角部重复接触、压力/连续性盖板条件。只修顶部局部策略，侧底及内部算法冻结。修复后重跑B0短检验+B1；不回头重写浮力。

### 情况3：A1、B1都不闭合，而相同装配方腔有回流

先比较曲面与角部几何/法向/邻域覆盖，以及三维与二维差异。可以在同一小入口加入矩形容器配置，保持三维和物性/驱动相近，用来区分“3D共同问题”与“圆柱曲面问题”；这仅在本分支触发，不首轮新增大矩阵。

### 情况4：几何检查正常，但密度重置不断制造压力跳变

只在A上做一次有额外重置/无额外重置的对照；保留Riemann连续性与耗散，关闭额外重置不等于关闭质量方程。若无重置发散，报告其不稳定，不直接判其物理正确。根据delta_p和力项定位是否需要更一致的密度控制。

### 情况5：局部尖峰随dt减半明显改变

先修时间步/步进/力更新。不能以升黏度掩盖异常。只有确认真实黏性步进正确后才讨论数值耗散。

### 情况6：两次修正仍无决定性改善

提交空间数据、控制体通量、实际算子和两次对照的结论，明确剩余候选。停止继续叠开关，下一决策可采用已验证的替代壁面或独立轴对称参考；不返回旧B2不相容解析速度测试。

## 11．通过后怎么回到法国主线

A1/B1形成可信回流后，保留该小入口作为回归。将成功的组件组合逐项接回原法国case：

1. 同样固定工程T、常mu；必须先得到与小入口一致的流量/流型。
2. 加载文献约束重建T，明确其温标/补全假设；不把固定T等同热方程验证。
3. 恢复mu(T)，用已有变黏度证据支持；不同时恢复全部热/电磁。
4. 开温度方程、固定非均匀Q，检查热流与温度自洽演化。
5. 最后接回电磁反馈。

若小入口正常而旧法国入口失败，以装配差异为根因线索，不再次修改同一个已通过小入口的浮力公式。

## 12．执行纪律与交付

### 必须提交

- `ASSEMBLY_DIFF.md`：C1/旧F1/新入口的实际调用与差别。
- `RUN_MATRIX.csv`：A0/A1/B0/B1及条件触发的补充运行。
- `MASS_FLUX_AND_STORAGE_REPORT.md`：截面/壁面通量、储存、覆盖率。
- `CLOSED_CIRCULATION_REPORT.md`：原始速度、方向、返回通道及不确定性。
- `FAILURE_LOCALIZATION.md`：若失败，哪项证据排除了什么，还不能排除什么。
- 初态、关键空间快照、诊断脚本、有效配置、源码diff和复现命令。

状态分开：WIRING_OK、QUIET_CONTROL_OK、BUOYANCY_RESPONSE、CLOSED_RETURN_CONFIRMED、STEADY、QUANTITATIVE_REFERENCE_VALIDATED。最后一项本实验没有独立参考时必须为NOT_COVERED。

不要交付只有passed=true的JSON；不要以“平均速度增长”“图像好看”或“文献输入仍缺”为结束标准。GPU上运行不了的项目单列工程阻塞，短CPU参考可继续，不默认切换整个项目。

### 本轮结束时必须回答的四句话

1. 全封闭无滑移圆柱是否形成可信闭合回流？证据在哪？
2. 只改顶部滑移，回流、泄漏和尖峰如何变化？
3. 若失败，边界、密度—压力、时间积分中哪个得到直接证据支持？
4. 哪个组件组合可以接回旧法国case，或下一次唯一需要改变的项目是什么？

**本轮成功不要求法国双涡、不要求电磁、不要求完美高阶精度；要求用受控的小问题证明基本回流，并将失败定位到可行动的模块。**


<!-- ===== END 205_OPHELIE_MINIMAL_CYLINDER_CLOSED_CIRCULATION_CURSOR_PLAN_2026-11-01.md ===== -->

---

<!-- ===== BEGIN 206_OPHELIE_BUOYANCY_CLOSURE_AND_THERMAL_HANDOFF_CURSOR_PLAN_2026-11-02.md ===== -->

<a id="206-ophelie-buoyancy-closure-and-thermal-handoff-cursor-plan-2026-11-02-md"></a>

# 206_OPHELIE_BUOYANCY_CLOSURE_AND_THERMAL_HANDOFF_CURSOR_PLAN_2026-11-02.md

# Doc206：浮力模块阶段收尾、剩余问题定位与热场主线交接

日期：2026-11-02。执行对象：Cursor。  
依据：Doc205_French_Reconnect_Handoff_2026-11-02.zip及本轮审阅。  
优先目标：在有限计算预算内，把“能产生基本自然对流”正式收尾，并把热场和法国主次环流对齐推进下去。

## 0．执行摘要：不要再重新开始

当前已观察到：最小圆柱回流、法国入口回流、无质量重置、变黏度及热/电磁接入后仍有大尺度返回通道。浮力公式和现有Riemann WCSPH路线不再作为首要怀疑对象。

本轮不更换压力算法、不新增delta-SPH、不重写GPU隐式黏性，不重新开发一个法国大case。保留已有入口、算子和有效结果。

本轮工作分两条有明确结束条件的线：

- **H线：收尾水动力。** 修正通量后处理；确认局部连续性、边界泄漏、时间步/空间敏感性；冻结可用的自然对流装配。
- **T线：继续主线。** 用解析热流检验校准边界面积及散热离散，再固定热源推进自洽温度场。

H线通过后立即交接T线，不因Fig.8尚未匹配而再次宣称浮力未完成。反过来，不能用回流图代替边界/连续性验证。

**“阶段收尾”只适用于明确配置范围，不是承诺所有浮力工况及法国模型都已全面验证。**

## 1．必须纠正的事实与状态

| 项目 | 已有证据 | 本轮要求 |
|---|---|---|
| 大尺度回流 | t50原始VTK显示轴侧下降、底部向外、靠外热区上升、顶部返回 | 保留，不再从零研究浮力公式 |
| 最外冷壁下降带 | 切片可见 | 与上升带分区，不混成一个“壁面平均” |
| 总质量 | t50 mass_init=mass_end≈97.230372kg，quiet-reset关闭 | 作为有效进展；另补局部控制体账本 |
| 通量 | slab选取全厚2b却除b | 绝对流量约多一倍，先修脚本重算旧CSV |
| EQ | 相同因子抵消 | 不因系数修正全盘否定旧EQ，但需检查采样宽度 |
| 热边界面积 | 顶−10.5%、侧−10.3%、底−6.4% | 解析热流验证，不先调h/功率 |
| t50热场 | 50kW输入，末态散热约108kW，平均降温33.3K | 是守恒的冷却瞬态，不是热稳态 |
| Fig.8“热点位置” | 部分参考取自人为补全场 | 不当作作者精确热点 |
| C1标准方腔 | Ra/Pr、Nu分歧及部分速度误差未完全关闭 | 列独立验证欠账，不默认为已通过 |

保留历史报告，在新状态表写明修订原因，不覆盖原始数据。

## 2．冻结实际成功组合

建立 `FROZEN_HYDRO_CONFIG.json` 与 `ASSEMBLY_VERIFIED.md`，从真实源码及运行命令提取：

- 几何、轴向、R/H/dp、粒子生成方式、物理壁面位置、法向来源。
- rho0、mu、beta、声速、精度、核定义与h/dp、邻域数。
- EOS、背景压、PerturbationMomentum、两种Riemann耗散、limiter、核修正。
- 密度求和、正则化及归一化的精确调用；不能把C1与F1不同的密度控制写成同一套。
- 顶部压力盖、滑移标签、黏性接触列表、角部规则。
- ForcePrior写入、清除及浮力刷新顺序；温度空间固定或随热方程演化的区别。
- quiet-reset-mass、完整GravityForce、均值吸收等历史诊断是否关闭。

当前t50候选：R=0.25m、H=0.185m、dp=0.015m、c0=20m/s、beta=1e-5；变黏度为借用曲线，不能作为期刊精确物性。固定温度/常mu的小圆柱保留为快速回归。

冻结不代表不准修bug；每次修正生成新配置标识与源码哈希，不静默改变旧基线。

## 3．H1：修通量统计，尽量零仿真成本完成

### 3.1 明确旧脚本错误

脚本 `section_mass_flux_from_field_csv.py` 令band=1.5dp，选择abs(z-z0)<=band，取样全厚2band，却使用rho*uz*Vol/band。

第一步改为全厚归一化：

\[
\dot m(z_0)\approx\frac{1}{2b}\sum_{|z_i-z_0|\le b}\rho_i u_{z,i}V_i.
\]

标记为SLAB_AVERAGE，不冒充无误差的数学截面。若取样层被物理域截断，按实际积分权重重新推导；本轮三个内部截面不应截断。

先用均匀人工速度场检验脚本：圆柱截面积A=pi*R²，解析流量rho*U*A。选用已知体积权重的规则点云；离散几何误差另列，不能误当公式系数错。再用具有正负分区的人工场检查两侧流量与EQ。

重算已有A1、B1、法国Phase1及t20/t50 CSV。旧数据若仅有VTK且缺rho/Vol，不猜常数冒充重算；从服务器导出原始CSV。**不用因后处理错误重跑长时仿真。**

### 3.2 减少宽层平均的误导

dp15mm旧全厚45mm，占H约24%。依次比较b=0.5dp、dp、1.5dp。小带宽采样层数/覆盖不足时标UNDER_RESOLVED，不硬挑EQ最小的带宽。

输出每个带宽的取样数、有效体积、平均速度、正负流量及EQ。若可行，追加圆盘网格插值积分作独立对照；必须包含全圆盘、覆盖率与缺失区，不能将空洞自动填0。

### 3.3 识别主流与外围下降带

至少输出z/H=0.25、0.5、0.75的径向分箱剖面：rho、T、mu、uz及流量贡献。初始分区可用r/R<0.4、0.4–0.7、0.7–0.9、0.9–1，并以实际分辨率检查最外层是否足够采样。

最终判断依据是径向符号变化，不强迫每个分区都满足预设方向。报告上升峰位置、外侧零交叉及近壁下降带宽；不把r/R>0.7合成单一“壁面上升”。

## 4．H2：局部连续性与边界是否真正闭合

### 4.1 控制体质量账本

取下半圆柱和下四分之一圆柱两组固定空间控制体。至少使用连续多帧：

\[
R_M=\frac{M(t_2)-M(t_1)}{t_2-t_1}+\overline{\int_{\partial CV}\rho\boldsymbol u\cdot\boldsymbol n\,dA}.
\]

包括内部水平截面及物理壁面的通量。不要用全体粒子的mass恒定代替局部守恒。粒子穿越CV边界会导致计数跳动，使用时间窗积分/平均，并给出离散估计，避免对单帧差分下结论。

用质量流量尺度归一化；流量接近beta0背景时标NO_RESOLVED_FLOW。不得用接近零的净通量作分母。

### 4.2 不穿透检查

顶、底、侧分别输出正/负法向质量通量、绝对值和净值。面值由一致的近壁插值获得；同时输出首层粒子法向速度，但两者不能混为一谈。

角部单列覆盖率与法向规则。Glyph伸出容器不是穿墙证据；以原始位置和真实壁面距离检查越界。

### 4.3 本项目阶段验收门槛

以下为项目工程筛查门槛，不是文献理论误差界：

- 准稳态截面EQ目标<=5%；下截面允许先<=10%，须同时说明带宽敏感性与储存项。
- 瞬态未稳态时不强求EQ为零，改用控制体储存闭合，归一化残差目标<=5%。
- 物理壁面绝对泄漏通量相对典型内部循环通量目标<=5%；面插值不可靠则标UNCERTAIN，不能通过净值相消验收。
- 运行中不重写mass；总质量变化目标<=1e-6或有记录的浮点求和误差包络。
- 明显越界不得超过0.5dp；严格面外点数量/位置也要报告。此门槛不等于允许真实泄漏。
- 若相同离散尺度的通量不确定性已大于门槛，应加密诊断或做独立面重建，不能只放宽门槛。

## 5．H3：最少补跑集，不再展开全排列

优先复用全部旧结果及服务器检查点。仅补缺少证据的项目：

| 任务 | 配置 | 作用 |
|---|---|---|
| H0 | 最终成功装配beta=0，5s | 确认不靠0.2s短静止掩盖慢漂移 |
| Hdt | 固定空间工程T、常mu，15mm，dt上限减半 | 确认回流与顶部尖峰对步长的敏感性 |
| Hdp | 同一固定T/常mu，12mm | 一次空间敏感性，不同时改温度/物性 |
| Hcv | 最终选定配置多帧导出 | 完成控制体账本；能从旧场导出则不新增运行 |

A/B已经有效的结果不全部重跑。Hdt/Hdp比较与基线完全相同的物理时刻；至少到20s，若用作稳态结论则需末段趋势依据。不能用5s加密结果比较旧20s。

固定T采用Doc205小入口的光滑r²场或现有法国工程场之一，必须统一选择并写出公式。旧小入口与法国R2线性r场并不严格相同，所以旧幅值“同量级”不是精确复现证明。本轮不要混用这两种输入来判断误差。

建议一致性目标：循环通量、主要回流方向及核心速度dt减半变化<=5%，dp15→12变化<=10%。未达到时只追踪产生变化最大的模块；不全面重写。

常物性显式黏性CFL同样生效，记录dt_adv/dt_ac/dt_visc及实际限制。现有Implicit不作为本轮前置；不叠显式与隐式黏性。

beta0背景相对受驱动相干流建议低于10%；若达不到，先定位背景而不是放大beta。对1e-5主线保持真实小驱动，不用强驱动代替正式验收。

## 6．H4：旧失败原因做有限归因，不再为“彻底解释历史”无限计算

交付 `ROOT_CAUSE_EVIDENCE.md`，使用四级标签：CONFIRMED、SUPPORTED、POSSIBLE、NOT_TESTED。

1. 已确认：旧初始化曾在均匀T时产生明显运动；旧quiet-reset改变质量；旧验收不足；通量厚度系数错误。
2. 已支持：简化装配、边界/压力配套、取消质量重置后能够回流。
3. 尚未单独证明：哪一条顶部/法向/密度改动独立导致回流恢复。

先比源码和有效调用，不盲目重跑。确需因果归因，最多两项对照：

- 在当前小入口中只比较已有的两种顶部实现，保持其余完全一致。
- 或只比较当前与旧密度控制，先短静止再受驱动短窗。

不要同时退回旧初态、改mass、改壁面及改T。旧错误初始化无需重新启用到生产；已有日志足够支持其风险结论。

如果两次对照仍不能唯一归因，按组记录“装配一致性修复有效”，保留未知项即可。阶段收尾不要求证明某一行是唯一根因，也不允许写成已找到唯一边界bug。

## 7．浮力模块何时正式收尾

输出 `BUOYANCY_STAGE_ACCEPTANCE.md`，逐项列结果与证据文件：

- 基础公式与更新接线：复用O0/浮力审计，确认版本覆盖。
- beta0静止控制：最终装配短中窗通过。
- 固定非均匀T：原始矢量、三截面正负流、空间返回通道一致。
- 局部连续性：控制体储存与面通量在误差范围闭合。
- 固体边界/滑移：无足以解释回流的泄漏，尖峰不主导相干运动。
- 时间/空间敏感性：一组有效受控对照，误差或限制明确。
- mass固定且无quiet-reset。

全部满足时授予：**ACCEPTED_FOR_FRENCH_ENGINEERING_COUPLING**。

同时明确不覆盖：任意高Ra自然对流、完整2024方案、法国双涡定量、所有变黏度曲线、EM精度。没有独立速度参考解时，不标ALL_PHYSICS_VERIFIED。

C1此前Ra/Pr、Nu问题仍是独立定量验证欠账：复用Doc204参数修正方案，按精确Ra/Pr检查一组标准方腔；可与T线排队运行。没有完成时写VALIDATION_DEBT，不假装标准定量验证已结束。若该测试揭示重大速度/热流错误，应重新评估受影响范围，而不是机械维持阶段PASS。

**若H线只有一项门失败：最多两轮单项修改。仍失败则明确阻塞项，不扩大到整套浮力研究。**

## 8．T1：热边界先用解析热流校准

本部分不返工水动力，优先使用现有T3/单面诊断入口，不新建多个大case。

### 8.1 已有面积问题

R=0.25、H=0.185时，几何目标：

\[
A_{side}=2\pi RH\approx0.290597\ \mathrm{m^2},\quad
A_{top}=A_{bottom}=\pi R^2\approx0.196350\ \mathrm{m^2}.
\]

旧离散面积分别约0.26055、0.175725、0.183825m²。面积偏小不能解释散热过强；补足面积可能让散热更多，因此不把area-normalize当降温调参。

### 8.2 最小热流矩阵

冻结位置，关闭流动和热源，给定已知均匀流体温度，读取瞬时边界功率；若温度被外部固定，只检验瞬时功率，不能拿外部覆盖后的DeltaE做能量守恒。

- 侧面对流单独开启，其他面关。
- 底面对流单独开启。
- 顶面对流单独开启、辐射关。
- 顶面辐射单独开启、对流关。

参考：Q=h*A*(T-Tinf)；辐射Q=eps*sigmaSB*A*(T^4-Trad^4)，温度全部Kelvin。至少两个温度检验辐射四次方关系；确认没有重复乘面积、体积或质量。

各项相对解析误差目标<=2%。逐面记录原始/归一化面积、角部归属、实际温度取值位置、符号和单位。不允许同一热边界单元因primary-face规则遗漏必要面，或角部双算；采用何种面积分配必须可追溯。

均匀场测试通过后，加一项有导热阻力的平板Robin测试，检查真实边界温度与壳层内部粒子温度的区别。否则单纯面积正确仍不能证明冷壁层散热准确。复用既有导热模块，参考用一维串联热阻解析解。

### 8.3 修复的回归范围

只修OPHELIE边界面积/通量适配，保留旧配置。短复跑一个固定热源热—流候选，不重跑电磁全套。检查新边界是否改变流型和温度趋势。

## 9．T2：冻结50kW非均匀热源，推进自洽热场

### 9.1 建立可追溯固定Q基线

从当前已生成的EM快照选择一份Q(r,z)，记录来源、时间、电导率、频率、电流、原始功率及映射；将离散玻璃吸收功率固定50kW。关闭live EM更新，其他成功水动力组件冻结。

固定Q是隔离反馈的诊断，不声称文献最终热源准确。QCSV与SPH采样归一化分别报告，不能仅凭总功率正确证明空间形状正确。

### 9.2 初态与物性

原t50是从特定初态演化出的冷却瞬态，不能直接当稳态。允许复用检查点，但比较时记录累计物理时间、T/mu/u以及初始焓。

期刊主线beta=1e-5、rho0=2750；mu/cp/k曲线来源逐项列清，借用CES/博士的曲线不改标签冒充期刊。先不人为提高黏度陡度来凑双涡。

使用已有mu(T)配方可继续推进工程线，但定量对齐需保留物性不确定性。冷区高黏导致显式步长过小时，复用已验证隐式或CPU参考，不先修改mu上限压低真实阻力。

### 9.3 热平衡判断

每个时间窗报告：Qin、Qside、Qbottom、Qtop_conv、Qtop_rad、dH/dt和残差。

变cp采用H(T)=integral(cp(T)dT)；温度重映射/质量或物性曲线改变导致的能量变化也计入，不以巨大交换量归一化掩盖小净热量误差。

阶段目标：累计能量残差<=输入能量的1%（冷却时另报交换与净储能尺度）；若要称热稳态，末段|Qin-Qout|/Qin<=5%，Tmean趋势与关键剖面也稳定。残差小但Qout约108kW仍是冷却，不能标稳态。

热边界核对正确后，若持续冷却：先检查真实冷壁温度、导热阻力、未解析冷层、物性和热源位置。不能立刻把Q调到108kW、把h改30或把eps关0来关闭文献工况问题。softcool只保留工程基线。

### 9.4 分辨率与运行预算

15mm无法保证解析mm级冷层；用实际alpha估算sqrt(alpha*t)，并检查边界层内粒子层数。

H线空间敏感性通过不等于热边界层也收敛。T线若出现明显近壁误差，先同一固定Q/物性下做一次12mm或10mm短对照；不立即启动百万粒子长跑。

只有热流、温度剖面、主次流型显示受益，才扩大长时/细网格预算。声速、dp、热边界参数不同时改变。

## 10．文献对照：只与真实读图约束比较

复用Doc204原始等温线CSV，不能把人为补全T场的热点或675K温差当作者精确数据。

必须区分：DIGITIZED_CONTOURS、RECONSTRUCTED_ASSUMED_INPUT、SIMULATED_FIELD。

- Fig.8：比较0.5–0.9等温线在归一化坐标中的位置、中央交点和热核包络；Tstar为T/Tmax。
- 当前0.21中位、0.48的p95温度残差照实保留，不能只用相关系数0.8盖过它。
- 增加空间等温线距离，逐等级等权或按弧长权重，避免同一条曲线上大量相关采样点虚增证据。
- 图上的缺口及未标定外曲线不补成观测；源图没有明确Tmax时不宣称绝对K定量匹配。
- Fig.7：重点比较内部上升带、主回流和最外壁下降带，不能仅用“壁升轴降”通过双涡。
- 若有psi，必须反算速度检验；外围下降带与主流之间是否形成独立闭合次涡，需要连续区域证据，不是单帧一排向下箭头。

不要求本轮一定达成法国双涡；要求将偏差定位到热场/近壁解析度/物性/热源中的可行动项目，不再泛化为浮力公式失效。

## 11．EM只维持可用接口，不抢占本轮收尾

当前spatial sigma反馈4次成功是接线证据；dp_em=50mm、电磁误差与更新周期仍有独立限制。

50kW归一化意味着固定吸收功率，不能同时宣称固定线圈电流下的总热功率反馈已验证。下一轮恢复live EM时，分别比较Q形状与原始功率，保存归一化系数。

本轮不改EM核心、不碰A-phi分支，不为TEAM7等未完成项阻塞固定Q热—流主线。

## 12．执行顺序与时间限制

以可获得GPU时段安排，以下是优先级/预算目标，不是未经测量的运行时间保证：

1. **第一工作块：优先当天完成。** 冻结配置、修通量、重算已有数据、空间分区、状态修订。此块原则上不需要长仿真。
2. **第二工作块：集中一个GPU队列。** 只补H0/Hdt/Hdp/Hcv缺失证据；同时复用短热边界解析测试。CPU只作短参考，不重复全矩阵。
3. **第三工作块：作收尾决策。** H满足则发布浮力工程接入验收，立即推进T；H失败只修一个主阻塞，最多两轮。
4. **后续主线：固定Q热—流。** 按检查点续跑和时间窗评估，不无条件再跑到某个秒数。

每次新增运行前写明：它要区分哪两个解释、改变哪个变量、使用哪个现有参考、什么结果会改变决策。回答不了这四点就不运行。

不再要求全部历史小case绿灯才推进；但新证据若揭示真实守恒/泄漏或步进错误，必须处理，不能因赶工掩盖。

## 13．运行脚本与交付规则

修复脚本将失败吞掉的问题：旧脚本存在保存bin_exit后仍继续、后处理`|| true`等。新的主验收脚本必须：

- 仿真退出非零时返回失败，不继续使用旧输出。
- 输入/源码/二进制/初态/参考哈希作为缓存键；输出目录独立。
- 后处理失败标诊断缺失，不产生总PASS。
- 不用“最新文件”盲选快照；按case_id和实际时间匹配。
- 缺少ref时不使用硬编码旧EQ替代实际参考。
- 控制台区分计算完成、诊断完成、工程验收、定量验证四种状态。

最终包至少包含：

1. FROZEN_HYDRO_CONFIG.json、ASSEMBLY_VERIFIED.md。
2. 修正后的通量脚本、人工场验证、旧新数值对照。
3. ROOT_CAUSE_EVIDENCE.md，不能伪称唯一根因。
4. BUOYANCY_STAGE_ACCEPTANCE.md：每项通过/失败/未覆盖、适用范围。
5. THERMAL_BOUNDARY_ANALYTIC_REPORT.md：面积及逐面功率。
6. 固定Q热—流的温度、热功率、径向速度分区与文献约束图。
7. VALIDATION_DEBT.md：C1定量、高Pr、黏度曲线、EM精度等独立欠账。
8. 关键原始CSV/空间快照、命令、配置、源码diff和哈希。

## 14．最后必须给用户的结论

- 当前成功组合是否可以冻结用于法国工程热—流计算？
- 是否仍存在足以影响回流的边界泄漏或压力—密度异常？
- 旧失败原因哪些已证实，哪些只能作为组合性解释？
- 剩余误差主要属于水动力、热边界、温度结构还是文献输入不确定性？
- 下一步唯一优先动作是什么，哪些测试不再重复？

**不把“所有未知都研究完”作为结束标准。用足够的守恒、边界和敏感性证据关闭本配置的浮力工程接入，把后续资源转向法国温度场和主次环流对齐。**


<!-- ===== END 206_OPHELIE_BUOYANCY_CLOSURE_AND_THERMAL_HANDOFF_CURSOR_PLAN_2026-11-02.md ===== -->

---

<!-- ===== BEGIN 207_OPHELIE_ROBIN_VALIDATION_CALIBRATION_AND_BUOYANCY_CLOSE_CURSOR_PLAN_2026-11-02.md ===== -->

<a id="207-ophelie-robin-validation-calibration-and-buoyancy-close-cursor-plan-2026-11-02-md"></a>

# 207_OPHELIE_ROBIN_VALIDATION_CALIBRATION_AND_BUOYANCY_CLOSE_CURSOR_PLAN_2026-11-02.md

# Doc207：近壁热传递验证、法国环流校准与浮力阶段收尾

日期：2026-11-02。接续Doc206结果及随后讨论。执行对象：Cursor。

## 0．最新授权与执行优先级

用户已明确：甲方当前需要法国自然对流主次环流的趋势与形态；文献未提供的物性、温标和边界参数允许校准，不必等待作者原始数据。

本轮主线：**冻结已能回流的水动力组合 → 独立验证近壁导热—散热 → 接回法国固定Q热—流 → 有限校准缺失参数，比较主次环流。**

后处理和C1严格定量作为有限补充任务，不再阻塞所有主线。除非发现真实穿墙、明显非物理压缩或数值发散，不重新研究浮力公式、改压力路线、上delta-SPH或重写GPU隐式。

本任务书更新Doc206的执行优先级：不是要求将所有CV指标做到完美才能推进。仍需如实保留未完成项，不能用调参掩盖已确认的实现错误。

## 1．现状及必须澄清的表述

- 大尺度回流已实现；此前看到的外围下降带不因缺少严格标签而被否定。
- 当前没有证据证明粒子明显穿墙或自由液面跳出。顶部是平坦不穿透滑移模型，并非本轮求解自由液面变形。
- H2的IDW/CV失败属于未确定诊断问题；不预先认定全是后处理错误，也不预先认定流体在漏。
- 法国双环流主要为内部主环流与靠冷壁外围环流，不是强制上下堆叠两个涡。
- C1已能形成自然对流，但严格Ra/Pr匹配、Nu定义及部分速度指标仍有欠账，不等于方腔完全失败。
- Doc206的T1a证明均匀温度下面积与公式接线正确；T1b壳层Robin恒等式没有独立证明温度层准确。
- 现有108kW散热与50kW输入是能量守恒的冷却瞬态；在近壁温度表示未验证前，不等于物理模型已完全正确。

## 2．交付目标与不再重复的工作

本轮最低交付：

1. 一维瞬态Robin独立参考与同热核数值比较，给出表面温度/剖面/热流/累计散热。
2. 根据结果保留当前热边界，或新增一个OPHELIE局部且验证过的候选，不改公共默认。
3. 法国固定Q热—流中生成主回流、外围流向与温度约束的对照。
4. 一张参数来源/可校准表及有限校准结果。
5. 修正H2带符号均值，验证测量方法；浮力工程可用结论与未覆盖项分开。

复用Doc205/206小圆柱、法国入口、O0证据、H0/Hdt/Hdp及等温热流结果；不重跑全矩阵，不再以面积双计/live EM为主要过冷解释，不重复文献图数字化。

## 3．任务A：真正独立的瞬态Robin验证【最高优先】

### A1．最小物理问题

用现有热验证入口增加局部子配置；不存在则新增小型热验证入口，不能复制法国大程序。

二维矩形用于承载相同SPH热核，温度仅沿x变化：0<=x<=0.12m，0<=y<=0.12m；x=0对流散热，其余三面绝热。流动冻结，重力/浮力/热源/辐射关闭。

固定常物性：rho=2750kg/m³，cp=1150J/(kg K)，k=4W/(m K)；初温Ti=1473K，Tinf=303K，h=300W/(m² K)。温度通过热方程演化，禁止每步覆盖温度后再宣称验证导热。

验证入口与生产热核保持同源。若生产实现依赖圆柱壳层面积分类，抽出“面积—代表体积—通量施加”的局部接口，平板使用明确面面积；不能在验证中另写一套热通量实现而不检查生产路径。

### A2．解析参考与稳定计算

在扰动未到达远端的时间内使用半无限平板解；x方向指向物体内部：

\[
\alpha=k/(\rho c_p),\quad \eta=x/(2\sqrt{\alpha t}),\quad b=h\sqrt{\alpha t}/k,
\]
\[
\theta(x,t)=\frac{T-T_\infty}{T_i-T_\infty}
=\operatorname{erf}(\eta)+e^{-\eta^2}\operatorname{erfcx}(\eta+b).
\]

其中erfcx(s)=exp(s²)*erfc(s)，使用scipy.special.erfcx避免溢出。t=0单独返回Ti，不除0。

表面温度和热流：

\[
T_w=T_\infty+(T_i-T_\infty)\operatorname{erfcx}(b),\qquad
q''=h(T_w-T_\infty).
\]

累计散热可由对q''独立高精度时间积分求得，或由温度剖面的能量亏损积分交叉核对。参考脚本不能调用被测SPH热核。

来源指引：B.S. Yilbas等，Convective Heating of Solid Surface: Entropy Generation Due to Temperature Field and Thermal Displacement, Entropy 5(5), 467起；https://www.mdpi.com/1099-4300/5/5/467 。实施时以原式和上述等价变换复核；不是法国文献的温度解。

按此常物性独立计算：alpha≈1.2648e-6m²/s；t=10s时Tw≈1189.94K；t=50s时Tw≈969.20K。它们仅用于本平板验证，**不能直接拿圆柱壳层平均温度与Tw相减，称为仿真误差**。

确认整个比较窗远端不受热扰动：例如50s时L/(2sqrt(alpha*t))约7.5，半无限近似适用。若延长到远端受影响，改用独立一维有限域参考，不能继续沿用半无限解。

### A3．分辨率与运行集

先复现生产等效粗层dp=15mm；增加局部细级7.5/3.75mm，必要时1.875mm用于确认趋向解析解。二维小域细化成本远低于法国3D全局细化。

不要求一次把所有级跑满。先粗/细两级判断误差方向；只有无法解释时补中间级。若法向与切向分辨率必须相同，采用规则二维粒子即可，不开发多分辨率框架。

输出t=1、5、10、20、50s温度剖面、真实壁面重建温度、瞬时通量、累计散热与能量。热扩散时间步随dp²限制，并做选定细级一次dt减半；不能把时间误差归到边界。

### A4．比较定义及项目门槛

- 解析T(x_i,t)与粒子代表点温度比较，不把第一层粒子温度当解析Tw。
- 数值Tw单独输出其重建方法；原法没有Tw自由度就标NO_INDEPENDENT_WALL_TEMPERATURE，不硬命名shell为wall。
- 比较q''、累计散热、以初始温差1170K归一的温度误差；远离边界的常温区域不能掩盖近壁误差。
- 细级目标：累计散热<=3%、t>=5s热流<=5%、近壁剖面归一RMS<=3%。这是项目目标，非文献误差界。
- 特别报告粗15mm误差及网格趋势；细级通过不意味着生产15mm自动准确。
- t趋零时通量/温度尚未解析，单独展示，不通过随意删点让结果好看。

## 4．任务B：必要时修近壁表示，不先修改文献h【由A触发】

### B1．决策

- 若A证明当前壳层方法在生产分辨率下足够准确：保留，不开发新边界；过冷转查物性/热源/状态。
- 若A显示粗层系统性高估散热：优先新增局部壁温/ghost候选，保留旧实现作对照。
- 若误差来自内部热扩散核：先修导热，不把内部错误吸收到等效h。

### B2．候选实现原则

一种可选局部表达为：

\[
q''=k(T_i-T_w)/d,
\]
\[
q''=h(T_w-T_\infty)+\varepsilon\sigma_{SB}(T_w^4-T_{rad}^4).
\]

无辐射时可解析消元；有辐射时用有界、单调的标量求根，输出残差和失败状态。Ti代表何处温度、d对应什么位置必须由离散几何确定。

**这只是候选表达，不能直接套上就判修复：**

- d不能凭经验固定2dp或取sqrt(alpha*t)来凑热流。
- 若现有热核已解析粒子到壁面的导热，这段阻力不能再重复串联。
- 分布式壳层体积源与面通量模型需要一致替换，不同时施加旧散热和新散热。
- h已包含的坩埚/接触热阻不能重复添加；这里只处理数值自由度到物理面的表示。
- d→0、网格加密时应恢复相应物理Robin条件；候选仍须通过任务A，而不是只满足自身恒等式。

新增例如 `--thermal-wall-model legacy-shell|wall-temperature`，实际旗名与现有配置兼容。仅OPHELIE局部；当前生产默认先不静默更改。

### B3．角部和能量

逐面面积使用几何目标；角部允许贡献多个物理面时必须给各自面积，不能漏面/重复面。所有q''乘面积后再换算到粒子热容，单位从W/m²到W到K/s逐步核对。

能量账本与实际施加源同一来源，但另保留独立解析功率检查。两边复制同一错误公式得到残差0不能算验证。

## 5．任务C：接回法国固定Q热—流【主线】

### C1．冻结输入

沿用成功水动力装配与Doc206固定Q快照，记录Q文件哈希、归一化规则和50kW玻璃吸收功率。先关闭live EM；不改压力、密度控制、beta和黏性算子。

保留R=0.25m、H=0.185m、beta=1e-5、rho0=2750等已确认配置。保持坐标映射；原自然对流case用z，其他Y轴STL不得混用。

选择一个初态与一致的检查点。旧热边界检查点可以用于探索，但比较时要说明旧模型塑造的初温场；热模型对照需同初态，不拿不同历史运行直接比较终点。

### C2．最少对照

1. 已有legacy-shell结果只读复用。
2. 新模型在相同初态做5s稳定性，再20s。
3. 只有温度剖面/热流仍明显演化且计算健康，才续跑较长窗；不默认再跑一批50/100s。

比较：Tmin/mean/max、壁面温度、壳层温度、各面功率、累计焓、主回流与外围下降带。新模型如果使热流增加也如实报告，不按是否“更暖”选模型。

### C3．热状态判断

瞬态应满足dH/dt=Qin-Qout，变cp使用焓积分。累计能量残差目标<=输入能量1%，并额外报告净储能尺度，避免大交换量掩盖误差。

热稳态需要末段净热量/输入量<=5%、温度剖面和流量趋势稳定。50kW输入、108kW输出且温度持续下降只能叫冷却瞬态。

不要把稳定在约1500K设为硬编码强制条件；通过合理热平衡和允许校准的缺失参数去接近目标。

## 6．任务D：有限参数校准，尽快对齐甲方关心的环流

### D1．参数分类与授权

用户允许校准文献缺失/含糊参数。建立 `PARAMETER_PROVENANCE_AND_CALIBRATION.csv`：名称、单位、当前值、来源页码、是否可校准、范围依据、每次使用值。

| 类别 | 处理 |
|---|---|
| 几何/期刊明确参数 | 优先保持，如直径500mm、50kW目标工况、beta主线及明确换热系数 |
| 完整mu(T)节点缺失 | 允许校准正值、温度升高黏度下降的曲线形状，保留参考点约束 |
| 绝对温标/初场缺失 | 允许用于初始化或重建的温标假设，不能强制运行中Tmax恒定 |
| 发射率、环境温度等借用值 | 来源不明确时允许校准，明确与博士/工程假设的区别 |
| 代码单位、面积、离散错误 | 先修正，不作为调参变量 |

若明确参数也需偏离以形成工程演示，单独保存 `ENGINEERING_CALIBRATED` 分支配置，不冒称全部采用期刊值。不要让过往“所有参数不得调整”的文字阻止当前已授权的缺失参数校准。

### D2．校准顺序与预算

先完成A/B，避免用mu/h补偿热边界实现错误。第一轮只选择对当前误差最相关的一类不确定输入：优先mu(T)冷区斜率或环境/发射率中的一项，依据温度与速度偏差决定。

每类先基准+低/高两个候选，不做全组合。参数范围有来源用来源；无来源记录为工程探索范围。曲线不能出现负值或不连续尖峰；不可通过降低mu到明显不合理来单纯放大速度。

最多两类参数的分阶段筛选、每类不超过3个候选；筛选用中等时窗，最终只续跑最好的一至两个。保存全部结果，不能只保留最像图的一次。

### D3．优化目标

主目标：内部主环流+靠冷壁外围反向环流、两者间上升通道、方向和相对范围。不是强制上下两个规则涡。

次目标：Fig.8已标0.5–0.9等温线的空间形状及中央交点。温度场是自洽演化产物，不能每步覆盖文献T后再声称热场复现。

不以速度大作为目标；文献无可靠速度标尺时只报告实际速度，不反推一个虚构参考值。

### D4．如何使用现有数字化

复用Doc204的原始等温线、方向及ROI，不重复提取。人为补全场的热点(0.85,0.70)、675K温差都不是作者原始精确数据。

用真实读图约束比较：各级等温线距离、主/次回流区域与方向。相关系数只作补充，不能掩盖等温线残差；ROI是搜索区域，不是精确涡心。

成果可标记 `CALIBRATED_CIRCULATION_TREND_MATCH`，表示经参数校准实现工程趋势对齐；与原始参数全面验证分开。当前目标无需先拿到作者完整数据。

## 7．任务E：后处理小修与可信度补充【并行且限时】

### E1．先修明确统计问题，不重跑

现 `wall_face_mass_flux_and_cv.py` 首层只输出带符号均值。增加体积加权mean(abs(un))、RMS、p95/max及原均值；防止正负抵消。

使用导出的mass计算控制体质量，并检查rho*Vol与mass在同一采样时刻是否同步。不同步先修导出/解释差异，不能用此差异直接判连续性失效。

截面流量延用已修正/2band，记录带宽。顶部真实面通量、首层速度、几何越界分开，不将箭头长度当穿墙证据。

### E2．先验证测量器，再评价求解器

用已知无散度、边界不穿透的解析圆柱速度场测试IDW/CV。可用轴对称流函数构造，无需跑流体：

设s=r/R、w=z/H，f=s²(1-s²)²，g=w²(1-w)²，psi=U R² f(s)g(w)，则

\[
u_r=-U(R/H)s(1-s^2)^2g'(w),
\qquad u_z=2U(1-s^2)(1-3s^2)g(w).
\]

其中g'=2w(1-w)(1-2w)，轴线上ur=0，uz=2Ug。转换成笛卡尔速度时显式处理r=0。所有物理面速度为0，内部有循环且无散度。

在与实际相似的点云采样上，用旧IDW重建面速度；如果仍报大泄漏，证明其不适于该验收。不得直接将面速度强制为0后宣称测量器通过；需测试不施加答案约束时的重建误差及网格趋势。

CV测试先稳态解析场、固定密度，后再处理真实多帧。粒子固定mass与速度连续性各自覆盖范围分开说明。

### E3．限时与退出

E线限定一轮脚本修正+一轮独立解析场测试，优先旧数据零仿真完成。纯测量误差未解决时保留 `DIAGNOSTIC_UNCERTAIN`，不阻塞A/B/C/D。

如果发现真实明显越界、局部密度持续异常或原生声学更新错误，才升级为水动力阻塞。不能因一个未验证的IDW指标失败重写整个求解器，也不能提前认定其全部是方法误差。

## 8．任务F：C1方腔的有限补账【非法国主线前置】

复用现有2D入口，严格Ra=1e5、Pr=0.71的参数生成见Doc204/206；先一个中级配置，必要时一个细级。保留旧表格输入作为另工况，不能将实际Ra=1.203e5与1e5参考直接验收。

报告两种Nu定义、两个中心线速度极值；旧水平速度偏差约15%、Nu约2倍分歧若仍存在，定位与本轮热边界/近壁梯度是否同源。

该任务单独排队，不重新启动2024全部算例。不因未完成就声称法国浮力完全不可用；如发现通用热/力算子重大错误，再界定实际受影响范围。

## 9．浮力阶段如何收口

本轮采用分层验收：

1. `BUOYANCY_ENGINEERING_USABLE`：已有基本公式、beta0、稳定回流、质量固定及时间/空间敏感性支持，冻结成功组合。
2. `WALL_CV_DIAGNOSTICS_VERIFIED`：测量方法经过独立检验才授予；未完成可以保留欠账。
3. `CALIBRATED_CIRCULATION_TREND_MATCH`：自洽热—流经记录的缺失参数校准后符合主次回流趋势。
4. `FULL_LITERATURE_VALIDATION`：本轮不默认要求/授予。

基础浮力工程可用可以收尾，不再因主次涡未完全吻合回到公式排查。热边界数值准确性和校准目标仍需分别报告。

## 10．软件结构与回归

- 全部改动隔离OPHELIE，不改全框架默认，不动A-phi分支。
- 旧成功小圆柱保持可运行，只在相关力学装配改变时短回归；纯热/后处理改动不全重跑水动力。
- 复用现有字段、GPU执行策略、邻域与热核，不复制多个法国主程序。
- 每个运行目录输出effective_config、source/input/binary哈希、实际物理时间、退出原因和完整日志。
- 仿真/后处理失败不能被 `|| true` 吞掉后产生PASS；缺参考不使用硬编码值替代。
- 不将“最新CSV”默认作为某工况结果，按运行ID、时间和配置匹配；初态和热历史相同才能作因果对照。

## 11．执行顺序与停止条件

1. 当天优先A解析参考与粗细热核对照；E脚本检查可同期进行，不要求多智能体或额外并行编程。
2. A失败才触发B；B最多两个有根据候选，每个必须回到解析问题验证。
3. 选定热边界后C接回固定Q；根据热流/温度是否仍演化决定续跑。
4. D有限校准直接服务甲方环流趋势；F独立排队。

若两轮局部修正无改善，交付温度/热流误差与具体失败项，不再展开十几个开关。可以用小CPU双精度热参考，不能让GPU性能重构拖延准确性。

本轮不继续：全局升功率抵消未知散热、盲目百万粒子、反复延长CV导出帧、重写浮力、同时扫描所有参数、把softcool偷偷改名成期刊原工况。

## 12．交付清单

- `ROUND207_DECISION.md`：本轮唯一主误差、修复与下一动作。
- `ROBIN_INDEPENDENT_REFERENCE_REPORT.md`：解析式、实际物性/时间、温度/热流/累计能量及收敛。
- `THERMAL_WALL_IMPLEMENTATION.md`：旧新自由度位置、面积、阻力和能量说明。
- `FRENCH_FIXED_Q_COMPARE.md`：同初态旧新热边界对照及主次回流。
- `PARAMETER_PROVENANCE_AND_CALIBRATION.csv`与全部候选结果。
- `FIG7_FIG8_COMPARISON.md`：真实读图约束与重建假设分开。
- `HYDRO_DIAGNOSTIC_STATUS.md`：绝对/RMS法向速度、解析场测试、CV限制。
- `BUOYANCY_ENGINEERING_CLOSEOUT.md`：冻结组合、适用范围、哪些事不再重复。
- `VALIDATION_DEBT.md`：C1/高Pr/EM精度等剩余项及是否影响当前目标。
- 原始关键CSV/空间快照、脚本、配置、源码diff和复现命令。

**最终应让用户明确知道：浮力计算已能稳定用于工程推进；近壁热传递是否已被独立验证；哪些参数经校准；法国主次环流形态还差在哪里。不要再用一个含糊的总PASS或总FAIL概括所有模块。**


<!-- ===== END 207_OPHELIE_ROBIN_VALIDATION_CALIBRATION_AND_BUOYANCY_CLOSE_CURSOR_PLAN_2026-11-02.md ===== -->

---

<!-- ===== BEGIN 208_OPHELIE_THERMAL_FIELD_AND_BUOYANCY_CLOSE_CURSOR_PLAN_2026-11-03.md ===== -->

<a id="208-ophelie-thermal-field-and-buoyancy-close-cursor-plan-2026-11-03-md"></a>

# 208_OPHELIE_THERMAL_FIELD_AND_BUOYANCY_CLOSE_CURSOR_PLAN_2026-11-03.md

# Doc208｜OPHELIE 浮力阶段收口、热边界修正与法国自然对流推进任务书

日期：2026-11-03  
适用对象：Cursor 开发执行、结果审查、甲方工程演示准备  
依据：`DOC207_ROUND_PACK_2026-11-03.tar.gz` 的代码快照、平板对照、法国固定 Q 原始 CSV/summary；Doc204 图7/8数字化资料；法国自然对流期刊原文。

> 本轮目标：保留已经能够产生回流的水动力装配，修正已查明的监测问题，补齐局部热边界证据，推进接近文献的非均匀温度场。不要再次把全部工作退回浮力公式、黎曼求解器或框架替换。
>
> 本文中的误差门槛是本项目工程验收建议，不是法国论文或 SPH 社区的统一标准。未达到时应报告原因和局限，不得通过删除早期数据、修改标签或改色标获得“通过”。

## 0. 本轮决策与范围

### 0.1 已接受的工程状态

- 浮力模块按 `BUOYANCY_ENGINEERING_USABLE` 阶段收口：已有 β=0 对照、质量守恒、固定温度驱动和法国接回的回流证据。
- 此标签表示可以继续做热—流耦合，不等于 C1 标准方腔、近壁连续性、所有流场精度均已完成验证。
- 现有 WCSPH／Riemann 装配保持。没有新的独立反证前，不新增 ISPH、δ-SPH 或全局压力离散重构任务。
- A-phi／Laplace 自研方法分支不在本轮修改范围；所有修改限制于 OPHELIE 局部模块、对应案例和后处理。
- 优先 GPU；若已有 CPU 隐式黏性或参考程序能可靠解决刚性，可作为局部备选。不得为追求全 GPU 重新开启无期限求解器开发。

### 0.2 必须纠正的判断

1. Doc207 已证明旧壳温 Robin 在平板粗网格上存在系统性过散热偏差，但没有证明法国全部热场差异都来自这一项。
2. 新 wall-temperature 方案明显改善误差；7.5 mm 下早期热流与累计散热仍未全部达标，不能写成严格 PASS。
3. 法国工况实际施加新壁温通量，但监测输出仍调用旧公式；先修这个明确错误，再讨论散热量和参数选择。
4. 从均匀 1473 K 初态开始的冷却，不能直接等同于错误。文献比较对象是耦合收敛后的非均匀温度场。
5. ε=0.4 只证明更暖、同一时间窗内回流方向相近，尚无依据认定它是“最佳文献匹配参数”。
6. `CLOSED_RETURN_OBSERVED` 目前主要是上下流量和方向证据，不能自动升级为两条闭合流线或两个独立涡已确认。

### 0.3 优先级与执行关系

| 优先级 | 任务 | 依赖 | 本轮完成标准 |
|---|---|---|---|
| P0 | A：统一新旧热边界的监测、账本、summary | 无 | 同一状态下通量口径一致；已有报告数值更正 |
| P0 | B：新壁温平板 3.75 mm 补证 | 可与 A 的整理工作独立进行 | 全时刻误差表、收敛判断和明确结论 |
| P0 | C：法国边界面积与新通量短对照 | A | 面积误差受控；保留回流；热账本一致 |
| P1 | D：二维轴对称热初猜及映射检查 | A；正式接回需 C | 一份不固定温度的热场接回结果 |
| P1 | E：图7/8直接比较及稳定性观察 | D；C 可先试用后处理 | 同坐标、同归一化的图与误差表 |
| P2 | F：缺失参数有限校准 | D/E 暴露明确误差后 | 只调有依据的少量参数，明确目标 |
| P2 | G：C1 与测量器独立欠账 | 不阻塞 A–F | 记清状态，不冒充已完成 |

执行是分阶段的，不要求一口气完成全部长时间算例。每阶段先提供可审查的结果；小算例失败只回退对应局部问题。

## 1. Doc207 冻结事实与参照数据

### 1.1 独立平板：t=50 s，相对半无限 Robin 解析解

| 模型 | dp | 瞬时 q 误差 | 累计 E 误差 | 近壁归一 RMS |
|---|---:|---:|---:|---:|
| legacy-shell | 15 mm | +32.8% | +27.3% | 约3.9% |
| wall-temperature | 15 mm | −6.1% | −14.3% | 约1.0% |
| wall-temperature | 7.5 mm | −1.4% | −4.5% | 约0.55% |

7.5 mm 新模型在 t=5/10 s 的 q 误差约 −9.0%/−5.8%，不能只引用 t=50 s 结果。

### 1.2 法国固定 Q：目前可以确认的回流方向

- 轴心中部下降；
- r/R=0.7–0.9 环带上升；
- r/R=0.9–1.0 最外侧下降。

这与文献主回流和外围回流的机制相容，但还需要直接矢量图、流线与分隔位置证据。

### 1.3 温度场主要差距

Doc207 t=20 s 三组 ε 的粒子温度范围均约为 1350–1488 K，Tmin/Tmax≈0.907。法国图8包含 T/Tmax=0.5–0.9 的等温线，低温区深入内部，并非只在物理壁面上。

注意：粒子温度范围不包含所有重建面温；仍须区分 Tw 与 Ti。但补出低温 Tw 不能替代内部温度场的比较。

### 1.4 当前散热数据必须重新解释

完整 0–20 s 的实际累计账本平均散热为：

| ε | E_boundary/t | E_side/t | 输入功率 |
|---:|---:|---:|---:|
| 0.9 | 119.8 kW | 77.3 kW | 50 kW |
| 0.6 | 110.7 kW | 77.3 kW | 50 kW |
| 0.4 | 104.1 kW | 77.3 kW | 50 kW |

旧报告 136.7/124.8/116.5 kW 来自混入旧公式的监测列，而且取样时间窗与上表稍有不同。不能把两组数当作同一模型的独立一致性证据。冷却瞬态判断仍成立。

## 2. 开发前冻结与代码入口

### 2.1 先建立可追溯基线

保存 Doc207 源码、编译参数、命令行、输入 Q 文件、参数表及其哈希。保持旧结果只读；新结果进入 `doc208_runs/`。不得覆盖 Doc207 数据并沿用旧标签。

不要预设服务器仓库绝对路径。实际定位项目根目录后，记录 `git rev-parse HEAD`、branch、相关修改的 diff 和二进制哈希。未提交源码也必须归档，单独 HEAD 不足以复现。

### 2.2 本轮实际生产执行入口

Doc207 运行脚本实际执行：

`tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_reduced_pressure/test_3d_ophelie_french_natural_reduced_pressure.cpp`

不要只因为报告写着 frozen-Q，就修改另一个同名相似的旧案例：

`test_3d_ophelie_french_natural_convection_frozen_q`

先以运行日志、summary 的 argv 和实际二进制为准核实入口。

局部热模块：

- `extra_src/shared/electromagnetic_ophelie/diagnostics/electromagnetic_ophelie_thermal_diffusion_one_way.h`
- `extra_src/shared/electromagnetic_ophelie/diagnostics/electromagnetic_ophelie_joule_to_heat_one_way.h`

小平板入口：`test_3d_ophelie_robin_plate_halfspace`。

### 2.3 固定的主线设定

| 项目 | 本轮基线 | 说明 |
|---|---|---|
| R、H | 0.25 m、0.185 m | 自然对流几何；不混入搅拌几何 |
| 竖直轴 | 当前自然对流代码 z | 显示可映射；不套用搅拌 STL 的 Y 坐标 |
| β | 1e−5 K−1 | 保留 |
| ρ | 2750 kg/m³ | 保留 |
| Tref | 1473 K | 物性锚点／浮力参考；不强制等于全域初温或稳态平均温度 |
| c0 | 20 m/s | 保留 |
| 玻璃吸收总功率 | 50 kW | 不因冷却瞬态提高到约100 kW |
| Q 的空间形状 | 现有固定 Q_rz 快照 | 明确是数值快照，尚非文献最终热源已验证 |
| live EM | 关闭 | 热—流阶段先分离定位 |
| 顶面流动边界 | 现有平顶滑移装配 | 不新增 Marangoni |
| side/base 流动边界 | 现有闭壁与无滑移装配 | 不重写 |
| 温度物性 | 现有可追溯 k(T)、cp(T)、μ(T) | 报告借用／拟合来源及范围 |
| wall-temperature | 新工程配置显式启用 | 旧模型保留为回归对照 |

稳定配置一经确认，应通过具名配置和启动脚本固定。保留旧默认不意味着新工程运行仍可误用旧模型：启动时打印、summary 记录，并校验选项确实生效。

## 3. 任务 A：修复监测与实际热边界不一致（P0）

### 3.1 已确认错误位置

在 reduced-pressure 案例 `write_monitor_row()` 内，目前 `bc_active` 分支无条件调用：

`hostOphelieThermalFrenchNaturalHeatLossPowers(...)`

但热更新路径在 `use_wall_T` 时使用新 wall-temperature 算子，并使用新通量累计 `ledger.e_boundary`。监测函数随后覆盖 `last_P_*`，进一步污染最终 summary 的瞬时功率。

### 3.2 必须修改

- 按实际模型选择监测通量。`use_wall_T` 时调用 `hostOphelieThermalWallTemperatureHeatLossPowers(...)`，否则才使用旧函数。
- 最好封装一个 OPHELIE 局部调度入口，供写监测、初始损失、最终 summary、分面审计共同调用，避免多处重复 if 后再次漏改。
- 面积、gap、k(T)、温度取值和采样相位必须与施加模型一致。
- 监测函数必须只读，不推进温度、不重新初始化字段、不修改累计账本。
- 为每条记录标注模型和采样相位。若源／导热／BC 以分步方式执行，说明功率是在 BC 前还是步后取值。
- host 分面功率函数目前用 `Tw_rep` 拆分顶面对流与辐射；当前独占主面时需检查无误。若未来支持角点同时多面，必须使用各面的 Tw，不能用混合平均 Tw 计算辐射。

### 3.3 增加三种互相独立的功率核对

1. 同一状态：实际模型 reduce 的总功率，与 host 按相同物理公式计算的总功率比较。
2. 分面和：Pside+Pbottom+Ptop_conv+Ptop_rad 与 Ptotal 比较。
3. 累计量：Eside+Ebottom+Etop 与 Eboundary 比较；同时独立核对 ΔH=Ein−Eout。

建议双精度 host 核对相对误差 ≤1e−4；实际 float/GPU 规约若超出，应给出精度分析，不能直接放宽到百分数级。

短跑能量账本建议保持既有 ≤1% 输入能量门槛；当前数据约1e−6量级，不应修监测后明显退化。

### 3.4 时间平均必须重写清楚

窗口 [t0,t1] 的平均实际损失优先取：

`P_loss_avg = (E_boundary(t1) − E_boundary(t0)) / (t1 − t0)`

分面同理。监测功率积分作为另一条检查，必须按真实时间加权；不要对不等间隔记录简单求算术平均。

“散热随时间变化”不能解释正确积分后的能量不一致。若粗采样积分不准，应加密采样或明确积分误差，不得把33%残差解释成正常漂移。

### 3.5 最小运行与交付

- 新、旧模型各跑一个相同初态的 1–5 s 短检查；不重跑原三组20 s矩阵。
- 在新模型同一时刻保存 Ti、Tw、gap、k、权重、分面 q、总 P 的审计样本。
- 只用原 Doc207 累计列重新计算已有三组平均值，标记 `REANALYSIS_OF_EXISTING_DATA`。
- 原有瞬时新模型功率如果没有保存足够字段，不伪造修正版历史；标记不可恢复，必要时仅重跑关键基线。

交付：`A_MONITOR_FIX_REPORT.md`、最小 diff、`POWER_AUDIT.json`、修正后的 `DOC207_REANALYSIS.csv`。

## 4. 任务 B：新壁温处理的最后一轮小算例补证（P0）

### 4.1 保持已定义的独立问题

半无限 x≥0，初始 Ti=1473 K，环境 T∞=303 K，h=300 W/(m² K)，ρ=2750、cp=1150、k=4。无流动、无热源、无辐射；远端与横向绝热，L=0.12 m。

使用生产相同的导热与 wall-temperature 算子，不另造只为测试通过的边界更新。

令 α=k/(ρcp)，η=x/(2√(αt))，b=h√(αt)/k：

`θ=(T−T∞)/(Ti−T∞)=erf(η)+exp(−η²) erfcx(η+b)`

`Tw=T∞+(Ti−T∞) erfcx(b)`，`q=h(Tw−T∞)`。

检查 Tw(10 s)≈1189.945 K、Tw(50 s)≈969.202 K。

Doc207 解析脚本实际采用的 `exp(2ηb+b²) erfc(η+b)` 与上述公式等价，数值结果可保留。清理其冗长且互相矛盾的注释，不要误以为需要把正确结果改成 `exp(+η²)erfcx(...)`。

### 4.2 本轮只新增一个主要分辨率

新增 wall-temperature，dp=3.75 mm，输出 t=0、1、5、10、20、50 s。复用已有15/7.5 mm结果形成收敛表，不新增全局三维法国细网格任务。

若3D平板规模不合适，可采用经核对的二维同物理平板；但维度、测度和生产内核差异必须说明，不能把二维新结果与三维旧结果混称同一网格序列。

### 4.3 必须检查的量

- q(t) 相对解析误差；Eout(t) 相对解析误差；粒子位置上的 T(x_i,t) 误差；Tw_recon 误差。
- Ti 与 Tw 不能混比。温度归一误差保留以1170 K为分母的定义，同时给出 K 单位误差，避免大归一分母掩盖壁温偏差。
- t=1 s 与 t=5 s 的欠散热不得删除。累计 E 不能从一个更晚时刻重新清零来掩盖早期误差。
- 检查时间步减半敏感性：仅在误差不随空间加密改善、或早期时间步可疑时触发一次。

### 4.4 判定与停止条件

沿用原目标：t≥5 s 的 |q误差|≤5%；各报告时刻 |E误差|≤3%；近壁归一 RMS≤3%。逐时刻列出，不以 t=50 单点替代全时段。

- 若全部满足：标记该分辨率的 `ROBIN_REFERENCE_PASS`；另列15 mm生产误差，不能跨分辨率授予精度。
- 若多数改善但早期仍未过：标记 `CONVERGENT_WITH_EARLY_TIME_DEBT`，保留明确误差范围；法国热初猜与形态工作可以继续，但不得声称热边界全面验证。
- 若不收敛或误差明显反向增大：只排查 gap、面积分配、第一层温度代表位置和导热边界截断，不回头改浮力。
- 本轮最多增加一次局部边界候选或一次必要时间步核对。不得形成无限加密或连续更换热模型。

## 5. 任务 C：法国工况面积和壁温短对照（P0）

### 5.1 已知问题

Doc207 法国工况 `area_normalize=0`：侧壁面积偏小约10.34%，顶面约10.50%，底面亦有偏差。平板面积正确不能替代圆柱面积检查。

不能把约10%的离散面积缺口留在基线里，再通过 ε 或 h 校准去补偿。

### 5.2 最小对照矩阵

保持同一均匀1473 K初态、dp=15 mm、固定Q=50 kW、h_side=300、h_bottom=35、h_top=20、ε=0.9、环境303 K、现有变物性和水动力装配：

| ID | 热壁面模型 | 面积归一化 | 用途 |
|---|---|---|---|
| C0 | wall-temperature | 关闭 | 修监测后的可复现旧配置 |
| C1 | wall-temperature | 开启 | 下一阶段具名工程基线 |

先 t=1–5 s；稳定后 C1 到20 s。C0仅在修监测短跑不能与Doc207对接时补20 s；不要默认重复全部矩阵。

ε=0.4 保留已有敏感性资料，本阶段不再扫。这里的 C1 是本任务局部编号，不是“标准方腔C1”；输出目录用完整名字，避免混淆。

### 5.3 面积归一化审计

- 分面目标：Aside=2πRH；Atop=Abottom=πR²。
- 面权重实际进入算子；不是只修 summary 的面积数值。
- 核对核内使用的 Vol 与构造权重时使用的 mass/ρ0 是否一致。WCSPH 若导致微小差别，报告实际施加面积，不能只报名义面积。
- 分面面积相对误差建议 ≤0.5%。角点当前采用最近主面独占分配，先保持，不顺手改成多面分配。
- gap 使用所属物理面的真实几何距离，不能使用2dp、热扩散长度或用于拟合的常数。
- 粒子排序后 gap、面权重、k、cp、ΔT 的对应关系保持正确；若会重标记，记录频率和重标记前后误差。

### 5.4 验收

面积和功率审计通过、质量不重置、无明显逃逸／失稳，保留合理的回流方向即可完成本阶段。不要求20 s达到热稳态，也不要求为了平均温度不下降而调整已知热参数。

新配置具名保存，例如 `french_fixedQ_wallT_areaNormalized_reference`。启动脚本显式写 `--thermal-wall-model=wall-temperature --area-normalize`；启动失败、字段缺失或审计失败必须非零退出，禁止 `|| true` 吞掉核心验收失败。

## 6. 任务 D：建立非均匀热初猜，释放温度后接回（P1，主线）

### 6.1 为什么需要这一阶段

文献的图7/8来自收敛的二维轴对称耦合。当前均匀1473 K初态只演化20 s，内部仍很接近初温，两者不是同一热状态。

约1500 K是论文对浴温的概述，不能当作平均温度必须保持1500 K的验收点。也不能用 Tmin 降低作为唯一优化目标：必须同时比较温度分布形态、热收支和回流。

### 6.2 首选：轻量轴对称导热参考生成初猜

在独立 Python 工具或独立参考案例中，解圆柱轴对称热问题：

`ρ cp(T) ∂T/∂t = (1/r)∂[r k(T) ∂T/∂r]/∂r + ∂[k(T)∂T/∂z]/∂z + Q(r,z)`

- 先不含对流，只服务热初猜和无流动热参考；不能称为法国最终解。
- 推荐守恒有限体积形式，使用圆环单元体积与各面面积，轴心 r=0 自然零径向通量，避免1/r奇点。
- 侧、底和顶面的边界使用与C阶段同一连续物理条件；顶面对流加绝对温度四次方辐射。
- 参考程序应在物理面施加通量，不复刻SPH壳层偏差。串联距离仅代表参考网格中心到面的数值半格，不附加一层未经定义的物理热阻。
- 使用现有固定Q快照，但在轴对称单元上重新积分归一到50 kW：积分包含2πr，不可将二维平面面积当体积。
- Q的网格、原始来源、覆盖区域、映射和缩放必须单独记录。该快照可能带有早期温度／电导率状态，不能称为文献收敛Fig6热源。
- 可用隐式时间推进或直接稳态迭代；若用伪时间加速，输出中明确没有真实瞬态时间意义。
- 至少两级廉价参考网格，建议径向/轴向约64×48、128×96，并在壁附近适当聚集。不能仅凭网格数认定精度，需检查边界通量与温度剖面的变化。
- 参考解总体热平衡建议 |Qin−Qout|/Qin≤1%，两级主要温度剖面差异与通量差异需报告。达不到则只作未收敛初猜，不冒充独立稳态参考。

无对流的热解与法国真实热—流状态可能差异明显。这是预期的用途限制，不是通过强制拟合文献来修饰的理由。

### 6.3 次选：文献数字化构造带不确定性的初猜

若轻量热工具预计开发成本过大，可复用 Doc204 数字化数据构造初猜，但必须遵守：

- 图8给的是T/Tmax，不是绝对K温度。Tmax单独标注为假设／校准量，不伪称作者给定值。
- 0.9等温线内热核未完全约束、凸包外区域无数据，必须标记外推区域和构造方法。
- 可以先使用物理热参考或平滑受约束延拓填缺，不直接把NaN当0 K、不使用任意热斑当作者数据。
- 该初猜来自文献，再用同一图在t=0评分属于构造一致性，不是独立验证。必须在释放温度后报告漂移，最好另用未参与构造的剖面或实验数据评价。

两个初猜路径选择一个先完成，不要求本轮同时开发两套。

### 6.4 映射到SPH粒子前的必要检查

1. 记录参考坐标 r/R、z/H 与自然对流代码坐标映射，明确 z_bottom 和几何中心；禁止上下翻转。
2. 采用稳定有界的插值；不产生负温度、局部尖峰或未报告外推。
3. 输出源网格与粒子场的 Tmin/Tmax/Tmean、总焓和等温线叠图。
4. 验证 Q 积分仍50 kW、质量不改变；不为匹配温度而改粒子质量或密度。
5. 使用程序实际热参考常数初始化 ΔT，保证 `Temperature = thermal_reference + DeltaT`；不得只写Temperature而留旧ΔT。热参考常数与Boussinesq Tref若分开，分别记录。
6. 依据映射后温度更新μ、k、cp及其他当前启用物性，排序后再次核对字段归属。
7. 压力／背景状态使用已验证装配的一致初始化，不能直接套用均匀温度初态的特殊平衡设置；不得每步重写密度来压制启动误差。

### 6.5 高黏度与启动保护

低温区形成后μ可能显著上升，现有热启动下μ_max只有约12.6 Pa·s，尚未覆盖文献高黏度冷区。

- 先输出预估 μ_min/max、各温区体积分数、黏性时间步限制和预计耗时；检查现有μ曲线有效范围及cap。
- 不得静默减小μ或降低cap来获得更强回流。文献以高黏度区表征近似固化，不要求本轮新增完整相变／FSI模型，但必须记录当前是否有对应冷区。
- 若GPU显式黏性变慢，优先调用已有可靠隐式／CPU局部备选，并做短时间一致性对照；不要为这一初猜任务重写GPU隐式求解器。
- 若初猜包含大体积低于材料曲线范围的温度，先检查初猜与物性适用性，不直接长跑。

### 6.6 最小接回运行

- D0：新初猜，β=0，短时间1–5 s；温度可演化，检查噪声、质量、边界和热账本。
- D1：同初猜，β=1e−5，其他与D0相同；先5 s，健康后到20 s。
- 初始化后解除固定温度：Q、导热、边界和物性正常更新。记录没有每步回写文献温度或初猜。
- 若回流已稳定、温度仍有持续漂移，则保留该事实，再决定分段延长或改进初猜；不能仅因到20 s还没稳态就再改浮力。
- β=0测试不要求完全无温度变化；它用于区分热演化与水动力伪流。

交付：初猜工具与输入、收敛历史、映射审计、D0/D1原始场、温度释放证明、耗时和步长统计。

## 7. 任务 E：直接比较图7/8，而不是只看标签（P1）

### 7.1 使用已有数字化资产

先定位 Doc204 数据包 `204_FIG7_FIG8_DATA_AND_CURSOR_PLAN_2026-11-01.zip`，其中：

- `reference/fig8_labelled_contour_points.csv`
- `reference/fig8_midline_contour_intersections.csv`
- `reference/fig7_direction_samples.csv`
- `reference/fig7_topology_regions.csv`
- `reference/digitization_provenance.json`
- 原图及qa叠图。

不要再次声称“论文图抓不出来”，也不要重新数字化已完成的轮廓。先检查坐标、归一化及不确定性；无法找到资产时明确列出缺失文件，保留后处理接口，不编造参考数据。

### 7.2 Python最小图件

每个关键状态至少输出：

1. 轴对称平均T/Tmax云图，与论文0.5/0.6/0.7/0.8/0.9轮廓叠加。
2. 子午面ur、uz等长度箭头，颜色显示速度大小；另给按幅值缩放图，避免微弱箭头全部看不见。
3. 使用直接速度插值的流线图，清楚标注轴对称平均、插值方法、无数据区域和物理边界。
4. z/H=0.25、0.5、0.75处T与uz径向剖面；图3对应的半高温度剖面可另作参考，不混用其不同h工况。
5. t=0、5、20 s及后续分段末态的同色标对照。

速度变换：`ur=(x ux+y uy)/r`，使用实际轴心坐标；r≈0单独处理。先在r-z分箱作体积加权平均，输出每箱样本数；不能把全三维粒子直接散点插值后误称轴对称平均。

温度使用 `T/Tmax`，其中T为绝对K温度，Tmax来源与时刻清楚。严禁改用 `(T−Tmin)/(Tmax−Tmin)` 获得看起来相似的轮廓。若使用轴对称平均场最大值，须同时报告原始粒子最大值和差异，避免噪声最大值控制全部归一化。

### 7.3 温度定量误差

在数字化等温线点上插值模拟T/Tmax，按每个标签输出：有效覆盖率、误差中位数、P95、带符号偏差。

- 不在无数据处默认为0，不只挑匹配较好的点。
- 数字化曲线上的密集像素不是独立实验样本。统计同时按等温线汇总，不能靠某条长曲线的点数支配总分。
- 给出坐标读图不确定性及±一个局部粒子间距的敏感性；粗网格下不要宣称毫米级涡中心精度。
- 若目标等温线在模拟场中不存在，直接标为 `CONTOUR_ABSENT`，不要自动改等值线。

### 7.4 环流形态验收

分别列出：

- 轴心下降是否存在；内部上升带的位置、宽度和峰值；最外侧下降带是否存在。
- 主回流与外围反向回流是否各自出现闭合流线；是否有明确分隔／停滞区域；跨至少两个输出时刻是否稳定。
- 外围涡是否依赖一层粒子或插值边界补零；对合理插值带宽变化是否消失。
- 图7人工方向样本只有方向和较大读图不确定性，没有可用于速度误差的m/s标尺。比较角度／方向一致性时剔除近零速度不确定点并报告覆盖率。

禁止仅根据“轴下降、环带上升、外层下降”三个符号授予“双涡已验证”。也不要求把文献两涡硬解释为上下两圆涡；其描述为内侧主胞与靠冷壁外围胞。

### 7.5 流函数和连续性诊断的使用限制

- 直接速度箭头与直接速度积分流线为首要形态证据，流函数作为辅助。
- 轴对称不可压关系可用 `ur=−(1/r)∂ψ/∂z`、`uz=(1/r)∂ψ/∂r`；但若两路径重建ψ明显不一致，报告测量／连续性残差，不能只选好看的路径。
- 未通过解析场验证的IDW壁通量和CV残差不能单独判定求解器漏流；同样不能用 signed-mean un 接近0证明不穿透。
- 壁面报告至少包括有符号均值、|un|均值、RMS、P95及粒子越界距离。固定零壁值可作带物理约束的可视化方案，但不能反过来证明求解器满足壁条件。

### 7.6 是否延长运行的依据

不从U≈0.1 mm/s直接断言必须跑多久，也不承诺20 s足够。按真实输出检查：

- 末段输入输出热量差和ΔH变化；
- T(r,z)、Urms、上升带位置、主次涡形态在相邻时间窗的变化；
- 是否仍受初猜显著控制。

建议分段5→20→50 s，之后按耗时和变化趋势决定是否到100 s或更久，不默认启动一天长跑。若从近稳态初猜启动仍迅速远离目标，优先检查Q形状／边界／物性，而不是仅加长时间。

热稳态工程候选可用：末窗 |平均Qin−平均Qout|/Qin≤5%，且主要温度和流型在连续两个窗口变化较小。形态变化建议≤5%作为观察阈值，窗口长度须与实际演化相符，不能用极短窗口制造稳定。

## 8. 任务 F：有目标的有限参数校准（P2）

### 8.1 校准前提

A修正完成、面积处理明确、参考与模拟归一化一致、至少有一个温度释放后的非均匀热场结果，才开启新扫描。

用户已允许对文献未给清楚的参数进行工程校准，不必继续等待作者数据。已给定的几何、50 kW、β及明确边界设定优先保持；若使用论文不同h工况，必须以独立工况命名，不伪称同一case。

### 8.2 参数角色

| 参数 | 可以解决的主要问题 | 不应期待解决的问题 |
|---|---|---|
| ε、辐射环境温度 | 顶部散热、表面冷区 | 不能自动修正侧壁离散偏差 |
| μ(T)冷区曲线 | 冷区流动抑制、回流分布、速度幅值 | 不能替代正确热源和热边界 |
| 初猜绝对温标 | 缩短启动、检验初态敏感性 | 不能固定回写作为最终热解 |
| Q的空间形状 | 热核位置、内部升流位置 | 不能在未披露情况下冒充已验证EM |
| h_side | 不同文献工况的侧壁传热差异 | 不能用来补偿已知面积／gap错误 |

文献列有h_side=100、200、300、400 W/(m² K)工况，但图7/8目标的具体工况来源仍需标明。不能把代码300直接称为图7/8唯一确定值，也不随意选最低h只为减小平均散热。

### 8.3 扫描预算与评分

- 每次只改变一类参数，最多“基线＋两个候选”；已有ε矩阵可复用，不重复做同样扫描。
- 本轮最多新增一类参数探索。优先由E的实际误差决定，而不是预设继续调ε。
- 评分至少包含图8温度轮廓误差、图7方向／流型证据和热量账本；已知物理值偏离程度单列。
- 不以Tmean最接近1473 K、温度最高、速度最大作为唯一目标。
- 最佳参数必须能在温度释放后维持或趋向相近状态。若仅t=0相似、随后漂移，标记初猜匹配，不标记稳态匹配。

未达到完全定量吻合不阻止工程演示；但交付必须同时给出误差、未稳态状态和参数来源。

## 9. 电磁接回的边界（本轮不扩散）

当前固定Q模式意味着：使用既有电磁计算得到的空间分布快照作为热源，并非均匀加热，也并非实时电磁—温度双向耦合。

本轮核对Q支持域、映射和50 kW积分，不新开TEAM7边界攻关。固定Q下热—流形态达到可解释状态后，下一阶段再接σ(T)和live EM，比较Q形状变化及热核迁移。

届时必须做“固定Q vs live EM”对照；不能把固定Q得到的环流当作电磁反馈已验证。也不将文献约20次耦合迭代解释为20 s。

## 10. 独立欠账的有限处理（G，非主线阻塞）

### 10.1 标准方腔

严格Ra、Pr、Nu、速度剖面和边界梯度的一致性仍未全部通过。旧SMOKE_PASS不能改写为基准定量验证完成。

若资源允许，按原独立计划补一组严格物性与中／细网格；若不能完成，完整保留欠账和对结论的限制。不要用法国图像匹配来替代标准基准。

### 10.2 壁面／CV测量器

按Doc207解析无散场先验证测量器，再分析生产数据。不要重复数轮仅换插值半径、没有解析参考的扫描。

若只剩后处理欠账且未发现实际异常，继续热—流主线；若发现显著真实壁通量、密度异常或大范围非物理散度，则提交最小复现并局部处理，不能以“非阻塞”忽略实质数值错误。

## 11. 统一收口与交付标签

| 标签 | 允许授予的证据 | 当前／本轮目标 |
|---|---|---|
| BUOYANCY_ENGINEERING_USABLE | 已有独立浮力、对照和回流工程证据 | 保持 |
| MONITOR_MATCHES_APPLIED_BC | 新旧路径监测与实际通量核对 | 本轮必须完成 |
| ROBIN_REFERENCE_PASS | 指定网格、指定时刻通过独立解析门槛 | B按结果判定，不预授 |
| FRENCH_THERMAL_BOUNDARY_ENGINEERING_READY | 面积、通量、映射受控且公开生产误差 | 本轮争取 |
| CIRCULATION_DIRECTION_OBSERVED | 轴下降／内部升流／外侧下降 | 已有，直接图件复核 |
| TWO_CELL_TOPOLOGY_SUPPORTED | 两胞及分隔的直接证据、插值敏感性检查 | 本轮争取 |
| CALIBRATED_CIRCULATION_TREND_MATCH | 与图7方向／形态直接比较＋校准来源公开 | 通过E/F后再授予 |
| THERMAL_FIELD_TREND_MATCH | 图8同归一化轮廓及温度释放后的证据 | 本轮争取 |
| FULL_LITERATURE_VALIDATION | 热场、流场、来源、收敛及耦合均有充分证据 | 不作为本轮强制目标，不提前授予 |

浮力阶段收口后保留回归短测，不再把每次法国热场不匹配都重新归咎于浮力公式。

## 12. 建议实施批次

### 批次1：小修改与现有数据重分析

完成A；生成纠正表；冻结新基线；修正Doc207过度标签。首先交付这部分，不等待大算例。

### 批次2：小算例与圆柱短跑

完成B新增一级、C面积短对照；明确热边界剩余误差。不扩大到全框架修改。

### 批次3：热初猜接回与直接图像

完成一条D初猜路线，运行D0/D1，执行E。先回答“温度场能否进入论文相近分布、流型是否随之改善”，再决定新增参数扫描或延长时间。

### 批次4：少量校准与阶段交付

只在E暴露具体形态差距时执行F。提供可复现的具名工程配置、甲方可用图件及明确限制。E2/C1独立欠账不自动扩展本轮范围。

## 13. 运行与打包规则

建议目录：

```text
docs/ophelie/reference/doc208_runs/
  00_entry/
  A_monitor_fix/
  B_robin_finer/
  C_area_wallT/
  D_thermal_initial_guess/
  E_fig7_fig8_compare/
  F_calibration/
  status/
  scripts/
```

每个运行至少提供：

- 完整argv、环境／设备、git与二进制标识、源码diff、输入文件哈希；
- 模型、坐标、面积处理、物性来源、是否固定温度、是否live EM；
- summary、连续时间序列、关键粒子场；
- 有独立参考的误差表、Python原图及生成脚本；
- wall time、物理时长、步数、限制时间步的算子；
- 失败信息及停止原因，不把缺失字段写成0或true。

最终打包为 `DOC208_ROUND_PACK_<实际日期>.tar.gz`。需包含所用Q和数字化参考的必要文件，或明确随包依赖及哈希，不依赖开发机绝对路径和未打包旧目录。

`00_entry/ROUND208_DECISION.md` 必须逐条回答：

1. 监测和实际施加热边界是否一致？旧报告哪些值已更正？
2. 新壁温在哪个分辨率、哪些时刻通过？15 mm生产误差多大？
3. 圆柱实际施加面积是否正确？是否出现新的真实边界异常？
4. 非均匀初猜如何生成、是否含文献拟合、释放温度后怎样变化？
5. 图7主次涡证据是什么，图8轮廓差距多大？
6. 当前是启动瞬态、工程趋势匹配，还是已达到可解释的稳定状态？
7. 浮力阶段可否维持工程收口？下一阶段最需要推进的是热场、Q形状还是EM反馈？

## 14. 本轮禁止的重复工作

- 不重复全面排查已验证的Boussinesq公式。
- 不直接切换整个流体框架来追逐两条好看的流线。
- 不因Pout>50 kW就认定实现错误、提高热源或把所有冷却调弱。
- 不把“温度更暖”当作文献匹配更好。
- 不以粗网格真实粒子场没有低温等温线为由，用色标重映射制造轮廓。
- 不把文献构造初猜在t=0的匹配当作预测验证。
- 不只展示t=50 s尾点掩盖早期热流和累计能量误差。
- 不让新旧通量、不同面积、不同时间窗的数据再次混入一个对比表。
- 不静默修改物性曲线、质量、热功率、温度参考或边界参数。
- 不让未验证的流函数／IDW诊断替代原始速度证据。

**本轮成功的核心：修复明确错误、把热边界误差说清楚，并用释放演化后的非均匀温度场检验法国主次回流；让后续判断建立在直接数据上，而不是继续重复“浮力能不能算”的讨论。**


<!-- ===== END 208_OPHELIE_THERMAL_FIELD_AND_BUOYANCY_CLOSE_CURSOR_PLAN_2026-11-03.md ===== -->

---

<!-- ===== BEGIN 209_OPHELIE_CLEAN_RIEMANN_CAVITY_ROOT_CAUSE_AND_CLOSURE_CURSOR_PLAN_2026-11-03.md ===== -->

<a id="209-ophelie-clean-riemann-cavity-root-cause-and-closure-cursor-plan-2026-11-03-md"></a>

# 209_OPHELIE_CLEAN_RIEMANN_CAVITY_ROOT_CAUSE_AND_CLOSURE_CURSOR_PLAN_2026-11-03.md

# Doc209｜干净低耗散黎曼方腔对照、根因定位与定量验收

日期：2026-11-03  
修订：R1（补充法国环流可信度边界、实际算子映射和分级接回验收）  
执行对象：Cursor／SPHinXsys开发者  
前置：Doc208结果包及本轮讨论  
主目标：确认原始低耗散黎曼流程能否完成严格差温方腔，定位当前实现的速度低估，并完成可复现的定量验证。

> 优先级调整：本文件覆盖Doc208中“方腔定量欠账不阻塞继续法国校准”的执行顺序。下一轮先解决方腔；法国当前结果冻结保存，暂停新增热源形状、μ(T)、ε或h的匹配扫描。
>
> 不预设原始求解器必然通过，也不预设某项自定义修改必然错误。通过同条件对照回答问题；不再以一个开关的最终指标变化替代调用链与状态一致性证据。

## 0. 本轮必须回答的问题

1. 当前方腔在严格Ra=1e5、Pr=0.71下，低估来自测量、未稳定、热场、还是水动力？
2. 原始低耗散黎曼／连续性／黏性／壁面流程，在最小热—浮力接入后能否通过？
3. 如果原始版通过，当前OPHELIE附加流程中哪一项使结果退化？
4. 如果原始版也失败，两版共同的热边界、浮力接入或测量是否有误？
5. 最终选择哪套可复现的流体配置，后续如何接回法国工况？

本轮不开发新的电磁算法，不推进TEAM7，不修改独立A–phi／Laplace分支，不新增GPU隐式黏性，不启动法国长时耦合。

## 1. 已知事实、确认错误和待验证假设

### 1.1 当前严格方腔结果

Doc208 n160复验：

| 指标 | de Vahl Davis参考 | 当前结果 | 偏差 |
|---|---:|---:|---:|
| 平均Nu | 4.519 | 4.2641 | −5.64% |
| 水平中心线竖向速度极值 v* | 68.59 | 59.6208 | −13.08% |
| 竖直中心线水平速度极值 u* | 34.73 | 28.4194 | −18.17% |
| v极值位置 x/L | 0.066 | 0.0721 | +0.0061 |
| u极值位置 y/L | 0.855 | 0.8633 | +0.0083 |

当前物性已经通过`--dvd-exact`对齐；不能再把历史Ra≈1.203e5工况的较好数字作为严格通过证据。

### 1.2 已确认的测量错误：导热L2多除一次粒子数

当前源码累计：

`se = Σ(Ti−Tref_i)^2`

`sn = N * ΔT^2`

随后计算`√[se/(sn*N)]`，重复除以N。均匀权重下正确公式是`√(se/sn)`。

一般粒子布局推荐统一使用体积加权：

`eT = sqrt[ Σ Vi(Ti−Tref_i)^2 / (ΔT^2 Σ Vi) ]`。

必须修正并重算旧报告。此错误不直接解释速度低估，不允许把“修L2”写成“修好自然对流”。

### 1.3 已确认的代码顺序，尚未确认其物理影响

当前每个对流步存在：

`density_sum → density_reg → rho=ρ0*(DensitySummation/DensSumSigma0) → acoustic`。

最终Density重写可覆盖前一步对Density的修改；若正则化另更新Vol、EOS或其他字段，开关又会混入额外影响。

因此，Doc208“关dens_reg更差，故已排除该项”的结论需要降级为未隔离的对照。不要先判定密度求和本身有错；要核实每个阶段实际读写内容。

### 1.4 其他尚未证实的假设

- 背景压力／扰动压力／壁面压力使用不一致；
- 壁面黏性或压力外推产生额外阻力；
- Riemann耗散压弱低速回流；
- 温度边界或导热使驱动温差分布偏差；
- 末态没有达到足够稳定，或中心线取样存在偏差。

目前不得把其中任一项写成唯一根因。原始低耗散黎曼是首要对照，不是预先指定的正确答案。

### 1.5 法国现有环流的证据范围（R1补充，必须落实到当前状态报告）

当前法国几何的回流与方腔共享部分关键压力、连续性和壁面实现，但配置不是完全相同：法国还涉及变物性、高黏度冷区、CPU隐式黏性及滑移顶面。共享范围需由实际调用链证明，不能笼统写成“同一个算法，所以所有结论互相适用”。

| 内容 | 当前允许结论 | 当前不允许结论 |
|---|---|---|
| β=0短跑零速度、β>0出现回流 | 有温差驱动机制的正面证据 | 整个自然对流求解已定量准确 |
| Doc208 D1/F1单时刻速度插值出现两个回流区域 | 主、外围回流形态已观察到，值得保存 | 稳定双涡及涡强度已验证 |
| 构造非均匀初温后释放温度 | 能完成该初态下的短时热—流推进 | 已摆脱初猜影响并预测出文献收敛状态 |
| 方腔速度仍偏低13%–18% | 共享计算路径有必须排查的精度风险 | 法国流场一定错误，或误差恰为同一百分比 |

现有法国环流既不能直接判为伪流，也不能作为求解器准确性的独立证明。它属于待验证的应用结果。

当前报告优先用`BUOYANCY_MECHANISM_SUPPORTED`与`FRENCH_CIRCULATION_OBSERVED_UNVALIDATED`。历史`BUOYANCY_ENGINEERING_USABLE`不删除原记录，但必须说明其只表示机制／推进可用，不表示自然对流定量收口。

法国旧结果冻结保存，仅允许无需新增仿真的数据整理；方腔通过之前，不用继续调Q、μ或边界来补偿尚未明确的流体误差。方腔通过后，也必须按第11节核实实际迁移与接回结果，不能自动将旧法国场升级为已验证。

## 2. 版本隔离：首先得到真的“原始版”

### 2.1 不覆盖当前工程

使用独立worktree或独立checkout。当前OPHELIE工作树、未提交代码、输入和法国结果全部保留。

禁止在当前目录执行`git reset --hard`、`git clean -fd`、批量恢复src、丢弃未提交修改，或为构造原始版强行切换当前分支。

允许在新目录建立只用于本轮验证的分支／工作树及独立build目录。使用实际绝对路径，不照抄开发机路径。

### 2.2 选择基线的顺序

1. 查清OPHELIE修改涉及的基础版本、提交和共享声学改动历史。
2. 优先选同一基础版本中、相关OPHELIE共享流体修改之前的干净提交。记录具体hash，不根据“原始”“main”等名字判断。
3. 若找不到对应干净版本，可选择可追溯的上游提交作版本对照；必须标记`UPSTREAM_VERSION_CONTROL`，列出额外版本差异，不能据此单独证明某一补丁是根因。
4. 不默认用最新上游替代历史基线，否则库升级和算法修改混在一起。需要拉取时先确认remote，避免更改当前分支。

Doc208包记录的HEAD为`67fef44ef8a2cd6f5f7712bf2e66038cc6de89d3`，它只是定位线索，不代表干净基础版本；当前服务器HEAD和未提交修改需重新确认。

### 2.3 必须验证原版身份

核实并记录以下实际编译来源：

- 声学两半步、Riemann solver／limiter；
- 密度求和及正则化；
- 黏性、壁面、ForcePrior更新；
- 状态方程、位置更新、排序和邻域配置；
- 核函数及平滑长度定义。

输出文件hash、compile_commands中的实际include路径及二进制hash。不能只有case源码干净，实际链接仍读当前修改过的共享头文件。

不要仅在已修改求解器上设置`PerturbationMomentum=0`就宣称恢复了原始代码。原始对照必须由原始共享实现编译。

## 3. 两个主版本＋一个后续修复版本

| ID | 名称 | 核心水动力 | 应用接入 | 用途 |
|---|---|---|---|---|
| O | ORIGINAL_MINIMAL | 干净原始低耗散黎曼流程 | 最小热输运＋Boussinesq | 检验原始路线 |
| M | MODIFIED_AUDITED | Doc208修改版 | 保持当前方腔物理设定 | 可追溯问题复现 |
| R | REPAIRED_SELECTED | 由证据选择O或局部修复M | 共用审计后的物理与后处理 | 最终定量验收 |

O与M在初始对照阶段均不按结果调参。修测量不算改变物理解，但需记录版本；修物理更新只能进入R或具名实验版本。

### 3.1 原版O允许新增的内容

- 独立方腔应用入口；
- Temperature、温度增量及必要热物性字段；
- 经审计的导热与热边界应用适配；
- 最小Boussinesq体力，严格遵循原版ForcePrior生命周期；
- 只读输出、监测、独立后处理接口。

禁止在建立O的过程中改共享声学方程、换limiter、增加当前特殊密度归一化、引入论文另一套壁面、加入未声明的压力投影或速度清零。

原版已有的标准密度处理、位置／输运更新、常规压力参考保留。原版流程中必须声明选择了哪个成熟示例的时间循环；不能把“尽量少代码”误解为删除原版必要步骤。

### 3.2 最小Boussinesq接入的压力约定

可采用减去恒密度静水项的压力表述：仅施加温差浮力，求解相应动态压力；不要同时再施加完整重力并重复扣背景。

若选择包含完整重力的等价总压表述，必须明确静水压、密度和初态如何匹配。O的首选是原版流程中最少修改、物理一致的表述，不为了与M同名而移植M的特殊压力开关。

总压与扰动压数值不要求逐点相同；比较物理压力梯度及运动。所需常量压力偏置、EOS参考和壁面压力处理必须明确，不能随时变化。

## 4. 两版共同的严格物理合同

### 4.1 几何与边界

二维正方腔，x,y∈[0,L]；y竖直向上。左壁303 K，右壁283 K，上下绝热；四壁无滑移、不穿透。初始流体293 K、速度0。无自由面，无EM，无热源，无Marangoni。

温度Dirichlet位于物理面x=0/L，不是壁粒子中心；角点热条件交界需单独标记。单位厚度二维测度，热量积分不能误用三维体积。

### 4.2 参数及定义

| 参数 | 值／定义 |
|---|---|
| L | 0.0386 m |
| ρ0 | 1.204 kg/m³ |
| cp | 1006 J/(kg·K) |
| β | 0.00343 K−1 |
| g | 9.81 m/s²，方向向下 |
| ΔT | 20 K |
| Tref | 293 K |
| Ra | 100000 |
| Pr | 0.71 |
| ν | 1.657702883504525e−5 m²/s |
| α | 2.334792793668345e−5 m²/s |
| k | 0.02827957066718147 W/(m·K) |
| μ | ρ0ν；输入动力黏度时不能直接填ν |
| 速度尺度 | α/L≈0.0006048685994 m/s |
| 无量纲时间 | t*=αt/L² |
| t*=0.5 | 约31.9077565 s物理时间 |

启动时独立计算Ra/Pr，不只打印目标标签。`--dvd-exact`不能因参数解析顺序被其他选项覆盖；将最终生效配置写入JSON。

### 4.3 数值设定

- O/M主比较先n=80，dp=L/80；正式选择R后n=160。
- 初始粒子位置、壁层、核函数、h/dp和实际物理壁位置尽量相同；不相同必须列为混杂因素。
- 声速主档共同采用Doc208约12.3071686 m/s或双方一致且事先声明的值。暂不以调声速凑结果。
- 原版低耗散limiter使用其原始实现和默认系数；输出实际类型与系数，不靠`scheme=riemann`字符串确认。
- 记录Real精度；O/M优先相同设备、精度和编译优化。先GPU，CPU仅作为有辨别力的参考，不要求所有组合都跑。
- 时间步规则及热子步稳定限制明确；不得在达到子步上限后静默丢失未推进的物理时间。
- 温度随粒子输运，导热叠加；若有transport velocity／shifting，必须说明温度与粒子更新一致性，不随手混入别篇文章修正。

## 5. 第一批：修测量，建立独立输出（必须先完成）

### 5.1 修复导热L2并验证误差程序

实现正确体积加权L2。同时用无需GPU的合成数组验证：

- 全部粒子误差为同一常数ε，输出应为|ε|/ΔT，与N无关；
- 将每个样本复制一次，误差不应改变；
- 精确场输出0；非均匀权重输出手算结果。

这是针对已发现的测量错误的必要验证，不要编写大量仅镜像实现的测试。

### 5.2 输出真实中心线，不再只取带状极值

主比较对象：

- u(x=L/2,y)，竖直中心线；
- v(x,y=L/2)，水平中心线；
- 同位置θ=(T−Tc)/ΔT及必要壁面温度剖面。

使用局部插值到几何中心线，保存带符号完整剖面；禁止在宽带内直接挑最大粒子作为中心线值。

至少两种合理插值／取样宽度作一次敏感性核对。输出覆盖率、最近样本距离和插值误差；不在无数据区任意外推。

参考定义：v正向峰在热壁附近x/L≈0.066；u正向峰在上半部y/L≈0.855。冷热侧／上下半部都输出，不使用“只要相近就优先选某半边”的90%逻辑选择有利位置。

### 5.3 Nu的测量合同

主Nu来自实际施加的热壁通量：先按接触面身份区分热壁／冷壁，再积算，不仅凭粒子距哪边更近分类。记录边界核完整支持域。

`Nu_hot = Qhot/(kΔT)`，`Nu_cold = −Qcold/(kΔT)`，这里单位厚度壁长L已与热流尺度中的L抵消；仍应保留完整面积定义以避免2D/3D混淆。

独立辅助：沿壁面分段重建流体温度法向梯度，采用物理Dirichlet约束并报告拟合窗口敏感性，再沿壁长积分。不要把所有高度混为一个未经解释的拟合，也不要将错误单侧差商乘2。

现有Contact≈Fit是有用证据，不单独等于边界准确；需与导热解析对照共同检查。

### 5.4 稳态不能只看全域Umax平台

固定t*=0.5是观察节点，不是自动稳态证据。

在末段连续记录中心线u/v极值、Nu_hot/cold、总动能、全域T和u的变化、壁热流差及内部热能变化。

建议窗口：Δt*=0.02，至少连续三个窗口。相邻窗口的平均速度剖面L2变化≤1%，Nu变化≤0.5%，温度场变化以ΔT归一≤0.5%，作为工程稳定候选；同时报告窗内波动，不能只靠均值掩盖持续大振荡。

若有缓慢稳定趋势，允许在同一配置续至t*=0.75，必要时最多1.0，再作判断；若始终波动或偏差平台，不通过无限延长替代根因定位。

### 5.5 记录频率，避免重跑后才缺数据

- 低成本标量：每Δt*≤0.002或更细；用实际时间戳计算时间平均。
- 中心线剖面：每Δt*=0.01；末段比较窗完整保存。
- 全粒子轻量CSV：t*=0、0.05、0.1、0.25、0.4、0.45、0.48、0.5；若续算补末窗。
- 因果诊断仅保存前数个对流步的各阶段，及稳定后一个代表步，不输出每个声学步全部粒子。
- VTP可低频或关闭，CSV必须含OriginalID、位置、速度、T、ρ、质量、Vol、压力及压力定义；不要依赖ParaView才能验收。

## 6. 第二批：两个廉价的公共接入检查

### 6.1 导热解析对照

复用一维热传导方腔：左热右冷、上下绝热、冻结速度，准确稳态为T=Th−ΔT·x/L、Nu=1。

O/M使用同一物理热边界和同源热更新适配。先n40或n80，不新增热算子路线。

检查修正后的温度L2、通量和独立梯度。目标：温度L2≤2%，梯度Nu≤3%，实际通量Nu≤5%，热量不平衡≤1%。保留角区和全域结果。

若失败，先修共同热接入；不得在该状态下解释自然对流速度偏低为水动力错误。

### 6.2 浮力接线与ForcePrior

在O/M上分别核对T=293、303、283 K的瞬时加速度：0、向上0.336483、向下0.336483 m/s²。该数值只用于浮力项，不是粒子总加速度。

核对黏性与浮力按原版生命周期分别贡献，前后更新不会重复累加或覆盖；不得每步清零总ForcePrior来掩盖问题。

均匀293 K、无热驱动的闭腔短检查应无明显自发流动。空间均匀非零浮力在闭腔会产生压力调整，不能用自由加速解析解验收容器整体运动。

瞬时力对照是必要的接线检查，不再把它当作完整自然对流准确性的证明。

## 7. 第三批：O/M同条件方腔主对照

### 7.1 固定执行顺序

1. 复用Doc208的M基线及原始数据；如缺独立后处理所需字段，补一组M n80，不无条件重跑所有历史n160。
2. O n80先短smoke，确认热、力、边界、时间推进健康，再完成至t*=0.5及末段观察。
3. 统一后处理得到O/M完整剖面与误差。
4. 根据第8节决策选择后续，不自动做所有候选×分辨率矩阵。

O的最小应用适配若不能编译，修适配，不把编译失败当作原始算法不适用。适配差异必须归档。

### 7.2 方腔工程验收门槛

这是项目验收目标，不是宣称与基准达到相同高精度：

- 实测Ra、Pr相对目标误差≤0.1%；
- 稳定性证据满足第5.4节；
- 平均Nu相对4.519误差≤5%；
- v、u中心线极值相对68.59、34.73各≤10%；
- 对应位置误差≤0.05L；
- 热壁／冷壁通量不平衡≤1%，并与内部热能变化相容；
- 无逃逸、无质量重置、无持续密度失稳。

所有条件共同判定。不得只通过Nu就宣布自然对流已验证；也不要求本轮追求参考表末位精度。

完整剖面与当前只有极值的参考数据分开：若没有可信全剖面参考，O/M可相互比较形状，不能捏造对文献剖面L2。

## 8. 决策分支：每一步都要能缩小问题

### 8.1 O通过、M失败

原始路径优先成为R，先验证n160，再考虑法国接回。不得因为已经开发很多特殊开关而强行保留M。

根因归因分组最多检查：

- Gρ：自定义初始排布归一化／密度重写；
- Gp：背景／扰动压力、声学密度模式及壁面压力配套；
- Gw：壁面声学、黏性和ForcePrior定制。

这些组之间有依赖，不机械逐行cherry-pick。每次移入必须形成物理一致且可编译的完整组；无法独立移入时将组扩大并明确降低归因粒度。

最多两次有依据的组对照。若只能证明“这一组使结果退化”，就保留组级结论，不虚构唯一代码行根因。原始版已过时，不必为了精确归责重跑数十组合。

### 8.2 O、M都失败且偏差接近

优先查共同项：物性、热边界物理位置、导热、浮力、同一测量误差、同一未稳态初态。

先核对第6节与末段剖面。此时不能说“原版也坏了，所以必须换ISPH／δ-SPH”。

若共同接入均有证据，再执行第10节受控诊断中的一个；不是启动所有诊断。

### 8.3 O、M都通过

核对旧／新测量、二进制、热边界、实际配置差异。只有旧场用新测量就恢复时，才能归为后处理；重新编译后恢复必须解释构建／代码差异。

### 8.4 O失败、M较好

保留M，同时确认O是否省略原版必要步骤、是否用了不同核／壁层、是否接入了不一致热路径。不把此结果预先解释为自定义方案更先进。

### 8.5 两次局部尝试仍无改善

停止堆开关，提交最小失败例、实际方程、阶段数据和已排除项。下一决策才考虑成熟替代壁面或独立参考；不在本轮自行扩展为框架重写。

## 9. 密度—压力一致性专项：不是再做一次无效开关

### 9.1 必须阅读实际实现

查清DensitySummation、DensityRegularization<Internal>、自定义归一化、EOS和声学两半步的真实读写字段。不要仅根据类名或旧文档猜测。

为每个阶段列出Density、DensitySummation、DensSumSigma0、Pressure、BackgroundPressure、VolumetricMeasure、Mass的读写表和host/device同步。

### 9.2 最小阶段快照

按稳定OriginalID记录：

- 对流步开始；
- density_sum后；
- density_reg后；
- 自定义归一化后；
- 第一声学步前／后；
- 对流步末。

每阶段检查：总质量、ΣVol、ΣρVol；ρVol与Mass的逐粒子差；Pressure与该压力定义下EOS是否一致；速度和动量增量。

若Pressure储存总压，EOS残差应按`p−pbg−pEOS(ρ)`等实际定义计算；不得用错误压力口径制造FAIL。

### 9.3 正确解释密度正则化消融

- 若reg只改Density且随后完全覆盖，之前对照不能用来衡量reg阻尼。
- 若reg另改Vol／EOS而归一化只改Density，先检查状态不一致的影响。
- 若有清楚且一致的设计理由，保留并量化；不因“有覆盖”就必然删除。
- 禁止只关密度求和／正则化后保留所有其他依赖项，跑崩后宣布原始方法无效。

修复应使用原版认可的一套密度流程，或完成一致的状态更新；不为让极值通过而直接乘压力或密度系数。

## 10. 若还失败：只选择一个有辨别力的诊断

### 10.1 首选：固定可信温度场，分开热与流动

若已有独立可信的方腔参考温度场，可映射后固定驱动O/M。比较同一温度下的速度，判断完整热演化是否是主要差异。

必须说明：仅有de Vahl Davis几个极值不等于已有完整温度场。不要凭温度云图随意补场当作独立权威参考。

若只有M自己的温度场，也可以用于O/M相同驱动的比较，但只能证明两条水动力路径是否一致，不能证明它们都准确。

该诊断不替代自由演化方腔最终验收。没有现成可信参考时，不为了本项另开发一套完整CFD求解器。

### 10.2 如怀疑黏性／壁面，复用相容的小测试

先检查实际μ／ν系数、内外壁贡献、ForcePrior生命周期。已有Couette／衰减测试只有在调用同一实际算子时才能复用结论。

需要新增解析驱动时，优先封闭方腔内满足无滑移、无穿透的制造解；不要重用有限闭容器中的水平单向正弦速度作为严格解析解，因为它与侧壁不相容。

可选流函数`ψ=UL sin²(πx/L)sin²(πy/L)`，定义`u=∂ψ/∂y`、`v=−∂ψ/∂x`，四壁速度均为0。用独立符号脚本生成稳态动量所需体力，求解器只接收体力，不能每步回写速度。此项仅在证据指向壁面／黏性时触发，不要求首轮实现。

### 10.3 如怀疑Riemann耗散，先量化而非换方案

分别记录声学两半步中压力平均、速度／压力跳跃耗散及壁面项的实际动量／动能变化。分裂算子的动能贡献可以分析，但不能将不同采样阶段的力直接相加当作精确预算。

只有证据支持耗散主导后，才允许一个OPHELIE局部、具名的强度对照，原版默认保持不变。

不能直接用普通估计`νnum≈c0h`归咎低耗散limiter，也不能用调小真实μ、调大β／ΔT或改Ra实现通过。

声速或精度敏感性也只在有证据时择一检查；不要同时形成c0×精度×网格×方案矩阵。

## 11. 最终选定R后的验收和法国接回

### 11.1 方腔正式收口

R至少完成n80、n160的同配置对照，使用独立后处理、稳定窗口和完整剖面。细级达到第7.2节工程门槛，中细级变化可解释，残余误差如实报告。

若只有n160偶然尾点通过、窗平均／其他时刻明显失败，不能收口。若门槛未过但根因已定位，报告`ROOT_CAUSE_LOCALIZED_NOT_CLOSED`，不要改门槛或授予完成标签。

### 11.2 两类结论分开

- O通过，证明干净原版＋本次最小接入能准确求解这一标准自然对流问题；不自动覆盖所有高Pr、强温变黏度、滑移顶面工况。
- R接回法国后能维持回流，证明修复没有破坏应用；不自动证明法国文献已定量复现。

### 11.3 先提交方腔—法国实际算子映射表（R1补充）

在运行法国回归前新增`CAVITY_TO_FRENCH_OPERATOR_MAP.csv`，按实际执行路径逐项列出：

| 项 | 方腔最终R | 法国当前M | 法国拟接回R | 需要的证据 |
|---|---|---|---|---|
| 声学两半步及limiter | 实际模板／源hash | 同左 | 同左 | 确认已迁移R，而非仅同名 |
| 密度求和、正则化、重写 | 实际执行及顺序 | 同左 | 同左 | 没有重新引入已排除的问题组合 |
| 压力、EOS、背景分离 | 场定义及更新相位 | 同左 | 同左 | 几何／初态改变后仍一致 |
| 位置／输运／排序 | 调用顺序及字段跟随 | 同左 | 同左 | 热、力、μ字段不错位 |
| 浮力与ForcePrior | 力项及生命周期 | 同左 | 同左 | 无重复施加、覆盖或过期值 |
| 黏性 | 空气常μ及实际算子 | 玻璃变μ／implicit_split | 明确选择 | 方腔通过不替代变μ隐式验证 |
| 壁面 | 四壁无滑移 | 侧底无滑移＋滑移顶 | 明确选择 | 不重新照搬未经验证的顶面补丁 |
| 热边界 | 左右定温、上下绝热 | Robin＋顶面辐射 | 保持同物理条件 | 独立热边界误差仍需披露 |

表中必须填写真实源文件、配置和执行顺序；不能用“共用SPHinXsys”替代。差异是预期存在的，但每个差异都应说明为什么不会使方腔证据失去适用范围，或明确还需哪项局部检查。

### 11.4 法国分级接回，不一次恢复全部修改（R1替换原短回归）

方腔完成正式验收后才启动。使用Doc208 D1未变形Q基线；F1/F2人为热源变形暂不作为验证基线。

| 阶段 | 固定内容 | 变化内容 | 本阶段回答的问题 |
|---|---|---|---|
| FR-A | 相同几何、同一给定T、常μ=4、同一壁面物理条件 | 法国旧M与接回R水动力对照 | 流体路径修正是否显著改变回流和伪流 |
| FR-B | 同一给定T、R水动力、同一壁面 | 恢复μ(T)和已有可靠隐式路径 | 冷区黏性处理是否稳定、是否正确改变流型 |
| FR-C | 同一初始T、R与变物性、原始固定Q=50 kW | 解除固定T，恢复导热及热边界 | 温度释放后回流能否维持、热场怎样演化 |

执行细节：

1. **初态一致**：优先复现Doc208原始映射后的t=0粒子温度。若改用其t=20场作初态，必须另命名并重置新运行物理时间；不得把不同热历史的“20 s”直接相减当修复效果。
2. **FR-A是诊断**：常μ和固定T不一定保留原来的两胞形态，不能因此直接判R失败；应比较同条件M/R、噪声、边界及回流变化，并解释原工况为何依赖变μ。固定T结果不能称为热方程预测。
3. **FR-B是新增适用范围**：复用既有变μ／隐式验证，只回归实际受影响路径。没有可复用证据或出现异常时，补最小相容测试；不因为方腔常μ已过就跳过这一层。CPU隐式可以保留，不要求本轮GPU化。
4. **顶面处理明确**：四壁无滑移方腔的成功，不能直接证明滑移顶面准确。记录真实壁面法向／切向响应；若出现明显异常，先做局部、相容的顶面检查，不把旧特殊开关全部恢复。
5. **FR-C释放真实温度**：热源、导热、损失及物性正常更新，不每步回写初猜。固定总功率、Q空间映射和面积归一化全部审计；live EM仍关闭。
6. **时长有限**：FR-A/FR-B先1–5 s，健康后FR-C到5 s、20 s；是否继续50 s依据变化趋势与预算，不默认长跑。β=0控制仅在迁移影响浮力／压力装配或出现伪流时重跑，其他已有证据可复用。
7. **统一后处理**：至少保留两个时刻的直接速度矢量、流线、T/Tmax以及主次胞位置。固定环带符号仅作辅助，`OTHER`或预设环带内速度变号不等于没有环流。
8. **修正图8评分**：沿用Doc208参考数据时，缺失等温线仍纳入误差与覆盖率报告，不因`CONTOUR_ABSENT`从总分删除；不以假定热点(0.85R,0.70H)作为作者精确目标。

本轮法国回归的最低验收：无新的显著非物理流、质量不重置、实际压力／密度／体积一致，热量账本受控，固定T和释放T结果分开，回流变化有直接图件及适用范围解释。它不要求本轮法国文献全部定量吻合。

如果R需要共享代码修改，只迁移已验证的最小补丁，并回归一个受影响的原版流体例子。不得用“OPHELIE局部”名义忽略共享修改范围。

### 11.5 迁移后的结论升级条件（R1补充）

- 方腔通过，但法国尚未接回：只授予方腔工程定量通过；旧法国结果仍为待验证观察。
- 法国FR-A/FR-B通过：可说明选定水动力与玻璃黏性／边界组合能够推进；不说明温度场已预测准确。
- FR-C通过并有跨时刻形态证据：可授予`FRENCH_SHORT_REGRESSION_ONLY`，明确实际时窗与剩余不稳定性。
- 若需要声称稳定双胞或文献趋势匹配，另需分隔／涡中心证据、插值敏感性和温度演化的稳定性；不能仅由一张好看的流线图授予。
- 绝不以参数调整抵消未解释的速度低估；也不以法国两胞存在反向证明方腔缺口已解决。

## 12. 任务预算、停止规则与交付节奏

### 批次A：无需完整长跑

冻结版本；确认O的源文件来源；修L2和中心线后处理；完成公共热／力小检查；输出密度读写审计。发现明确实现错误即可局部修复，不等到全部长跑结束。

### 批次B：唯一主对照

O/M n80同条件结果＋末段窗口。每次启动前估计粒子数、声学步数、GPU耗时及字段大小。优先前台可见，失败非零退出；不得把缺失数据写成0或PASS。

### 批次C：按决策树最多两项局部定位

选择R，完成n160正式验收。时间延长、分辨率加密和耗散参数变化不同时进行。

### 批次D：成功后最小法国回归

只有通过方腔工程门槛才继续应用级校准。若未过，本轮仍需给出可审查失败例和下一唯一动作，不能仅回复“建议再跑”。

## 13. 输出结构与运行清单

建议目录：

```text
doc209_runs/
  00_entry/
  versions/
  measurement_audit/
  thermal_and_force_checks/
  original_O/
  modified_M/
  density_pressure_audit/
  targeted_experiments/
  repaired_R/
  french_short_regression/
  scripts/
```

必须交付：

1. `ROUND209_DECISION.md`：原版通过与否、当前偏差、已证实根因／未证实假设、最终R及下一步。
2. `VERSION_AND_SOURCE_PROVENANCE.json`：O/M/R提交、未提交diff、关键源文件hash、实际include路径、设备、精度、二进制hash。
3. `PHYSICS_CONTRACT.json`：最终实测Ra/Pr、μ/ν/α/k、声速、几何、核、壁位置、热条件、压力定义和时间推进。
4. `MEASUREMENT_CORRECTION.md`：L2修正、旧数据可否重算、中心线插值、Nu定义和误差程序的小验证。
5. `DENSITY_PRESSURE_STAGE_AUDIT.md`＋阶段CSV：真实读写和归一化影响，不能仅报开关名称。
6. `CAVITY_COMPARISON.csv`：同定义的O/M/R n80/n160、稳定窗口、Nu/u/v/位置和误差。
7. 中心线带符号剖面CSV、温度／速度图、必要原始粒子场及可复现脚本。
8. `CAUSE_EVIDENCE.md`：每个因果结论对应哪一项受控改动；不能只用相关性排序嫌疑。
9. `REMAINING_LIMITATIONS.md`：高Pr、变物性、滑移顶面、EM和法国文献仍有哪些未覆盖。
10. `CAVITY_TO_FRENCH_OPERATOR_MAP.csv`与`FRENCH_TRANSFER_REVIEW.md`：实际路径差异、FR-A/B/C结果及允许升级的结论；若方腔未过，法国运行标NOT_RUN并保留旧结果状态。

打包`DOC209_ROUND_PACK_<实际日期>.tar.gz`。不要嵌套历史大包；包含复现必要输入、少量旧基线及其来源索引。

## 14. 标签和禁止事项

允许标签：

- `MEASUREMENT_VERIFIED`
- `ORIGINAL_SOURCE_VERIFIED`
- `CAVITY_ENGINEERING_QUANTITATIVE_PASS`
- `ROOT_CAUSE_LOCALIZED_NOT_CLOSED`
- `FRENCH_SHORT_REGRESSION_ONLY`
- `BUOYANCY_MECHANISM_SUPPORTED`
- `FRENCH_CIRCULATION_OBSERVED_UNVALIDATED`

本轮方腔通过不授予`FULL_FRENCH_VALIDATION`。若仅SMOKE_PASS、形态正确或浮力公式正确，不授予方腔定量通过。

禁止：

- 覆盖当前工作树或丢弃未提交内容来“恢复原版”；
- 使用修改过的共享库却把新case叫原版；
- 同时更改物性、核、壁面、密度和limiter后宣称定位单项根因；
- 为匹配基准调真实ν、β、ΔT；
- 用当前自己的温度场冒充独立参考；
- 在不相容边界下做解析速度验证；
- 继续把已知L2错误、密度覆盖疑点当作已通过证据；
- 法国图像匹配掩盖方腔速度低估；
- 原版通过后仍无依据坚持保留特殊修改；
- 因原版一次失败直接转ISPH／δ-SPH或重写框架。

**本轮收口标准：得到一条版本清楚、状态更新一致、独立测量可信、标准方腔定量达标的流体路径，并明确当前OPHELIE附加修改的必要性和影响。**


<!-- ===== END 209_OPHELIE_CLEAN_RIEMANN_CAVITY_ROOT_CAUSE_AND_CLOSURE_CURSOR_PLAN_2026-11-03.md ===== -->

---

<!-- ===== BEGIN 210_OPHELIE_PARTICLE_MOTION_NATIVE_WORKFLOW_AND_CAVITY_TEST_PLAN_2026-11-04.md ===== -->

<a id="210-ophelie-particle-motion-native-workflow-and-cavity-test-plan-2026-11-04-md"></a>

# 210_OPHELIE_PARTICLE_MOTION_NATIVE_WORKFLOW_AND_CAVITY_TEST_PLAN_2026-11-04.md

# Doc210｜物理位置更新、原生推进流程与方腔定量测试任务书

日期：2026-11-04  
修订：R1（transport velocity实际启用审计与ON/OFF受控对照列为必做）  
执行对象：Cursor／SPHinXsys开发者  
依据：`DOC209_ROUND_PACK_2026-11-04.tar.gz`、Doc209 R1及后续审查。  
本轮主线：核对并恢复完整、相容的原生低耗散黎曼推进流程，定位粒子坐标变化后误差增长，再完成方腔定量验收。

> 本文接续Doc209，不重复全部历史验证。保留干净声学核心、已有受控对照及法国结果。方腔未通过前，法国匹配扫描和live EM仍暂停。
>
> “移动粒子”在Doc209冻结位置对照中指`update_pos.exec()`控制的物理坐标更新，不等于particle shifting，也不等于transport velocity。用户已确认仓库存在transport velocity实现：本轮必须核实方腔是否真正执行，并完成原生TVF开启／关闭对照；仍不预设它是唯一根因，也不增加第二套shifting。

## 0. 本轮必须改变的工作方式

不再只更换声学头文件、切换若干开关，然后把结果称为“原生求解器完整验证”。必须核对实际时间循环、状态变量与粒子几何更新的整体流程。

本轮结束至少回答：

1. 当前O入口与同版本成熟CK流体示例相比，缺少或改变了哪些必要步骤？
2. 物理坐标、位移变量、物理速度／输运速度是否按同一约定推进？
3. 误差首先来自状态同步、粒子排布、邻域／核几何、力更新，还是仍未定位？
4. 修复是否改善移动粒子制造解，并同时改善自由演化方腔？
5. 方腔是否真正达到定量与稳态门槛，还是仅完成了局部定位？

本轮不开发新δ-SPH／ISPH框架、不改A–phi分支、不新增电磁／TEAM7任务、不新增GPU隐式黏性，不把冻结坐标作为生产修复。

## 1. Doc209保留事实与需要纠正的结论

### 1.1 n80可比较结果

以下为同一新中心线测量口径下的结果，不能与Doc208旧带状极值直接混比：

| 配置 | V* | 相对68.59 | U* | 相对34.73 | Nu_contact |
|---|---:|---:|---:|---:|---:|
| M全套 | 53.22 | −22.4% | 26.43 | −23.9% | 4.231 |
| O干净声学＋当前driver | 61.55 | −10.3% | 29.44 | −15.2% | 4.894 |
| M关dens_reg和dens_norm | 63.86 | −6.9% | 29.66 | −14.6% | 4.894 |
| O单独加dens_norm | 54.65 | −20.3% | 24.85 | −28.5% | 4.244 |

自定义packing dens_norm解释了O/M一部分速度差；不要再以相同方式重复该矩阵。它不解释关闭后仍有的残余误差。

### 1.2 移动粒子制造解

- 移动位置，t=0.05 s：内部速度L2约0.12%。
- 移动位置，t=0.5 s：内部速度L2约167%。
- 冻结位置，t=0.5 s：约1.2%。

这是定位线索，不是`UpdateParticlePosition`实现有错的证明。冻结位置改变了输运和离散几何，不能把两者看作物理方程完全相同的精度对照。

### 1.3 必须降级的排除项

| 原报告说法 | 本轮正确处理 |
|---|---|
| 净声学半步动能变化很小，排除Riemann耗散 | 半步包含多个贡献；净值可能抵消。保留“未分离，不能排除” |
| 黏性瞬时内部L2≈3.8%，排除黏性问题 | 只支持内部瞬时算子基本合理；近壁、状态变化和长期推进仍未覆盖 |
| dens_reg单独运行发散，证明原生正则化有害 | 还需核对其初态、求和尺度、EOS／Vol配套，不推广到原生方法 |
| n160更差，所以网格加密无效 | 先确认实际终时、稳态窗口、配置及粒子状态；当前不足以证明收敛趋势 |
| 干净声学文件=干净完整流程 | 声学核心与应用时间循环分别验明 |

Doc209的`ROOT_CAUSE_LOCALIZED_NOT_CLOSED`只表示部分退化因素已定位，不表示剩余方腔误差的唯一根因已经找到。

## 2. 明确三个不同操作

| 操作 | 含义 | 本轮要求 |
|---|---|---|
| 物理位置更新 | 按积分约定推进实际Position，基本材料运动为dx/dt=u | 必须存在且与原生流程一致 |
| Transport velocity／TVF | 用原生输运速度及配套修正维护粒子排布、推进位置／方程 | 确认方腔实际执行状态，必须做原生配套ON/OFF对照 |
| Particle shifting | 额外重定位及必要场重映射 | 本轮不自行再加第二套 |

`UpdateParticlePosition`的名字不能证明其内部采用u还是输运速度，也不能证明它只做简单Euler更新。必须阅读同版本实际实现。

冻结位置是诊断；不能以“冻结后误差小”改成固定粒子生产求解器。若另走Eulerian SPH需完整方程和独立验证，不属于本轮。

## 3. 版本、入口与修改边界

### 3.1 保留三个版本身份

- `O209`：Doc209原始声学＋当前应用driver，保持可复现。
- `N210`：按同版本成熟原生CK示例恢复完整必要流程的候选。
- `R210`：经过证据选择的最终配置，初始可以等于N210，不预授通过。

不覆盖O209和M，不重置当前工程、不丢未提交代码。使用独立worktree／build或已有隔离目录。禁止将新build链接回另一版本共享头文件。

### 3.2 本轮定位入口

包内主要参照：

- `src_snapshot_O/validation/test_2d_ophelie_cavity_original.cpp`
- `src_snapshot_M_cavity.cpp`
- `src_snapshot_O/shared/electromagnetic_ophelie_laplace.h`
- `src_snapshot_O/shared/electromagnetic_ophelie_boussinesq.h`

包中没有完整核心推进源码，Cursor需从实际同版本仓库补齐；不得凭这些应用快照猜测原生UpdateParticlePosition／AdvectionStepSetup的内部实现。

确认真实编译入口，不因目录名natural_convection或original而认定算法身份。归档HEAD、未提交diff、关键共享源码、编译命令、Real精度与二进制hash。

## 4. P0：先核对完整原生流程，不启动新参数矩阵

### 4.1 选择一个参考示例

在同一基础提交中选择成熟、可运行的CK单相低速黏性流动示例，尽量是有壁面的受限流动。记录选取理由和调用链。

- 优先复用同版本lid-driven cavity或已验证的有壁黏性示例；先确认存在，不编造目标名。
- 若选Poiseuille，必须保留其真实周期／入口边界，不把该几何改成封闭容器后仍沿用解析解。
- 不将自由面密度处理直接移植成闭腔方案；区分应用边界选择与通用核心步骤。
- 没有完全相同的原生自然对流示例也可继续：使用原生水动力循环，只做最小热和浮力接入。

### 4.2 提交逐阶段差异表

对参考示例、O209、N210并列填写：

1. 初始粒子／壁粒子生成与松弛；初始核求和和压力参考。
2. DensitySummation的实际输出；DensityRegularization是否写Density、Vol、EOS或其它字段。
3. 每步Density、Mass、Vol、Pressure更新时刻。
4. AdvectionStepSetup读写的Position／Displacement／Velocity等实际变量。
5. TVF／transport correction／粒子正则化是否存在、作用变量和时机。
6. 黏性、浮力和其他ForcePrior的生命周期。
7. 两个声学半步分别更新什么；位移在哪里累积。
8. UpdateParticlePosition到底覆盖、累加或重建Position，使用什么速度与dt。
9. 何时更新cell list、邻域、核几何、修正矩阵和壁面信息。
10. 排序及字段跟随、host/device同步。
11. 热子步在哪一套粒子位置和邻域上执行。

输出`NATIVE_WORKFLOW_DIFF.md`，列真实类名、源文件及实际代码片段，不只列构造过的对象。

### 4.3 密度流程专项澄清

O209的`dens_sum=1`不是完整密度流程已经保留的证据。从源码确认求和后是否真正消费结果。

关闭reg和norm之后，若求和结果仅存入DensitySummation而不进入密度或其它更新，报告为“声学连续性推进＋未消费的求和诊断”，不能称作保留全部原生密度修正。

若恢复原生正则化，需要同时检查初始排布、核归一和EOS／Vol配套。Doc209初始求和有约1% packing偏差线索，不能在不核对这些因素时直接开reg，以其发散否定原版方法。

自定义`rho=ρ0*(当前sum/初始sum)`保持关闭，不作为默认精度修复。

### 4.4 选择候选的规则

- 找到缺失／不一致的原生必要步骤：形成最小完整补丁N210，记录变更组。
- 原生示例确实使用TVF：恢复其配套步骤、系数及更新约定，不只增加一个位置修正调用。
- 所选参考示例不使用TVF：继续定位同版本仓库中实际使用TVF的成熟入口，核对其适用条件；按第4.5节完成用户要求的受控测试，不把该示例的选择当作取消测试的理由。
- 发现多个独立疑点：先修可确定的生命周期／同步／积分错误。只有相互依赖的原生步骤才能作为一组恢复，禁止把所有疑点一起改掉。

### 4.5 transport velocity专项（R1新增，必做）

#### A. 确认实现与启用状态

先在同版本仓库定位已有transport velocity／transport correction实现和成熟调用入口。术语可能不同，不只按类名猜测。记录：

- 类与源文件、实际作用字段、配套动量／密度／热输运处理；
- 原生默认系数和必要前置字段／修正矩阵；
- 方腔O209是否构造、是否实际执行、执行频率与阶段；
- 是否只修改位移，还是维护独立输运速度；Position最终是否消费了该修正。

“仓库有实现”“case include了头文件”“对象已构造”都不能作为已开启的证据。新增或导出执行计数，并输出实际修正位移／速度的RMS和最大值。

均匀静止排布下修正接近零可能正常；应在移动MMS／方腔的非均匀排布状态核对实际效果，不能为了显示非零而修改默认系数。

#### B. 定义可比较的ON/OFF配置

- `TVF_OFF`：正常物理坐标更新，关闭TVF及其专属配套项；不是冻结Position。
- `TVF_ON_NATIVE`：启用同版本原生TVF及保持离散一致性所需的完整配套步骤，使用原生默认系数，不额外加particle shifting。
- 两者保持相同几何、物性、声速、limiter、核、分辨率、初态、壁面物理条件、精度、设备和后处理。
- ON若需要不同的密度／场修正规则，将这些明确列入“原生TVF组合”，不能声称仅改了一行坐标更新；自定义packing dens_norm仍不作为稳定化捷径恢复。
- 配置开关先核实已有CLI；不存在则新增OPHELIE局部清晰选项，本文不假定现有参数名。summary必须记录请求状态、实际生效状态、调用次数和参数。

若某一现有TVF实现无法用于当前CK／设备／边界，交付具体缺失接口和最小适配方案。复用原生数值公式，不自行设计第三种方法，也不静默改成shifting充数。

#### C. 防止测试体力覆盖TVF配套项

Doc209移动MMS里有直接写`ForcePrior = ViscousForce + manufactured_source`的诊断代码。必须核对TVF配套项是否也使用ForcePrior、速度／位移或相关历史字段；若是，不能让这条写入覆盖它们。

制造源仍来自解析物理场；以原生认可的独立体力机制或一致组合接入。保留黏性和必要输运项，防止重复累计；不能为了让ON通过而用数值解反算抵消源。

温度及其它材料字段也必须按原生输运约定处理。空气方腔若ON只修坐标、温度仍按不相容路径更新，不能用它评价TVF对自然对流的真实影响。

#### D. 测试顺序与预算

1. **移动MMS配对**：OFF和ON同初态，先0.05 s，健康后各到0.5 s；两者都正常移动坐标。原OFF结果在代码／配置和测量相同、字段足够时可复用，但需补一段执行状态审计。
2. **n80自由方腔配对**：OFF和ON使用相同Ra/Pr和温度自由演化，先短smoke，再完成相同物理时长和稳态观察。已有合格OFF基线可复用，不强制重跑历史M。
3. 比较速度L2、排布／核矩、压力波动、Nu两口径、中心线速度、稳态变化和实际耗时。不能只比较粒子看起来是否整齐。
4. ON明显改善MMS而方腔无改善：检查热输运及边界配套，不能直接宣布自然对流已修好。
5. ON/OFF均失败：记录首个异常阶段和修正是否真实进入Position，不开启TVF系数扫描、shifting叠加或新的limiter矩阵。

本轮TVF测试固定一套原生参数，不做系数调优。n160只对选定且在n80有可信改善的配置执行。

## 5. P0：建立不会误判的记录和小检查

### 5.1 时间与完成状态

新增并实际计算：

- `requested_end_time_s`、`actual_end_time_s`；
- `requested_tstar_end`、`actual_tstar_end`；
- 对流／声学／热步数，真实dt累计；
- `termination_reason`、是否达到目标、是否达到稳态；
- 触发步数上限、NaN、无进展、墙钟中断时明确非正常完成。

记录频率按物理时间阈值触发，禁止用`steps % 200==0`作为唯一观测条件。JSON目标时间不能冒充实际推进时间。

可用double影子计时核对Real时钟累计误差，先只读诊断，不同时改变所有计算精度。子步未覆盖完整对流步时不得静默丢弃时间。

### 5.2 物理位置更新一致性

在前几个对流步、MMS误差开始上升附近、末态各保存一组：

- Position更新前后；
- 实际位移变量、物理速度及存在时的输运速度；
- 所用dt与速度采样相位；
- 原生更新公式预测的位移与实际位移差。

应按真实分裂积分公式核对，不强制把所有更新拟合成`Δx=u_end*dt`。检查重复累加、使用过期位移、错dt、host/device不同步。

建议一致性残差用double离线重算，并按float舍入规模解释；不是事先规定任意固定小数门槛。出现相差一倍等量级错误先修，不能通过增加shifting掩盖。

### 5.3 状态一致性

每阶段记录Mass、Density、Vol、Pressure及实际参考压力，检查：

- Mass是否保持；逐粒子ρVol与Mass是否匹配；
- 压力与该阶段EOS／背景定义是否一致；
- ViscousForce、PreviousViscousForce、BuoyancyForce及ForcePrior是否对应同一粒子；
- 排序后OriginalID、T、温度增量、物性和力历史是否跟随。

若某阶段故意延迟更新，注明设计和使用范围，不对尚未进入正式更新阶段的字段错误套EOS判据。真正被下一算子消费前必须满足其状态约定。

不能每步清零速度、重置质量、强制密度或解析压力来获得稳定。

### 5.4 邻域与离散几何

少量固定OriginalID样本：内部、四面近壁、角区，以及误差增长最快区域。分别记录：

- 最近邻距离／dp、邻居数；
- 核求和S_i=ΣV_jW_ij；
- 一阶核矩／线性重现矩阵的误差与条件数；
- 壁面法向、物理距离，最小／最大粒子间距；
- 邻域列表与当前Position、参考位置或位移重建坐标的关系。

对固定种子抽样（例如64个粒子）使用CPU暴力几何查询，核对实际核支持半径内候选是否遗漏、重复，以及核几何究竟在哪套坐标计算。

重要限制：原生双时间步方法可能有意在对流区间冻结邻域／参考几何。不要因为未每个声学步重建就判错；按原生推导核对其参考坐标、位移与允许重建时机。

“每个声学步都重建邻域”不是本轮默认修复，只能在发现真实遗漏或过期几何时作短对照。

### 5.5 输出量与频率

MMS物理时刻：0、0.01、0.025、0.05、0.075、0.1、0.2、0.5 s。

- 轻量标量每0.005 s或更细；误差明显增长阶段自适应加密。
- 上述关键时刻全粒子CSV；至少包含ID、位置、u、ρ、Mass、Vol、p及误差。
- 阶段快照仅前几个对流步和误差首次超过1%／5%／20%附近的少量步；避免输出每个声学步全部粒子。
- 保留场误差空间分布图，不能只报一个最终L2。

## 6. P1：移动粒子制造解——主定位测试

### 6.1 保留已有制造场，不增加一批新变体

几何L=0.0386 m，U=1e−3 m/s，μ／ρ沿用严格空气方腔，四壁无滑移，T=Tref、浮力关闭。

`ψ=UL sin²(πx/L) sin²(πy/L)`

`u=∂ψ/∂y`，`v=−∂ψ/∂x`。

解析速度幅值不是U；最大速度为πU量级。`Umax/U≈π`不能判为发散。

采用稳态连续Navier–Stokes的源：

`f = ρ (u_ref·∇)u_ref − μ ∇²u_ref`，制造压力为常数。

连续方程可写为`ρ(u_ref·∇)u_ref = −∇p + μ∇²u_ref + f`。注意修正文档中将含对流源错误写成Stokes零左端的表述。

### 6.2 先验证参考生成器

独立脚本／符号计算检查div u_ref=0、四壁速度0、解析Laplace与对流导数，随机内点代回方程。

源f始终从解析场与当前物理坐标计算，不从数值速度构造“抵消误差”的源；求解器不能每步回写解析u。

移动粒子的参考速度应在其当前物理位置评价，不按初始位置比较。参考轨迹积分只在位置更新约定有疑点时触发，不默认另建轨迹矩阵。

### 6.3 必须处理TVF／shifting适用性

如果N210使用输运速度推进坐标，物理速度、密度及温度等变量必须使用原生相容的输运／修正方程。不能只把`dx/dt=u`替换成`dx/dt=u_transport`，其余不变。

若原生方法采用额外位置重定位，严格按其场修正／历史字段规则执行，记录修正位移和相关动量变化。不能声称坐标修正天然不影响温度、速度或ForcePrior历史。

本轮必须完成第4.5节的原生TVF ON/OFF配对；只用一种原生排布维护方法，不叠加新particle shifting，也不通过调大修正系数拟合MMS。

### 6.4 最小运行矩阵

| ID | 配置 | 时间 | 必要性 |
|---|---|---:|---|
| M0 / TVF_OFF | 已审计的无TVF流程＋正常移动＋只读观测 | 0.05→0.5 s | 既有O209确属OFF时可复用末态，定位首次增长 |
| M1 / TVF_ON_NATIVE | 同条件启用完整原生TVF流程，正常移动 | 0.05→0.5 s | 必做ON/OFF对照，检查修正真实生效 |
| M2 | 若必要，一项有依据局部对照 | 与M1同窗 | 仅用于区分尚存候选 |

M0仪器化若明显改变原结果，先解释二进制／配置／同步差异，不能把变化当作修复效果。

不重复全套冻坐标测试。旧冻结结果保留为定位辅助，标签改`FROZEN_GEOMETRY_DIAGNOSTIC_ONLY`，不授予拉格朗日NS准确性。

### 6.5 误差与判定

体积加权速度L2：`sqrt(ΣV_i|u_i−u_ref(x_i)|² / ΣV_i|u_ref(x_i)|²)`；全域与去3dp近壁的内部同时输出。

另报绝对速度误差、L∞、密度范围、散度及核矩。零参考速度点不单独作相对除法。

工程目标：t=0.5 s内部L2≤5%、全域≤10%，且误差没有明显失控增长；这是本轮定位目标，不是文献统一标准。

- 由167%显著下降，但尚未达门槛：记`MMS_IMPROVED_NOT_CLOSED`，不得仅因好转授予PASS。
- 满足：记`MOVING_MMS_ENGINEERING_PASS`，随后做方腔。
- 不满足：按首次异常阶段定位；最多再做M2，不继续扩展制造场、冻结模式、limiter组合。

MMS改善不自动证明方腔改善；二者必须分别比较。

## 7. P1：若发现原生流程差异，如何形成修复

### 7.1 位置／状态更新错误

发现重复位移、旧dt、力数组错位、EOS／Vol不同步等确定错误时，先修这些实现问题。保持物理参数和原始limiter不变。

提交修复前后的同一步数据以及MMS误差曲线，避免只提交最终PASS。

### 7.2 排布维护步骤遗漏

若参考示例和调用链证明某项原生步骤是完整流程的一部分，恢复整套配合。源文件／系数／作用时间保持可追溯；不要临时自制弱化密度重写替代。

必须报告修正前后最近邻、核矩、压力和速度误差的变化，不能仅以粒子更整齐作精度证据。

### 7.3 邻域／核几何不相容

修正实际坐标、参考坐标与位移的读取约定，或恢复原生正确更新时机。局部CPU参考证实后再改GPU路径。

不在没有遗漏证据时强制每声学步重建，避免把不同离散算法和额外耗时混入对照。

### 7.4 仍无明确差异

如果完整原生循环与当前O一致，或M1未改善，明确写出事实；不继续声称“肯定缺少transport velocity”。

重新检查MMS源、边界压力相容性和近壁误差。必要时只选一项相容的既有黏性／壁面回归。不能从冻结坐标较好直接推出新shifting方案。

## 8. P2：回到严格自由自然对流方腔

### 8.1 保留物理合同

L=0.0386 m，左303 K、右283 K、上下绝热、四壁无滑移；ρ0=1.204、cp=1006、β=0.00343、g=9.81。

Ra=1e5，Pr=0.71：

- ν=1.657702883504525e−5 m²/s；
- α=2.334792793668345e−5 m²/s；
- k=0.02827957066718147 W/(m·K)，μ=ρ0ν；
- 速度尺度α/L，t*=αt/L²；t*=0.5约31.9077565 s。

保持Doc209声速主档约12.307 m/s与原始低耗散limiter。输出实际Ra/Pr和源hash，不按希望的结果调整ν、β、ΔT或功率。

### 8.2 先n80，不同时重跑所有版本

按第4.5节完成TVF_OFF与TVF_ON_NATIVE的n80配对，真实温度演化；不得固定T或冻结位置。使用Doc209修正后的同一独立后处理，与O209和M的已有n80结果对比。若O209已核实为有效OFF且数据完整，可作为OFF基线。

先检查MMS与方腔的改善是否同向。如果只改善MMS而方腔不改善，不宣布共同根因已解决；下一疑点可能在热边界或耦合，而非继续叠加排布修正。

### 8.3 完成稳态观测

每Δt*≤0.002记录Nu、中心线u/v极值、总动能与平均温度；每0.01保存中心线剖面；末窗保存全场。

至少连续三个Δt*=0.02窗口，检查：速度剖面窗口均值变化≤1%、Nu≤0.5%、T场变化以ΔT归一≤0.5%，同时报告窗内波动。

若到0.5仍有清楚的缓慢演化，可分段续至0.75，最多1.0；若偏差已平台或持续不稳定，不靠无限加时解决。

n160只在n80与流程验证有明确改善后执行。必须记录实际终时和相同稳态条件；不把两个不同阶段的尾点当空间收敛。

### 8.4 Nu两口径不一致仍须保留

Doc209 O n80的Nu_contact≈4.89，而Nu_fit≈4.09；不能说“真值必然被两者夹住”，也不能取平均凑4.519。

主值继续由实际施加边界通量定义，辅助给沿壁分段梯度重建、窗口敏感性及全域能量账本。纯导热Nu≈1不自动保证对流热边界层准确。

若恢复输运步骤后Nu口径仍明显不一致，单列热边界／重建问题；不要把它混为MMS动量误差。

### 8.5 方腔定量验收

- 实测Ra/Pr相对目标误差≤0.1%；
- 满足稳定窗口，而非仅运行到目标时间；
- 平均Nu相对4.519误差≤5%；
- 中心线v*相对68.59、u*相对34.73，各误差≤10%；
- 极值位置相对0.066、0.855，各≤0.05L；
- 热流收支与内部能量变化相容，无质量重置和非物理逃逸；
- n80／n160变化有一致解释。

位置和速度必须按带符号真实中心线取样，不能宽带挑最大值或使用90%半区偏好。仅Nu通过不能授予自然对流定量通过。

## 9. 动能分析的边界：本轮不再错误排除

现有整半步ΔKE含压力、外力、黏性等贡献，不是纯Riemann耗散。每步很小、总和略正均不能排除大项抵消；第二半步不直接改速度也可能间接影响下一步。

本轮不把详细能量分解作为首要开发任务。若已有明确异常且需要定位，才记录同一相位的分项速度增量／功；无法分离时标`DISSIPATION_NOT_ISOLATED`。

禁止因为没有纯耗散证据就调limiter凑基准，也禁止把整步净值写成“已证明无数值耗散”。

## 10. 法国案例与电磁暂时如何处理

保持Doc209 R1原则：法国当前结果为`FRENCH_CIRCULATION_OBSERVED_UNVALIDATED`，不是已验证的预测结果。原场和原配置冻结保存。

方腔达到验收后，按实际算子映射，固定T／常μ→变μ／隐式→释放T／固定Q逐步短回归。高Pr、变μ、滑移顶面与空气方腔适用范围分开。

本轮不再调法国热源上移／径向压缩，不用类似双涡反向证明方腔正确；不将固定Q结果算作电磁已验证。

## 11. 批次、预算与停止条件

| 批次 | 必须内容 | 不自动展开 |
|---|---|---|
| A | 原生流程审计、TVF执行状态、时间与阶段观测、参考生成器核对 | 新TVF／shifting方法开发 |
| B | 移动MMS的TVF_OFF／ON_NATIVE配对与阶段定位 | 冻结模式／制造场大矩阵 |
| C | 必要时最多M2一次局部辨别 | dens×limiter×c0×网格扫描 |
| D | n80自由方腔TVF_OFF／ON_NATIVE对照，改善后n160 | 全部历史案例重跑 |
| E | 方腔通过后按Doc209 R1法国短回归 | 法国长时校准与EM攻关 |

若A发现确定的代码错误，修复后直接验证，不为了凑批次继续运行已知错误配置。

若B/C后仍无法定位，提交最小失败例及第一异常阶段，写`ROOT_CAUSE_NOT_YET_IDENTIFIED`。禁止把“与粒子移动相关”升级成“输运算法错误已定位”。

所有候选只在隔离目录中形成，旧结果保留。没有达标就如实报告，不降低门槛、不仅靠内部裁剪区域通过、不用修正系数补物理误差。

## 12. 必须交付的数据与文档

建议根目录：`docs/ophelie/reference/doc210_runs/`。

```text
00_entry/
native_workflow_audit/
time_and_state_audit/
moving_mms/
particle_geometry/
selected_fix/
cavity_n80/
cavity_n160/
french_regression/
scripts/
```

必须文件：

1. `ROUND210_DECISION.md`：O209实际缺失／不一致项、首个误差来源、修复因果证据、MMS和方腔是否达标。
2. `NATIVE_WORKFLOW_DIFF.md`：参考示例和O209/N210逐阶段差异，真实源文件／公式／调用顺序。
3. `POSITION_UPDATE_CONTRACT.md`：物理速度／输运速度／位移／Position的关系，TVF与shifting是否实际执行。
4. `TIME_AND_TERMINATION.json`：目标与实际终时、步数、Σdt、结束原因、完成状态。
5. `STATE_STAGE_AUDIT.csv`：ID、时间、阶段、ρ、Vol、Mass、p、位置／位移、各力项及同步信息。
6. `PARTICLE_GEOMETRY_HISTORY.csv`：邻居数、最近邻、核矩／求和和参考几何，误差首增时刻。
7. `NEIGHBOR_SAMPLE_CHECK.json`：CPU抽样查询方法、核支持定义、遗漏／重复及允许的冻结配置解释。
8. `MMS_REFERENCE_CHECK.json`：独立参考方程、边界、源项代回；不每步回写解。
9. `MOVING_MMS_COMPARISON.csv`：M0/M1必要M2的全域／内部误差时间史，配套误差空间图。
10. `CAVITY_COMPARISON.csv`：真实终时、稳态、Nu两口径、中心线极值与误差；旧数据与新数据标清测量口径。
11. 选定源码、最小diff、完整运行命令、输入、设备精度、编译与二进制hash；实际依赖的核心推进头文件必须随包。
12. `REMAINING_LIMITATIONS.md`：仍未覆盖的热边界、近壁、变物性、法国及EM内容。
13. `TVF_ACTIVATION_AUDIT.json`与`TVF_ON_OFF_COMPARISON.csv`：真实调用、修正量、原生配套项、MMS／n80方腔同条件结果和误差；未完成须给具体阻碍，不得只写“可选下一步”。

标签按证据分别使用：`NATIVE_WORKFLOW_AUDITED`、`POSITION_UPDATE_CONSISTENT`、`MOVING_MMS_ENGINEERING_PASS`、`CAVITY_ENGINEERING_QUANTITATIVE_PASS`。不得用任何单项PASS替代整体。

最终归档`DOC210_ROUND_PACK_<实际日期>.tar.gz`，不递归塞旧大包。完成报告必须明确区分：已查明、仅相关、尚未排除、未运行。

## 13. 最后检查：避免Cursor执行歧义

- 不把物理位置更新叫particle shifting。
- 不在未看实现前认定UpdateParticlePosition使用哪一种速度。
- 不把“原版示例使用TVF”推广成所有SPH都必须用TVF；但本轮已明确要求的原生ON/OFF对照不能省略。
- 不只恢复位移修正、遗漏配套动量／密度／温度处理。
- 不因冻结几何较好就把粒子固定作为修复。
- 不以原生正则化在不匹配初态中发散为由否定原生流程。
- 不把未消费的DensitySummation输出写成密度已经修正。
- 不用步数取模打印来证明稳态或实际终时。
- 不把u_max/U≈π的正常制造场幅值当作速度暴涨。
- 不因内部黏性快照较好就排除近壁和长期推进。
- 不把两个Nu口径的中间数当作真实答案。
- 不调整真实黏度／β／温差，也不继续法国调参掩盖方腔误差。

**本轮成功应表现为：完整推进链有可核对依据，移动粒子误差的来源被缩小或修复，且同一配置在严格方腔中取得可重复的定量改善／通过；不是又增加一批开关和PASS标签。**


<!-- ===== END 210_OPHELIE_PARTICLE_MOTION_NATIVE_WORKFLOW_AND_CAVITY_TEST_PLAN_2026-11-04.md ===== -->

---

<!-- ===== BEGIN 211_OPHELIE_CAVITY_CAUSAL_CLOSURE_AND_EM_VALIDATION_CURSOR_PLAN_2026-11-06.md ===== -->

<a id="211-ophelie-cavity-causal-closure-and-em-validation-cursor-plan-2026-11-06-md"></a>

# 211_OPHELIE_CAVITY_CAUSAL_CLOSURE_AND_EM_VALIDATION_CURSOR_PLAN_2026-11-06.md

# Doc211：方腔归因收口与电磁独立验证并行执行任务书

日期：2026-11-06  
执行对象：Cursor；项目：冷坩埚招投标／法国 OPHELIE 复现  
输入基线：DOC210_ROUND_PACK_2026-11-06；电磁以当前源码审计为准，Doc201 结果仅作历史证据。  
本轮目标：把方腔的进展变成可复核基线，同时补齐电磁验证链。不是重新选一套浮力算法，也不是调整法国工况参数追图。

## 0. 执行原则与本轮边界

1. 保留 R210 的源码、结果和完整命令。此前 n80 末态工程定量通过有效，不因补测而删除或改写。
2. 本轮分 H（流体归因及热边界）与 E（电磁独立验证）两条线，可并行排程；修改同一个共享文件时必须串行集成、重新构建并确认版本。
3. 先完成归因与验证，再启动法国耦合正式计算。本轮允许整理法国驱动的算子映射，不启动法国参数扫描、长时计算或 live EM—热—流耦合。
4. 不修改独立 A-phi／Laplace 论文分支，不以本轮测试为由覆盖其算法。OPHELIE 局部算子也应放在明确的局部目录或受控分支中。
5. 不叠加 density normalization、limiter、shifting、人工阻尼或新的压力方案。没有定量证据，不继续更换浮力公式。
6. 下文新给出的容差为本项目工程验收目标，不冒充文献标准。历史测试已有明确、合理且更严格的标准时保留原标准。
7. 失败是有效结果。不得为获得 PASS 修改参考值、放宽标准、挑选时刻或切换测量口径；不得把 NOT_RUN 写成 FAIL。
8. 不要求补齐旧计划全部形式化表格。只补本轮判断所需数据，避免重复工作。

## 1. 已确认事实与不能扩大的结论

### 1.1 Doc210 方腔与移动 MMS

| 项目 | 已确认结果 | 适用范围 |
|---|---|---|
| MMS，B OFF／TVF OFF，0.5 s | 内部速度 L2 约 167% | 该移动粒子配置失败 |
| MMS，B OFF／TVF ON，0.5 s | 内部速度 L2 约 167% | 单独 TVF 未解决该失败 |
| MMS，B ON／TVF ON，0.5 s | 内部 0.893%，全域 0.895% | 明确改善；尚缺 B ON／TVF OFF 分离对照 |
| 上述成功 MMS 的 TVF 位移修正 | RMS 峰值约 1.7e-11 m | 支持 B 修正路径是主要原因，不等于完整因果证明 |
| R210 方腔 n80，t*=0.5 | Nu_contact=4.50861；V*=65.25797；U*=31.43398 | 末态满足已定工程误差门槛 |
| R210 一维导热 | Nu_contact≈0.9448；Nu_fit≈0.9866；T_L2≈0.379% | 接触热流仍比解析值低约 5.5% |
| 已有 n160 | ghost OFF，但未采用最终 fluid-k 配置 | 不构成 R210 分辨率验证 |

R210 命令要素：`--dvd-exact --tvf-native --no-thermal-dirichlet-ghost --thermal-dirichlet-fluid-k --n 80 --tstar-end 0.5`。

当前 `quantitative_pass` 主要比较末态误差，未合并连续稳态窗口判断。最终配置的热流已较稳定，但仅凭全场最大速度历史不足以确认主流场稳态。

### 1.2 电磁历史证据

Doc201 的命名候选 `doc201_unified_legacy_flux_norm_off` 已有相位、梯度势组合、变电导率制造解及圆柱自感的通过项；同时仍有下列缺口：

- U1 重构电流散度指标约 0.143，内部边通量指标约 0.140，未过当时标准。
- 某 DivSigmaA 对照可使一种重构散度很小，但边通量误差仍很大，不能拼接成统一通过。
- 一个场／功率相关指标约 5.87%，其具体定义须从原脚本确认，不能只看指标名推断含义。
- 圆柱自感已有总功率误差约 0.61%、分布相关性较好等结果；仍需核查复数场误差和迭代自洽性。
- 尚无本次审查确认的 TEAM7 完整定量通过证据。

这些是历史状态，不自动等于今天代码状态。E0 先核查后续变更与已有测试，能复用就复用；不得把上述历史失败当作当前仍失败的既定事实。

## 2. 代码与数据组织

建议建立隔离工作目录 `doc211_runs/`，至少包含：

```text
doc211_runs/
  00_entry/                 # 总结、版本、执行入口
  H_causality/              # B×TVF 对照
  H_cavity/                 # n80 稳态与 n160
  H_thermal_boundary/       # 公式审计、导热和热流复核
  E_inventory/              # 当前生产路径与候选路径映射
  E_consistency/            # 复数场、电流和热源相容性
  E_cylinder/               # 解析／半解析圆柱
  E_team7/                  # 复用原 case 的边界定位
  scripts/                  # 少量统一后处理脚本
  source_snapshot/          # 本轮实际相关源码及 diff
```

- 优先扩展已有验证入口，不复制整套方腔、圆柱或 TEAM7 项目。不同 run 使用配置与独立输出目录。
- 如需工作树，从实际 R210 本地修改导出并确认；不能只检出其 HEAD，因为 R210 包含未提交改动。
- 不使用 `reset --hard` 或清理未跟踪文件来建立“干净基线”。
- 每个 run 保存源码提交、未提交 diff 哈希、二进制哈希、设备、精度、完整参数、实际开始／结束时间与终止原因。
- 同一 GPU 的精度测试尽量顺序运行；独立 GPU 可分配 H/E。不要为了并行同时改写输出或重编译运行中的二进制。

## 3. H0：拆开声学 B 修正与 TVF 开关

### 3.1 先审计，后修改

当前 `--tvf-native` 联动开启声学 `LinearCorrectionCK`。需要独立控制：

| 开关概念 | 控制对象 |
|---|---|
| acoustic_correction | 声学两半步 Inner 和 Wall 的修正模板选择 |
| tvf_enabled | 原生 KGI／TVF 位移修正流程是否执行 |
| correction_matrix_required | 根据声学或 TVF 的需求决定是否计算 B 矩阵 |

可以新增局部 CLI；名称由代码习惯决定，但任务记录必须写实际参数。不得声称旧程序原本已有这些独立参数。

关键要求：

1. B ON 必须同时落实到声学两半步、Inner 与 Wall，不能只执行矩阵计算却仍实例化 NoKernelCorrection 声学算子。
2. B OFF／TVF ON 若原生 TVF 需要 B 矩阵，可为 KGI／TVF 计算 B，但声学仍走未修正路径。“B OFF”在矩阵表中专指声学修正 OFF。
3. TVF OFF 仍正常移动物理粒子，不等于冻结 Position。
4. 不添加额外 particle shifting；TVF 系数沿用原生 0.2，不扫系数。
5. 保留旧 `--tvf-native` 的兼容行为，若与新显式参数冲突，程序应报错或清楚记录解析优先级，不能静默切换。
6. 导出实际模板／路径标识、矩阵执行次数、TVF 执行次数、位移 RMS／最大值及其相对 dp 的幅度。
7. MMS 体力仍按当前空间坐标独立解析计算，不覆盖速度；审计 ForcePrior 生命周期，防止覆盖外力或跨步累积。

## 4. H1：最小因果矩阵

### 4.1 移动 MMS

使用 Doc210 完全相同的制造解、参数、初态与壁面条件，n80，先 0.05 s，再 0.5 s。物理时间与方腔 t* 不混淆。

| ID | 声学 B 修正 | TVF | 执行要求 |
|---|---:|---:|---|
| H1-00 | OFF | OFF | 旧结果满足复用条件则不重跑 |
| H1-01 | OFF | ON | 旧结果满足复用条件则不重跑 |
| H1-10 | ON | OFF | 必做的新关键对照 |
| H1-11 | ON | ON | 优先复用 Doc210；独立开关改造后先短跑确认路径等价 |

复用条件：实际算子、参数、精度、设备执行路径、初态、终时与误差定义一致，且有可追溯结果；仅名称相同不够。改造触及共同执行逻辑时，成功配置至少补一次正式复现。

MMS 工程目标：0.5 s 内部 L2≤5%、全域 L2≤10%，不发散、不冻结物理运动。输出误差随时间曲线，不能只保存末点。

归因判定：

- H1-10 与 H1-11 都过且明显优于未修正组：可写“声学修正足以消除本组主要失败，TVF 非必要条件”。
- 两者误差差异不超过 max(0.2 个百分点，较大误差的 10%)：可记为本测试中影响较小，不能写成所有工况 TVF 无用。
- 只有 H1-11 过：检查 TVF／KGI 相容性、位移和粒子质量，报告组合依赖，不强行归因于 B。
- 若 H1-10 异常：先检查新开关有没有漏掉矩阵更新或错误声学模板；只允许一个由证据支持的修复，再复测。

### 4.2 方腔对应对照

仅在 MMS 完成后运行最终热边界设置下的 B ON／TVF OFF n80 方腔，与 B ON／TVF ON 比较。

不要求重新跑全部四组合长时方腔。保持实际 Ra、Pr、温差、声速、黏性、时间推进与测量方法一致。

比较均以稳态窗口统计为准：中心线 U/V、峰值位置、Nu、固定网格速度与温度、近壁粒子分布。若 TVF 对法国后续高黏度工况的必要性尚未验证，在结论中保留该适用范围。

## 5. H2：最终方腔配置的稳态与分辨率

### 5.1 物理条件锁定

使用严格 DVD 方腔：Ra=1e5，Pr=0.71，L=0.0386 m；左热右冷，Th=303 K、Tc=283 K；上下绝热；四壁无滑移；初始 293 K。沿用严格模式的物性计算，不恢复旧近似 Ra／Pr 默认值。

参考值：Nu=4.519；水平中心线 Vmax*=68.59，位置 x/L=0.066；竖直中心线 Umax*=34.73，位置 y/L=0.855。

参考量无量纲速度使用 alpha/L；t*=t alpha/L²。输出实际 Ra／Pr 和实际终时，不能只写请求值。

### 5.2 稳态测量

1. n80 从 t*=0.5 附近开始检查；没有兼容 restart 则重跑，禁止伪造续算。
2. 保存固定空间采样网格的速度／温度，配合中心线数据；不得按粒子序号直接比较空间稳态。
3. 每个 Δt*=0.02 的窗口内至少有 10 个样本，连续检查三个窗口。原先每窗口一个样本不足以计算窗口平均和波动。
4. 采用固定内部采样点和固定插值方法；记录邻域覆盖率。无法覆盖处标记缺失，不补零，不沿壁强制制造速度回流。
5. 相邻窗口的空间速度均值变化 L2≤1%，温度均值变化 L2/ΔT≤0.5%，Nu 均值变化≤0.5%；速度归一化采用固定参考量或明确的全场 RMS 下限，防止零速区除零。
6. 同时报告窗口内波动、中心线峰值均值／范围和全场最大速度。均值稳定但波动仍大时，只可称统计稳定，不能直接称稳态解。
7. 无量纲时间先到 0.56 检查；必要时延至 0.7，再到最多 1.0。仍未满足则停止并定位，不无限延长。

量化目标沿用：Nu_contact 相对误差≤5%，中心线 U/V 各≤10%，峰值位置绝对误差≤0.05L。通过记录需同时注明稳态状态，不能只输出一个综合 PASS。

### 5.3 同配置 n160

- 必做最终 R210 热边界设置的 n160；与 n80 对齐稳态标准，不仅对齐末时刻。
- 全部测量使用相同无量纲坐标和空间口径，检查中心线 IDW 支持半径敏感性。
- 报告两分辨率各自相对参考的误差和相互差异。两分辨率不能给出可靠收敛阶，也不称严格网格无关。
- 两档都达工程目标、结果无显著恶化，方可写“本项目两档分辨率验证通过”。
- 若 n160 某项越界或明显变差：优先分开检查稳态、空间采样、热边界层和粒子质量；本轮不自动追加 n320。

## 6. H3：热边界公式和一维导热收尾

### 6.1 必须解释的公式

把实际热接触离散完整写出，标明：核梯度权重、负号、体积因子、距离与正则项、流体 k、壁粒子 k、ghost 温度跳跃以及物理壁面位置。

重点区分两类边界：

- 指定温度的数学边界：虚粒子参数用于实现边界，不代表一块具有真实导热热阻的固体。
- 真实共轭导热界面：两侧 k、距离和接触热阻都具有物理含义，应依据串联热阻处理。

因此不能把 `harmonic(kf,kw)≈2kf` 一概称为错误，也不能只因改成 kf 后接近 DVD 就称其普遍正确。ghost 倍增和导热系数处理必须连同几何一起解释，避免重复补偿。

### 6.2 最小导热测试

1. 在 R210 固定粒子一维导热上复用 n80，补 n160；左右定温、上下绝热、无流动，初态与终时明确。
2. 解析稳态温度线性，Nu=1。给出体积加权温度误差、冷热端接触通量、独立梯度通量和第一层粒子几何。
3. 温度误差定义：sqrt[ΣVi(Ti−Texact)²/(ΔT²ΣVi)]，不得额外除粒子数。
4. 独立梯度复核采用有固定壁面值约束的局部拟合，比较 2dp/3dp/4dp 支持宽度，不选择最接近参考的一档作为唯一结果。
5. 检查不同层粒子到物理壁面的距离、积分面积、角点归属、两壁功率符号和单位。
6. 工程目标：T_L2≤2%，Nu_contact 误差≤5%，独立梯度 Nu 误差≤3%。当前 5.5% 如未解决，准确记为未过，不能改门槛到 6%。

### 6.3 修改与停止条件

先诊断再改代码。仅在发现明确的距离、重复系数、热阻或积分错误时实施一次针对性修复；修复后回归一维导热和 n80 方腔，确认不是用一种误差抵消另一种误差。

若没有公式错误而主要是有限分辨率：保留 R210、报告加密趋势和适用范围，不继续试验 face／rmin／Shepard 参数组合。允许“速度验证通过、热边界部分通过”的分项结论，不强求整体标签。

附加检查：`q_imbalance` 小只说明两端瞬时热流接近平衡。应独立检查 ΔΣmcpT 与边界热输入积分的差，并核查内传热是否全局守恒；变物性阶段需改用一致焓定义。本轮常物性即可，不扩展变 cp。

## 7. E0：冻结统一电磁路径与证据清单

先完成当前源码映射，输出 `EM_IMPLEMENTATION_MAP.md`：

| 必填项目 | 要写清什么 |
|---|---|
| 时谐约定 | exp(+iωt) 或 exp(−iωt)，峰值或 RMS，相量单位 |
| 激励 | A_ext／线圈电流来源、频率、几何、实际电流幅值 |
| 电势方程 | 实部／虚部各自实际算子、右端、约束和电流边界 |
| 电流 | 求解使用的边通量、可视化重构 J、是否为同一离散定义 |
| JouleHeat | 依据哪一种 J 或 E、平均功率系数、负值处理、归一化位置 |
| 自感 | A_ind 计算、迭代次序、松弛、停止条件和最终状态一致性 |
| 生产路径 | 法国驱动实际调用的类、模板和参数，与候选验证路径的差异 |

不得把实部使用一种离散、虚部使用另一种离散的默认实现，自动继承统一候选的测试结果。

冻结一个候选；历史实现仅作为对照，不并列大规模扫描。全局总功率归一化在基础验证中关闭。先记录未归一化误差，再单独检查归一化操作是否只改变整体幅度。

源文件或测试数据缺失时，先查现有仓库与已归档资料；确实缺失则列具体项，不用“缺数据”泛泛阻断整条线。

## 8. E1：复用基础测试并定位电流／热源一致性

### 8.1 基础烟测

复用已有 U1–U5：相位旋转、纯梯度势抵消、A0+gradχ、空间变 σ、非零复数电势。

- 保留原解析／制造解、边界与误差归一化。当前统一候选未变且证据足够的测试不全量重跑。
- 公共算子发生修改后，基础烟测必须重跑，防止修散度破坏复数相位或可积梯度抵消。
- 纯梯度抵消参考电流为零，误差分母使用预先定义的激励尺度，不能用近零数归一化。

### 8.2 同一快照同时计算的四组量

1. 求解器残差：实际离散电势方程及对应 RHS 归一化。
2. 边通量守恒：逐粒子通量和、内部边配对反对称性、边界法向电流。
3. 重构 J：与解析 J 的复数场误差，以及独立散度、边界分布。不要把后处理梯度算子当成求解器本身。
4. 焦耳热：局部 Q 分布、体积分总功率、与当前相量约定相符的 J／E 功率关系。

峰值相量下常用平均功率为 Q=0.5 Re(J·E*)；若使用 RMS 则不带 0.5。实数各向同性 σ>0 时可与 |J|²/(2σ) 对照；不得跨不同定义的重构 J、边功率直接假定逐点相等。

明确哪些关系应严格离散满足，哪些只能随分辨率趋近。只有当输入场和离散定义相容时，才将功率恒等式作为严格测试。

必须分区输出内部、近壁和全域误差，至少给 RMS／最大值／分位数及空间图；不能用外围大量近零粒子稀释内部误差。散度尺度采用固定 Jchar/L，并记录 Jchar 公式。

### 8.3 有限分支定位

| 观察 | 优先动作 |
|---|---|
| 原生边通量已不守恒 | 检查边配对、边界通量、复数 RHS 和线性收敛 |
| 边守恒通过，但重构 J 误差大 | 检查重构的一致性、局部条件数和近壁支撑；优先局部重构修复 |
| J 合理，但 Q 分布／功率错 | 检查 0.5、σ、共轭、单位、体积及归一化位置 |
| 仅近壁失败 | 分开检查法向／切向误差和边界覆盖，不同时重写体积算子 |

每条失败先只做一个有证据的修复，不建立多种算子排列组合。不得通过投影后处理生成漂亮 J，然后仍用不一致的旧 J 计算 Q 而宣称整体修复。

验收：保留历史测试容差并注明来源。若旧指标定义不合理，必须同时输出原指标与修正指标，解释变化，不追溯改写旧 PASS／FAIL。

## 9. E2：圆柱趋肤／自感定量验证

### 9.1 复用现有 U6

优先采用已有轴对称圆柱解析／半解析测试，不新建平行项目。先把参考模型写清：无限长理想圆柱还是有限长模型，外部场方向、边界、σ、μ、ω、圆柱半径和采样区域。

理想无限长解不能直接与有限高度圆柱端部全域误差混算。允许使用中部剖面作有限长度近似，但必须说明端部影响并用长度或区域敏感性检查支持。

注意 μ 表示磁导率，与流体动力黏度符号区分。趋肤深度 δ=sqrt[2/(ωμσ)]。

### 9.2 最小测试矩阵

- 主验证点：已有 δ/R≈1，保持其激励与物性，做两档空间分辨率；复用有效旧档。
- 补充点：δ/R≈3，作为较弱趋肤对照；沿用已选离散，不调系数。
- δ/R≈0.5 仅在前两点通过且边界层得到足够分辨时追加，本轮不强制。
- 不要求用法国玻璃参数替代标准测试参数，也不在本轮启动 σ(T) 热反馈。

### 9.3 输出与目标

输出复数 J 或相应解析可比场的径向实／虚部、幅值和相位；低幅值区域的相位误差应掩膜并报告阈值。还需 Q(r)、未归一化总功率、最大局部误差、分辨率和 δ/dp。

工程目标：可比区域复数场 L2≤5%，总功率误差≤3%；分辨率加密不出现无解释恶化。相关系数只能辅助，不能替代幅值误差。

自感迭代在同一最终状态重算：φ、J、A_ind、Q。保存更新前后相对变化及最终方程残差，避免输出不同迭代轮的混合场。

增加一次容差收紧对照：线性与外迭代收紧后，关键场／功率变化目标≤0.5%，或显著小于离散误差。不能仅凭迭代增量小证明方程求解正确。

## 10. E3：TEAM7 边界问题定位

### 10.1 本轮定位，而非孤注一掷

TEAM7 继续使用原 case；现有标准几何、线圈、测量点和实验数据应整理进同一入口。优先完成一个已有频率的定量对照，通常先原有 50 Hz；其他频率待其解释清楚再扩展。

必须核对并记录原始基准来源：板的开孔和厚度、线圈位置／截面／匝数、激励安匝峰值或 RMS、坐标、探针位置、所比较磁场分量、实验／参考曲线和单位。禁止凭记忆填标准参数。

### 10.2 一次只改变一个边界因素

1. 用当前配置复现异常，并画空间误差，区分导体表面、电势约束、自感迭代或外域截断。
2. 若磁场采用 Green 函数／Biot–Savart 的开放域积分，不添加无意义的空气盒敏感性；若是有限外域 PDE 才做外域尺寸对照。
3. 若怀疑导体绝缘边界，测 n·J 与局部守恒，固定其他条件做一种边界修复。
4. 若怀疑线圈激励，先在无导体状态对照磁场与几何／积分收敛，避免把激励误差归给板边界。
5. 若 σ 被错误扩展到空气或支撑区，优先修正导电域定义；不要用总功率归一化掩盖。

全程保留未归一化结果。能对齐实验测点时输出逐点表和曲线，误差分母使用参考峰值等明确尺度，不对过零点计算失真的逐点相对误差。

目标：复用已有 TEAM7 验收标准；若历史未设标准，运行前登记建议场曲线归一化 L2≤5%、峰值误差≤10%，并说明实验误差与数据来源。这是工程目标，不能替代原基准的定义。

若缺测点数据，只能完成场与边界一致性检查，状态写“基准定量未完成”；不能以残差小替代 TEAM7 实验验证。

## 11. 排程与停止规则

### 11.1 建议执行顺序

| 阶段 | H 线 | E 线 |
|---|---|---|
| 第一批 | 冻结 R210，独立开关，H1-10 MMS | E0 路径审计、复用清单、已有基础烟测 |
| 第二批 | n80 稳态与 TVF 对照 | E1 统一快照的电流／热源一致性 |
| 第三批 | 一维导热审计与 n160 | E2 圆柱主点、加密和迭代收紧 |
| 第四批 | 最终配置 n160 方腔 | E3 TEAM7 单频率边界定位 |
| 汇总 | 归因和验证分项结论 | 当前候选可用范围与剩余阻断项 |

先短测再长跑；正式运行前确认输出监测已生效，避免跑完后发现只有末点。

### 11.2 防止再次循环

- 某失败首次出现：检查配置与指标定义；找到证据后最多做一个针对性修改。
- 同一假说经两次有控制的测试仍无支持：停止该方向，写清已排除内容，不继续换近义开关。
- 不因缺少漂亮的综合 PASS 而反复调系数；也不因一个末点 PASS 就关闭不一致性问题。
- 本轮达不到全部目标时，必须产出具体的“哪条关系不满足／在哪个区域／误差多大／已排除什么”，不能只写“仍需优化”。

## 12. 收口与下一阶段准入

分别给出以下结论，不把所有模块压成一个标签：

1. **MMS 归因**：是否确认声学 B 修正为主要因素，TVF 是否必要。
2. **方腔速度**：末态、稳态、两档分辨率各自是否通过。
3. **热边界**：公式解释、导热解析对照和方腔热流是否一致；仍有多少误差。
4. **电磁基础离散**：边通量、重构 J、复数约定、Q 是否相容。
5. **电磁自感**：圆柱解析对照及最终状态是否通过。
6. **TEAM7**：定量通过／边界已定位／基准数据不足，明确是哪一种。

法国固定温度场短回归的前提是 H 线证明流体配置可信；它不必等待 TEAM7 全部收口，但依据用户本轮优先级，本轮先完成验证汇总再启动。

法国采用冻结电磁热源推进温度演化前，还需确认 E1/E2 对该热源路径有支撑，并审计法国热边界；不能直接复制方腔 Dirichlet 处理到 Robin／辐射边界。

TEAM7 尚未通过不自动证明轴对称法国工况不可用；须判断失败是否涉及共用电流／热源算法。若是共性缺陷，阻断相应电磁耦合；若仅限 TEAM7 特有几何或外域，记录限制后单独决策，不泛化结论。

本轮不自动授予“法国自然对流已验证”或“电磁全部正确”。

## 13. 必交付文件与图

### 13.1 入口

- `ROUND211_REVIEW.md`：先给实际结论、证据和剩余问题，不能复制计划充当完成总结。
- `RUN_INDEX.csv`：每个 run 的实际配置、版本、终时、终止原因、指标、分项状态、结果路径。
- `REPRODUCE_COMMANDS.md`：真实可执行命令；区分建议命令与已执行命令。
- `SOURCE_AND_BINARY_PROVENANCE.json`：源码、diff、二进制、环境与测试对应关系。

### 13.2 H 线

- `B_TVF_CAUSALITY.csv` 与结论：四组 MMS 对照、复用来源、实际声学路径及 TVF 统计。
- `CAVITY_STEADY_WINDOWS.csv`：窗口内样本数、均值、波动、相邻变化、覆盖率。
- `CAVITY_N80_N160.csv`：同配置参考误差与分辨率差异。
- `THERMAL_BOUNDARY_FORMULATION.md`：实际公式、几何解释、一维导热与修复依据。
- 图：MMS 误差历史；方腔中心线 U/V 对照；Nu 历史；导热温度与通量对照。

### 13.3 E 线

- `EM_IMPLEMENTATION_MAP.md` 与 `EM_VALIDATION_LEDGER.csv`。
- 电流守恒／重构／热源一致性表和空间误差图。
- 圆柱复数场与 Q 的径向曲线、分辨率表、迭代收紧对照。
- TEAM7 标准参数来源、测点对照、边界诊断与结果限制。

图必须由打包数据和脚本重现。参考曲线、计算曲线、误差定义、时间／分辨率／单位明确；不只提交截图或终端 PASS 文本。

### 13.4 打包验收

保留必要源码、diff、参数、原始 CSV/JSON、关键日志和少量图；无需打包 build、缓存、全部临时脚本或海量 VTP。

README 必须回答：

1. 为什么 Doc210 改善，哪项证据补全了因果？
2. 最终方腔配置在稳态和加密后是否仍可信？
3. 热边界 5.5% 导热缺口是否解决，未解决如何限定结论？
4. 当前法国电磁路径与已验证候选是否一致？
5. 电磁剩余问题属于求解、重构、边界还是热源转换？
6. 下一轮可以启动哪一种法国短回归，还有哪个具体阻断项？

> 本轮成功标准是获得可复核、能指导下一步的结论；不是新增一批开关，也不是让每张表最后一列都变成 PASS。


<!-- ===== END 211_OPHELIE_CAVITY_CAUSAL_CLOSURE_AND_EM_VALIDATION_CURSOR_PLAN_2026-11-06.md ===== -->

---

<!-- ===== BEGIN 212_OPHELIE_TARGETED_HYDRO_EM_REPAIR_CURSOR_PLAN_2026-11-06.md ===== -->

<a id="212-ophelie-targeted-hydro-em-repair-cursor-plan-2026-11-06-md"></a>

# 212_OPHELIE_TARGETED_HYDRO_EM_REPAIR_CURSOR_PLAN_2026-11-06.md

# Doc212：方腔近壁稳定性与电磁一致性定向修复任务书

日期：2026-11-06  
基线：`DOC211_ROUND_PACK_2026-11-06.tar.gz`；保留 Doc210／211 原始结果。  
执行对象：Cursor。本文是下一轮实施任务，不是已完成结果报告。

## 0. 本轮需要交付什么

本轮从“重复测试并登记失败”转向“定位具体不一致并实施有依据的修复”。限定四个任务：

1. **H：解释并修复 n160 方腔速度波动与热流偏差。** 保留 n80 已通过基线，不再换浮力公式。
2. **E-A：统一生产电磁实／虚部离散。** 首先解决相位不变性与电流守恒，不再用另一套候选的 PASS 背书默认实现。
3. **E-B：拆开圆柱功率误差。** 分离解析模型、几何积分、场误差，并补完整复数场比较。
4. **E-C：在简单几何上验证绝缘边界。** 通过后再迁回 TEAM7，停止直接在复杂薄板上盲调。

H 与 E 可在独立设备／隔离工作树并行；E 内按 E-A→E-B／E-C 推进，几何积分离线分析可提前。最终只允许一个经过验证的 OPHELIE 电磁候选进入集成，不保留多套“默认生产配置”并列使用。

本轮默认不运行法国长时计算、参数拟合、文献流线扫描或新的变物性全耦合。法国部分先完成调用路径映射，待阶段准入通过后只允许受控短回归。

## 1. 必须保留与修正的已有结论

### 1.1 流体事实

| 项目 | Doc211 证据 | 本轮解释 |
|---|---|---|
| 移动 MMS，B ON／TVF OFF | 0.5 s 内部 L2≈0.893% | B 修正足以解决此制造解的主要失败 |
| 方腔，B ON／TVF OFF | 短时发散 | 此方腔不能直接取消原生 TVF 配套 |
| 方腔 n80，B＋TVF，R210 热边界 | t*=0.70：Nu≈4.4506；V*≈66.062；U*≈31.660；窗口速度变化≈0.50% | 该配置工程通过，有明确适用范围 |
| 方腔 n160，同热边界 | t*=1：Nu≈4.8178（+6.61%）；窗口速度变化≈3.3% | 热流与稳定性未过，不能选择 n80 来注销问题 |
| 一维导热 n80／n160 | 接触 Nu 约−5.5%／−4.6% | 与对流方腔偏差方向不同，不支持“接触算子总偏热” |

不重新跑 MMS 全矩阵。仅共同声学实现被修改时，复测成功 MMS 作为回归。

### 1.2 电磁事实

- 默认生产路径存在实虚 RHS 不一致；相位旋转 J≈9%、Q≈11% 误差。
- 统一候选的相位与部分制造解通过，但不代表其所有守恒检查都通过。
- 圆柱功率误差约 6%～9%，不能直接归因于电磁算法或有限高度。
- 圆柱中部约 1% 指标由幅值 CSV 可复算；现有数据不能据此称完整复数矢量场误差约 1%。
- TEAM7 的单次 RHS Neumann＋LevelSet 试验失败，且整板 7700 粒子都进入边界掩膜。排除的是这套实现有效性，不能排除绝缘边界是原因。

### 1.3 法国短跑结论修正

Doc211 法国 5 s 日志中 `kernel_correction=0`，不能视作 R210 声学修正迁移成功；`validation_passed=0`，只能称部分管线完成烟测。

`dT≈35 K` 是 Tmax−Tmin，不是整体升温。末态平均温度约 1468.5 K，需与真实初态核对；热源约 49 kW、末态损失约 153 kW，属于显著非平衡冷却。后续报告分别写 Tmax、Tmin、Tmean、空间温差及平均温变。

## 2. P0：先修复源码—二进制—结果对应关系

Doc211 包内方腔源快照未找到稳态监测实现，哈希清单也没有方腔二进制。先确认是漏打包还是使用了不同源码；不得自行推断数值数据造假或全部无效。

必须完成：

1. 找到实际生成 `steady_samples` 和 `windows.csv` 的代码／脚本；保存源码及调用方式。
2. 对实际方腔、圆柱、TEAM7、法国驱动分别记录二进制 SHA256、源码树与未提交 diff；单个 HEAD 不足以标识本项目。
3. 检查是否存在旧二进制、错误工作目录、同名入口编译目标错配。
4. 只在真实源码缺失／无法复现时重建最小监测，不重复编写已有工具。
5. 未识别 CLI 必须报错，运行横幅导出实际解析后的配置；不要默默忽略拼错参数。
6. 结果字段分开：完成运行、数值稳定、解析／基准对照、稳态、守恒、完整验证。综合 smoke PASS 不覆盖分项失败。

建议使用 `doc212_runs/`：`00_entry/`、`H_wall/`、`E_unified/`、`E_power/`、`E_boundary/`、`integration/`、`scripts/`、`source_snapshot/`。不复制整套 case；扩展现有验证入口。

禁止清理用户未提交工作。创建隔离工作树须包含 R210／211 实际本地改动，而非只 checkout 历史 HEAD。独立 A-phi／Laplace 研究分支不修改、不合并。

## 3. H0：锁定方腔基线与诊断窗口

保持严格 DVD 模式：Ra=1e5，Pr=0.71，L=0.0386 m，Th=303 K，Tc=283 K，初温293 K，四壁无滑移、上下绝热。物性、核函数、声速、热时间步和积分算法沿用 R210。

参考 Nu=4.519；中心线 Vmax*=68.59@x/L=0.066，Umax*=34.73@y/L=0.855；速度尺度 alpha/L，时间 t*=t alpha/L²。

保持 B＋原生 TVF＋ghost OFF＋fluid-k。TVF 关闭的发散只作历史对照，不重跑长时。不能额外加 shifting、density normalization、人工阻尼或改 β。

优先复用可信 restart；无兼容 restart 时重跑并标明成本，不伪造续算。n80/n160 采样同一无量纲晚期区间，稳定性比较还要覆盖各自已知平台窗口。

## 4. H1：一次运行收集足够定位的数据

### 4.1 分区与空间观测

按粒子到物理壁面的真实距离分区：0～dp、dp～2dp、2dp～4dp、内部；角点另标。分区由当前 Position 决定，保留粒子唯一 ID，避免排序后错配。

固定空间网格速度／温度仍用于稳态；粒子数据用于定位局部原因，不能按粒子序号直接比较空间稳态。

各区每个采样时刻输出 RMS、P95、最大值；重要异常输出粒子 ID 与位置：

| 诊断组 | 必需量 |
|---|---|
| 运动与密度 | u、rho/rho0、EOS 压力、到壁距离/dp、越界／穿透独立统计 |
| 粒子质量 | 邻居数、最近间距/dp、核矩特征值或奇异值、B 条件数和最大范数 |
| TVF | 纯修正位移/dp，与物理位移的比值；修正前后近壁粒子分布 |
| 动量 | 压力／声学、黏性、浮力各自加速度，按实际使用字段取值 |
| 温度与热流 | T、近壁温度跳跃、接触热流、各距离层对总热流贡献 |

B 条件数只作诊断；未定位前不启用截断或钳制。热流分别按热壁／冷壁实际接触对象累计，不能仅靠宽带位置近似归属角点。

采样：连续三个 Δt*=0.02 窗口，每窗至少10点；仅保存少量完整快照，例如异常前、异常峰、恢复后。额外高频诊断可只对选定异常粒子和一个短窗口执行，避免巨大输出。

### 4.2 必做的调用相容性审计

画出实际一整个对流步的调用顺序，逐一核对：

- Position／Displacement 更新时间，声学两半步，KGI／B 更新与 ForcePrior 生命周期。
- TVF 是否只改变位置、是否需要与温度输运配套；按本仓库原生定义判断，不凭“TVF≈shifting”推断。
- 更新位置后，热模块使用的邻居索引、核距离与核梯度分别属于哪个时刻；当前坐标与旧邻域是否混用。
- 双时间步框架有意冻结声学邻域不自动是错误；只检查是否违反所选原生流程及其几何假设。
- 用少量粒子暴力搜索核支撑，检查真正漏邻居；不要看到声学子步未重建就强制每子步重建。
- 排序后的温度、黏性力、浮力历史和 B 等字段是否同步排列；诊断是否读取 GPU 最新数据。
- 热推进实际子步累计时间与流体对流步是否一致；是否触及子步上限而静默漏推进。

## 5. H2：按证据选择一个修复，不扫矩阵

| 证据 | 优先修复方向 | 验证方式 |
|---|---|---|
| 热算子使用的位置／邻域不相容且存在实际漏邻居 | 调整必要的热前邻域更新或几何刷新 | 同快照热算子前后对照＋导热＋方腔回归 |
| B 极端值与局部速度爆点同步 | 核查近壁支撑、矩阵构造和壁面贡献 | 固定快照检查矩再现；修几何／构造后回归，不先钳制 B |
| TVF 导致近壁堆积／壁面方向错误 | 与成熟原生闭腔流程核对其适用范围和墙贡献 | 只修发现的实现差异，不换算法／扫系数 |
| 压力或密度异常先于速度噪声 | 核查 EOS、Vol=m/rho、壁声学项及历史外力 | 阶段快照与原生实现对照；必要时短窗口减半时间步 |
| 速度较稳但热流由少数过近粒子主导 | 检查接触几何／积分及温度跳跃一致性 | 独立热阻／壁距复核，修复后导热和方腔共同验收 |
| 高频变化与时间步有关，空间异常不明显 | 只做一次 dt/2 短窗口对照 | 同物理时间对齐，不能同时改声速与核参数 |

若当前快照不足以区分原因，先增加一次针对性的测量，不直接同时改多个算子。报告必须给“异常首先出现在哪个阶段、哪个区域、哪个量”，不能只说近壁可能有问题。

如果怀疑单精度消减或矩阵精度，只针对已定位的局部算子做 CPU 双精度同快照复算；不要未经证据把整个工程迁成双精度。

### H 验收

每个有依据的修复先过一维导热和 n80，随后过 n160；声学公共逻辑变更时补成功 MMS。基准物理参数不变。

- Nu_contact 误差≤5%；中心线 U/V 误差各≤10%；峰位置误差≤0.05L。
- 连续三个窗口的速度均值空间 L2 变化≤1%，温度变化/ΔT≤0.5%，Nu均值变化≤0.5%；同时报告窗内波动，不能仅报告平均抵消。
- 若 n160 不过但有清晰修复方向，继续完成对应修复与回归；“本轮停止延时”不等于“停止解决”。
- 不自动加 n320，不因 n80 过就授予分辨率通过；两个分辨率也不能声称严格收敛阶。

## 6. E-A0：统一实虚部装配，而不是只切一个标志

以当前 `exp(+iωt)` 峰值相量约定为准：

```text
E = -grad(phi) - i*omega*A
Er = -grad(phi_r) + omega*Ai
Ei = -grad(phi_i) - omega*Ar
Q = 0.5*Re(J dot conj(E))
```

对实、虚部使用同一套：电导权重、LHS、RHS 激励离散、边界、规范约束、归一化策略和求解容差。允许符号因实虚分裂不同，不允许一个链用边源而另一个链用不相容散度。

复用现有统一候选作为起点，不复制第二套“新统一候选”。对照旧生产路径保留原始失败证据。

实现要求：

1. 优先抽取共享的复数分量装配／求解函数，由 active A 字段与符号参数区分实虚；避免两份函数继续漂移。
2. 两个组件的真实实例与参数在运行时输出，不能只输出 `unified=true`。
3. 弄清 gauge penalty、零均值、固定参考点如何消除常数零空间。特别检查罚项是否改变物理解，而非仅规范。
4. 规范约束采用同一方式。若临时保留小罚项，输出其量纲、相对算子尺度及对守恒残差的贡献，做一次收紧对照；不能盲用1或1e-8。
5. 保持零均值源项投影的前后统计；若投影移除了显著物理源，先查兼容条件，不用投影掩盖装配错误。
6. 基础验收关闭功率标定、场幅缩放和后处理 J 投影； solver-local 缩放如为数值条件化必须严格逆变换并记录。
7. 不在未通过前直接覆盖法国默认路径；通过后以共享生产 profile 接入，不让测试入口独享特殊补丁。

## 7. E-A1：用离散恒等式定位守恒失败

在同一小规模现有测试上保存 LHS、RHS、规范项、边通量、边源、重构 J 与 Q。同一个快照、同一个 A、同一个 phi，避免跨轮混合。

首先写清：边通量和与电势线性方程残差应有什么精确关系。Vi、边电导单位、正负号、边界源项和 gauge 项全部纳入。

特别注意：带 penalty 的方程收敛后，裸边通量和可能剩下 penalty*phi 或源项投影量。不能一边保留这些项，一边要求裸通量和机械等于机器零；也不能把物理泄漏自动称为规范误差而不量化。

先验证“离散恒等式自身相等”，再判断物理电流守恒是否足够好：

- 若恒等式不成立：装配／审计／边权不一致，修这一处。
- 若恒等式成立但 RHS 与边源不兼容：修复 RHS，不先改 J 重构。
- 若守恒通量可信但 J 重构不准：单独检查重构矩阵与近壁支撑，保留原通量作依据。
- 若 Q 用另一套 J／E：先确定唯一热源路径；不同离散功率表达式的比较要注明预期等价范围。

对内部／近壁／全域输出误差；避免大量近零区域稀释。零电流测试的误差分母采用固定激励 Jchar，不用近零计算结果自身归一化。

## 8. E-A2：最小回归与验收

沿用现有 U2–U5 和 A1 入口，不另建重复 case：

| 测试 | 必须比较 |
|---|---|
| 相位旋转 A→iA | J→iJ；Q不变；phi按同规范变换 |
| 纯梯度 A=gradχ | 对应 phi 抵消后 J≈0，包含边界而非只测内部 |
| A0+gradχ | 物理电流／热源不变，输入变换遵循固定时谐规范 |
| 变σ与非零复phi | 两组件在相同空间变系数下仍一致 |
| 非零电流测试 | 边恒等式、物理守恒、重构场与Q分别验收 |

项目目标：相位旋转 J 相对误差≤1%、Q≤2%，其余保留原合理标准。若先前候选已达1e-6量级，新实现不能把明显退化合理化为“仍在1%内”。

边恒等式误差应与线性残差及实际浮点精度同量级，不能预定机器零；物理边守恒沿用已定义的≤1%工程目标，并单独报告规范／边界贡献。

最终选择必须是一套同时通过相位与守恒的路径。禁止用“候选相位通过＋旧生产重构散度通过”拼成整体通过。

## 9. E-B：圆柱误差拆分，先离线再改算法

### 9.1 已复算的定位线索

以下数据由 Doc211 中部幅值 CSV 复算，按其轴向重复布局外推全高，供定位而非正式验收：

| 工况 | 总功率相对连续解析误差 | 同粒子解析幅值积分的几何偏差 |
|---|---:|---:|
| δ/R=1，dp=0.025 | +9.26% | 约+5.95% |
| δ/R=1，dp=0.020 | +6.53% | 约+3.51% |
| δ/R=3，dp=0.025 | +5.98% | 约+5.94% |

弱趋肤功率偏差几乎可以由几何积分解释，但正式结论必须导出全域数据验证，不将中部外推当作已测全域。

### 9.2 正式计算三个功率

在实际全体粒子位置 xi、体积 Vi 上计算同一解析参考：

```text
P_cont   = 连续理想圆柱上对解析 Q 的高精度积分
P_sample = sum_i Vi * Q_exact(xi)
P_num    = sum_i Vi * Q_num(xi)

P_num-P_cont = (P_sample-P_cont) + (P_num-P_sample)
```

分别报告：几何积分偏差、同采样场功率偏差、总偏差。恒等式须闭合。

如果数值求解是有限高度而参考为无限长，最后一项仍含模型差异，不能全部命名为纯数值误差。分中段／端部报告，必要时只增加一次较长圆柱对照；不预先声称有限高已经解释差值。

同步输出 ΣVi、πR²H、ΣVi r² 及其连续精确矩。弱趋肤下 Q近似正比r²，有助于解释为什么体积误差不足以单独说明功率误差。

不得把 P_num 乘 P_cont/P_sample 后作为“原始功率通过”。可作误差分析辅助，但原始数据和总误差必须保留。若确需修正边界粒子体积，必须一致用于源积分、算子和场重构并回归，不只改功率表。

### 9.3 补完整复数场

全域 CSV 至少导出 xi、Vi、σ、Jr/Ji 三分量、Er/Ei、Q，以及同位置解析复数场。

完整复数向量误差：

```text
eJ = sqrt[sum Vi*(|Jr-Jr_ref|²+|Ji-Ji_ref|²)
          /sum Vi*(|Jr_ref|²+|Ji_ref|²)]
```

另报幅值误差、强场区相位误差、径向／轴向杂散分量。低幅值区域相位掩膜阈值事先固定并输出。禁止拟合全局相位后隐藏符号错误；若展示相位对齐诊断，须同时保留未对齐主结果。

保持主点δ/R=1两分辨率及弱点δ/R=3，复用有效既有结果；不扩大频率矩阵。共享电磁装配修复后至少重跑上述组。

工程验收：可比区域完整复数场 L2≤5%；同采样功率偏差目标≤3%；连续总功率≤3%单独标记。只有后者也达标才写“连续功率验证通过”。场与总功率不混为一个 PASS。

## 10. E-C0：简单绝缘边界测试

新建一个轻量边界验证模式或小 case 是合理的，因其物理目标与现有圆柱不同；优先放入已有制造解验证入口，不复制 TEAM7。

### 10.1 主测试：封闭长方体、常量矢势、零电流解析解

采用三维长方体导体，常σ、给定空间常量复数 A0、无自感、无真实线圈积分。解析：

```text
phi_exact(x) = -i*omega*A0 dot (x-xc) + C
E_exact = 0
J_exact = 0
Q_exact = 0
```

固定同一规范常数 C。物理边界条件为绝缘 n·J=0，即 dphi/dn=-iωA0·n。

必须通过待验证的自然／Neumann 边界路径求解，不能在每个粒子直接写入解析 phi 或用全边界 Dirichlet 替代后宣称绝缘边界通过。解析注入可作算子单测，但与求解验证分开。

先取 A 沿一个坐标方向，再做一次整体90°相位旋转。法向直接由解析面给出，用于先验证公式和装配；随后以数值几何法向替换，分离几何误差。

### 10.2 面贡献必须正确表示

薄体的同一粒子可能同时受两个相对表面的影响，不能只存“最近面一个法向”然后将整粒子当单侧边界。

- 明确是面片积分、虚粒子、弱自然边界或其他离散表示。
- 如果现有边通量形式已隐式包含绝缘条件，证明其离散机制，再判断显式 RHS 是否重复施加。
- 如采用表面到体积映射，应给各面独立权重与法向；不能仅缩小 shell 因子凑到 PASS。
- 首轮厚长方体保留足够内部粒子，避免一开始全域都是边界；随后保持面内解析度，减小厚度测试双面重叠。

量化目标：以 Jchar=σω|A0| 为参考，J的归一化 RMS≤1%，φ梯度误差≤3%；几何面通量与各面法向电流分别检查。由于真 J=0，不使用 |Jn|/|Jt| 作为主门槛。

### 10.3 第二测试：非零电流圆柱

使用 E-B 无自感轴向磁场圆柱的周向电流，理论对侧壁与端面均无穿透。沿用统一实现检查 n·J、切向电流、Q和场误差，防止“把电流全部压成零”也通过第一测试。

法向以几何面定义，角／棱线分面处理；单个体粒子两侧贡献不混成平均零法向。

## 11. E-C1：迁回 TEAM7 的条件与有限任务

仅当简单厚体、薄体双面及非零电流测试均有可信结果后，迁回 TEAM7 50 Hz。

先确认当前 TEAM7 是否调用与 E-A 相同的装配及边界 profile，不能继续测试旧实虚混合路径。

TEAM7 顺序：

1. 固定真实几何、材料、线圈安匝与峰值／RMS、探针位置和单位；从原始基准数据核对，不凭记忆补。
2. 纯线圈磁场与独立 Biot–Savart／细分收敛源场比较。实验带导体的总场不能直接当作纯线圈准确性参考。
3. 边界修复前后在同配置板上比较，输出各物理表面电流和薄板双面处理；不要只输出 mask数量。
4. 一向感应测试若感应场过大，先检查单位、ω、σ、体积、源电流和算子守恒。Bind/Bcoil大本身不是所有物理情形都非法，但当前实验严重不符必须解释。
5. 一向近似不是总场全耦合解；不能要求其直接等同完整实验，也不能因其失败永久禁止进入合理的自洽求解。只有激励、离散与边界基本正确后才按实验目标运行受控 Picard。

保留原 TEAM7 容差，给逐点曲线与固定尺度 L2；不得做 j_post_scale、最佳幅值拟合或移动探针位置来授予物理通过。

若简单边界仍未过，本轮 TEAM7 状态写“前置边界未过”；停止重复复杂算例，但继续解决已定位的简单边界装配问题。

## 12. 集成：将验证通过的路径接入法国入口

本轮先做映射表，不默认开跑：

| 模块 | 验证入口 | 法国入口需要证明 |
|---|---|---|
| 声学 | 方腔 B＋TVF | 同类 Inner／Wall 修正真实执行，而非只存在类定义 |
| 顶面 | 方腔没有自由滑移顶面 | 法国顶面另行确认，不能继承四壁无滑移结论 |
| 热边界 | 方腔定温 | 法国 Robin／辐射独立路径，不复制 fluid-k CLI |
| 黏性 | 方腔常物性 | 法国变μ及隐式路径另标范围，不自动继承 |
| 电磁 | E-A／E-B共享profile | 实虚部、J/Q和约束参数与验证一致 |
| 热源传递 | EM粒子→流体 | 守恒、空间映射和 freshness真实测量 |

满足 H 的明确可信范围及 E-A/E-B 后，可做法国短回归；TEAM7未过是否阻断，按共用缺陷还是专属薄板问题说明。若 H n160仍不过，不以“法国看着有流动”覆盖，只能做标明限制的探索运行。

短回归先固定热源／物性／初温，用同状态β=0对照排除数值漂移；不通过周期EM自动打开导热来绕过接口缺陷。冻结Q、热扩散和流体开关应能独立控制并输出。

正式运行前证明声学修正真实启用，且传热与浮力共用预期温度字段。后续若启动周期EM，单独核查更新次数与热源传递；`transfer_unmeasured` 不得写作验证通过。

## 13. 排程、停止与继续修复的规则

| 优先级 | H线 | E线 |
|---|---|---|
| P0 | 补实际监测源码／二进制身份，布置诊断 | 冻结统一装配与历史对照；离线功率分解 |
| P1 | n80/n160同窗口定位异常 | 相位／离散恒等式／守恒修复 |
| P2 | 一个有依据的修复及回归 | 全域圆柱复数场与功率拆分；简单边界 |
| P3 | 稳态和分辨率分项验收 | 边界迁回TEAM7；生产profile集成 |

停止的是无信息增益的路线，不是停止解决核心失败：

- 禁止再次只延长方腔时间、只收紧 Picard、只缩小 TEAM7 shell 因子。
- 若一个假说失败，保存证据并转到它排除后剩余的具体问题；不要重复同一开关，也不要把失败简单归档后宣布收口。
- 若基础离散恒等式不成立，优先从函数级装配重构；无需再跑法国长时或复杂TEAM7。
- 在现有范围内允许完成有依据的代码修复与必要回归，不因“只做诊断”的旧措辞停止实现。
- 无法继续的真实阻断须具体写明：缺哪份源码／哪项数据／哪条数学关系未闭合。不能笼统写“需进一步研究”。

## 14. 必须提交的结果

1. `ROUND212_REVIEW.md`：实际修复了什么、证据、未过项、与Doc211相比的数字变化。
2. `ROOT_CAUSE_AND_FIX.md`：每项假说→观测→定位→代码改动→回归；区分已证明与推测。
3. `RUN_INDEX.csv`：真实命令、参数、版本、设备、实际终时、终止原因、分项结果；旧结果标明复用来源。
4. H近壁分区统计、少量完整快照、B／TVF／热流关联图与同配置n80/n160曲线。
5. 电磁两组件装配映射、离散恒等式检查、相位与守恒回归表。
6. 圆柱 `P_cont/P_sample/P_num` 全域拆分、复数场CSV和曲线；提供可复算脚本。
7. 简单厚／薄体绝缘边界结果与非零电流对照；TEAM7是否达到迁移条件。
8. 实际相关源码、diff、二进制哈希及监测实现；不打包build和海量中间输出。

所有工程门槛均按原定义记录，不更改参考值／挑时刻／选择分辨率来变绿。烟测、解析验证、稳态和耦合验证分开。

### 最终需要明确回答

- n160方腔的异常首先发生在哪里，修复后是否消失？
- 电磁是否已经用同一套实虚部离散同时通过相位与守恒？
- 圆柱6%～9%功率偏差中，几何积分、场误差和模型差异各占多少？
- 绝缘边界在简单几何上是否可信，薄板双面贡献是否处理正确？
- 法国驱动是否真正接入验证通过的算子，而不是仍沿用旧默认？

本轮进展应体现为上述问题的答案变得明确、代码缺陷得到修复；不能只增加若干失败日志和新的“可收口”标签。


<!-- ===== END 212_OPHELIE_TARGETED_HYDRO_EM_REPAIR_CURSOR_PLAN_2026-11-06.md ===== -->

---

<!-- ===== BEGIN 213_OPHELIE_LEGACYFLUX_EMF_SKIN_CLOSURE_CURSOR_PLAN_2026-11-07.md ===== -->

<a id="213-ophelie-legacyflux-emf-skin-closure-cursor-plan-2026-11-07-md"></a>

# 213_OPHELIE_LEGACYFLUX_EMF_SKIN_CLOSURE_CURSOR_PLAN_2026-11-07.md

# Doc213 — LegacyFlux／edge EMF 趋肤场闭合（Cursor 计划）

日期：2026-11-07  
前置：Doc212（相位／边门 PASS；场 L2 ~5–6% PARTIAL；F5／F6 证明换 RHS 救场砸边）

## 1. 目标

在**生产路径不变**的前提下（shared complex `phi_rhs=legacy-flux` EMF；LHS=`legacy-pairwise`）：

1. 把圆柱弱趋肤中部 `L2_rel_J_mid` 从 ~5–6% **明显压下**（目标 ≤2%，理想 ≤1%），且  
2. **边门** `edge_rel_l2≤1%` 与 **相位 A1** 不回退。

禁止：把 `--doc212-complex-rhs=div-sigma-a|legacy-pairwise` 升为默认；禁止回退实虚分裂。

## 2. 已钉死事实（继承 Doc212）

| 事实 | 证据 |
|---|---|
| 场债在 LegacyFlux／EMF 离散，非 Picard／ghost／phi_rhs_ghost | `doc212_runs/E_power/EB_FIELD_ROOT_CAUSE.md` |
| 统一 DivΣA：场好、边 ~100% FAIL | F5 |
| 统一 pairwise：同交易 | F6 |
| 生产 shared EMF：边 ~0.23% PASS、场 ~5.5% | EA1 |

## 3. 工作顺序

| 步 | 内容 | 门 |
|---|---|---|
| **G0** | 基线三元组：no-self 生产／DivΣA；+ 生产 Picard+e1；粒子 CSV 含 `rhs_leg−rhs_div` | 复现 Doc212 数字 |
| **G1** | 径向／区域分解场债（rim vs interior；`rhs_leg_minus_div` vs \|J\| 误差相关） | 定位 EMF 哪一段 |
| **G2** | **只改**共享 EMF 装配（权重／A 插值／近壁通量），每刀同时报 edge+L2+A1 | 候选修复 |
| **G3** | 若 EMF 内无法兼得：设计**边一致**新 RHS（须过 e1），再考虑替换生产 | 仅当 G2 穷尽 |
| **P1 并行可选** | E-C 薄体绝缘短测（不阻塞 G0–G2） | EC 薄板 |
| **暂缓** | 法国长回归宣称、TEAM7、SPH-ε 再压 1% | — |

## 4. 运行约定

- GPU 命令由用户执行；**短期跑不进 tmux**。  
- 结果落 `docs/ophelie/reference/doc213_runs/`（gitignored `doc*_runs/`）。  
- 每批跑完后评审再改代码；不自动开法国／TEAM7。

## 5. 验收

- 生产 profile 不变旗下：`L2_rel_J_mid` 改善且 `edge_rel≤1%`、A1 PASS。  
- `ROUND213_REVIEW.md` + `ROOT_CAUSE_AND_FIX.md` 更新链。  
- 未达标不得改法国默认。


<!-- ===== END 213_OPHELIE_LEGACYFLUX_EMF_SKIN_CLOSURE_CURSOR_PLAN_2026-11-07.md ===== -->

---

<!-- ===== BEGIN 214_OPHELIE_CONSISTENT_EM_REPAIR_AND_CONVECTION_VALIDATION_CURSOR_PLAN_2026-11-08.md ===== -->

<a id="214-ophelie-consistent-em-repair-and-convection-validation-cursor-plan-2026-11-08-md"></a>

# 214_OPHELIE_CONSISTENT_EM_REPAIR_AND_CONVECTION_VALIDATION_CURSOR_PLAN_2026-11-08.md

# Doc214：电磁一致性定向修复与自然对流验证收口开发方案

日期：2026-11-08  
面向：Cursor／SPHinXsys OPHELIE 开发分支  
最新电磁输入：`DOC213_ROUND_PACK_2026-11-08.tar.gz`  
最新已审查流体结果：Doc211；Doc213 包未提供新的方腔或法国流体结果。  
目的：用有限、具有明确判别力的试验定位缺陷，完成针对性修复；不再扫描互相补偿的开关。

## 0. 必须遵守的任务定义

本轮包含两条主线，可在独立工作树／设备并行：

- **E：电磁电势、边界、电流和热源相容性修复。** 先做矢势常量平移不变性测试，再依据离散关系修代码。替代离散仅在现有实现明确无法相容时进入。
- **H：恢复方腔 n160 定位。** 不重做已完成的 B／TVF 全矩阵，不再换浮力公式，定位热流偏高与速度波动的共同或独立来源。

功率后处理为离线修正任务，可与两线并行。TEAM7 和法国长时计算暂缓；不把大量复杂工况运行作为本轮成果。

不要把“停止盲扫参数”解释成“只登记失败、停止修复”。在用户授权范围内，完成必要诊断、针对性代码修改和回归，不必为每个局部修复另求确认。

保留所有历史结果与可工作的基线。不得改动独立 A-phi／Laplace 论文分支；本方案中的局部投影候选属于 OPHELIE 诊断路径，不自动合并到独立研究算法。

## 1. 当前状态：哪些成立，哪些尚未成立

### 1.1 电磁

1. 共享实虚部 DivΣA 已通过整体相位旋转，圆柱中部周向电流幅值误差约1%，这是真实进展。
2. 该配置尚未同时通过常量A绝缘测试；开启另一种边界修正后EC0／EC1可过，但圆柱及TEAM7变差。
3. 不同配置分别通过，不能合并成统一方案验证通过。
4. 旧EMF边残差不一定适用于新生产离散；但撤销旧门槛后必须建立与实际生产电流相容的守恒检查，不能只换指标名。
5. 圆柱约1%是幅值指标，完整复数矢量场精度仍未确认。
6. 功率几何误差已得到部分解释；强趋肤正式分解与中段积分口径仍需修正。

### 1.2 浮力与自然对流

| 已审查项目 | 结果 | 能支持什么 |
|---|---|---|
| 移动MMS，B ON，TVF OFF | 0.5 s内部速度L2≈0.893% | 该移动制造解的主要错误得到修复 |
| 方腔B ON、TVF OFF | 发散 | 当前方腔仍需原生TVF配套，不取消 |
| n80方腔，B＋TVF，R210热边界 | t*=0.70定量与稳态窗口通过 | 已有有限范围的工程验证 |
| n160同配置方腔 | t*=1.0，Nu约+6.61%，速度窗口变化约3.3% | 加密与稳定性尚未通过 |
| 一维导热n80／n160 | 接触Nu约−5.5%／−4.6% | 温度误差小，但壁面热流离散仍有偏差 |
| 法国自然对流 | 尚无经审查的完整定量复现 | 不可声明法国环流已验证 |

“浮力公式能工作”“一个分辨率方腔通过”“自然对流完整收敛验证通过”是三个不同结论，不混用。

## 2. 代码身份与最小组织

建立 `doc214_runs/`，子目录：`00_entry`、`E_invariance`、`E_operator`、`E_candidate`、`E_power`、`H_cavity`、`integration`、`scripts`、`source_snapshot`。

开始前：

- 锁定实际Doc213生产profile和回滚版本，保存未提交diff；不要仅用HEAD标识本地修改。
- 核对方腔实际稳态采样代码、后处理脚本和二进制哈希。Doc211包内曾缺少对应监测实现，先补齐真实来源。
- 各测试导出解析后参数、实际模板／算子标识、程序路径和哈希；未知参数报错，不静默忽略。
- 旧数据可复用，但必须与此次代码、配置、误差定义可对应。不要重建已有case或重复写监测脚本。
- 不清理用户未跟踪工作、不强制reset；共享文件的合并和构建串行完成。

以下测试ID是任务编号，不代表仓库已有CLI。新增局部参数时记录真实名称，不编造可直接运行的旧命令。

## 3. E0：先做同配置矢势常量平移测试

### 3.1 固定条件

复用现有圆柱无自感验证，使用同一粒子布局、σ、频率、LHS、RHS、边界、规范、容差和J/Q路径。关闭自感、功率标定、幅值恢复类后处理和J投影。

首先使用未做moment体积修正的基线，避免同时改变几何；后续验证moment模式时单独记录。

三组输入：

| ID | 输入 | 目的 |
|---|---|---|
| E0-R | 原旋转矢势Arot | 原非零周向电流基线 |
| E0-C | 常量A0 | 同配置零电流绝缘检查 |
| E0-RC | Arot＋A0 | 判断物理解是否依赖矢势选择 |

首组A0沿横向x，幅度预先取B0R/2，与旋转矢势边缘量级相当。不是优化参数，不扫描幅度。需要分离实虚问题时，再整体乘i做一次对照。

对于 exp(+iωt)：

```text
E = -grad(phi) - i*omega*A
A' = A + A0
phi' = phi - i*omega*A0 dot (x-xc) + C
E' = E; J' = J; Q' = Q
```

常量A0不改变B。测试属于同一OPHELIE实现的内部不变性检查，不调用独立A-phi分支。若求解器规范固定，应在比较phi时只移除允许的常数C，不拟合线性项来隐藏误差。

### 3.2 禁止的“通过方式”

- 不直接给每个粒子注入解析phi作为求解结果。
- 不在检测到EC0／EC1或常量A时自动换RHS或边界模式。
- 不对E0-RC自动先减去A0后再把结果称为一般离散已验证。解析分离可作诊断，但不能代替修复。
- 不用J后处理投影强制绝缘，同时保留另一套未经修正的Q。
- 同一profile必须完成三组；另一候选可作显式独立对照。

### 3.3 输出和工程判据

保存位置、体积、Ar/Ai、phi_r/phi_i、梯度、Er/Ei、Jr/Ji、Q、边界标签和求解残差。

- E0-R与E0-RC完整复数J相对L2差≤1%，Q场相对L2差≤2%。这些是本项目目标，不是通用文献标准。
- E0-C使用Jchar=σω|A0|归一化，RMS电流≤1%；分母不取应为零的参考电流。
- 分内部、侧壁、端部报告；不能只报全域平均。
- 报告phi变换误差与梯度变换误差；同时保留不做任何对齐的物理场比较。
- 所有残差都输出原始值和尺度，确认失败不是求解未收敛或读错字段。

### 3.4 跑完必须执行的分支

| 观察 | 下一动作 |
|---|---|
| 三组都过，但旧EC0失败 | 查入口／参数／规范／后处理差异，先统一代码路径 |
| 内部好、近壁差 | 进入E1，重点检查边界贡献和支撑截断 |
| 全域差，线性phi也不能再现 | 进入E1，先修梯度／装配，不碰TEAM7 |
| 不变性过、绝对解析场仍不准 | 检查几何、参考模型和空间离散，不能继续归因相位 |

## 4. E1：把当前实现写成可检查的离散关系

### 4.1 小规模算子审计

在几百粒子、固定几何问题上，用当前代码实际算子导出或通过基向量施加得到：

- G：phi到实际计算E所用梯度；
- K：实际电势LHS，不含规范项与含规范项分别保存；
- b(A)：实际激励RHS，体积项与边界项分开；
- J/Q所用重构或边关系；
- 规范约束、源项投影、数值缩放。

矩阵审计可在CPU双精度进行，先与同布局GPU算子作用对照。大矩阵只作稀疏或matrix-free处理；禁止对实际大算例逐基向量生成密集矩阵。

### 4.2 必查关系

1. 常数phi的梯度应接近零；线性phi在内部及近壁能否再现常量梯度，单独报告。
2. E0-C解析phi施加到K后，Kphi与b的差在哪里；先剔除规范常数差，再分析边界贡献。
3. 解析注入仅用于定位：随后必须从正常初值求解，再验证真实输出。
4. 若用边通量，逐边反对称性、节点汇总与方程残差关系写清Vi和符号。
5. 若用重构电流，明确对应弱散度或离散散度。不能把不相关的旧边残差当硬门，也不能取消所有物理守恒检查。
6. penalty*phi、零均值源项投影和边界通量全部入账。求解残差小不代表原始无泄漏电流问题被正确求解。
7. 对称性／正定性／零空间按真实矩阵检验；不能未检查就因为旧代码用PCG而继续强用PCG。

### 4.3 先尝试最小修复

优先修复已证明的：边界遗漏或重复、符号、体积权重、错误active A字段、实虚参数漂移、规范泄漏、近壁梯度不一致。

每次只针对一个明确缺陷修改；但一个数学关系需要共同修改LHS与RHS时，应作为一套一致修复，不机械限制为只能改一个文件或一个算子。

禁止再扫描Neumann amp、shell厚度、罚项任意值以凑三个门同时通过。边界几何厚度可做有物理解释的离散研究，不是禁止一切几何调整。

原配置与修复配置各保留一个明确ID。局部修复后立刻回归E0三组、整体相位和原圆柱场，不等TEAM7跑完才发现退化。

## 5. E2：现有算子无法相容时的条件性替代

### 5.1 启动条件

只有E1已定位到“现有K、b和生产电流构造不相容，局部接线修正不足以解决”时，启动OPHELIE局部一致变分投影原型。

不把已测失败的 `DivSigmaGrad` 标志再次当新方案。必须明确新原型的梯度、伴随、权重、约束和边界，不能只是换名字。

### 5.2 原型数学定义

在固定导电域、给定A、σ>0下，令G为同一个离散梯度，M为体积×电导率权重，在复数空间求：

```text
min_phi (G phi + i omega A)^* M (G phi + i omega A)

K = G^* M G
b = -G^* M (i omega A)
E = -(G phi + i omega A)
J = sigma E
Q = 0.5*Re(J dot conj(E))
```

星号为共轭转置；G为实矩阵时对应普通转置。当前M中Vi与σ必须为正，记录其范围。

先处理常数零空间，采用明确参考点或零均值约束；检测额外零模态／近零模态，不盲加大罚项压住。

这一定义保证的是相应离散弱残差关系，不自动保证逐粒子点值divJ、逐边局部守恒或所有表面点Jn精确为零。自然绝缘边界的离散表达仍需要依据自由度和积分方案检查，不能由“变分”两个字直接宣布完成。

若G不能再现线性场，E0常量平移仍可能失败；先解决梯度支撑和秩，不用最小二乘形式掩盖坏梯度。

### 5.3 实施顺序与退出条件

1. CPU双精度小问题建立明确参考，检查秩、对称性、残差、能量和E0三组。
2. 两档粒子布局检查近壁Jn及圆柱场，确认不是零电流过而非零电流被抹掉。
3. 参考原型成功后，再实现SYCL matrix-free `G` 与严格配对的 `G^*`；检查伴随内积恒等式，GPU与CPU一致。
4. 优先复用现有线性求解基础设施；方法选择按真实矩阵性质，不要求无依据地全部GPU隐式重写。
5. 如果原型失败，输出具体失败关系或零模态，回到G／边界构造；不并行再开第三、第四种候选。

这是有明确前提的备选路线，不承诺天然比现有方法更准，也不等价于用户独立A-phi方法的迁移或修改。

## 6. E3：统一profile验收，而非各选各的通过配置

通过的同一套profile必须完成：

| 检查 | 目的 |
|---|---|
| 常量A封闭导体 | 零电流与绝缘抵消 |
| 无自感圆柱旋转A | 非零周向电流保真 |
| 圆柱旋转A＋常量A | 物理场对矢势选择不敏感 |
| 整体相位旋转 | 实虚部一致，Q不变 |
| 已有纯梯度／变σ制造解 | 复用旧测试，公共代码改后回归 |
| 同离散关系守恒 | 方程、实际J与边界相容，而非只看新残差小 |

模式可因真实边界类型变化而改变，但相同物理绝缘边界不能为了不同测试A而切算法。若某模式只适用受限轴对称物理，须明确限定，不称通用电磁方案。

工程容差沿用Doc212／214定义；先前已达1e-6的测试不可明显退化后仅凭仍在1%内称无回退。必要时以相同求解容差比较。

只有同一profile通过上述检查，才继续完整圆柱自感及TEAM7。TEAM7先固定一个频率，不再大矩阵扫边界模式。

## 7. P：功率测量修正与完整复数场验证

### 7.1 先修复中段积分

当前源码按粒子中心选中段并累计整Vi，却以名义带厚计算连续参考；H=0.6、dp=0.02时，中段名义厚度0.09 m，选中层的完整厚度约0.08 m，不能把约9%偏低直接解释成端部物理。

正式口径二选一并固定：

- 同批粒子上积分解析Q，数值／参考均使用同一Vi和同一掩膜；
- 对切片边缘控制体按实际交叠体积积分，数值／参考一致。

同时保留连续名义区间参考作为几何偏差分析，不混为同采样误差。不要只乘一个“有效高度系数”就称严格修复。

### 7.2 正式功率三分解

```text
P_cont   = 对应解析模型的连续积分
P_sample = sum Vi * Q_exact(xi)，所有实际粒子
P_num    = sum Vi * Q_num(xi)
P_num-P_cont = (P_sample-P_cont)+(P_num-P_sample)
```

无自感使用规定A参考；自感强趋肤使用full_skin参考。不能用无自感参考分解强趋肤，再把场差称为数值误差。

有限高度数值模型对无限长参考，最后一项仍含模型差异。高度扫描可提供证据，但不是所有余项都已被定量解释。

moment体积修正保留opt-in；无自感Q∝r²时匹配ΣVr²后功率接近零误差是预期积分性质，不能单独作为场算法精度证明。不直接推广到法国流体质量或任意复杂几何。

### 7.3 完整场

保存Jr/Ji、Er/Ei各分量以及同点解析值，计算体积加权完整复数矢量L2。幅值、相位、径向／轴向杂散分量分别报。

不拟合整体相位后隐藏主误差；低场相位阈值预先固定。先离线复用已有全域字段，没有字段再最小重跑。

目标：可比区域完整复数J误差≤5%；同采样功率误差≤3%；连续总功率≤3%独立记账，不能前者通过就宣称后者通过。

## 8. H：方腔尚欠哪些工作，怎样有控制地修复

### 8.1 本轮范围

保留n80已通过的R210：B＋原生TVF、ghost OFF、fluid-k。物理严格Ra=1e5、Pr=0.71；不改变β、μ、ΔT或边界来凑参考。

已有n80和n160速度峰值都进入10%误差范围，但n160 Nu与速度稳定性未过。因此本轮主问不是“有没有浮力”，而是“加密后为何还有近壁／运动／热流不一致”。

### 8.2 一次定位运行

优先兼容restart，在实际波动窗口采样；无restart则从相同初态重跑并说明。每Δt*=0.02至少10点，连续3窗；保留3～5个完整快照，不输出海量VTP。

按物理壁距分0～dp、dp～2dp、2dp～4dp、内部与角点：

- 速度、密度、EOS压力及各力加速度；
- 最近粒子距、邻居数、核矩／B条件数；
- TVF纯位移/dp及物理位移，修正前后几何；
- 温度、各层接触功率及异常粒子对总热流贡献；
- 固定空间网格速度／温度用于稳态，粒子唯一ID用于追踪。

必须核查热算子在更新Position后使用的邻域和核几何时刻是否相容，先取少量粒子暴力核支撑检查；不机械要求每声学子步重建邻居。

### 8.3 根据最先异常的量修复

| 首个异常 | 动作 |
|---|---|
| 几何／邻居错误先出现 | 修时刻、更新或字段同步；同快照算子复核 |
| B极端／近壁支撑退化 | 查构造与墙贡献，不先钳制B |
| TVF堆积／壁向偏移 | 对照原生闭腔流程，修实际差异，不扫TVF系数 |
| rho／p先出现高频异常 | 核查EOS、Vol、壁声学和历史外力；必要时一次dt/2短窗 |
| 热流由过近少数粒子主导 | 检查壁距、跳跃、接触权重与积分面积，导热／方腔共同回归 |

固定快照CPU双精度复算仅在证据指向消减或条件数时使用，不自动把全程序改双精度。

### 8.4 H验收

- 常物性一维导热先过，声学改变时复测成功MMS。
- 方腔n80／n160均报完整中心线与温度剖面，不只峰值。
- 工程目标：Nu≤5%误差、中心线U/V各≤10%、峰位置差≤0.05L。
- 连续三个窗口均值变化：空间速度L2≤1%、温度/ΔT≤0.5%、Nu≤0.5%；同时报告波动幅度，不用均值抵消噪声。
- 必须说明两分辨率差异，不授予严格收敛阶；需要论文收敛研究时，在缺陷解决后再加有依据的第三档，不现在盲跑n320。
- 不再单靠延时追稳态。修复必须有“异常→代码原因→前后对照”的证据链。

## 9. 论文验证范围：当前可以写到哪里

### 9.1 当前可支持的表述

在给定n80配置上，标准差分加热方腔的中心线速度、换热量及稳态监测达到了本项目工程门槛；这证明该配置能产生合理的浮力驱动自然对流，并有定量基准支持。

### 9.2 当前不支持的表述

- 自然对流求解器已完成网格／粒子收敛验证；
- 所有分辨率或所有热边界均准确；
- 法国高Pr玻璃、变μ隐式路径、自由滑移顶面已经由空气方腔一并验证；
- 法国双环流或完整电磁—热—流耦合已经定量复现；
- 方腔通过说明任何误差都只来自法国物性或文献边界不确定性。

### 9.3 写文章前最小补齐项

1. 锁定最终算法和实际源码，解决或清楚界定n160热流与稳定性问题。
2. 同配置分辨率研究，完整速度与温度曲线、Nu、稳态和热能量平衡，不能仅挑n80。
3. 热边界解析验证及热流定义解释；一维n80接触Nu−5.5%缺口真实报告。
4. 作为法国应用文章，增加与其高Pr／高黏性路径相近的一项针对性验证，例如同Pr的可信参考解或同物性制造解；不能套用空气方腔参考值。
5. 若保留变μ／隐式黏性，补相应算子与时间步准确性证据；若顶面自由滑移对结果关键，补独立边界检验。
6. 法国流线拓扑、温度分布和适用参考数据独立对照；方腔验证与工程验证在文章中分开。

不为了发文章立刻新增十个case。优先完成上述已知缺口，再按论文实际使用的物理路径决定附加验证。

可选英文写法（目前阶段，用于内部草稿）：

> At the tested n80 resolution, the cavity calculation agrees with the benchmark within the adopted engineering tolerances. Resolution sensitivity and the remaining wall heat-transfer discrepancy require further assessment.

正式投稿时按最新结果重写，不把这段限制性描述当作永久免验证条款。

## 10. 集成与停止条件

### 10.1 不再无休止测试的执行规则

每次运行前填写：唯一问题、预期能区分的两种结果、结果各自触发的代码动作。无法填写则暂不运行。

一轮只维护原基线和一个修复候选；旧实验参数可留作明确诊断，但不形成并行生产配方。

若E1证明结构不相容，进入E2，不重复壳幅扫描；若E2发现坏梯度或额外零模态，修对应构造，不开第四套RHS。

有限预算优先用于修复后回归，不用于重复失败的法国长算或TEAM7参数矩阵。

### 10.2 法国短回归准入

先完成E0/E3统一profile与H可信范围。法国入口必须输出实际修正声学、TVF、热边界、黏性和电磁profile；`kernel_correction=0`的旧短跑不继承R210验证。

Robin／辐射与方腔Dirichlet分开。冻结Q、热扩散与流体开关独立控制，不能靠开周期EM绕过冻结热源无导热的接口问题。

短跑分别报告Tmean变化与Tmax−Tmin；初态非平衡冷却不称整体升温。缺失热源传递／freshness测量不得记为耦合验证通过。

若n160仍未解决，只能说明受限探索的范围，不以法国出现环流覆盖基础问题。TEAM7特有薄板难点不自动阻断轴对称探索，但共用守恒缺陷必须阻断可信电磁耦合结论。

## 11. 必交付文件

- `ROUND214_REVIEW.md`：实际修复与数字变化，主线分开结论。
- `E0_INVARIANCE.csv`：三输入同profile及误差分区，保留失败基线。
- `DISCRETE_OPERATOR_RELATIONS.md`：K／G／b／边界／规范／J的实际关系，定位证据和修复。
- `CANDIDATE_DECISION.md`：现有实现局部修复或条件性E2原型的选择依据；未触发E2则明确NOT_RUN。
- `POWER_SPLIT.csv`：正确参考下连续／同采样／数值功率及切片体积；完整复数场数据与脚本。
- `H_ROOT_CAUSE_AND_FIX.md`、近壁分区与稳态CSV、n80/n160基准曲线。
- `VALIDATION_SCOPE.md`：能支持的论文表述、尚未验证的法国路径。
- `RUN_INDEX.csv`、真实命令、源码diff和二进制身份；图与结论能由打包数据复算。

验收最后回答：同一配置是否同时解决零电流、非零电流、常量平移不变性和守恒？n160方腔异常是否找到并修复？如果尚未解决，下一步是否已经缩小到具体矩阵关系或局部代码，而不是又增加一批待优化标签？


<!-- ===== END 214_OPHELIE_CONSISTENT_EM_REPAIR_AND_CONVECTION_VALIDATION_CURSOR_PLAN_2026-11-08.md ===== -->

---

<!-- ===== BEGIN 215_OPHELIE_CAVITY_NEIGHBOR_AND_EM_ADJOINT_REPAIR_CURSOR_PLAN_2026-11-08.md ===== -->

<a id="215-ophelie-cavity-neighbor-and-em-adjoint-repair-cursor-plan-2026-11-08-md"></a>

# 215_OPHELIE_CAVITY_NEIGHBOR_AND_EM_ADJOINT_REPAIR_CURSOR_PLAN_2026-11-08.md

# Doc215：方腔邻域／传热修复与电磁真实伴随验证开发任务书

日期：2026-11-08  
执行对象：Cursor；项目：冷坩埚招投标／OPHELIE。  
电磁基线：`DOC214_ROUND_PACK_2026-11-08`。  
流体基线：实际包含 R210／Doc211 改动和稳态监测的工作树；Doc213／214 没有新的方腔结果。  
本文是下一轮实施要求，不是已完成结果。

## 0. 本轮的两项核心交付

**H线必须实际启动：**查明位置更新后热计算是否使用不相容邻域，完成有依据的修复；在同一算法下比较n80／n120／n160，解释热流偏高与速度波动。

**E线先纠正实现：**修复校正梯度伴随中的邻居矩阵与体积权重问题，先完成CPU显式矩阵基准，再做GPU算子对照和物理回归。不要以错误伴随的失败否定整个变分路线。

两线可在独立设备／工作树并行，不相互等待。共享文件集成、构建和生产profile变更串行进行。P功率后处理只作离线附属任务，不抢占H/E主线。

禁止默认开展：新的Neumann幅值／壳厚扫描、第三套散度族、TEAM7长跑、法国长时耦合、为吻合基准修改β／物性／边界温度。

不修改独立A-phi／Laplace论文分支。本轮局部电磁原型不等同于该研究分支的修改或迁移。

## 1. 执行边界与历史结论修正

### 1.1 需要保留的成果

- n80方腔在R210配置下主要特征量和窗口稳态达标；保留结果，不改写。
- 共享DivΣA实虚路径已解决整体相位不一致，仍作为失败定位对照。
- 生产E0常量A测试失败、解析线性phi注入后抵消成功，是确定的装配缺口证据。
- LegacyFlux／全矩壳能过常量平移但损害无自感圆柱场，保留作对照，不再升默认。

### 1.2 不能沿用的过强结论

- G2“场差”不能用于否定正确的校正变分离散：当前伴随公式存在实现缺陷。
- G2b用校正散度替代伴随后爆炸，不能自动证明真正的转置组合不稳定。
- 线性phi注入成功只证明该测试的重构性质，不证明任意粒子梯度、曲面或零空间全部正确。
- Doc214表中的P0近机器零误差来自`no-self-induction / prescribed_A_no_self`，不是趋肤求解达到机器精度。
- H线、P功率口径在Doc214未执行，不能记成已处理或已排除。

## 2. P0：版本、程序与目录

建立 `doc215_runs/`：`00_entry/`、`H_snapshot/`、`H_cavity/`、`E_cpu/`、`E_gpu/`、`E_physical/`、`P_postprocess/`、`scripts/`、`source_snapshot/`。

1. 记录真实HEAD、分支、工作树、未提交diff、二进制路径与哈希；不能只记录同一HEAD。
2. 找到实际方腔稳态CSV生成代码和脚本，并随包提供；不要把旧快照当成最终运行源码。
3. 每次运行保存完整参数、解析后实际配置、设备、精度、请求／实际终时和终止原因。
4. 旧结果能复用则复用；改造公共执行逻辑后，成功基线至少做一次必要回归。
5. 不清除用户未提交工作，不使用强制reset；不因本任务自动push。
6. 延用已有入口和监测，新增参数只控制本轮必要诊断；未知参数必须报错。

## 3. H0：方腔固定条件与问题定义

严格二维差温方腔：L=0.0386 m，Ra=1e5，Pr=0.71；左壁303 K、右壁283 K；上下绝热；四壁无滑移；初始293 K。

数值基线：声学B修正＋原生TVF，关闭自定义density normalization；`thermal_dirichlet_ghost=false`、`thermal_dirichlet_fluid_k=true`。不增加particle shifting、不调整TVF系数。

**不引入额外ghost／镜像粒子。**继续使用现有壁面粒子。旧温度ghost开关只是虚拟温度跳跃，不生成新粒子；本轮保持关闭，不把它作为调参工具。

当前需要解释的数据：

| 量 | n80，t*=0.70 | n160，t*=1.0 |
|---|---:|---:|
| Nu_contact | 4.4506 | 4.8178 |
| Nu_fit，近壁4dp线性拟合 | 4.1016 | 4.1042 |
| 中心线V峰值误差 | −3.69% | −0.92% |
| 中心线U峰值误差 | −8.84% | −5.82% |
| 窗口速度变化 | 约0.50% | 约3.3% |

接触热流变化、拟合梯度近乎不变是诊断线索，不表示拟合Nu必然正确。温度剖面曲率和拟合区间变化必须考虑。

## 4. H1：移动后邻域的同快照对照——先做这一项

### 4.1 代码疑点

已审查快照顺序：

```text
声学推进 → UpdateParticlePosition → 热扩散子步
→ cell更新 → relation更新 → 下一轮力学
```

热计算可能使用新Position与旧邻居索引。是否产生误差要由实际核距离、缓存和邻居覆盖检查决定，不预先认定为根因。

### 4.2 选取快照

优先从可复现的n160波动窗口取一个完成位置更新、尚未热推进的状态；n80取同类状态作对照。无兼容restart时重跑并标明，不能伪造续算。

保存Position、粒子ID、Velocity、T／温度增量、rho、Mass、Vol、墙字段、时间步、B、原邻域信息以及本步TVF位移。确保CPU/GPU同步，不按排序后数组下标混配。

### 4.3 A/B计算，状态必须相同

- A：保持当前旧邻域，计算Inner／Wall热算子，不推进温度。
- B：只重建cell／relation及必要几何缓存，计算同一热算子，不改变T／rho／Vol。
- 保存两个分支的原始算子输出；避免B覆盖A或被下一次计算清零。
- 对64个分层选择的粒子，CPU暴力核支撑检查应有邻居、旧列表和新列表；包含侧壁、角点、内部与异常速度附近。
- 检查关系更新是否重排粒子、重分配缓冲区、使旧InteractKernel指针失效；按框架规范更新kernel，不只调用表面API。

### 4.4 必交指标

| 比较 | 定义／要求 |
|---|---|
| 邻居覆盖 | 漏邻居／多余邻居、是否仍在实际核支撑内 |
| 热算子差 | 全域／0～2dp近壁的加权L2及最大差，明确归一化尺度 |
| 壁面功率差 | 热壁、冷壁独立积分；不能只看冷热差接近零 |
| 一步温度差 | 用相同热时间步预测，不启动整段正式运行 |
| 位置变化 | 对流步物理位移/dp、TVF位移/dp，判断旧支撑适用范围 |

对近零算子用固定kΔT/L²或其他事先声明的尺度，不用近零值作分母。相对差约1%可作为需要关注的工程信号，但不是唯一硬门；同时看绝对量和误差空间分布。

### 4.5 跑完后的确定动作

| 结果 | 下一步 |
|---|---|
| 旧列表实际漏邻居，热流／局部温度差显著 | 修正位置后热前更新顺序，进入H2 |
| 列表相同但缓存几何不同 | 修缓存刷新／算子生命周期，再做同快照对照 |
| A/B在数值误差范围相同 | 明确排除本快照的邻域解释，保留数据；转H3，不把更改顺序包装成修复 |
| 仅个别状态有差异 | 再取一个正常／峰值快照，定位发生条件；不展开长时参数扫描 |

## 5. H2：有证据后实施局部顺序修复

候选顺序：

```text
声学推进 → 更新Position → 更新cell和relation／几何
→ 热扩散子步 → 下一轮力学
```

1. 保持原双时间步物理和声学子步逻辑，不强制每声学子步重建邻域。
2. 固定壁面未移动时，不机械重复不必要的壁cell构建；以框架实际依赖为准。
3. 热子步内若Position不再改变，可复用一次新邻域，不每个热子步重建。
4. 不同时改变B更新、黏性、浮力、壁温公式、粒子排序或热时间步系数。
5. 确认下一轮力学使用预期的新几何；若新顺序已更新，避免再无意义重复重建。
6. 热子步累计时间必须覆盖实际流体推进时间；触及迭代上限要报错／记录未完成，不能静默跳过。

## 6. H3：近壁热流与速度波动定位

邻域修复后仍异常或H1已排除邻域时，收集一个短窗口，按当前物理壁距分：0～dp、dp～2dp、2dp～4dp、内部、角点。

每区输出速度／压力／密度RMS和P95／最大值，邻居数与最小间距，B最小奇异值／条件数，TVF纯位移/dp，近壁温度跳跃和热流贡献。

### 6.1 壁面热流

- 按实际接触壁对象分冷热壁，检查角点归属和积分面积。
- 输出前1%、5%高贡献粒子占总功率比例，检查是否少数过近粒子主导。
- 核查壁粒子中心、物理壁面和流体粒子位置，不把壁粒子温度赋值自动等同于物理面上的精确定温。
- 独立拟合Nu比较2dp／3dp／4dp，并说明是否约束壁温；不能挑最接近参考的一档。
- 不将harmonic或fluid-k单独称为普遍正确；温差、距离、虚拟温度、体积和权重需整体相容。

### 6.2 TVF与温度

测量真实方腔TVF位移，不能用MMS中极小位移来代替。比较纯TVF整理前后固定空间温度采样，确认是否引入额外温度搬运。

若要加入温度输运补偿，必须先依据原生TVF定义推导；不能直接复制其他shifting方法公式。本轮不因名字近似就视为等价。

### 6.3 动量波动

记录异常前后压力／声学、黏性、浮力加速度，确定最先异常的是粒子几何、B、压力还是速度。只根据证据修对应实现。

必要时只做一次相同状态起步的dt/2短窗口对照；不同时调整声速与核支撑。怀疑精度时先做局部算子CPU双精度复算，不迁移全程序。

## 7. H4：修复后n80／n120／n160统一检查

| 配置 | dp=L/n | 二维流体粒子约数 |
|---|---:|---:|
| n80 | 0.4825 mm | 6400 |
| n120 | 0.32167 mm | 14400 |
| n160 | 0.24125 mm | 25600 |

壁粒子另计。三档物理条件、核函数、h/dp、TVF系数、边界公式一致；时间步按同一稳定性规则生成，不要求数值dt完全相同。

先导热和n80短回归，随后n160检查改动是否针对异常；再补n120形成趋势。保留旧失败结果，不要求重跑所有旧算法长时组合。

初始计划t*=0.70，检查三个Δt*=0.02窗口，每窗≥10样本；必要时到1.0，但不再靠无限延时解决已成平台的误差。已有终时不同的旧结果只作历史参考，新的比较需对齐稳态定义。

验收：

- 参考Nu=4.519；Vmax*=68.59@x/L=0.066；Umax*=34.73@y/L=0.855。
- 工程误差目标：Nu≤5%，中心线U/V各≤10%，位置差≤0.05L；不能修改门槛来变绿。
- 窗间变化：固定空间速度L2≤1%，温度/ΔT≤0.5%，Nu≤0.5%；同时报告窗内波动与覆盖率。
- 导出完整中心线速度／温度曲线和Nu随时间；未取得参考曲线的量标记REFERENCE_MISSING，不编造对齐结论。
- n120用于判断趋势，不用于寻找可替代n160的“通过分辨率”。三档结果不单调时解释原因，不强行拟合收敛阶。
- 分开报告完成运行、特征量误差、稳态、能量平衡和分辨率趋势。任意一项失败不被综合smoke标签覆盖。

能量：按实际热子步累计施加功率，与ΣmcpΔT比较。核对真实符号：现有实现为`dT += -dt*Lap/(rho0*cp)`，审计应使用实际算子定义，不能按旧文档中相反符号推导。

## 8. E0：修复前先保存伴随错误的最小反例

目标文件：`electromagnetic_ophelie_doc214_e2_variational.h`。核对`ComputeOphelieE2ScalarGradientCK`与`ComputeOphelieAdjointGradientDotCK`。

当前问题：G在粒子i使用Ci；“伴随”却把来自vi和vj的贡献都使用Ci。Ci≠Cj时一般不等于真实转置，近壁尤其需要检查。

最小代数反例（等体积、两节点、一维分量，Ci=2、Cj=1）：

```text
G       = [[ 2, -2], [ 1, -1]]
H_old   = [[ 2,  2], [-1, -1]]
G^T     = [[ 2,  1], [-2, -1]]
H_old G = [[ 6, -6], [-3,  3]]
G^T G   = [[ 5, -5], [-5,  5]]
```

取x=(0.7,1)，旧组合xᵀKx=−0.36，正确组合=0.45。这是公式反例，不是圆柱实际仿真，也不能保证纠正转置后所有物理测试一定通过。

本轮应把该反例做成一个极小、独立的代数测试，连同真实几何矩阵测试。不要继续仅依赖圆柱场失败来判断伴随正确与否。

## 9. E1：CPU双精度显式参考，先固定数学定义

### 9.1 实际G的构建

从同一真实小粒子布局导出Position、Vol、σ、邻居、核梯度和Ci。先在几百粒子以内：

```text
(G phi)_i = sum_j a_ij * (phi_i-phi_j)
a_ij = 实际粒子i那一行使用的校正核系数向量
```

矩阵G为3N×N，明确行列、单位、正负号。每行自己的Ci只在构建该行时使用。

优先直接导出系数；若通过逐基向量调用实际算子构建G，须确认系数、邻域、规范和临时字段在每次调用中不变。仅对小问题这样做，不用于大算例。

### 9.2 两种合法内积形式只能选一种作为主实现

定义Mv=diag(Vi)，Mσ=diag(Viσi)，向量场对应的体积权重扩展到三分量。

**首选CPU参考：普通欧氏内积。**

```text
K = G^T Mσ G
b = -G^T Mσ F
F = i omega A              # 复数整体写法
E = -(G phi + F)
J = sigma E
```

实虚分裂沿用exp(+iωt)：Fi=+ωAr，Fr=−ωAi。不能因为代码有a_sign就省略独立符号复核。

**若GPU沿用体积加权PCG：**

```text
Kv = Mv^(-1) G^T Mσ G
bv = -Mv^(-1) G^T Mσ F
```

Kv在Mv内积下自伴，而非一般欧氏对称。需检查`Mv Kv = Kv^T Mv`。或者进行对称变量变换后用普通PCG。

不得把欧氏K直接塞进体积加权PCG，尤其非均匀Vol时。不得重复或遗漏Vi／σ。

### 9.3 转置的安全实现原则

CPU参考直接用矩阵转置，不手写另一套“类似散度”。matrix-free转置按原G行系数散射：

```text
for each directed row i and neighbor j:
    s = dot(a_ij, vector_i)
    out_i += s
    out_j -= s
```

体积和σ权重可先乘入vector，但只乘一次。若改为gather，必须正确使用反向行的a_ji，包含Cj、Vi和实际支撑，不能简单把Ci换成Cj就认为所有权重已正确。

若邻域非互易，gather还需要真实反向邻接。散射原G非零项天然给出转置，是小参考的优先方式。

### 9.4 规范与零空间

先检查G常数零模态及额外零模态。采用明确固定参考点或零均值约束消除常数自由度；固定点时做一致消元，不能只覆盖一行破坏对称性。

不把全域`penalty*phi`称为纯常数规范，它会改变非恒定电势。参考验证优先不靠大罚项；若保留正则化，单独报告其影响，不能用它掩盖秩缺陷。

## 10. E2：矩阵验收通过前，不跑物理候选比较

CPU双精度小问题至少覆盖：等体积、不等正体积、不同Ci、常σ、正的空间变σ。用固定随机种子生成测试向量。

| 测试 | 内容 | 工程目标 |
|---|---|---|
| 伴随恒等式 | dot(Gx,y)与dot(x,Gᵀy)，同权重定义 | 相对归一化差≤1e-10，病态时报告条件数与绝对量 |
| 对称性 | K−Kᵀ；或MvKv−KvᵀMv | 归一化Frobenius≤1e-10 |
| 非负性 | xᵀKx及小矩阵谱 | 负值不得超过矩阵尺度×舍入容差 |
| 常数模态 | G1 | 与浮点误差一致 |
| 额外零空间 | 约束后最小奇异值／特征值 | 无未解释零模态，不能盲加罚项 |
| matvec | 显式Kx与分步G／权重／转置 | 与数值精度相符 |

随机二次型只能辅助，不能代替小矩阵谱检查。CPU代数成立只证明实现内部结构正确，不代表SPH边界／物理精度已自动正确。

如果上述失败，直接修系数、内积、约束或转置；禁止先跑更多圆柱／TEAM7寻找“较好配置”。

## 11. E3：CPU物理参考，再修GPU

### 11.1 同一profile三组输入

复用无自感圆柱：R旋转A、C常量A0、RC=R+C。A0预定为(B0R/2,0,0)，不扫描幅值。

关闭功率归一化、J后处理投影、自动切边界模式、moment体积修正。φ必须真正求解；解析phi注入仅作单独算子检查。

工程目标沿用：C的Jrms/Jchar≤1%；RC与R的完整复数J差≤1%、Q场差≤2%；非零圆柱完整复数J对同位置解析参考≤5%，不是只看相关系数。

必须导出全部Jr/Ji、Er/Ei、phi和Q。无自感解析场幅值准确不等于趋肤／自感通过。

### 11.2 若代数通过而场仍差

将问题分开：

- G是否在全域再现线性场；
- 正确K是否存在额外近零模态／棋盘模态；
- 旋转A的GᵀMσF主要集中在哪里，是否是曲面采样／积分支撑误差；
- 电流是否满足相应离散弱守恒，物理表面Jn是否随加密改善。

不能要求正确的转置恰好把任意离散旋转A映射成零；也不能仅因连续divA=0就删去所有边界贡献。

只做一个由诊断支持的梯度／几何修复，伴随继续从同一G生成。若仍失败，给出秩／再现／积分的具体结果，不再开第三种手写伴随或退回壳幅扫描。

### 11.3 GPU实现

CPU参考的代数和物理都通过后，再修SYCL路径：

1. 按同一行系数实现G与转置；并行散射需原子或安全归约，gather需真实反向系数。
2. 每次调用正确清零输出，核依赖、GPU同步和字段覆盖明确。
3. 先对照同快照Gx、Gᵀy和Kx，再对照完整解；单精度工程目标可先取1e-4相对差，误差尺度与求解容差一起记录。
4. 使用和矩阵性质匹配的求解器／内积；旧LegacyPairwise对角可作经检验的正预条件器，但不能假装是新K的真实对角。
5. 保持CPU参考，GPU失败先定位算子差而非再改数学公式。

不预先承诺GPU速度。必要时CPU参考先作为准确性依据；最终优先GPU，但不以速度牺牲本轮修复正确性。

## 12. E4：共同回归与生产限制

同一修复profile还需整体相位旋转、已有纯梯度／变σ测试、离散守恒。不得按测试名自动选择不同算法。

常量phi注入通过不能注销曲面边界、额外模态或局部电流检查。变分弱残差通过也不自动等于逐边守恒与逐点Jn为零，必须清楚定义实际检验量。

上述通过后才恢复圆柱Picard，保留full_skin参考，做一次迭代收紧确认最终φ、A、J、Q属于一致状态。

TEAM7暂不跑长时或开关矩阵。只有统一profile通过基础非零／零电流检查后，才安排一次固定50Hz的受控回归。法国生产默认不在本轮修复过程中来回切换；新profile未过前，旧DivΣA只能称受限基线，不能称已验证正确。

## 13. P：附属功率口径修复

利用已有字段和脚本离线完成，不另开算法：

```text
P_cont   = 正确解析参考在连续区域的积分
P_sample = 同一实际粒子／体积上解析Q的积分
P_num    = 同一实际粒子／体积上数值Q的积分
P_num-P_cont = (P_sample-P_cont)+(P_num-P_sample)
```

无自感／强趋肤参考分开。中段按实际相同体积比较，不把中心筛选整粒子与名义0.09m带厚的差异称为物理端部效应。有限高对无限长参考的模型差异仍单列。

保留原始功率和所有归一化前数据；moment模式保持诊断opt-in，不用改功率表的方式获得通过。

## 14. 工作顺序与防止再次跑偏

| 批次 | H线 | E线 |
|---|---|---|
| 第一批必须完成 | H1旧／新邻域同快照对照 | 两节点反例＋实际CPU G矩阵与真转置 |
| 第二批 | 证据支持的局部修复，导热／n80短回归 | E2矩阵性质与内积验收 |
| 第三批 | n160异常窗口定位／修复验证 | CPU三输入物理验证 |
| 第四批 | n80／n120／n160同配置比较 | GPU逐算子对照、基础回归，必要时Picard |

如果H线未启动，不能把电磁报告当作本轮全部交付。不能因某一线未过而让另一线无限等待。

每次新运行先写：唯一问题、两种可能结果、各自触发的代码动作。代数未过不跑物理；局部快照能回答的不跑长时。

本轮允许修复及必要回归，不把“停止错误实现”解释成“禁止改正确实现”。如果已确认结构不对，继续按数学定义修，不新增类似开关。

不强行要求全部绿色；但未过必须缩小到具体算子、矩阵或局部几何关系，不能继续泛称“需要优化”。

## 15. 结果包必须提供

### 15.1 入口

- `ROUND215_REVIEW.md`：两线各自完成了什么、实际改动、误差变化和剩余阻断。
- `RUN_INDEX.csv`：包含复用／执行／NOT_RUN、真实参数和分项验收。
- 实际源码diff、二进制哈希、监测代码、真实执行命令。

### 15.2 H线

- `H1_NEIGHBOR_SNAPSHOT_COMPARISON.csv`：旧／新邻域热流、算子和覆盖差。
- `H_CALL_ORDER_BEFORE_AFTER.md`：真实执行顺序及缓存／kernel生命周期。
- 近壁分区统计、异常时序、接触功率贡献与独立梯度对照。
- `CAVITY_N80_N120_N160.csv`、完整中心线与稳态窗口、热能量账本。
- `H_ROOT_CAUSE_AND_FIX.md`：异常→定位→代码修复→对照；未证实内容标假说。

### 15.3 E线

- 两节点反例脚本与真实几何G、权重、约束及CPU参考。
- `ADJOINT_MATRIX_CHECKS.json`：伴随、对称性、谱、零空间、matvec误差；不只一个PASS标签。
- `CPU_GPU_OPERATOR_COMPARISON.csv`：同快照逐算子误差。
- 三输入同profile的完整复数场、相位与Q结果及独立解析对照。
- `E_ROOT_CAUSE_AND_FIX.md`：Ci／Cj、体积权重和内积如何修复，是否真正达到同一离散关系。

### 15.4 结论约束

当前自然对流只能称部分基准验证，n160未过前不宣称完成分辨率验证。最终论文需要完整剖面、稳态、热边界和分辨率证据；空气方腔不能自动覆盖法国高Pr、变μ隐式和自由滑移顶面。

本轮不以新标签数量衡量进展。成功标志是：方腔完成邻域因果对照并解释异常，电磁得到真实伴随和可复核CPU参考，再判断物理准确性与GPU实现，而不是又提交一批错误组合的失败日志。


<!-- ===== END 215_OPHELIE_CAVITY_NEIGHBOR_AND_EM_ADJOINT_REPAIR_CURSOR_PLAN_2026-11-08.md ===== -->

---

<!-- ===== BEGIN 216_OPHELIE_3D_EM_REFERENCE_AND_CAVITY_TIMESTEP_CURSOR_PLAN_2026-11-09.md ===== -->

<a id="216-ophelie-3d-em-reference-and-cavity-timestep-cursor-plan-2026-11-09-md"></a>

# 216_OPHELIE_3D_EM_REFERENCE_AND_CAVITY_TIMESTEP_CURSOR_PLAN_2026-11-09.md

# Doc216：实际三维电磁参考与方腔有限时间步复核执行计划

日期：2026-11-09  
执行对象：Cursor；项目：冷坩埚招投标／OPHELIE。  
基线：`DOC215_ROUND_PACK_2026-11-09.tar.gz`及对应真实工作树。

## 0. 本轮决策与目标

用户已明确：采用更小时间步获得稳定性是合理的，不因此推翻当前方腔实现。保留Doc215三档通过成果，补有限时间步敏感性即可；主开发精力转向电磁。

本轮必须完成三项：

1. **H：方腔有限复核。** 保留现有CFL策略，最多先加两组代表性减步测试，不继续反复修改浮力／热边界。
2. **E：真实圆柱CPU参考与GPU逐算子对照。** 不再用两节点／四节点玩具矩阵替代真实三维算子验证。
3. **E：解释并修复旋转场被非物理投影的原因。** 根据实际矩阵、谱、弱残差、几何与积分证据选择一项修复，不再扫描det限幅或Neumann系数。

不自动开展TEAM7长跑、法国长时耦合、额外分辨率大扫描、新的RHS／伴随族。P功率后处理可离线补齐，不抢占电磁主线。

继续使用既有壁面粒子；不新增需要额外布置的ghost／镜像粒子。不修改独立A-phi／Laplace论文分支。

## 1. 已确认事实与本轮不再重复的工作

### 1.1 方腔成果

| 量 | n80 | n120 | n160 |
|---|---:|---:|---:|
| Nu误差 | −2.9% | +1.9% | +1.4% |
| 中心线V峰值误差 | −5.4% | −1.7% | −1.4% |
| 中心线U峰值误差 | −8.3% | −7.6% | −6.0% |
| 窗口速度变化 | 0.41% | 0.26% | 0.21% |

三档主要工程指标和窗口稳态已经通过。Nu非严格单调，不强行计算收敛阶；不能说所有局部量都已验证。

H1同快照排查未发现足以解释旧Nu偏差的邻域变化；不重跑整套旧／新邻域矩阵。保留位置更新后热前重建的相容顺序。

### 1.2 方腔时间步事实

新源码实际使用对流CFL：默认0.2；当ν/h相对设定速度尺度较大时自动收紧为0.125。Doc215中n80=0.2，n160=0.125，n120以日志实际值为准。

这是一种数值时间步选择，允许保留。报告需如实说明，不称所有分辨率采用同一数值CFL，也不能单独将改善归因于邻域更新。

### 1.3 电磁成果与缺口

- Ci／Cj伴随公式已修复，是有效进展；不重写第三种伴随。
- 全校正G路径通过常量A与平移不变性，但旋转场精度未过。
- 当前CPU脚本只覆盖2～4节点标量图；未验证实际2400粒子三维圆柱。
- 原脚本M只用σ，未完成体积权重、真实规范、真实谱和真实GPUmatvec对照。
- `det(C)`大不等于矩阵病态。禁止沿用“det到几十所以校正爆炸”的结论。
- 无自感规定A的场对照与自感趋肤验证分开。不得因参数打印δ/R=1就称关闭自感的运行验证了趋肤。

## 2. P0：冻结、复用与组织

建立 `doc216_runs/`：`00_entry/`、`H_dt/`、`E_export/`、`E_cpu/`、`E_gpu/`、`E_weak/`、`E_fix/`、`P_post/`、`scripts/`、`source_snapshot/`。

1. 保存实际HEAD、工作树和所有相关未提交diff，尤其方腔CFL选择、电磁伴随、G校正与约束。
2. 记录二进制哈希和每个run最终解析参数；日志横幅若仍打印旧LHS名称，应新增实际执行路径标识，不能只靠文字解释。
3. 新测试复用现有入口；若扩展CLI，文档明确哪些是新增参数。未知参数报错。
4. 已有数据满足同源码／配置／误差定义条件就复用；不重新运行全部Doc213–215失败组合。
5. H/E可在独立GPU／CPU排程并行；同一共享文件合并、构建和输出写入需串行。
6. 不清理用户未提交工作，不自动push；失败实验代码清楚标识，不混入默认profile。

## 3. H1：最小时间步敏感性矩阵

保持严格DVD物理参数、B＋TVF、ghost OFF、fluid-k以及已通过热前邻域更新顺序。不给β、μ、声速、核函数、h/dp、温度边界或TVF系数做任何同步调整。

### 3.1 先做两组减步运行

| ID | 分辨率 | 既有基线 | 新测试 |
|---|---:|---|---|
| H160-half | 160 | 对流CFL=0.125，t*=1.0 | 显式固定对流CFL=0.0625；同终时与窗口 |
| H80-half | 80 | 对流CFL=0.2，t*=0.70 | 显式固定对流CFL=0.1；同终时与窗口 |

现有`--advection-cfl`可显式覆盖自动策略，使用前核对实际解析。不要同时乘`advect_dt_scale=0.5`造成重复减半。

这里只改变对流时间步系数，声学与热子步仍按相同既有稳定规则计算。实际声学步数、热子步数和累积时间都记录，解释它不是所有dt统一减半。

n120先复用Doc215。只有两端显示明显敏感且需要判断阈值时才补一组n120；不默认扫描更多CFL。

### 3.2 启动、统计和判定

- 可以使用完整且兼容的同状态restart；不能只恢复Position／Velocity而丢失T、rho、历史力或时间状态。
- 若用晚期restart，减步后先经过充分适应再统计，不能立即测原步长留下的状态；没有证据支持restart完整时，从相同初态运行。
- 每Δt*=0.02窗口≥10点，连续3窗；记录Nu、完整中心线U/V与温度、速度场RMS变化和波动。
- 比较相同物理坐标、同一窗口定义；不要按步号对齐。
- 两步长差异目标：Nu均值≤2%，中心线速度曲线相对L2≤3%，无明显稳态恶化；是本项目敏感性目标，不是文献规定。
- 参考误差仍单独报告：Nu≤5%、U/V峰值≤10%、峰位差≤0.05L。

若减步差异小，保留当前较经济策略，H阶段收口。若差异较大，先查实际dt、更新频率、热子步累计和TVF步长依赖；不要按“哪组更接近参考”挑选配置。

只有差异确实较大且仍有收敛趋势时，最多在受影响分辨率补一次中间步长，例如原CFL的0.75倍；不无限加密时间步。出现非单调或更小步长反而大幅改变结果，应报告时间分裂／实现敏感性并定位，而不是继续凑参考。

H的这些补测不阻塞E线，也不撤销Doc215已得到的三档工程验证成果。累计热能量账本如仍缺失，补采样，但不藉此重写流体求解器。

## 4. E1：导出实际三维圆柱算子——本轮第一优先级

### 4.1 固定参考问题

复用Doc215圆柱：R=0.1 m、H=0.6 m、dp=0.02 m、约2400粒子、常σ=16、规定旋转A。主验证先关闭自感、功率标定、J投影、det限幅、moment体积修正。

起点为全校正G＋真伴随，不扩展多套候选。无校正数据只作已有对照。

保留全部实际粒子，不能用去掉边界的内部子集代替同一个圆柱问题。可另用几百粒子小圆柱做全谱检查，但实际2400粒子算子仍必须导出并对照。

### 4.2 导出内容

- 唯一粒子ID、Position、Vol、Mass、σ、核支撑、h／dp、邻接关系；
- 每个有向邻接的真实向量系数a_ij，包含行i自己的Ci和Vol_j；
- 每个粒子完整3×3校正矩阵C及求逆前矩阵，不只导出det；
- 实际Ar/Ai、phi_r/phi_i、Er/Ei、Jr/Ji、Q；
- 实际规范项、源项投影、缩放、PCG内积、预条件器；
- 实际GPU对若干固定向量施加的Gx、伴随y、Kx、b。

### 4.3 明确G的排列

```text
(G phi)_i = sum_j a_ij * (phi_i-phi_j)
G : 3N × N
```

明确向量按粒子三分量交错或按分量分块存储，所有M和A排列一致。重复边合并，保存真实方向；不得在CPU“顺手对称化”邻接后再称同一算子。

优先直接导出稀疏系数，不对实际2400粒子反复调用GPU逐基向量构造密集矩阵。保存CSR／COO及一个能独立读取复算的脚本。

## 5. E2：两层CPU参考，区分实现与离散问题

### 5.1 第一层：忠实复现当前GPU方程

使用导出的同一组系数，用CPU双精度计算当前GPU真正采用的G、权重、规范、RHS和LHS。

当前代码的伴随为普通转置，向量预乘σ，必须核对有无Vi以及求解器使用的内积。先“按实际做了什么”复现，不能直接换成计划中的数学式再比较后说GPU错误。

若Vi全部相等，某些遗漏只表现为整体比例，但规范项的相对尺度可能变化；若Vi不相等，这种等价一般不成立。分别记录，不能把加上Vi保证能修好当前39%场误差。

### 5.2 第二层：明确目标变分方程

```text
Mv = diag(Vi)
Mσ = diag(Vi*σi)，扩展至向量场三分量
F = i*omega*A

K = G^T Mσ G
b = -G^T Mσ F
E = -(G phi+F)
J = sigma*E
Q = 0.5*Re(J dot conj(E))
```

CPU采用普通内积和明确的参考点／零均值约束。固定参考点用一致消元，避免仅改一行破坏对称性。

如果GPU保留体积加权PCG，应对应：

```text
Kv = Mv^(-1)*K
bv = Mv^(-1)*b
```

并验证Mv自伴关系。对称变换到普通PCG也可，但必须明确。禁止混用两种内积。

当前全域`penalty*phi`不只是固定常数，忠实复现中保留并记录；目标参考先采用精确规范。两者有差异时比较其场影响，不把换约束后的变化混叫转置修复。

### 5.3 真转置与检查

CPU直接用G.T，或按原G行系数散射生成转置。不要再手写第三套“校正散度近似伴随”。

检查实际几何的：

- 常数／线性再现；
- 伴随内积恒等式；
- K对称性或MvKv自伴关系；
- 小几何完整谱、实际几何约束后最小特征值／近零模态；
- 当前CPU方程与目标方程的权重和规范差异。

归一化CPU双精度代数目标先取1e-10；若条件数或导出精度受限，报告绝对量、矩阵尺度和误差来源，不靠放宽门槛隐藏实现差异。

只做随机非负二次型不能排除额外零模态。K按G.T构造后对称是结构性质，不是物理准确性证明。

## 6. E3：CPU/GPU逐算子与完整解对照

固定相同粒子、相同导出系数、相同场和约束，依次比较：

1. Gx；
2. G.Ty或相同定义的加权伴随；
3. Kx和b；
4. 实／虚部phi去除允许常数后的差；
5. E、J、Q和同采样功率。

先用固定随机向量，再用线性phi和实际旋转A。CPU双精度参考与GPU单精度的算子相对差目标1e-4量级，最终解的差需结合条件数与线性容差说明。

若算子不同，定位到具体系数、反向邻接、清零／累加、内积、临时字段或GPU同步并修复；不要先改几何或校正矩阵。

若算子相同但解不同，查求解器和约束；预条件器必须适合真实K。不要把LegacyPairwise对角当作新K的准确对角。

若CPU/GPU解一致且场同样错误，明确结束“GPU实现错误”方向，进入E4；不继续调迭代容差。

## 7. E4：旋转场为何被投影——重点诊断

### 7.1 已有原始数据的支持线索

Doc215的无自感R支路，使用全域完整复数J复算：

| 配置 | 全域复数J相对L2 | 同粒子参考功率下降 |
|---|---:|---:|
| 全校正G＋真伴随 | 约39.16% | 约15.33% |
| 无校正G＋真伴随 | 约29.59% | 约8.76% |

误差平方与功率下降近乎一致，符合离散最小二乘移除部分场的行为。它支持“旋转A含非物理离散梯度投影”的假说，但不是单独证明全部装配正确。

### 7.2 能量与投影恒等式

在同一个Mσ内积下计算：

```text
E_num = -(F+Gphi)
r = G^T Mσ (F+Gphi)
2P = ||F+Gphi||²_Mσ
```

无额外罚项、无强制非物理约束且方程收敛时，应满足相应正交条件。若有规范消元，只在允许试探空间检查；若有penalty／约束反力，完整入账。

比较`||Gphi||²_Mσ`、`||F||²_Mσ-||F+Gphi||²_Mσ`及交叉项，判断是在准确投影一个错误的离散空间，还是仍有装配／求解不一致。

禁止仅将功率缩回解析值、删除Gphi或把旋转A特判为phi=0来获得通过。

### 7.3 空间与模态分解

- 输出b_rot=G.T Mσ F在内部、侧壁、端盖、角／棱交界的分布与尺度；
- 与C／RC支路b作对照，不能把常量A表面驱动一起删掉；
- 检查phi、Gphi以及错误电流的空间模态；
- 找出对b有显著响应的近零／低能量模态，区分棋盘模态、边界振荡和光滑物理模态；
- 对线性／二次多项式试探函数检查体积分与边界积分的离散关系，独立表面法向和面积来自几何，避免循环定义。

连续divA=0不意味着任意离散节点b都应被手工置零；需要考察离散分部积分、物理曲面和试探空间的一致性。

## 8. E5：用奇异值与几何证据选择一项修复

### 8.1 先纠正“det大=病态”的判断

导出求逆前矩阵及C的奇异值、条件数κ2、主方向、线性再现残差、支撑数和壁距。示例C=2I，det=8但κ2=1，不是病态。

删除本轮默认候选中的det上限回退；保留旧cap结果作失败证据，不继续扫描cap。

### 8.2 只按以下分支行动

| 诊断证据 | 一项优先修复 | 必须回归 |
|---|---|---|
| CPU/GPU系数／作用不同 | 修索引、转置、权重、同步或约束 | 实际几何逐算子＋三输入 |
| 真正高条件数／近秩缺失 | 查支撑与矩阵构造，处理秩而非det幅值 | 奇异值、线性再现、非零场 |
| 额外零模态主导 | 修已确认的离散稳定性缺陷；先给出模态证据和一致性说明 | 谱、平移不变性、场／能量 |
| 非病态C，主要是曲面求积／离散分部积分缺陷 | 在CPU上设计一项几何一致的局部积分／系数修正 | 几何矩、常量／线性再现、旋转场及变σ |
| 权重／规范影响显著 | 修相应定义，解释量纲和等价变换 | 同快照功率、约束与守恒 |

几何／积分修正必须独立于输入A确定，不能为Arot单独拟合消掉b；验证时保留不参与构造的场作为检查。

如果选择几何积分修正，先在CPU明确其定义、正体积／面积、边界法向、矩条件与适用域；不直接再写新GPU核。不能只把后处理功率权重改好而求解仍用旧权重。

不新增ghost粒子。必要的几何面片或积分点属于积分表示，须清楚说明且不作为新增流体自由度；若不需要则不引入。

一项修复后，先验证原缺陷是否被消除，再做完整物理回归。若证据不支持，返回已缩小的问题，不再平行扩展多套候选。

## 9. E6：同一profile的统一物理回归

同一几何／约束／权重／边界实现下运行：

| 组 | 物理目的 |
|---|---|
| R：旋转A | 正确非零环向电流 |
| C：常量A0 | 封闭绝缘域零电流 |
| RC：旋转A＋A0 | J／Q对矢势选择不变 |
| 整体乘i | 复数实虚路径一致 |
| 原有纯梯度、正的变σ | 基础抵消与变系数回归 |

沿用A0=B0R/2横向输入，不扫描幅值。禁止按case名称自动切路径、解析注入代替求解、J后处理投影和功率归一化。

工程目标：C归一化电流≤1%；R／RC完整复数J差≤1%，Q场差≤2%；R对同位置解析参考完整复数J≤5%。内部／侧壁／端部／全域均报告，不能只用中段幅值。

量化守恒需区分实际方程弱残差、物理边界Jn和旧边通量残差。撤销不相容旧指标时，写明新指标的物理／离散定义；不能取消所有守恒要求。

只有统一profile通过上述检查才恢复圆柱自感：φ、A_ind、J、Q最终状态一致，full_skin参考及一次迭代收紧对照。不把无自感场通过说成自感已通过。

## 10. P：补齐功率与切片口径

优先离线读取已有全域数据：

```text
P_cont   = 与该物理模式匹配的连续解析积分
P_sample = sum Vi*Q_exact(xi)
P_num    = sum Vi*Q_num(xi)
```

无自感使用规定A参考，强趋肤使用full_skin；有限高对无限长差异单列。切片数值和解析使用同一实际体积，名义带厚参考仅用于几何偏差分析。

moment体积修正继续opt-in诊断。不得用恰好匹配r²功率矩来替代电流空间精度，也不直接修改法国流体质量。

## 11. 法国与TEAM7的安排

本轮不自动跑法国长时或TEAM7多配置。保留现有生产DivΣA为受限对照，不继承候选的验证标签。

统一profile通过后，先映射到实际法国入口并确认所有实虚、J/Q、规范和字段相同，再安排一次受控短回归。Robin／辐射、顶面滑移、高Pr与变μ分别标记覆盖范围。

TEAM7若恢复，只固定一个频率和同一个已验证profile；纯源场与导体实验总场参考分开，不做幅值／探针拟合。

方腔有限时间步复核完成后，不继续被电磁拖回去反复改流体。后续法国流体单独验证可在冻结可信热源下开展，但不能把未经验证电磁热源的结果宣称文献复现准确。

## 12. 排程与交付约束

| 优先级 | 执行内容 | 完成条件 |
|---|---|---|
| 首批 | E1实际2400粒子算子导出；H两组有限减步排程 | 有真实G、C、权重及运行日志 |
| 第二批 | 忠实CPU复现＋目标CPU参考＋GPUmatvec | 明确实现差异还是离散误差 |
| 第三批 | 旋转场弱残差、能量、奇异值和模态 | 确定一种有证据的修复 |
| 第四批 | CPU修复＋GPU对照＋统一物理回归 | 同一profile达标或缺口进一步具体化 |

每次运行前写明唯一问题、可能结果和各自下一动作。能用离线矩阵／快照回答的，不长跑。

若仅提交玩具矩阵脚本，本轮E1未完成。若仅继续导出det并尝试cap，本轮E5未完成。不得把“候选暂不升生产”当作替代诊断和修复的交付。

必交文件：

- `ROUND216_REVIEW.md`：H有限复核、E实际参考、修复与剩余问题。
- `H_TIMESTEP_COMPARISON.csv`：实际CFL／dt统计、窗口、曲线差、基准误差。
- `REAL_CYLINDER_OPERATOR.npz`或等价稀疏格式、元数据和CPU复算脚本。
- `CPU_GPU_OPERATOR_AND_SOLUTION.csv`：G／转置／K／b／phi／J／Q比较。
- `CORRECTION_SVD_AND_WEAK_RESIDUAL.csv`：奇异值、条件数、空间与模态结果。
- `ROOT_CAUSE_AND_FIX.md`：证据→定位→修复→回归，假说明确标记。
- 三输入完整复数场与功率表、必要曲线；附真实相关源码、diff、二进制哈希。
- `RUN_INDEX.csv`：实际执行／复用／NOT_RUN分开，不能将建议命令当完成记录。

保留所有原始失败基线；不为得到绿色标签调整参考、指标或分辨率。最终结论应说明电磁缺陷是否已经从实际CPU参考中定位并修复，而不是再次以“表面校正不兼容”笼统收口。


<!-- ===== END 216_OPHELIE_3D_EM_REFERENCE_AND_CAVITY_TIMESTEP_CURSOR_PLAN_2026-11-09.md ===== -->

---

<!-- ===== BEGIN 217_OPHELIE_GENERAL_EM_REPAIR_AND_VALIDATION_CURSOR_PLAN_2026-11-09.md ===== -->

<a id="217-ophelie-general-em-repair-and-validation-cursor-plan-2026-11-09-md"></a>

# 217_OPHELIE_GENERAL_EM_REPAIR_AND_VALIDATION_CURSOR_PLAN_2026-11-09.md

# Doc217｜OPHELIE 通用电磁修复、独立验证与法国工况接入任务书

日期：2026-11-09  
执行对象：Cursor / SPHinXsys OPHELIE 开发分支  
依据：DOC216_ROUND_PACK_2026-11-09 真实结果、源快照与独立算子检查  
优先级：真实残差 → 完整解空间 CPU 修复 → 独立场精度验证 → GPU → 法国入口映射。

## 0. 本轮必须达成的目标

本轮不再把“某个开关让圆柱通过”当作通用电磁模块完成。本轮需要回答：

1. 当前求解残差到底测量了什么，是否真实计算？
2. 线性电势以外的场能否正确求解？
3. 能否在真实三维圆柱上，保持完整粒子电势自由度，同时正确处理旋转场、常量场、二次电势和变电导率？
4. CPU 已验证的同一离散方程，能否在 GPU 上复现？
5. 通过验证的路径是否确实能被法国 case 调用？

方腔保留已有结果，阶段性收口。下一轮不默认新增方腔 CFL 扫描，不重新修改浮力、压力、黏性和热壁面。A-phi/Laplace 独立论文分支不在本轮修改范围。

### 0.1 本轮最低交付与完整交付

| 层级 | 必须交付 | 允许的结论 |
|---|---|---|
| 最低交付 | 残差报告修复、P1 能力边界、完整空间 CPU 候选的真实成败、独立场测试数据 | 可以明确定位尚未修复的环节；不得写通用电磁 PASS |
| 目标交付 | 完整空间 CPU 候选通过全部基本场、变 σ 制造解；GPU 对照通过 | 可称选定算例范围内的通用电势离散已验证 |
| 条件性交付 | 自感与空间分辨率检查、法国入口离线/短时映射 | 可准备法国耦合；不自动称法国文献已复现 |

任务允许出现真实 FAIL，但不允许通过缩小解空间、改宽门限、填零残差把 FAIL 改成 PASS。

## 1. Doc216 已知事实与本轮不重复项

### 1.1 方腔

| 工况 | 原 CFL / Nu | 减半 CFL / Nu | 结论 |
|---|---|---|---|
| n80 | 0.2 / 4.389 | 0.1 / 4.260 | Nu 变化约 −2.95%；半步相对基准 −5.7% |
| n160 | 0.125 / 4.583 | 0.0625 / 4.442 | Nu 变化约 −3.09%；半步相对基准 −1.7% |

Doc215 n80/n120/n160 的既有工程门结果保留。Doc216 表明约 3% 时间步敏感性；不能声称时间独立性已完全证明，也不需要因此撤销所有浮力成果。

创建 `H_STATUS_FROZEN.md`，记录物理参数、CFL、稳态窗口、Nu/速度指标、热边界/TVF/核校正状态及限制。复用原始数据，不重跑。不得以“原 CFL 更贴近参考值”为理由选择时间步；选择依据是稳定、误差可接受和成本可控。

### 1.2 电磁

- 真实圆柱：2400 粒子，R=0.1 m，H=0.6 m，dp=0.02 m，σ=16 S/m，f≈1.5831435 MHz。
- 真实校正 G 已导出，边数 107360；CPU 同一离散方程能够复现原错误。
- CPU/GPU 电流差约 0.045%；CPU 全场电流误差约 39.16%；错误投影移除约 15.33% 功率。
- 校正矩阵最大条件数约 2.10。行列式最大约 31.6，不等于矩阵病态；停止按 determinant cap 扫描。
- P1 限制 φ∈span{1,x,y,z}，本轮 R/C/RC 测试表现好，但不能表达一般非线性电势。
- 变 σ E0 通过的是常量附加场不变性，不是变 σ 的独立场精度验证。
- 自感记录中电势残差为代码填零；不能引用为完整方程收敛证据。
- 生产默认仍为 DivΣA + LegacyPairwise + boundary=none；诊断 P1 尚非生产默认。

### 1.3 已尝试但没有完成修复的路线

| 已有尝试 | Doc216 表现 | 本轮处理 |
|---|---|---|
| 原 GᵀσG 全空间投影 | R 全场误差约39% | 保留失败基线，禁止重复跑一轮后再称新发现 |
| 保留原 K，仅把 RHS 换成几何表面项 | R 好，C 误差约31%，RC 不变性失败 | 不再原样重复；证明只换 RHS 不足 |
| `E5_neumann_divsym_cpu_trial.json` | C/RC 结果有希望；R 求解器 info=800 | 优先读原代码、核对方程与绝对残差，再决定复用 |
| P1 能量投影 | 对本轮低阶输入有效 | 留作对照，不升通用默认 |

`info=800` 不能忽略，但若 RHS 接近零，可能涉及相对残差归一与绝对容差。必须检查 ||b||、||r||、初值与迭代停止条件，不能直接判定物理场错，也不能凭场图好就判求解成功。

## 2. 工作区、源码与路径约束

### 2.1 起点审计

Doc216 记录：

- EM：`feature/EM-Biot-Savart`，HEAD `bc411c43718ca7ecaf4c1ed3e29aefeb4e6090c4`，但有大量未提交改动。
- 流体：`.worktrees/doc210-native`，HEAD `67fef44ef8a2cd6f5f7712bf2e66038cc6de89d3`。
- 最终圆柱二进制 SHA256：`d33fae869f98e432ce9c29225e428471771755cbe304812c753a40764ac7dde5`。

执行前核对当前状态。这些值是历史身份，不是保证当前工作树仍相同。禁止 reset/clean 丢弃改动；新建 worktree 也不能假设未提交代码会自动带过去。

先存 tracked diff、untracked 文件清单及必要源快照；保留原任务书。沿用当前可编译工作区，在 OPHELIE 诊断头、测试入口和独立 CPU 脚本中隔离本轮改动。不要向 main 或远端自动推送，不改独立 A-phi 工作区。

### 2.2 优先检查的文件

实际仓库路径由现有文件定位，不凭本任务书虚构目录：

- `electromagnetic_ophelie_phi_component.h`
- `electromagnetic_ophelie_self_induction.h`
- `electromagnetic_ophelie_doc216_e5_p1_projection.h`
- `electromagnetic_ophelie_doc214_e2_variational.h`
- `test_3d_ophelie_analytic_cylinder_skin.cpp`
- `e2_cpu_real_cylinder_operator_checks.py`
- `e5_cpu_geom_div_neumann.py`
- Doc216 中 `neumann_divsym` 试验脚本；若未打包，定位现有工作树原文件并一并归档。

CPU 结果必须注明输入 NPZ/source hash，GPU 结果必须注明源码和二进制 hash。避免旧脚本配新二进制。

## 3. E-R：先修报告与验收逻辑

### 3.1 删除“未经测量的残差零值”

P1 分支当前将 `phi_*_solver_rel_residual`、`phi_eq_res_vol_*` 直接写成 0。替换成真实指标或明确 `NOT_MEASURED` 状态。

新增明确的状态枚举/字段，不用 0 表示“本方案不提供该量”。如框架需要浮点值，使用状态+数值配对；未测量时不得参加 PASS 判断。不要只填 NaN 而让比较逻辑偶然绕过门限。

### 3.2 分清四种指标

令 F 为某一复数分量的规定驱动，E=−(Gφ+F)，J=σE。

| 指标 | 定义/意义 | 不允许的替代 |
|---|---|---|
| 线性求解残差 | r=Kφ−b；记录绝对量和相对量 | 不用剩余电场大小代替 |
| P1 驻点残差 | BᵀGᵀMσ(GBa+F)，B 为基函数矩阵 | 不能称完整空间守恒残差 |
| 完整离散弱缺陷 | GᵀMσ(Gφ+F)，注明所用 G 与度量 | 该 G 自身有缺陷时也不是独立物理真值 |
| 独立物理缺陷 | ∇·J、边界 J·n、闭合区域通量，独立重构 | 不用生产同一残差自证全部守恒 |

P1 现有 `residual_rel = ||Gφ+F||σ/||F||σ` 改名为 `remaining_field_ratio`。正确旋转电场中可接近1；它不是求解残差。旧日志字段保留映射说明，避免跨轮误读。

对零/近零 RHS，采用 `||r|| ≤ atol + rtol*||b||`；atol 依据无量纲尺度和浮点精度设置并记录，不为单个测试临时调宽。所有求解器 info、迭代数、停止原因完整输出。

### 3.3 统一体积权重与规范约束

连续目标为 ∫σ|∇φ+F|² dV。粒子求积必须明确 Mσ=diag(Viσi)；当前等体积会掩盖漏 Vi 问题。P1 主机循环当前仅使用 σ，修正或显式限制为等体积专用。

先以“固定一个参考自由度”或“体积加权零均值”精确处理常数零空间。不要用全域 `+penalty*φ` 后仍称只是固定常数规范；它会改变非常数模态。若为复现旧 GPU 方程保留 penalty，必须标为 legacy reproduction。

P1 只需解三个斜率，常数由规范确定，避免用极小常数钉扎混淆条件数。通用非对称 Neumann 系统不得未经确认就使用 CG；按矩阵性质选稀疏直接解/GMRES等，并记录真实残差。

### 3.4 E-R 验收

重放一组 R、一组 C 和一组自感，确认：无硬编码残差；近零 RHS 正常结束；报告能区别 reduced-space/full-space；旧 `passed=1` 不会自动继承新“完整验证”状态。

## 4. E-M：用固定测试集约束修复，禁止只看圆柱旋转场

### 4.1 统一符号和误差口径

对实部、虚部分别定义 F，符号必须从现有 E=−∇φ−∂tA 推导并在代码注释中固定；不得为了匹配某张图改变符号。

复数全场电流误差：

`eJ = sqrt(sum_i Vi (|Jr−Jr*|²+|Ji−Ji*|²) / sum_i Vi (|Jr*|²+|Ji*|²))`。

零参考电流时以规定的 `Jchar=σref Fchar` 归一，不除以接近零的真值。φ 先去同一规范，再比较；主要物理量为梯度、J、Q。

同时输出全域/内部/侧壁/端面分区；径向平均、|Jθ|相关系数只作辅助，不替代完整复数向量误差。JouleHeat 延续现有峰值/RMS约定，若使用峰值相量 Q=|J|²/(2σ)，参考也使用同一约定。

### 4.2 固定的基本测试

| ID | 输入与精确性质 | 检验目的 |
|---|---|---|
| M0 | φ常数、线性坐标 | 常数零空间、线性重现 |
| M1/C | F常量，φ*=−F·x，J*=0 | 线性电势、绝缘边界、零电流 |
| M2/R | 既有轴对称旋转F，φ*=0 | 不错误抵消真实环向场 |
| M3/RC | 同一R叠加同一C | 常量附加场不变性 |
| M4 | ψ=ψ0[(x−cx)²−(y−cy)²]/R²，F=−∇ψ，φ*=ψ，J*=0 | 二次调和电势、非线性规范抵消 |
| M5 | 非零环向J* + 上述二次ψ，构造F=−∇ψ−J*/σ | 同时检验非线性φ与真实电流 |
| M6 | M5改为正的空间变σ | 变σ独立场精度 |

M4 是纯梯度附加场；不能因它不是原线圈输入就删除。它测试电势算法必须满足的梯度抵消能力。它与常量附加场测试是同一逻辑的更高阶检查。

M5/M6 采用圆柱轴 z：`J*=J0*(-(y−cy)/R, (x−cx)/R, 0)`。该场散度为零，且对圆柱侧面及端面无穿透。σ可取 Doc216 的 `σ0*(1.5+0.5*(x−cx)/R)`，在域内保持正值。于是 φ*=ψ、J=J* 是独立已知解，且 ψ 不在线性 P1 空间内。

ψ0、J0 通过固定尺度选取，使 |∇ψ| 与 |J*/σ| 同量级；写入配置，所有候选复用。禁止某项过小而让错误被主场掩盖。

### 4.3 两层制造解不得混淆

1. **离散代数测试**：F=−Gψ，检查算法是否能在其完整解空间内抵消已知离散梯度。用于找实现/子空间限制，不能代替连续精度验证。
2. **连续制造解测试**：在真实坐标上用解析 ∇ψ、σ、J* 构造 F。生产重构近似会产生空间误差，应随分辨率检查。

独立复核已发现：当前 P1 对 M4 的连续与离散输入均几乎保留全部伪电流，剩余比例约100%。本轮复现并归档，不把它作为新问题再做参数扫描。

### 4.4 本轮门限（固定后不得追着结果改）

- CPU代数：常数/伴随/矩阵装配等目标1e-10量级，说明无量纲归一；单精度导出系数的连续线性重现另报，不混为双精度误差。
- M1、M3：C电流≤1% Jchar；RC−R电流≤1%，Q≤2%。
- M2：全场复数向量电流误差≤5%。
- M4：零电流误差≤5% Jchar，且细化不恶化；离散精确梯度版本需接近线性求解精度。
- M5/M6：全场eJ≤5%，同粒子采样总Q误差≤5%；φ梯度误差单列，不能只比总功率。
- 若dp=.02不能满足空间门，允许一次dp=.01细化；必须报告两档并区分“空间未收敛”和“离散结构错误”。

## 5. E-N：完整电势空间 CPU 修复——本轮主任务

### 5.1 先复用现有 Neumann 试验，避免重复建设

读取 Doc216 `neumann_divsym` 原代码，写一页方程审计，明确：

- 未知量是否仍每粒子一个φ；
- 内部算子到底是哪个散度/梯度/拉普拉斯组合；
- RHS与边界是否从同一方程推导；
- 表面法向、面积、边界评估位置如何给出；
- 是否缺失常量零空间约束；
- `info=800` 是否因近零 RHS 的相对容差；
- 是否只输出场值、未计算边界残差。

先修审计确定的实现问题，跑 M1/M2/M3/M4。若结构成立再扩展 M5/M6，不能只修到R通过就停止。

### 5.2 明确目标方程

求解：

`div[σ(grad φ + F)] = 0` in Ω；

`n·σ(grad φ + F) = 0` on ∂Ω。

等价强形式为 `div(σ gradφ)=−div(σF)`，边界 `n·gradφ=−n·F`（σ正）。

变σ时包含 `gradσ·gradφ` 与 `gradσ·F`。禁止沿用常σ时 `divF=0` 就把变σ体源设零。禁止假设任意 coil A 都满足当前离散 divA=0。

### 5.3 如果现有试验没有一致构造，允许一个明确替代候选

CPU 实现局部多项式重构的完整空间 Neumann 配点原型，作为本轮唯一备用离散候选：

1. 内部点使用尺度归一的二次局部拟合，输出一阶/二阶导数；邻域秩、奇异值及条件数必须记录。
2. 拟合得到 `σ Δφ + gradσ·gradφ` 的线性行，RHS为 `−div(σF)`。
3. 对制造解先用解析RHS隔离LHS误差，再用和真实A相同的数值RHS构造。只有数值RHS版本通过才允许迁移生产。
4. 边界行在真实圆柱侧面/端面位置评估法向导数。粒子中心通常距真实壁面半个dp；不得把近壁粒子中心当作真实边界而不说明误差。
5. 可用现有近壁粒子的局部拟合在真实边界点外推导数；边界评估点不是新增流体ghost粒子，不参与动力学。记录其位置、法向和与未知量的映射。
6. 几何转角处避免用一个随意平均法向代替两个面的约束；先对侧面/端面各自评估，再明确组装方式和方程计数。
7. 求解器与矩阵性质匹配；严格处理常数零空间及Neumann相容性。兼容性修正大小必须输出，不能大幅投影RHS后隐藏物理缺陷。
8. 用同一φ的局部梯度重构E/J/Q；如果还用旧G重构，必须把“求解梯度”和“输出梯度”的差异单独验证，不能混用后声称同一离散通过。

这是待验证候选，不是已证明的修复。完整空间强式配点也不自动保证局部守恒，必须通过第6节独立通量检查。不要为了按期交付跳过该限制。

### 5.4 推进分支与停止条件

顺序固定：已有Neumann试验审计/修复 → 基本场 → 制造解 → 独立守恒 → 一次细化。

若首选在基本场失败，最多修复两个有明确证据的结构问题；不得新增多个边界开关组合扫描。如果备用完整空间候选也失败，提交真实失败模式和最小复现，停止GPU迁移。

此时可制作一个独立CPU有限元/有限体积参考解，求解同一标量电势方程，用于辨别粒子离散与几何/参考数据的问题；它属于验证工具，不迁入A-phi论文分支，也不意味着立即更换整个工程框架。不得临时造一个低阶全域电势空间冒充独立高精度参考。

## 6. E-C：独立连续性、边界与热源检查

### 6.1 必须补测的量

- 圆柱真实侧面及上下端面的J·n，面积加权L2与最大值分列。
- 完整外边界净电流，归一为 `|∫J·n dS|/(Jchar*S)`。
- 内部至少8个固定子区域的闭合通量，避免只看全域正负抵消。
- 内部独立重构散度，使用不同于求解残差的评估；局部拟合秩不足点标记无效，覆盖率必须输出。
- 独立结果先用解析J*验证测量本身；否则不得用测量噪声判求解器错误。

本轮工程目标：边界法向L2/Jchar≤5%，全域净通量≤1%；局部闭合通量以子区域表面积和Jchar归一，目标≤5%。当参考测量本身未达精度时，标为MEASUREMENT_UNRESOLVED，不自动PASS，也不据此立即改求解器。

### 6.2 功率拆分，避免几何采样误差混入电磁误差

同时输出：

- `P_cont`：连续解析/独立参考积分；
- `P_sample`：在同一粒子、同一Vi、同一分区上采样解析Q后求和；
- `P_num`：当前数值Q积分。

分别报告 P_sample/P_cont 的几何求积误差、P_num/P_sample 的电磁场误差。总功率接近不代表空间分布准确。

M4/零电流以参考Jchar构造功率尺度，禁止用接近零的参考功率作分母。

## 7. E-G：CPU通过后才迁GPU

### 7.1 迁移要求

保留同一未知量空间、边界行、RHS、规范约束和重构方法。CPU/GPU都输出Gx/算子作用、RHS、φ、J、Q抽样，不只比较最终总功率。

如果使用非对称算子，不得为复用PCG把它未经说明改成正规方程；这会改变条件数和残差含义。

### 7.2 测试集合与容差

固定M1–M6及复数相位旋转。CPU/GPU全场J相对差目标≤0.1%，Q≤0.2%；近零场用Jchar归一。困难矩阵若达不到，必须输出真实条件/迭代证据，不仅提高最大迭代数后给PASS。

对同一复数输入整体相位旋转，J应同相位旋转、Q保持不变。禁止只测虚部，再认为复数两分量都已验证。

### 7.3 GPU性能定位

Doc216 P1会同步到主机并用Eigen解小系统，不应宣传为全GPU求解。本轮CPU原型可以保留；若GPU通用迁移未完成，明确标出CPU/GPU数据同步与求解时间。用户接受准确优先，不能为速度保留已知错误。

## 8. E-S：自感结果的真实收口

### 8.1 先修残差，再复核原P1圆柱作为对照

原P1自感不是本轮主修复路径，但其积极结果应保留。记录原始数据：中部径向J指标约0.96%，总功率误差约6.7%，固定点残差约3.4%。这些不是完整三维场误差的替代。

原门限 `rel_P<0.40` 仅是宽松烟雾门，不用于“高精度电磁验证”结论。不要修改旧历史记录，新增严格门和状态字段。

### 8.2 只做两次有目的的迭代收紧

在选定算法通过E-G后，同一dp=.02、自感松弛0.3：

1. j_tol=.01，最多40步；
2. j_tol=.003，最多80步。

比较未松弛固定点残差、最终场一致性、J/Q/功率变化。若两次J差≤1%、功率差≤1%，可认为该范围内迭代误差已受控；否则报告未收敛，不默认继续无限加步。正确电势驻点/方程指标同时输出。

Picard残差达标只说明迭代自洽，不能代替解析场精度。

### 8.3 局部热源比异常必须解释

Doc216 `Q_o/c_sph=39.3326`、`Q_o/c_an=25.3471`，差约55%。先定位定义：点值还是bin均值，中心是否接近Q=0，是否同坐标/同bin/同体积。

若指标对中心小分母敏感，保留原值并新增固定有限体积区域平均比，不删除原异常。用同采样解析Q比较，排除粒子位置与bin口径误差后，再决定是否需细化。主场应同时比较全场Q误差和边界加热层厚度，不只看相关系数。

仅当证据指向空间分辨率，增加一次dp=.01；同时记录实际delta/dp。当前delta/R=1、dp=.02仅约5个粒距跨越delta，不能代表更强趋肤工况已验证。

## 9. F-M：法国入口映射，先离线验证再短时耦合

只有E-G通过，才把新路径接入法国case；P1不可仅凭圆柱通过就升默认。

### 9.1 显式profile，防止路径错配

建立可打印的profile，命名沿用仓库约定即可，至少区分：

- 历史production（保留复现）；
- P1受限对照；
- 本轮验证的完整空间候选。

打印φ空间、梯度/散度/边界、规范、复数分量、σ来源、自感开关、相量约定、功率归一和代码hash。开关冲突直接报错，不静默改开关。

### 9.2 法国几何离线电磁快照

沿用已确认玻璃直径500mm及现有高度/线圈几何，逐项打印来源。自然对流与搅拌轴向不得混用。

先固定温度、常σ=16，导出A/φ/J/Q和守恒指标；再使用明确的非均匀温度快照映射σ(T)，做同样检查。温度快照可为诊断场，必须标注，不作为文献真实温度数据。

均匀σ轴对称算例下P1可能接近正确，可作为对照；若P1与完整空间结果接近，只能说明本快照近似可用，不能推及非对称搅拌和任意σ(T)。

### 9.3 热源归一前后都输出

文献自然对流对照优先沿用50kW；如当前甲方配置另有60kW，单列配置，不能混作同一文献工况。禁止靠目标功率归一掩盖局部Q错误。

输出原始P、归一系数、归一后P、相同归一规则下的空间Q分布。常量附加场/相位旋转前后用同一归一协议，优先先比较未归一结果。

### 9.4 短时耦合允许做什么

仅验证EM→Q→热模块的数据传递、更新频率、变量身份和能量收支，不以几秒环流图宣布文献复现。默认不启动新的长时法国自然对流，也不把调物性凑流线插入本轮通用EM修复。

TEAM7先保留既有状态；完整空间/边界修复通过后再做针对性迁移，不另开大规模TEAM7参数扫描。

## 10. 执行顺序与退出条件

| 阶段 | 工作 | 进入下一阶段条件 |
|---|---|---|
| A | 起点审计、H冻结、E-R报告修复 | 残差真实，源码身份明确 |
| B | 固定M0–M6，复现P1能力边界 | 测量/解析构造可独立检查 |
| C | 审计并修复已有Neumann CPU试验 | 基本场与非线性制造解通过 |
| D | 必要时一个完整空间备用候选 | 全场精度与独立边界/通量通过 |
| E | CPU→GPU同方程对照 | M1–M6及复数相位通过 |
| F | 自感收紧、局部Q及一次必要细化 | 迭代误差与空间误差分清 |
| G | 法国离线快照/短时映射 | 实际入口确实使用验证路径 |

优先完成A–D，避免还没确认方程就急着移GPU。H冻结文件可与CPU工作独立准备，但不派生新的流体修复任务。

每个失败必须给出：最小输入、期望解、实际误差、哪条证据支持哪类问题、只修改了什么。最多保留首选与一个备用候选；不增加第三/第四组经验组合。

## 11. 交付目录与报告格式

建议根目录：`docs/ophelie/reference/doc217_runs/`。

| 目录 | 必须内容 |
|---|---|
| `00_entry/` | 代码身份、源/二进制hash、完整运行命令、状态表 |
| `H_frozen/` | H_STATUS_FROZEN.md，原结果引用与限制 |
| `E_reporting/` | 残差语义、填零修复diff、重放日志 |
| `E_manufactured/` | M0–M6配置、解析定义、粒子CSV、误差表 |
| `E_cpu/` | 已有试验审计、完整空间矩阵/系数、CPU脚本、失败基线 |
| `E_conservation/` | 真实边界采样、独立测量自检、局部通量与覆盖率 |
| `E_gpu/` | 同算子对照、相位测试、性能与同步说明 |
| `E_selfind/` | 迭代收紧、Q异常复核、必要细化 |
| `F_mapping/` | 法国profile映射、离线快照、可选短时能量收支 |

根目录报告必须包括：

- `ROUND217_CLOSURE.md`：分别说实现完成、代数通过、物理精度通过、生产接入完成；禁止合成一个笼统PASS。
- `RESULTS_TABLE.csv`：测试ID/算法/空间/σ/自感/分辨率/误差/阈值/状态/hash。
- `ROOT_CAUSE_AND_FIX.md`：哪些旧判断被修正，哪些是真正证据，哪些仍是假设。
- `REMAINING_LIMITATIONS.md`：不能把未做的复杂几何、强趋肤、TEAM7写成通过。
- `NEXT_DECISION.md`：下一步只给一个主方向，附失败时明确分支。

数据较大可以压缩稀疏矩阵，但不能只留漂亮汇总，删除能够独立重放的真实算子、几何与配置。脚本/配置相对路径可重放，不依赖某个旧工作树的绝对路径。

## 12. Cursor完成前的检查清单

- [ ] 没有在P1/其他新路径硬编码方程残差为0。
- [ ] `remaining_field_ratio` 与线性求解残差已分开。
- [ ] 规范约束没有被全域惩罚项偷偷替代。
- [ ] 体积权重、相量约定、复数符号已统一。
- [ ] M4/M5/M6的正确解来自独立解析构造，而非候选算子的输出。
- [ ] 原样失败的几何RHS尝试没有重复包装成新修复。
- [ ] 完整空间方案没有在内部退化为全域P1。
- [ ] 近零RHS停止条件同时检查绝对量。
- [ ] 独立边界测量先通过解析场自检。
- [ ] CPU通过后才迁移GPU，迁移的是同一个方程和边界。
- [ ] 径向平均好看没有替代全场复数误差。
- [ ] 功率比较区分连续积分、同粒子采样和数值误差。
- [ ] 自感PASS不依赖填零φ残差或40%宽松功率门。
- [ ] 法国case入口打印的profile与验证profile一致。
- [ ] 方腔未被无证据重新打开，独立A-phi分支未被修改。

本轮收口原则：保留方腔与对称圆柱的真实成果，优先补上“非线性电势 + 变σ + 真实边界 + 完整空间”这条证据链，再向法国自然对流继续推进。


<!-- ===== END 217_OPHELIE_GENERAL_EM_REPAIR_AND_VALIDATION_CURSOR_PLAN_2026-11-09.md ===== -->

---

<!-- ===== BEGIN 218_OPHELIE_D2_BOUNDARY_COMPATIBILITY_AND_FRENCH_THERMAL_CURSOR_PLAN_2026-11-10.md ===== -->

<a id="218-ophelie-d2-boundary-compatibility-and-french-thermal-cursor-plan-2026-11-10-md"></a>

# 218_OPHELIE_D2_BOUNDARY_COMPATIBILITY_AND_FRENCH_THERMAL_CURSOR_PLAN_2026-11-10.md

# Doc218｜D2 边界与相容性修复、C++ 验证闭环及法国热边界接入

日期：2026-11-10  
执行对象：Cursor / SPHinXsys OPHELIE 分支  
基线：DOC217_ROUND_PACK_2026-11-10  
本轮原则：修两个已定位问题，补齐实际 C++ 验证，通过后回到法国热边界；不再扩展候选算法。

## 0. 请先读：本轮究竟要完成什么

Doc217 的 D2 局部二次 Neumann 配点方案已经在 Python 制造解上取得进展，自感迭代也基本稳定。本轮继续这一方案，不重写电磁框架，不回退到全域 P1，不重开方腔修复。

当前优先修复的两项是：

1. **C++ 边界方程左右两侧评估位置不一致**：电势法向导数在真实壁面，驱动场 F 却取近壁粒子中心值。
2. **边界 RHS 算术减均值不成立**：它可能改变有非零体源的正确 Neumann 问题。

必须把 Python 成果落实到实际 C++ D2 路径，不能继续用“同一离散家族”代替实际等价性验证。

### 0.1 执行顺序

| 顺序 | 任务 | 完成标准 |
|---|---|---|
| A | 固化起点、重现两个缺陷 | 旧行为有最小反例与完整命令 |
| B | 统一壁面 F、取消无依据相容性处理 | 基本场和新增非零体源测试通过 |
| C | 实际 C++ M1–M8 与 Python 对照 | 同输入、同矩阵/重构、真实残差可复核 |
| D | 数值 RHS 下独立通量与功率口径 | 不再用解析 RHS 结果覆盖实际路径 |
| E | 自感与法国离线回归 | 修复没有破坏已有有效结果 |
| F | 法国固定热源 + Robin/辐射短时验证 | 能量收支、温度响应、归一协议明确 |

A–D 是本轮必须完成的主任务。E 是必要回归；F 在 A–D 通过后执行。不得为了完成 F 而绕过前面的失败。

### 0.2 明确不做的事

- 不重新扫描方腔 CFL、不修改流体黎曼/压力/黏性/TVF。
- 不开发第三种电势离散，不做 determinant cap/边界权重组合扫描。
- 不默认启动法国长时自然对流或大规模 TEAM7 迁移。
- 不修改独立 A-phi/Laplace 论文分支。
- 不要求本轮必须完成全 GPU Krylov；CPU 准确性优先，但须真实报告成本。
- 不新增流体 ghost 粒子；真实壁面评估点只用于电磁边界与积分。

## 1. 基线事实与结论边界

### 1.1 可以继承的成果

- 方腔 Doc215/216 工程验证状态保留，约3%时间步敏感性记录保留。
- D2 Python 数值 RHS：M6 电流误差约0.76%、Q误差约0.97%。
- Doc217 D2 自感：j_tol .01→.003，径向J差约0.051%、功率差约0.091%。
- 法国 Ø500mm 离线入口、一次无扩散集总加热交接已跑通。
- 残差填零问题已有修复基础，不能重新引入。

以上均按已完成范围继承，不能推导成通用三维电磁、全设备求解或法国自然对流复现完成。

### 1.2 独立复核发现（不是实际 C++ 重跑结果）

使用 Doc216 真实2400粒子圆柱 NPZ、Doc217 D2 Python 装配进行代数复核：

| 测试 | 设置 | 电流相对误差 |
|---|---|---:|
| 调和二次电势 | F取真实壁面 | 约2.6e-10 |
| 同一问题 | F改取粒子中心，模拟当前C++边界 | 约9.9% |
| 非调和二次电势 ψ=x² | 关闭边界减均值 | 约2.6e-9 |
| 同一问题 | 开启边界减均值 | 约48.4（无量纲，约4840%） |

这些数值用于说明缺陷，不是本轮要求精确复制的浮点指纹。必须在本地相同输入上重现趋势，并明确 Python 诊断与实际 C++ 验证的区别。

### 1.3 热源比疑点的已知解释

Doc217 无自感 R：数值外侧/中心Q比约39.22215；在相同粒子与体积权重上采样解析 Q∝r²，约39.22239。两者几乎一致。

因此，无自感情形的39对连续环带比25，主要是采样口径差，不应据此继续修改求解器。自感仍需用同样的同粒子解析对照补齐证据。

## 2. 起点、文件与实现隔离

### 2.1 起点保护

历史主树 `feature/EM-Biot-Savart` HEAD 为 `bc411c43718ca7ecaf4c1ed3e29aefeb4e6090c4`，但 Doc214–217 有大量未提交改动；该 HEAD 本身不是本轮完整源代码身份。

开始前输出当前 git status、tracked diff、必要 untracked 源文件清单和源/二进制 hash。不能 reset/clean 丢弃改动；新建 worktree 不会自动继承未提交代码。不得自动 commit/push。

保留 Doc217 原始数据和任务书；本轮写 `doc218_runs/`。历史 production/P1 只作对照，D2 保持显式 opt-in，未通过本轮验收前不改默认。

### 2.2 必须检查的代码

在实际仓库定位，以下均为文件名而非承诺的绝对路径：

- `electromagnetic_ophelie_doc217_d2_local_poly.h`
- `electromagnetic_ophelie_phi_component.h`
- `electromagnetic_ophelie_self_induction.h`
- `test_3d_ophelie_analytic_cylinder_skin.cpp`
- 法国 `natural_em`、`complex_joule_to_heat_one_way` 入口
- `d2_cpu_local_poly_deepfix.py`
- `ec_cpu_local_poly_conservation.py`

优先修改 D2 隔离实现和对应测试，不重构不相关公共流体模块。新增测试放入现有圆柱入口的独立模式；若必须新建 executable，使用既有 `tests/extra_source_and_tests/3d_examples` 结构，不复制整套生产代码。

## 3. A：先保存最小失败，再修复

### 3.1 保留两条旧行为诊断

只在诊断测试中允许切换：

- `wall_F_mode=center_legacy / wall_reconstructed`；
- `compat_mode=boundary_mean_legacy / none`。

名称为建议，实际实现需写清CLI及默认。legacy模式不得流入正常法国生产入口。正式回归固定 `wall_reconstructed` + `none`。

先记录：输入坐标、σ、F中心值、真实壁点、法向、F壁值、装配行、缩放、规范约束、解和完整电流误差。禁止只提交截图。

### 3.2 重现后不继续扫参数

沿用 dp=.02、圆柱R=.1/H=.6、h0=2.5dp、hmax=4dp；不同时改变粒子云、邻域、精度、正则化和边界处理。先单因素证明问题，再修复。

## 4. B1：统一真实壁面上的驱动场评估

### 4.1 方程与位置必须统一

对一个复数分量，约定：

`F = a_sign*omega*A`；`E = -(gradφ + F)`；`J=σE`。

真实边界方程：

`n(xw)·gradφ(xw) = -n(xw)·F(xw)`。

当前C++是左侧xw、右侧xi，必须修正。实部与虚部分别沿用现有符号，不因某个测试改符号。

### 4.2 推荐实现：共用局部拟合的值评估

给定粒子i的实际成功拟合状态（ids、h、系数映射），对每个驱动分量求局部多项式系数c，在真实壁点使用 `B((xw-xi)/h)·c` 评估F。

- 壁面驱动值与电势导数使用同一个真实位置和法向。
- F的值评估与导数评估使用统一的自适应拟合状态，避免一处用h0、一处用hmax却不记录。
- 常量/线性F应能重现；一般非多项式F报告外推误差，不能调用制造解解析函数偷偷替代生产数值评估。
- 允许测试模式额外使用解析F(xw)，用于隔离壁面重构误差；最终通过必须来自数值壁值路径。
- 规范去均值只作用φ，不作用F。

### 4.3 统一装配、∇σ与重构的失败策略

当前装配会逐级扩邻域；∇σ/输出重构却可能只按邻居数量决定是否扩邻域，并在拟合失败时默认为0。统一为同一 `adaptiveFit()` 策略，并缓存实际成功的状态。

以下情况必须显式失败退出或返回失败状态，不能继续输出PASS：

- 内部/边界装配行拟合失败；
- ∇σ拟合失败；
- F壁值拟合失败；
- E/J/Q重构拟合失败；
- SparseLU分解失败、求解失败或解非有限。

删除 `1e-12*φ_i=0` 的物理行兜底；它仍等价于钉住该点电势，不是无害的软Neumann。常数规范的一个参考自由度单独标识。

输出真实邻域最小数、最大/分位条件数、扩邻域计数、失败位置；不能只报中位数和 `rank_fail=0`。

## 5. B2：取消错误减均值，正确检查相容性

### 5.1 为什么不能减边界算术平均

求解：

`div(σ gradφ) = f = -div(σF)`；

`σ n·gradφ = g = -σ n·F`。

连续相容性为 `∫Ω f dV = ∫∂Ω g dS`，不是 `mean(boundary_scaled_rhs)=0`。即使σ恒定，只要f非零，边界通量也可以非零。

当前边界行经过h缩放，节点分布不等面积。对这些行减算术平均会改变边界数据。正式路径默认不进行该操作，保留旧值仅用于诊断日志。

### 5.2 先在不改RHS的情况下求解

保留原始物理K、b，再加显式规范约束。记录：

- K·1（右常数零空间缺陷）；
- 求解后的原始物理行残差，包括被pin替换的那一行；
- pin前后矩阵/方程身份；
- 未缩放的内部残差与边界残差，分别归一；
- 物理积分相容性测量及其测量误差。

稀疏配点K一般非对称，`K·1≈0`不意味着`1ᵀK≈0`。不得把右零向量当作左零向量。若需要诊断离散相容性，可计算/估计左零向量y并检查yᵀb，但不默认投影或修改物理数据。

先修造成不相容的源项/边界位置/重构问题。若仍需离散兼容处理，先证明其合理性和修正大小，单列为未解决项；本轮禁止用一个新的经验系数代替旧减均值。

### 5.3 规范与残差报告

固定参考φ或体积零均值只用于常数规范。若先pin再demean，代数残差在pin解上计算，物理残差在最终输出φ上计算；二者都报告。

真实输出 `||r||`、`||b||`、`atol`、`rtol`、停止判据。使用绝对+相对容差处理近零RHS。当前带分母下限的 `linear_rel` 应另名为有效归一残差，不能让用户误认为始终等于 `||r||/||b||`。

CPU主机局部矩阵和SparseLU建议使用double，GPU粒子字段仍可Real；若做此调整，先分开记录“边界修复”和“精度变化”的影响。不要边改精度边宣称所有改善只来自边界。

## 6. C：实际 C++ 统一回归，不再只在Python验证

### 6.1 已有M1–M6全部沿用

| ID | 精确性质 | 关键验收 |
|---|---|---|
| M1/C | 常量F，线性φ抵消，J=0 | J/Jchar≤1% |
| M2/R | 环向旋转F，φ=0 | 完整向量电流误差≤5% |
| M3/RC | R+C，J与R相同 | ΔJ≤1%，ΔQ≤2% |
| M4 | ψ∝x²−y²，F=−∇ψ，J=0 | 数值壁值路径J/Jchar≤5% |
| M5 | 二次ψ + 散度零环向J* | eJ≤5%，Q场误差≤5% |
| M6 | M5改正的非均匀σ | 数值体源、数值壁值、数值∇σ全部启用后同上 |

坐标均相对圆柱中心。保留Doc217固定幅值，不通过让梯度项过小来掩盖边界错误。

### 6.2 新增且必须执行的两个非零体源测试

新增M7：`ψ=ψ0*((x-cx)/R)^2`；新增M8：`ψ=ψ0*((z-cz)/H)^2`，H为半高。

两者均取 `φ*=ψ`、`F=−∇ψ`、σ=16、J*=0。体源由方程得到，不强制设0：

- M7：`f=σ*2ψ0/R²`；
- M8：`f=σ*2ψ0/H²`。

M7暴露非零体源与侧面通量，M8检查上下端面边界。全部使用数值体源与数值壁值作最终门；解析版本仅作隔离对照。零电流误差目标≤5% Jchar，真实壁面法向误差≤5%。

侧壁/端面交界不能只报“side优先”就认为端面也已验证。边缘附近按每个物理面分别评估独立残差。如果端面漏约束导致失败，针对该缺陷修复，不扩大成新算法。

### 6.3 测试必须走真实调用链

在现有圆柱入口加入制造解模式，调用生产D2同一个装配、求解和重构函数。不得另写一个仅供测试的“正确C++求解器”。

每个测试至少导出pos/Vol/σ/F/φ/J/Q/真实壁点与法向。给出模式冲突检查，制造解模式不得同时保留线圈驱动或自感迭代。

Python与C++使用同一输入CSV/NPZ；比较K作用、b、规范后的φ、J/Q和壁值。旧Python即使写“numeric RHS”，壁值仍可能使用解析函数；本轮必须显式区分 `volume_rhs_mode`、`wall_F_mode`、`grad_sigma_mode`。

同为double时算子/解差目标1e-8量级；涉及单精度字段时全场J差≤0.1%、Q≤0.2%作为工程对照门。超标给出来源，不把松门限当修复。

### 6.4 复数相位回归

对非零复数输入执行一次整体相位旋转：J应同相位旋转，Q应不变。零电流用固定Jchar归一。不能只测虚部就称复数实现完成。

### 6.5 分辨率与测试预算

先完成同一dp=.02回归。若误差有明确空间特征，最多增加一次dp=.01；二次多项式重现本应很好，若仍大误差，优先查位置/符号/装配而不是立即细化。

所有PASS必须带字段：测试ID、真实算法、数值/解析源项、壁值模式、σ模式、精度、dp、hash、各指标及阈值。

## 7. D：补独立通量与功率证据

### 7.1 先改报告范围

Doc217守恒结果仅对应解析RHS、同类拟合测量；保留原数据，改写为“解析制造解下的初步守恒检查”。本轮至少对M2/M5/M6/M7/M8的数值路径重新评估。

### 7.2 外边界积分

使用与粒子分布独立的圆柱几何积分点：侧面按角度×高度布点，端面按极坐标或三角面元积分，面积权重之和核对解析面积。

在这些点从数值φ/F/σ重构J；制造解测量另用解析J*直接评估，检查测量自身误差与覆盖率。不能在测量时用解析F替换数值解应有的F而不标注。

输出每个面 `L2(Jn)/Jchar`、max、净电流和面积覆盖率。净通量用Jchar×实际总面积归一；目标L2≤5%、净通量≤1%。覆盖率不足95%标记MEASUREMENT_UNRESOLVED，不能默认为0。

### 7.3 子区域必须真正做表面积分

沿用8个固定子区域，每个子区域包含其真实外侧面、内部切面及端面部分；逐面输出法向/面积/通量。

分别计算 `∮J·n dS` 与独立 `∫divJ dV`，比较两者。不能再把后者命名成已独立测得的表面通量。子区域归一用其自身表面积，不用整个容器表面积稀释误差。目标闭合通量≤5%，测量先在解析J*上自检。

### 7.4 功率与热源比

明确三种量：

- P_cont：独立连续积分/高精度几何积分；
- P_sample：同粒子、同Vi采样解析Q；
- P_num：数值Q积分。

例如均匀σ、J*=J0*(-y/R,x/R,0)时，峰值相量约定下 `P_cont=π*H*R²*J0²/(2σ)`，H为半高。变σ时用独立高精度积分并收紧积分点检查，不用P_sample冒充P_cont。

无自感与自感都输出 `Q_o/c_num`、`Q_o/c_sample_ref`、`Q_o/c_cont_ref`。固定区域r>.75R与r<.25R，统一z范围和权重。主判据看同采样Q场误差，不再追着39对25改求解器。

## 8. E：必要回归与法国离线映射

### 8.1 自感只保留有意义的两档

修复后沿用Doc217：j_tol=.01/max40 与 .003/max80，relax=.3；若同一离散已缓存可共用几何。比较完整终态J/Q，径向均值作为辅助。

目标两档J差≤1%、功率差≤1%，真实φ残差达标。保留全皮参考精度和总功率约6%历史误差，不把Picard收敛解释成物理误差为0。

不得通过40%功率宽松门宣称严格精度验证通过。CPU求解耗时、装配、重构、host/device同步及峰值内存分开记录；不做性能优化大项目。

### 8.2 法国离线三快照

复用Doc217 dp=.015、10474粒子reload，校验hash和几何；不重新生成不同粒子云混淆回归。

1. 常σ=16、固定电流，修复D2与Doc217结果对照。
2. 既有诊断空间变σ，明确不是文献真实温度。
3. 若前两项通过，用明确来源的非均匀T快照得到σ(T)；σ(T)表不完整时不自行宣称文献一致，记录采用值和可调项。

打印R=.25、总高=.185、轴向、f、线圈参数及来源；线圈细节未确认的仍保留不确定标签。不要把历史代码值自动写成文献值。

诊断快照比较边界Jn、Q空间分布、原始功率，不能只因production与D2总功率差小就判所有场正确。

### 8.3 修正功率模式口径

Doc217 FM约1.64W、SC约2.14W且scaling状态不同。检查所有入口最终有效参数，保证打印和实际执行一致。

同温度/σ/线圈/电流、同自感/归一配置下，离线与热交接入口应得到一致Q。目标P差≤0.1%；若不同，列出唯一原因并固定，不保留模糊的“历史默认可能生效”说明。

增加一次50kW标定检查：输出原始P、比例s、实际电流变化、归一后P和Q。固定σ且电磁系统线性时I按sqrt(s)缩放、Q按s缩放；若代码只缩Q，明确是热源标定而非重新求解场。

必须区分“初始标定后固定I”和“每次EM更新保持恒功率”。两者是不同工况，显式选一，按既有法国复现配置说明；本轮不悄悄互换。甲方60kW另列，不混进50kW文献对照。

## 9. F：通过后开始法国固定热源 + Robin/辐射

### 9.1 目的与范围

本阶段检验真实热边界与Q传递，先关闭流体运动，固定粒子几何。它不是法国自然对流结果验收，不要求短时间形成双环流。

使用同一份通过验证的D2热源快照，初始温度沿用现有明确配置；绝对温度使用K，辐射采用K⁴。表面、侧面、底面参数逐项写出数值/来源/是否校准。不能把上轮的 `thermal_bc=none` 算例当作Robin已验证。

### 9.2 三个有限测试

| ID | 热源 | 边界 | 检查目的 |
|---|---|---|---|
| T0 | 固定Q，50kW诊断标定 | 绝热，无流动 | 完整热模块中的Q积分与升温 |
| T1 | Q=0 | 侧/底Robin、顶部对流+辐射 | 冷却符号、面积、量纲和失热分项 |
| T2 | 同一固定Q | 同T1 | 加热-散热共同能量收支与空间温度响应 |

不能用简单手工ΔT更新替代真实热扩散/边界算子。先单步核对，再推进到1s物理时间；若稳定性限制使代价异常，报告稳定dt与子步数，不为了跑完而改边界参数。

T2增加一次热时间步减半，其他不变。温度均值/热点及累计边界热量变化目标≤2%，只做这一次，不扩展成无限时间步扫描。

### 9.3 能量记录

逐步/累计记录：输入EM热量、侧面对流、底面对流、顶部对流、顶部辐射、流体内能变化、未解释残差。

统一向外失热为正：`ΔU = E_in - E_side - E_bottom - E_top_conv - E_top_rad + residual`。

使用实际热模块的ρ/cp/Vi约定；cp(T)变化时用对应内能/焓积分，不拿常cp公式冒充精确累计能量。残差以输入热量、失热绝对量与内能变化中的稳定尺度归一，避免加热和冷却接近抵消时分母趋零。

门限：单步能量闭合≤1%、累计≤2%；如边界面积/离散积分误差主导，先给独立面积自检，不能通过降低失热系数强行过门。

输出Python r-z温度剖面、Q剖面、顶部/侧面温度及失热曲线，使用固定范围和物理坐标。此阶段无流动，不画“环流已形成”的图。

### 9.4 如何衔接下一阶段

当T0–T2通过，提供接回已验证浮力路径的配置清单和一次短程运行预算。若时间允许，可做≤1s热-流入口烟雾测试，只验证变量更新/力方向/能量身份；不得据此宣布文献环流对齐。

法国长时、σ(T)持续反馈、μ(T)以及Fig7/8比较列为下一轮任务，不在本轮未完成核心修复时抢跑。此前Fig7/8数字化若可复用，记录文件身份与精度限制，不重新凭空生成参照流线。

## 10. 退出规则：避免再次漫无目的测试

- 壁值修复后M4/M7仍失败：检查同一输入下Python/C++装配行、F壁值、规范和符号；不先改邻域参数。
- 去减均值后出现不相容：检查原始物理行、体源与表面通量；不恢复经验减均值掩盖。
- 数值RHS失败、解析RHS通过：只定位F/σ重构和微分，不改热边界。
- C++通过但独立通量测量失败：先做测量自检；不能自动归因求解器。
- 只有法国入口不同：查profile、归一、数据同步与边界几何；不要重新打开方腔。
- A–D未过：停止F，只交付失败最小复现和唯一后续方向。

本轮最多一次必要分辨率细化、一次T2热时间步减半；禁止增加新的候选家族来避开已定位错误。

## 11. 交付结构与验收清单

根目录建议 `docs/ophelie/reference/doc218_runs/`。

| 子目录 | 必交内容 |
|---|---|
| `00_entry/` | source/bin hash、参数、命令、状态、依赖版本 |
| `A_repro/` | 两个缺陷的旧行为与最小输入 |
| `B_fix/` | 壁面值/相容性/失败策略的精简diff和解释 |
| `C_cpp/` | M1–M8完整CSV、Python/C++对照、相位测试 |
| `D_measure/` | 外边界/8子区真实积分点和通量、解析自检、功率口径 |
| `E_regression/` | 自感两档与法国离线快照 |
| `F_thermal/` | T0/T1/T2能量与温度剖面、减半热步比较 |

根报告：`ROUND218_CLOSURE.md`、`RESULTS_TABLE.csv`、`ROOT_CAUSE_AND_FIX.md`、`REMAINING_LIMITATIONS.md`、`NEXT_DECISION.md`。

完成前逐项检查：

- [ ] 边界F与gradφ在同一真实壁面位置评估。
- [ ] 数值壁值通过，而非只用制造解解析壁值。
- [ ] 正常路径取消边界RHS算术减均值。
- [ ] M7/M8验证非零体源及侧/端面通量。
- [ ] 不再存在拟合失败后φ=0、gradσ=0或gradφ=0的静默兜底。
- [ ] SparseLU求解状态检查完整，失败阻止PASS。
- [ ] 装配、重构和测量的模式/精度已明确。
- [ ] C++实际调用链通过M1–M8，不继承Python PASS。
- [ ] 独立通量覆盖数值RHS，不把∫divJ直接改名表面积分。
- [ ] P_cont不是P_sample的别名。
- [ ] 自感Q比有同粒子解析对照，历史差异不被删除。
- [ ] 法国不同入口的最终功率配置一致且可追溯。
- [ ] Robin/辐射实际运行，所有能量分项有符号与单位。
- [ ] CPU主机求解被真实标注，未称全GPU完成。
- [ ] 方腔状态保留，独立A-phi分支没有被改动。

最终结论分开写：实现修复 / 场精度 / 独立守恒 / 自感迭代 / 法国热边界。不要用一个 `passed=1` 合并所有层次。


<!-- ===== END 218_OPHELIE_D2_BOUNDARY_COMPATIBILITY_AND_FRENCH_THERMAL_CURSOR_PLAN_2026-11-10.md ===== -->

---

<!-- ===== BEGIN 219_OPHELIE_VALIDATION_CLOSURE_AND_FRENCH_CONVECTION_RESTART_CURSOR_PLAN_2026-11-10.md ===== -->

<a id="219-ophelie-validation-closure-and-french-convection-restart-cursor-plan-2026-11-10-md"></a>

# 219_OPHELIE_VALIDATION_CLOSURE_AND_FRENCH_CONVECTION_RESTART_CURSOR_PLAN_2026-11-10.md

# Doc219｜验证收口与法国自然对流重启：详细执行任务书

日期：2026-11-10  
执行对象：Cursor / SPHinXsys OPHELIE 开发工作区  
基线：DOC218_ROUND_PACK_2026-11-10  
主方向：有限补齐验收缺口 → 固定电磁空间热源 → 接回已验证浮力路径 → 联合检查温度和环流。

## 0. 必须遵守的范围与结果定义

Doc218 已在实际 C++ D2 路径通过 M1–M8，壁面驱动场取值和错误边界减均值已修复。本轮继续该方案，不开发新电磁离散，不重新打开方腔，不改独立 A-phi/Laplace 论文分支。

现在的重点不是继续寻找电磁路线，而是修正有限的参考计算、测量自检、能量积分、失败传播和入口激励，然后真正恢复法国自然对流推进。

### 0.1 阶段与进入条件

| 阶段 | 工作 | 进入下一阶段条件 |
|---|---|---|
| A | 复用CSV修正功率参考/完整测量链 | 参考身份正确，原结果重新计分 |
| B | 失败传播、严格不变性/相位、最终二进制回归 | 同一实际调用链稳定通过 |
| C | 统一入口激励，完成逐步能量积分 | 固定热源T0/T1/T2真实闭合 |
| D | 已验证流体算法映射到法国case | 明确没有悄悄换回旧流体配置 |
| E | 固定空间Q下短程自然对流 | 温度/流向/能量/测量可解释 |
| F | 下一阶段决策 | 只给一个主方向，不再列多个可选算法 |

A的离线复算与B/C的开发可独立安排，但所有修改、依赖运行按先后执行。不要求增加并行代理或额外执行环境。

### 0.2 本轮不做

- 不扫描CFL、核校正阈值或多套边界权重，不重新修方腔。
- 不升级D2为全仓库默认；法国验证入口使用显式D2 profile即可。
- 不要求本轮完成全设备Krylov；主机double SparseLU可以保留，需如实报告。
- 不大规模迁移TEAM7，不开发搅拌STL工况。
- 不直接启动数百/数千秒计算，不在本轮扫描μ(T)、h、σ来凑双环流。
- 不把短时出现一个回流，称为文献复现或双涡结构已验证。

### 0.3 允许的收口结论

1. 圆柱测试范围内D2电势验证已阶段性收口；尚非任意几何、强趋肤、TEAM7均已验证。
2. 法国固定热源与热边界的能量传递已真实验证；不能用一个 `passed=1`覆盖未测量量。
3. 法国自然对流已恢复可诊断推进；是否对齐文献，必须同时看温度和流动。

## 1. Doc218 已知数据：保留成果，修正过度结论

### 1.1 可继承的有效结果

| 项 | Doc218结果 | 本轮处理 |
|---|---:|---|
| M1常量场抵消 | eJ≈0.37% | 保留 |
| M2旋转场 | eJ≈0.009% | 保留 |
| M4二次电势 | eJ≈0.18% | 保留 |
| M5非线性电势+电流 | eJ≈0.22% | 保留 |
| M6变σ | eJ≈2.35%、Q场误差≈3.05% | 保留，不误读为参考功率脚本正确 |
| M7/M8非零体源 | eJ≈0.22%/0.62% | 保留 |
| 自感两档 | 径向J差≈0.054%、功率差≈0.089% | 迭代误差受控，物理精度另列 |
| 50kW热入口主机double | φ代数残差约4e-11 | 保留，但不是完整物理行/重构误差 |

### 1.2 独立复核已发现的验收问题

- M3−M2严格门未在C++实现；直接CSV复算电流差约0.544%、Q差约0.768%，实际满足门限。
- D脚本混用F的RMS与环向电流壁面幅值，导致部分P_cont/P_sample错误。
- IDW测量自检绕过插值；解析粒子电流经同一IDW后，M2侧壁法向L2约3.252%，基本解释原测量的3.25%。
- 热能量用末态失热率×1s估计，未逐步累计；带符号的 `energy_cap_rel_err <= positive_tol` 不能验证Robin能量闭合。
- ∇σ/重构失败静默留零、SparseLU失败传播还不完整。
- 法国两个入口功率比近似(8/7)²，强烈指向初始安匝数/每匝电流差；必须用源码和同激励重跑确认，不能只根据比值宣布根因完成。

## 2. 起点保护、身份与文件范围

历史EM主树为 `feature/EM-Biot-Savart`，HEAD `bc411c43718ca7ecaf4c1ed3e29aefeb4e6090c4`；Doc214–218含未提交代码。流体验证工作树历史HEAD为 `67fef44ef8a2cd6f5f7712bf2e66038cc6de89d3`。两者都必须重新核查实际状态，不能只靠HEAD推断源码相同。

开始前保存当前tracked diff、必要untracked源文件清单、源hash、二进制hash和编译精度。禁止reset/clean丢弃工作，不自动commit/push。文档和历史结果保留；本轮写 `docs/ophelie/reference/doc219_runs/`。

优先定位并修改：

- `d_independent_flux_power.py`
- `electromagnetic_ophelie_doc217_d2_local_poly.h`
- `electromagnetic_ophelie_doc218_d2_manufactured.h`
- `electromagnetic_ophelie_phi_component.h`
- 法国 `natural_em` 与 `complex_joule_to_heat_one_way` 入口
- 实际提供热扩散、Robin/辐射更新的头文件，不能只在入口拼能量估计。
- 已有法国热-流入口；先审计其流体算法身份，再复用。

保留原函数/文件名称没有问题，新增日志标Doc219修复身份。不要复制第二套求解器来专供测试。

## 3. A1：修复功率参考，不重新运行电磁

### 3.1 输入与公式

先复用Doc218 `particles_M2/M5/M6.csv`。从制造解配置读取R、半高H、σ0、B0、ω、Fchar与J0；独立参考的幅值来自输入定义，不能从数值解反推，也不能用含多个物理分量的F_RMS重建J0。

统一峰值相量约定：`Q*=|J*|²/(2σ)`。

常σ、`J*=J0*(-y/R,x/R,0)` 时：

`P_cont = π*H*R²*J0²/(2σ)`，其中H是半高，总高为2H。该公式本来就适用于随r线性变化的环向电流；不能写成“只适用于恒定幅值方位电流”后解释一倍差异。

M2用 `J*=-σF_R` 的规定旋转场；M5/M6用制造解规定J0。M6连续功率用正的σ(x)做独立高精度几何积分，不把σ均值代入常σ公式。

### 3.2 必须回收的三个功率

| 名称 | 定义 |
|---|---|
| P_num | ΣVi Q_num |
| P_sample | 同粒子同Vi采样独立解析J*与σ后的ΣVi Q* |
| P_cont | 连续解析公式或独立高精度几何积分 |

报告P_num/P_sample的场误差、P_sample/P_cont的几何求积误差；两者不混合。非均匀σ的积分将点数加倍一次，目标P_cont变化≤0.1%。

### 3.3 用于复核的参考数值

Doc218 CSV独立复算约为：

| 项 | 正确复核值 | 原脚本问题 |
|---|---:|---|
| M2 P_cont | 66974W | 原约34023W，幅值混淆 |
| M2 P_sample | 69310W | 与P_num几乎一致 |
| M5 P_sample | 69310W | 不应以驱动F_RMS生成J0 |
| M6 P_sample | 48058W | 原约82411W，参考电流幅值错误 |
| M6 P_num | 48074W | 与正确P_sample差约0.033% |

这些是数据复算检查值，容许CSV有效位数造成小差异。不得通过硬编码以上数值让脚本通过；公式和配置改变后应自然更新。

### 3.4 A1验收

输出 `POWER_REFERENCE_CORRECTION.md`、修正脚本、逐case机器可读数据、原/新参考并列表。明确撤回旧“连续与采样相差一倍是预期”的解释，不撤销正确的数值场结果。

## 4. A2：完整测量链自检，避免把插值误差当漏电

保留真实几何表面积分和8个子区域。新增三条不同路径：

1. **解析直评**：在积分点直接计算J*，检查几何/法向/求积；
2. **解析粒子→IDW→积分**：把J*采样到实际粒子，经与数值J完全相同的kNN、权重和覆盖规则；
3. **数值粒子→IDW→积分**：实际待测结果。

必须同时报告三条结果、以及2与3的差；第2条是测量误差底限，不是可以任意相减后宣称物理误差为零的校正项。

M2第2条侧壁Jn/Jchar约3.252%，预期接近第3条3.25%。确认后说明该量主要受测量插值限制，不再为此改电磁方程。

覆盖判据以实际dp/h和最大最近点距离等物理尺度定义，不能只用 `distance <= 1.5*median(distance)` 自动声称有效。分别报告各面与每个子区域。

门限沿用：净通量≤1%，子区域闭合通量≤5%，完整测量覆盖≥95%；局部法向L2目标≤5%。若测量链自身未达精度，标 `MEASUREMENT_LIMITED`，同时保留解析场精度验证结论，不把它偷偷改成PASS或FAIL求解器。

只做一次现有积分点密度加倍对照；不扫描kNN/IDW权重以调到最好结果。输出全流程自检，而非仅截图。

## 5. B：完成失败传播与最终C++回归

### 5.1 统一自适应拟合

∇σ、装配、F壁值、φ重构与报告用同一成功拟合策略。缓存实际成功的ids/h/系数映射；不得装配扩邻域成功后，重构又返回不足秩h0并默认为0。

任何必要拟合失败都返回明确失败状态，含粒子位置、阶段、邻居数与条件数。禁止∇σ=0、gradφ=0的静默替代。结果变量不能被热模块继续消费。

### 5.2 SparseLU状态和原始物理残差

分解后与solve后分别检查状态、解是否有限；失败不得只打印WARN后继续。失败状态逐级传到 `exec...SolveReconAndPower`、Picard、法国入口，终止该次有效推进并返回非零状态。

保存pin前物理方程K_phys/b_phys；保留pin解及最终demean解。输出：

- pin系统求解残差；
- 最终φ的原始物理行残差，包含被pin替换的那一行；
- 内部与边界未缩放残差分项；
- 常数零空间K_phys·1缺陷；
- 解析/数值/实际输出精度及规范。

不再将pin前解的约4e-11代数残差宣传为单精度输出梯度或物理边界也达到4e-11。

### 5.3 double范围与失败测试

源码中SparseLU已是double，但部分局部矩阵/导数系数仍可能是Real。准确打印 `fit_precision / assembly_precision / solve_precision / output_precision`。本轮不默认全面升级；只有残差或重现性证据需要时才调整，并单独归因。

增加两项小型失败传播测试：人为让局部邻域不满足秩条件；受控触发求解失败状态。可以在测试层注入失败，不污染正式输入。验收是入口确实FAIL、不会产出可供后续使用的Q，而不是强求真实工程case自然失败。

### 5.4 补齐实际验收门

- 在C++正式实现M3−M2完整向量差≤1%、Q场差≤2%，不能只测M3对参考的5%误差。现CSV复算约0.544%/0.768%，应能通过。
- 一次非零复数整体相位旋转：J同相位旋转、Q不变，目标J差≤0.1%、Q差≤0.2%。明确实/虚符号和峰值/RMS约定。
- 最终二进制跑一次M1–M8，用同一输入粒子，沿用Doc218原阈值；包括数值体源、数值壁面值、数值∇σ。
- 只在新修复导致失败时定位；不新增算法候选、不盲目细化。

自感本轮不默认再跑两套长回归：若φ求解/重构输出改变，复跑j_tol=.003那一档；若仅报告/失败处理改变，保留Doc218两档结果并提供最终版本一次短回归。

## 6. C1：统一法国入口的有效电流

### 6.1 当前线索

`2.14263 / 1.64045 ≈ (8/7)²`，相应电流比约8/7；热入口日志的标定后电流174.584A与放大系数152.761之比也约8/7。

首要假设是一个入口保留默认8安匝、另一个显式7匝×1A。必须沿解析/CLI覆盖/参数同步/线圈构造逐步核对；它目前是强线索，不是无需核查的最终结论。

### 6.2 最终激励清单

统一输出：N_turn、I_peak_per_loop、I_rms_per_loop、NI、频率、相位、线圈中心列表、R、segments、固定电流或目标功率、是否自感、σ、几何与reload hash。

若同时给NI与I_per_loop且冲突，报错或执行明文优先级；不允许两次sync把用户设置静默覆盖。

### 6.3 有限验证

同一reload、σ16、无自感、同D2、每匝峰值1A、N=7，分别运行natural_em与热入口（热步为0或独立EM-only快照模式）。原始P差≤0.1%，粒子Q相对L2差≤0.2%。

若差异消失，归档根因与修复；若仍不同，比较A源和实际线圈列表，不能再仅标为“跨入口口径不同”跳过。

再做一次50kW初始标定：输出原始P、sqrt(Ptarget/Praw)、最终电流、重新求解P和Q。固定温度下应服从线性场/平方功率关系。

本轮法国主配置选择“初始标定50kW后固定I”；不每次更新自动恒功率，不混入甲方60kW。以后若开展恒功率对照，另列工况身份。

## 7. C2：真实累计热量与验收逻辑

### 7.1 在实际更新处累计，不在入口估计

每个热子步记录实际使用的能量：

- E_EM；
- E_side、E_bottom、E_top_conv、E_top_rad，向外失热为正；
- ΔU；
- 内部导热总和/非守恒残差；
- 其他明确施加的温度约束能量，如无则为0。

必须在施加热通量的实际时间层积累。若算子用T^n，累计对应P(T^n)*dt；若用隐式或分裂子步，按实际更新积分。不能用最终P_loss乘总时长，也不能为了闭合而从ΔU倒推失热。

常物性时 `U=Σm_i cp(T_i−T_ref)`；移动粒子时质量定义必须与流体一致，不把密度/体积重估带来的变化错计为热量。cp(T)未来开启时用对应u(T)积分，本轮先固定Table1物性。

### 7.2 正确的能量门

定义：

`R_E = ΔU − (E_EM − E_side − E_bottom − E_top_conv − E_top_rad + E_other)`。

采用非负归一尺度，例如：

`S_E = max(|ΔU|, |E_EM|+|E_side|+|E_bottom|+|E_top_conv|+|E_top_rad|+|E_other|, E_floor)`。

E_floor以固定物理能量尺度乘数值容差定义并打印；不能用过大常数稀释误差。

验收 `abs(R_E)/S_E`，不是带符号量单边比较。逐步≤1%，累计≤2%；同时保留绝对J值。

Q=0冷却、加热小于散热时仍按同一公式验收。失热正值和温度有限只是合理性检查，不代替能量闭合。

### 7.3 有限重跑T0/T1/T2

复用Doc218粒子与参数：T0绝热加热、T1无热源Robin/辐射冷却、T2加热+散热，先各1s。T2增加一次dt=.1→.05对照（稳定性若要求更小，记录真实限制）。

要求：能量真实闭合、各失热分项符号正确、dt减半温度与累计热量变化≤2%。温度/能量门必须进入最终passed，而非只放在MD人工估算。

旧 `energy_balance_rel_err`、`closure_mismatch_vol_frac` 如果仍按“仅Joule升温”定义，在Robin时标NOT_APPLICABLE；不要同时打印巨大旧错误和passed=1却不解释。

## 8. D：接回已验证流体路径，避免算法身份倒退

方腔通过不自动代表法国入口使用同一求解配置。建立 `HYDRO_PROFILE_MAPPING.md`，逐项核对：

- 压力/连续性与黎曼求解器；
- 黏性离散及隐式/显式路径；
- B矩阵与校正开关；
- TVF/transport velocity、粒子位置更新；
- 更新位置→cell/relation→热/力计算顺序；
- 浮力β、重力、Tref与静水压处理；
- 壁面不可穿透/无滑移、自由液面机械边界；
- 时间步限制、声速与温度耦合时序。

不要求把方腔几何和壁面物理照搬到法国；要求数值实现身份明确，法国特有边界有出处。旧法国入口若缺少验证后修复，只迁必要内容，禁止把整套旧实验开关叠加回来。

若两工作树难以直接共用，先构建一个隔离法国自然对流验证入口，复用已验证dynamics；列出最小迁移diff，不盲目合并整个流体分支。

## 9. E：法国固定空间热源下恢复自然对流

### 9.1 物理与运行身份

沿用已确认玻璃D=.5m、总高=.185m、dp=.015m。竖直轴明确打印；本阶段无搅拌。输出全部热边界参数、绝对温度K、重力、ρ/cp/k/μ/β、热源和功率模式。

先保持σ16、μ4Pa·s等当前Table1常物性配置，明确是阶段基线；不把它写成完整σ(T)/μ(T)文献模型。

本阶段仅由温差浮力驱动，明确关闭Lorentz体力、搅拌和Marangoni；除非已有文献与单独任务要求，否则不附加其他驱动力。侧/底机械边界与顶部自由表面/自由滑移处理分别打印，避免用方腔四壁无滑移设定代替法国顶部条件。温度随流体运动的输运应沿已验证的粒子热方程处理，不再额外重复加一遍网格对流项。

当前h_side=300、h_bottom=35、h_free=20、ε=.8作为Doc218延续配置，分别标来源与校准状态。ε不能声明是唯一文献规定值。初始T=1473K延续时，必须保留“实际初始温场尚未对齐文献”的限制。

### 9.2 固定Q必须是固定空间分布

从验证后的D2在固定几何/温度快照生成Q_ref(x)，标定总功率50kW后冻结该空间形状。

粒子移动后不能让初始Q值仅绑定粒子ID随粒子漂移，那会把电磁空间热源变成物质标签。每次热更新按当前位置从参考空间场取Q，或用等价已审计的空间映射。

输出映射误差、覆盖率、ΣViQ。映射先在原位置自检，Q相对误差目标≤1%、总功率差≤1%；不得出现未覆盖点静默Q=0。若使用每步离散归一维持50kW，明确这是固定空间形状/固定总功率诊断协议，报告归一系数；不声称等价于全部电磁随温度反馈。

本阶段不每步调用D2；先排除EM更新、σ反馈与流体之间的混淆。Q是W/m³，粒子功率QVi、升温率Q/(ρcp)，不要混用。

### 9.3 只做两个主工况

| ID | 设定 | 目的 |
|---|---|---|
| N0 | 同Q/边界，β=0，1s | 无热浮力基线，量化数值背景速度 |
| N1 | 同Q/边界，β为既有文献基线 | 温度驱动流动，逐级推进 |

N1从同一初态先1s，再续算至5s；若状态正常、成本允许，续至20s。是同一个case的检查点，不从头重复三次。上限20s，不自动转长时任务。

打印实际walltime/物理秒和剩余预算；若20s代价明显超出本轮可用资源，保留5s可续算checkpoint，并明确还不能据此判断稳态/双涡。不要为了完成20s改变物理黏度、温差或声速约束。

N0不是要求速度严格等于0，判定基于其相对浮力工况的量级、压强平衡及先前静止基线。若N1信号与N0同量级，先解释驱动力与时间尺度，不直接改算法。

### 9.4 必须同步输出温度与流动

每个检查点输出：

- Tmin/Tmax/体积平均T、冷热区位置、侧/底/顶温度；
- Q(r,z)与T(r,z)，统一固定色标和坐标；
- 原始三维速度、r-z环向平均速度、非轴对称波动强度；
- 速度RMS/最大值、径向/竖向分量；
- 边界法向速度诊断与覆盖率；
- 累计热量、内能与机械动能，明确这里只是热收支还是含机械能全收支。

以实际ν=μ/ρ、α=k/(ρcp)、H和实际温差计算Pr/Ra与黏性/热扩散时间。20s不等于热稳态；不能只因涡结构未形成就认定浮力失败。

本配置1473K下历史失热约140kW、加热50kW，平均温度下降是可能的正常结果。局部温差仍可驱动环流。不要把整体冷却误判成热源没加载，也不要通过调低h强行改成全域升温。

## 10. Python剖面、矢量与环流诊断规范

必须交付Python生成的温度剖面+速度矢量，ParaView仅作辅助。

1. 从原始粒子字段出发，输出插值网格、有效域mask和覆盖率；禁止在无数据区补成大片零速度后画闭合流线。
2. 用固定等长箭头显示方向，再另图用速度幅值着色；避免弱回流被按幅值缩放隐藏。
3. 输出至少一个真实子午面切片与环向平均r-z图，两者标注不同；不能把强非轴对称流动平均后当作真实二维流线。
4. 若计算轴对称流函数，明确采用 `u_r=−(1/r)∂ψ/∂z`、`u_z=(1/r)∂ψ/∂r` 的约定及拟合/积分残差；轴线特殊处理。
5. 对三维SPH场若连续性近似不够，流函数重构应标测量限制；矢量图仍保留，不凭一张ψ等值线判定物理解对错。
6. 输出主涡中心、旋向、近壁上升/下降带；外围第二涡只有出现独立闭合结构、方向一致且非插值伪影时才标识。

法国Fig7主/次环流不预设为上下两层。先读取已有图像/数字化参照，沿其真实结构比较。缺少作者速度数值时，不杜撰绝对速度误差。

Fig8若采用T/Tmax，使用绝对K与相同参考定义，不能用每张图min-max归一掩盖温差不足。图形数字化误差单列，先比温度热点/梯度，再谈流线是否接近。

## 11. 是否加入σ(T)：本轮有明确停止点

只有A–C通过且N1已可解释，才准备反馈接口清单，不默认把σ(T)/μ(T)与流体同时打开。

本轮优先产出固定Q自然对流结果；若其温度/流动异常，先在这一简化问题中定位。下一轮再按相同几何和热边界增加一次温度相关EM更新，比较Q与环流改变。

若当前资料不够确定σ(T)，保留明确可校准模型及来源；不能将博士论文或期刊中缺失的参数写成已确认。用户允许校准不等于可以省略数值验证或隐藏参数改动。

## 12. 验收、停止规则与测试预算

| 情况 | 唯一优先动作 |
|---|---|
| 修正参考后原CSV仍符合场门 | 保留成果，不重跑物理 |
| 测量自检解释3.25%法向误差 | 标清测量底限，不改方程 |
| 失败处理补丁改变场值 | 检查缓存/重构差异，最终M1–M8一次回归 |
| 两入口同电流仍不同 | 比较A源、线圈列表、配置，不用归一掩盖 |
| T0/T1/T2能量未闭合 | 查实际通量更新、单位/质量、边界能量；停止N1 |
| N1有流动但温场离文献远 | 优先温度/边界/初态与时间尺度，不继续泛化怀疑浮力 |
| 流线与矢量冲突 | 核对后处理与连续性测量，不凭图改算法 |

预算：离线参考复算；最终M1–M8一次；一次复数相位；两个小型失败传播；同激励两入口；T0/T1/T2各1s及一次T2半步；N0 1s、N1检查点1/5/最多20s。除新失败有明确证据，不增加测试矩阵。

## 13. 交付清单

根目录建议 `doc219_runs/`，保留可复放脚本和配置，输出：

| 目录 | 内容 |
|---|---|
| `00_entry/` | 当前代码/二进制身份、精度、命令、参数摘要 |
| `A_reference/` | 正确功率公式/配置、CSV复算、完整IDW自检 |
| `B_closure/` | 失败传播、原始物理残差、M1–M8、M3不变性、相位 |
| `C_excitation/` | 安匝数追踪、同I入口对照、50kW标定 |
| `C_energy/` | 每热步与累计收支、T0/T1/T2、严格门结果 |
| `D_mapping/` | 已验证流体→法国入口逐项映射与最小diff |
| `E_convection/` | N0/N1配置、checkpoint、原始数据、Python图和诊断 |

根报告必须有：

- `ROUND219_CLOSURE.md`：参考正确性/实现可靠性/电磁精度/热收支/自然对流分别结论。
- `RESULTS_TABLE.csv`：真实数值+阈值+PASS/FAIL/NOT_RUN/MEASUREMENT_LIMITED，不能用“≤5%”代替测得数值。
- `ROOT_CAUSE_AND_FIX.md`：区分电磁物理解错误、后处理参考错误、测量误差和验收逻辑错误。
- `REMAINING_LIMITATIONS.md`：几何范围、CPU主机求解、强趋肤/TEAM7、反馈/稳态等限制。
- `NEXT_DECISION.md`：根据N1温度与环流结果只给一个下一主方向。

最终检查：

- [ ] M2/M5/M6参考幅值来自制造解配置，而非F_RMS。
- [ ] IDW自检真正经过同一插值链。
- [ ] 失败传播阻止无效Q进入热模块。
- [ ] 原物理行和pin系统残差分开。
- [ ] M3严格门和相位测试已实际执行。
- [ ] 同I下法国两个入口的A/Q可比较。
- [ ] 热边界分项逐步累计，绝对能量残差参与passed。
- [ ] 旧仅Joule升温诊断在Robin模式不再误导。
- [ ] 流体profile可追溯到已验证配置。
- [ ] 固定Q是固定空间分布，粒子移动后正确重采样。
- [ ] Python温度/矢量图、覆盖率和原始字段一并交付。
- [ ] 未形成第二涡时先结合温场与时间尺度解释，未盲目调参数。
- [ ] 方腔和A-phi独立路线没有被重新改动。

本轮结束应当回答一个具体问题：在已核查的空间热源与热边界下，我们的法国几何里出现了怎样的温度分布和自然对流；若尚未接近文献，证据主要指向温度条件、演化时间还是流体实现。


<!-- ===== END 219_OPHELIE_VALIDATION_CLOSURE_AND_FRENCH_CONVECTION_RESTART_CURSOR_PLAN_2026-11-10.md ===== -->

---

<!-- ===== BEGIN 220_OPHELIE_SORTING_THERMAL_IDENTITY_AND_HYDROSTATIC_REPAIR_CURSOR_PLAN_2026-11-12.md ===== -->

<a id="220-ophelie-sorting-thermal-identity-and-hydrostatic-repair-cursor-plan-2026-11-12-md"></a>

# 220_OPHELIE_SORTING_THERMAL_IDENTITY_AND_HYDROSTATIC_REPAIR_CURSOR_PLAN_2026-11-12.md

# Doc220｜温度排序完整性、能量初态与法国静水平衡修复任务书

日期：2026-11-12  
执行对象：Cursor / SPHinXsys OPHELIE 工作区  
依据：DOC219_ROUND_PACK_2026-11-12 的日志、分区数据和入口源码  
主线：排序与温场身份 → 能量初始化 → 静水平衡定位及修复 → 短程浮力验证。

## 0. 本轮目标与边界

本轮不是继续增加自由液面“静音”开关。先解决一个高度可疑的状态/测量错误，再确定静水压与壁面离散的具体缺陷。

必须回答：

1. 第50个平流步排序是否把温度或测量索引打乱？
2. 预设热核的初始能量是否正确计入账本？
3. 法国几何中，无浮力仍有0.02m/s流动，主要由哪个实际算子/状态更新产生？
4. 在温场可信、背景流动受控后，冻结温场与释放温场的浮力响应是否一致可解释？

保留既有电磁D2、50kW热源传递和方腔成果；本轮不重开电磁算法、TEAM7、方腔、搅拌和独立A-phi/Laplace分支。

### 0.1 执行顺序及门

| 阶段 | 工作 | 通过后才能做 |
|---|---|---|
| A | 当前源码/字段/排序身份审计 | 排序最小实验 |
| B | 零时间排序审计、首次排序前后对照、修复 | 可信非均匀温场演化 |
| C | 能量初态/约束能量/验收身份修复 | 热核能量验收 |
| D | 法国静水分阶段力平衡测量 | 选择唯一修复方向 |
| E | 有证据的静水平衡修复及短静水保持 | 浮力信号检查 |
| F | 冻结温场与释放温场的成对短程计算 | 决定是否续算法国自然对流 |

B/C失败时不再跑5s/20s热核来讨论环流。静水未过时不通过增大β或调低μ来掩盖背景误差。

### 0.2 不扩展范围

- 不再次扫描密度上限、力抵消、背景压、多种TVF组合。
- 不用逐步清零速度/力的方式制造“静水通过”。
- 不以永久关闭排序作为最终修复；不新增ghost粒子布置。
- 不推倒重写整个SPHinXsys流体框架；只隔离修复本项目必要内容。
- 不默认启用δ-SPH、变μ/变σ、Marangoni或Lorentz体力。
- 不默认启动长时计算；主结果先限制在0～2s，符合门后最多续到5s。

## 1. Doc219证据：哪些是事实，哪些是假设

### 1.1 可以保留的成果

- 法国两个EM入口同为7匝、每匝峰值1A时，原始P≈1.64045W，入口电流差异已纠正。
- 初始标定50kW后固定电流的激励路径已经接通。
- D2 M1–M8、常量附加场不变性、相位与失败传播有验证记录。
- 固定几何T0/T1/T2累计热收支已经通过。
- 这些成果不等于法国移动流体的字段排序、静水平衡、完整环流已通过。

### 1.2 最高优先级异常

G1释放热核，β=1e-5：

| t(s) | 平流步N | T中心(K) | T外区(K) | T均值(K) |
|---:|---:|---:|---:|---:|
| 0.138649 | 40 | 1639.81 | 1497.00 | 1540.19 |
| 0.208561 | 60 | 1528.94 | 1537.30 | 1540.11 |

β=0的G0在同一区间也出现约111K中心降温和40K外区升温。最高温度和均值几乎不变。

入口每50个平流步调用 `fluid_particle_sorting.exec()`。最大温度空间跳变恰好跨过首次排序。这是**排序或测量索引错位的强证据**，但当前包缺少最新字段注册/排序实现头文件，尚不能宣布具体数组错误已经被证实。

该0.07s区间热扩散长度约0.00030m，不能在没有进一步证据时把中心约111K突变解释为正常导热/外围加热。不要继续引用“5s外围Q把热核反转”作为已确认结论。

### 1.3 第二异常：初始能量基准顺序

当前入口先 `ledger.setH0(prev_energy)`，后施加预设温场。G0/G1末态残差约7.49MJ，按当前粒子质量×cp折算约67K平均温升，与施加热核的初始额外能量相符。

这支持“初始能量漏记”的解释，不证明热扩散实际丢失7.49MJ。必须重建账本后再判守恒。

### 1.4 静水异常不能只称液面噪声

Doc219关闭浮力后的分区数据：

| 区域 | 速度RMS(m/s) | 诊断净加速度z RMS(m/s²) |
|---|---:|---:|
| 全域 | 约0.0192 | — |
| 侧壁带 | 约0.0367 | 约3.40 |
| 液面带 | 约0.0141 | 约7.62 |
| 深部 | 约0.0126 | 约1.88 |

侧壁带速度最强，深部也有显著流动。诊断加速度来自额外力评估，不等同每步持续物理加速度；但无浮力持续速度说明不是纯后处理伪影。

β=1e-5、ΔT=20～200K的浮力加速度仅0.001962～0.01962m/s²。应检查压力—重力—壁面平衡，不能只围绕液面密度加新开关。

## 2. 起点保护与完整源快照

历史主树 `feature/EM-Biot-Savart`，HEAD `bc411c43718ca7ecaf4c1ed3e29aefeb4e6090c4`，含未提交Doc214–219改动。以当前实际工作树为准，先保存必要diff/清单/hash；禁止reset/clean丢工作，不自动commit/push。

Doc219主入口：

`tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/test_3d_ophelie_french_natural_convection_frozen_q.cpp`

历史入口SHA256：`2c3c53cbd7d823ae383fd34b9d7a54788665de757b2b9b8b94b724019dab6237`。不能假定执行时仍相同。

本轮必须额外归档实际使用的依赖源，解决上一包只有入口CPP的问题：

- 温度/ΔT注册、热更新和物性注册函数；
- 粒子排序、evolving/sortable字段机制、稳定ID机制；
- Boussinesq与GravityForceCK，以及历史力更新机制；
- `AcousticStep1stHalfWithWallRiemannCorrectionCK`、2ndHalf及实际依赖；
- 密度求和、FreeSurface密度重整、EOS、B矩阵；
- 温度空间统计、子午面统计、热量账本；
- 本轮直接修改的测试、头文件、Python脚本及CMake。

不要用Doc211等旧快照证明当前注册实现一定相同。新计划和结果放 `docs/ophelie/reference/doc220_runs/`。

## 3. A：字段身份审计

### 3.1 逐字段清单

建立 `FIELD_OWNERSHIP_AND_SORTING.csv`，至少包含：

| 字段组 | 需要确认 |
|---|---|
| Position/Velocity/Mass/Density/VolumetricMeasure | 框架默认排序、变量注册时机、实际数组身份 |
| Temperature/ThermalDeltaT | 是否共同随粒子排列；更新后的恒等关系 |
| 粒子cp/k/μ、相关历史量 | 哪些是物质状态，哪些每步按T重算 |
| ForcePrior及Gravity/Boussinesq/Viscous的当前与Previous字段 | 是否排序，增量替换是否只发生一次 |
| B矩阵/梯度/热Laplace临时量 | 是否在排序后无条件重算再使用 |
| JouleHeat | 当前是空间重采样，不应把旧粒子Q误作物质状态 |
| 统计缓冲、host指针、邻域索引、区域标签 | 排序后重新获取/重算规则 |

分类不能全部简单标为evolving：物质状态应随粒子排列；几何统计与临时系数可以重算，但必须在首次使用前重算；空间Q按当前位置取值。严格沿当前SPHinXsys API注册，不凭旧API名称照抄。

### 3.2 稳定粒子ID

复用框架中确实随粒子移动且不随排序重新编号的ID；确认它不是“当前排序索引”。如果没有，新增仅用于审计的整数稳定ID，并按正式机制随粒子排序。

验证ID唯一、数量不变。不得用数组下标、易碰撞坐标hash或排序后重新生成的ID做对照。

## 4. B：排序最小实验、分支判断和修复

### 4.1 S0：零物理时间、单独排序

同一法国粒子云，生成非均匀温场与可识别物质标记；确保排序确实导致非平凡排列。暂停全部力学/热推进和位置更新，执行一次排序，按稳定ID匹配前后记录。

检查：

- 每个ID的Position、T、ΔT、物性及持久历史状态没有改变，只排列改变；
- `T = Tref + ThermalDeltaT`（若该字段定义确为此）排序前后都成立；
- 总质量、总热量、Tmin/Tmax/Tmean不变；
- 从位置重新分区得到的T中心/T外区不变；
- 清理/重建邻域后，读回实际设备字段再次核对，不只比host缓存。

纯重排不含算术更新，持久标量/向量期望逐ID逐位相同；若框架确有格式转换，说明原因并给出ULP级容差，不能接受K级变化。

输出排序前/后完整CSV与 `S0_COMPARE.json`。仅总和相同远远不够——错误排列也能守恒。

### 4.2 S1：首次排序前后原始运行

复用G0/G1热核输入，先只跑β=0到0.3s，排序频率保持50。密集记录第48/49/50/51/52步：排序前、排序后、邻域更新后、Q后、扩散后、Robin后。

输出各阶段稳定ID对应的状态，T中心/外区、能量、字段恒等式以及热/力调用时间。统计函数重新获取字段，不沿用旧host指针。

### 4.3 S2：仅关闭排序的短对照

同一初态、相同物理时间步策略，排序临时关闭，跑到0.3s。此项用于定位，不作为生产修复。

按分支处理：

| 观察 | 优先修复 |
|---|---|
| S0逐ID温度错位 | 注册/排序字段列表，修持久热状态 |
| S0温度正确，空间统计跳变 | 统计索引、区域映射、host缓存/指针 |
| 温度正确，ΔT错位或反之 | 统一两字段排序及热更新权威状态 |
| S0正确，推进后某一算子突变 | 缓存绑定、排序后临时量/邻域重建与调用顺序 |
| 无排序仍出现同样突变 | 不再归因排序；按Q/扩散/边界分阶段快照定位首次变坏处 |

### 4.4 修复原则

只改证据指向的注册/缓存/重建，不用每步按初始坐标重设温度掩盖错误。禁止保留一个“影子温度数组”只让测量好看。

如果多个更新算子分别写T与ΔT，明确权威状态和同步点，避免一个算子用旧ΔT覆盖新的T。排序后Boussinesq当前/Previous字段也必须正确，防止ForcePrior的增量替换失效。

### 4.5 B阶段验收

- S0持久字段逐ID正确；
- 修复后排序开/关到0.3s，温度剖面不再出现111K级突变；无排序本身不授PASS；
- 按稳定ID比较两次短程的T与速度，允许浮点/演化差异，但归一温度RMS差目标≤0.1%、全域最大差≤1%初始温差；超限定位，不扫参数；
- 修复版本开启正常排序，通过至少3次排序的短回归。将“全局热量守恒但空间错位”的反例保留为回归用例。

## 5. C：重建热初态能量与约束记账

### 5.1 初始化顺序

正式开始时的顺序应是：完成几何/静水准备 → 最终设置温场及一致ΔT → 初始化对应物性/力 → 计算当前真实内能H0 → 清空本次物理演化累计量 → 记录t=0快照 → 推进。

同步重置 `prev_energy`、`prev_energy_time`、采样窗口和相关账本时间，避免只修改H0但首步dE/dt仍从旧温场计算。

从checkpoint续算时应继承原账本，或者明确开启新统计区间并记录偏移；不能每次输出/重启就重置H0消除残差。

### 5.2 两类温场模式要分开

- **释放温场**：初始热核属于H0；此后只计实际Q、边界和其他已声明热量。
- **冻结温场**：每次重新施加温场相当于外部温度约束。记录施加前后内能差 `E_prescribe`；它不是零热量演化，也不能按普通自然对流稳态验证。

冻结模式应在实际浮力首次使用前保证温场正确，包括排序后。冻结模式仅用于力学诊断，不能用它掩盖释放温场的状态错误。

### 5.3 能量验收

`R = H(t)-H0-[E_Q-E_boundary+E_prescribe+E_other]`。

使用固定质量和一致cp/u(T)，给出绝对J值与非负尺度归一残差。初始t=0残差应为数值零；释放热核0.3s累计相对残差目标≤1%，随后2s≤2%。严禁以最终温度反推一个“修正H0”。

复核原7.49MJ是否等于漏记初始热核能量，给出差值。如果重建H0后仍残留明显误差，再查移动粒子、导热和边界，不把初始化解释当万能豁免。

`smoke_passed`、`thermal_energy_ok`、`signal_resolved`、`literature_comparison_ready`分列。旧 `passed=1`可能仅意味着能跑完，不得在MD写成法国环流已验证。

## 6. D：最小静水问题与算子分阶段测量

### 6.1 基线输入

保留法国D=.5m、总高=.185m、dp=.015m、侧/底壁面、实际自由液面。Q=0、β=0、温度常值、无热扩散/Robin、μ=4、ρ0=2750、c0沿当前值；TVF先关，所有cap/cancel/经验归一实验关闭。

不要生成新粒子云；保存原reload hash、壁面生成参数/法向及两者几何距离分布。任何wall gap异常先纠正，不用额外力去压住。

初始化ρ/p必须与实际EOS匹配。当前线性静水密度是弱可压近似，应记录其与实际EOS/静水方程的误差，不将任意温压初值称作精确平衡。

### 6.2 H0：固定位置的一次算子诊断

在克隆/可完整恢复的状态中测量，不能连续调用会更新速度/历史力的算子后把结果相减当作同一时刻。

依次记录以下同源状态的结果：

1. 刚施加解析静水ρ/p；
2. 密度求和前/后；
3. FreeSurface重整前/后；
4. B矩阵更新前/后；
5. inner压力项；
6. wall压力项；
7. GravityForce与其Previous/ForcePrior；
8. 真正声学子步使用的总力和一次速度增量。

分开输出压力梯度、黎曼耗散、黏性与外力；零速时不应有虚构的黏性驱动。检查wall Riemann的重力/参考加速度约定、法向方向、inner/contact是否被重复累计。

输出按稳定ID的CSV，至少包括位置、深度、ρ、ρ_sum、p_EOS、p_bg、p_total、p_h参考、B条件数、邻域数、各力、质量和实际Δu。

### 6.3 分区与量纲

分别报告侧壁/底壁/角点/液面/深部；定义互斥区与重叠区，避免重复统计。输出径向、竖向、有符号均值和RMS，不只全域RMS。

重点检查：密度重整是否覆盖了正确静水密度；邻域不完整是否把侧壁粒子当成自由液面；wall外推是否包含正确重力；B是否改善线性场却破坏常压力零力；接触项是否与inner校正一致。

同时报告a/g与a/(gβΔT_design)。ΔT_design使用明确诊断温差（例如200K），不能用异常Tmax−Tmin放大信号尺度让门限变宽。

### 6.4 H1：正常静水保持

在H0定位前不做长保持。先0.3s，再在修复后续到1s；静水保持期间不得清零速度、力或钉住流体粒子。原settle若结束后清零速度，必须另报清零前状态和清零后的持续驱动，不能以清零作为稳态验收。

用真实每步Δu或一致时刻力统计验证额外诊断快照，避免把一个额外acoustic_1st的状态当作生产平均加速度。

## 7. E：按证据选择一项静水修复，不增加开关族

### 7.1 优先顺序

| 已定位缺陷 | 对应修复 |
|---|---|
| 当前/Previous力或状态排序、ForcePrior更新错 | 修字段和增量替换顺序，再跑H0/H1 |
| 密度重整打掉静水分布 | 核查原生自由面密度策略与初始化/连续性衔接；修单一一致策略，不全域clip到ρ0 |
| wall静水外推/重力参考符号错 | 修wall状态/Riemann输入，保留总压物理 |
| inner/contact校正不一致 | 对齐同一算子家族，查常压力零力与线性p梯度 |
| 自由面识别混入侧壁 | 用实际邻域与几何证据修分类，不取消整带物理力 |

先做实施修复，不先引入新SPH方法。所有改动在OPHELIE隔离路径中，不改变其他case默认。

### 7.2 原生路径作为一个基准对照

如果当前带B/专用wall-correction路径不能保持静水，允许在**同一法国几何**中跑一次当前主库成熟的原生自由面+wall Riemann路径，去掉本项目额外修正；这是一项有目的的基准，不是组合扫描。

先读实际源码确认原生实现，不因类名含Riemann就认为一样。比较H0各项残差和H1速度。若原生更好，迁移其一致实现；若两者都差，优先查共用初始化、壁面、密度与状态生命周期。

### 7.3 扰动压力不是禁区，但不直接复用闭盖结果

若前述实施问题排除后，证据指向总压静水抵消精度，允许提出一套保持边界一致的扰动压力/静水平衡方案，先写方程和离散差异再实现。

`p_total=p_h+p'`本身不等于闭盖。若自由面总压为0，必须满足 `p'=-p_h` 在真实自由面上；wall Riemann与密度/EOS也必须一致。不能仅删GravityForce而继续沿用原总压wall输入。

历史flat-lid结果可作算法诊断，不能当开口复现。不得无依据禁止所有p'路线，也不得把此前闭盖脚本直接复制成当前答案。

本轮最多一个证据支持的修复主方案加上述一次原生对照；若仍失败，交付精确最小反例，不再追加cap/cancel/速度清零配方。

### 7.4 E阶段门：按实际信号量级验收

记录修复前后同一静水测试的全域/侧壁/深部速度、持续偏置、真实Δu及压力误差。不能用一条任意固定“U<5e-3”就宣布所有弱对流可分辨。

暂定工程目标：无驱动1s速度RMS降至≤1e-3m/s，并较原0.02m/s显著下降；这只是进入浮力诊断的门，不是文献验证门。若实际目标环流更弱，后续信号门需更严格。

静水没有下降却仅增大β得到明显流动，不授E通过。瞬时力RMS未完全消失不一定否决，但必须证明没有持续偏置、没有人工耗散静音，且对目标信号的污染可测。

## 8. F：重做可信温场下的浮力对照

### 8.1 冻结温场：只检验水动力响应

同一静水准备状态，规定同一径向热核ΔT=200K，两组β=0/1e-5；Q、扩散和Robin关闭，冻结约束明确记录。

先0.3s，再到2s。输出轴心/外围有符号速度和空间向量差，检查力方向、增长过程、压力响应。冻结是外部维持的诊断温场，不宣称文献物理解。

### 8.2 释放温场：只改变是否冻结

两组相同初始热核，β=0/1e-5，正常排序、50kW空间Q、扩散和Robin/辐射恢复；H0在预设温场后建立。

先0.3s核对首次排序，再2s；温场/背景流动/能量都可信时最多续到5s。不得沿用Doc219错位温场的checkpoint续算。

空间Q按当前粒子位置重采样，不能随粒子ID漂移；保证所有流体状态和热量账本采用实际同一时间间隔。

### 8.3 比较向量场，不做“两个RMS相减即无信号”

把两组速度映射到共同固定空间网格，先比较有符号矢量：

`δu(x,t)=u_beta_on(x,t)-u_beta_off(x,t)`。

输出δu的空间L2、方向、轴心/外围剖面、动能与涡量诊断。两个速度RMS相近不能证明向量场相同。

如果信号仍弱，先评估相同基线重复或空间非轴对称波动给出的不确定度。只有信号大于测量/背景不确定度才能授 `BUOYANCY_SIGNAL_RESOLVED`；“force_ok=1”只说明施力公式接通。

## 9. Python图、数据与时间尺度

必须输出排序前后T(r,z)、热核短程T图、β0/βon速度矢量及δu。使用统一坐标、同一插值方法/色标，保留有效域mask和覆盖率。

至少保留一个真实子午切片与环向平均图；非轴对称扰动单列。流函数只是辅助重构，写清轴对称质量守恒约定与重构误差，不能用插值回流区强行标第二涡。

α=4/(2750*1150)≈1.265e-6m²/s。按实际局部尺度估计扩散、速度和温度变化上界，核对突变是否物理合理。总高度.185m的热扩散时间约2.7e4s，黏性时间约23.5s；2/5s只用于修复与响应验证，不是热稳态或文献涡结构终验。

自然对流期刊参照只在温度场与算法身份可信后恢复。Fig7结构不预设上下两层；Fig8温标不能用每帧min-max归一掩盖差异。

## 10. 结果状态、测试预算与停止规则

| 状态 | 含义 |
|---|---|
| SORT_IDENTITY_PASS | 零时间重排与推进中字段身份均正确 |
| THERMAL_INITIALIZATION_PASS | 初始热核能量入账，约束能量分列 |
| HYDROSTATIC_DIAGNOSED | 找到具体更新/算子导致的静水偏差，有同状态证据 |
| HYDROSTATIC_IMPROVED | 无人工静音，保持测试达到约定工程门 |
| BUOYANCY_FORCE_CONNECTED | 施力公式匹配，不能替代流场验证 |
| BUOYANCY_SIGNAL_RESOLVED | 向量场差分方向/量级可分辨 |
| FRENCH_LITERATURE_VALIDATED | 本轮不默认授予 |

预算：S0一次；首次排序前后0.3s与关闭排序0.3s各一次；修复后至少3次排序短回归；能量同程核对；H0单状态分解；H1原/修复各0.3～1s；必要时一次原生对照；冻结与释放温场各β成对0.3～2s，符合门后释放温场最多5s。

若S0证明排序没有问题，必须转向实际首次变坏算子，不继续为了支持假设修改注册。若温场未可信，暂停静水以外的对流判断。若静水仍有原量级流动，停止继续跑热核/调物性凑结果。

## 11. 交付要求

建议根目录 `doc220_runs/`：

| 目录 | 内容 |
|---|---|
| `00_entry/` | 当前源码/依赖/二进制hash、命令、模式、精度、状态表 |
| `A_fields/` | 字段归属/排序清单、稳定ID说明 |
| `B_sort/` | S0及首次排序多阶段CSV、开/关排序对照、修复diff |
| `C_energy/` | 正确H0、约束能量、初始7.49MJ解释和重跑收支 |
| `D_hydro/` | 同状态力分解、密度前后、分区压力与加速度 |
| `E_fix/` | 唯一修复与必要原生对照、静水保持数据 |
| `F_buoyancy/` | 冻结/释放成对结果、共同网格δu与Python图 |

根报告：`ROUND220_CLOSURE.md`、`RESULTS_TABLE.csv`、`ROOT_CAUSE_AND_FIX.md`、`REMAINING_LIMITATIONS.md`、`NEXT_DECISION.md`。

必须回答：温场跳变是否由排序导致、具体哪个字段/缓存；能量误差有多少是初态漏记；静水误差在哪一条实际更新中产生；修复后β差分是否可分辨。不能只写“新策略改善噪声，建议下一轮继续”。

打包必须含当前依赖头文件、原始字段和可复放脚本；历史实验数据可引用但不能用旧头文件冒充当前实现。未运行、未测量、测量受限都显式标注，不将 `smoke_passed` 写成物理验证通过。

## 12. Cursor最终自检

- [ ] 排序问题被验证或排除，没有把假设直接写成已证实。
- [ ] 稳定ID确实稳定，零时间排序对照通过。
- [ ] T/ΔT/物性/历史力状态按同一粒子身份处理。
- [ ] 临时量、邻域、host指针在排序后正确重算/重取。
- [ ] 没有永久关排序或反复重设温度掩盖释放温场错误。
- [ ] H0建立在最终初温后，prev_energy等同步初始化。
- [ ] 冻结温度的外加能量被记录，不能冒充自由热演化。
- [ ] 静水诊断同状态且可恢复，额外子步没有污染生产数据。
- [ ] 侧壁、深部、液面问题分开，没有笼统归为FS。
- [ ] 静水保持不清零速度/力，不加经验抵消获得PASS。
- [ ] 历史闭盖方案未冒充自由面；p'路线没有被无依据一概禁止。
- [ ] β差分比较的是共同空间矢量场，而非只减两个RMS。
- [ ] 原始状态、完整依赖源和Python图一并交付。
- [ ] 电磁、方腔与独立A-phi成果未被无依据重开或修改。

本轮的终点是获得可信的温度空间分布和可诊断的弱浮力响应，而不是再增加一组静音开关或再生成一张无法归因的流线图。


<!-- ===== END 220_OPHELIE_SORTING_THERMAL_IDENTITY_AND_HYDROSTATIC_REPAIR_CURSOR_PLAN_2026-11-12.md ===== -->

---

<!-- ===== BEGIN 221_OPHELIE_FORCEPRIOR_LIFECYCLE_REPAIR_AND_BUOYANCY_RECOVERY_CURSOR_PLAN_2026-11-12.md ===== -->

<a id="221-ophelie-forceprior-lifecycle-repair-and-buoyancy-recovery-cursor-plan-2026-11-12-md"></a>

# 221_OPHELIE_FORCEPRIOR_LIFECYCLE_REPAIR_AND_BUOYANCY_RECOVERY_CURSOR_PLAN_2026-11-12.md

# Doc221｜ForcePrior 生命周期修复与持续浮力恢复：详细开发任务书

日期：2026-11-12  
执行对象：Cursor / SPHinXsys OPHELIE 工作区  
基线：DOC220_ROUND_PACK_2026-11-12  
唯一主线：修复总非压力力被清空的问题 → 验证持续施力及黏性保留 → 重新验收静水 → 恢复短程自然对流。

## 0. 本轮必须先完成的事

Doc220已经确认并修复温度/ΔT排序错误，预设温场后的H0初始化也已修正。扰动压力路径的无浮力速度明显降低，但同一新路径在每个平流步清空 `ForcePrior`，破坏了框架的增量力更新规则。

本轮先修这个实施问题。暂不将“环流没有增长”归因于高Pr、μ=4、运行时间短；不先增加计算时长或调物性。

### 0.1 本轮完成标准

1. `ForcePrior`在生产推进中保持有效非压力力分量之和，不被重力关闭逻辑整体清空。
2. 恒定浮力在每个实际动量更新中持续存在；连续调用同一力更新不会重复累计。
3. 黏性力在后续处理后仍保留，排序和检查点不破坏当前/历史力状态。
4. 最终修复版本重新跑静水及冻结/释放温场短程对照，输出真实空间速度差。
5. 结论区分力公式正确、总力正确、时间积分正确、静水可用和文献环流对齐。

### 0.2 不扩展的内容

- 保留Doc220扰动压力候选，保持opt-in；本轮不换压力算法，不新增静音开关族。
- 不重新验证方腔，不重做电磁D2、TEAM7和线圈几何。
- 不改独立A-phi/Laplace论文分支，不全局改动其他case默认。
- 不新增ghost粒子，不引入Marangoni/Lorentz/搅拌驱动。
- 不通过加大β、减小μ或清零速度获得可见流线。
- 不先续跑Doc220有缺陷路径的5～20s；从正确初态重新开始。

## 1. 已定位问题与证据

### 1.1 当前调用顺序

在法国入口平流步中，源码包含：

```cpp
fluid_viscous_force.exec();
apply_gravity_force_prior();
if (enable_boussinesq)
    boussinesq_force.exec();
```

`omit_gravity_force`分支当前又包含：

```cpp
hostZeroVecFieldFrench(particles, "GravityForceCK");
hostZeroVecFieldFrench(particles, "PreviousGravityForceCK");
hostZeroVecFieldFrench(particles, "ForcePrior");
```

最后一行不是“只去掉重力”，而是把所有已经登记的非压力力一并抹掉。黏性力恰好在该调用之前刚完成更新。

### 1.2 框架规则

在基线提交的 `force_prior_ck.h` 中：

```cpp
force_prior_[i] += current_force_[i] - previous_force_[i];
previous_force_[i] = current_force_[i];
```

核对位置：`src/shared/shared_ck/particle_dynamics/general_dynamics/force_prior_ck.h`。历史基线原文：

https://github.com/Xiangyu-Hu/SPHinXsys/blob/bc411c43718ca7ecaf4c1ed3e29aefeb4e6090c4/src/shared/shared_ck/particle_dynamics/general_dynamics/force_prior_ck.h

实际工作树可能有本地改动，执行前必须阅读并归档当前版本，不能只引用远端证明本地实现相同。

### 1.3 为什么重复调用浮力不能自动补回

恒定温场下，浮力新值≈旧值。总力被清零后，浮力更新添加的是 `F_new-F_previous≈0`，而不是完整F。因此，浮力专属字段非零、解析公式审计通过，总力仍可能几乎没有浮力。

`force_ok=1`只审计浮力分量，不能证明总力/时间积分正确。黏性分量的生命周期同样要检查。

### 1.4 时序现象

Doc220 2s冻结温场成对数据：轴心δUz在0.068s约1.834e-5m/s，到2s约1.833e-5m/s，几乎不增长。它与“初始受到一次推动、后续持续贡献丢失”的解释相符。

这是重要辅助证据，但不能仅根据曲线外观宣布所有异常都由此造成。必须以本轮逐阶段总力账本和真实积分验证收口。

## 2. 起点、代码范围与历史结论处理

工作区历史分支 `feature/EM-Biot-Savart`，HEAD `bc411c43718ca7ecaf4c1ed3e29aefeb4e6090c4`，Doc220记录有86行dirty。开始前保存当前status、必要diff和源/二进制hash；禁止reset/clean，不自动commit/push。

主入口：

`tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/test_3d_ophelie_french_natural_convection_frozen_q.cpp`

Doc220记录的入口SHA256：`47384f74c48b38775eeb57df3c24629024a5b4fb3c3bd020a308c5df9ac77830`；二进制SHA256：`cd1628a39fec23c676fc754209898168fe98cfe1beaa9524f32f6201bfa84f83`。这些用于定位历史，不保证当前相同。

主要修改入口的重力关闭/初始化/诊断逻辑、必要测试和日志；除非实际发现框架缺陷，不修改 `ForcePriorCK` 的全局增量机制。

### 2.1 继承与暂缓继承

| Doc220成果 | 本轮状态 |
|---|---|
| Temperature/ΔT加入evolving，排序身份通过 | 保留，做一次短回归 |
| H0建立于预设温场后 | 保留，继续检查推进热收支 |
| 电磁/空间Q/功率身份 | 保留，不重开算法 |
| pprime无浮力U_rms约1e-4 | 保留为旧版本测量；修力后必须重新测 |
| “高黏度/时间短导致不增长” | 降为未证实假设 |
| 浮力公式匹配 | 保留，不等同总力正确 |
| 文献自然对流/双涡 | 仍未通过 |

## 3. A：建立总力来源清单和最小旧行为复现

### 3.1 列出全部有效非压力力

建立 `FORCE_COMPONENT_OWNERSHIP.csv`：注册名、当前字段、Previous字段、更新时间、是否贡献ForcePrior、是否随排序、是否在checkpoint恢复。

至少覆盖：Gravity、Boussinesq、黏性inner/contact以及该入口实际启用的其他分量。不要假设黏性一定独立存储一个名为ViscousForceCK的字段；读实际CK实现确认它的注册与更新方式。

压力 `Force` 与非压力 `ForcePrior`分开，不把压力项加入ForcePrior账本。相同物理力不能既由框架登记又由入口手工加一次。

### 3.2 在原有路径采集一次证据

冻结热核、β开启、正常pprime路径，仅推进到至少两个平流步，记录：

1. 黏性更新前；
2. 黏性更新后；
3. 重力处理后；
4. 浮力更新前/后；
5. 每个声学动量更新前/后。

按稳定ID输出ForcePrior、有效各分量当前值和Previous值、Mass、T、Velocity、压力Force、dt及相位。捕获后不在诊断中额外调用浮力/声学更新来“刷新”观测值。

必须证明旧路径在哪一步破坏账本，而不只是做一个完全独立的Python示意。

## 4. B：最小修复设计

### 4.1 首选：启动前固定压力/重力模式

Doc220 pprime工况在构造/初始化时就确定不登记GravityForce。正常每步推进中：

```cpp
if (!omit_gravity_force)
    constant_gravity.exec();
// omit=true: do nothing to ForcePrior or other force histories.
```

这是行为示意，不要求逐字照抄。初始化时确认Gravity当前/Previous均为0、总力账本一致；构造完后不再为了“保险”每步清空任何共享总力。

可以重命名 `apply_gravity_force_prior` 使职责明确；命名修复不是重点，不能只把清零搬进另一个函数。

### 4.2 若必须从已登记重力状态切换

本轮不要求运行中切换。若现有初始化确实登记过g，不能只将Gravity和Previous同时写0却不从总力中撤销原贡献。

使用框架同一增量机制将该力的新值更新为0，或一次性明确从总力撤销已登记的重力并同步该分量历史。不得清空其他分量。记录切换前后总力和每个分量。

最简单可靠的本轮方案是让pprime模式在任何force更新前已固定，避免中途切换。

### 4.3 禁止的临时补丁

- 每平流步清空ForcePrior，然后只重加浮力。
- 每声学步直接 `ForcePrior += F_buoyancy`，造成重复累计。
- 清空所有Previous字段来让增量“重新起效”。
- 将专属浮力字段直接写进Velocity来绕开动量更新。
- 修改全局ForcePriorCK为覆盖赋值，破坏多力组合。
- 用更大β、更小μ或速度清零掩盖生命周期错误。

### 4.4 检索全部路径

搜索该入口及直接helper中所有对ForcePrior、各Previous字段的写操作，标明属于：一次性初始化、模式切换、正常推进、可恢复诊断或废弃实验。

生产路径不得触发诊断清零。若诊断需要临时改状态，必须完整恢复总力、所有相关当前/Previous、速度、ρ/压强/密度变化率及时间状态；不能只恢复Velocity。

## 5. C：总力账本验收

### 5.1 正式不变量

在所有当前分量都已按规定时序更新的时刻：

`F_prior(i) = Σ_k F_registered,k(i)`。

这里的和包括实际已登记的非压力力，不包括压力Force。未更新分量允许保持上次值，但不得凭空漏掉它。

定义逐粒子残差 `R_i=F_prior_i−ΣF_k,i`，同时输出绝对N值和归一值。归一分母使用各分量范数尺度与固定数值floor，不能用可能抵消为0的总力作唯一分母。

float路径目标归一RMS≤1e-5、最大值≤1e-4；double-only代数测试目标≤1e-10。接近零力用明确绝对尺度；不为达到门限增大floor。

### 5.2 三个时序断言

1. omit-gravity处理前后：除实际撤销的Gravity分量外，ForcePrior不变；正常固定模式下应完全不变。
2. 温场不变时再次调用浮力更新：总力不再增加，但原有浮力仍在。
3. 黏性更新后，经重力处理/浮力处理直至声学更新前：该黏性贡献仍包含在总力中。

分别检查至少两个平流周期中的所有声学子步，不只检查初始一次。

### 5.3 真实动量更新检查

声学更新前后直接核对：

`Δv_actual` 与该实际更新使用的 `(ForcePrior + Force)/Mass * dt_update`。

使用源码真实的更新系数/时间间隔，不能因为名称1stHalf就擅自用dt/2，也不能把一整步时间重复算两遍。若还有其他速度更新，逐阶段单列。

记录全域动量增量、外力冲量与壁面/压力作用，但不把壁面容器总动量当作无壁自由加速解析解。

## 6. D：两个小型、实际CK路径的持续施力测试

本节是高影响力生命周期修复的必要测试，不增加新的大规模物理benchmark。使用真实ForcePriorCK及实际动量UpdateKernel；纯Python复刻规则只能辅助，不能替代它们。

### 6.1 D1：恒定力持续积分

最小粒子体，测试模式不执行压力/壁面/黏性Interaction，只使用实际动量速度更新。设质量m、恒定非压力力F0，初速0。沿真实更新时序重复至少100个子步，累计实际dt。

预期：`v=F0/m * Σdt_update`。每次重复exec同一力时，ForcePrior保持F0，不变成0或nF0。

速度相对误差float目标≤1e-3（同时给绝对值）；dt采用一个固定值和一次减半对照，物理总时间相同。这个测试不能推广为法国流场应当自由加速，只验证持续力的积分。

### 6.2 D2：两个分量及一次改变/关闭

使用两个独立ForcePrior通道，先A+B，然后只改变A，再只关闭A，重复更新B。检查总力依次为对应A+B、A_new+B和B。

包含omit-gravity helper的实际调用，证明它不会删除B。不要仅测框架基类而绕过造成错误的入口helper。

### 6.3 D3：实际黏性保留（不另建长时验证case）

在法国同几何设置小幅非均匀初速、β=0/Q=0，执行实际黏性更新及后续helper。比较已计算的黏性分量与最终总力中的贡献，确认没有被抹掉。

计算该阶段离散黏性功率 `ΣF_visc·v`，与当前黏性模型预期的耗散方向核对。压力/壁面/体力功分开，不因总动能短时变化就直接判黏性错误。

只需少量子步；若黏性算子本身的空间精度需要研究，明确另列，不在此扩大成Poiseuille/方腔重验。

## 7. E：排序、静水与能量回归

### 7.1 排序与历史力

复用Doc220稳定ID排序审计，扩展当前/Previous力及ForcePrior。相同ID纯排序前后均不变。临时计算量如果重算，必须在使用前完成。

跑一次含≥3次正常排序的短程（必要时将诊断排序频率明确设小，但正式场仍用既有频率）。关闭排序只允许故障定位，不能作为正式配置。

checkpoint若保留持久力，必须同时保留对应Previous；若重建，明确采用一致初始化协议。至少做一次短resume身份检查，或明确本轮不提供可用resume，禁止错误继承旧checkpoint。

### 7.2 重新验收pprime静水

同Doc220法国几何和物性、Q=0/β=0，正常黏性、排序、密度和边界处理。先0.3s，再到1s；保持过程不清零速度/力。

记录是否仍有settle结束速度清零，并把它与后续1s无人工清零保持分开。不能用末态清零获得静水门。

原约1e-4m/s只是待复核基线。目标全域RMS≤1e-3m/s，侧壁/深部单列，并证明相对原0.02m/s明显改善。修复后结果改变不应被强行拉回原值。

若静水恶化，先查当前/Previous与黏性恢复、步长和边界相容性，不恢复总力清零来“修好静水”。

### 7.3 热量与冻结约束

保留正确H0初始化。释放热核0.3s原残差约1.03%，严格略超≤1%目标，不能继续写成完全PASS。

核对移动粒子的能量计算是否使用一致固定质量、cp与热更新体积；逐项记录Q、Robin/辐射、扩散内部残差。不要用放宽门限消除这0.03个百分点差异。

冻结温场每次重设前后计算约束能量差，累计E_prescribe；初始化热核能量与推进中约束能量分开。若本轮仅做力学冻结诊断、不验热收支，明确标NOT_APPLICABLE，不能打印虚假的热守恒PASS。

## 8. F：修力后的冻结/释放温场对照

### 8.1 固定物理条件

沿用D=.5m、总高=.185m、dp=.015m、μ=4Pa·s、β=1e-5、原声速/稳定步长策略。pprime保持显式开关，机械开口/自由面和侧底壁面身份不变。

本轮温场采用既有径向热核ΔT=200K作为诊断，不称论文真实温场。β0与βon从同一可信静水准备状态开始；冻结与释放工况初态一致。

| 工况 | β | 热处理 | 时长 |
|---|---:|---|---|
| F0 | 0 | 冻结温场，Q/扩散/Robin关 | 0.3s→2s |
| F1 | 1e-5 | 同F0 | 0.3s→2s |
| F2 | 0 | 释放温场，空间Q50kW+扩散+Robin | 0.3s→2s |
| F3 | 1e-5 | 同F2 | 0.3s→2s |

由同一次运行保存检查点，避免从头反复计算。主验收通过且趋势有意义时，F2/F3最多续到5s；不自动20s，不沿用错误总力路径checkpoint。

### 8.2 总力审计嵌入实际法国推进

在前两个平流周期逐声学步审计，此后每个平流步审计摘要，排序前后额外审计。保证声学动量Update实际读取的ForcePrior包含正确浮力与黏性。

除原 `Fz_analytic` 外，输出：

- 浮力分量范数、Previous范数；
- ForcePrior与有效分量和的误差；
- 声学更新前浮力占比与实际外力冲量；
- 按阶段的force丢失/重复计入计数（应为0）。

不允许为输出审计再次调用force.exec()改变被测状态。

### 8.3 空间比较与合理性

必须在共同固定空间网格上输出 `δu=u_beta_on−u_beta_off` 的有符号向量、空间L2、轴心/外围剖面及时间曲线。两个速度RMS相近不能证明无浮力响应。

只要求先证明持续驱动与方向可解释，不预设速度必须单调增长：压力和黏性可使局部速度饱和或改变。不得要求法国容器速度等于无壁自由加速公式。

若修复后δUz仍几乎恒定，先检查声学前力账本、空间δu和压力平衡；这些通过后才讨论黏性时间尺度/温场结构。不能直接回到“时间短，再跑长些”。

### 8.4 温度、Q与机械边界同步输出

- T(r,z)、Q(r,z)、真实子午切片及环向平均速度；
- 温度字段与ΔT一致性、排序前后状态；
- 空间Q随当前位置重采样，不能让Q跟着粒子ID漂移；
- 边界法向速度、液面高度/形变、侧壁滑移诊断；
- 总热收支与动量/力冲量身份。

保持固定坐标和色标，输出有效域mask、覆盖率和插值误差自检。Fig7/8对齐、双涡拓扑不在本轮自动授予范围。

## 9. pprime物理一致性需要保留的检查

低静水速度不自动证明自由面受力正确。记录 `p_total=p_background+p'`，真实自由面总压条件与p'边界对应关系、wall压力外推和密度/EOS约定。

当前新路径以p'推进并省去基准g，只允许省去与静水背景抵消的重力部分；温差浮力与黏性不能一起省去。不要只因无g静水安静，就认定对流物理都正确。

本轮不再调整pprime方程；如果力生命周期修好后暴露边界恢复力问题，提交独立证据和明确未通过状态，不增加隐藏盖板或静音处理。

## 10. 执行顺序与停止规则

| 顺序 | 工作 | 进入下一项条件 |
|---|---|---|
| 1 | 当前源码审计+旧行为2平流步复现 | 明确总力在哪一调用丢失 |
| 2 | 最小修复及账本断言 | 不变量通过，无重复累计 |
| 3 | D1/D2/D3短测试 | 持续力积分和黏性保留通过 |
| 4 | 排序/静水回归 | 状态身份正确，静水达到工程门 |
| 5 | F0/F1冻结温场 | 力学持续驱动可解释 |
| 6 | F2/F3释放温场+热收支 | 温场/能量/空间差分可信 |
| 7 | 有条件续到5s | 不再有实现缺陷主导的迹象 |

若账本失败：只修生命周期，不改μ/β/压力。若账本通过但真实积分失败：检查具体动量更新时序/时间系数。若静水失败：检查恢复的黏性/壁面，不允许清零补丁。若温度跳变复发：停止对流判断，恢复排序审计。

## 11. 测试预算与不再重复的工作

预算固定：旧行为两个平流周期；最小持续力/双通道/黏性保留；一次排序回归；静水0.3/1s；F0–F3各0.3/2s检查点；条件满足时释放工况最多5s。

不新增EM/方腔/TEAM7全套回归。若修改范围只在法国入口，不以“保险”为由重跑无关模块。主机/GPU同步审计仅短程开启，正式性能不能把大量CSV/同步成本算成算法本体。

## 12. 交付文件与结论层次

根目录建议 `docs/ophelie/reference/doc221_runs/`。

| 目录 | 内容 |
|---|---|
| `00_entry/` | 当前源/二进制hash、精度、命令、配置及变化范围 |
| `A_repro/` | 旧总力丢失的分阶段CSV、ForcePrior来源清单 |
| `B_fix/` | 精简diff、初始化/省略重力逻辑、写操作审计 |
| `C_force/` | 总力账本与实际Δv检查、D1/D2/D3 |
| `E_regression/` | 排序、静水、原/新状态对照 |
| `F_convection/` | F0–F3、空间δu数据、Python温度/矢量图和时序 |
| `thermal/` | 初态/冻结约束/释放热收支分项 |

源快照必须包含当前实际 `force_prior_ck.h/.hpp`、Boussinesq、ViscousForceCK、声学更新、温度注册和排序相关实现，不能再次只打包入口CPP。

根报告必须有：

- `ROUND221_CLOSURE.md`：分别说明force分量、总力账本、真实积分、静水、环流响应的结论。
- `RESULTS_TABLE.csv`：真实测量值、阈值、PASS/FAIL/NOT_RUN；不以“低于门限”替代数值。
- `ROOT_CAUSE_AND_FIX.md`：旧问题是否复现、实际修复、哪些历史结论受影响。
- `REMAINING_LIMITATIONS.md`：pprime边界、物性/时间尺度、文献环流仍未验证项。
- `NEXT_DECISION.md`：根据真实空间响应，只给一个下一主方向。

### 12.1 本轮不得混用的状态

| 状态 | 意义 |
|---|---|
| FORCE_COMPONENT_VALID | 浮力公式等分量正确 |
| FORCE_LEDGER_VALID | 总力正确包含所有有效分量 |
| FORCE_TIME_INTEGRATION_VALID | 实际速度更新正确使用持续力 |
| HYDROSTATIC_RECHECK_PASS | 修力后的静水满足工程门 |
| BUOYANCY_RESPONSE_RESOLVED | 共同空间场中响应可分辨且方向可解释 |
| FRENCH_LITERATURE_VALIDATED | 本轮不预设授予 |

## 13. Cursor最终检查清单

- [ ] 正常省略重力分支不再清空ForcePrior。
- [ ] 不手工重复加浮力，不清空所有Previous“重新起效”。
- [ ] 初态/切换/正常更新/诊断四种生命周期分清。
- [ ] 总力账本在第二个及后续平流周期依旧通过。
- [ ] 黏性更新后的贡献到动量更新前仍存在。
- [ ] 实际动量dt/半步系数经过源码核查。
- [ ] 恒定外力产生持续冲量，而非只第一步有效。
- [ ] 两个分量独立改变/关闭不会破坏另一个。
- [ ] 排序和checkpoint历史状态一致。
- [ ] 修力后静水重新验证，未强行继承旧低速度PASS。
- [ ] 释放热残差1.03%未被写成严格≤1%通过。
- [ ] 冻结约束能量与初始热核能量分开。
- [ ] F0–F3有共同空间δu矢量图，不只减两个RMS。
- [ ] 没有通过调μ/β/静音开关掩盖实现问题。
- [ ] 电磁、方腔和独立A-phi工作未被无依据改动。

本轮收口应给出一句可核查的结论：持续浮力和黏性是否已正确保留并参与每个动量更新；只有它成立，才继续讨论法国环流的物性、时间和文献边界差异。


<!-- ===== END 221_OPHELIE_FORCEPRIOR_LIFECYCLE_REPAIR_AND_BUOYANCY_RECOVERY_CURSOR_PLAN_2026-11-12.md ===== -->

---

<!-- ===== BEGIN 222_OPHELIE_FREE_SURFACE_PRESSURE_AND_FORCE_LIFECYCLE_CURSOR_PLAN_2026-11-15.md ===== -->

<a id="222-ophelie-free-surface-pressure-and-force-lifecycle-cursor-plan-2026-11-15-md"></a>

# 222_OPHELIE_FREE_SURFACE_PRESSURE_AND_FORCE_LIFECYCLE_CURSOR_PLAN_2026-11-15.md

# 第222轮：自由面压力相容性与外力生命周期修复执行计划

日期：2026-11-15  
适用项目：冷坩埚招投标／法国 OPHELIE 自然对流复现  
执行对象：Cursor；基于 DOC221_ROUND_PACK_2026-11-15 的实际代码继续开发

## 0. 本轮必须完成的事

本轮目的：弄清恢复持续浮力以后，为什么出现强烈上涌和液面抬高，却没有可信的自然回流；修复已识别的尖端力重复累加，并通过独立测试判断扰动压力自由面处理是否缺少重力恢复作用。

本轮不是以获得某张好看的环流图作为成功标准。必须得到以下一种明确结果：

1. 找到并修复压力／自由面不相容，未经人为速度约束的冻结温场计算产生可信回流；或
2. 将问题缩小到一个可复现的算子或边界状态错误，提供同一状态的分项证据和最小复现，说明下一次只需改什么。

禁止以新增阻尼、改物性或改判据替代这两种结果。

执行顺序：A 固定底座 → B 外力生命周期 → C 压力路径审计 → D 无热自由面测试 → E 只修被证实的缺陷 → F 冻结温场回流 → G 短程热耦合回归。

**B、C可以独立开展；D使用不含尖端干预的底座。E必须依据C/D证据，F必须在D通过后进行。不要同时改多套求解公式。**

## 1. 已有证据、待证假设与不能继承的结论

### 1.1 已有有效进展

DOC221 修复了 omit-gravity 分支每个平流步清空共享 ForcePrior 的错误；重力通道、浮力通道和黏性通道应保留这个修复。

| 证据 | DOC221 结果 | 本轮处理 |
|---|---:|---|
| 基础外力账本相对 RMS 残差 | 约2.8e-6 | 作为回归参考，不重写力框架 |
| 冻结热核、β=0，2s速度RMS | 约2.79e-5 m/s | 保留静止背景对照 |
| 同一热核、β=1e-5，2s速度RMS | 约1.16e-2 m/s | 持续浮力已恢复 |
| 2s轴心速度的β开关差 | 约+9.62e-3 m/s | 已非“浮力完全不起作用” |
| 短程释放温场能量残差 | 约0.85%～0.87% | 是短程回归证据，不等于稳态验证 |

### 1.2 当前真实异常

加高墙沿后，5s不再越缘，但最高粒子位置仍相对静置后增加约145mm；3s增加约49mm。容器玻璃名义高度约185mm。必须区分局部尖峰、平均液面、主体体积变化，不能用z_max单独宣布液体整体膨胀或体积守恒失败。

同时轴心上升、外围平均竖直速度仍向上，未经人为干预尚未得到清晰回流。

### 1.3 必须降级的旧结论

- N4/N5依赖全局平均Uz扣除和尖端外力，不能作为物理自由面成功基线。
- N4/N5数值相同只说明当时的筛选集合相同，不能称为“物理等价验证”。
- J中的自由面预算是近似运动学解释，不是独立的质量／体积守恒证明。
- 加高墙避免越缘是几何遏制措施，不是液面动力学正确性证明。
- Ur上下符号改变、轴升缘降、流函数出现正负值，都不能单独证明闭合双环流。
- 方腔验证的有效结论保留，但不能自动覆盖三维可变形自由面与本case的p′处理。

### 1.4 本轮两个优先问题

**问题一，代码可直接定位：** `hostApplyTipSoftForcePrior`直接执行`ForcePrior += m*a_tip`，在声学步、平流步内多次调用，包括用于打印活动情况的调用；未见匹配的历史贡献扣除。需要以当前完整调用链确认并修复。

**问题二，强假设、尚未证实：** p′路径去掉完整重力后，固定参考静水压力及自由面密度处理，可能没有正确保留液面变形的恢复作用。不得在未读声学实现和做D测试前宣称“根因已经完全确认”。

## 2. 范围与禁止事项

### 2.1 本轮不做

- 不修改独立A-phi／Laplace分支。
- 不重开TEAM7、电磁离散、趋肤、屏蔽、线圈和功率标定工作。
- 不更换整个流体框架，不引入新δ-SPH／ghost布局。
- 不同时扫描μ、β、温差、分辨率、声速、密度模式、阻尼系数。
- 不把减全局平均Uz、杀死上升速度、限位z、尖端阻尼默认为正式物理条件。
- 不将ρ_sum<ρ0直接当成可靠的真实自由面识别器；固壁附近也可能支撑不完整。
- 不直接启动20s或更长法国耦合计算。
- 不为了通过改拓扑阈值、删掉不利区域或只展示β差场。

### 2.2 代码隔离

沿用现有法国入口进行小范围修复；新建一个轻量验证入口或复用现有静水／液面入口，仅用于D测试。不要为每个开关复制一个大CPP。

建议验证入口名称：`test_ophelie_free_surface_pressure_consistency`。名称仅为建议；若已有等价入口则直接复用。

验证入口应调用与法国case相同的压力、密度、壁面、力更新和位置推进实现。若简化成完全不同的求解器，不能解释法国case。

优先CPU双精度小规模完成算子审计和快速判别，再用相同路径做GPU短程回归。最终法国入口优先GPU；CPU慢但正确的对照可保留。报告清楚实际在哪个设备执行。

## 3. A：固定一个干净、可追溯的底座

### A1. 保存版本与调用链

记录git SHA、dirty diff、编译命令、精度、执行设备、线程／GPU信息、输入文件哈希、所有最终生效参数。禁止只保存命令行而不保存别名开关展开后的配置。

本轮输入包只包含法国入口CPP和patch，缺少部分共享算子正文。因此必须从执行服务器保存实际编译使用的以下实现及相关diff：

- 两个声学半步及压力／密度状态更新；
- 内部与壁面Riemann状态构造；
- ForcePrior增量更新及Gravity／Boussinesq／Viscous通道；
- FreeSurface密度正则化、BackgroundPressure及PerturbationMomentum；
- 接触壁面法向和参考加速度处理；
- transport velocity、排序／evolving字段注册及位置更新。

明确指出：哪些是当前分支自定义实现，哪些沿用上游。不要仅根据变量名称推断物理含义。

### A2. 固定参考状态

采用本轮已有相同几何、粒子生成和静置流程；保留DOC221基础力修复、热字段排序修复、能量H0修复。建立无尖端干预的基线R0。

R0强制关闭：zero-net-Uz、hard/soft lid、tip kill、tip damping、incomplete-force cancellation，以及其他直接修改Velocity／Force来压制表面运动的诊断开关。若原底座依赖其中某项，应明确列出，不能隐式继承。

墙沿加高可以用于避免短程越缘，但两条压力路径必须相同；记录液体以上墙体是否被纳入接触与密度支持。不要一边加高墙，一边把它当成完全不影响求解的变化。

相同布局、扰动种子、温场、μ、β、c0、时间步策略、TVF开关和输出时刻。压力路径对照允许各自使用物理一致的静水密度初始化，不能强行让本来不同的状态量数值相同。

## 4. B：修复尖端外力生命周期；禁止继续借它凑环流

### B1. 先给出完整调用证据

定位`hostApplyTipSoftForcePrior`全部调用，逐处标明：调用频率、是否用于日志、在声学力积分之前还是之后、是否同一时刻重复执行、下一步是否有旧力撤销。

核对当前源码中每声学步调用、平流分支记录调用、黏性／浮力更新后的调用。不要只修一处。

以一个固定状态执行函数1、2、10次，输出调用前后ForcePrior和尖端贡献，证实是否重复累加。旧行为只保留在独立复现模式中。

### B2. 正确更新约定

若保留该诊断功能，必须按框架外力通道规则实现：

```text
new_tip_force = m * a_tip(current_state)
ForcePrior += new_tip_force - previous_tip_force
previous_tip_force = new_tip_force
```

可以复用标准ForcePrior dynamics；不要重建一套共享总力清零规则。对应字段必须随粒子排序，并按框架要求正确同步host/device。

所有粒子都应写入新尖端贡献，离开作用区或Uz≤0时写零并撤销旧贡献。禁用该功能、切换阶段、重启时也必须有一致行为。

日志函数只能读状态，禁止通过再调用施力函数获取n_touch。触碰数量可在施力时记录，或由无副作用的统计函数计算。

**不要给力额外乘dt。** ForcePrior存的是力，时间积分由声学推进处理；除非实际框架定义不同，此时先写清单位与调用契约。

### B3. 必过检查

| 检查 | 预期 |
|---|---|
| 同一状态重复执行10次 | 与执行1次得到相同总力 |
| 作用区内 → 作用区外 | 原尖端贡献撤销 |
| 上升 → 静止／下降 | 单向阻尼贡献归零 |
| 开启 → 关闭 | 不留下历史力 |
| 调用日志前后 | 任意粒子物理状态不变 |
| 排序前后按OriginalID比较 | 贡献对应同一粒子 |
| 重力＋浮力＋黏性＋尖端账本 | 和等于ForcePrior，误差符合精度 |

双精度小测试目标相对误差1e-10量级，单精度使用与运算规模相符的容差；近零量必须同时使用绝对误差，不能用接近零的分母制造大相对残差。

另做一个常量加速度的实际积分检查，确认更新次数变化不会使加速度翻倍。只测最终字段相等还不足以覆盖调用时序。

B通过后，N4/N5仍为诊断功能。本轮主线继续关闭；**不重扫c、k、tip_r，不尝试恢复旧的PASS图。**

## 5. C：逐公式核对扰动压力与自由面条件

### C1. 写出实际求解方程

交付`PRESSURE_CONVENTION.md`，至少明确：

1. Density、Pressure、BackgroundPressure各自定义；EOS输入及输出是什么。
2. 动量内部相互作用到底用p、p′、p_h还是组合；每个半步分别说明。
3. 连续性使用的Riemann速度／压力状态；是否与动量约定相同。
4. 壁面p、ρ、法向速度和参考加速度如何构造。
5. 自由面是显式边界状态，还是只依赖不完整支持和密度截断。
6. 重力在哪项抵消、是否仍残留、是否被重复计入。
7. 固定参考高度与实际液面高度分别在哪里使用。
8. TVF或其他正则化是否额外使用背景压力；不能混淆输运背景压力和静水压力。

为每条公式给出实际文件、函数和调用位置。只列开关表不算完成。

### C2. 数学核对基准

取z向上，参考密度ρ0，忽略表面张力，以大气压为表压零点：

\[
p_{total}=p_h+p',\qquad p_h=\rho_0 g(z_{ref}-z).
\]

若在实际液面z_s上要求p_total=0，则：

\[
p'(z_s)=\rho_0g(z_s-z_{ref}).
\]

因此，“自由面保持ρ≈ρ0，即p′≈0”只在参考平面附近适用，不能直接推广到任意变形液面。

当前`hostApplyFrenchBackgroundHydroPressure`使用`max(0,z_free-z)`。若截断导致p_h在参考液面以上梯度变零，而完整重力仍被全域移除，必须逐区域说明抵消关系还是否成立。

**上式是核对约定的依据，不是要求把一个公式直接粘贴到表面密度上。** 还必须考虑实际EOS、负p′可表示性、截断／密度下限、Riemann状态和壁面项。否则可能修了边界值，却破坏连续性或制造拉伸不稳定。

若代码把p_h仅作为记账字段，从未以相容形式进入实际积分，也必须明确指出。不能因为字段存在，就认为恢复作用存在。

### C3. 同一状态的力分解

取无热的小幅液面变形状态，在完全相同粒子几何上，输出：内部压力力、壁面压力力、重力、黏性、浮力、其他正则化力、净加速度。

注意有些力只存在于声学Force而非ForcePrior。只审查ForcePrior不能判断总加速度是否正确。若半步有Riemann耗散，单列或说明无法拆分的组合项。

对总压力路径与p′路径比较它们代表的**总物理加速度**；不要比较两个不同定义的Pressure字段数值是否相等。

## 6. D：先做无热自由面测试，不先跑法国热核

### D0. 公共设置和测量

β=0，关闭热源、传热和EM，不施加任何尖端约束，不扣全局平均速度。固定黏度、几何和位置更新逻辑。

优先用同一三维小圆柱／现有法国粗粒子布局，以免测试成功却没有覆盖实际边界。如果现成小水槽入口更快，可作快速预检，但至少一组必须回到法国圆柱的相同实现。

所有测试记录：U_rms、U_max、动能、质心z、p_total及p′、ρ误差、z_max、表面高度分位数、按半径分箱的表面高度。表面高度必须说明重建方法、有效bin数和空洞，不把最高一个粒子当成整个液面。

### D1. 平液面静水保持

两个物理一致路径：

- T：总压力＋完整重力＋相容静水密度／压力初始化；
- P：p′＋背景静水压力＋对应自由面／壁面约定，保留当前实现作修前对照。

首跑0.3s，正常后到1s；若仍在静置衰减，报告衰减趋势及实际时间尺度，不用任意时刻阈值宣布不稳定。

通过要求：无持续动能增长、无大规模定向漂移、压力随深度的趋势正确、内部净加速度小；以DOC221静水速度量级作为参考，而非跨设置硬性通用阈值。

### D2. 小幅液面变形的重力恢复

这是本轮最关键测试。静水安静不能替代它。

建立可追溯、总体积不变的小幅轴对称液面变形，例如中心略高、外围略低：

\[
\eta(r)=\epsilon[1-2(r/R)^2],\quad \int_0^R2\pi r\eta(r)dr=0.
\]

用解析变形域重新生成粒子或采用经过验证的初始化；记录离散体积偏差与初始密度。不要随意把整列粒子拉伸而仍假定原密度／体积完全一致。

ε应明显高于该分辨率的表面识别噪声，同时ε/H较小；先用一个可分辨幅值，记录ε/dp和ε/H。若粗粒子无法兼顾这两点，只对这个小测试提高分辨率，不启动法国全域分辨率扫描。

对照T与P，初始形状、物性相同，压力各按自身定义初始化。记录早期总加速度和随后形状／动能演化。必要时增加g=0对照，识别纯压力松弛或粒子正则化造成的表面变平。

**高黏玻璃可能表现为过阻尼恢复，不要求一定出现振荡。** 本轮不套用不相容的无黏液面波周期作硬指标。

判断重点：

- 是否存在与重力相关、推动形状恢复的响应；
- T、P描述同一物理问题时是否给出相容趋势；
- 原P是否在越过参考高度时出现异常失衡；
- 两条路径的差异能否对应C中的某个实际项。

不能要求每一个表面粒子瞬时都向下；恢复运动由非局部压力和连续性共同决定。以模态幅值、径向表面形状及力分解判断。

### D3. 条件式时间步复核

只对D2的决定性差异做一次dt和dt/2对照，物理输出时刻相同。若差异随步长明显变化，先查积分与调用频率；不要立即把它归为自由面物理错误。

### D结果决策表

| 结果 | 下一步 |
|---|---|
| T静水和恢复正常，P恢复异常 | 优先修P的压力分解／表面／壁面一致性 |
| T、P静水均失败 | 停止热核计算；查共同的压力、壁面支持、法向和初始化 |
| 静水均正常，恢复均异常 | 查自由面实现、几何初始化及共同位置／连续性更新 |
| T/P均通过，法国热核仍异常 | 转F的受热动力学与控制体收支，不能继续认定缺重力恢复 |
| 步长减半改变结论 | 优先查调用与推进时序；不扫描更多边界参数 |

## 7. E：只修证据指向的那一条路径

### E1. 首选：恢复与总压力物理等价的p′路径

若C/D证实自由面转换不相容：修复动量、连续性与边界状态的整套约定。可能涉及固定参考p_h的有符号延拓、实际表面p′条件、壁面压力外推、密度正则化的覆盖关系；按审计结果选必要项。

不能把`max(0,depth)`删掉就宣布完成，也不能全域冻结ρ来掩盖连续性。修后必须重跑D1/D2及原力账本。

### E2. 有证据的备选：总压力＋完整重力

若T路径通过而P修复涉及大面积自定义代码，可采用T作为本轮法国case底座，先保证准确。必须记录静水残余与热浮力信号之比；若静水噪声与热浮力同量级，不得因为出现大运动而通过。

不要在本轮同时维持三四套压力算法开发。

### E3. 文献平顶近似只作为明确建模选择

如果核对原文确认其自然对流采用固定几何上的平坦、滑移、不穿透顶面，可保留这种工程模型作为后续正式路线候选。不能仅凭本计划假定原文细节已确认。

实现必须在Riemann／壁面状态中施加相容法向条件，并区分其压力条件与可变形大气自由面。不可用每步Uz清零、z截断、全局减平均速度替代。

DOC221硬顶实验出现极性反转，已经说明“加顶就能解决”不成立。不要重扫dens_mode；如走此路线，先通过该顶面下的静水和简单定向受力测试，再进入热核。

E3不是默认并行分支。仅当文献依据明确且E1/E2不合适时，写出选择理由再实施；已有授权范围内无需为常规实现另行等待确认。

## 8. F：恢复冻结温场，验证真正的回流

### F1. 严格对照

使用DOC221同一冻结热核、ΔT=200K、同一μ和β；它是诊断温场，不是法国论文温度场。

至少两组：β=0与β=1e-5。采用D已通过的同一压力路径、同一液面模型，不含任何尖端／全局速度干预。

计算先0.3s，再2s；结果正常但尚未形成完整响应才到5s。时间不是越长越好，以动能、液面及压力诊断判断是否值得继续。

瞬时没有闭合回流不能直接判失败；但若已出现持续大幅液面增长、控制体收支失配或压力异常，应立即停止，保存首次异常前后快照。

### F2. 质量、体积与通量必须独立测量

固定粒子质量之和恒定只证明没有丢粒子，不证明速度场满足近似不可压缩约束。

同时输出：

- 总质量和Σm/ρ；后者是密度推算体积，不能单独代表几何占据体积。
- 从粒子占据域／表面重建得到的体积及不确定性。
- 三个预先固定高度截面的向上、向下和净体积通量。
- 截面上方的控制体储液变化，与净通量作比较。
- 固壁相对法向速度／近壁通量估计，以及测量层厚度。

对变形自由面，截面净通量可以非零，应满足上方储液变化。对固定不穿透顶面且近似不可压缩的封闭模型，截面净通量应趋近零。两种模型不能套用同一“净通量必为零”条件。

轴对称通量采用`Q=2π∫r u_z dr`时，应说明环带有效覆盖、方位平均方法和空bin。另用原始3D粒子控制体收支交叉检查，不让同一插值结果同时充当计算值与验证值。

禁止通过减平均Uz使通量过关。全域平均速度不是局部散度，也不是任意截面的净流量。

### F3. Python图与流函数检查

不依赖ParaView人工判断，自动生成：

1. 同一r-z剖面的绝对温度、u_r、u_z、速度矢量和覆盖mask。
2. β=0、βon绝对速度，以及两者差场；三者分别标注。
3. 早期／2s／最多5s的液面与速度时间序列。
4. 轴心上升带、内部上升带、外缘下降带的固定采样区统计。
5. 流函数等值线和闭合轮廓；只在数据覆盖与连续性允许时解释。

轴对称流函数约定可取`u_r=-(1/r)∂ψ/∂z`、`u_z=(1/r)∂ψ/∂r`，轴线使用相应极限处理。用两条积分路径或离散最小二乘重建，报告重建速度残差和路径不一致程度。

若连续性不满足，流函数只是近似投影；不得靠它证明原始速度场已经守恒。三维明显不轴对称时，应报告不同方位剖面差异，不能仅用方位平均掩盖异常。

环流通过至少需要：绝对速度中的上升与下降路径、上下径向返回、有效区域内闭合等值线，以及合理的控制体收支。Ur_flip仅作辅助量。

只有主环流通过后，再讨论法国文献的主次涡位置、旋向和强弱。冻结热核不是原文温场，本轮不强制它复现Fig7双涡。

## 9. G：有条件进行短程热耦合回归

仅在D通过、F未显示明显非物理漂移后执行。

沿用现有空间Q、扩散与边界换热，不改EM求解和功率。先0.3s，必要时1s；核对温度随粒子排序、μ(T)更新、H0和边界能量账本是否仍正确。

冻结热核与释放温场必须分目录，不能混用结果。记录Q实际总功率、是否固定Q、是否更新σ(T)、实际热边界与温度单位；未执行的耦合项明确写未启用。

这个测试只授予短程耦合回归通过，不授予法国稳态热场或Fig7/8复现通过。

## 10. 最小执行矩阵与停止规则

| 编号 | 内容 | 前提 | 规模限制 |
|---|---|---|---|
| B0 | 旧尖端力重复调用复现 | A | 固定状态，不长跑 |
| B1 | 修后力通道与实际积分检查 | B0 | 小规模CPU，随后GPU短检 |
| D1-T/P | 平液面静水 | C | 各0.3→1s |
| D2-T/P | 无热液面变形恢复 | D1可解释 | 按恢复时间尺度的短窗 |
| D3 | 决定性结果dt减半 | D2有差异 | 仅复核关键一对 |
| E-R | 一个修复候选回归 | D定位 | 不并行新增算法 |
| F0/F1 | 冻结温场β开关 | D通过 | 0.3→2→最多5s |
| G0/G1 | 释放温场短回归 | F可信 | 0.3→必要时1s |

若某测试无法在计划短窗获得判别，先依据g、尺寸、黏度估算响应时间并报告；允许合理调整该测试时长，不允许用统一20s替代分析。

出现NaN、密度显著失真、动能异常增长、粒子越域、力账本失配、输出覆盖严重不足时停止该组，保存首次发生时刻和前后快照。阈值在启动前写进配置，不根据最终图形修改。

**每次代码修复最多对应一个明确假设。若连续两个候选修复没有改变它所针对的独立测试，应停止该方向，提交最小复现与分项力证据，不继续加开关。**

## 11. 验收状态必须分开

| 状态 | 授予条件 |
|---|---|
| FORCE_LIFECYCLE_FIXED | 重复调用／退出作用区／日志无副作用／排序／实际积分均通过 |
| HYDROSTATIC_REST_PASS | 无热平液面静水通过 |
| FREE_SURFACE_RESTORING_PASS | 可变形模型的D2有可信重力恢复与步长复核 |
| FIXED_TOP_BOUNDARY_PASS | 若选择固定顶面，边界法向、压力与静水独立通过；不冒充上一项 |
| BUOYANCY_RESPONSE_PASS | β对照中出现方向合理且持续的物理响应 |
| RETURN_FLOW_VERIFIED | 未人为改速度的绝对速度场中回流与控制体收支可信 |
| THERMAL_SHORT_REGRESSION_PASS | 释放温场短程温度／能量／排序回归通过 |
| FRENCH_LITERATURE_VALIDATED | 本轮默认不授予；需原文温场、边界及Fig7/8定量对比 |

旧N4/N5报告保留以供追溯，但新增说明：其结果受尖端力生命周期和全局速度干预影响，不能作为正式物理验证。

## 12. 交付目录与必须回答的问题

建议输出：

```text
doc222_runs/
  A_baseline/       最终配置、版本、源码清单、调用链
  B_force/          修前复现、修后通道与积分检查
  C_pressure/       实际方程、压力字段流向、边界状态、分项力
  D_free_surface/   静水、变形恢复、关键时间步复核
  E_repair/         仅一个被证据支持的修复及回归
  F_convection/     冻结温场β对照、Python图、控制体收支
  G_thermal/        满足前提才产生
  ROUND222_RESULT.md
  NEXT_DECISION.md
```

每组保存完整启动命令、最终生效配置、stdout/stderr、原始CSV、必要粒子快照、Python脚本、有效覆盖信息和通过／失败理由。统一SI单位，温度明确K，ρ、μ、ν不得混用。

打包包括实际修改过的共享源文件和patch，不能再只有入口CPP。保存编译目标和依赖路径，便于检查是否编译了新代码。

最终报告必须直接回答：

1. 尖端力是否曾随调用次数累加，退出作用区是否留下残留？修复证据是什么？
2. 总压力和p′路径各自在哪里计入重力与自由面压力？
3. 当前p′是否缺少变形液面的恢复作用？证实、否定，还是证据不足？
4. 哪个测试首次区分两条路径？对应哪一项实际力／状态？
5. 不扣平均Uz、不加尖端阻尼，是否已经出现真实返回路径？
6. 液面抬高是局部尖峰、整体储液变化，还是测量／连续性异常？
7. 现有结论覆盖了什么、没有覆盖什么？下一步唯一最高优先问题是什么？

禁止只写“已收口／PASS”。每一结论对应原始结果和定义。若失败，失败本身有价值，但必须比本轮开始时更精确地定位问题。

## 13. 给Cursor的最终执行约束

本轮已经授权在OPHELIE隔离范围内修复并运行上述必要验证。遇到常规实现选择，依据测试证据推进，不需要反复等待确认。

但不要擅自把一个不同物理模型包装成法国自由面复现，不要因为修复后环流图变差就恢复已知有问题的旧阻尼，不要覆盖已有有效验证。

**先证明力更新正确，再证明边界能恢复；然后才判断热浮力能否形成回流。只修被独立测试证实的问题。这是本轮结束反复试错的核心。**


<!-- ===== END 222_OPHELIE_FREE_SURFACE_PRESSURE_AND_FORCE_LIFECYCLE_CURSOR_PLAN_2026-11-15.md ===== -->

---

<!-- ===== BEGIN 223_OPHELIE_FORCE_RESET_REPAIR_AND_PAPER_TOP_BOUNDARY_CURSOR_PLAN_2026-11-17.md ===== -->

<a id="223-ophelie-force-reset-repair-and-paper-top-boundary-cursor-plan-2026-11-17-md"></a>

# 223_OPHELIE_FORCE_RESET_REPAIR_AND_PAPER_TOP_BOUNDARY_CURSOR_PLAN_2026-11-17.md

# 第223轮整改与推进任务书：统一外力重置，对齐论文顶面，恢复可信自然回流

日期：2026-11-17  
项目：冷坩埚招投标／法国OPHELIE自然对流复现  
依据：DOC222_ROUND_PACK_2026-11-17及本轮代码审查  
执行对象：Cursor；本文为下一轮唯一主任务书

## 0. 给Cursor的执行摘要

用户临近验收。禁止再以大量压力开关、局部阻尼和物性扫描代替定位。本轮只允许两条有明确边界的工作：

1. **必须整改：** 修复阶段切换／变形／诊断中共享总外力与历史通道不一致的问题，重新判定受影响的旧测试。
2. **验收主线：** 建立固定平坦、滑移、不穿透顶面的论文对应模型；先用冻结温场证明无人工速度修正的回流，再进入接近论文的温场和空间热源。

可变形自由面的全面修复降为后续问题，**不得作为固定顶面路线的通用前置门槛**。固定顶面仍须通过自己的压力、连续性、滑移和不穿透测试。

固定顶面是对论文固定计算域和顶部slip描述的工程建模解释，不得声称作者公开了全部SPH实现细节。

本轮不承诺一定完成文献双涡；必须至少完成根因整改和一个可信模型的可解释结果。失败要定位到实际状态／算子，不能只交付“还需新一轮研究”。

## 1. 与222任务书的关系及授权边界

本任务书覆盖222中以下执行优先级：

- 不再要求先完成可变形自由面的所有恢复测试，才能开始论文平顶模型。
- 不继续E1/E2a/E2b局部液面／BackgroundPressure变体扫描。
- 不接受“p′方法结构不可能”“总压力路线已证伪”作为已证实结论。
- 已修好的温度排序、能量H0、omit-gravity不清共享总力、tip增量更新保留。
- 现有方腔、电磁有效验证保留，不因法国case的问题推倒重做。

允许直接在OPHELIE隔离范围内修复、编译和运行必要测试，无需为普通实现选择反复询问用户。不得改独立A-phi／Laplace研究分支，不推送或修改上游main。

若修改共享CK算子，新行为必须由清晰参数／独立策略限定，默认保持其他case行为；避免把法国case的特殊处理写成全局默认。

## 2. 本轮审查确认的事实与待证判断

### 2.1 有效事实

DOC222的tip自检证明旧函数存在重复累加，修后重复执行幂等、离开作用区撤销、日志无副作用。它只验证该诊断力通道，不代表所有阶段切换都正确。

DOC222平液面p′路径速度RMS约4.8e-5 m/s；裸总压力静水约0.019 m/s，仍有真实噪声欠账。无热变形实验异常，F/G未运行。

### 2.2 新定位的P0问题

入口CPP在变形之后执行：

```cpp
hostZeroVecFieldFrench(glass_particles, "ForcePrior");
// 重建几何、处理密度等
if (!local_cli.omit_gravity_force)
    constant_gravity.exec();
```

该分支未同步重置`PreviousGravityForceCK`。静置阶段已执行重力，故按框架增量规则：

\[
F_{prior}\leftarrow F_{prior}+F_{g,new}-F_{g,previous}.
\]

恒定重力新旧相等，总力清零后再次exec不能自动恢复它。现有代码构成明确的状态不一致；**其对具体运行丢失多少力，必须以新账本直接验证，不只根据运动图倒推。**

实际基线规则来源：
`src/shared/shared_ck/particle_dynamics/general_dynamics/force_prior_ck.h`，基线SHA `bc411c43718ca7ecaf4c1ed3e29aefeb4e6090c4`。

参考链接：
https://github.com/Xiangyu-Hu/SPHinXsys/blob/bc411c43718ca7ecaf4c1ed3e29aefeb4e6090c4/src/shared/shared_ck/particle_dynamics/general_dynamics/force_prior_ck.h

执行服务器若有不同实现，以实际编译源码为准，记录差异。

### 2.3 旧变形实验的初始化污染

`hostApplyDoc222AxisymmetricFsDeform`采用整列竖直伸缩；ε=0.02m、H=0.185m，局部Jacobian约为1±0.108。总体积积分为零不代表局部密度／粒子体积一致。

固定质量、沿用体积字段，再冻结ρ或密度正则，可能产生不相容初态。E2a首个输出动能约0.93J，必须先解释，不能只看1s结果。

### 2.4 p′结果不能过度解释

声学初始化`Pressure=EOS(ρ)+BackgroundPressure`；PM=1时动量使用`Pressure-BackgroundPressure=EOS(ρ)`。

所以只改p_h形状不改变p′动量，是该分解的代数性质，不能直接证明扰动压力方法不可用。可变形表面的总压力条件尚未相容实现，这与固定不穿透顶面是不同模型。

模式反号、动能从零增长也不是单独的失稳判据；应看总能量、边界功、耗散、初始误差和形状演化。

## 3. 文献目标：什么是原文明示，什么是本轮选择

文献：Numerical modelling of natural convection in molten glass heated by induction，Chemical Engineering and Processing 47 (2008) 449–455，DOI 10.1016/j.cep.2007.01.013。

| 项目 | 原文或现有依据 | 本轮执行约定 |
|---|---|---|
| 流动 | 层流、稳态、Boussinesq | 用瞬态SPH推进逼近稳定响应，不把20次耦合迭代解释为20s |
| 固壁 | no-slip | 侧壁和底面保持无滑移 |
| 顶部 | slip | 固定平坦顶面、不穿透、切向无黏性牵引作为本轮工程近似 |
| 顶部热边界 | 对流＋辐射 | 流体顶面虽然固定，热损失仍保留 |
| 主环流 | 沿轴心向下闭合 | 不把中心热核上升图当作论文目标 |
| 外围环流 | 冷壁下降形成外围涡 | 在论文驱动条件下检查外侧下降和内部上升带 |
| 温场／热源 | 上部外围加热、表面和侧壁冷却共同驱动 | 先做明确标注的诊断温场，再接近论文条件 |

模型对应说明必须引用原文2.1.1、2.1.2及Fig7附近文字。禁止宣称原文给出了本轮的镜像速度、EOS密度边界或压力偏置实现。

## 4. 必须遵守的范围限制

1. 不增加全局平均Uz扣除、tip阻尼、Uz kill、粒子z限位来让主线通过。
2. 不将`DensitySummation<ρ0`视为自由面的唯一可靠识别；靠固壁同样可能支撑不完整。
3. 不复制DOC221的closed-lid/dens-mode/ContSilent开关包作为新主线。
4. 不通过同时设置顶部p=0和u_n=0来过度约束固定顶面。
5. 不改μ、β、几何、功率来掩盖边界或力账本错误。
6. 不重开TEAM7、趋肤、屏蔽、线圈实体细节和新EM离散。
7. 不把冻结温场验证叫作完整电磁加热自然对流复现。
8. 不把smoke_passed、进程退出0或JSON旧字段passed=1当成物理通过。

## 5. 代码结构：只新增一个清楚的模型入口

建议在现有法国入口加一个显式模型参数，例如：

```text
hydro_model = paper_flat_slip_top
pressure_model = perturbation_fixed_top
```

名字由仓库规范决定，但必须是单一配置入口，不要求用户组合十余个历史开关。所有最终生效参数写到JSON。

旧可变形模型保留独立标识用于复现，不能静默更改其含义。非法组合启动即报错：平顶模型＋tip、zero-net、z-clamp、ContSilent等诊断开关不得同时启用。

尽量复用同一求解核心和已验证的方腔压力／连续性路径。不要复制当前近7000行CPP制造另一套无关联入口；可把通道重置、状态审计、模型配置和Python诊断拆成小模块。

优先GPU；允许小规模CPU双精度验证。任何CPU回退明确记录，不能称全部GPU已通过。

## 6. A：冻结底座和证据身份

### A1. 必须保存

git SHA、dirty diff、实际编译目标、编译器、Real精度、设备、运行命令、输入reload哈希、最终参数、修改过的共享头文件及源码。

不要擅自清理用户工作区或丢弃临时文件。保存任务书、有效结果和最小复现。能创建隔离工作目录／分支则创建，不能覆盖旧结果目录。

### A2. 本轮关键源码定位

以N222快照为参考，具体行号随修改变化：

- `hostApplyDoc222AxisymmetricFsDeform`及变形后初始化块（约5090–5140）；
- 全部`hostZeroVecFieldFrench(...,"ForcePrior")`调用；
- 诊断snapshot/restore是否保存各Previous字段；
- `AcousticStep1stHalf::InitializeKernel`及inner/contact的压力约定；
- `FreeSurface`密度正则、位置增量、梯度修正矩阵更新；
- 新旧顶面接触及黏性接触的body/relation配置。

A交付：`BASELINE_AND_CALL_GRAPH.md`，明确哪些运行用到了哪些初始化分支。

## 7. B：统一重置与恢复，先证明重力没有消失

### B1. 修前复现必须短、小、直接

沿用相同粒子状态：静置后、变形／重置前、清零后、gravity.exec后、首个声学步前，输出每个阶段：

- ForcePrior；
- GravityForceCK与PreviousGravityForceCK；
- Boussinesq当前／Previous；
- 黏性当前／Previous；
- tip当前／Previous；
- 其他实际注册的外力通道。

至少记录质量加权均值、RMS、最大残差，并保留固定OriginalID粒子的原始数值。只有`GravityForceCK=mg`不够，必须证实mg已进入实际使用的共享ForcePrior。

全重力、无其他外力时应有`ForcePrior/m=(0,0,-g)`；若有黏性残余，用分项相加核对，不能要求总力严格等于mg。

本复现无需把旧爆炸实验重新跑满1s。确认首步力丢失即可停止旧路径。

### B2. 只允许两种明确操作

**操作I：保持生命周期。** 普通位置／温度／邻域更新不要清ForcePrior；各力通道依标准增量规则更新。如果位置或速度改变，重新计算受影响的通道即可。

**操作II：完整重新初始化。** 确需重置时，将共享总力及所有已注册通道的current/previous作为一个事务处理，再重新执行全部当前启用的力通道。

```text
确认外力通道清单完整
清零共享ForcePrior
清零每个通道current与previous
按约定重新计算重力、黏性、浮力及其他启用通道
验证ForcePrior = Σcurrent_force
```

禁止只清总力、只清current、仅补一次重力、或用不完整的硬编码三个字段列表。若框架无通道注册枚举，用一个集中维护的明确清单，并给出缺字段即报错的校验。

重新初始化helper只用于阶段切换，不能每声学步清空重算来绕过框架设计。不要在保留旧Previous时把总力再叠加一遍。

### B3. 同时审计其他暂态字段

位置变形、速度清零和时间归零后，核对：DensityChangeRate、位置增量dpos、声学Force、黏性历史、压力、粒子体积、梯度修正矩阵、时间触发器及host/device缓存。

不应盲目全清所有字段；逐项写出“保留／重算／清零”的理由。在第一步声学积分之前必须完成必要的邻域与修正矩阵更新，不能等到第一个平流步之后才更新。

诊断若执行实际声学半步，snapshot/restore必须覆盖该半步改变的完整状态和所有力历史；只比较恢复前后U_rms相等不够。优先使用无副作用统计。

### B4. 验收检查

| 检查 | 预期结果 |
|---|---|
| 恒定重力重复执行 | 不增倍，不消失 |
| 完整重置→重新启用 | 总力正确恢复 |
| 重置→关闭一个通道 | 无旧贡献残留 |
| 相同状态重复重置初始化 | 最终力一致 |
| 排序前后按OriginalID | current、previous与总力对应一致 |
| 诊断前后完整状态 | 无改变或精确恢复 |
| 实际常量受力积分 | 固定物理时间内速度增量符合aΔt |

误差同时用绝对值与相对值；小规模双精度目标1e-10量级、单精度目标1e-5量级，具体容差按量纲／规模说明，禁止对零量只算相对误差。

B通过后，全局搜索剩余ForcePrior清零点，逐个标记为已处理／仅旧复现／不相关。不得只修变形分支而遗留诊断分支。

## 8. C：旧可变形测试只做纠错，不再成为验收主线

### C1. 修改结论标签

保留DOC222原结果，新增勘误，不覆盖原始log：

- 含总力清零不一致的D2总压力结果：`INVALID_FOR_ROUTE_REJECTION`。
- 列拉伸＋密度不一致的恢复实验：`INITIALIZATION_CONTAMINATED`。
- “p′结构不可能”改为当前可变形自由面实现尚未验证。
- D1未经过变形重置的约0.02m/s静水噪声仍保留为真实欠账。

### C2. 有限复核范围

最多重跑：修后一个总压力平液面静水，以及一个无热变形的早期短窗，用于确认B修复前后差异。变形初始化仍不干净时只报告诊断，不授予恢复通过。

不得继续扫local_eta bins、p_h形状、dens_mode或顶端阻尼。即使可变形恢复仍失败，也转入下述固定顶面验证；不要自动阻断全部主线。

如需将来正式验证自由面恢复，应重新生成一致变形域粒子、设置匹配密度／体积／压力并核对总机械能；本轮不展开这项完整开发。

## 9. D：论文平坦滑移顶面的确定实现要求

### D1. 几何与物理定义

玻璃R=0.25m、H=0.185m，竖直z。流体域固定；侧、底无滑移；顶面位于名义高度z_top，不穿透且切向无黏性牵引。平顶下名义流体体积约0.03632m³，记录实际离散体积偏差。

增加顶面边界的接触支持可以使用现有固体边界粒子，不要求新增多层动态ghost粒子系统。边界算法内部虚拟／镜像状态不等于需要用户手工布置ghost粒子。

顶面壁粒子的中心平面与实际物理边界位置区别必须说明，按现有粒子壁面惯例设置。液体和顶面不得重叠。检查顶面／侧壁交线是否重复计算支持、压力和黏性。

### D2. 首选压力路径锁定

首选复用已验证方腔的**固定域扰动压力／约化重力处理**，明确其密度正则、压力／连续性和壁面构造；仅将顶面改成滑移。

完整重力被静水分解消去时，只施加Boussinesq剩余体力。GravityForceCK关闭，实际浮力以`b=-β(T-T_ref)g_vec`定义并说明符号；温度高于参考时向上。

固定不穿透顶面不要求变形自由面的`p′=ρgη`条件。压力响应应由连续性和法向动量／壁面状态决定。**不要在顶面同时硬设p_total=0。** 若闭域使用压力基准或偏置，记录其对EOS与密度的影响；不可逐步减平均压力而不更新密度。

不能简单把当前FreeSurface密度路径原封不动罩在封闭顶面上。采用方腔已验证的固定壁支持与密度处理；禁止全域冻结ρ来压住噪声。

如果已验证方腔实际不是上述压力路径，应先写出真实方程，选用已验证的那条作为替代，并记录一次选择；不要凭名称复制开关。禁止同时开发多套压力主线。

若首选路径在本轮的独立边界测试失败，经定位仍不能修复，可仅增加一个总压力＋完整重力的受控对照。严禁一边保留全重力、一边再次施加已包含完整重力的浮力项。

### D3. 法向不穿透与切向滑移

通过边界Riemann／镜像状态实现。固定顶面n=(0,0,1)，其概念速度关系：

\[
u_{ghost}=u_i-2[(u_i-u_w)\cdot n]n,\qquad u_w=0.
\]

即反射法向分量、保留切向分量。实际关系若使用其他符号／局部坐标，必须推导相同物理效果。

两个声学半步和连续性都必须相容。禁止只在一半步设置滑移，另一半步仍用无滑移速度镜像。

压力壁面状态使用当前路径的实际外加速度。p′路径不补回已消去的全重力，但应核对Boussinesq法向体力的影响。不得为了安静直接关闭顶面连续性或压力跳跃。

不要仅把每对接触压力力强行投影到n就认为滑移已正确；优先复用经过验证的边界构造，检查动量与连续性兼容。

黏性接触只保留侧壁／底面无滑移；顶面不得施加无滑移切向阻力。实现可用分离contact或明确边界类型标记，但顶面仍参与必要的密度支持与压力作用。

禁止每步覆盖流体Uz或z坐标。出现法向运动要在边界状态和离散中解决。

### D4. 热边界独立设置

固定顶面是流动力学边界，不是冷却金属盖。顶面对流／辐射作用于原玻璃表面；不自动新增顶部固体导热、冷水温度Dirichlet或壁体热容量。

侧、底换热沿用既有文献参数；这一阶段不改换热系数。

## 10. E：三个小测试验证平顶，不直接上热核

### E0. 公共配置

主线诊断干预全部关闭。使用同一平顶实现、相同μ、c0、时间步和TVF策略。以现有方腔已验证策略为起点，记录TVF是否开启；不得一边改边界一边试多个TVF开关。

### E1. 等温静止与零浮力

β=0、Q=0、不传热。0.3s预检，再1s。记录U_rms、U_max、密度偏差、力账本、壁面法向速度、粒子位置与控制体收支。

工程目标：U_rms约1e-4m/s或更低，且无持续增长。该值是本项目的候选门槛，不是通用SPH标准。若最终热驱动较弱，还必须满足背景噪声小于驱动信号的5%。

### E2. 顶面切向无阻力的局部算子检查

构造顶面附近均匀切向速度、零法向速度的固定状态；只统计顶面接触贡献，侧／底贡献分开。预期顶面黏性切向牵引接近数值零，压力法向行为不因切向速度凭空变化。

这是无副作用或snapshot完整恢复的算子检查，不要求整个有无滑移底壁的容器保持均匀运动。不得因底壁阻力造成动能衰减而误判顶面滑移。

### E3. 法向扰动短程检查

给定小幅、平滑且初始质量加权平均速度为零的法向扰动，幅度取远低于c0并记录Ma；仅在初始化构造一次，不逐步减平均速度。

通过边界实际推进检查顶面法向响应、越界、密度与压力。无外部持续驱动时，不应发生长时能量自激或集体向上平移。

固定域中局部法向速度不可能在离散壁面附近处处精确为零，故同时报告：近顶层u_n RMS、穿过物理顶面的粒子数／质量、最大超越距离及核支持层厚度。不得把壁粒子中心位置误当物理顶面。

E1/E3关键结果各做一次dt/2复核。出现差异时先查声学／平流／邻域更新时间，不增加新的壁面模式。

如现有顶面实现无法通过E2或E3，停在该测试定位，不能直接开热核并以环流好看宣布修好了。

## 11. F：两种冻结温场，分清求解器与论文目标

### F1. 回归温场：沿用中心热核

同DOC221，ΔT=200K，β=0与β=1e-5各一组。只检查持续浮力和回流闭合：中心热核上升是预期，**不与法国Fig7的轴心下降作通过／失败对比**。

先0.3s，再2s；正常后到5s。记录绝对速度、β差场、控制体收支。无需出现双涡才能通过这个诊断。

### F2. 接近论文的诊断温场

先查仓库已有Fig8数字化数据、坐标映射及之前的温场构造，避免重新做一套相互矛盾的数据。列出来源文件和不确定性。

若已有可靠数字化温场：构造轴对称T(r,z)，保存原图坐标、插值、有效范围和误差；不能把稀疏点插值场当作精确原文数据。

若缺乏可靠数据：优先使用现有热模块在固定域上求得的温场；使用文献空间Q或已经验证的现有Q及原热边界。可先求无流动导热场作为流体初始条件，但明确它不是最终耦合温场。

只在以上两者都不可用时，构造一个**明确标注的定性温场**：上部内部环带偏热、轴心／外壁／表面相对冷。构造公式、温度幅度和宽度一次选定并记录，不为出双涡扫描。它只能证明给定驱动下流体能形成相应趋势，不能授予文献复现通过。

为避免自行选择参数，该最后备选固定为：`T(r,z)=1473+200*exp(-((r-0.15)/0.05)^2-((z-z_bottom-0.13)/0.035)^2)`，长度单位m，温度K，T_ref=1473K。只适用于本计划R=0.25m、H=0.185m的诊断；所有数字都是本轮选定的测试值，**不是论文提取值**。不因没出现双涡而修改这些数字；该组仅检查热环带下的返回路径。

如冻结温场来自文献Fig8，温度恢复使用绝对K和T/Tmax定义；不得随意做min-max归一化或把摄氏温度直接乘比值。Tmax来源缺失则保留温度尺度不确定性。

冻结温场是固定在空间的T(r,z)，粒子移动后按空间重新取值；不要混成粒子携带常温。每次更新后按同一时序更新浮力。

### F3. 物性锁定

流体验证先固定μ=4Pa·s、ρ0=2750kg/m³、β=1e-5K^-1，或使用当前已确认的同一常物性底座，并说明任何差异。不要在回流验证时同时启用未审定的μ(T)。

对应ν约1.45e-3m²/s、H²/ν约23.5s，2–5s只能是早期响应。**不能因为5s尚未双涡就宣布算法失败。**

若静水、边界、收支均通过，且动能稳定、回流正在建立，允许同一配置续算到10s，必要时20s。只延长这一组，不开参数扫描；仍不宣布热场已经稳态。

## 12. G：Python后处理与独立收支，统一通过标准

### G1. 必须输出的图

1. 温度T、u_r、u_z、速度矢量；所有剖面坐标、单位一致。
2. β=0与βon绝对速度图，另画差场；不能只展示差场。
3. 三个以上时刻或时间窗，注明瞬时／平均及窗口样本量。
4. 流函数等值线、闭合轮廓、有效覆盖mask和重建误差。
5. 轴心、内部环带、外缘速度与动能随时间曲线。
6. 三个固定高度的向上／向下／净通量与控制体储量变化。

### G2. 流函数重建

可采用轴对称约定：

\[
u_r=-\frac1r\partial_z\psi,\qquad u_z=\frac1r\partial_r\psi.
\]

处理轴线极限；使用双路径积分或最小二乘并回算速度误差。路径不一致明显时，说明原速度尚不满足相应连续性，不能把投影后的流函数当独立通过证据。

若继续使用表面高度／模态投影统计，必须同时拟合平均高度与形状幅值，避免用离散非正交基把整体平移误当模态变化；报告分箱面积权重、空洞和离群点处理。

检查至少两个方位剖面与方位平均的差异；若明显非轴对称，报告三维性，不用平均消掉异常。

### G3. 固定顶面收支

对近似不可压缩固定域，水平截面的净通量应小。定义参考：

\[
e_Q=\frac{|Q_{up}+Q_{down}|}{Q_{up}+|Q_{down}|+Q_{floor}}.
\]

只在正负流量已明显高于静止噪声时使用该比值；Q_floor按β=0背景确定并预先记录。

工程诊断目标e_Q≤5%，5%–10%标注待改进，>10%不能直接宣布收支通过。这些是项目暂定门槛，不是文献标准；若插值测量误差本身大，应先改测量或与原始3D控制体交叉检查，而非调阈值。

固定域WCSPH仍允许小幅压缩储量变化，因此也比较上方控制体质量与密度储量，不能只凭一个净通量比值判断连续性。

近顶法向RMS目标小于代表性环流速度的5%，越界粒子／最大越界距离完整报告。没有颗粒穿墙不等于速度连续性已经合格；反之一个插值矢量朝外也不证明粒子实际穿墙。

### G4. 回流与文献对比

诊断热核：中心上升＋返回路径＋控制体闭合即可，不强求两个涡。

论文近似温场：检查轴心下降、内部上升带、最外缘下降以及两个独立闭合涡。主次涡不预设为简单上下叠放；依据原Fig7的r-z位置辨识。

不能以Ur符号翻转、流函数正负或circ_proxy非零单独授予双涡通过。

若冻结温场来自Fig8，同一数据不能同时被称为独立温场预测验证。此阶段验证的是给定温场下的流动响应，真正热场预测留给H。

## 13. H：仅在回流可信后接空间热源

沿用现有已验证EM空间Q、归一和映射，不改EM求解器。先冻结Q，演化温度与流动；后续再恢复现有σ(T)更新策略。

首跑0.3→1s检查温度排序、μ更新时序、能量H0和边界损失。之后仅在热场与流动正常时延长。

自然对流文献功率使用50kW；60kW是不同客户工况，分目录。不要把历史均匀55kW混入论文对照。若当前Q未满足50kW配置，先记录并修正输入，不能只在后处理改标签。

流体顶面热损失用玻璃对流／辐射，固定边界粒子不得新增导热散热通道。温度单位K；辐射四次方必须用K。

热初猜可以来自已有热计算或文献近似，以缩短热启动；必须标明来源。H²/α约2.7e4s，不能用几秒温度演化无变化或未达稳态来否定流体模块，也不能把热初猜直接当成预测稳态。

甲方允许未知物性／换热参数校准，但本轮只有在B–G可信后才讨论。校准应单列参数、依据、范围和校准目标；不允许改数值错误来拟合图形。

## 14. 执行矩阵与停止条件

| 阶段 | 必须运行 | 是否阻塞平顶主线 |
|---|---|---|
| B | 修前重力丢失短复现、统一重置、实际积分、排序与日志检查 | 是 |
| C | 旧结果勘误、最多两组有限复核 | 不要求可变形恢复成功 |
| D | 平顶单一配置及完整压力／速度／热边界说明 | 是 |
| E | 静水、切向算子、法向短程、关键dt/2 | 是 |
| F1 | 冻结中心热核β对照 | 是；验证基本回流 |
| F2/G | 论文近似温场、Python图、独立收支 | 文献趋势主线 |
| H | 空间Q温度演化短回归 | F/G可信才启动 |

按执行依赖推进，不把本轮拆成若干“下次再做”的文档。遇到已定位的小缺陷直接修复并回归。

如静水失败，先输出区域分项压力／外力和壁面状态；如滑移失败，只修top边界；如基本回流成功但论文趋势失败，先比较温场，不重写浮力。

连续两个修复候选未改变其针对的独立测试，停止该假设，提交最小复现和分项证据。禁止第三次继续叠开关。

出现NaN、力账本失配、显著越界、密度失真或无驱动动能持续增长时自动停该组，保存第一次异常前后状态；禁止把明显爆炸跑满预设时长。

速度／密度／边界误差阈值在运行前写入配置，结果出来以后不得根据图形修改。

## 15. 验收状态与对外交付边界

| 标签 | 含义 |
|---|---|
| FORCE_RESET_VERIFIED | 所有阶段切换的外力状态一致，实际积分通过 |
| PAPER_TOP_BOUNDARY_VERIFIED | 平顶的滑移、不穿透、压力／连续性通过 |
| DIAGNOSTIC_RETURN_FLOW_VERIFIED | 冻结诊断温场中真实回流与收支通过 |
| PAPER_TREND_REPRODUCED | 论文近似条件下主次环流位置与旋向有明确对照 |
| THERMAL_COUPLING_SHORT_PASS | 空间Q下短程温度／能量／流动回归通过 |
| FRENCH_QUANTITATIVE_VALIDATED | 需要进一步物性、温场、速度及分辨率／时间步证据；不默认授予 |

验收图必须标明模型：固定平顶／可变形自由面、冻结T／演化T、固定Q／更新EM、已校准参数、实际物理时间。

“指定温场下有双涡”和“完整EM—热—流自洽预测双涡”是两个交付层次，不能混写。应向甲方展示已完成层次及未完成项，不把诊断功能包装成最终模型。

若本轮只到可信大回流，也应输出完整证据和差距，不宣布文献双涡完成。

## 16. 结果目录、源码与最终报告

```text
doc223_runs/
  A_baseline/
  B_force_reset/
  C_old_test_errata/
  D_paper_top_model/
  E_boundary_tests/
  F_frozen_temperature/
  G_python_and_budgets/
  H_spatial_q_thermal/        满足前提才生成
  RESULTS_TABLE.md
  ROUND223_REVIEW.md
  NEXT_DECISION.md
```

每个运行目录独立保存：命令、最终JSON配置、退出码、log、原始CSV、选定粒子快照和Python图。不要多个并发进程共用一个output目录。

打包所有修改过的共享源文件、入口、构建配置、patch、reload身份哈希、图与原始测量，不再只有主入口CPP。二进制可仅给哈希，实际编译路径和构建时间必须有。

最终报告必须回答：

1. 重置是否确实丢失重力？修前／修后每单位质量贡献各是多少？
2. 还发现了哪些残留current/previous或暂态状态错误？
3. 固定平顶采用了什么真实压力、连续性、法向与切向状态？
4. 原有硬顶失败与本次实现具体差在哪里，而非只改开关名？
5. β=0背景、βon响应、时间步复核是否支持真实回流？
6. 论文近似温场与中心热核是否被明确区分？
7. 是否出现主次涡，位置／旋向／收支依据是什么？
8. 当前可以交付哪一级证据，距离完整文献复现还缺哪一项？

## 17. 给Cursor的最终硬约束

**先修状态错误，再用论文对应的固定顶面模型验证回流；可变形自由面全面开发不再占据验收关键路径。**

本轮允许必要的实现判断，但不允许自选新算法、新参数扫描或新研究主题。若配置／资料缺失，先检索项目已有证据并记录缺口；不能编造已数字化数据、精确物性曲线或文献边界细节。

每个PASS对应一个定义、原始数据和反例检查。最终目标是可解释、可重复的物理结果，不是通过开关把图形调成环流。


<!-- ===== END 223_OPHELIE_FORCE_RESET_REPAIR_AND_PAPER_TOP_BOUNDARY_CURSOR_PLAN_2026-11-17.md ===== -->

---

<!-- ===== BEGIN 224_OPHELIE_PAPER_TREND_THERMAL_ALIGNMENT_AND_ACCEPTANCE_CURSOR_PLAN_2026-11-17.md ===== -->

<a id="224-ophelie-paper-trend-thermal-alignment-and-acceptance-cursor-plan-2026-11-17-md"></a>

# 224_OPHELIE_PAPER_TREND_THERMAL_ALIGNMENT_AND_ACCEPTANCE_CURSOR_PLAN_2026-11-17.md

# 第224轮开发任务书：锁定有效求解配置，完成后期双环流与文献热场对齐

日期：2026-11-17  
项目：冷坩埚招投标／法国OPHELIE自然对流  
执行依据：DOC223_ROUND_PACK_2026-11-17及其原始数据、代码审查  
执行对象：Cursor。本轮以本文为主，不继续叠加221—223中已经弃用的试验。

## 0. 本轮目标、顺序和停止边界

本轮不是继续证明“有没有浮力”，也不更换求解器。当前已经有无人工速度干预的回流，并在指定环带温场下出现双回流形态。下一步必须回答：

1. 这个结构是否在完整后期时间窗内持续存在，顶部和质量收支是否可信？
2. 给定论文近似温场时，主次环流的位置和旋向是否与Fig7相符？
3. 使用50kW现有空间热源，温场与流动能否向论文条件发展？

执行顺序：A 固定底座 → B 修正测量与判定 → C 延长原F2并检查数值可信度 → D 文献参考资料与固定温场比较 → E 空间Q热场推进 → F 有条件恢复电磁反馈 → G 验收结果表。

B与D资料整理可独立开展；C、D中的流体运行使用同一锁定模型。E可在短程基线正确且数值无红线时推进，不必等待缺失的图像资料。F必须在E解释清楚后才启动。

**禁止将缺少文献数字化资料解释为“流体失败”；也禁止将指定温场双涡解释为完整电磁—热—流复现通过。**

本轮必须交付可重复结果及明确差距，不接受仅生成另一份计划、等待下一轮继续。

## 1. 基线事实与当前完成层次

### 1.1 已有证据，不重做大范围验证

| 项目 | DOC223证据 | 本轮处理 |
|---|---|---|
| 外力重置 | 清总力但不清Previous会丢失重力；完整重置后恢复约−g | 保留修复，短程回归即可 |
| 平顶静水 | U_rms约3.4e-7m/s；dt/2约6.4e-6m/s | 有效背景参考，不把两者当严格时间收敛 |
| 顶面滑移 | 顶部切向黏性贡献0，侧底仍有贡献 | 保留实现 |
| 法向扰动 | 动能衰减至约1.1%，短程无越界 | 保留，长期边界另查 |
| 中心热核 | 轴升缘降至5s | 基本回流通过，非论文目标 |
| 环带温场F2 | 轴降—环升—缘降，原始速度场可绘出两个回流单元 | 指定温场的双胞证据，尚非论文验证 |
| H1空间Q | 50kW，温度演化至1s，能量账本稳定 | 短程耦合通过，尚非稳态 |

### 1.2 本轮必须补齐的具体欠账

- 现F2图平均从约0.05s启动直到10s，不能作为后期稳态场。
- 某些时间窗只有一个样本、实际覆盖不足1ms却标记ok=1，必须修。
- F1/F2长期最高粒子可超名义顶面约0.1—0.15dp，未量化质量占比和持续性。
- 原截面净流量不平衡约3%、5.6%、6.4%；只有网格估计，独立控制体证据不足。
- 原ψ重建采用单径向积分，RMSE约1.5e-4／2.5e-4m/s，不能直接证明原速度无散。
- 未使用可信Fig7/8参考做位置与温度比较。
- H1损失约156kW、输入50kW，整体缓冷；其初态与热边界解释尚需核对。
- H是冻结Q，命令包含no-self-induction，未恢复σ(T)周期刷新；不能称完整论文电磁模型。

## 2. 冻结求解配置，防止重复试错

### 2.1 本轮默认保持

| 配置 | 锁定值／约定 |
|---|---|
| 玻璃几何 | R=0.25m，H=0.185m，z竖直 |
| 基准粒距 | dp=0.015m |
| 流动边界 | paper_flat_slip_top：固定平顶滑移、不穿透；侧底无滑移 |
| 压力 | DOC223已通过的p′路径，完整重力不重复加入 |
| 密度 | 静置阶段原处理；正式推进不周期性dens重初始化，密度仍由连续方程更新 |
| 核梯度修正 | OFF，沿用已通过配置 |
| TVF | OFF，沿用DOC223配置 |
| 诊断人工干预 | tip、zero-net-Uz、Uz kill、z-clamp、ContSilent全部OFF |
| 基准ρ0、μ、β | 2750kg/m³、4Pa·s、1e-5K^-1；如实际基线不同须显式核对 |
| 基准c0、时间步 | 沿用DOC223最终生效值；另做一个关键dt/2对照 |
| 初始阶段 | 保留外力完整重置、温度字段排序与H0修复 |

温度依赖物性应按运行分层启用。C诊断工况保持原来的物性实现；先记录DOC223实际是否恒μ，禁止仅凭日志mu=4就认定没有μ(T)。D/E中启用哪些曲线必须在配置表明确列出。

### 2.2 只有以下情况可改底座

B/C证实某具体实现错误，例如窗口权重、体积权重、字段排序或边界标记错误。先提交修前证据、单一修改和对应回归。不能因为结果不像论文就开密度、KC或压力模式扫描。

修改共享代码仍限定OPHELIE策略，保留其他case默认行为；独立A-phi／Laplace分支不改。

### 2.3 基准、重启和输出身份

保存SHA、dirty diff、实际编译源、Real精度、设备、完整命令、最终生效参数、reload与Q场哈希。

若已有可靠checkpoint，可续算F2；必须保存全部current/previous力、密度／压力、温度、位置增量、时间与累积统计状态。先用短区间证明续算与不中断运行一致。

若checkpoint不完整，直接按相同初态重跑到20s，不能只导入位置和速度伪称从10s续算。重启不得重新静置、清零速度或重新注入初始温场而不说明。

## 3. 所有标准的性质和判定规则

下列数值门槛是**本项目预先制定的工程判定标准**，不是法国论文公布的误差标准，也不自动构成发表级验证。

每项使用：`PASS / WARN / FAIL / NOT_EVALUATED`，附reason和数据路径。

- PASS：证据覆盖该指标与对应工况。
- WARN：结果有趋势但达不到严格门槛，不得被顶层程序自动提升为PASS。
- FAIL：可信数据明确违反门槛，或出现实现／运行红线。
- NOT_EVALUATED：缺参考、缺样本、未运行，不能写成零误差或默认PASS。

总标签由列出的必需项推导。禁止直接读取旧`passed=1`、`smoke_passed=1`或退出码0授予物理通过。

缺少参考与数值失败应采用不同reason。不得为了结束任务修改阈值或裁掉不利结果。

## 4. A：一次性基线回归

只做必要短回归：零浮力平顶静水1s，检查新测量代码不改变原物理结果；已有外力重置／tip自检若未改相关实现，不重复长跑。

静水要求无定向持续增长、U_rms≤1e-4m/s；后续有浮力信号应至少达到对应背景噪声的20倍。该要求覆盖本轮工程比较，不是通用公式精度证明。

自检：相同物理状态调用测量前后，Velocity、Density、Position、ForcePrior及所有历史通道不改变。需要执行动力学的诊断，必须完整snapshot/restore并逐字段检查。

## 5. B：先修后处理，再用它作判定

### B1. 完整时间窗

统一窗口：0—5、5—10、10—15、15—20s；延长后使用20—25、25—30s等同长窗口。

每窗记录requested_duration、实际覆盖时长、coverage_ratio、样本数、最大采样间隔、是否跨重启。权重按实际物理时间，不按样本个数。

有效窗要求：覆盖≥99%，至少50个时间采样，最大间隔≤0.2s；不得用插值向未计算区间外推补齐。

只覆盖约0.0007s的10—15s窗口，必须NOT_EVALUATED。终点刚越过10s的一个样本不能证明10—15s存在。

瞬时场、单窗平均场、启动全过程平均场分目录和标题。主报告以完整后期窗口为主。

### B2. 采样区域和权重固定

保留全域网格及有效mask；额外区域统计统一采用r/R与z/H。

- 轴区：r<0.08m；环区：0.12≤r≤0.18m；缘区：r>0.18m。
- 同时提供全高统计和排除上下各一粒距的内部统计，不能只选有利版本。
- 区域均值使用对应体积／环带面积权重；不能把各r-bin简单算术平均冒充体积平均。
- 不隐藏空bin；流函数积分不能把未知bin默认为真实零速度。

截面取z/H=0.25、0.50、0.75；如取最近网格层，报告实际高度与偏差。

### B3. 流函数与闭合轮廓

约定`u_r=-(1/r)∂ψ/∂z`、`u_z=(1/r)∂ψ/∂r`。处理轴线极限与积分常数。

改进原单向累积方法：采用两条独立积分路径作一致性检查，或在有效网格上最小二乘重建；必须回算u_r/u_z并报告误差。

对流场内部有效区定义`Eψ = ||u_reconstructed-u_original||_V / max(||u_original||_V, U_noise*sqrt(V))`。PASS≤10%，WARN为10%—20%，>20%不能单靠ψ授予闭合双涡；不自动判流体失败，先核对重建／采样误差。

同时绘制原始速度插值的streamlines与ψ等值线，两者的旋向与涡心应一致。增加采样bin数只能改善重建，不能称粒子分辨率提高。

### B4. 质量、体积与截面收支

必须分别输出：总质量、Σm/ρ、固定域体积、控制体质量及几何占据偏差；各有不同含义。

对水平截面报告Q_up、Q_down、Q_net，先输出不带大floor的原始量，再判定：

`e_Q=|Q_net|/(Q_up+|Q_down|)`。

只有总交换量高于β=0噪声通量20倍时评估相对误差，否则NOT_EVALUATED。不要用大Q_floor人为压低e_Q。

工程标准：≤5% PASS；5%—10% WARN；>10% FAIL于该测量项。至少一条原始3D粒子控制体路径独立检查，不能让同一r-z插值既算流量又证明自身正确。

WCSPH固定域仍有小幅压缩储量变化。对上方控制体同时比较质量流率与dM/dt；实际密度不恒定时，质量收支比Q净值更直接。原始Q不平衡而质量储量能解释时，报告压缩效应，不能冒充严格不可压缩。

无丢粒子／相变时总质量相对变化目标≤1e-6；检查Σm/ρ相对静置后漂移：≤1% PASS、1%—2% WARN、>2%需停止延长并定位。重建体积误差单列，不能用Σm/ρ代替几何占据体积。

## 6. C：原环带工况延长，先回答是否稳定且边界可信

### C1. 温场与配置不变

沿用原F2：

`T=1473+200*exp(-((r-0.15)/0.05)^2-((z-z_bottom-0.13)/0.035)^2)`。

此式为人工诊断温场，非Fig8。保持空间冻结T、Q=0、β=1e-5及其余锁定配置；不移动热点或调整温差以维持双涡。

先计算到20s。无红线、仍在稳定发展但未满足后期不变性时，可到30s；本轮该组默认不超过30s。若仍非稳态，标记NONSTEADY并报告趋势，不改物性强行加快。

### C2. 后期不变性的工程标准

比较相邻完整5s窗口：

- U_rms、动能及主次环流强度变化≤5%；
- 归一化速度场加权L2差≤10%；
- 两个涡心位置变化≤一个dp；
- 主次旋向和分隔结构不变；
- 顶部及密度误差无持续累积。

所有条件满足，可授予`LATE_WINDOW_STABLE`，只是该冻结T工况的后期稳定证据；不是热稳态或法国文献验证。任一不满足，保持NONSTEADY，不能只挑一条速度曲线。

### C3. 顶面检查必须定量

定义物理顶面z_lid，区分壁粒子中心和接触边界平面。

记录每个采样时刻：z_max、超顶粒子数、超顶质量比例f_over、超越深度分布、近顶u_n RMS；近顶层固定为顶面下2dp，并单列已超顶粒子。

本轮暂定：

| 指标 | PASS | WARN | FAIL／停止延长 |
|---|---|---|---|
| 最大超顶距离 | ≤0.1dp | 0.1—0.25dp | >0.25dp |
| 超顶质量比例 | ≤0.1% | 0.1%—1% | >1% |
| 近顶法向RMS／内部代表流速 | ≤5% | 5%—10% | >10% |
| 发展趋势 | 无持续累积 | 证据不足 | 后期持续恶化 |

这些容差承认粒子边界有限分辨率，不能描述为“精确零穿透”。流速未超过噪声20倍时，法向比值不评估，改报绝对值。

WARN不得删掉粒子或覆盖z修复；如WARN稳定、其余指标可靠，可继续D/E并如实标注边界精度限制。FAIL先定位物理平面、壁面支持、法向及接触状态，禁止直接恢复旧hard-lid包。

### C4. 一个针对性的时间步对照

在同一配置下dt/2计算到与基准相同物理时刻，至少覆盖一对完整后期窗口。不要只比较静水dt/2。

代表速度／环流强度差≤10%，速度场L2差≤15%，涡心位置差≤一个dp且旋向相同，作为工程时间步鲁棒性PASS。超过则先查时间推进与边界；不通过修改μ或β解决。

### C5. 双胞判定

同时具备：原速度场两侧反向回流路径、各自内部涡心、独立闭合轮廓、共同上升带、不同方向的环流；并满足B的有效场覆盖。

涡不能只有一个孤立bin极值：每个涡的核心邻域应至少覆盖约2×2个有效采样单元，且跨时间窗可追踪。若窄外围涡只有一列bin，标记未充分分辨。

此外须报告涡的实际宽度／dp；小于约2dp的结构不得授予充分分辨，不能靠增加采样bin数量通过此项。

轴降—环升—缘降仅为必要的辅助证据，不单独授予双胞；流函数颜色正负也不单独授予。

## 7. D：Fig7/8参考提取和固定温场验证

### D1. 资料复用优先，不能声称已拿到不存在的数据

先在项目文档／reference目录查找已有Fig7/8数字化、原PDF、对应博士论文补充。列出命中的文件和版本，检查是否真包含坐标、图例和数据，而不只是旧任务书。

主参考为《Numerical modelling of natural convection in molten glass heated by induction》，Chemical Engineering and Processing 47 (2008) 449–455，DOI：10.1016/j.cep.2007.01.013。博士论文仅补期刊缺项，不能静默把不同工况混在一张比较表里。

当前DOC223包未提供可信Fig7/8参考数据，所以本轮不得默认已完成。若执行环境缺原PDF，写REFERENCE_MISSING并继续C/E；不要把人工环带公式改名Fig8。

### D2. 提取物理信息，不对每条流线伪造精确值

统一坐标r/R、z/H。记录图号、页码、原图像素尺寸、裁切框、轴向方向、镜像关系和图例。

Fig7：旋向、涡心区域、内部上升带、外围下降带、分隔线和主次相对大小；只有图例明确给出数值／单位时才提取ψ或速度大小。流线间距不自动等于速度大小。

Fig8：可识别等温线／色阶及其值、热点区域、轴心／侧壁／顶部相对温度。T/Tmax必须保留绝对K定义；不得做min-max替换。

若Tmax缺失，输出归一化温场及温标不确定性，不用1473K表格物性参考温度冒充图中Tmax。缺失绝对温标可做单列假设的诊断，不能授予绝对温度定量通过。

所有手工／像素读数带误差范围，重复选点或以线宽估计不确定性。无数值图例时只授予拓扑参考。

交付`fig7_features.csv`、`fig8_samples.csv`、有效mask和`REFERENCE_PROVENANCE.md`；无法提取的列保留缺失，不填猜测数字。

### D3. 固定论文近似温场的运行

构造受参考约束的T(r,z)，保存插值方法、数据范围和不确定区域；不加超出图中证据的局部热点来人为制造第二个涡。

固定温场下运行同一流体模型。先0.3→2→5s作检查，再根据C得到的响应时间选择20s，必要时30s；执行完整后期窗口分析。

本阶段输入温场来自Fig8，则温度场本身不是独立预测验证；验证对象是该温场下的流动与Fig7是否相容。

### D4. 文献趋势通过标准

必须全部满足：

1. 轴心下降、内部环带上升、外缘下降的空间关系与原图一致。
2. 两个环流旋向、主次位置和分隔关系与原图一致，不预设上下叠放。
3. 涡心、上升带位置偏差不超过`max(0.10L, 2dp, 文献读图位置不确定度)`，径向L=R、竖向L=H；报告误差预算，不将此粗分辨率标准称精确匹配。
4. 至少两完整后期窗口保持趋势；数值可信度无FAIL，WARN完整披露。
5. 结果来源、温场假设与是否校准清楚。

满足则`PAPER_TREND_FROZEN_T=PASS`。绝对流速和ψ若无可靠文献值，标记NOT_EVALUATED；禁止授予定量通过。

## 8. E：现有50kW空间Q下推进温场，先解释热边界

### E1. 固定一条热源基线

复用H1的空间Q、50kW及映射，先不改σ或电磁求解器；每次空间采样记录实际总功率。禁止换均匀55kW，也不把客户60kW工况混入论文对照。

明确记录：Q来自哪次EM计算、σ值／温度、频率、线圈、电磁自感是否开启、是否按目标功率归一。冻结Q不是实时温度反馈。

### E2. 热边界审计表

每个表面列出面积、温度采样方式、h、环境／冷却温度、辐射率、单位和来源；侧／底／顶部完全分开。

固定平顶只改变流动边界，不能额外成为冷金属导热盖。检查无重复顶面热损、无重复侧壁热交换、无把玻璃体积当表面积、无摄氏温度进入T^4。

与原文逐项比对；原文给范围则保留范围，缺值明确标为假设／可校准。不得将当前默认值自动标记为文献精确值。

H1约156kW损失、50kW输入的缓冷首先用账本解释：`dH/dt=P_in-P_loss`。匹配不等于边界物理已正确，不匹配也不能靠校准h掩盖代码错误。

### E3. 物理启动和热初猜分开

运行E-cold：沿用1473K等温初态，1→5s，正常后最多10s，观察热点、边界冷却、ΔT和速度。它是启动诊断，不要求短时热平衡。

运行E-warm仅在有可信热初猜时启动：使用文献近似温场或现有热模块生成的温场，随流体移动正常演化T；不每步重新指定T。

若通过无流动热计算预热／预冷获得初猜，保存独立热计算时间、边界、能量账本与结果，明确不是耦合流动稳态。可使用现有轴对称或隐式热模块加快获得初猜，但不另开发新求解器。

不得只放大Q、压低cp或每步重设温度来冒充物理加热时间。稳态求解／伪时间若已有成熟实现，必须说明方程和残差，不能直接称真实瞬态。

### E4. 热与耦合标准

| 指标 | 通过标准／解释 |
|---|---|
| 实际P_heat对50kW | 相对误差≤0.1%；超过先查映射／归一 |
| 累积能量闭合 | `|ΔH-E_in+E_out|/max(|E_in|,|E_out|,能量噪声下限)≤1%` |
| 能量噪声下限 | 由β=0短回归／浮点量级确定并固定，不能用总初始焓掩盖增量误差 |
| 无红线 | 无NaN、异常密度、力账本错误或持续越界 |
| 短程耦合 | T、浮力、速度更新与能量账本对应即可；不要求P_in=P_loss |
| 热准稳态候选 | 完整窗口平均`|dH/dt|/P_in≤5%`，且`|P_in-P_loss|/P_in≤5%`，两者一致 |
| 热场稳定 | 相邻完整窗温度场变化RMS≤2%温差尺度，热点位移≤dp；近等温时用绝对误差说明 |

热稳态标签只在候选条件持续两个完整窗口且物性／边界稳定时授予；1s能量残差很小不能替代热稳态。

耦合推进中至少输出T_mean/T_min/T_max、ΔT、轴心／热环带／侧壁／顶面温度、三项表面损失、速度与主次涡指标。先对比温场，再解释流向。

### E5. 校准严格后置且有限

用户允许未知参数校准，但只有边界审计与代码收支正确后才允许。已明确的几何、功率、浮力符号不得改。

每次只选一个缺失／不确定参数，以原文范围或明确物理范围为依据，最多三个预先选定值；记录改变前后的温度与流动。优先选热边界不确定项，而非改β凑流线。

如缺文献范围，先列为待定而不启动任意扫描。校准结果必须标CALIBRATED；用来拟合的图不能再称独立验证。不得通过校准掩盖明显数值不收敛。

## 9. F：电磁反馈为条件任务，不重新开EM算法研究

E已得到可信温场／流动趋势，且σ(T)曲线与控制方式来源明确时，才复用现有周期刷新恢复EM反馈。

固定一个刷新策略：按温度导致的电导率变化控制，候选目标每次刷新间最大相对σ变化≤5%；初次测量变化率后据此确定时间间隔。再用间隔减半做一次对照，不扫描多个周期。

必须明确固定总功率还是固定电流；论文50kW功率对照下，若采用保持50kW，应说明这是功率控制，不是EM空间分布已经独立验证。刷新后不能静默改变驱动含义。

记录每次刷新时间、σ范围、EM残差、Q场变化与总功率、映射误差。自感／趋肤近似单列；本轮不强制重写自感算法，也不将no-self-induction模型称完整文献电磁解。

如果E仍处于初始热发展，F标NOT_RUN，不能为了“全耦合”标签提前增加一个难以归因的变量。

## 10. 有限分辨率复核：只检验一个已经有意义的结果

C或D已得到稳定双胞／可信趋势后，选一组关键工况将dp由0.015m改为0.010m，保持物理参数和边界位置，重新生成匹配粒子布局与必要壁面支持。

不得用同一粗粒子reload改dp标签冒充加密。空间Q采样保持物理场一致；需要时重新映射并核对功率。

比较涡的数量／旋向、涡心／上升带、U_rms、环流强度、顶部偏差和收支。位置差≤0.1L且强度差≤15%可记工程鲁棒；否则报告分辨率敏感，先定位薄外围涡／壁层分辨问题。

两档分辨率只提供敏感性证据，不能声称已完成严格网格无关性或收敛阶验证。不要在热场未形成前启动三四档长时耦合。

## 11. 自动判定与人工图审必须同时存在

至少实现一个独立评价脚本，读取原始CSV／配置，输出每项value、unit、threshold、status、reason、source_path。

用四类小数据验证评价脚本：完整／截断时间窗、零信号、已知闭合流场、故意引入净通量或缺失bin。确保一个样本不能通过完整窗、未知不能自动变零、零分母不产生虚假PASS。

已知流场可用解析轴对称ψ构造速度测试重建与通量，但它只验证后处理，不代替流体物理验证。

图审须展示：原始绝对速度、温度、mask、ψ、流线、后期窗口对比和边界误差。不只交一张光滑流线。

## 12. 最终状态分层，防止验收误述

| 标签 | 必要证据 |
|---|---|
| DIAGNOSTIC_TWO_CELL_OBSERVED | 人工温场下原速度有双胞形态，未必稳定 |
| DIAGNOSTIC_TWO_CELL_STABLE | C的后期不变性＋无数值FAIL＋闭合路径证据 |
| PAPER_TREND_FROZEN_T | D参考有出处、位置旋向符合标准；温度是输入 |
| SPATIAL_Q_COUPLING_VALID | E空间Q下T和流动实际演化、账本闭合、无数值红线 |
| PAPER_TREND_COUPLED | 耦合预测的T与主次环流趋势同时符合原图，披露校准与EM近似 |
| PAPER_QUANTITATIVE_VALIDATED | 可靠数值参考、误差预算、时间步／分辨率证据齐全；本轮不默认授予 |

存在WARN可给出对应观察结果，但不能省略WARN写成“全部通过”。严格标签所需项缺失则NOT_EVALUATED；数值红线则不能授予稳定／预测通过。

本轮成功不要求一次取得所有标签，但至少完成B修复、C后期证据和E热边界／启动解释；D/F受参考和阶段条件影响时，写明阻塞原因和唯一下一动作。

## 13. 运行预算与失败后的固定分支

| 情况 | 必须执行的下一动作 | 禁止动作 |
|---|---|---|
| 原F2仍稳定发展 | 同配置到20s，条件满足才到30s | 改热点／β保双涡 |
| 顶部WARN但不累积 | 保留诊断并推进文献／热场，做关键dt或dp对照 | 清零Uz或裁切粒子 |
| 顶部FAIL | 查物理平面、壁层、法向与边界状态 | 引入tip封口 |
| ψ不一致但原速度有双胞 | 修重建、mask和时间窗，比较独立控制体 | 投影后替换原速度 |
| 论文温场缺绝对温标 | 做归一化参考，绝对项NOT_EVALUATED | 直接取1473作Tmax |
| 空间Q下流弱、ΔT刚发展 | 按实际热时间／初猜解释，继续受控热启动 | 判浮力失效或改μβ |
| 能量守恒但损失大 | 审计温度与热边界来源，区分初态和稳态 | 直接调h使损失等于50kW |
| D趋势过、E趋势不过 | 优先比较Q、T和热边界 | 再换流体求解器 |

已有算法错误的修复最多测试两个针对同一假设的候选；未改变独立证据则停止该假设，输出最小复现，不继续叠开关。

## 14. 交付目录与报告要求

```text
doc224_runs/
  A_baseline/
  B_metrics_and_tests/
  C_f2_late_windows/
  D_fig7_fig8_reference/
  D_frozen_paper_temperature/
  E_thermal_boundary_audit/
  E_spatial_q_cold_start/
  E_spatial_q_warm_start/       有初猜才运行
  F_em_feedback/               满足条件才运行
  R_dt_dp_checks/
  ACCEPTANCE_TABLE.csv
  ACCEPTANCE_TABLE.json
  ROUND224_REVIEW.md
  NEXT_DECISION.md
```

不同进程不得共用output。保存配置、代码身份、原始CSV、必要VTP、图、评价脚本和所有修改过的共享源码／patch。

最终报告开头直接回答：

1. 后期双胞是否稳定，依据哪两个完整窗口？
2. 顶面、密度、质量收支是否有FAIL或WARN？
3. Fig7/8参考究竟拿到了什么、缺什么？
4. 给定论文近似温场后，流动趋势是否匹配？
5. 50kW空间Q预测的温场为何接近或偏离论文？
6. 当前双胞是指定温场响应，还是温度演化预测？
7. 电磁反馈是否实际更新、采用哪些近似？
8. 现在能交付哪个层次，尚缺的最高优先项是什么？

**最终执行原则：保留DOC223的突破，把工作重心转到可信时间窗、文献温场和热源驱动；每个结论有原始证据，每个缺口有固定处理方向，不再凭图形好看判断成功。**


<!-- ===== END 224_OPHELIE_PAPER_TREND_THERMAL_ALIGNMENT_AND_ACCEPTANCE_CURSOR_PLAN_2026-11-17.md ===== -->

---

<!-- ===== BEGIN 225_OPHELIE_ABSOLUTE_TEMPERATURE_SLIP_AUDIT_AND_WARM_RESTART_CURSOR_PLAN_2026-11-18.md ===== -->

<a id="225-ophelie-absolute-temperature-slip-audit-and-warm-restart-cursor-plan-2026-11-18-md"></a>

# 225_OPHELIE_ABSOLUTE_TEMPERATURE_SLIP_AUDIT_AND_WARM_RESTART_CURSOR_PLAN_2026-11-18.md

# 第225轮整改与推进任务书：绝对温度保真、浮力审计、顶面总滑移与暖启动复核

日期：2026-11-18  
项目：冷坩埚招投标／法国OPHELIE自然对流复现  
依据：DOC224_ROUND_PACK_2026-11-18及其源码、日志、输入CSV审查  
执行对象：Cursor。本轮以本文为唯一主任务书。

## 0. 本轮目标与执行顺序

本轮只解决四项，不新增电磁算法，不重新更换流体框架：

1. 保留CSV中的绝对温度，分离热量、物性与浮力参考温度。
2. 修复近零总浮力的错误相对误差判据，恢复可信的逐粒子力审计。
3. 补齐O3顶面的声学＋黏性＋连续性总滑移检查。
4. 核对Fig8重建温场，使用修正后的输入重跑冻结温场及50kW暖启动。

执行顺序：A 固定底座 → B 温度与参考量整改 → C 浮力审计 → D 顶面总滑移 → E Fig8重建核查 → F 修正后的冻结／释放温场 → G 有条件后期延长与交付。

B/C通过后可开展短程F；D与E可独立整理，但长期物理结论必须包含其检查结果。

禁止再用“热损失大”直接决定调h，也不能以“30s没达到稳定”直接否定求解器。先排除输入与边界解释问题。

本轮应完成必要代码、短测试和有判别力的运行，不要只输出另一份工作计划。

## 1. 本轮事实、必须更正的结论

### 1.1 保留的有效进展

- 外力重置丢失重力已定位并修复。
- 人工环带F2在30s仍有轴降—环升—缘降趋势；不能再说当前完全算不出回流。
- O3顶面长时没有测得超顶，值得保留为候选。
- F2在20—25s窗口的主要速度标量对dt/2差异小于0.1%。这只是标量鲁棒，不是完整场与拓扑已验证。
- 网格截面净通量误差改善至约2.6%—3.2%，独立3D收支仍未补齐。
- 暖启动dt/2可运行至20s且能量账本闭合，是短程耦合证据；未热稳态不是耦合全部失效。

### 1.2 必须整改的温度路径

函数`hostApplyDoc224RzTemperatureGrid`在demean时执行：

```cpp
delta_t[i] -= demean_shift;
temp[i] = t_ref + delta_t[i];
```

本轮日志：输入CSV范围750—1425K；demean_shift约−394.272K；温度整体被抬高，运行中最高温度达到约1831K，平均被推向1473K。

因此，“Tmax1500_ASSUMED输入”的实际绝对温度已被改变。这直接影响Robin换热、T^4辐射和温度依赖物性，旧暖启动不能直接用于判断原CSV的热损失。

这项问题不必然解释恒物性冻结温场的流向反转。均匀温度偏移对温度梯度没有直接影响；不能只修demean就宣布D已解决。

### 1.3 必须修正的浮力判定

日志总浮力约1.59e-8N，解析约−3.40e-8N，绝对差约5e-8N，却因分母近零得到rel_err≈1.47、force_ok=0。不能以此判断浮力错误，也不能一直保留一个不解释的WARN。

### 1.4 O3标签范围

O3目前是“声学关闭AcousticWallFreeSlip＋黏性OpenCup”。E2只证明顶部黏性切向贡献为零，不覆盖声学压力切向贡献和连续性镜像状态。

本轮补验证，不预先认定O3错误，不直接恢复旧P2开关包。

## 2. 锁定范围与禁止事项

沿用R=0.25m、H=0.185m、dp=0.015m、z竖直、O3候选顶面、已通过的p′路径、无完整重力重复项、正式推进不做周期性dens重初始化、KC OFF、TVF OFF。

基准ρ0=2750kg/m³、μ=4Pa·s、β=1e-5K^-1；明确实际恒物性／μ(T)配置。本轮因暖启动默认步长已发生失稳，修后暖启动从`dt_scale=0.5`开始，另做一个关键dt_scale=0.25短对照。

人工速度修正全部OFF：tip、zero-net-Uz、Uz kill、z-clamp、ContSilent。不得更改独立A-phi／Laplace分支。

禁止事项：

1. 不为避免总浮力非零而改绝对Temperature。
2. 不把浮力参考温度同时传给辐射、物性或热焓公式。
3. 不为让损失等于50kW而直接调低h或ε。
4. 不修改β、μ、热点位置来强行保留双涡。
5. 不继续扫p_h、dens_mode、顶面阻尼和核修正开关。
6. 不在修正后的热状态尚不明确时开启σ(T)周期EM刷新。
7. 不以程序退出0／旧passed=1代替本轮物理验收。

## 3. A：身份、勘误和配置约束

保存SHA、dirty diff、实际编译源、精度、设备、reload与温场CSV哈希、完整命令及最终生效配置。不要覆盖DOC224原结果。

新增`DOC224_ERRATA.md`：

- 带demean的暖启动标记`ABSOLUTE_T_SHIFTED_LEGACY`。
- O3顶面“不穿透通过”和“总滑移尚待补证”分开。
- C4标量时间步一致性有效，但完整场L2／涡心比较未覆盖。
- 近零总力rel_err失效，撤回对应force_ok结论。
- 缺导热的冷启动保留为配置对照，不与启用导热的物理工况直接比较。

建立`warm_spatial_q`运行预检：freeze_T=0、thermal_diffusion=1、Q_actual≈50kW、顶面无额外冷盖、无legacy温度平移、无人工速度干预。任一不符直接报错，不得运行后才发现漏开导热。

## 4. B：绝对温度与三类参考量分离

### B1. 明确字段语义

| 字段／参数 | 唯一含义 | 允许使用方 |
|---|---|---|
| Temperature | 实际绝对温度K | 传热、辐射、物性、EM、浮力、输出 |
| T_thermal_ref | 热焓／热增量参考 | 既有热账本及其DeltaT |
| T_property_ref | 文献物性参考点，例如1473K | 物性定律定义；不替代实际T |
| T_force_ref | Boussinesq剩余体力参考 | 浮力通道与一致的压力／边界约定 |

允许这些参考值数值相同，但禁止代码依靠它们必须相同。

`OphelieThermalDeltaT`保留现有热模块要求的`Temperature-T_thermal_ref`语义。浮力目前读取Temperature及独立设备标量参考，不应借改热DeltaT实现均值消除。

### B2. 温度加载函数整改

普通路径只执行：

```text
T_abs[i] = CSV.sample(r_i,z_i)
thermal_deltaT[i] = T_abs[i] - T_thermal_ref
```

不加偏移，不把Tmean钉在1473K，不因浮力零均值覆盖Temperature。

冻结温场：粒子移动后按同一个空间CSV取值；允许重采样，不允许根据当前粒子温度均值再平移场。

释放温场：只在初始阶段加载一次，随后由传热与对流演化；不得每次输出或更新邻域时重新加载CSV。

旧`--init-T-rz-demean`不得静默改成新含义。正式运行遇到旧参数应报错并说明迁移；若必须复现旧结果，使用显式legacy模式，并禁止用于验收主线。

### B3. 本轮浮力参考温度的确定选择

现有人工环带F2保持原T_force_ref，避免改变旧对照。

新Fig8输入组采用：**初始CSV加载到粒子后，计算一次质量加权平均温度作为T_force_ref，并在该次运行期间固定。**

\[
T_{force,ref}=\frac{\sum_i m_iT_i(0)}{\sum_i m_i},\quad
F_{b,i}=-m_i\beta[T_i-T_{force,ref}]g_{vec}.
\]

这样初始净浮力接近零而绝对T不变。该值是数值参考选择，不是文献温度测量。

在初始粒子质量相等时，旧版本`T_legacy=T_csv+C、T_force_ref=1473`与新版本`T=T_csv、T_force_ref=mean(T_csv)`应给出相同的初始逐粒子浮力，C=1473−mean(T_csv)。将此作为直接核对，证明热损失变化不是任意改变了初始浮力。若粒子质量不等，旧算术均值与新质量均值的差异另列。冻结温场后旧版本每次重新去均值，而新版本参考固定，不能据此要求整段历史逐点完全相同。

暖启动降温后，不逐步更新T_force_ref；产生的均匀平均体力部分由当前闭域压力响应承担。若后续要用动态参考，需另行给出压力／边界一致性及能量解释，本轮不加此方案。

修改参考标量后按框架增量规则重算浮力通道，保留ForcePrior历史一致性；不得直接清共享总力。

热初值确定以后重建H0，并清零新的热演化累积项；不要在演化途中不断重建H0以消除残差。

### B4. 必过的温度检查

1. 在解析温场`T=1000+100*r/R+50*z/H`上检查双线性采样与host/device一致性。
2. 对真实CSV，导出原始采样T、实际Temperature及其差值；最大差≤1e-3K为本工程单精度目标，双精度使用相应更严容差。
3. CSV范围与粒子采样范围分别报告；粒子没落在极值节点不等于错误。不得要求粒子最大值恰等于网格最大值。
4. 同一位置重复执行加载10次，T与热DeltaT不发生累积偏移。
5. 只改变T_force_ref，Temperature、瞬时热损失、k/cp/μ/σ在同一物理状态不变。
6. 冻结／释放两种路径，排序后按OriginalID核对字段，重采样规则正确。

## 5. C：浮力审计修正，不再除以接近零的总力

### C1. 使用真实参考温度

审计必须读取求解器实际使用的`getForceReferenceTemperature()`或等价共享标量；不再默认t0=1473K。

逐粒子参考力取实际m、β、T、g、T_force_ref。比较专用浮力字段，同时核对其确实进入共享ForcePrior；专用字段正确不等于实际积分总力正确。

### C2. 明确评价量

\[
E_{particle}=\frac{\sqrt{\sum_i|F_{b,i}-F_{ref,i}|^2}}
{\max(\sqrt{\sum_i|F_{ref,i}|^2},F_{floor})},
\quad E_{net,abs}=\left|\sum_i(F_{b,i}-F_{ref,i})\right|.
\]

非零信号：单精度E_particle≤1e-5作为工程目标；总力绝对残差≤`max(1e-5N,1e-5*Σ|F_ref|)`。同时报告相对于Σ|F_ref|的净残差，不只报告相对于|ΣF_ref|。

近零逐粒子参考力时不算虚假相对误差，报告绝对最大残差和总残差。小规模零力测试目标每粒子≤1e-6N、总≤1e-5N；实际规模变化应说明浮点预算，不能结果出来后随意放宽。

旧的约5e-8N总残差应在这种判据下被正确识别，不再显示成相对误差146%的物理失败。

### C3. 最小检查矩阵

| 输入 | 预期 |
|---|---|
| 所有T=T_force_ref | 浮力为零 |
| 均匀T高于参考 | 浮力向上，幅值正确 |
| 正负温差抵消 | 逐粒子正确，净力近零，不误判FAIL |
| T整体+C，T_force_ref也+C | 浮力逐点不变；这是代数检查，不代表传热不变 |
| 只变T_force_ref | 浮力按公式改变，绝对T与热损失不变 |
| 同状态重复exec | 幂等，无力累加 |
| 修改参考后exec | 总力更新正确，无旧通道残留 |

至少一个实际声学积分短测试，确认审计通过的浮力进入运动。若原有测试已覆盖未改部分，复用证据，不重建整套case。

## 6. D：O3顶面总滑移检查，保留有效不穿透证据

### D1. 必须保存真实共享源

打包实际编译使用的两个声学半步及壁面状态实现，而非只提供1st half。核对：顶面标记、压力接触、速度镜像、连续性、Riemann耗散、黏性关系及边角粒子。

现O3声学默认无滑移状态不一定必然错误；连续层面或离散对称性可能使算子等价。但不能仅凭OpenCup黏性为零就宣布总边界滑移通过。

### D2. 固定状态的切向不变性检查

使用相同平面顶面补丁／实际顶面中心区域，排除侧壁、底面和交线影响。比较零速度与均匀切向速度U=(0.01,0,0)、法向速度0的两组相同几何状态。

保持相同密度、压力、外力和相容EOS；测量只执行必要单步／算子，不让不同位置演化污染比较。

分项输出顶部引起的：

- 1st half压力力的法向／切向分量；
- 2nd half连续性增量和耗散力；
- 顶部黏性力；
- 合计切向加速度与机械功。

均匀切向滑动不应凭空产生明显法向密度源或无滑移摩擦。比较完整壁面聚合贡献，不能只用单粒子对的切向分量判失败；SPH各对贡献可能通过邻域对称相消。

顶面压力状态／虚拟速度关系若逐对不是纯法向镜像，应说明聚合后为何满足同样的物理条件，并用下述测试验证，而不强制某一个旧开关。

### D3. 工程判据

采用同一状态的顶部无滑移黏性对照作为有信号参考。顶部残余总切向加速度RMS应≤`max(1e-6m/s², 1%*a_noslip_reference)`；超过1%但≤5%为WARN，>5%需定位。

切向运动导致的连续性增量差RMS，以`ρ0*U/h`归一：≤1%为PASS，1%—5%为WARN，>5%需定位；同时报告绝对值及实际h。

阈值是本项目诊断标准，不是论文公布误差。若边界补丁受侧壁影响无法隔离，先修测试布局或分项统计，不能调阈值。

再做小幅法向扰动短检，确认仍无长期越界／动能自激。只有切向测试通过但法向失效，也不能替换O3主线。

### D4. 失败后的唯一允许路线

如证实某声学项引入不合理切向作用，只修对应的压力／速度状态或连续性构造；先回归切向、法向和静水，再运行F。

禁止为消切向力直接关闭连续性，禁止用Uz kill或z-clamp补法向，禁止重新扫描多个dens_mode。

最多两个针对同一假设的候选修复；若无判别，提交最小算子复现和证据，停止该假设。不能直接把全部流动衰减归因于顶面。

## 7. E：Fig8重建场核查，分开“图像数据”和“补全假设”

### E1. 必须复用的资料链

核查原Fig7/8、Doc204样点／等温线、部分线性重建、Doc205完整ASSUMED温场和DOC224输入CSV。主参考为Chemical Engineering and Processing 47 (2008) 449–455，DOI 10.1016/j.cep.2007.01.013。

DOC224中的完整CSV不是10201个独立文献观测。需要说明从稀疏样点到101×101网格的填充、外推、平滑和温标假设。

### E2. 交付三张核对图及一张表

1. 原Fig8＋数字化点／等温线叠加，标出坐标原点和r/R、z/H方向。
2. 重建T/Tmax＋观测凸包／有效mask＋未观测补全区域。
3. CSV采样到粒子后的实际T/Tmax，与重建网格同图比较。
4. 每个可识别参考点：原值、重建值、误差、像素／线宽不确定度及是否用于拟合。

在有明确标值的参考点上，归一化温度误差目标：RMS≤0.02、最大≤`max(0.05,读图不确定度)`。这是重建一致性检查，不是热求解预测验证；参与拟合的点不能当独立验证点。

不能可靠读值时NOT_EVALUATED，不补造精确值。观测凸包外单独报告面积／体积占比，不当成有证据区域。

### E3. 重点检查会改变环流的几何特征

热点是否位于原图对应区域；轴心、内部环带、冷壁、表面之间的相对温度；径向温度梯度的符号；是否把空白／无数据区填成750K冷块；是否错误镜像、上下翻转或把直径坐标当半径。

未知区域不得新增专门制造第二个涡的热点。若旧重建不合格，只修复有出处的映射／插值问题，保存v1/v2差异和原始样点，不重新拟合到流线结果。

### E4. 温标限制

本轮默认继续保留Tmax=1500K作为明确ASSUMED温标，但Temperature不再平移。CSV最高1425K可能只对应0.95等值线／补全范围，不自动强行抬到1500K。

若无法确认绝对Tmax，绝对温度定量保持NOT_EVALUATED。不要用1473K物性参考点冒充论文Tmax。

本轮不扫描多个Tmax；先弄清同一温标下新旧加载的差异。确需温标敏感性，列为后续一个受控问题。

## 8. F：有限运行矩阵，先隔离温度平移影响

### F0. 相同初态的一次热损失探针

在相同粒子位置和相同边界下比较：

- legacy平移温度（只读复现，不作为验收）；
- 保真CSV温度＋独立T_force_ref。

输出T_min/T_mean/T_max、H0、侧／底／顶部对流和辐射损失、物性范围。该探针无需推进20s，也不要求两者热损失接近。

可增加一个严格代数对照：legacy整体平移量为C时，原T与原参考同时加C，浮力应一致；传热和物性按绝对T变化。用它分清力学参考变化与热学真实变化。

### F1. 修正CSV的冻结温场

Q=0，冻结温场、传热关闭，保持实际T为原CSV，采用B3一次确定的T_force_ref。先0.3→2→5s，若无数值红线到20s，必要时30s。

D验证若需要边界修复，则运行修复后主线；同时保留最少一组O3短程对照，避免同时改温场与边界却无法归因。

主要判断：方向、涡心、环带上升及后期趋势，而不是要求改demean后一定修复轴升。

如果绝对温度与参考一并平移后、恒物性冻结流动差异很大，先查壁面压力参考／实际施力／初始化，不先判温标改变导致真实流动差异。

### F2. 修正CSV的50kW暖启动

配置必须满足：freeze_T=0、thermal_diffusion=1、空间Q=50kW、固定原物性基线、绝对温度保真、参考温度固定、无诊断速度修正。

先0.3→1→5s；正常后到20s。使用dt_scale=0.5；选前5s以dt_scale=0.25复核T、速度、能量和峰值发展。

比较同一物理时刻；若高温初猜引起μ(T)变化而自动收紧黏性时间步，必须输出实际限制项。不要仅凭减半后能跑完就宣布时间步独立。

重点输出：修正前后实际T与损失差异、热点演化、轴／环／缘温度与速度、密度、顶部、分项力、能量。

### F3. 同样冷启动短回归

1473K等温、导热开启、同一Q与边界，先1s，必要时5s。用于证明新参考量分离未意外破坏原有热模块。

不重跑缺导热的旧E1作为正式物理结果。若B/D修改并未涉及该路径且已有完全相同可信结果，可复用并记录身份。

### F4. 后期延长只允许一组

修正后的F1或F2在30s仍发展，但无数值红线且有可解释趋势时，最多选择一组延长至60s；以10s检查点输出。

选择优先级：能直接区分“论文温场能否驱动目标流动”的F1；若F1已明确而F2热演化更有判别力，可选择F2。不得两组都无条件长跑。

H²/ν约23.5s，30s上限不是物理失败边界。物理热时间更长；60s也不能默认热稳态。

## 9. 测量要求：保留224中有效工具，补缺项而不重写

### 9.1 后期专用场

输出完整5s窗：0—5、5—10、10—15、15—20等。覆盖≥99%、样本≥50、最大间隔≤0.2s；时间加权、不向未计算时段外推。

每个完整窗单独保存r-z速度场、T场、mask和涡心／ψ。不能只输出全程平均场，再用窗口标量声称后期拓扑稳定。

轴／环／缘区域沿用224定义与体积权重，另外输出两个方位剖面，检查方位平均是否掩盖三维性。

### 9.2 收支

保留网格Q_up/Q_down/Q_net及ψ回算；至少一条独立3D粒子控制体的质量储量／质量流率路径。

截面收支相对误差只有信号高于背景20倍时评估。≤5% PASS、5%—10% WARN；>10%先定位采样／压缩储量，不能直接改求解器。

Σm/ρ相对静置后漂移≤1% PASS、1%—2% WARN、>2%停止延长并定位。不能以总质量常数代替连续性验证。

ψ重建Eψ≤10% PASS、10%—20% WARN；>20%不得单靠ψ授闭合涡。原速度流线与ψ等值线应一致。

### 9.3 顶面

沿用224预定义的物理顶面与采样层。报告最大超顶距离、超顶质量比例及近顶法向速度。

最大超顶≤0.1dp、超顶质量≤0.1%、近顶法向RMS≤内部代表速度5%为PASS目标；0.1—0.25dp／0.1%—1%／5%—10%为WARN；超过上限或持续累积应停止延长。

边界WARN不等于禁止全部热诊断，但不能授予无保留的预测稳定标签。不允许删粒子或限位使指标通过。

## 10. 验收标准：分开输入正确、数值可跑与物理对齐

以下是预先规定的项目工程门槛，不是文献公布标准。

| 标签 | 必需证据 |
|---|---|
| ABSOLUTE_T_PRESERVED | B检查通过，实际T与CSV一致，无隐式平移 |
| BUOYANCY_AUDIT_VALID | C逐粒子与绝对残差通过，真实参考正确 |
| TOP_TOTAL_SLIP_CHECKED | D声学／连续性／黏性总检查及不穿透复核 |
| FIG8_RECONSTRUCTION_CHECKED | E坐标、样点、mask、补全和温标清楚；缺失项可NOT_EVALUATED |
| SPATIAL_Q_SHORT_NUMERICAL_PASS | T实际演化、Q正确、能量闭合、无数值红线 |
| LATE_FLOW_STABLE | 两相邻完整窗主要强度变化≤5%，场L2差≤10%，涡心变化≤dp，拓扑保持，无误差累积 |
| PAPER_TREND_FROZEN_T | 可信输入下，主次涡位置／旋向与Fig7相符；温度是输入 |
| PAPER_TREND_COUPLED | 释放温场后预测T与环流趋势均有对应证据 |

暖启动短程数值通过不要求P_in=P_loss，也不要求达到LATE_FLOW_STABLE。后期没稳住不能把这个标签一并撤回。

能量闭合：`|ΔH-Ein+Eout|/max(|Ein|,|Eout|,合理能量噪声下限)≤1%`；不得以初始总焓作为大分母掩盖增量错误。实际P_heat对50kW误差≤0.1%。

热准稳态另需完整窗平均|dH/dt|/Pin≤5%、|Pin−Ploss|/Pin≤5%，并与温场稳定对应；本轮不强制获得。

文献涡心／上升带位置的工程比较容差仍用`max(0.10L,2dp,读图不确定度)`，L分别为R/H；必须报告这是粗分辨率趋势标准。

每个指标使用PASS/WARN/FAIL/NOT_EVALUATED，并给value、unit、threshold、reason、原始结果路径。不得把缺参考记为数值失败，也不得把未知记为零误差。

## 11. 固定失败分支，限制Cursor自行扩展

| 结果 | 下一动作 | 禁止动作 |
|---|---|---|
| 修改参考仍改变Temperature／热损失 | 修字段／函数依赖，停暖启动 | 校准h补偿 |
| 总浮力近零但旧rel_err大 | 改逐粒子／绝对判据 | 关闭浮力或改β |
| O3切向测试FAIL | 只查被证实的声学／连续性项 | 全面更换流体框架 |
| Fig8重建与原点不符 | 修坐标／插值／mask并保存差异 | 调热点使流线好看 |
| 保真T后损失仍大 | 查原文边界、初态与实际表面T | 立刻压低h使等于50kW |
| F1符合而F2不符合 | 比较Q、温场演化和热边界 | 重新怀疑已验证的浮力公式 |
| F1仍轴升 | 核对重建梯度与总边界，报告证据 | 只因Tmax假定就任意扫描 |
| 暖启动dt/2仍爆 | 停止，保存首个异常前后分项状态，查黏性／热／声学步限制与排序 | 再改温度均值压住计算 |
| 30s仍局部发展 | 按F4选一组有条件到60s，或标NONSTEADY | 将30s作为算法失败证明 |

热边界参数校准和σ(T)刷新本轮默认不启动。只有以上整改及输入解释清楚、且结果确实指向参数缺口时，列出唯一最高优先后续任务，不扩展成多参数扫描。

## 12. 交付目录与最终报告

```text
doc225_runs/
  A_baseline_and_errata/
  B_temperature_reference/
  C_buoyancy_audit/
  D_top_total_slip/
  E_fig8_input_audit/
  F0_heat_loss_probe/
  F1_frozen_csv/
  F2_warm_spatial_q/
  F3_cold_regression/
  G_late_windows_and_3d_budget/
  ACCEPTANCE_TABLE.csv
  ACCEPTANCE_TABLE.json
  ROUND225_REVIEW.md
  NEXT_DECISION.md
```

保存实际修改的入口、两声学半步、温度／浮力相关模块、patch、构建配置、输入CSV及来源、采样核对、原始CSV、后期专用场、必要VTP和Python脚本。不同运行输出独立，禁止并发共用output。

最终报告必须先回答：

1. 旧demean使实际温度改变多少？新温度是否逐点保真？
2. 三种参考量现在分别是多少，使用者是谁？
3. force_ok误判是否消除，逐粒子和真实总力是否正确？
4. O3总滑移是否通过，是否修改了某个实际声学项？
5. Fig8旧重建有多少真实观测、多少补全，坐标／梯度是否可靠？
6. 相同几何下，去掉绝对温度平移使三项热损失变化多少？
7. 冻结T与释放T结果分别是什么，是否有完整后期窗口？
8. 当前通过的是短程数值耦合、稳定流动，还是论文趋势？下一最高优先项是什么？

**本轮核心：不再把参考温度操作混入绝对温度，不再把近零分母误差当浮力失败，不再只凭黏性测试宣布总滑移成立。先使输入与判定可信，再判断法国文献复现差距。**


<!-- ===== END 225_OPHELIE_ABSOLUTE_TEMPERATURE_SLIP_AUDIT_AND_WARM_RESTART_CURSOR_PLAN_2026-11-18.md ===== -->

---

<!-- ===== BEGIN 226_OPHELIE_TOP_SLIP_FIG8_AND_LOCAL_VISCOSITY_CURSOR_PLAN_2026-11-19.md ===== -->

<a id="226-ophelie-top-slip-fig8-and-local-viscosity-cursor-plan-2026-11-19-md"></a>

# 226_OPHELIE_TOP_SLIP_FIG8_AND_LOCAL_VISCOSITY_CURSOR_PLAN_2026-11-19.md

# 第226轮执行任务书：补齐顶面与温场欠项，验证局部温变黏度，推进法国自然对流

日期：2026-11-19  
项目：冷坩埚招投标／法国OPHELIE自然对流  
输入基线：DOC225_ROUND_PACK_2026-11-19  
执行对象：Cursor。本文是本轮唯一主任务书。

## 0. 本轮必须完成什么

优先解决三个有明确证据的缺口：

1. **补完顶面总滑移检查**：不再只用“OpenCup黏性为零”证明声学、连续性和总边界都正确。
2. **验证Fig8重建一致性**：不再把“来源、补全数量说清楚”当作原始等温线已经核对通过。
3. **审计并接通局部μ(T)**：确认实际黏性相互作用读取逐粒子黏度，而非全域平均值；通过独立解析算例后，比较恒μ与局部μ(T)的冻结温场响应。

这些检查完成后，才用可信温场和已验证黏性路径推进50kW暖启动。

**本轮不能以“已归档、PARTIAL、下一文档再做”替代顶面探针与温场误差核查。** 如确有阻塞，交付最小复现、实际缺项和对其他工作影响，继续完成不依赖该项的任务。

## 1. 为什么现在转向这些问题

### 1.1 保留的有效修复

- Temperature已不再被demean整体抬高约394K。
- T_thermal_ref与T_property_ref保留1473K，Fig8组T_force_ref取初始质量均值并固定。
- 逐粒子浮力审计通过，旧近零总力相对误差只作信息。
- 保真暖启动dt_scale=0.5可运行20s，能量闭合；前5s与0.25标量对照基本一致。
- O3不穿透已有长时证据，保留为本轮候选边界。

不得回退这些修复，不重开浮力符号、共享ForcePrior清零或温度平移方案。

### 1.2 恒μ结果不等于论文物性结果

DOC225冻结与暖启动日志仍为`constant_thermal_rheology=1`、μ=4Pa·s，而输入温场约750—1425K。

法国期刊明确采用温度相关的分段线性物性。表1的4Pa·s是1473K参考值，不是全温区恒定值；原文以μ>1e4Pa·s识别近似凝固冷区，给出的黏度变化范围跨多个数量级。

因此，恒μ=4是求解器诊断底座，不能作为温变物性论文复现的最终判据。它的价值保留，下一步是有控制地增加局部μ(T)，不是推倒原算法。

主参考：Numerical modelling of natural convection in molten glass heated by induction，Chemical Engineering and Processing 47 (2008) 449–455，DOI 10.1016/j.cep.2007.01.013。

## 2. 执行顺序与范围

A 身份与欠项表 → B 顶面总算子探针 → C Fig8重建核查 → D 现有局部黏度与物性曲线审计 → E 黏性独立验证与高黏度推进 → F 冻结温场对照 → G 空间Q暖启动。

B、C、D可分别准备；修改后的边界必须通过B，局部黏性必须通过E，才用于F/G的正式结论。

不修改独立A-phi／Laplace路线；不重新开发电磁算法、线圈、自感或可变形自由面；不扫描μ、β、h与热点位置凑图。

允许在OPHELIE隔离范围内修复、编译、运行。优先GPU；已有CPU隐式路径若更可靠，可先作为准确性对照。计算较慢可以接受，但运行前必须估算成本，不启动明显不可完成的显式高黏度长跑。

## 3. A：锁定配置与结果身份

| 项 | 本轮约定 |
|---|---|
| 几何 | R=0.25m，H=0.185m，z竖直 |
| 基准dp | 0.015m |
| 压力／连续性 | DOC225的固定顶面p′路径，完整重力不重复加入 |
| 顶面 | O3候选，接受B的定向修复，不恢复旧开关包 |
| 密度 | 正式推进不做周期性dens重初始化；保留连续方程 |
| KC／TVF | OFF／OFF，沿用有效底座 |
| 人工运动干预 | tip、zero-net-Uz、Uz kill、z-clamp、ContSilent全部OFF |
| 温度 | 原CSV绝对K保真，不做demean平移 |
| 浮力参考 | Fig8组初始质量平均温度，一次设置后固定 |
| 其余基准物性 | ρ0=2750kg/m³，β=1e-5K^-1；cp/k先保留基准，单独列差异 |
| 物理对照功率 | 50kW空间Q；本轮不混入60kW客户工况 |
| 暖启动时间步 | 从dt_scale=0.5开始；黏性实际限制另查 |

保存git SHA、dirty diff、实际编译源、精度、设备、输入／reload／Q哈希。旧结果不可覆盖。

修改共享算子时用独立黏度类型／策略或受控参数隔离，其他case默认行为不变。

新增`DOC225_OPEN_ITEMS.md`：顶面总滑移未测、温场参考点误差未交、局部黏性读取链未证明；这三项不能因本轮改名而消失。

## 4. B：顶面总滑移探针——本轮必须执行

### B1. 实际源码与状态链

归档两个声学半步、顶面法向、镜像速度、压力外推、连续性、Riemann耗散与黏性接触源码。列出O3下每一个实际启用的分支。

不得仅检查`AcousticWallFreeSlip`标志名称，也不能仅凭声学采用旧无滑移镜像就断言聚合算子必错。以物理不变性测试判别。

### B2. 固定状态对照

顶面中心补丁或可隔离的平面边界，排除侧／底及交线。相同位置、密度、压力、外力，比较：

- U=(0,0,0)；
- U=(0.01,0,0)m/s，法向速度0。

保持EOS与压力初始化相容。固定状态执行必要算子或短步，分项记录：1st-half压力力、2nd-half耗散力、连续性增量、黏性力，以及总切向加速度和功。

完整snapshot/restore或无副作用测量；不得为测试修改主线粒子速度后忘记恢复。

比较完整壁面聚合贡献。SPH单对粒子力可能含切向分量，不以单对非零直接判总边界失败。

### B3. 判据与失败处理

沿用225预设工程标准：总切向加速度RMS≤`max(1e-6m/s²,0.01*a_noslip_reference)`为PASS，1%—5%为WARN，>5%需定位。无滑移参考用同状态下的明确有信号对照。

切向滑动引起的连续性差RMS，以ρ0U/h归一，≤1% PASS、1%—5% WARN、>5%需定位。输出实际h和绝对值。

再做小幅法向扰动与零浮力静水短回归，不穿透不得因切向修复被破坏。

若失败，只修证实存在问题的压力／连续性／镜像项；禁止静音连续性、直接裁切速度或引入局部阻尼。最多两个单假设候选；均无判别则提供最小复现与分项证据。

交付`TOP_OPERATOR_PROBE.csv/json`和`B_RESULT.md`，必须有实际数据。源码归档本身不构成B通过。

## 5. C：Fig8重建核查——必须交付误差，不只归档来源

### C1. 输入链

复用原Fig8、Doc204样点／等温线与部分重建、Doc205完整CSV及meta。当前10201节点的6326插值、1308热点填补、2567凸包外NN数量保留，但这些数量不等于重建正确性。

温标仍为Tmax=1500K ASSUMED，CSV实际最高1425K；不自动把最高值抬至1500K，不用1473K作Tmax。

### C2. 必交文件

1. `fig8_reference_vs_reconstruction.csv`：参考坐标、原T/Tmax、重建值、误差、读图不确定度、是否参与拟合。
2. 原图＋样点／等温线＋重建等温线叠加图。
3. 观测凸包／热点填补／NN区域mask图，给出面积与轴对称体积权重占比，不能只数网格点。
4. 同一物理坐标上的粒子采样T与CSV插值差异图。
5. 轴心、热环带、侧壁与表面的温度剖面，检查径向／竖向梯度符号。

在可靠标值点上：归一化误差RMS≤0.02、最大≤max(0.05,读图不确定度)为工程一致性目标。参与拟合点是重建检查，不是独立热预测验证。

缺少可靠值时NOT_EVALUATED，不猜测补点。原数据无效区与补全假设完整展示。

### C3. 输入修复规则

仅修有证据的坐标、单位、镜像、插值和mask错误；不得根据流线结果调整热点使第二个涡出现。

检查是否将无观测区填成冷块、热点填补是否改变关键梯度、是否把直径当半径、上下方向是否一致。

若旧场不合格，保存v1与修后v2及差异；若没有明确错误，则保留旧场并标出补全不确定性。不能因后期轴升就认定温场必错。

## 6. D：局部黏度读取链审计——先找已有实现，再决定改什么

### D1. 现有资产优先复用

历史快照中已存在`electromagnetic_ophelie_french_material_laws.h`、`ParticleViscosity`、`OphelieDynamicViscosity`、温变物性表、raw/used黏度统计及mu_cap处理。

必须在当前执行仓库查找最新版本，不将历史快照直接覆盖当前代码。先给出文件、类、调用位置和是否真正用于法国入口。

当前入口可见`glass_viscosity.setReferenceViscosity(mu_now)`，这只能证明参考值被更新，不能证明逐粒子μ进入了相互作用；也不能仅根据这一句认定当前完全没有局部黏度支持。

### D2. 把完整链画清楚

`Temperature → μ_raw(T) → μ_used → 注册字段／排序／同步 → 材料闭包类型 → inner μ_ij → wall μ_i → 实际黏性力 → ForcePrior或隐式更新 → 时间步限制`。

对每一环给实际源码与字段名称。特别检查：

- `ViscousForceCK<...,Viscosity,...>`还是`ParticleViscosity`等实际类型；
- `getInterParticleViscosity`与`getOneSideViscosity`读常数还是设备字段；
- 粒子场存在但模板仍绑定常黏度的情况；
- CK内核缓存是否按值保存参考值，字段变化是否能立即被读取；
- 壁面是否仍用全域平均黏度；
- 更新是否发生在本次黏性力计算之前；冻结T阶段是否也按T更新μ；
- 密度、体积、μ字段在排序后是否仍对应同一粒子。

上游基线接口核查参考：
https://github.com/Xiangyu-Hu/SPHinXsys/blob/bc411c43718ca7ecaf4c1ed3e29aefeb4e6090c4/src/shared/shared_ck/particle_dynamics/fluid_dynamics/viscous_force.hpp

该接口通过材料黏度类型取得粒子对和单侧黏度；具体执行以服务器实际源码为准。不要改整个黏性框架来替代正确绑定已有材料类型。

### D3. 物性表来源

导出`mu_T_source.csv`：T(K)、μ(Pa·s)、原出处、原图／页码、插值、端点外推、是否缩放锚定1473K、可靠温区。

优先期刊明确数据；缺项用博士论文／已有数字化补足并标注。核对历史Fig编号，不能把不同文献中的Fig2混用。旧代码注释不作为唯一文献依据。

表1范围1—1e8Pa·s不能直接当完整μ(T)曲线。无法取得可信曲线时，仍可完成E的算法验证；法国温变物性组标ASSUMED／REFERENCE_INCOMPLETE，不编造文献节点。

线性μ插值与logμ插值会不同；按原图坐标和数据来源确定一次并记录，不为得到环流改变插值方式。

### D4. raw、used、cap分开

每个时刻输出μ_raw与μ_used的最小、最大、分位数、体积加权均值；黏度比、截断粒子／体积比例，以及温度超出数据温区的比例。

μ_cap或端点保持必须显式显示。不能把截断后结果标成完全采用文献物性。

## 7. E：局部黏性必须通过的最小验证

### E0. 物理与离散契约

目标是空间变黏度牛顿应力`∇·[μ(∇u+∇uᵀ)]`的适用离散，不能简单把常ν乘以Laplace(u)而忽略μ空间变化。

优先复用现有稳定的对称粒子对黏性。给出μ_ij混合方式、壁面μ取值、量纲、对称性与耗散性。不要未经验证就把几何平均替换成算术／调和平均；算法是否合适由下述界面剪切测试判别。

### E1. 常数退化、局部读取与排序

1. 所有μ_i=4Pa·s时，局部实现与原常μ实现的力／短程速度应一致：单精度相对差目标≤1e-5，双精度使用对应更严容差；近零量同时给绝对误差。
2. 在有已知剪切的同一状态，交换两块区域的μ，保持全域平均μ不变；局部黏性力应按区域变化。若完全不变，说明仍可能使用平均值。
3. 只修改μ字段、不重建几何，重新调用后立即可见更新；按需要正确使CK缓存失效。
4. 排序前后按OriginalID比较μ及力；host/device一致。
5. inner对称黏性对的净内力接近零；无外功的孤立黏性子步不增加动能。壁面运动的功单列，不能把移动壁输入能量误判自激。

### E2. 两层Couette解析验证（必须，不用法国流线替代）

二维平行板，z∈[0,H]，底板静止、顶板沿x以U移动；x周期；无重力、无热、无压梯度。两层等密度，仅黏度不同，无表面张力／界面附加物理。

推荐小规模：H=Lx=0.02m、U=0.01m/s、ρ=2750kg/m³；界面a=H/2；下层μ1=4，上层μ2=40Pa·s。先dp=H/40，必要时H/80。

稳态解析：

\[
\tau=\frac{U}{a/\mu_1+(H-a)/\mu_2},\quad
u(z)=\begin{cases}
\tau z/\mu_1,&z\le a,\\
\tau a/\mu_1+\tau(z-a)/\mu_2,&z>a.
\end{cases}
\]

速度在界面连续，μdu/dz恒定。比较整个速度剖面、界面速度、壁面剪切和两侧界面剪切；不能只比较平均流速。

工程目标：速度L2≤5%、界面速度误差≤5%；远离数值界面各约2h处的剪切误差≤10%；界面邻域完整展示，不因误差大删掉整段。

若粗分辨率不通过，做一次加密确认误差是否减小；不能只调混合系数拟合解析。另用μ1=μ2常黏度退化确认。

ρ只用于瞬态与时间步，稳态解析使用动态黏度μ，禁止混用μ与ν。

### E3. 为何这些测试足够进入工程对照

它们检查常数兼容、局部读取、排序／缓存、应力连续、耗散及实际积分，不是仅复述实现。通过后可进入F，但不代表全部变黏度SPH问题已得到发表级证明。

复用已有等价且身份清楚的验证，补缺项即可，不重新创建十余个测试入口。

## 8. 高黏度如何推进：先算成本，不能盲目显式长跑

### 8.1 必须先输出成本估计

根据实际μ_used,max、ρ、h和实现的黏性稳定条件估算限制时间步及20s所需步数。采用实际公式，典型显式限制量级为`dt_visc∝ρh²/μ_max`，不能按平均μ估计。

论文冷区μ可很高，直接全部投入显式计算可能不可行。若估计需要超过1e7个黏性子步，或预计墙钟成本明显超出本轮预算，不得启动长跑；必须进入已有隐式路径检查。

这个步数是任务成本保护阈值，不是物理稳定性标准。

### 8.2 路线优先级

1. 搜索并复用项目已有的隐式黏性及之前验证结果，确认能接受局部μ。
2. 已有GPU隐式实现可用，则GPU优先。
3. GPU尚不可靠但CPU隐式已有正确实现，可先以CPU完成准确性对照，并明确数据传输／性能代价。
4. 如现有隐式实现只支持常μ，优先扩展其系数读取，先通过E1/E2，不另建第二套压力／流体框架。

不能因用户希望GPU就用严重截断μ冒充正确物理。准确性优先，GPU优化后置。

### 8.3 隐式接入的强约束

如采用后向Euler黏性更新，应写清实际分裂时序和方程，例如`(M+dt*K(μ))u_new=M*u_star+wall_rhs`。K来自已验证黏性离散，不直接照抄符号；检查其对称／耗散条件。

不得同时通过隐式速度更新和显式ForcePrior重复加入同一黏性作用。模式切换必须撤销旧显式黏性通道贡献，保持重力／浮力及其他通道不变。

记录线性残差、迭代次数、收敛失败及每步动能／壁面功。若对称正定条件不满足，不能直接用CG并把“迭代结束”当收敛；优先修已识别离散／边界问题。

时间步不再受显式黏性限制不等于时间精度独立。E2及F关键结果仍做一个dt/2对照。

### 8.4 黏度截断只作诊断

允许沿用已有mu_cap做明确标注的成本／敏感性检查，但必须输出raw/used差异、截断体积比例和对冷区速度的影响。

本轮不扫描多组cap。确需判断截断影响时，只选两档有物理／成本依据的cap，且保留未截断的小规模隐式参考。不能用低cap结果宣称论文冷区已再现。

不自动把μ>1e4区粒子速度清零、不改成固体、不加入新FSI／相变；若后续采用静止冷层近似，需要另列物理建模假设。

## 9. F：同一冻结温场的恒μ／局部μ(T)对照

### F1. 必须隔离变化

选择C已核查的同一温场版本、同一顶面、同一ρ/β/c0、同一初态与时间统计。基准恒μ=4，候选为已通过E的局部μ(T)。cp/k保持原诊断值，冻结温场下不传热、Q=0。

若B修改边界或C修正温场，先补一个短恒μ基准，不能让“旧边界旧温场恒μ”与“新边界新温场变μ”直接比较后归因黏度。

初始T_force_ref一次质量平均后固定；Temperature保真。μ按相同空间冻结T求值，粒子移动后及时更新，不把μ固定在初始粒子ID上。

### F2. 运行顺序

0.3→2→5s；力、密度、顶面和黏性求解正常才到20s。若仍有可解释发展，最多选一组到60s，10s检查点；不同时延长所有组。

选择统一物理时刻及完整后期窗比较。温变μ冷区响应可能很慢，不能以和恒μ完全相同的启动时间要求稳态。

### F3. 必须比较的量

- 原始T(r,z)、μ_raw/μ_used(r,z)与冷区范围；
- 轴、内环、外缘的u_z；
- 两个涡心、旋向、分隔线与上升带；
- 按温度／黏度分层的速度和黏性耗散，而不只看全域U_rms；
- 高黏区是否表现为低运动，不能预设某区必须精确静止；
- 顶面法向、密度／体积漂移、控制体收支；
- CPU/GPU路径、dt与隐式残差。

若局部μ改变了主次涡位置或消除了轴向反转，记录因果证据；若没有，也保留结果，转查温场补全与边界，而非改曲线强行匹配。

## 10. G：空间Q暖启动仅在F解释清楚后运行

沿用50kW冻结空间Q，导热明确开启，温场只在初始化加载一次，随后演化；本轮不恢复σ(T)周期EM反馈。

只在已有正确更新路径上先启用局部μ(T)，cp/k仍保持基准用于隔离黏度影响。不得仅取消`constant_thermal_rheology`就同时改变μ、cp、k、热焓而不说明。

如当前开关捆绑多物性，应拆成清楚的独立策略配置，保持原默认可复现。实际用了哪些物性写进JSON。

μ在温度更新后、下一次黏性求解前更新；隐式矩阵／缓存系数同步，排序字段齐全。高黏区仍按8节成本策略处理。

0.3→1→5s，正常到20s；明确较大散热与热初态的关系。默认不校准h、ε，不把减小损失作为黏度修复验收条件。

能量用与实际cp策略一致的焓；H0在最终初态后设置，绝不反复重建消除残差。实际功率50kW误差≤0.1%，增量能量闭合≤1%作为工程目标。

如果G趋势偏离而F符合，先比较Q与温场演化，不再重写浮力。完整温变cp/k与EM反馈属于后续独立阶段。

## 11. 测量与标准：沿用有效工具，补齐后期与3D证据

### 11.1 后期专用场

完整5s窗覆盖≥99%、样本≥50、最大间隔≤0.2s；每窗输出T、μ、u_r/u_z和mask。不能只输出全程平均。

区域统计体积加权；原始绝对速度、两个方位剖面、方位平均、ψ等值线与流线同时输出。

ψ回算相对误差≤10% PASS、10%—20% WARN；更大时不能单靠ψ授闭合涡。图形插值不等于粒子分辨率提高。

### 11.2 收支与边界

补独立3D控制体质量收支；网格截面Q误差≤5% PASS、5%—10% WARN，低信号不计算虚假相对值。WCSPH储量变化单列。

Σm/ρ相对静置后漂移≤1% PASS、1%—2% WARN、>2%停止延长。质量恒定不替代连续性验证。

顶面仍按225标准：超顶≤0.1dp、超顶质量≤0.1%、近顶法向RMS≤内部代表速度5%为目标；超过0.25dp／1%／10%或持续累积，停止该组并定位。不得裁切粒子。

### 11.3 稳定与论文趋势

相邻完整后期窗：主要强度变化≤5%、速度场L2变化≤10%、涡心变化≤dp、旋向与分隔结构保持、误差不累积，才授予该工况后期稳定。

论文趋势须在可信温场和明确物性来源下检查轴降—内环升—外缘降及两独立回流。位置工程容差沿用`max(0.10L,2dp,读图不确定度)`，L分别为R/H；不能称精确定量验证。

所有门槛为项目工程标准，不是论文公布值。缺可信μ(T)或温标时，算法验证可通过，论文物性对齐保持NOT_EVALUATED。

## 12. 标签与自动判定

每项使用PASS/WARN/FAIL/NOT_EVALUATED，含值、单位、阈值、原因和原始路径；不得用进程退出0或旧passed=1替代。

| 标签 | 必要条件 |
|---|---|
| TOP_TOTAL_OPERATOR_CHECKED | B实际切向与法向测试有数据 |
| FIG8_INPUT_CONSISTENCY_CHECKED | C误差表、叠加、mask及补全风险明确 |
| LOCAL_MU_READ_PATH_VERIFIED | D读取链＋E1字段扰动／排序／常数退化通过 |
| VARIABLE_MU_OPERATOR_VERIFIED | E2解析剖面和剪切通过，耗散／推进无异常 |
| TEMPERATURE_DEPENDENT_MU_COMPARISON | F同温场、同边界、有来源曲线与明确cap对照 |
| SPATIAL_Q_SHORT_NUMERICAL_PASS | G温度演化、能量与数值门通过；不要求热稳态 |
| PAPER_TREND_FROZEN_T | 给定可信温场下论文主次环流对应 |
| PAPER_TREND_COUPLED | 演化预测的温场和流场同时对应，披露物性与EM近似 |

未知不能记零误差；数据说明完毕不能替代一致性验证；使用平均μ不能授予LOCAL_MU标签。

## 13. 固定失败分支，避免扩展任务

| 情况 | 下一动作 | 不允许 |
|---|---|---|
| 顶面探针仍未实施 | 直接补完B后再做新长跑 | 再写PARTIAL并称本轮全部完成 |
| μ场存在但力不变 | 修材料类型、字段读取／缓存 | 只刷新全域平均μ |
| 常μ退化不一致 | 查系数、单位、wall与ForcePrior | 直接上法国高黏度 |
| 两层Couette失败 | 查界面应力、混合方式和空间分辨 | 调μ曲线拟合速度 |
| 真实μ显式成本不可行 | 复用／扩展已有隐式，CPU准确性可先做 | 用低cap宣称文献物性已匹配 |
| Fig8参考不一致 | 修可证明的输入问题并保存版本 | 改热点追流线 |
| 恒μ／局部μ都后期轴升 | 检查温场梯度、边界与3D性 | 宣称浮力公式不可用 |
| F对而G不对 | 查空间Q与热演化 | 同时改μ、β、h和Q |

本轮禁止新开TEAM7、自感、屏蔽、机械搅拌或相变／FSI主题。若出现新缺陷，仅修影响本轮证据的必要项。

## 14. 交付与最终报告

```text
doc226_runs/
  A_baseline/
  B_top_operator_probe/
  C_fig8_consistency/
  D_local_mu_audit/
  E_viscosity_validation/
  F_frozen_mu_comparison/
  G_warm_spatial_q/
  H_late_windows_and_budgets/
  ACCEPTANCE_TABLE.csv
  ACCEPTANCE_TABLE.json
  ROUND226_REVIEW.md
  NEXT_DECISION.md
```

保存实际修改的共享源、入口、材料类型、时间步／隐式实现、patch、物性表来源、原始温场、参考误差表、测试原始曲线、后期专用场与评价脚本。不要只有入口CPP。

每个运行目录独立输出，保存完整命令与退出码。重启若不具备全部力历史／黏度状态，则重跑同初态，不能仅导入位置速度伪称续算。

最终报告开头必须回答：

1. 225欠交的顶面探针这次是否真的执行，哪一项通过或失败？
2. Fig8样点误差是多少、哪些区域仍依赖补全？
3. μ(T)曲线出处和可信温区是什么，有无锚定／外推／截断？
4. 实际黏性算子是否逐粒子读取μ，如何证明不是平均值？
5. 两层Couette与常μ退化是否通过，界面剪切如何？
6. 高黏度用了显式还是隐式，是否GPU，是否发生双计黏性？
7. 同一冻结温场下，恒μ与局部μ(T)的涡心、旋向和冷区速度差多少？
8. 修正后暖启动完成哪一级证据，当前唯一最高优先缺口是什么？

**本轮核心：补完已经明确的欠项，复用已有局部物性能力，以独立解析验证确认真实黏性作用，再解释法国环流。不要再把开关开启、字段存在或数据归档当作实现与验证完成。**


<!-- ===== END 226_OPHELIE_TOP_SLIP_FIG8_AND_LOCAL_VISCOSITY_CURSOR_PLAN_2026-11-19.md ===== -->

---

<!-- ===== BEGIN 227_OPHELIE_ACTUAL_TIMESTEP_PERIODIC_NEIGHBORS_AND_VISCOSITY_REPAIR_CURSOR_PLAN_2026-11-19.md ===== -->

<a id="227-ophelie-actual-timestep-periodic-neighbors-and-viscosity-repair-cursor-plan-2026-11-19-md"></a>

# 227_OPHELIE_ACTUAL_TIMESTEP_PERIODIC_NEIGHBORS_AND_VISCOSITY_REPAIR_CURSOR_PLAN_2026-11-19.md

# 第227轮修正任务书：真实时间步约束、周期邻域与局部黏性验证闭合

日期：2026-11-19  
项目：冷坩埚招投标／法国OPHELIE自然对流  
输入基线：DOC226_ROUND_PACK_2026-11-19  
执行对象：Cursor；本轮以本文为唯一主任务书。

## 0. 本轮唯一主线

先修两个已定位的问题，再判断局部温变黏性：

1. 法国入口的黏性间隔只是触发阈值，没有硬性约束实际积分跨度。
2. 新Couette算例仅回绕坐标，没有完整跨周期边界邻域。

必须补齐常`Viscosity`与均匀`ParticleViscosity`退化对照；真实周期的常μ Couette先通过，再做两层4/40。之后才重新运行法国局部μ组。

本轮不把“仍然爆炸”直接解释成显式路线不可用，也不为了继续推进而降低cap、删除异常粒子或改变误差门槛。

**执行顺序：A 固定证据 → B 真实时间推进修复 → C 均匀μ退化与耗散检查 → D 真周期Couette → E 法国局部μ短程 → F 条件式隐式推进 → G 验收。**

B与D的边界准备可以独立进行；E必须在B/C/D通过后正式运行。已有用户暂停的旧失败任务保持暂停，不自动恢复；使用修正代码和独立目录新跑必要对照。

## 1. 对DOC226结论的修正

### 1.1 保留的有效证据

- 顶面O3总算子探针通过：切向力小，连续性差约0.323%，OpenCup黏性为零。
- Fig8拟合一致性RMS约0.0131，最大约0.0636，保留WARN与补全不确定性。
- 绝对温度保真、浮力参考独立、逐粒子力审计及ForcePrior重置修复保留。
- `ParticleViscosity`已接通，场改变会影响力；壁面黏性接入也有修复。

不要重复开顶面模式扫描，不重做全部Fig8数字化，不退回demean温度平移。

### 1.2 当前不能成立的推断

DOC226 F1c虽然输出`dt_adv≈5.11e-4s`，主循环仍是：

```cpp
const Real acoustic_dt = time_stepper.incrementPhysicalTime(acoustic_time_step);
run_acoustic(acoustic_dt);
if (advection_step(advection_viscous_time_step)) { ... }
```

`TriggerByInterval`只在累计时间大于interval后触发，不截断声学步。本包日志声学步量级约7.04e-4s，已经大于上述黏性刷新间隔。

因此，修“刷新公式和缓存”不等于修“真实积分跨度”。不得写“CFL全部排除，已证明显式局部μ不能用”。

Couette使用普通Inner/Contact关系，hostWrapPeriodicXY只将坐标拉回域，未构造跨边界支持。常μ也失稳，不能用它直接否定变μ离散。

新增`DOC226_ERRATA.md`，保留原始失败数据，但标记`TIME_SCHEDULING_NOT_FULLY_CONSTRAINED`与`PERIODIC_NEIGHBORS_INCOMPLETE`。这是修正解释，不是删除失败。

## 2. 范围与禁止项

固定R=0.25m、H=0.185m、dp=0.015m、O3候选顶面、当前p′／连续性路径、KC OFF、TVF OFF。

保留绝对T、独立T_force_ref、物性来源与raw/used/cap统计。初始Fig8参考温度一次质量加权确定后固定。

不修改独立A-phi／Laplace分支。不做新TEAM7、自感、屏蔽、线圈、搅拌或可变形自由面工作。

禁止tip、zero-net-Uz、Uz kill、z-clamp、ContSilent；禁止通过改变μ(T)表、β、h、Q或热点位置修数值错误。

优先GPU，但准确性优先。周期Couette可先在成熟CPU实现完成；不要为了全GPU另造不完整周期条件。框架内部周期镜像／索引映射不等于用户必须手工布置工程ghost壁粒子。

本轮共享库改动应最小化。法国入口的时间组织先在本case隔离修正；不得直接全局改变所有TimeStepper触发语义。

## 3. A：固定代码、调用与运行身份

记录SHA、dirty diff、实际编译入口和共享源、Real精度、设备、reload／T／Q哈希及最终生效配置。

必须提供：

- 法国入口实际时间循环、温度／μ更新、黏性力更新、两个声学半步、位置／邻域更新时间；
- `TimeStepper::incrementPhysicalTime`、TriggerByInterval及黏性时间步类；
- 黏性inner/contact/ForcePrior更新链；
- Couette周期处理、壁面速度字段、邻域关系和评价函数。

用一个表列出“哪些函数改变物理时间、哪些只是测量、哪些改变力”，不能再把诊断调用当无副作用。

旧结果保留，所有新运行独立目录；并发运行不得共用output。

## 4. B：把黏性限制落实到真实推进

### B1. 首先增加实际调度审计

新增`actual_step_audit.csv`，至少包含：

```text
t_start,t_end,dt_actual,dt_acoustic_candidate,
outer_start,outer_budget,outer_elapsed,outer_remaining,
dt_viscous_candidate,mu_raw_max,mu_used_max,
viscous_force_last_update_time,viscous_force_age,
position_update_count,neighbor_update_count,viscous_update_count
```

关键是记录同一份显式黏性力在实际积分中保持了多久，而不只是声学步有多小。

若黏性力外层更新、声学内层重复使用，实际黏性作用跨度是外层区间总长。只让每个声学步≤dt_visc而整个力刷新间隔仍超限，也不算修复。

先对旧F1c做短复现／调度复现，输出首次跨越限制的证据，不必把旧爆炸重新跑满。

### B2. 首选实现：匹配外层区间，截断声学末步

复用SPHinXsys成熟双时间步组织。伪代码用于表达约束，实际调用须保持本case既有物理时序：

```text
while t < end_time:
    使用当前粒子状态、T、μ和邻域
    更新必要的黏性系数及实际限制
    Dt = min(inertial_limit, viscous_limit, end_time - t)
    保存本外层区间的起点和限制
    恰好一次构造／更新当前外力和黏性作用

    integrated = 0
    while integrated < Dt:
        remaining = Dt - integrated
        dt = min(acoustic_candidate, remaining, end_time - t)
        校验 dt > 0 且没有超过限制
        只增加一次物理时间
        执行两个声学半步
        按原约定更新浮力等声学阶段量
        integrated += dt

    恰好一次完成位置／排序／周期／邻域更新
    按该区间实际时间完成一次热更新（若开启）
    冻结空间T则在新位置重采样，再更新下一段需要的μ
```

重要约束：

1. 禁止新的显式外层循环与旧advection_trigger分支同时执行相同物理操作。
2. 不能既手动更新时间又调用会增加时间的旧接口。
3. 不能每声学步和每外层区间重复累积同一段热源／边界损失。
4. 更新位置后必须更新cell／邻域，再计算依赖新位置的力；保持evolving字段排序一致。
5. 外层阶段的μ/T保持与力系数一致；若必须中途修改μ，则先结束／重建该区间预算，不沿用旧上限。
6. 需要提前重建邻域时可提前结束外层区间，但按真实已积分时间处理，不能假定完成Dt。
7. 零剩余量用与时间精度相符的容差结束；不得通过人为跳过一段非零物理时间消除小步。

先只在法国OPHELIE模式中启用受控实现。已验证方腔等case不需要一起重写。

### B3. 统一触发器的职责

物理强约束由实际积分预算实现；输出、统计、EM刷新等事件可以使用触发器。

禁止通过`present_time > interval`改成`>=`就宣称解决：跨越间隔的声学大步仍然存在。

如果采用框架已有matched-interval帮助函数，必须确认它确实对最后一个子步取min，并输出实际累计长度；不能仅凭函数名推断。

### B4. 明确黏性限制的来源

先使用实际μ_used,max、ρ和h计算已有黏性限制，不能用平均μ。核对wall项的系数是否已包含在经验CFL安全系数中。

对冻结几何的无外功黏性子步，检查动能不增加；如果满足经验h²/ν限制仍自激，进一步检查实际离散的谱／行和稳定界，而不是盲目再减cap。

对已证明为正系数对称扩散的离散，可按`a_i^visc=Σ_j c_ij(u_j-u_i)+wall`的实际系数估计速率上界；需包含壁面，给出符号与量纲推导。不可未经检查照搬某个常数。

### B5. 必过自动检查

| 条件 | 判定 |
|---|---|
| 每个dt_actual≤声学候选与外层剩余 | 必须满足，仅允许浮点容差 |
| 黏性力实际使用区间≤其有效限制 | 必须满足 |
| sum(dt_actual)=实际物理时间增量 | 必须满足 |
| end_time不被一个完整大步跨越 | 必须满足 |
| 一次外层只更新一次对应位置／热账本 | 必须满足 |
| μ突变后下一段预算实时收紧 | 必须满足 |

设置一个确定性调度小测试：声学候选大于、小于、非整除外层Dt三种情形，核对最后一步截断；它验证调度，不代替流体验证。

B通过后保存旧／新审计对比，才允许写`ACTUAL_TIMESTEP_CONSTRAINED=PASS`。

## 5. C：先完成均匀ParticleViscosity退化与耗散检查

### C1. 两路径同状态比较

固定相同粒子、邻域、速度、压力、边界和μ=4Pa·s，对比：

- 原常`Viscosity`路径；
- 全部μ_i=4的`ParticleViscosity`路径。

分别比较inner、wall、更新后的ForcePrior；必须让两路径壁面贡献采用相同更新次序，不能一条漏wall、一条包含wall后比较。

单精度非零力相对差目标≤1e-5；双精度更严。近零力同时用绝对误差与参考尺度，输出逐粒子最大残差。

随后用修后的同一调度进行0.3s实际速度对照，目标主要速度量差≤1%；若不符，先查材料类型、系数、缓存、接触顺序及更新频率。

### C2. ForcePrior与壁面贡献

确认inner→wall→当前黏性通道更新的实际顺序，保证wall剪切进入总力，且不会执行两次相同贡献。

重复同状态exec应幂等；排序后按OriginalID比较；从旧黏性模式切换时旧current/previous贡献必须一致撤销／重建，不直接清共享总力。

### C3. 纯黏性耗散

无移动壁外功、无浮力的冻结几何上给平滑初速度，单独检查黏性推进的动能与内力。使用相同预算，动能不得出现超出浮点误差的持续增长。

对移动壁工况，壁面功单列，不能把能量输入误认为自激。Couette真实驱动由实际壁面速度字段进入，不以日志U非零为证据。

C没过时不进入法国μ(T)，也不开展高黏度长跑。

## 6. D：修复Couette周期邻域与评价

### D1. 真周期的最低要求

跨x、y周期面的流体粒子必须互相进入核支持范围，距离用正确周期位移或等效镜像邻居；包括角部同时跨两方向的关系。

只wrap坐标不满足要求。边界附近缺邻居不能靠缩窄采样区域、增大横向宽度或删除离群粒子掩盖。

优先复用现有成熟周期算例／CPU周期实现；可把验证改为二维x周期平行板，避免新增不完整3D薄层。二维验证与法国三维使用同一黏性核心或明确对照。

若采用框架内部周期镜像支持，记录虚拟粒子同步、排序、去重和反力累计；不要求用户手工布置工程ghost边界。

### D2. 先做周期支持预检

规则粒子、均匀密度、常速度场，比较周期面两侧与内部的：邻居数／核和、密度重初始化误差、压力残差、黏性力。

对一个粒子整体平移一个周期长度，物理邻域与力应等价。检查同一邻居是否被重复加入，周期镜像质量不能进入真实总质量统计。

先通过这个预检，再跑Couette。否则标`INVALID_PERIODIC_TEST`，不允许写黏性算子失败。

### D3. 壁面驱动核对

确认实际黏性算子读取的壁面平均速度／速度字段正是顶板设置的U；压力壁、黏性壁和法向位置一致。底板0、顶板U，周期方向无侧墙。

初始条件与对比对象一致。可从静止推进到稳态，也可从解析速度作保持／扰动检查；若仍是瞬态，不与稳态解析直接判误差超标。

### D4. 常μ先过，两层随后

沿用226的小规模几何：H=Lx=0.02m，U=0.01m/s，ρ=2750kg/m³，界面a=H/2；先μ1=μ2=4，再μ1=4、μ2=40Pa·s。

稳态两层解析：

\[
\tau=\frac{U}{a/\mu_1+(H-a)/\mu_2},\qquad
u(z)=\begin{cases}
\tau z/\mu_1,&z\le a,\\
\tau a/\mu_1+\tau(z-a)/\mu_2,&z>a.
\end{cases}
\]

常μ退化应为线性速度。先H/40，只有周期正确且粗粒子误差合理时才用H/80复核。

后期稳定需至少两段窗口的剖面变化小；按实际ν和H估计响应时间。不能因某个中间时刻恰好符合稳态解析就选它通过。

### D5. 原始数据为主，禁止异常值过滤决定PASS

主评价使用预先固定的物理采样区域，保留所有真实粒子。原始L2、最大速度、横向速度、密度、异常粒子数量都输出。

可额外输出稳健误差用于定位，但绝不能由删除高速度粒子后的误差授予通过。若同一区域存在明显异常速度，先判数值状态不合格。

工程目标：原始速度L2≤5%，界面速度误差≤5%，远离界面约2h的两侧剪切误差≤10%；界面邻域曲线仍完整展示。

检查μdu/dz连续、壁面剪切与解析τ对应。不能只以界面速度一项通过宣布局部黏性通过。

粗细两档应呈一致物理解与误差改善趋势；两档不构成严格收敛阶证明。

## 7. E：法国局部黏性只做有限、有判别力的短程

### E1. 前置门

必须已满足B真实步约束、C常μ退化、D真实周期与解析验证。否则本节只允许用于定位的短探针，不得形成物理对照结论。

冻结同一Fig8温场、Q=0、同一O3、同一物理参数和初态。温度保真，T_force_ref固定，cp/k保持基准以隔离μ。

### E2. 最小运行矩阵

| 组 | 黏度路径 | 目的 | 时长 |
|---|---|---|---|
| E0 | 原常μ=4 | 修后调度基准 | 0.3s |
| E1 | ParticleViscosity均匀4 | 材料路径退化 | 0.3s |
| E2 | 局部μ(T)，沿用诊断cap400 | 检查原F1c失败是否仍存在 | 先0.05s，通过到0.3s |
| E3条件式 | ParticleViscosity均匀400 | 仅在E2仍失败时区分整体刚度与空间对比 | 0.05s |

不默认增加cap40或十几个cap值。cap400只是历史可追溯的诊断，不是论文物性。

如果E2仍失败，保存首个异常前后实际dt、黏性力年龄、μ、ρ、压力、分项力与动能，定位异常最先出现在墙、冷区、界面还是内部。

E3也失败：优先查高黏度推进／谱稳定界；E3稳定E2失败：优先查变系数离散、黏度跳变、字段更新与壁面取值。不可只根据一组爆炸决定换整个流体框架。

### E3. 修后时间步复核

仅对有意义的E2做一次外层预算安全系数减半对照。真实步长／力保持区间同时记录，不能只改命令行数字。

短程主要速度与能量趋势一致、无异常增长，是继续的必要条件；0.3s运行正常仍不能授予论文环流通过。

## 8. F：何时转隐式，以及如何复用

### F1. 转隐式条件

满足以下任一项可以转：

- B/C/D已通过，实际约束满足但高黏度显式仍因已解释的刚度成本不可行；
- 完整文献μ(T)所需步数／墙钟成本明显过高。

不能用“未修正时间调度的爆炸”作为充分证据。即使显式cap400修后可用，也不意味着完整冷区高μ显式可行。

### F2. 复用顺序

先检查已接通局部μ的`test_3d_ophelie_french_natural_reduced_pressure`及已有`ParticleViscosityHarmonicDampingRate`／InteractionSplit隐式阻尼。

GPU隐式可靠则优先GPU；CPU已有可靠路径可先作为正确性对照。移植黏性推进策略，不复制第二套压力框架、不回退温度或边界修复。

隐式路径必须先通过相同的常μ退化和Couette测试，再用法国冻结温场。不能仅以残差小或不爆炸判断准确。

同一黏性作用只积分一次：隐式速度更新启用时，旧显式黏性ForcePrior贡献应正确撤销／关闭，不允许叠加。

记录实际分裂时序、壁面项、守恒／耗散、迭代残差与时间步精度。只针对黏性子步，不额外迁移全套流体力学。

### F3. 完整μ与cap披露

保留μ_raw/μ_used、截断体积、曲线来源／锚定／外推。当前cap400截断了超过一半粒子，不能称论文冷区已对应。

真实曲线不完备时算法可验证，但论文物性仍ASSUMED／REFERENCE_INCOMPLETE。不编造μ(T)节点，不通过cap拟合流线。

## 9. 后续热与电磁位置：本轮先不扩大

本轮优先完成黏性闭合。通过后若有余量，可做冻结T到2—5s的有条件对照；不要求马上启动20s暖启动。

下阶段顺序为：可信局部μ冻结温场 → 50kW冻结电磁空间Q下温度／流动演化 → σ(T)及EM刷新。

明确配置标签：

- `fixed_temperature`：指定T，Q=0；不是电磁加热自然对流。
- `frozen_em_q`：电磁来源空间Q固定，温度演化；不是无电磁，也不是完整反馈。
- `temperature_em_feedback`：T→σ(T)→EM→Q→T实际刷新。

论文采用交替耦合与σ欠松弛，未给一个能直接照搬的物理秒数刷新间隔。282kHz与原文300kHz、no-self-induction与论文模型差异单列，暂不同时修改。

## 10. 验收标签与标准

所有数值门槛为本项目预先规定的工程标准，不是论文原始误差标准。

| 标签 | 必要证据 |
|---|---|
| ACTUAL_TIMESTEP_CONSTRAINED | B实际每步与力保持区间均受限，调度测试通过 |
| CONSTANT_MU_PATH_EQUIVALENT | C同状态力及实际推进退化通过 |
| PERIODIC_NEIGHBORS_VERIFIED | D跨边界支持、平移不变性与无重复质量通过 |
| COUETTE_CONSTANT_MU_PASS | 常μ原始剖面、稳态与无异常通过 |
| COUETTE_VARIABLE_MU_PASS | 两层原始速度、界面与剪切通过 |
| FRENCH_LOCAL_MU_SHORT_PASS | E2真实受限短程通过，dt复核有证据 |
| IMPLICIT_LOCAL_MU_VERIFIED | 若走F，独立验证与不双计通过 |

不授予论文双涡、长期热平衡或完整EM反馈标签。实际结果用PASS/WARN/FAIL/NOT_EVALUATED，附值、单位、阈值、原因与数据路径。

旧稳健／过滤误差只能作为辅助列，不参与主PASS。缺数据不能默认PASS。

## 11. 运行与停止规则

1. 常μ退化未过，禁止继续变μ法国长跑。
2. 周期预检未过，禁止继续加密Couette。
3. 首次明显数值异常时保存前后状态并停止，不能把爆炸跑满预设时长。
4. 持续动能增长不能仅以“不穿墙”判为正常；无驱动黏性子步应耗散。
5. 分项力账本失配、负／非有限密度、时间预算越界时直接FAIL。
6. 每个缺陷最多两个有明确假设的候选；未改善证据则停该假设，提交最小复现，不叠新开关。
7. 不重做已经通过的顶面与Fig8全部工作，只做受本轮修改影响的必要回归。
8. 保留方腔等已验证路径；共享库有改动时做其关键短回归，不开展无关全库长测。

## 12. 交付目录与必须回答的问题

```text
doc227_runs/
  A_baseline_and_errata/
  B_actual_timestep/
  C_uniform_mu_equivalence/
  D_periodic_support/
  D_couette_validation/
  E_french_local_mu_short/
  F_implicit_if_needed/
  ACCEPTANCE_TABLE.csv
  ACCEPTANCE_TABLE.json
  ROUND227_REVIEW.md
  NEXT_DECISION.md
```

保存实际源码、patch、输入身份、编译命令、退出码、调度原始CSV、周期支持数据、原始速度剖面、未过滤误差、异常快照及评价脚本。

最终报告开头直接回答：

1. 原F1c是否实际跨过黏性预算？旧／新实际跨度分别多少？
2. 现在黏性力在两次刷新之间使用了多久，是否满足其稳定限制？
3. 均匀ParticleViscosity是否与原常黏度一致？壁面贡献是否只更新一次？
4. Couette是否真有跨周期邻域？常μ是否先通过？
5. 两层Couette原始误差与剪切是否通过，有无删除异常值？
6. 修后cap400失败是否消失；若仍失败，首个异常位于何处、哪项先增长？
7. 是否需要隐式，依据是正确性问题还是计算成本？是否复用已有实现？
8. 当前可以授予什么标签，下一最高优先项是什么？

**本轮结束条件：真实时间约束、材料路径退化、周期验证和局部黏性解析验证形成闭合证据；不能再把未合格的小算例失败当成法国自然对流路线失败。**


<!-- ===== END 227_OPHELIE_ACTUAL_TIMESTEP_PERIODIC_NEIGHBORS_AND_VISCOSITY_REPAIR_CURSOR_PLAN_2026-11-19.md ===== -->

---
---


<!-- ===== BEGIN 228_OPHELIE_LINUX_ELASTIC_PADDLE_SOFTWARE_ADDENDUM_2026-11-20.md ===== -->

# 228｜OPHELIE Linux 软件开发补充任务书：不锈钢弹性搅拌桨与接入前测试

日期：2026-11-20  
执行对象：Cursor / OPHELIE 开发人员  
配套文件：`OPHELIE_DESKTOP_IMPLEMENTATION_CURSOR_PLAN_2026-11-19(1).md`、甲方几何与物性原包、`OPHELIE_BRANCH_STATUS_2026-11-20.md`。

## 0. 本补充件的优先级和本轮目标

本文件补充原软件任务书；冲突处以本文件及用户最新指示为准。不要重新创建一套软件架构。

用户已明确：

1. 最近只交付 **Linux 软件**。暂停 Windows 移植、SYCL CPU Windows 后端和 Windows 安装包任务。
2. **先将甲方机械搅拌项目中的桨改为真正参与力学求解的不锈钢弹性体，完成粗分辨率测试，再接入该项目的软件运行功能。** 可以提前搭建不依赖求解结果的 GUI、装配预览和文件接口。
3. 软件仍为一个程序、两个示例项目；自然对流不含桨，机械搅拌含弹性桨。
4. 当前按用户与甲方讨论口径，固体力学和传热学分别独立验证。不要在本轮追加疲劳寿命、完整热应力或全面强耦合验收任务。
5. 原始运行证据已由用户在本机保存各轮 round pack。整理证据索引即可；不要因此要求重跑全部历史测试。
6. TEAM7 尚无完整通过证据，仍标记待完成。软件可开发，不得写成电磁模块已全面达到合同 10% 指标。

本轮完成目标：**甲方真实桨几何 → 弹性体粒子 → 不锈钢材料 → 轴端驱动 → 流体载荷双向交换 → 真实固体积分 → 可复跑短测 → Linux 软件项目接入**。

不得以“材料字段改成 stainless steel，但仍整桨刚体旋转”作为完成。

## 1. 开工检查与版本固定

先检查实际工作树，不根据状态摘要推断代码已实现。

记录：分支、Git SHA、工作区改动、编译器和构建目标、SYCL 设备、固体/流体精度、实际配置及几何校验和。保留用户已有改动；不硬重置、不强推、不删除旧结果。

优先沿用 `feature/EM-Biot-Savart`，可在隔离工作树/功能分支实施后集成。保持独立 A-phi/Laplace 研究分支不受影响。

本轮需要查明并记录：

| 内容 | 必须确认的实际实现 |
|---|---|
| 甲方项目入口 | 当前甲方装配在哪个目标；若尚无，沿现有工程入口接入，不假定法国尺寸就是甲方尺寸 |
| 原桨运动 | 哪个类覆盖 Position / Velocity / Normal；Simbody 或运动学约束在哪里执行 |
| 弹性固体 | 已有材料、本构、应力积分、参考构形修正、固体时间步和约束类 |
| FSI | 压力和黏性载荷传递、壁面平均速度/加速度、接触邻域更新的已有成熟算例 |
| GPU 支持 | 上述每个环节实际是 GPU、CPU 还是混合；不能仅由 SYCL 主程序推断 |
| 近期修复 | ForcePrior 生命周期、真实时间步约束在机械搅拌入口是否生效 |

提交 `BASELINE_AND_CAPABILITY.md`，列出采用的成熟参考算例和源码路径。先看可复用实现，再作局部改造，不重写通用 FSI 框架。

## 2. 改造范围：轴端规定运动，桨体弹性响应

### 2.1 必须实现的物理含义

- 桨采用三维弹性固体，内部粒子的运动来自应力积分、流体作用力及已启用的外力。
- 电机驱动通过**轴端夹持区/明确的驱动区域**施加；该区域作给定轴旋转。
- 非夹持区域允许真实弹性变形，流体边界使用实际变形后的桨位置、速度和法向。
- 旋转角度由角速度时间积分得到；启动阶段使用平滑斜坡，不能把变转速过程写成 `theta = omega(t) * t`。
- 区分材料点位移与相对刚体旋转的弹性偏移。旋转产生的大位移不能直接显示为“大变形”。

不要同时用整桨刚体运动类和固体应力积分更新同一批自由粒子。约束只覆盖选定夹持区域；禁止每步把所有固体位置重置到解析旋转位置。

### 2.2 推荐实现路线

首选已有经过验证的 SPHinXsys 弹性 FSI 与规定运动约束组合。复用能正确处理有限转动的本构与应力更新，避免使用仅适合全局小转角的线性应变实现进行整圈旋转。

若采用旋转参考系，必须完整处理所需惯性项和载荷/速度变换，并给出成熟依据。本轮不为赶工引入未经测试的新旋转坐标算法。

保留旧刚性桨为显式的回归配置，例如 `paddle.mechanics=rigid_prescribed`。新甲方搅拌预设必须选择 `elastic_driven`。该字段是建议接口名，可适配现有命名；必须只有一个生效来源。

### 2.3 轴、夹持区和装配

沿用原任务书的统一 STL 变换。机械轴以甲方装配为准，不能将偏心桨强行移至坩埚中心。历史估计偏心约 -67 mm 仅供核查，不是未经检查的硬编码。

驱动区优先选轴顶部明确的一段实体截面。输出驱动区粒子数、包围盒、轴向范围、轴点和单位轴向量；导出与装配叠加截图。驱动区不能为空，也不能误选整个桨。若 STL 不含完整夹持段，选择现有轴上可解释的驱动区域并记录近似。

各计算环节统一坐标，不要让 CAD +Y 轴、求解器 +Z 轴和扭矩取矩轴混用。

## 3. 不锈钢材料默认值与配置

用户尚未指定牌号时，采用 **304/1.4301 室温常物性开发预设**，允许后续替换。不是高温寿命或实际桨材料认证。

| 参数 | 开发默认值 | 单位与说明 |
|---|---:|---|
| 密度 rho_s | 7900 | kg/m³ |
| 杨氏模量 E | 2.00e11 | Pa；输入 200 GPa 只转换一次 |
| 泊松比 nu | 0.30 | 本轮工程假设，需在来源中注明 |
| 导热系数 k_s | 15 | W/(m·K)；仅在开启固体热模块时读取 |
| 比热 cp_s | 500 | J/(kg·K)；同上 |
| 线膨胀系数 alpha_s | 16e-6 | 1/K；预留，不代表已实现热应变 |
| 弹性应力参考态 | 初始无应力构形 | 在求解坐标中构建，不能把初始装配旋转当作应变 |

rho、E、k、cp、alpha 的参考：Outokumpu Core range，Core 304/4301，物性表。
https://www.outokumpu.com/-/media/files/products/core/outokumpu-core-range-datasheet.pdf

材料标记：`source_status=development_assumption`；保存原始值、SI 值、参考温度和常量/温变模式。

第一轮采用等温力学，固体与流体可统一为 293.15 K，浮力与热源关闭。流体采用当前已有可用常物性预设并记录，不能误读法国高温玻璃参数为室温实验实物性质；这是功能短测。

后续暖态测试如沿用室温 E，必须注明“常 E 功能测试，非高温材料预测”。不自动外推室温钢材参数到 1400 K 以上作为工程强度结论。

不要因为钢有导电性就自动把桨加入 EM 导体域；首版保持原有电磁模型范围，是否计算钢桨涡流另列任务。导电性、热膨胀、塑性、疲劳不因填写材料名而默认启用。

## 4. 粒子生成与粗分辨率约束

使用用户提供的真实桨 STL；采用实体粒子，不能用仅供展示的线圈或表面采样粒子替代。

1. 检查连通性、封闭性、最薄叶片厚度、轴径、叶片与轴连接、初始浸没和整圈扫掠。
2. 尽量使最薄承载厚度至少具有约 3–4 层实体粒子。这是本轮可解析性准入检查，不代表应力已经收敛。
3. 粗分辨率不能让叶片断开、消失或只剩一层却仍宣称弹性计算有效。若不满足，优先局部/固体加密，使用框架已有多分辨率接触能力；没有成熟支持则采用可接受的统一间距，不临时发明耦合算法。
4. 原任务书约 9.3 mm 的桨壁间隙需要按实际装配重新计算。约 5 mm 流体间距只能用于受限运动冒烟，不能据此授予近壁流场或扭矩精度通过。
5. 流体初始体积应扣除桨的真实浸入体积；自然对流项目仍不生成桨求解体。
6. 输出固体/流体/壁面粒子数、质量体积、厚度层数、最小间隙与 dp 比、夹持区及流固接触支持诊断。

若原始薄片在可承受粒子数下不可解析，标记 `BLOCKED_UNDERRESOLVED_SOLID` 并交最小修正建议；不得擅自加厚叶片、换成壳模型、降低 E 或退回刚体并标 PASS。

## 5. 流固耦合和时间推进：必须核对的实现点

### 5.1 真实双向载荷

复用已有压力/黏性流体作用力到固体的路径，并用于固体动量推进；自由固体粒子的实际运动返回流体接触边界。

- 压力与黏性载荷分别可追踪；禁止只输出“估计载荷”而未作用于固体。
- 说明固体采用的应力张量，以及流固力是否已经进入总力，避免 ForcePrior 或外力重复计入/被清零。
- 固体亚步完成后，按成熟 FSI 算例更新给流体使用的平均壁速度/加速度和法向。
- 粒子移动后，在需要使用新邻域之前更新对应 cell 与接触邻域；固体参考构形邻域与当前流固邻域分开管理，不能全部随意重建或全部冻结。
- 排序必须保持材料、应力、参考位置、夹持标记及热字段的对应关系。

### 5.2 固体亚步与计算成本

真实不锈钢 E 会引入远小于流体步长的弹性时间步。不能把流体时间步直接用于固体，也不能通过降低 E 冒充真实不锈钢来提速。

按已有 FSI 调度，在每个流体耦合区间内对固体作亚循环：实际固体步长不超过固体稳定限制和该区间剩余时间，最后一步截断，累计时间精确覆盖耦合区间。约束、载荷施加、应力半步和速度/位置更新顺序沿用已验证实现。

记录流体实际 dt、固体 dt 最小/最大、每个区间固体亚步数、步长越界次数、真实物理时间和墙钟耗时。第 227 轮实际步长修复不能只在自然对流打开；机械搅拌需单独确认。

GPU 优先；若成熟固体/FSI 环节目前只能在 CPU 执行，允许先完成准确的 CPU/GPU 混合短测，明确显示各段后端与同步成本。不得假称全 GPU，也不在本轮为此重写固体 GPU 后端。

## 6. 测试顺序和判定：先单模块，再甲方几何

以下阈值为本轮开发验收，不替代合同中双方确认的正式指标。每个测试保存命令、配置、生效参数、日志和原始结果。

### S0：固体力学独立验证

先复跑一个已有成熟弹性基准，例如小变形悬臂梁静挠度。优先复用现有算例和参考，不重复造测试框架。

若采用端部合力 F 的细长矩形梁，参考挠度 delta=F*L³/(3*E*I)，I=b*t³/12。载荷实际分布、有效长度及边界须与解析假设匹配；不能将不同加载方式混用。施加平滑载荷，确认已稳定后比较，不用任意瞬时值代替静挠度。

本轮要求：挠度相对误差 ≤10%，结果有限、约束明确、非零加载产生非零响应。应力若比较，应远离夹持和加载奇异区；不要求奇异点峰值满足 10%。已有同版本证据可复用，但须核对确实使用本轮所选材料/算子路径。

### T0：传热独立验证

按用户口径独立验证，不要求先把热应变并入桨。优先复用已通过的导热解析基准。可选均匀材料一维平板瞬态导热、其他面绝热：

`T(x,t)=T0+A*sin(pi*x/L)*exp[-alpha*(pi/L)^2*t]`，alpha=k/(rho*cp)，两端固定 T0，初态同式 t=0。

采用同一独立热模块路径、匹配边界与材料，比较非零时刻温度扰动的相对 L2，而非用巨大绝对温度掩盖误差。至少一个有明显衰减的时刻误差 ≤10%，附合理分辨率检查。已有对应证据可复用。

T0 可与桨机械测试独立开展，不阻塞无热力学的 E0/E1。GUI 不得据此显示桨热应力或热疲劳已完成。

### E0：真实桨静止保持，先隔离固体初始化

真实桨、选定夹持区、等温、零转速、无流体载荷、重力等外力关闭。跑至少 100 个有效固体步，检查无 NaN、无翻转/负 Jacobian、初始应力近零、无自发显著运动。记录无量纲应力 sigma/E 与位移/dp 的范围；不要把任意单精度残差当结构失效。

### E1：真实桨非零载荷测试，证明没有被锁成刚体

保持夹持区，给叶片指定小幅平滑已知载荷，运行卸载或响应窗口。至少两个非零载荷级，报告载荷和弹性偏移响应；小变形范围内应体现合理的近似比例关系。

要求：自由区确实由固体求解器积分、存在可分辨的非零弹性响应、应力和内力实际更新、无整体位置覆盖。若单精度下响应过小，使用独立加载测试调整测试力至可分辨但仍小应变；不偷偷降低 E。该项不是实际运行强度工况。

### E2：真实桨规定旋转，无流体或关闭流体作用

对夹持区施加平滑转速斜坡，先作 0.05–0.1 s 冒烟。检查旋转轴、符号、约束及固体亚步正确后，再作至多一轮有明确时长的延长测试。

整圈测试时长按 `integral omega(t) dt` 达到 2*pi 判定，不能只用含斜坡情况下的 `60/rpm` 代替。可先报告短角度通过，再补整圈；未跑满不能标整圈 PASS。不得为了缩短墙钟时间任意提高转速。

检查弹性偏移（相对驱动刚体参考轨迹）、应变、应力、粒子支持和稳定性。实际旋转有惯性载荷，不能要求有转速时所有应力严格为零。若大转角阶段出现伪应变，先修本构/约束客观性，不用大阻尼掩盖。

### F0：甲方搅拌项目，等温弹性 FSI 粗分辨率短跑

使用真实坩埚与弹性桨装配、开发常物性流体、已确认或明确标记的开发转速、平滑启动；关闭 EM、热源和浮力，先隔离机械耦合。

先跑 0.05–0.1 s；通过后按实测墙钟预算选择一个延长窗口。要求：

- 流体实际受到桨运动驱动，固体实际接收压力和黏性载荷；
- 固体内部积分次数和有效自由粒子数非零，夹持区确实仅占一部分；
- 无 NaN、无应变/Jacobian 异常，质量和推进时间有效；
- 检查几何穿透与数值近壁重叠的区别，输出最大越界量及数量；明确穿透容差的几何定义，不直接沿用任意 0.5dp 放行；
- 扭矩绕真实偏心轴计算，压力、黏性、驱动反力分项区分；驱动反力包含固体惯性，不能直接当作流体扭矩；
- 同一配置能重复运行，结果在相同精度及并行容差内一致。

可对照一次同装配、同流体条件的刚性桨短跑。若弹性偏移较小，应有相近整体流动响应；这是辅助一致性检查，不能强制要求两者完全相同，也不能用相似云图替代 E1。

### F1：恢复热/电磁路径的集成短测

F0 通过后才恢复现有稳定热源路径。先用常流体黏度和既有空间 Q，核对温度、入热和输出，再考虑温变黏度；不要同时开启未验证的局部高黏度、自感和 sigma(T) 反馈排错。

明确选择：真实 EM 计算得到的固定空间 Q、规定热源、固定温场诊断三者之一。只有前者能描述为“已计算 EM 热源”；固定 Q 也不等于 sigma(T) 动态反馈。

弹性体首轮不强制开启热应变。若实际启用固体导热/流固换热，须检查界面热交换与能量收支；若未启用，结果和界面明确说明，不伪造桨温度场。

## 7. 输出字段、文件与通过标签

固体至少输出：Position、Velocity、参考位置或可恢复参考构形、夹持区标记、实际可用应力张量、von Mises、弹性偏移及其定义。若现有实现不提供应变，明确字段不可用，不随意从速度替代。

von Mises 由物理 Cauchy 应力的偏应力计算；若求解器存储 Piola 应力，需用现有正确转换后再显示。不能把任意应力分量或应力范数命名为 von Mises。

驱动 CSV：物理时间、theta、omega、转轴、流体压力/黏性作用力和绕轴扭矩、可用驱动反力、固体 dt 与亚步数、最大应变/弹性偏移、固体求解耗时。

建议交付目录：

- `elastic_paddle_runs/S0_solid_reference/`
- `elastic_paddle_runs/T0_thermal_reference/`
- `elastic_paddle_runs/E0_hold/`
- `elastic_paddle_runs/E1_load_response/`
- `elastic_paddle_runs/E2_driven_rotation/`
- `elastic_paddle_runs/F0_isothermal_fsi/`
- `elastic_paddle_runs/F1_heated_integration/`

每项含 `effective-case.json`、`summary.json`、日志、原始 CSV、少量实际 VTP、复跑命令。事件/索引与原软件接口保持一致。

| 标签 | 授予条件 |
|---|---|
| SOLID_REFERENCE_PASS | S0 定量通过 |
| THERMAL_REFERENCE_PASS | T0 定量通过或同路径有效历史证据 |
| ELASTIC_PADDLE_ACTIVE | E0/E1 通过，实际本构和自由区积分接通 |
| DRIVEN_ELASTIC_ROTATION_SHORT_PASS | E2 短窗通过；整圈单列 |
| ELASTIC_FSI_SHORT_PASS | F0 真实双向短跑通过 |
| HEATED_ELASTIC_INTEGRATION_PASS | F1 通过，明确已开启的热/EM 范围 |
| LINUX_GUI_WORKFLOW_PASS | 安装包实际启动求解、显示并重开真实结果 |

严禁使用一个 `passed=true` 替代以上范围。粗分辨率功能通过不等于应力精度、扭矩精度、论文复现或合同全部通过。

## 8. 对原软件任务书的具体修改

| 原任务 | 本补充件后的执行方式 |
|---|---|
| M0 双平台 | 改为 Linux 基线及 Qt/VTK 试编；Windows 暂停 |
| M2 规定运动刚性桨 | 改为轴端驱动弹性桨；先通过 E0/E1/E2/F0，再接软件机械项目 |
| M1 自然对流 | 保留，可独立先接入；不得因无桨项目改动破坏回归 |
| M3 四环源 | 优先使用已有可用源，保证与 STL 元数据一致；新的有限截面算法不阻塞弹性桨和 GUI |
| M4 桌面工作流 | 沿用一个软件两个项目，加入材料及 mechanics 类型和真实固体字段 |
| M5 Linux 打包 | 近期发行主线，必须在目标机器实际试测 |
| M6 Windows | 本轮取消执行，保留后续规划 |

项目粒子缓存指纹增加：桨力学模式、固体分辨率、参考构形/预处理版本、夹持区选择。材料 E 改变通常不要求重生成几何粒子，但必须更新运行快照及固体初始算子，不得复用过时运行状态。

装配 STL 仅作初始预览。结果页用实际弹性粒子或其正确重建表面显示变形；不能把初始 STL 刚体动画当成弹性结果。显示放大变形时明确倍率。

机械项目缺少有效桨、材料、夹持区或可用固体后端时，校验阻止启动；不静默退回刚体。自然对流开关关闭搅拌时，按原任务书移除求解桨。

## 9. 本轮暂缓事项与遗留问题

- 暂缓 Windows、疲劳寿命、塑性、损伤、大规模细分辨率应力计算及钢桨 EM 涡流。
- 暂缓为对齐法国 Fig7/8 开新一轮大范围参数扫描。
- 显式局部 mu 在法国几何约 2.69 s 的尖峰仍未定位；隐式在另一入口越过该时刻，不作为所有路径已解决的证明。
- TEAM7 尚未完整通过，解析柱场精度亦需依据原始结果逐量签字。守恒、残差和总功率归一不替代场精度。
- 合同验证本轮按固体力学、传热分别整理；记录用户此口径，不自行扩写成疲劳或完整热应力验收。

## 10. 执行顺序、停止条件及最终交付

顺序：版本固定/能力盘点 → S0 与 E0/E1 → E2 → F0 → F1 → 机械项目软件接入 → Linux 包现场复测。T0 独立整理/补测；自然对流 GUI 和通用接口可在依赖满足时同步推进。

一项失败先给具体数值、首次异常时间、粒子/区域、实际步长及最小复现命令，再修一个主因后复测。禁止无限增开材料/阻尼/分辨率/转速扫描。

若遇到薄片不可解析、已有本构不能处理整圈转动、缺失真实流固载荷接口等阻塞，提交明确 BLOCKED 与最小修复范围。不得为了“赶完成”把全桨锁住、降 E、关闭载荷或仅放动画。

最终提交：

1. 可复建提交 SHA、依赖及构建命令；
2. 两套甲方项目配置和几何校验和；
3. 固体力学与传热独立验证结果/有效历史证据索引；
4. 弹性桨 E0/E1/E2、FSI F0、热集成 F1 的实际结果及标签；
5. Linux 软件包、启动方式、依赖、实际设备与已知限制；
6. 一页中文状态表，明确“弹性体已接通”“短跑通过”“整圈通过”“合同精度通过”各自状态。

**本轮最重要的判据：桨的自由区域确实在求解弹性力学、流固作用确实生效，然后由同一 Linux 软件调用真实项目完成计算与结果回放。**


<!-- ===== END 228_OPHELIE_LINUX_ELASTIC_PADDLE_SOFTWARE_ADDENDUM_2026-11-20.md ===== -->


<!-- ===== BEGIN 229_CLIENT_MATERIAL_AUDIT_AND_CURSOR_INPUT_RULES_2026-11-20.md ===== -->

# 229｜甲方回复后的物性、自然对流温度目标与双机开发补充规则

更新日期：2026-11-20  
原编号保留。本文件是第 228 轮及原软件计划的补充，不是新的独立开发轮次。

## 0. 资料与优先级

依据：甲方原文件 `实验中粘度和电导率(1).docx`、最新 `回复.docx`（包括两张内嵌图片），以及用户转述的补充消息。

**本版替代旧 229 的物性解释和待确认列表。** 与 228 或原软件任务书存在冲突时，按本版更新的甲方数据和排期执行。Linux 软件、真实弹性桨及独立固体/传热验证要求仍保留。

本版主要变化：

- 黏度原图错误已由甲方确认，公式为准，温区 1050–1800 K。
- 回复中的新图片明确电导率为自然指数，旧版按 log10 计算的甲方电导率及相关倍率不再作为当前配置依据。
- 转速 10 rpm，开发频率 100 kHz，30 kW 为电源输入功率。
- 当前定量重点为甲方自然对流热电偶温度约 1000℃；机械搅拌定性，约 1400℃仅为待补全时间与测点定义的实验趋势信息。
- 两机同时开展软件和物理验证；自然对流温度工作不等待弹性桨全部完成。

## 1. 已收到的甲方答复及配置状态

| 条目 | 甲方答复 | 当前处理 |
|---|---|---|
| 黏度 | 原图有问题，以公式和数值表为准，提供修正版图 | 以公式直接计算；修正版图趋势一致，旧图作废用于拟合 |
| 温度范围 | 1050 K–1800 K，无其他数据 | 作为甲方本轮指定使用窗口；不是已提供原始测点的声明 |
| 电导率底数 | 新图写 lnσ=4.5−8500/T，提醒使用 exp | 当前甲方公式按自然指数实现 |
| 电导率可信度 | 可能有问题，问题较大可用法国文献数值 | 允许有来源的法国替代预设；不伪称实测 |
| 转速 | 10 转/分钟 | 10 rpm，角速度约 1.04719755 rad/s；旋向仍待确认 |
| 频率 | 先设 100 kHz | 两个甲方项目使用 100000 Hz；保留“暂定”来源 |
| 功率 | 电源输入 30 kW | 记录 generator_input_power_W=30000；不是 absorbed_power_W=30000 |
| 电流、电压 | 千安量级，明天提供准确数字 | 未确认；不把 1000 A 当作已确认值；峰值/RMS、各环馈电也待确认 |
| 密度 | “2.7” | 暂按 2.7 g/cm³=2700 kg/m³解释，单位需补确认；不得输入成 2.7 kg/m³ |
| cp、k、beta | 参考法国文献 | 可先采用已追溯的法国表1常值，见第 4 节 |
| 自然对流温度 | 未开搅拌，热电偶处约 1000℃ | 1273.15 K；测点目标，不是全域平均值/最大值或初温 |
| 搅拌温度 | 开搅拌后最高可到 1400℃，搅拌定性即可 | 1673.15 K；峰值位置/时刻待确认，不强制做定量验收 |

## 2. 黏度与电导率：修正后的公式

内部温度统一为 K。

动态黏度（Pa·s）：

`mu(T)=exp(-4.0+8500.0/T)`

电导率（暂按 S/m）：

`sigma_client(T)=exp(4.5-8500.0/T)`

新图给出的等价十进对数为：

`log10(sigma) = (4.5-8500/T)/ln(10)`

**不能写成 `pow(10,4.5-8500/T)`。** 旧文件中该写法是基于旧曲线的暂时推断，现被甲方新图替代。旧版 1473 K 的 0.0536 S/m、1800 K 的 0.599 S/m 及相应倍率不得继续用作当前材料参数。

公式复算值（不是实验原始点）：

| T/K | 对应℃ | mu/(Pa·s) | sigma_client/(S/m) |
|---:|---:|---:|---:|
| 1050 | 776.85 | 60.05363594 | 0.02745415 |
| 1273.15 | 1000 | 14.53201178 | 0.11345444 |
| 1473 | 1199.85 | 5.87400286 | 0.28068105 |
| 1673.15 | 1400 | 2.94538010 | 0.55976520 |
| 1800 | 1526.85 | 2.05900369 | 0.80073740 |

两公式按数值计算满足 `mu*sigma=exp(0.5)`。这只能作为底数和实现一致性的快检，不能证明材料关联式物理正确。

动态/运动黏度、S/m 与 S/cm 的区别继续检查；电导率单位仍需补确认。T<=0、非有限值必须拒绝。排序和温度刷新后材料字段需正确同步。

窗口外禁止无标记外推。开发用端点保持必须明确开启并记录受影响数量/体积分数；冷区大范围越出窗口时不能称该材料已适用。室温弹性 FSI 功能测试继续用独立常黏度开发流体，不把熔体关联式外推到 300 K。

## 3. 法国电导率替代的选择规则

甲方已允许在电导率存在较大问题时采用法国文献数值，但未指定是哪一条。保持三个可区分的来源：

1. `client_ln_sigma_provisional`：本版修正的 exp 公式，甲方声明仍可能有问题。
2. `french_table1_sigma16_borrowed`：此前使用的法国表1常量 16 S/m。作为可追溯的工程替代候选，不声称是甲方实测。
3. `french_sigma_T_verified_source`：只有查清现有曲线对应哪篇文献、单位与温标后才启用；不能把表1常量16、约43的关联式值及图片随意拼接成新曲线。

建议先保留现有稳定的法国 sigma16 工程替代路径，另以修正甲方公式作一次同条件对照，再记录选定的软件开发预设。每次切换均输出 material_id、来源和使用值，不按“哪个更贴近1000℃”静默选材料。

给定电流模式必须重新求 EM 和原始吸收功率；固定吸收功率模式必须保留归一前功率、归一系数及实际流体入热。总功率归一不能验证局部 Q 和峰值。

本轮不强制同时增加 sigma(T) 反馈、自感新算法和新线圈有限截面算法；优先沿用已验证路径，增量测试。

## 4. 其他材料与目标温度的意义

密度暂定 2700 kg/m³；cp=1150 J/(kg·K)、k=4 W/(m·K)、beta=1e-5 K⁻¹可作为法国表1借用常值开发基线。与法国原 rho2750 区分；不得悄悄修改法国回归配置。

这些常值是已授权参考来源的工程借用，不代表甲方材料已测量，也不代表已启用全温区温变物性。温变 cp/k/beta 后续仅使用来源明确、实际算子确实读取的路径。

自然对流测点 1000℃时黏度约14.53 Pa·s，而不是历史常值4；1400℃时约2.95 Pa·s。同一几何下它们会改变流动响应，不能继续把全部熔体视为 mu4 并声称采用甲方物性。

测点从约1000℃升到约1400℃可能涉及热量输运/温场重新分布，也可能伴随功率和运行阶段变化；现有答复不足以确定原因。不得把搅拌本身当作凭空新增400℃的热源，也不得要求自然对流全域最高温必须为1000℃。

## 5. 功率口径：30 kW 不能直接替代历史吸收功率

配置必须区分：

- `generator_input_power_W=30000`：甲方提供。
- `coil_current_A`、`current_convention`、各环相位：待确认。
- `P_em_absorbed_raw_W`：给定电流计算所得熔体吸收功率。
- `absorbed_power_target_W`：如采用规定吸收功率开发模式，独立设定并标明假设。
- `coupling_efficiency`：如需估计，定义为熔体吸收/电源输入，不是现有测量值。

不能把历史法国50kW或搅拌60kW吸收预设直接复制到30kW输入的甲方项目。

若开发时用 P_abs=eta*30000，eta 在0至1之间，且明确说明是同一口径下对直接电磁入热的估计；上游功率测量、其他散热与其他热输入分别记账。没有其他能量来源时，不能让直接电磁吸收平均功率超过电源输入。

准确电流到达前，可以独立进行界面、弹性桨、规定空间热源和热边界测试；不得以随意试一个千安电流后归一功率，声称已得到真实电源预测。电流若用于敏感性计算，需少量、明确标记的取值，不扫描至“凑中温度”。

## 6. 当前验收重点与排期修正

### A. 自然对流：近期首要定量任务

目标是在甲方几何、测点和已确认实验条件下计算热电偶位置温度，与约1273.15 K对照。法国Fig7/8保留研究/模块回归价值，不再作为甲方新几何交接前必须完全一致的条件。

- 核对真实感温点，不能仅把热电偶STL末端当成已确认位置。
- 采用稳定的邻域插值/采样；输出采样坐标、覆盖、温度和物理时间，不用单个最近粒子编号长期代替固定测点。
- 输出该点温度时间序列、截面温度、全域平均/极值、热输入、各边界损失和能量账本。
- 约1000℃目前缺少测量时间、稳定窗口、误差条和初态。只能作为目标量级；不能仅凭某一时刻碰到1000℃授予定量PASS。
- 不把1273.15 K直接钉在测点当边界或设为全域初温后，用初始时刻声称达标。
- 若用目标温度校准未知等效换热或吸收效率，明确CALIBRATED，参数来源与允许范围固定；同一数据不能再称独立预测验证。
- 误差口径由双方固定。先报告摄氏度差/开尔文差（数值相同）及温度曲线误差，不自行选择绝对K分母或℃分母宣布合同10%通过。
- 短跑只验证功能和能量路径。长期/稳态对比按温度、能量与流动收敛判断，不能把几秒短跑直接当实际实验稳态。

### B. 搅拌：近期定性，弹性桨任务保留

甲方转速10 rpm：匀速一圈需6 s，含启动斜坡总时间更长。保留228的小时间窗冒烟，再按实测成本延长；不要提高转速冒充目标工况，也不要求现场短跑必须等一圈算完。

真实弹性桨、轴端驱动、力学短测与双向流固作用仍需实现。输出桨旁流动、热量输运和测点趋势；约1400℃暂不作为强制精度门槛。

不强制为了达到1400℃提高输入功率或替换材料。不因“定性即可”而用刚体动画冒充弹性响应。

### C. 两个项目的独立交接

自然对流不含桨，不应等待弹性桨完成才推进其温度测试或GUI接入。机械项目按228完成弹性/FSI门控后交接。两项目分别就绪、分别接入同一Linux软件，软件框架开发可同步开展。

## 7. Cursor 下一步有限执行清单

1. 替换旧甲方pow10电导率解释为新ln/exp公式，保留资料版本记录；不要修改法国材料ID。
2. 更新甲方预设：rho暂定2700、10rpm、100kHz、输入30kW、公式温区1050–1800K。分别标记确认、暂定、借用和未知，输出effective配置。
3. 对第2节各点做公式复算检查；1473K应为 mu≈5.874、sigma≈0.280681，浮点相对误差目标1e-4。检查用户输入℃只转换一次。
4. 先跑无桨自然对流的静态几何/测点与热账本检查；初态、冷却边界如未知采用明确开发配置，不将测量温度强制写为场约束。
5. 电流到达前，做一次法国sigma16与修正甲方sigma的同电流、同频率、同几何原始Q对照，前提是已有EM路径可用。保存未归一功率，不同时扫描多个物性/电流/边界。
6. 新实验温度主线先用可靠流体路径，常mu取实际测试参考温度下的公式值并明确“冻结mu”；局部mu(T)随后单独接入并检查实际步长/稳定性，不能宣称常mu已经是温变材料。
7. 弹性桨负责方继续228测试；温度目标主线不因该项等待。机械项目以10rpm作最终配置，暖态材料与室温力学测试配置分别保存。
8. 补全实验时刻与电流后，固定一套自然对流定量配置，再决定必要的物理时长或成熟稳态加速路线。不要先启动无边界依据的超长GPU计算。
9. 提交更新参数表、公式复算CSV、测点定义、功率口径、短测日志、通过范围及尚缺信息；无需重跑全部历史验证。

## 8. 本版仍需甲方补充的最少信息

1. 电流/电压准确值、峰值还是RMS、频率是否最终100kHz、线圈连接与各环电流关系；30kW实际测量位置/读数方式。
2. 自然对流1000℃是哪个时刻/是否稳定值，初温、加热持续时间、液深及准确感温点。
3. 搅拌1400℃是否同一热电偶峰值，开启搅拌前后的功率是否相同、转向及达到该值的时间。
4. 密度2.7的单位是否为g/cm³；电导率是否仍按S/m；冷却水/环境温度及能提供的热边界信息。

不再询问“黏度旧图还是公式为准”；该项已经答复。电导率底数按新图片采用ln，不再默认log10。

## 9. 两台电脑并行：Linux 软件与物理验证的分工及同步

本节为用户最新安排，替代此前“整体暂停物理验证”的排期建议：**软件机优先交付，物理机可继续验证；未经回归的新物理实现不自动进入软件审查版。** 第 228 号中弹性桨作为机械项目正式接入前置的要求仍然保留。

### 9.1 Git 对另一台电脑的实际影响

- 两台电脑各自 clone，拥有独立工作区、HEAD 和本地分支。分支名称相同不代表共享本地目录。
- A 机修改、commit：只改变 A 的本地仓库。push：更新远程分支，也不自动改 B 的源文件。
- B 执行 fetch：更新远程跟踪信息，不改变当前检出的源文件。
- B 执行 pull、merge、rebase、cherry-pick、switch/checkout 等，才可能改变当前源文件；B 的 IDE 若配置自动拉取或自动切换也要核查。自动 fetch 本身可以保留。
- 上述隔离前提是两个独立本地目录，不是共用网络盘、双向自动同步源码目录或两个编辑器连接同一个服务器工作区。
- Git 不会主动改已启动进程的内存，但重建/覆盖正在使用的程序、运行中读取的材料或输入文件仍会干扰计算。每次运行固定可执行文件、配置和输入副本。

### 9.2 推荐分支安排

从同一已记录且可构建的提交开始，继续使用同一个 OPHELIE 仓库，无需迁入 SPHinXsim 或主分支。

| 工作区 | 推荐本地/远程工作分支 | 职责 |
|---|---|---|
| 软件机 | `ophelie-software-linux` | Qt/VTK、CLI/JSON/事件、项目/装配、粒子和结果显示、Linux 打包；使用固定求解器版本 |
| 物理机 | `ophelie-physics-validation` | 本文件材料审核、场精度、TEAM7 等已有缺口、必要的物理修复与验证 |
| 原集成分支 | `feature/EM-Biot-Savart` | 集成经过验证的改动；记录审查版 SHA/标签 |

各工作分支在各自独立工作区提交。跨机器需要传递已提交代码时，再按项目已有授权的私有远程或 Git bundle 流程同步；本任务书不要求发布到公开仓库。

弹性桨只指定一个实现负责方。建议物理机承担核心弹性/FSI实现并尽早交付最小通过提交，软件机同时完成通用界面和自然对流项目；不要两台分别重写同一弹性桨入口。若已有软件机在做该实现，则沿用其所有权，物理机只协助测试，不重复开发。

采用不同工作分支后，双方仍可独立 commit 和 push。软件机只接收明确选定且通过测试的物理提交，不自动跟随物理分支最新状态。

### 9.3 如果坚持两台都使用同名集成分支

也可以独立本地 commit，但两边都向同一个远程分支 push 会出现分叉和非快进拒绝，需要整合；不能保证无冲突。不要用强推解决。

推荐同名集成分支只由一方维护，另一方有源码改动时建立工作分支再提交。仅运行验证不修改源码时可以保持固定 SHA，但不要把测试输入和必要脚本也全部丢弃为临时文件。

### 9.4 文件所有权和集成顺序

软件机主要修改 UI、配置解析/校验、适配层和打包目录；物理机主要修改算子、材料和验证项目。公共 CMake、共享头、机械项目入口、参数默认值和输出字段是潜在冲突区，修改前由任务记录明确所有者。

合并没有文本冲突也不等于兼容：物理代码可能改变材料默认值、字段单位、输出目录或数值调度。软件必须通过有版本的契约读取生效能力。

物理机交接每项修复至少包括：

1. 基于哪个 SHA、对应提交及所需前置提交；
2. 变动文件和是否改变配置/输出语义；
3. 一条最小复现命令和原始结果；
4. 实际设备、求解器构建配置和通过范围。

软件机在独立集成/试测目录合并或 cherry-pick（仅适用于可独立提取的提交），完成构建、受影响算例及 GUI 接口短回归后，才更新软件发行基线。不要在主构建或现场演示中途 pull 后继续沿用旧可执行文件。

两机独立 build、运行、预处理和安装目录；材料文件按版本随项目存放，不用网络共享文件实时覆盖。证据包由用户本机继续保存，新增结果记录 SHA 与配置校验和。

明后天审查版固定一个 SHA/版本号。物理机可继续推进，完成的新结果先进入物理分支；只有明确需要的修复经复测才纳入审查版。


## 10. 给甲方的简短后续确认文字

老师您好，收到。黏度我们按公式采用，温区1050–1800K；电导率按您新图中的自然指数exp实现，并保留法国文献参数作为可追溯替代。转速按10rpm、频率先按100kHz。30kW我们记为电源输入功率，不直接等同于熔体吸收功率。

为优先对齐自然对流测点约1000℃，还想确认：该温度是加热多久后的稳定值，还是过程中的某一时刻？初温、液深和热电偶实际感温点分别是什么？搅拌后1400℃是否为同一测点的峰值，开搅拌前后功率是否保持不变？密度2.7我们先按g/cm³换算为2700kg/m³，电导率单位先按S/m，请您确认。准确电流电压提供时，也请注明峰值/RMS和线圈连接方式。


<!-- ===== END 229_CLIENT_MATERIAL_AUDIT_AND_CURSOR_INPUT_RULES_2026-11-20.md ===== -->


<!-- ===== BEGIN 230_CLIENT_CASE_PHYSICS_LANDING_2026-11-20.md ===== -->

# 230｜甲方算例物理落地（材料 / 几何 / 预设）

日期：2026-11-20  
依据：Doc228、Doc229、软件任务书、甲方几何（已放入对应 relax `data/`）与回复。  
本轮**不包含** Linux GUI、Windows、弹性桨 E0–F1 完成声明。

## 1. 已落地代码

| 文件 | 作用 |
|------|------|
| `electromagnetic_ophelie_client_material_laws.h` | 甲方 `mu=exp(-4+8500/T)`、`sigma=exp(4.5-8500/T)`；冻结 μ；法国 σ=16 借用路径；公式复算 CSV |
| `electromagnetic_ophelie_client_geometry.h` | 两套 STL 角色、AABB、线圈四环、名义热电偶、CAD→求解器变换说明 |
| `electromagnetic_ophelie_client_case_params.h` | NC/搅拌 development 预设；30 kW **电源输入**账本；effective JSON |
| `test_3d_ophelie_client_case_audit` | 主机审计入口：公式门控 + STL 清单 + geometry/effective JSON |
| `OphelieClientNcGlassCylinderShape` / `OphelieClientNcCrucibleWallShape` | 甲方 NC **参数化**液腔与数值壁（不用壳 STL 填流体） |
| `test_3d_ophelie_client_nc_glass_relax` | 甲方 NC 液腔粒子生成 + SYCL-CK relax + Reload；写出测点/材料 meta |
| `applyClientNcMeltToFrenchReduced` / `applyClientNcDevelopmentToLiterature` | 把甲方熔体/材料映射进现有 French-reduced knobs |
| `frozen_q` CLI `--client-nc-geometry` / `--client-nc-material` / `--zero-q` skip-EM | 甲方 Reload 等温静水入口 |

## 2. 材料口径（已编码）

- **作废**：`pow(10, 4.5-8500/T)` 旧电导率读法及 1473 K≈0.0536 / 1800 K≈0.599 旧表。
- **现行甲方公式**：自然指数；1473 K → μ≈5.874、σ≈0.280681（相对误差门 1e-4）。
- **NC 开发默认**：`french_table1_sigma16_borrowed` + **冻结** μ(1000℃)≈14.532 Pa·s。
- **对照预设**：同几何 / 同频率的 `client_ln_sigma_provisional`（供后续同电流未归一功率对照；电流仍未知，本轮不跑假 kA）。
- ρ=2700 暂定；cp/k/β 借用法国表1并标记 borrowed。
- **不改**法国回归材料 ID / 法国算例默认值。

## 3. 几何口径（已编码）

- STL 单位 mm；CAD 竖直 +Y；推荐变换：坩埚外底轴心→原点后 `(x,y,z)→(x,-z,y)/1000`。
- 自然对流无桨；搅拌桨偏心历史值约 −67 mm（须复核后才硬编码驱动轴）。
- 热电偶尖端 = AABB 下端名义点，`status=unconfirmed_sensor_location`。
- 水冷管 / 线圈展示；零面积三角形保留原件并记入说明。

## 4. 功率口径（已编码）

- `generator_input_power_W=30000`（confirmed）。
- `coil_current_A` / peak-vs-RMS = **unknown**（不填 1000 A）。
- 禁止把法国 50/60 kW 吸收功率直接复制为甲方预设。

## 5. 如何复跑审计

已跑通一次，产物在：

`docs/ophelie/reference/client_case_audit_2026-11-20/`

```bash
# 配置并编译 test_3d_ophelie_client_case_audit 后：
./test_3d_ophelie_client_case_audit \
  --repo-root /path/to/SPHinXsys \
  --output-dir docs/ophelie/reference/client_case_audit_2026-11-20
```

期望：`client_formula_check.csv`、两套 `*_geometry.json` / `*_stl_manifest.csv`、三份 `effective_*.json`、`AUDIT_SUMMARY.md`；进程退出码 0。

## 6. 编译与运行（请本机执行）

### 6.1 材料/几何审计（无需 GPU 长跑）

```bash
cd /home/ubuntu/SPHinXsys/build
cmake --build . --target test_3d_ophelie_client_case_audit -j$(nproc)

./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_case_audit/bin/test_3d_ophelie_client_case_audit \
  --repo-root /home/ubuntu/SPHinXsys \
  --output-dir /home/ubuntu/SPHinXsys/docs/ophelie/reference/client_case_audit_2026-11-20
```

### 6.2 甲方 NC 液腔 relax（需 SYCL/GPU 构建）

```bash
cd /home/ubuntu/SPHinXsys/build
cmake --build . --target test_3d_ophelie_client_nc_glass_relax -j$(nproc)

cd tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_nc_glass_relax/bin
./test_3d_ophelie_client_nc_glass_relax \
  --dp=0.010 \
  --melt-height=0.180 \
  --relax-steps=1000 \
  --relax-vtp-every=100 \
  --output-meta-dir=./client_nc_relax_meta
```

可选：`--client-ln-sigma` 切换甲方 ln-σ 对照材料；默认法国 σ=16 借用 + 冻结 μ(1000℃)。

检查：`client_nc_relax_meta/{melt_geometry.json,effective-case.json,probe_nominal.csv,README.md}`，以及 Reload / VTP。

液深 `0.180 m` 为开发暂定，甲方确认后改 `--melt-height=`。

## 7. 等温静水短跑（Reload → frozen_q）

`test_3d_ophelie_french_natural_convection_frozen_q` 已接：

- `--client-nc-geometry`：甲方 NC 参数化液腔（R≈0.1415、H 默认 0.180、底 z≈0.001、4 匝线圈、100 kHz）
- `--client-nc-material`：冻结 μ(1000℃) + σ16 借用；关闭浮力；**不要**用法 50 kW 标定
- `--zero-q`：**跳过 EM Phase A**，`frozen_q=0`（等温静水）
- `--melt-height=`：与 relax 同口径（仅配合 `--client-nc-geometry`）

不要单独用 `--glass-height=`（会把 center 重置为 H/2，与甲方底面 z=0.001 冲突）；需要改高度用 `--melt-height=` 或同时设 `--glass-center=`。

```bash
cd /home/ubuntu/SPHinXsys/build
# 若刚改过源码，先 reconfigure 再编（未知 target 时）：cmake .
cmake --build . --target test_3d_ophelie_french_natural_convection_frozen_q -j$(nproc)

FROZEN_BIN=tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/bin
RELOAD=../test_3d_ophelie_client_nc_glass_relax/bin/reload   # 相对 FROZEN_BIN

cd "$FROZEN_BIN"
./test_3d_ophelie_french_natural_convection_frozen_q \
  --client-nc-geometry \
  --client-nc-material \
  --zero-q \
  --zero-buoyancy \
  --dp=0.010 \
  --melt-height=0.180 \
  --reload-dir="$RELOAD" \
  --end-time=0.50 \
  --hydro-settle-time=0.15 \
  --state-record-interval=0.05 \
  --summary-json=./client_nc_isothermal_hydro_summary.json
```

绝对路径也可：`--reload-dir=/home/ubuntu/SPHinXsys/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_nc_glass_relax/bin/reload`（目录内须有 `Reload.xml`）。

期望日志含：`[ophelie][client-nc] geometry`、`[ophelie][client-nc] material`、`[ophelie][stage4.1] zero_q: skip EM`；静水后 `escaped≈0`、`U_max ≪ c0`。

短冒烟可用 `--end-time=0.05`；正式静水验收用上面的 `0.50`。未显式 `--reload-dir=` 且带 `--client-nc-geometry` 时，会优先找 `../test_3d_ophelie_client_nc_glass_relax/bin/reload`。

**已知坑（已修）**：
1. `--zero-q` 跳过 Doc219 spatial-q coverage abort（否则 `n_uncovered=N`）。
2. `--zero-q` 在 settle 后补清 `DensityChangeRate/Force/Displacement`；Doc224 blowup 门改为 `0.05·c0`（约 0.66 m/s），避免静水残速 ~0.05 误杀。
3. `material_check passed=0` 在甲方冻结 μ 路径上是预期日志，不中止。

### 7.1 本机结果（2026-11-21，exit 0）

| 量 | 值 |
|----|----|
| `passed` / `smoke_passed` | 1 / 1 |
| `escaped` | 0/11232 |
| end `U_max` / `U_rms` | 0.0135 / 0.0038（settle 后衰减） |
| `P_joule` / `em_updates` | 0 / 0 |
| `boussinesq` | 0（`--zero-buoyancy`） |
| 备注 | Robin 冷却使 `dT` 缓增到 ~3.7 K（无扩散、无 Q）；等温静水门看 escaped + U 衰减，不钉绝对 T |

产物：`bin/client_nc_isothermal_hydro_summary.json`、`output/GlassBody_ite_*.vtp`。

## 8. 受控温差短冒烟（已跑）

冻结 μ≈14.5 下 Ra 被埋；本步只验：**温差场可施加 + 有/无浮力对照不炸**，**不**宣称可观测自然对流。

### 8.1 本机结果（2026-11-21，两跑均 exit 0）

| 跑 | boussinesq | dT | U_max | U_rms | escaped | passed |
|----|------------|-----|-------|-------|---------|--------|
| A `client_nc_dT20_boussinesq_summary.json` | 1（β=1e-5） | 19.82 | 0.01350 | 0.00381 | 0 | 1 |
| B `client_nc_dT20_zero_buoyancy_summary.json` | 0 | 19.82 | 0.01350 | 0.00379 | 0 | 1 |

结论：冻结径向 ΔT=20 K 稳定；A/B 速度几乎相同 → **在冻结 μ 开发预设下浮力信号不可分辨**（与 `buried_defer_finer_mesh` 一致）。等温静水 + 受控温差冒烟均 PASS。

## 9. Doc228 S0/T0：用仓库现成算例写报告（不自造通用力学 case）

用户口径（2026-11-21）：**纯模块独立验证**不必新写悬臂/导热 case；用 SPHinXsys 自带固体/流体/传热算例即可。自写的 `test_3d_ophelie_doc228_s0_cantilever` **作废**（且 `DistributingPointForces` 要壳 `Thickness`，体积梁会 `Thickness is nullptr`）。

S0/T0 报告证据（本机为 **SYCL 构建**：上游 `tests/3d_examples` 在 SYCL 下不编入；用 `extra_source` 里已有算例）：

```bash
cd /home/ubuntu/SPHinXsys/build

# S0 固体：现成 LG 被动悬臂（extra，可 SYCL 编）
cmake --build . --target test_3d_passive_cantilever_LG -j$(nproc)
cd tests/extra_source_and_tests/3d_examples/test_3d_passive_cantilever_LG/bin
./test_3d_passive_cantilever_LG

# T0 传热：OPHELIE 已有 1D 导热解析门（不是新造框架）
cd /home/ubuntu/SPHinXsys/build
cmake --build . --target test_3d_ophelie_thermal_1d_conduction -j$(nproc)
cd tests/extra_source_and_tests/3d_examples/test_3d_ophelie_thermal_1d_conduction/bin
./test_3d_ophelie_thermal_1d_conduction --summary-json=./t0_thermal_summary.json
```

非 SYCL 构建还可补跑上游 `test_3d_mr_cantilever_beam` / `test_3d_passive_cantilever`。  
不锈钢 304 预设头 `electromagnetic_ophelie_client_stainless_steel.h` 留给 **甲方桨 E0+**。

## 10. 甲方搅拌路线（对齐法国搅拌 case 结构）

用户口径（2026-11-21）：**按法国搅拌 case 装配甲方电磁+搅拌**；桨改为弹性力学；材料用开发说明书（Doc228/229）。不是另起一套通用力学玩具。

| 层 | 法国参考 | 甲方目标 |
|----|----------|----------|
| Relax | `french_stirring_glass_relax` → Glass/Rotor/Wall Reload | 同结构；熔体参数化液腔−桨；桨=甲方 STL |
| EM+流动 | `french_stirring_em`（host-spin 刚桨） | 同 Phase A/B；桨=`elastic_driven` + 304 |
| 材料 | 法国表1 / 60 kW | Doc229 冻结μ/σ预设；30 kW **电源输入**；电流 unknown |

头文件：`electromagnetic_ophelie_client_stirring_geometry.h`（脚手架）。  
E0（桨单独 hold）是接入全案前的 Doc228 固体隔离步，**不是**最终交付形态。

### 10.1 E0 桨 hold（源码树路径；先修 BoundingBox 后再跑）

```bash
SRC=/home/ubuntu/SPHinXsys
cd "$SRC/build"
cmake .
cmake --build . --target test_3d_ophelie_doc228_e0_paddle_hold -j$(nproc)

ST_RELAX=$SRC/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_glass_relax
# CAD + 预处理桨在 stirring relax/data/（CMake 会拷到 E0 bin/data）
python3 "$ST_RELAX/preprocess_client_paddle_stl.py"
BINE0=$SRC/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_doc228_e0_paddle_hold/bin

OUT=$SRC/docs/ophelie/reference/elastic_paddle_runs_2026-11-21/E0_hold
mkdir -p "$OUT"
cd "$BINE0"
./test_3d_ophelie_doc228_e0_paddle_hold \
  --stl=./data/client_stir_paddle_solver_m.stl \
  --dp=0.004 --solid-steps=100 \
  --out-dir="$OUT" --summary-json=summary.json
```

甲方 CAD：`client_nc_glass_relax/data/`、`client_stirring_glass_relax/data/`（已从 `docs/ophelie/software/*stl` 迁入）。  
其后：client stirring relax → F0 等温弹性 FSI。

### 10.2 E0 本机结果（2026-11-21，exit 0）

| 量 | 值 |
|----|-----|
| `PADDLE_HOLD_PASS` | 1 |
| n / n_clamp | 10119 / 649（夹持约 6.4%，非整桨） |
| E | 2e11 Pa（304，未降 E） |
| max_disp / vonMises | 0 / 0（无载等温 hold） |
| WARN | `dp=0.004` 叶片可能欠分辨；E1 应力前加密 |

产物：`docs/ophelie/reference/elastic_paddle_runs_2026-11-21/E0_hold/summary.json`

### 10.3 F0 等温弹性 FSI 短跑（编译 + 运行；本机 GPU 跑）

前置：`test_3d_ophelie_client_stirring_glass_relax --bodies=all` 已写出 Reload（`n_rotor≈4332`，`dp≈0.005`）。

```bash
SRC=/home/ubuntu/SPHinXsys
cd "$SRC/build"
cmake .
cmake --build . --target test_3d_ophelie_client_stirring_f0_elastic_fsi -j$(nproc)

RELAX_RELOAD=$SRC/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_glass_relax/bin/reload
BINF0=$SRC/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_f0_elastic_fsi/bin
# CMake 已从 stirring relax/data 拷贝 client_stir_paddle_solver_m.stl → bin/data
mkdir -p "$BINF0/output"

OUT=$SRC/docs/ophelie/reference/elastic_paddle_runs_2026-11-21/F0_isothermal_fsi
mkdir -p "$OUT"
cd "$BINF0"
./test_3d_ophelie_client_stirring_f0_elastic_fsi \
  --reload-dir="$RELAX_RELOAD" \
  --rotor-stl=./data/client_stir_paddle_solver_m.stl \
  --end-time=0.05 --omega-ramp=0.05 --c0=8 \
  --clamp-top-band=0.05 \
  --max-wall-hours=2 \
  --out-dir="$OUT" --summary-json=summary.json
```

要点：桨=`elastic_driven` + 304（E=2e11，不降 E）；关 EM / Q / Boussinesq；重力保留；夹持=桨顶 50 mm 暂定带；监控 CSV=`output/client_f0_elastic_fsi_monitor.csv`。不锈钢固体声速高 → 固体亚步多，墙钟可能较长。

### 10.4 F0 本机结果（2026-11-21，exit 0）

| 量 | 值 |
|----|-----|
| `ELASTIC_FSI_SHORT_PASS` | 1 |
| t / wall | 0.050 s / ≈41 s |
| n_glass / n_rotor / n_clamp | 43917 / 4332 / 278（夹持≠整桨） |
| E | 2e11 Pa（未降） |
| f_pres / f_visc / Tq | ≈1.79 N / 0.44 N / −0.012 Nm（双向载荷非零） |
| tip / e_off | ≈0.071 m/s / ≈0.22 mm |
| solid_steps | 70759 |

产物：`docs/ophelie/reference/elastic_paddle_runs_2026-11-21/F0_isothermal_fsi/summary.json`

### 10.5 F1 加热弹性 FSI（规定体积热源；非 EM）

Doc228：F0 后恢复稳定热源路径；本轮选 **规定体积热**（电流未知，不得标成“已算 EM”）。常黏度、无浮力、无桨热应变；热电偶 Observer 采样 T。

```bash
SRC=/home/ubuntu/SPHinXsys
cd "$SRC/build"
cmake .
cmake --build . --target test_3d_ophelie_client_stirring_f1_heated_elastic -j$(nproc)

RELAX_RELOAD=$SRC/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_glass_relax/bin/reload
BINF1=$SRC/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_f1_heated_elastic/bin
OUT=$SRC/docs/ophelie/reference/elastic_paddle_runs_2026-11-21/F1_heated_integration
mkdir -p "$OUT" "$BINF1/output"
cd "$BINF1"
./test_3d_ophelie_client_stirring_f1_heated_elastic \
  --reload-dir="$RELAX_RELOAD" \
  --rotor-stl=./data/client_stir_paddle_solver_m.stl \
  --end-time=0.2 --omega-ramp=0.05 --c0=8 \
  --heat-power=1e5 \
  --clamp-top-band=0.05 \
  --max-wall-hours=2 \
  --out-dir="$OUT" --summary-json=summary.json
```

默认 `P=1e5 W` 为短窗可见 ΔT 的开发热功率（非 30 kW 电源输入）。看 `HEATED_ELASTIC_INTEGRATION_PASS`、`energy_rel_err`、`probe_T_K`。

### 10.6 E1 非零载荷响应（无流体；不锈钢 E）

Doc228：夹持保留，对自由区施加两级已知总力（质量加权写入 `ForcePrior`），检查可分辨弹性偏移、应力更新、夹持静止及近似比例关系。不降低 E。

```bash
SRC=/home/ubuntu/SPHinXsys
cd "$SRC/build"
cmake .
cmake --build . --target test_3d_ophelie_doc228_e1_paddle_load -j$(nproc)

BINE1=$SRC/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_doc228_e1_paddle_load/bin
OUT=$SRC/docs/ophelie/reference/elastic_paddle_runs_2026-11-21/E1_load_response
mkdir -p "$OUT" "$BINE1/output"
cd "$BINE1"
./test_3d_ophelie_doc228_e1_paddle_load \
  --dp=0.004 --solid-steps=5000 --clamp-top-band=0.05 \
  --force1=1.5e4 --force2=3e4 \
  --out-dir="$OUT" --summary-json=summary.json
```

看 `LOAD_RESPONSE_PASS`、`ELASTIC_PADDLE_ACTIVE`、`e_off_ratio` / `ratio_rel_err`。粗 `dp` 下不宣称应力精度。

### 10.7 E1 本机结果（2026-11-21，exit 0）

| 字段 | 值 |
|------|-----|
| `LOAD_RESPONSE_PASS` | 1 |
| `ELASTIC_PADDLE_ACTIVE` | 1 |
| F1 / F2 | 1.5e4 / 3.0e4 N（+Y，自由区质量加权） |
| e_off_free 1/2 | ≈17.3 / 34.5 mm |
| e_off_ratio / force_ratio | ≈1.991 / 2（rel_err≈0.45%） |
| sigma/E 1/2 | ≈0.0043 / 0.0086 |
| clamp_e_off | 0 |
| wall_s | ≈48 |

产物：`docs/ophelie/reference/elastic_paddle_runs_2026-11-21/E1_load_response/summary.json`

### 10.8 E2 规定旋转短窗（无流体；10 rpm）

Doc228：夹持区平滑转速斜坡，先 0.05–0.1 s 冒烟。检查轴、符号、约束与固体亚步；弹性偏移相对刚体参考轨迹。整圈 2π 不在本短窗宣称。

```bash
SRC=/home/ubuntu/SPHinXsys
cd "$SRC/build"
cmake .
cmake --build . --target test_3d_ophelie_doc228_e2_paddle_spin -j$(nproc)

BINE2=$SRC/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_doc228_e2_paddle_spin/bin
OUT=$SRC/docs/ophelie/reference/elastic_paddle_runs_2026-11-21/E2_driven_rotation
mkdir -p "$OUT" "$BINE2/output"
cd "$BINE2"
./test_3d_ophelie_doc228_e2_paddle_spin \
  --dp=0.004 --end-time=0.05 --omega-ramp=0.05 --clamp-top-band=0.05 \
  --out-dir="$OUT" --summary-json=summary.json
```

看 `DRIVEN_ELASTIC_ROTATION_SHORT_PASS`、`theta_rad`、`sign_ok`、`max_clamp_err_m`。

### 10.9 E2 本机结果（2026-11-21，exit 0）

| 字段 | 值 |
|------|-----|
| `DRIVEN_ELASTIC_ROTATION_SHORT_PASS` | 1 |
| `full_revolution_pass` | 0（短窗不宣称） |
| t / θ / ω | 0.05 s / ≈0.026 rad / ≈1.047 rad/s（10 rpm） |
| max_e_off / clamp_err | ≈0.33 mm / 0 |
| tip_u / sign_ok | ≈0.068 m/s / 1 |
| solid_steps / wall_s | 73645 / ≈346 |

产物：`docs/ophelie/reference/elastic_paddle_runs_2026-11-21/E2_driven_rotation/summary.json`

## 11. 仍未做

1. 电流到位后：σ16 vs ln-σ 未归一 Q 对照；再换真实 EM Q。  
2. Doc228 E2 整圈（∫ω dt = 2π；本短窗不宣称）。  
3. 测点时间序列 + 1000℃ 窗口（不得钉测点温度冒充 PASS）。  
4. ~~两案软件交接包~~ → 见 `docs/ophelie/reference/client_software_handoff_2026-11-21/`（HANDOFF.md + COMPILE_AND_RUN.sh）。整圈 / EM Q / 1000℃ 定量仍开放。

## 12. 给软件机的契约提示

读取 `effective_*.json` 的 `material_id`、`power.generator_input_power_W`、`mechanics_mode`。  
自然对流短跑与搅拌 Doc228 门控（E0/E1/E2 短窗、F0、F1）本机已过；软件机可接两案预设，但不得把短窗标成整圈 PASS，也不得把规定体积热标成已算 EM。


<!-- ===== END 230_CLIENT_CASE_PHYSICS_LANDING_2026-11-20.md ===== -->


<!-- ===== BEGIN 231_OPHELIE_PHYSICS_NC_THERMAL_BALANCE_EM_AND_ELASTIC_PADDLE_CURSOR_PLAN_2026-11-21.md ===== -->

# 231｜物理机下一轮：甲方自然对流测温与热平衡、电磁缺口及弹性桨延长测试

日期：2026-11-21  
执行对象：原 OPHELIE 物理开发电脑上的 Cursor。  
软件机执行230，使用已交接版本。本轮物理修改不自动进入软件审查版。

## 0. 本轮目标及范围

第一优先：在甲方无桨几何上，建立可信的“实际入热—散热—储热—测点温度”链，判断约1000℃目标是否可达到，以及还缺哪些实验条件。

第二优先：从现有解析柱、TEAM7证据中明确电磁尚未通过的量，围绕一个具体问题推进。

第三优先：对已接通的10rpm不锈钢弹性桨延长测试，补真实旋转与流固稳定性证据。

禁止把法国Fig7/8完全对齐、Windows、重做GUI或迁入SPHinXsim作为本轮前置。不要重新实现已经交接的弹性桨，也不重新开展全框架改造。

本文件基于用户转述的交接摘要：NC静水/ΔT20、弹性桨E0–E2短窗/F0/F1已有首轮通过。**开始前必须核对实际交接包与源码；这些不是本任务书重新授予的PASS。**

本轮交付至少应形成：实际代码身份、输入状态表、无桨热账本及测点检查结果、有限暖态试跑结论、电磁缺口表、软件机可选择接收的修复提交。若条件不足，给明确BLOCKED/NOT_RUN及最小缺项，不用“还需继续长跑”代替分析。

## 1. 版本隔离与开工盘点 A

保留软件已接收的基线SHA/配置/binary哈希。物理机用独立工作分支，例如 `ophelie-physics-validation`，与软件机独立构建和运行目录。

先列明：

- 当前SHA、dirty状态、目标与源码入口、实际后端与精度；
- 两个甲方案例预处理与运行命令、effective及summary；
- 新弹性桨与旧刚性路径如何区分；
- 真实时间步约束、ForcePrior及材料字段绑定在哪些入口生效；
- 交接F1实际采用固定温场、规定热源、冻结EM空间Q还是在线EM。

保存 `A_BASELINE_AND_OPEN_ITEMS.md`。对照源码核查，不凭旧轮号猜测；若交接已解决旧缺口，登记证据而不重复开发。

用户本机已有历史round pack，不要求另建重复备份体系。验收表索引到具体包/文件/配置；新结果仍完整保存。

## 2. 输入条件：沿用最新版229，禁止旧值回流

| 输入 | 本轮依据与约束 |
|---|---|
| 几何 | 甲方NC STL及现有装配；不可复制法国R0.25/H0.185代替；液深沿实际已记录值并标来源 |
| 黏度 | mu=exp(-4+8500/T)，Pa·s，甲方指定1050–1800K |
| 电导率 | 新回复为sigma=exp(4.5−8500/T)，暂按S/m；不是旧pow10公式 |
| 法国sigma替代 | 允许借用；现有sigma16路径可作工程基线，来源独立，不伪称实测 |
| 密度 | 暂按2700kg/m³解释甲方“2.7”，单位待确认 |
| cp/k/beta | 可用已追溯法国表1常值1150/4/1e-5，标借用 |
| 供电 | 暂定100kHz；30kW是电源输入；电流/电压如有最新答复，先核对幅值/RMS、单位、各环关系 |
| 温度目标 | NC热电偶约1000℃=1273.15K，初态/时刻/稳态待确认 |
| 搅拌 | 10rpm，旋向按已确认值或明确假设；约1400℃为定性参考 |

快检：1273.15K下mu≈14.532Pa·s，sigma≈0.113454S/m；1473K下mu≈5.874、sigma≈0.280681。旧甲方sigma0.0536不是最新公式值。

运行前分别给出“已确认、甲方暂定、法国借用、开发假设、未知”。未知不填成看似真实的默认实验值。

## 3. B：先验证固定测点和热账本，再加长计算

### B1 固定热电偶采样

核对CAD到求解坐标变换、单位、液深和实际感温点。热电偶STL末端仅是候选位置，未经确认不能当实测点。若缺准确位置，暂用明确命名的开发测点；不在事后移动测点凑温度。

采用已有稳定插值在固定空间位置采样，记录邻居数、覆盖权重及是否落在流体有效域。不要长期跟踪某一个最近粒子编号作为固定测点。

做两个小检查，可采用现有后处理脚本，不强迫写新框架：

1. 均匀温场，测点应返回该值；
2. 已知线性温场，检查内部测点插值误差。开发目标：误差≤温度变化幅值的1%；界面附近若不满足，记录支持不足，先修采样而非改热源。

提供坐标表、与装配/粒子叠加剖面和 `probe_sampling_check.csv`。未确认物理感温点与插值算法通过分开记录。

### B2 能量收支：统计实际施加的量

对封闭玻璃域定义：

`DeltaH = E_Q - E_side - E_bottom - E_top_conv - E_top_rad + E_other + R_E`

各边界向外为正；反向热通量允许为负，不能截零。若包含机械耗散等其他能量，明确是否真的进入热方程，不能只由力学功率推断已加热。

H使用实际cp模式：常cp用同一参考态；cp(T)用焓积分。初场加载/重置后建立热账本基线；不能把初始化温度替换造成的DeltaH当时间推进能量误差。

- Q对实际流体体积积分，并与热推进中实际能量增量一致；
- 侧、底、顶面通量分项累计，对角区不重复计面积；
- 辐射使用K与正确T^4，参数来源明确；
- 冷却管显示、等效侧底边界不重复施加损失；
- 不额外扣除一项“维持玻璃温度能耗”；内部导热和对流是输运，不是全域额外散热；
- 固定温场诊断不能参与热平衡验收。

每个输出区间用累计能量差除以真实Delta t求平均功率，不能算不等时间间隔日志行的算术平均。

定义开发检查：`eps_E=abs(R_E)/max(abs(E_Q)+sum(abs(E_boundary_terms))+abs(E_other),1J)`。同时报告原始J、DeltaH及参考区间，防止小输入/大抵消时比例失真。

### B3 有限热模块检查

先检索同一算子/参数路径的既有证据。证据有效则引用，不重新跑全套。

若接入甲方新配置后缺必要证据，仅补：

- 均匀绝热、规定体热源：短窗平均升温与DeltaT=P*dt/(M*cp)一致；本测试关闭流动，检查能量入口；
- 各启用边界的静态温场瞬时功率：与独立面面积积分比较；侧/底/顶面对流/顶面辐射分项，不用初始诊断结果冒充长时验证。

开发目标：短窗能量eps_E≤1%，分项功率与独立积分差≤2%。超过时定位单位、面积、角区、施加/监控不一致或时间累计问题，先修复再跑耦合。阈值是本轮开发检查，不是合同温度误差定义。

## 4. C：无桨暖态热流动的有限推进

### C0 模式与起点

基于交接已稳定流体路径，只改变本项需要的开关。必须解除固定温场，让热扩散、Q、热边界与浮力实际生效。

起点可用均匀1273.15K作为**接近目标的数值初猜**，不可把它当已知实验初温；也可用交接已有暖态场。始终记录initial_guess，不用初始测点接近1000℃作为通过。

先采用冻结mu=14.532Pa·s及明确来源的常cp/k；标记冻结物性。浮力参考温度按现有一致方案固定，不随时重置绝对温标或人为消除真实温差。

### C1 热源的两种条件，不能混用

**准确电流已知且EM路径可靠时：**使用确认电流/频率及选定材料重算原始Q。记录原始EM功率、映射功率与实际施加功率，检查给定电流模式未被目标功率归一悄悄覆盖。

**电流仍未知或EM精度未闭合时：**使用单独命名的规定空间热源开发模式，明确Q来源。若借用EM形状，注明生成条件不同，不当作已确认实验Q。

不得自动把30kW电源输入当熔体吸收功率，也不得复制法国50kW/旧搅拌60kW。若缺已记录P_abs基线，用`P_abs=eta*30000`及显式声明的eta作诊断；选择依据与未知状态保留，不能为“看起来像实验”隐去假设。若无选择依据，先以B项单位功率归一检查推进，实验温度标签保持BLOCKED_INPUT，不无休止等待也不伪造物理输入。

### C2 时间窗口与运行预算

先0.1–1s功能检查；通过后5s；再决定是否到20s。已有同配置有效证据可跳过重复短窗。每次延长必须在报告说明将解决什么不确定性。

记录物理秒/墙钟秒、最小dt、限制来源、Urms/Umax、测点温度、温度范围、mu使用值、越界/数值异常。先按实测性能估算加长成本，再启动。

20s不是稳态证明。若短窗已经显著温漂，先用功率账本解释方向，不盲目加长。总功率净差给出近似热容温漂预测，与实际平均温漂比较。

### C3 稳态/统计稳定判据

输出三个相互独立的判断：全域热量、固定空间温场、流动状态。

开发用建议门限（固定在运行前，不能看到结果后放宽）：

- 窗口平均`abs(P_in-P_out)/max(abs(P_in),abs(P_out),1W)≤2%`；
- 窗口内独立能量闭合eps_E≤1%；
- 测点及体积平均温度趋势绝对值≤0.01K/s，窗口首末温度差≤1K；
- 固定空间网格上温度场变化的归一L2≤1%，分母采用约定温差尺度并给下限，不能用绝对温度掩盖温差误差；
- Urms和环流/截面统计不再持续漂移，允许真实波动但应报告统计量。

至少连续三个非重叠窗口通过。初始窗口可用10s作诊断；若热/流动特征时间更长，必须延长窗口并说明依据。仅初始几秒满足上述小漂移，仍不足以证明长期稳态。

若初猜恰好接近热平衡，条件允许时从相差50K的另一暖态初猜作有限交叉检查，验证向相同候选场收敛；不要求本轮强行完成两条超长计算。未完成则标候选稳态/尚未检查初态依赖。

短跑完成、数值稳定、热平衡、实验对齐四个状态分别记录。

### C4 不长跑也要给出有用判断

稳态未达时，至少说明：净功率正/负、温漂是否与热容匹配、主要散热项、当前瓶颈是物理时间不足还是数值不稳定/输入不全。

若成本过大，可提出复用已有成熟稳态热或分区耦合方法；新伪时间加速若改cp/惯性等，必须另作算法验证并以原物理方程复核，不能把加速时间当真实实验时间。本轮不默认开发新稳态求解器。

## 5. D：物性与边界的一次受控升级，不做无目的调参

C路径稳定后才接局部mu(T)，保持其他输入不变。检查每粒子读值、温度刷新/排序、真实黏性步长、显式/隐式路径及截断比例。

1050–1800K内mu约60.05至2.059Pa·s，不直接照搬法国cap400配置。窗口外的端点保持或其他政策必须显式记录，不把截断当真实冷区模型。

出现此前类似局部尖峰时，导出首次异常粒子的坐标、T、mu、密度/压力、邻域支持及各力项；不只保留Umax超过门槛。先用同配置常mu对照隔离，不立即扫多组cap和阻尼。

热边界按甲方水冷条件建立。法国h只能作借用假设，不能因主线不匹配随意调所有面。需要校准时，先固定输入功率口径，只动一个有物理意义的未知参数；最多一个基线加两个有界对照，依据灵敏度决定下一步。

单个测点不足以同时辨识吸收效率、侧壁h、顶面发射率和sigma。用1000℃校准的结果标CALIBRATED；独立预测验收需另有数据。合同误差分母/窗口由双方确认，本轮不自行宣布10%通过。

## 6. E：电磁验证的有限任务

### E0 先清账，不沿用模糊PASS

整理现有解析柱、TEAM7的原始参考、配置、场误差和未通过项。最近交接可能已有新结果，以实际证据为准。

表中分开：源场、复电流、磁场、边界法向电流、守恒、积分功率、局部Q/峰值、迭代残差、分辨率与运行后端。内部残差/功率归一不是外部场精度验证。

### E1 解析柱优先

核对参考解对应无限/有限长、外场方向、频率、sigma、自感和边界，不拿不一致模型强比。沿已有门控找出最主要失败量，仅修一个原因并用原指标复测。

源场幅值、相位/RMS约定必须先一致。保留未归一Q与损耗。解析解附近零值用有物理尺度的范数，不用逐点相对误差除零。

### E2 TEAM7只攻一个已有明确缺口

先选择当前失败最集中的一项，如边界邻域/法向电流/电流闭合，附失败位置图与误差分布。构造最小对照区分边界、源场和自感迭代，不同时更换所有phi算法。

修复后同输入同指标比较；若模块误差目标为≤10%，预先说明比较的量、位置、范数、频率与采样，不删去不利点再授PASS。未达标继续FAIL/PARTIAL，不阻塞软件现有功能。

低电导玻璃测试不能替代高导铝板TEAM7；TEAM7未过也不自动证明全部玻璃加热不可用，需按模型范围和误差类型判断。

### E3 甲方电流到达后

固定100kHz及实际四环数据，确认各环电流幅相、峰值/RMS。比较修正甲方sigma与已选法国替代的一次同条件结果，输出原始P_abs和Q；不能各自归一到相同功率后只比较总功率。

如P_abs超过同口径电源输入，或非常小，先排查安匝/电流分配、单位、幅相、重复积分和模型适用性，不直接归一消掉问题。

## 7. F：10rpm弹性桨延长测试

先沿交接已通过的等温F0延长，暂不叠加新的mu(T)、EM反馈或热应变。保持真实不锈钢E和密度、轴端驱动与自由区积分；不得为速度降E或锁死整桨。

顺序：已有短窗 → 0.5s → 根据性能决定一次1–2s或整圈。达到整圈必须满足积分theta≥2pi；10rpm含斜坡时不等同固定6s。记录actual_angle、full_revolution_pass及最大弹性偏移/应变/应力、Jacobian、实际solid dt和亚步数。

检查压力/黏性载荷与驱动反力区别、流体实际响应和桨壁间隙。固体真实后端CPU或混合如实说明，不强制本轮重写GPU固体。

稳定后复用原F1热集成作一次短测。搅拌1400℃暂不设硬性门槛，不能为达标人为追加热源。

## 8. 优先级、预算及停止条件

| 顺序 | 执行内容 | 转下一项条件 |
|---|---|---|
| P0 | A版本/输入，B测点与热账本 | 数据可追溯、实际能量链可解释 |
| P1 | C无桨暖态有限试跑 | 可给出温漂/热平衡判断，不强制本轮达到最终稳态 |
| P2 | D局部mu接入或一个边界敏感性 | 仅在C有必要且输入足够时做；不同时调两类因素 |
| P3 | E电磁缺口表与一项修复 | 有具体前后误差证据；失败也说明原因 |
| P4 | F弹性桨延长 | 资源允许，优先保留可靠性；不阻塞P1 |

同机可用空闲GPU独立跑已固定的测试，但明确设备、显存及输出目录，避免互抢影响判断；不要因有多GPU就无控制发起扫描。

任何阶段：NaN/越界步长/能量异常/明显非物理失稳即停止加长。保存首次异常状态。输入未知则标BLOCKED_INPUT，继续独立基础检查。超预算先汇报成本与将解决的问题，不用无限长跑代替定位。

## 9. 目录、状态与交付

建议 `doc231_runs/` 下按A/B/C/D/E/F分组；大结果不必进Git，但必要配置、分析脚本、参考索引和任务书应可追踪。

必须产出：

- `ROUND231_REVIEW.md`：已解决的问题、证据、未解决原因、下一最高优先项；
- `ACCEPTANCE_TABLE.csv`：PASS/FAIL/BLOCKED/NOT_RUN及配置/文件索引；
- `INPUT_STATUS.csv`：确认、借用、假设、未知；
- `thermal_ledger.csv`、`probes.csv`、关键截面图和有限原始场；
- `EM_VALIDATION_GAPS.md`及本轮选定缺口的修复前后结果；
- 弹性桨测试如执行，保留角度、固体步长及真实力学状态；
- `SOFTWARE_HANDOFF_DELTA.md`：是否建议软件机接收、提交SHA/前置依赖、接口/参数变化和短回归命令。

标签至少区分：PROBE_SAMPLING_PASS、THERMAL_LEDGER_PASS、NC_WARM_SHORT_PASS、NC_THERMAL_BALANCE_CANDIDATE、NC_STEADY_STATE_CHECKED、NC_EXPERIMENT_COMPARISON、EM_FIELD_VALIDATION、ELASTIC_ROTATION_EXTENDED。

未提供实验时刻/稳态/测点等条件时，NC_EXPERIMENT_COMPARISON不得标正式PASS。软件机只接收确有必要且复测通过的修复，不自动更新整条物理分支。

## 10. 本轮必须回答的六个问题

1. 测点取值是否可信，位置是否已由甲方确认？
2. 实际吸收/施加热量是多少，30kW输入是否被错误当成吸收功率？
3. 当前温度升降能否由实际能量收支解释，主要损失来自哪里？
4. 距热平衡是物理时间不足、边界/功率未知，还是数值实现有错？
5. 电磁剩余失败具体在哪个量/区域，本轮是否减少了该误差？
6. 哪些修复可交给软件机，哪些仅限诊断，哪些仍未通过？

**本轮以能解释的测温和能量证据推进交付，不以云图像不像、短窗没崩或任意一时刻碰到1000℃作为完成。**


<!-- ===== END 231_OPHELIE_PHYSICS_NC_THERMAL_BALANCE_EM_AND_ELASTIC_PADDLE_CURSOR_PLAN_2026-11-21.md ===== -->


<!-- ===== BEGIN 232_OPHELIE_PRODUCTION_THERMAL_PROBES_BOUNDARY_AUDIT_AND_EM_REPAIR_PLAN_2026-11-22.md ===== -->

# 232｜物理机研究计划：生产测点、真实热算子、冷却层演化及电磁单项修复

日期：2026-11-22  
执行对象：OPHELIE物理开发机上的Cursor。  
基线：`client_physics_handoff_doc231_2026-11-22.tar.gz`及其源码/配置；软件机仍执行230，固定已经接收的版本。

## 0. 本轮要解决的问题与结论边界

本轮不换计算框架，不重新实现弹性桨，不用调参代替验证。先把231的独立检查接到实际生产路径，再判断是否需要延长计算或调整边界。

231已经证明：开发均匀10kW热源能驱动温度演化；本组初态1秒内平均降温约1.36K，与累计热收支量级一致。

231尚未证明：

1. 生产求解中固定热电偶的温度曲线可信——交付的`probes.csv`是独立均匀/线性场测试，不是暖态时间序列；
2. 生产GPU热源算子通过绝热验证——B3使用自写CPU加热循环；
3. 实际散热边界的面积、温度采样与通量正确——分项独立检查未完成；
4. 当前启动运动已可解释为自然对流——存在明显初始声学/静水残余运动；
5. 测点1000℃稳态不可达到——57kW为近均匀1000℃初场下的瞬态散热，不是已知最终稳态损失；
6. TEAM7已通过——本轮仅历史清账，无新修复证据。

**本轮最小成功成果：真实测点时间序列 + 生产热源检查 + 生产边界独立核对 + 一组受控暖态延长/浮力对照 + 一项电磁误差的前后证据。**

物理精度、数值稳定、短跑完成和实验对齐继续分开记录。软件机不自动拉取本轮代码。

## 1. A｜版本与231结论更正

### A1 固定身份

记录实际SHA、dirty文件、补丁/新文件、编译目标、实际GPU/CPU路径及精度。231包的基线为`832dbc782...`，但部分修改尚未提交，不能只检出该SHA就声称复现了231。

保留源码、配置和小型必要输入；大型Reload/结果可单独保存并记录校验和。禁止覆盖用户旧成果或使用整个源文件复制破坏软件机修改。

### A2 修订报告口径

提交`A_ERRATA_231.md`，明确：

- B1仅为独立插值小测试通过，生产测点尚未接通；
- B3为CPU显式加热与账本自洽检查，不是生产GPU热算子已验证；
- 热账本闭合不等于边界模型/离散正确；
- 57kW瞬态散热不足以推出内部测点1000℃稳态不可能；
- 甲方电流未知只阻塞甲方给定电流预测，不阻塞有指定激励的标准电磁验证。

保留原始结果，不篡改旧日志。新报告引用本更正，不延续过宽PASS标签。

## 2. 全轮固定输入与有限变更原则

主热流动试验继续使用231甲方几何：R约0.1415m，开发液深0.180m，dp0.010m，以实际Reload元数据为准。液深和热电偶感温点未确认状态保留。

固定暖态初猜1273.15K，rho暂2700、cp1150、k4、beta1e-5；冻结mu约14.532Pa·s。它们的来源沿229，不修改法国回归配置。

规定10kW均匀吸收热只作为开发热源，**不是真实EM，不是甲方确认吸收功率**。30kW仍为电源输入。甲方电流若新到达，单独版本接收，不中途替换本轮对照输入。

热边界先保留231值：侧300、底35、顶面对流20 W/(m²·K)，发射率0.8，参考环境/水温300K。全部标借用/假设。本轮先核对实现，不先调h。

暂不叠加局部mu(T)、sigma(T)反馈、自感新算法、skull新模型、Marangoni或新壁面方案。各独立制造场测试可采用专门温场，但必须复用生产算子和几何，测试后不污染正式初态。

## 3. B｜把固定测点接进真实生产运行（最高优先）

### B1 生产接入

复用231独立采样函数或现有成熟Observer，提取为共享调用路径；生产和验证调用同一采样实现。若改用另一种插值，不可直接沿用旧B1的PASS。

固定空间坐标采样当前粒子位置与温度，随流体运动和排序正确刷新。不能追踪固定粒子ID代替固定空间测点。

至少包含：名义热电偶点、轴线中部、内部半径点及靠壁诊断点。热电偶候选约(0,0,0.130)m，位置仍标`unconfirmed_sensor_location`；坐标变换/基准按实际CAD确认。不得为了贴近目标温度事后移动。

输出真实`probes_timeseries.csv`：time_s、probe_id、solver坐标、T_K、T_C、有效邻居数、支持半径、权重和、覆盖/无效状态、配置/位置确认标记。

区分界面最近粒子数据、带宽平均和严格表面温度。支持不足时输出null/invalid及原因，不回填目标温度或上一时刻伪造有效数据。

采样应在明确的完整热/流体更新阶段执行，文件记录真实时间。无需每个声学子步强制同步；使用固定采样间隔即可。

### B2 同路径验证

均匀场、内部线性温场各测一次，直接调用生产采样入口。均匀误差目标≤0.001K，线性场误差≤场变化幅值1%；记录实际浮点容差和覆盖。

再在真实暖态短跑中保存一个时刻的当前粒子Position/Temperature，由离线脚本独立复算至少一个测点，与在线结果比较。确认不读初始温度或过期host数据。

**通过条件：** 同路径制造场与在线/离线一致性通过，并产生至少两个不同非零物理时刻的有效样本。物理感温点未确认不阻塞算法通过，但实验标签不授PASS。

## 4. C｜真实热源算子的绝热验证

### C1 必须复用生产路径

在甲方NC入口加入明确诊断配置，或复用已有入口：关闭热边界、流体运动和EM，开启与正式运行相同的热源算子、温度同步、热子步调度和字段存储。

不能再调用独立的`hostApplyUniformHeat`作为生产热算子验证。若正式路径本来就是CPU，必须如实说明；不以程序编译启用SYCL等同热算子在GPU执行。

优先复用已有源接口，新增边界关闭开关时限定诊断范围，默认工况不受影响。

### C2 独立参考

参考质量来自真实粒子Mass求和，参考cp来自实际材料配置，规定功率使用用户输入P；不能只拿源内同一个归一因子回算作为独立参考。

常cp绝热参考：`DeltaT_exact=P*t/(M*cp)`。

独立计算：`DeltaH=sum(m_i*cp*DeltaT_i)`，用增量温度避免约1273K绝对温标的float32消减；同时检查可视化Temperature与参考温度+DeltaT一致。

先10kW、0.2s；必要时一次0.4s或半热步长检查，不无限扫描。验证实际P、dt累计、终点截断和热源是否重复施加。

**开发门限：** 质量加权平均升温相对误差≤1%，累计能量误差≤1%，时间与实际执行一致、无重复入热。低输入导致近零量时同时输出绝对误差，不用不合理分母授PASS。

## 5. D｜生产热边界的分项独立核对

### D1 面积、分区与符号

使用真实NC粒子和实际边界算子。输出各面作用粒子、离壁距离、所分配面积、角区归属及总面积；侧/底/顶不要重复计算。

解析圆柱参考：A_side=2*pi*R*H；A_bottom=A_top=pi*R²。参考R/H来自本轮实际液域，不从旧法国几何取值。

面积归一后总面积相等只是检查之一。还需看面积分布和温度取样：20mm厚边界带的粒子温度不自动等于真实壁面温度，须说明采用的近似。

每个边界热流以向外为正；若T低于环境应可为负。对流和辐射不应无理由截零。辐射使用K与正确T⁴。

### D2 同状态两个温场

不先做长时计算。选择：

1. 均匀1273.15K；
2. 一个明确可复算的非均匀温场，例如仅沿z线性变化，温区选1100–1300K。

依次只开启侧、底、顶面对流、顶面辐射；计算生产算子实际能量变化或实际热率。对显式算子可用单步冻结初场诊断，若使用隐式边界，按其实际离散时刻匹配参考，不混用步前/步后值。

独立脚本对解析表面作独立积分或解析计算，不复用生产粒子面积权重生成“独立参考”。

**开发门限：** 各分项总功率差≤2%，实际热方程变化与对应边界能量记账差≤1%；失败时同时给原始面积、功率及误差位置。总面积正确但非均匀温场失败，应先定位离壁取样和权重，不改h凑总功率。

### D3 通过与失败后的动作

通过：只能说明该分辨率下这组边界离散检查通过，不能说明借用h是甲方真实h。

失败：先修一个直接原因，再跑同一对照。若仅靠壁温度分辨率不足，记录最小必要分辨率对照方案，不立即重建框架或开启多参数扫描。

完成B/C/D后，才恢复有流动、有入热和全部边界的主测试。

## 6. E｜启动扰动与浮力的匹配对照

231速度中止门从0.05m/s改为0.05*c0，约0.664m/s；日志早期Umax约0.188m/s，不能只因没触新门而称浮力准确。

### E1 两组主对照

使用同一初态、粒子、压力/重力路径、热源、边界、步长和输出：

- E_ON：beta=1e-5；
- E_OFF：只关闭Boussinesq增量，保留真实重力和静水支撑。

不得把E_OFF写成关闭所有重力，也不要同时换平顶/自由面/压力方案。

两组先1s；如果仍主要是启动衰减，可将有效组延长至后续F窗口，不另开大量分辨率测试。

输出Umax/Urms、动能、密度偏差、液面变化、测点温度、热功率及截面速度，使用固定时间窗口比较。信号与噪声无法区分时标`BUOYANCY_SIGNAL_UNRESOLVED`，不凭流线或颜色授PASS。

### E2 实际步长检查

日志显示231主跑`doc227_actual_timestep=0`。本轮先输出真实声学dt、外层实际跨度、黏性力刷新年龄、热子步dt及其限制；常mu不代表可以省略约束检查。

若发现跨预算，按已有231之前的成熟修复接入，所有匹配对照统一使用修后的调度；不要同时改压力和边界。未越界则保留当前路径并登记证据，不只为打开某个轮号开关而改代码。

数值硬中止和物理可信度分开：保留非有限值、实际步长越界、显著粒子失效及Ma检查。速度超过0.05m/s可作诊断警告，不自动等同崩溃；低于0.05*c0也不自动等同可信。

## 7. F｜冷却层演化：一组有限延长，不先调边界

### F1 固定配置

在前述生产验证通过后，用E_ON配置从相同初态跑到5s；可用1s结果续算仅当restart包含所有必要状态且已验证，否则重跑到5s。

保持10kW开发源、原h和冻结mu；输出每0.1–0.5s测点/热账本，每1s左右少量截面。具体频率按实际成本选择，记录真实保存时间。

必须输出：

- 名义热电偶、中心及靠壁温度曲线；
- 各面边界带面积加权温度，注明它是离散取样温度；
- z中面径向温度剖面、轴向剖面与冷区变化；
- 各边界损失功率及累计能量；
- 全域质量加权温度、极值、动能/速度、粒子有效性。

### F2 5s后只作一次延长决策

若边界损失开始下降、冷区形成且数值/能量检查正常，说明初始57kW不是稳态常量，可进一步到20s观察；不要求本轮20s内稳态。

若散热几乎不变，检查墙带温度是否实际更新、热扩散/输运是否正确、时间尺度与边界带厚度，不能仅凭“不变”宣布物理错误。

若出现非物理失稳/步长越界/热账本失败，停止加长，导出首次异常状态。若5s成本已不合理，报告性能与主要同步/亚步成本，提供最小改进，不盲跑。

### F3 最终判断

本轮允许结论为“冷却层演化中，尚未稳态”。候选热平衡需输入/输出功率近似平衡、固定空间温度与流动统计不再持续漂移；不能只看总能量或一个点。

1000℃是测点目标，不是必须维持整个熔体均匀1000℃。禁止用初场瞬时损失比较电源功率后证明目标不可能；也禁止为了达标提高P_abs超过30kW输入或任意降低h。

窗口外mu(T)、skull或真实电磁形状会影响以后稳态，但本轮先完成冻结mu/开发源下的可信演化，不同时引入。

## 8. G｜电磁独立支线：本轮必须推进一项，而非只登记

优先复现现有解析柱中部复电流场约5%–6%误差。使用标准算例自带激励，甲方电流未知不阻塞。

1. 找到真实原始参考、条件、指标与运行命令；无限/有限长、频率、sigma、自感及幅值/RMS一致。
2. 固定一个基线复跑，输出复数J分量/幅相、区域定义、误差及未归一损耗。
3. 依据空间误差形态与残差，选一个主假设：源输入、离散/边界、场重构或自感收敛。最多一个受控修复/对照，不同时更换所有phi方案。
4. 用同一参考和指标比较前后值；无改善就报告假设被否定，不另造新的PASS口径。

若解析柱所需输入在当前工作树缺失，先恢复已知输入。确有外部阻塞可选TEAM7最小边界/诱导场对照，但必须说明原因、运行证据与具体变化；禁止再次只交历史缺口表。

没有修复成功也可以形成有效研究结果：明确误差不随某因素变化、排除一个原因。但不能称精度缺口已关闭。

独立参考场误差、边残差、守恒和功率归一继续分开。解析柱某区域低于10%不等于TEAM7或局部Q峰值整体通过。

## 9. H｜软件可接收代码的必要防护

规定热源接口可作为开发功能，但接收前必须：

- 用显式“参数已出现”标记处理模式，不以P>0唯一判定来源；
- 拒绝负值、NaN、Inf和非法字符串；零功率支持与否明确，不允许输入无效后自动退回未知电流EM；
- 与zero-q、固定温场、给定电流/目标功率等冲突在完整解析后统一检查，交换命令顺序仍给相同结果；
- 规定源模式中EM状态为NOT_RUN/disabled，不能用phi_residual=0或em_solver_ok=1伪装EM收敛；
- 日志/summary/effective统一显示规定入热、实际入热、generator30kW及实际源模式，清理遗留target50kW等误导性最终参数输出；
- 规定源与在线EM更新互斥，除非存在已定义且验证过的组合；禁止运行中又被旧EM刷新覆盖。

最小软件回归：合法10kW、非法值、冲突参数两种顺序、源模式输出一致性。不要将放宽速度中止门独立包装成“物理修复”。

仅交小提交/补丁及必要前置依赖，不整文件覆盖软件工作区。物理新诊断默认不改变软件示例。

## 10. 本轮不作为前置的事项

- 弹性桨整圈：本轮可以暂缓。先解决无桨温度主线；上一轮0.05s约346s成本已记录，不无条件启动长跑。
- 全部法国Fig7/8对齐、Windows、GUI、疲劳、钢桨EM涡流、全温变强耦合继续留在其他任务。
- 不为已有有效短测重复造验证；本轮新增检查必须对应已识别的生产接入缺口。

## 11. 验收标签及停止条件

| 标签 | 条件 |
|---|---|
| PRODUCTION_PROBE_CONNECTED | B生产时间序列及在线/离线检查通过；位置确认独立 |
| PRODUCTION_HEAT_SOURCE_VERIFIED | C使用真实生产算子通过，不是CPU替代式 |
| BOUNDARY_FLUX_AUDITED | D四类边界的均匀/非均匀积分通过 |
| ACTUAL_SCHEDULE_AUDITED | E实际步长与刷新时间检查完成，无违规或修复后通过 |
| BUOYANCY_RESPONSE_STATUS | 记录可分辨/未分辨/异常及ON-OFF数据，不强制授PASS |
| COOLING_LAYER_EVOLUTION_OBSERVED | F真实壁带/中心/测点与分项散热时间证据 |
| NC_STEADY_STATE | 本轮默认OPEN，只有实际充分收敛证据才更新 |
| NC_EXPERIMENT_MATCH | 缺实验位置/时间/吸收条件时保持BLOCKED_INPUT |
| EM_ONE_HYPOTHESIS_TESTED | G真实基线与一次受控对照，给出排除或改善结论 |
| SOFTWARE_DELTA_READY | H参数/状态回归通过，提交身份完整 |

失稳/非有限值/调度越界立即停止加长，保存首次异常位置与字段。边界检查未过，不用调h补偿。输入缺失不阻塞无关标准验证。所有NOT_RUN给原因，不把它写成已解决。

## 12. 产物与本轮必须回答的问题

使用`doc232_runs/A...H`独立目录，避免覆盖231。交付：

- `ROUND232_REVIEW.md`、`ACCEPTANCE_TABLE.csv`、版本与输入身份；
- 231结论更正；真实测点CSV与采样核对；
- 生产热源绝热结果、边界独立积分原始值/误差；
- 真实累计热账本，标明区间能量与瞬时功率的区别；
- ON/OFF匹配对照、5s及条件性20s温度/速度/壁带曲线；
- 电磁基线及一个假设的修复前后数据；
- `SOFTWARE_HANDOFF_DELTA.md`，列源码提交、前置条件、复跑命令及不能声明的能力。

必须回答：

1. 真实生产测点温度能否被离线数据复算？
2. 生产热源与各边界的误差分别是多少，独立参考来自哪里？
3. 壁带温度和57kW散热如何随时间变化，而非仅全域均温如何变化？
4. 开/关浮力差异是否可从启动残余运动中分辨？
5. 原“1000℃不可达”推断应如何修正，当前还缺哪一条证据？
6. 电磁本轮排除了哪个原因或改善了哪个量？
7. 哪些改动值得交给软件机，哪些仅为诊断？

**完成标准是补齐真实生产链的证据，并收窄一个电磁问题；不是增加更多测试名称，也不是把短跑完成改名为稳态或实验通过。**


<!-- ===== END 232_OPHELIE_PRODUCTION_THERMAL_PROBES_BOUNDARY_AUDIT_AND_EM_REPAIR_PLAN_2026-11-22.md ===== -->


<!-- ===== BEGIN 233_OPHELIE_BUOYANCY_CAUSALITY_BOUNDARY_SCOPE_TIMESTEP_AND_EM_CONVERGENCE_CURSOR_PLAN_2026-11-22.md ===== -->

# 233｜物理机任务书：浮力归因、边界适用性、实际步长与电磁收敛

日期：2026-11-22  
执行对象：原OPHELIE物理开发机Cursor。  
依据：第232号任务书、`client_physics_handoff_doc232_2026-11-22.tar.gz`及本轮审核。

## 0. 总目标：收窄四个问题，不继续泛化排障

本轮只集中推进：

1. 同配置20秒开/关浮力对照，确认生产工况中浮力是否真正影响流场；
2. 检验232新增的全域竖向温度梯度外推是否适用于非线性/径向温场；
3. 用真实执行时间记录审计流体、黏性力刷新与热子步，不能用“未中止”授PASS；
4. 解析柱先查自感收敛与误差定义，再决定是否需要一次加密。

保留232已完成的生产测点、真实热源绝热验证和20秒瞬态证据，不再重造这几项。暂不启动几百秒热平衡计算，不改框架、不引入新自然对流论文方案、不在本轮重做弹性桨整圈或GUI。

本轮结束应能回答：**甲方自然对流入口的浮力链是否有效，当前热边界修正是否可用于工程，时间调度是否满足实际限制，电磁误差是否由迭代未充分收敛或参考比较口径造成。**

## 1. A｜身份、已有结论与本轮输入固定

### A1 代码身份先固定

232包记录HEAD=`832dbc782...`，源码dirty；主要运行结果与最新H接口修改未必来自同一次构建。

先检查实际工作区，整理必要源码与配置，形成可重建的本地提交或明确父SHA+补丁哈希。用户已授权持续开发，常规本地提交属于必要版本管理；不再仅以“未要求commit”长期保留不可追溯状态。不自行公开推送、不强推、不覆盖用户改动。

记录可执行文件SHA256、构建命令、实际设备/精度、Reload校验和、完整argv、effective配置和结果目录。若必须保留未提交文件，明确列出，不能只填HEAD代表实际执行版本。

物理机继续使用独立工作分支，软件机保持固定版本。优先把不同功能拆为可择收提交：测点/CLI、边界诊断、调度审计、实际物理修复。

### A2 已确认与需收窄的结论

| 232成果 | 本轮采用的准确口径 |
|---|---|
| 生产测点与离线复算 | 算法/数据路径已通过；名义感温点仍未确认 |
| 生产SYCL热源 | 10kW绝热0.2s升温/能量误差约0.5%，满足原1%开发门限 |
| 边界均匀/竖向线性场 | 该两类制造场通过，不代表任意边界层已验证 |
| 20s热演化 | 均温约1247K、名义热电偶约1257K、散热约53.7kW；未稳态，非真实EM |
| 1s浮力ON/OFF | 信号未分辨，不能授浮力PASS，也不能直接证明浮力公式错 |
| 实际步长审计 | 仍缺真实力年龄/外层跨度证据；估计步长和未中止不足 |
| 电磁柱长对照 | H翻倍无改善，有限柱长不是当前主要误差来源；精度缺口未关闭 |

### A3 固定主配置

采用232甲方NC几何、dp=0.010m、开发液深0.180m、同一Reload；初始均匀1273.15K，rho暂2700、cp1150、k4、冻结mu约14.532Pa·s。

规定10kW均匀开发入热；侧300/底35/顶对流20 W/(m²·K)、发射率0.8、环境/水温300K，来源均按229记录。30kW只是电源输入。

不同时开启mu(T)、sigma(T)反馈、新skull模型或修改h。准确甲方电流如已到达，登记为后续配置版本，不中途污染本轮匹配对照。

## 2. B｜真实时间调度审计（长对照前完成）

### B1 必须记录的实际量

审计写到CSV，至少包含：

- 声学子步开始/结束物理时间及实际dt；
- dt声学候选、惯性/黏性限制、当前mu及估计采用的尺度；
- 最近一次黏性力刷新时间、当前力保持年龄、位置/邻域更新次数；
- 外层实际跨度及相关限制；
- 热步实际dt、最近温度刷新时间和累计热推进时间；
- 结束时间截断及声学/外层/热步累计时间的对应关系。

需覆盖**当前实际生效的旧调度路径**，不能只在`doc227_actual_timestep=1`时生成审计然后给另一条路径授PASS。允许先在0.2–1s检查；若后续速度/限制变化明显，20s保留稀疏审计和最差值。

### B2 判断与修复

逐项检查声学dt是否超自身限制、黏性力保持区间是否超采用的显式稳定限制、热积分是否漏/重复、实际时长是否越过目标。浮点误差给明确小容差，不能以大容差掩盖一个完整子步越界。

若存在越界，复用已有成熟匹配外层预算方案，最后子步截断；修复后ON/OFF都用相同新调度。若无越界，保留当前路径，不为打开某个历史开关而强改代码。

不得同时改变压力、壁面和时间调度来“修复稳定性”。输出审计最差值与违规计数；没有记录时标签为NOT_AUDITED，而非PASS。

## 3. C｜20秒浮力匹配对照与归因（本轮主任务）

### C1 对照身份

两组只改变Boussinesq增量：

- C_ON：beta=1e-5；
- C_OFF：关闭浮力增量，保留重力与原静水支撑。

几何、粒子、初始化、settle、边界修正模式、Q、材料、时间推进、输出阶段完全相同。不能同时切换总压力/扰动压力、自由面/平顶或kernel correction。

优先复用232的20s ON数据，但只有在二进制、源码、参数、边界和调度均能对应时才可复用；否则在固定版本重跑两组，不拿最新OFF比较无法确认版本的ON。

先做0.2s一致性检查，再到20s；已有相同固定版本有效短测可省略重复。热/流动差异比较至少使用后段窗口[5,10]、[10,20]秒，不能把启动平均掩盖后期差异。

### C2 比较场，而不是只相减最大速度

输出两个运行的实际速度、温度、压力/密度及固定空间剖面。粒子运动不同后，按共同有效域映射到同一固定网格；提供插值覆盖掩码。不能按排序后的数组下标直接逐粒子相减。

至少给出：

- 场差`R_u=||u_ON-u_OFF||₂ / max(||u_ON||₂,||u_OFF||₂,u_floor)`，同时报告绝对场差和参考速度下限；
- Urms/Umax、竖向/径向分量、动能和可用环流统计；
- 固定测点温差、截面温差、全域温度/各边界功率差；
- 后段矢量图，统一坐标/色标，区分瞬时与时间平均。

ON/OFF的RMS数值相近不代表两流场相同。差异存在也不能仅凭很小相对差自动认为具有物理意义。

先用确定性快照重放或一组必要的重复短测估计实现/并行波动尺度；本轮不要求无限重复。以高于该尺度、时间上持续且方向合理的响应作为可分辨证据，不预先规定“必须超过10%”以免漏掉真实弱响应。

### C3 不论长窗结果如何，确认实际浮力力场

在本轮生产入口选一个有明确非均匀温场的状态，导出：T、实际T_force_ref、beta、重力方向、粒子质量、浮力专用力、ForcePrior及前一力历史。

独立参考`a_b=-beta*(T-T_force_ref)*g_vector`，具体正负约定与实际g一致；检查热上升/冷下降的方向。T_force_ref固定且可追溯，不得通过移动绝对温标制造温差。

在**同一个冻结状态**比较开启/关闭浮力时的总力差，保存/恢复所有有关字段和历史，不让诊断污染后续运行。检查专用力正确、确实进入总力、重复执行不累加、下一次黏性/重力更新不清掉它、声学动量步读取了有效总力。

优先复用已存在ForcePrior诊断，沿本轮入口验证；不重新造另一套力学框架。独立参考与力增量在合理浮点容差内一致，原有1e-5相对门限可在有足够非零尺度时采用；接近零值用绝对门限并报告，不除零。

### C4 结果分支与停止条件

1. **力链正确，20s场响应可分辨：** 该工况浮力机制通过归因检查；尚不代表实验温度或法国流线对齐。
2. **力链错误：** 先修具体丢力/温度/参数/调度问题；短测和同配置复测后再延长，不另扫材料或h。
3. **力链正确，但响应仍弱或未分辨：** 补一个有限诊断，使用已有非均匀温场（含水平温差）驱动同一入口，关闭热演化以隔离热场生成；只说明该机制能否响应，不算甲方热过程已验证。不能仅增大beta凑环流。

若需要减小启动噪声，只选一个有理由的静水初始化或settle方案做对照；改变后ON/OFF共同采用。不要同时更换整个求解器。

## 4. D｜全域温度梯度外推的适用范围与隔离

232当前实现对全域粒子拟合单一dT/dz，并默认用于顶/底面外推。它对竖向线性场成立，但不保证适用于局部边界层；侧壁也缺少对应径向处理。

### D1 先变为明确模式，保留旧行为

给边界表面温度估计一个有来源的显式模式：旧边界带取样、232全域竖向拟合、其他已验证局部方案。具体字段名适配现有CLI，不新增多重冲突开关。

232全域拟合在完成本节前视为实验性，不静默改变软件示例。无需立即开发新局部重构，先测现状和旧方案。

### D2 两个有判别力的制造场

复用真实生产边界与独立几何积分，固定粒子状态和材料：

- 径向场：`T=1200+100*(r/R)^2` K；检验侧壁表面温度与内部带温度差，以及顶底径向非均匀辐射积分。
- 局部竖向非线性场：`T=1200+100*(2*(z-zb)/H-1)^2` K；对称时全域拟合梯度接近零，但顶/底附近真实梯度非零，检验全域外推是否失效。

两者都在甲方指定温区内，测试不涉及低温物性外推。参考直接在解析几何面计算T并积分，不用生产粒子权重生成参考。

逐项给出四类边界功率与单步真实能量变化。沿用原开发门限：各分项功率差≤2%，实际能量对账≤1%；同时报告原始误差，不能只改能量口径以获得PASS。

### D3 若失败，按原因做最小修正

- 全域拟合无法表示局部梯度：限定其用途，不直接在工程中默认启用。
- 若已有成熟局部重构/壁面温度方案，复用并对同一测试复测。
- 若确为粗边界带空间分辨率误差，可做一次有依据的厚度/dp检查，不通过任意修改h补偿。

不因旧模式也不够准就默认新模式正确。工程模式使用范围、误差和临时近似应明确登记。

**边界模式一旦改变，不能把新20s结果直接与旧232曲线混为同一配置。** 浮力对照在同一选定模式内进行；不强制把每个实验模式都跑20s。

## 5. E｜质量焓与算子能量口径并列检查

232绝热生产源约0.5%误差已满足原开发阈值，本轮先解释，不强迫无依据追零。

同一实际状态记录：

- 固定粒子质量`sum(m_i)`；当前`sum(V_i)`、密度范围及`rho0*sum(V_i)`；
- 质量焓`DeltaH_m=sum(m_i*cp*DeltaT_i)`；
- 算子口径`DeltaH_rhoV=sum(rho0*V_i*cp*DeltaT_i)`；
- 实际热源/边界累计能量、热子步跨度。

温变cp未启用时先按常cp。初场/参考态初始化后建立基线，不把换初场计为能量不闭合。

正式运动结果以固定质量焓为主要物理账本，rho0V作辅助离散诊断。如果因WCSPH密度/体积变化出现差异，要量化其规模与时间增长趋势；不能只把门改为更容易闭合的口径并宣布问题消失。

若差异在约定≤1%且不增长，记录限制即可；若持续增长或解释不了，定位质量、体积、源归一和温度积分链。不要顺带改变cp/rho/h去掩盖。

## 6. F｜解析柱：先收敛与指标，后分辨率

### F1 更正指标含义

232的`L2_rel_J_mid≈0.00042`与平均电流误差4.18%、中部功率误差17.8%同时存在。先检查实际指标源码/后处理，确认是否经过幅值归一、最优比例拟合、分箱或区域重定义。

新增/明确未经幅值拟合的比较：

`E_J = sqrt(sum(V_i*|J_i-J_ref(x_i)|²) / sum(V_i*|J_ref(x_i)|²))`

复数场比较真实复分量；若只能比较幅值，另命名，不冒充复数相位精度。参考在相同粒子坐标求值，区域与体积权重固定。

径向分箱时解析参考也对同一采样集合/同一权重平均，不能把bin中心解析值直接当稀疏粒子的均值参考。

功率并列报告：数值Q的离散积分、相同粒子采样上的解析Q积分、解析几何体积分。这样区分场误差和体积积分误差，不能删除不利区域。

### F2 一次严格收敛对照

保持232基线dp=0.025m、R0.1m、H0.6m、sigma16、300kHz及相同源输入。日志原固定点残差约0.0396、最终一致性残差约0.0276；它们与真实场误差不同，但需要检查其影响。

先仅收紧迭代容差，使真实固定点/一致性残差达到约1e-3或当前可可靠达到的更低值，迭代上限设有限预算；记录未松弛残差、松弛步长及各物理指标随迭代变化。不要只看较小的relaxed_step宣称收敛。

若达不到容差，报告停滞原因和残差历史，不降低目标后标严格收敛PASS。必要时一个中间容差对照可帮助判断，但不展开多维扫描。

比较未经归一的J误差、径向误差、全柱/中部功率与232基线；不得用目标功率重新标定电流掩盖误差。

### F3 是否加密的决策

- 严格收敛后误差明显下降：先固定收敛设置，不同时加密。
- 误差基本不变、指标采样已澄清：允许一次径向分辨率加密，其他条件固定；检查内核稀疏采样和源/自感离散。
- 不再继续增加柱长；232已显示该因素不是主因。

本轮不以甲方电流未知为由暂停标准算例。TEAM7不要求同时重跑全套，保留明确未通过状态；解析柱一项通过不能替代TEAM7/峰值Q的验证。

## 7. G｜提交后的真实短跑与软件交接

232最新H源码包含`em_status=NOT_RUN`、规定热源身份和空EM结果字段，但部分20s日志仍是旧`em_solver_ok=1`口径。属于源码/运行证据版本差异，不能只用CLI自检盖过真实运行。

在整理后的提交上至少做一次0.2s实际规定热源短跑，检查：

- effective/运行事件/终端最终摘要/summary对热源模式含义一致；
- EM未执行时明确NOT_RUN，`em_solver_ok=null`或等价不可用状态，不是假成功；
- 生成器输入与规定吸收功率分开，实际热输入与配置一致；
- 固定测点在实际输出中可读；
- 非法输入和冲突顺序防护沿用232有效回归，不重复造接口。

软件机优先接收测点、CLI防护、状态一致性等通过提交。边界实验性外推和未确认物理修复暂不默认合入。提供父SHA、提交、必要前置依赖、配置/格式变化和复跑命令；禁止整文件覆盖软件修改。

## 8. 执行顺序、预算和停止条件

| 优先级 | 内容 | 依赖及限制 |
|---|---|---|
| P0 | A身份、B实际步长短审计 | 长运行前；若越界先修 |
| P1 | C力链快照和20s浮力对照 | 同版本同边界；尽量复用可核验ON |
| P2 | D边界制造场、E双口径能量 | 静态/短步即可，不用每项都长跑 |
| P3 | F电磁指标+严格收敛对照 | 可用空闲设备独立推进，目录/设备隔离 |
| P4 | G固定提交真实短跑与交接 | 只交已测功能，不整分支自动同步 |

可先快速做静态D判定边界模式风险，再安排C长对照；保持每组对照内部一致。不要同时测试多种压力、h、beta、dp、阻尼和材料。

出现非有限值、实际步长违规、粒子明显失效或能量持续异常时停止延长，保存首次状态。输入不足标BLOCKED，继续不依赖该输入的任务。

本轮无必要不跑超过20s，不追求立即热平衡或1000℃误差达标。现在的初态就是约1000℃，短时接近目标不能作为预测验证。

## 9. 文件与状态交付

建议`doc233_runs/A...G`，不得覆盖232结果。必要源码、分析脚本和小型配置可追踪；大场/Reload单独保存校验和。

交付：

1. `ROUND233_REVIEW.md`、`ACCEPTANCE_TABLE.csv`、版本与运行身份；
2. 实际步长/力刷新年龄审计CSV及最大违规表；
3. 20s ON/OFF速度/温度场差、后段统计、固定测点曲线与力链快照验证；
4. 两个新边界制造场的独立参考、原始误差及模式选择；
5. 固定质量焓与rho0V口径差异；
6. 电磁指标定义、严格收敛前后值、如有则一次加密对照；
7. 提交后真实短跑的schema/EM身份核查和软件择收说明。

状态明确区分：ACTUAL_TIMESTEP_AUDITED、BUOYANCY_FORCE_CHAIN_VERIFIED、BUOYANCY_FLOW_RESPONSE、BOUNDARY_RECONSTRUCTION_SCOPE、MASS_ENTHALPY_ACCOUNTING、EM_ITERATION_SENSITIVITY、SOFTWARE_RUNTIME_IDENTITY。

本轮必须回答：

- 浮力是否实际进入生产动量推进，20s场响应是否可分辨？
- 232全域梯度外推在哪些温场有效、哪些无效，软件是否可默认使用？
- 当前实际推进是否跨过黏性/声学/热预算？
- 约0.5%的绝热能量差来自何处，运动后是否增长？
- 电磁4%左右幅值和17.8%中部功率误差中，有多少来自收敛或参考采样？
- 哪些提交适合软件接收，哪些仍是诊断？

**禁止用“跑完”“曲线相似”“内部残差小”替代以上因果与误差证据；也不因结果未通过而临时更换指标或删除不利区域。**


<!-- ===== END 233_OPHELIE_BUOYANCY_CAUSALITY_BOUNDARY_SCOPE_TIMESTEP_AND_EM_CONVERGENCE_CURSOR_PLAN_2026-11-22.md ===== -->


<!-- ===== BEGIN 234_OPHELIE_FORCE_PRECISION_FIXED_T_RESPONSE_AND_EM_SPATIAL_REPAIR_CURSOR_PLAN_2026-11-22.md ===== -->

# 234｜OPHELIE 浮力归因收口、热边界修复与电磁空间误差定位执行任务书

- 日期：2026-11-22
- 执行对象：物理开发机上的 Cursor。
- 上轮基线：第233轮交接，分支 `ophelie-physics-doc233`，报告 HEAD `e84f3c3bc`，父 tip `832dbc782`。
- 本轮建议分支：从实际核实的233物理提交新建 `ophelie-physics-doc234`。
- 软件机继续独立开发。物理改动经明确验证和选择性移交后接收，不要求软件机立即整体合并。
- 本文是执行计划，下面的候选原因均不是已查明的代码缺陷。

## 0. 本轮必须回答的问题

本轮目标不是增加更多 PASS/FAIL 标签，而是用少量有区分力的试验查明：

1. 0.163%的浮力总力差，是单精度消减误差、诊断位置不对，还是浮力被覆盖/重复/没有被实际动量推进读取？
2. 黏性力“年龄违规”是真正超出了更新预算，还是时间戳、比较尺度或计数逻辑有误？
3. 在固定、非均匀温度场下，当前工程流体路径能否产生与浮力对应的速度场变化？
4. 上下热边界对非线性温度场的误差，能否用有限范围的局部处理降低？
5. 电磁解析柱的中部功率17.8%误差，多少来自场离散，多少来自采样积分及参考定义？

完成这些工作后，再决定是否恢复客户温度演化以及法国自然对流复现。不要在未回答前重新铺开长时计算或多参数扫描。

### 0.1 第233轮应保留的事实

| 项目 | 已有结果 | 本轮正确处理 |
|---|---|---|
| 专用浮力通道 | 公式对比误差约4.06e-8 | 保留，不重做公式推导 |
| 总力 ON−OFF | 相对专用浮力约1.63e-3 | 未达到旧门限，但原因未定 |
| 20 s ON/OFF | 速度统计量相近 | 不能等价为真实矢量场差小 |
| 步长审计 | 声学未报告越界；黏性年龄计数与比值待解释 | 先审计测量代码，再改调度 |
| 热时间累计 | 为0 | 说明记录缺失，不能说没有热推进 |
| 全域梯度外推 | 抛物竖向温场失败 | 保持实验模式 |
| 质量焓 | 20 s残差约0.4%，按累计散热归一 | 不能标成1e-4量级 |
| 电磁严格收敛 | 固定点约7.1e-4，但物理误差几乎不降 | 查指标和空间误差，不继续主要调迭代 |
| 软件规定热身份 | EM NOT_RUN已明确 | 保留，避免规定热伪装成真实EM |

### 0.2 不在本轮扩展的内容

- 不替换整个流体框架，不重新混合2024论文算子、黎曼压力、密度稳定等多个方案。
- 不修改独立 A-phi/Laplace 研究分支。
- 不引入新壁面 ghost 粒子布置框架。
- 不同时开启 μ(T)、σ(T)、新水冷模型、导电坩埚或底板。
- 不通过增大β、降低μ、重调散热系数制造“更好看的环流”。
- 不重跑标准方腔全套；若修改公共流体内核，才补与修改直接相关的既有回归。
- 不要求搅拌整圈或重做弹性桨；本轮不改变已交接的不锈钢弹性桨设置。
- 不把1000℃作为当前工况已通过实验匹配的证据。

## 1. A｜版本、文件及结论身份

### A1 开工核对

1. 记录实际分支、HEAD、工作树状态、编译配置、Real类型、GPU和二进制SHA256。
2. 如果工作树有未提交代码，先明确哪些构成233结果，哪些是后续改动；不得覆盖或清理用户文件。
3. 本轮所有匹配对照使用相同二进制、Reload和有效参数。需要修改代码后，重新生成版本记录。
4. 不使用未记录的环境变量静默改变边界、精度、调度或材料参数。
5. 确认软件机当前固定版本；本轮只准备可选择接收的补丁，不直接影响其工作树。

### A2 本轮目录与提交建议

代码优先放在既有 OPHELIE case/diagnostics 内。公共库修改必须有直接证据支持。

产物统一放到 `docs/ophelie/reference/doc234_runs/`：

| 子目录 | 内容 |
|---|---|
| A_identity | 版本、构建、参数、运行清单 |
| B_force_precision | 同一快照力链、舍入分析、实际动量读取 |
| C_schedule | 审计定义、真实时间跨度、必要修复前后 |
| D_fixed_T | 冻结温场、短时速度场、Python矢量图 |
| E_boundary_energy | 制造温场、局部边界验证、双能量口径 |
| F_em | 指标定义、解析参考采样、一次加密 |
| G_handoff | 有效配置、验收表、软件接收建议 |

按功能拆分提交：诊断、必要物理修复、后处理、报告分开。不得把大量运行产物全部提交进Git。

## 2. B｜先查浮力力链：区分浮点误差与真实覆盖

### B1 旧判定不能直接继承

233的总力差相对误差约0.163%，构建为Float。在含重力的总力上做ON−OFF，可能丢失小浮力的有效数字。

估算只用于选诊断：β=1e-5、ΔT=10 K时，浮力约为重力的1e-4。单精度大力相减的相对浮力误差可能显著放大。但最终是否由此造成，必须用实际每粒子力尺度验证。

禁止两种做法：
- 未检查舍入与通道，就断言“浮力没有进入动量”；
- 仅把1e-5门限改宽，然后宣布已解决。

### B2 同一粒子快照，先做不积分的检查

固定相同的Position、Temperature、Mass、Density、Pressure、邻域、时间和所有既有力通道；从同一快照分支计算ON/OFF，禁止先积分ON再拿变化后的状态算OFF。

独立参考使用：
`a_b_ref = -beta * (T - T_force_ref) * gravity_vector`
`F_b_ref = mass * a_b_ref`

在代码中明确gravity向量的符号以及Force/Acceleration单位。T_force_ref必须读实际有效值，不能一边用当前质量均温、一边用初始温度。

每粒子导出至少：
- 稳定粒子ID、位置、T、T_force_ref、mass；
- gravity、独立a_b_ref、专用F_b；
- 加浮力前ForcePrior、加后ForcePrior；
- ON/OFF差；
- 浮力之后其他动力学修改的ForcePrior；
- 进入声学动量更新前真正读取的力或加速度；
- host/device同步位置与数据类型。

如果实际API的力通道名称不同，使用真实名称，报告映射关系；不要创建一个无人读取的新变量来“验证通过”。

### B3 舍入分析与必要的双精度对照

1. 在host用double独立计算同一已导出的Float快照参考。double重算参考本身不能证明Float总力相减准确。
2. 对生产Float的加法和减法，给出与该位置实际大力大小相关的ULP或舍入尺度估计。不要仅用全局平均力。
3. 双精度对照优先局限到诊断算术/独立诊断目标；若条件允许可用小规模CPU double参考，不要求重编整套软件或全量FP64 GPU。
4. 确认专用Fb是否被加入一次、是否被后续初始化/重置覆盖、是否有前一子步残留。
5. 记录“诊断前后读取到host”的同步是否正确，排除只比较了旧设备数据。

误差同时报告绝对量与相对量。对接近零的粒子，禁止直接除以接近零的Fb后取最大相对误差。报告全场体积/质量权重L2以及非零浮力粒子分组统计；零温差粒子另报绝对误差。

### B4 最少三种制造温场

只做同一几何快照，无须完整流动：
- 均匀温度：若T_force_ref取同一均温，应得到近零浮力增量；
- 原有制造非均匀温场：复现233的0.163%；
- 有正负温差的非均匀温场：检查冷热粒子的方向与符号。

这些是力链单元诊断，不是自然对流拓扑验收。

### B5 实际动量通道检查

至少跟踪一个真实声学子步，证明Boussinesq执行后的有效力被生产动量更新读取，并记录其后是否被其他动力学改写。

不能要求完成压力响应后的速度逐粒子严格等于a_b*dt：压力、黏性和约束会共同作用。要检查的是实际离散更新中浮力项被读入的贡献，以及必要的力/质量换算是否正确。

### B6 结论分支

| 结果 | 后续 |
|---|---|
| 专用Fb正确、动量确实读取、Float差异被舍入尺度解释且double诊断显著改善 | 记为精度受限的诊断，不修改物理力链 |
| 存在覆盖、漏读、重复累加、错单位或错误参考温度 | 只修具体链路，重做同一快照与短时响应 |
| 误差明显超出舍入解释、来源仍不清 | 保持OPEN，报告残差粒子分布与具体调用位置，禁止继续堆长跑 |

验收不依赖单个硬相对门限；必须同时具备独立公式、实际读入力、绝对误差和精度解释。不能仅做host数学演示。

## 3. C｜修正实际调度审计，查明“年龄违规”

### C1 先核对定义

233记录：max_viscous_age约0.00238 s，max_age_over_visc约2.1e-4，却几乎每步计数违规。

必须在报告中列出：
- force_age的起止事件；
- age_over_visc的确切分子、分母、单位；
- 允许的年龄预算来自哪个生产调度公式；
- 违规计数的确切布尔表达式；
- 每次计数发生在刷新之前还是之后。

如果比值不是“年龄/允许黏性跨度”，更名并说明，不沿用误导名称。

### C2 每个外层/力刷新区间记录

至少记录：
- 物理t起止、实际外层跨度；
- 每个声学dt及实际累计；
- 黏性力刷新时刻、当前年龄、该刷新所对应的允许保持预算；
- 实际使用的声学、平流、黏性、热限制及其来源；
- 粒子位置、cell、邻域更新的事件次序和计数；
- 热模块实际执行的子步dt和累计热时间。

黏性力年龄不能直接与一个无关的声学dt比较。显式黏性、隐式黏性与冻结力外层保持的限制也不能混为同一公式。先识别本case实际采用的路径。

### C3 修复热累计记录

在真正调用温度更新的路径中累积其物理dt，记录是否分裂积分、是否多次调用。若采用分裂子步，按完整时间区间解释累计，不能因重复计数误报为推进两倍时间。

验收：
- 实际热推进覆盖与物理时间一致，误差在端点对齐及可解释舍入范围内；
- 没有热路径时明确NOT_RUN，不用0伪装通过；
- 记录值可由短跑原始行重新累计。

### C4 运行与修复预算

先0.2 s；如果事件样本不足，最多延至1 s。
- 如果只是审计错误：修审计，保持物理调度不变。
- 如果证实真实预算违规：按现有框架兼容方式修正外层跨度/刷新位置，重做相同短跑。
- 未证实前禁止直接翻转旧CFL开关、全面改变声学顺序或更换黏性算子。

输出至少一段可读事件轨迹，解释一次完整外层及子步。不要只给累计PASS/FAIL。

## 4. D｜固定温度场短测：直接确认浮力驱动

### D1 先复用已有数据，不先重跑20 s

检查233原始产物是否含指定时刻速度矢量、粒子坐标或足够重建的网格场。只有速度RMS、单个环流代理或PNG不算完整矢量场证据。

如果缺场，明确写MISSING_FIELD_DATA，不继续声称“场差只有1%”。233窗口统计差只能作为统计量差。

### D2 固定物理设置

仍用客户NC圆柱、dp=0.010 m、H=0.180 m、相同Reload；冻结μ约14.532 Pa·s、ρ2700、β1e-5，保留同一重力、压力、壁面、自由面设置。

本项关闭热源推进、热扩散和热边界温度更新，直接施加固定Eulerian非均匀温场：
`T(r,z) = T0 + 20 K * (r^2/R^2 - 1/2)`
其中T0=1273.15 K，R为本case实际液体半径。

该函数具有径向温差；不同高度在同一半径温度相同。面积均值为T0，离散均值仍按生产T_force_ref机制记录。它不是法国论文温度重建，仅用于验证浮力响应。

每次生产浮力读取之前按当前粒子位置更新规定T，避免粒子携带初始标签后温场已变。ON/OFF共用完全相同的温度定义。

为什么不用纯竖向线性温场作为唯一测试：它可能对应静力平衡或稳定分层，不能保证产生目标环流。径向温差用于避免把“无环流”误读成力链错误。

### D3 最小运行矩阵

| 运行 | 温度 | 浮力 | 时间 |
|---|---|---|---|
| D_ON | 上述固定场 | β=1e-5 | 先1 s |
| D_OFF | 相同固定场 | 仅关闭Boussinesq | 同1 s |

关闭浮力只关闭浮力增量，保留重力和所有其他设置。若需要静水预处理，两次使用相同经过核实的预处理状态及温度施加时刻，预处理时长另记。

先通过B/C，再执行本项。如果1 s响应太弱，结合实际力尺度与数值噪声确定是否延至5 s；不得无条件再跑20 s。温差不放大，μ不降低。

### D4 场保存与Python检查

至少保存温场施加起点及0.1、0.5、1 s；若延长，补2、5 s。输出粒子ID、位置、速度、T和必要权重，不必每步写VTP。

生产速度是原始证据。后处理可将ON/OFF插值到相同Eulerian网格，并导出：
1. ON、OFF和差值的同色标速度剖面；
2. ON/OFF原始矢量图及差值矢量图；
3. 若计算流线，说明插值、有效区域、积分设置；
4. 质量/体积加权的矢量差L2与其绝对m/s数值；
5. 竖向速度径向剖面，以及内部上升/下降区域；
6. 末时刻及时间演化，不能只挑一张最佳图。

粒子移动后不能把同一粒子ID直接当作同一空间坐标。可附Lagrangian差，但与Eulerian场差分开命名。

圆柱轴对称平均时，用真实r与z、旋转投影后的u_r与u_z；轴附近u_r方向不适定，需正确处理。报告周向离散/非轴对称程度，避免平均抹去真实结构。

有效网格必须ON/OFF都有支持，记录覆盖率；容器外不补零。流函数是后处理派生量，若连续性不足，不以闭合流函数图替代原始速度证据。

### D5 判据与分支

- B证实浮力进入动量，且D出现与冷热位置相符、时间上可解释的ON/OFF矢量响应：可记“固定温场浮力响应验证通过”；不等于法国双环流通过。
- 差异接近插值或运行波动尺度：只增加一次有针对性的重复或更直接的内部粒子比较，估计误差底；不能仅按ON/OFF统计差给噪声标签。
- 力链正确而流动响应弱：检查实际温差、黏性、压力平衡、边界约束及静水背景，量化竞争尺度，再决定延长至5 s是否有意义。
- 若OFF背景持续主导，先定位OFF的静水/自由面/边界残余流；不要把它命名为温度驱动自然对流。

不强制要求短时间内达到闭合大环流，也不要求速度严格随时间线性增长。

## 5. E｜上下热边界的有限修复及能量口径纠正

### E1 更正既有能量结论

233 ON 20 s摘要：
- E_Q=199976 J；
- E_boundary=1105680 J；
- 输入减散热=-905704 J；
- H_mass=-910047 J；
- 固定质量焓残差约-4343 J，按累计散热归一约0.393%。

原报告约1e-4的残差更接近rho0V口径，不能称作固定质量焓残差。本轮必须由原始日志完整精度重新计算；上述数值只用于指出口径问题。

新增明确字段：
`deltaH_mass_J`、`deltaH_rho0V_J`、`E_source_J`、`E_boundary_out_J`、
`residual_mass_J`、`residual_rho0V_J`、各自归一分母及相对值。

固定质量焓计算用起止温度差：
`DeltaH_mass = sum(m_i * cp * (T_i(t)-T_i(0)))`
本轮cp固定。使用绝对焓时参考温度必须固定；不得把每步变化的质量均温当作焓零点，消去真实能量变化。

相对残差给统一、公开的分母，例如max(|E_source|,|E_boundary|,能量数值下限)。低输入低散热时同时给绝对J，不用极小分母制造巨大相对值。

### E2 保留两种边界模式身份

- band：旧带取样基线；
- global_dt_dz：实验全域梯度外推。

不能只修改日志名称而实际仍默认使用实验外推。effective配置、命令、summary和算子实际路径必须一致。

233的抛物场失败说明全域梯度不能代表局部边界梯度；不证明所有局部边界都需重写。

### E3 局部修复的实施顺序

1. 审查当前换热面面积权重、作用带厚度、表面温度估计、热通量和能量注入是否一致。
2. 在固定粒子快照上复现均匀、linear_z、radial_quad、parabolic_z四种制造场。
3. 若主要误差来自表面T取样，优先在上下表面使用已有邻域支持的局部重构；沿局部法向取内侧粒子，检查支持量与拟合条件。
4. 若没有可复用的稳定局部实现，先做一次dp减半的制造场对照，区分带取样阶数误差和代码错误；不要立刻创建新通用边界框架。
5. 局部重构增加明确实验模式，报告支持不足时的回退。拟合失败不能输出未定义温度或静默大幅限幅。
6. 若为通过抛物场选择二次局部拟合，必须额外保留一个不同形式的光滑温场或独立交叉验证；不能用制造场解析公式直接填回生产表面温度。

### E4 边界验收

本轮沿用制造场功率误差≤2%、能量注入误差≤1%的内部目标，并分别列每一面，不能只报总散热抵消后的误差。

解析参考必须对实际几何面温度独立积分。辐射项用T^4及一致绝对温标，不把平均T的四次方冒充四次方的平均值。明确辐射环境温度和对流环境温度。

区分：
- 表面温度重构误差；
- 面积/通量离散误差；
- 已算出的通量是否正确更新温度；
- 全局固定质量能量残差。

通过并不证明强边界层已收敛；失败也不能直接归因浮力公式。

本项以静态制造温场为主，不要求每一种边界模式再跑20 s。

## 6. F｜电磁：先把误差定义查清，再做一次空间加密

### F1 基线固定

继承233解析圆柱：
R=0.10 m、总高0.60 m、dp=0.025 m、B0=0.01 T、
sigma=16 S/m、f=300000 Hz，源定义及相位约定不变。

使用已达到约1e-3真实固定点残差的严格收敛配置，记录真实残差、松弛步长、phi残差以及最终一致性检查。不得只凭converged字段。

本轮不继续增加柱长；232柱长对照已显示该调整没有明显改善主要指标，但不要将其写成排除了所有三维边缘效应。

### F2 指标必须追溯源码与参考

逐项列出L2_rel_J_mid、mean_rel_J、rel_P_mid、rel_P的：
- 精确数学定义与实际代码位置；
- 数据区域；
- 粒子/网格/径向bin支持；
- 权重与归一；
- 是否幅值拟合；
- 是复数场、向量幅值还是某个分量；
- peak/RMS约定。

新增未经幅值拟合的复数场误差：
`E_J_complex = sqrt(sum(V_i*|J_i-J_ref(x_i)|^2) / sum(V_i*|J_ref(x_i)|^2))`

如果只能比较幅值，另报E_J_amplitude，不冒充相位验证。参考在相同粒子坐标上计算。用同一体积权重及区域，保留绝对量。

径向bin参考也在同一组粒子坐标采样并加权，不能把bin中心解析值当整组粒子参考。

### F3 功率拆成三个量

针对全柱和明确中部切片分别计算：
1. P_num_discrete：数值J/Q在粒子体积上的积分；
2. P_ref_discrete：同一粒子集合、同一权重上的解析Q积分；
3. P_ref_geometric：同一几何区域的解析积分。

给出：
- 数值场相对同采样参考的功率差；
- 解析离散采样相对解析几何积分的差；
- 最终数值功率相对几何参考的差。

若粒子代表体积与切片边界相交，明确采用点选权重还是截断权重；解析几何区域必须与定义匹配。检查有效切片体积及其相对几何体积误差。

Joule项明确峰值复幅值下的1/2约定或RMS约定。禁止用目标总功率归一掩盖标准例误差，也不能删掉内圈不利数据。

### F4 一次加密，禁止再停在“迭代误差不敏感”

指标采样已澄清后，从dp=0.025 m加密到dp=0.0125 m，其他物理设置不变。先检查显存、相互作用复杂度及最大迭代成本；若预计超出可用预算，选择一个更温和的中间分辨率并提前记录理由，不改物理输入。

这是一次空间误差对照，不是要求证明严格收敛阶。记录：
- 粒子数、区域体积；
- 严格迭代的最终残差；
- 时间与显存；
- 同采样复场误差、径向J/Q、三种功率；
- 内圈支持数与自项/近邻处理。

若加密后迭代残差未达到相同量级，先注明无法干净归因，不把结果直接当空间收敛。

### F5 依据结果决定修复对象

| 结果 | 下一步 |
|---|---|
| P_ref_discrete与几何参考差明显，加密后缩小 | 优先查体积、切片和积分采样 |
| 同采样J/Q误差随加密明显下降 | 空间离散为重要因素，记录工程分辨率限制 |
| 体积和采样正确，但同采样场误差不降 | 定位自项、近邻核、源幅相、实际边界设定；只选证据最强的一处修 |
| 功率约定有固定因子错误 | 修peak/RMS或1/2约定并重做基线，不与分辨率收益混算 |
| 指标定义仍不清 | 保持OPEN，不用较小L2宣称电磁精度通过 |

本轮无需等待甲方电流参数才能跑解析标准例。TEAM7、局部峰值Q仍是独立验证缺口；解析柱结果不能替它们背书。

## 7. G｜恢复温度演化与法国复现的进入条件

只有以下条件具备，才恢复一个短的客户热演化回归：
- 浮力进入生产动量的证据已明确；
- 固定温场ON/OFF矢量响应可解释；
- 实际调度审计的口径矛盾已解决；
- 本次使用的热边界模式已明确，残差字段口径正确。

恢复时只跑一个已选定配置，先1 s，必要时5 s；不展开所有模式的20 s矩阵。规定10 kW保持为诊断热源，记录EM NOT_RUN。

客户温度定量与法国文献的自然对流拓扑不是同一个验收：
- 客户几何、输入电流、测点位置/时间尚需对应实验；
- 法国NC几何与热条件需另用明确配置；
- 不要求客户短跑复制法国论文的双环流；
- 不因短跑均温接近1000℃就标实验通过。

当前输入10 kW、短窗散热约54 kW，说明此时净冷却。它不是整个模型在所有时间都无法达到稳态的证明；也不能把当前瞬态散热直接当最终所需吸收功率。

本轮不重开μ(T)/σ(T)。待上述链路明确后，再单独安排真实EM热源与温度反馈。客户暂借sigma=16的状态仍需显式标注，不能称甲方材料已经核定。

## 8. 执行顺序、运行预算和停止规则

### 8.1 推荐顺序

1. A：身份及233原始数据核对。
2. B/C：源码级诊断审查与0.2 s短跑，优先纠正误判。
3. D：冻结温场1 s匹配对照，保存原始场并出Python图。
4. E：静态边界制造场、必要局部修复及能量账本纠正。
5. F：电磁指标和采样核对，再执行一次加密。
6. G：满足条件后做一次短回归，整理可交接内容。

电磁与静态制造场可在独立目录执行；共享GPU要避免资源争抢。共用bin/output的流体运行必须串行或显式独立输出目录，不能互相覆盖。

### 8.2 本轮运行上限

- 力链：静态制造场及至多一个真实子步轨迹；
- 调度：0.2 s，必要时1 s；
- 固定温场：ON/OFF各1 s，仅有必要时各延至5 s；
- 边界：静态制造场；必要时一次分辨率对照；
- 电磁：严格基线证据充分则复用，加一次分辨率；
- 恢复热演化：有条件的单配置1–5 s。

这些是默认预算，出现明确错误时先修错误，不机械跑完全部矩阵。不要以还没跑完矩阵为由拖延已经定位的修复。

### 8.3 停止与转向规则

- 诊断误差被舍入解释：停止追求不合理Float相对门限，转到实际动量响应。
- 调度只是计数错误：停止改流体顺序。
- 固定温场无清晰响应：停止耦合长跑，保留场与实际力证据定位。
- 严格EM容差无改善：停止继续主要调容差，转采样/空间离散。
- 局部边界支持不足：保持实验模式，报告支持问题，不静默回退后称全部通过。
- 本轮未触及的TEAM7和整圈搅拌，列遗留，不写成新近失败或已解决。

## 9. 交付文件与报告要求

这次不能只回传摘要。交接至少包含：

1. ROUND234_REVIEW.md：按“问题→证据→归因→修复→剩余影响”撰写。
2. ACCEPTANCE_TABLE.csv：区分PASS、OPEN、FAIL、NOT_RUN和受限适用。
3. RUN_COMMANDS.md：实际运行命令、执行顺序、工作目录、产物路径及退出码；示例命令与已执行命令分开。
4. effective配置及版本、二进制、Reload校验信息。
5. 本轮源码补丁/提交差异，包含力链审计、调度计数、边界和EM指标定义。
6. 力链每粒子关键CSV与实际读入力事件轨迹。
7. 实际调度CSV原始行、热dt累计、统计脚本。
8. 固定温场指定时刻粒子场/剖面CSV、共同网格、Python脚本及ON/OFF/差值矢量图。
9. 每个边界制造场逐面功率、能量、独立参考计算，以及局部支持失败统计。
10. 双口径焓账本CSV与明确的残差分母。
11. 电磁原始日志、径向数据、复场误差、三种功率及一次加密对照。
12. SOFTWARE_HANDOFF_DELTA.md：可以择收什么、什么仍不能进入工程默认。

完整大场可外置，但包内必须有稳定索引及SHA256；关键稀疏场、CSV、脚本、补丁和足够核对的日志片段应直接打包。不只写一个审核者访问不到的本机路径。

### 9.1 验收表建议字段

`item,status,scope,metric,absolute_error,relative_error,normalization,config_sha,evidence,limitation`

建议项目：
- FORCE_FORMULA_VERIFIED；
- FORCE_CHANNEL_AND_MOMENTUM_READ_VERIFIED；
- FLOAT_CANCELLATION_EXPLAINED；
- ACTUAL_SCHEDULE_AUDITED；
- THERMAL_TIME_ACCOUNTED；
- FIXED_T_BUOYANCY_RESPONSE；
- BOUNDARY_FACE_POWER_SCOPE；
- MASS_ENTHALPY_ACCOUNTING；
- EM_METRIC_DEFINITION_VERIFIED；
- EM_SPATIAL_COMPARISON_COMPLETED；
- SOFTWARE_SELECTIVE_DELTA_READY。

不要让一个顶层passed=1覆盖尚未通过的物理项。既有passed/smoke_passed如仅表示程序结束，要明确其作用域。

### 9.2 最终必须回答的六句话

1. 0.163%力链差异到底由什么造成，是否影响真实动量输入？
2. 黏性年龄是否真的违规，热子步是否完整覆盖物理时间？
3. 固定温场的实际速度场是否显示可解释的浮力贡献？
4. 非线性表面温度下的热通量是否改进，在哪些范围仍有限制？
5. 电磁17.8%中部功率缺口主要落在场、采样还是尚未定位？
6. 下一轮应进入温度驱动耦合，还是只剩哪个明确缺陷需要修？

## 10. 发给Cursor的执行约束

请按本任务书自行完成已明确授权的诊断、局部修复和有限测试，不要在允许的一次电磁加密前再次等待一般性确认。

涉及真实预算不可承受、缺关键源文件、或与本机用户明确约束冲突时，记录具体阻塞及已完成证据。不要凭猜测填补数据或把OPEN改PASS。

允许准确的CPU诊断和小规模双精度参考；生产GPU路径尽量保留。不要为了全GPU名义放弃必要验证，也不要因此将整个项目切换求解框架。

本轮成功标准：用可信测量查明最关键的因果链，完成最小必要修复，并使下一步只面对明确、有限的物理问题。


<!-- ===== END 234_OPHELIE_FORCE_PRECISION_FIXED_T_RESPONSE_AND_EM_SPATIAL_REPAIR_CURSOR_PLAN_2026-11-22.md ===== -->


<!-- ===== BEGIN 235_OPHELIE_OFF_FLOW_EM_GEOMETRY_AND_DIAGNOSTIC_CLOSURE_CURSOR_PLAN_2026-11-22.md ===== -->

# 235｜OPHELIE 无浮力背景流归因、电磁积分域修正与验证收口任务书

- 日期：2026-11-22。
- 执行对象：物理开发机上的 Cursor。
- 基线：第234轮交接包，报告提交 `15a4f7c6b`，实际工作分支为 `feature/EM-Biot-Savart`。
- 本轮范围：静水初态与重置、固定温场浮力响应、电磁几何积分、诊断一致性，以及一个有限的热边界补充检查。
- 本文给出授权执行的开发和测试内容，不要求再就常规局部修复及明确列出的短测逐项确认。

## 0. 本轮目标与执行原则

本轮集中解决两个会直接影响后续复现的障碍：

1. 关闭浮力后仍有明显背景流动，使真实浮力响应难以判别；
2. 电磁解析柱的粒子积分域与连续几何域不一致，导致同采样场精度很好，却仍存在较大的几何总功率差。

另修正力链误差输出和验收标签，避免继续根据不准确的PASS/FAIL决定开发方向。

本轮不追求把所有误差消到机器精度，不增加新的大范围参数扫描。已有结果应复用，代码修复必须与具体证据对应。

### 0.1 第234轮必须保留的有效进展

| 项目 | 已有证据 | 正确结论 |
|---|---|---|
| 浮力公式 | 专用Fb与独立参考误差约4.06e-8 | 公式实现有较好证据 |
| 固定温场 | 冷轴/热缘，ON−OFF出现轴相对下沉、缘相对上升 | 浮力已产生方向合理的响应 |
| 电磁复场 | dp减半后误差约0.462%降到0.131% | 同坐标场误差呈下降趋势 |
| 匹配外层调度 | 短跑无年龄违规，G热累计1 s | 该配置的审计通过 |
| 上下边界 | local_dt_dz改善竖向抛物制造场 | 有限范围有效，不能泛化到任意二维温场 |
| 固定质量焓 | G残差约0.392% | 当前短跑规模可接受，不是首要阻塞 |

### 0.2 必须纠正的234结论

- CSV重算的力链相对误差约2.88e-4，即0.0288%，不是零。
- 233的0.163%与234的0.0288%不是同一数值；不得未经复现把历史差异全部归为浮点消减。
- `50 × expected_rel_floor`本轮约等于10.1%的相对容许误差，不能作为可靠的自动解释/通过门。
- 实际动量读取在分项笔记为PARTIAL，顶层表却为PASS，必须统一。
- 电磁几何功率差不能因被归为采样误差而直接关闭。
- 旧功率误差按解析参考归一，新指标按两者最大值归一；17.8%和15.1%不能作为改善前后比较。

## 1. A｜固定身份、参数与改动范围

### A1 开工检查

记录实际HEAD、工作树、编译器与选项、Real类型、二进制SHA256、Reload校验、运行配置。检查233/234补丁是否完整进入本机代码。

沿用用户当前分支工作安排；需要隔离时可使用独立worktree或任务分支，但不得自动切换有未提交改动的工作树，不得清理用户文件。软件机保持其已选版本。

代码提交至少分成：
1. 诊断与后处理；
2. 有证据支持的初态/重置修复；
3. EM积分定义及必要实现；
4. 边界补充检查与报告。

### A2 浮力测试固定参数

保持234客户NC几何：
- dp=0.010 m，液体高度0.180 m，液体半径按effective记录，当前参考约0.1415 m；
- rho=2700 kg/m³，mu约14.532 Pa·s，beta=1e-5 K^-1；
- 同一Reload、重力、压力路径、壁面和自由面模型；
- 使用234已验证的匹配外层调度，显式记录开关及effective路径；
- 固定温场测试关闭热源、扩散和热边界推进，仅保留规定温度对浮力的作用。

不得通过改beta、mu、声速、重力或强加速度抑制来获得更安静的图。若具体证据指向某一数值配置，先报告原因，再做单因素诊断，不展开多维组合。

### A3 不在本轮进行

- 不更换整套黎曼流体求解器，不重开δ-SPH/WCSPH方案混搭。
- 不引入新的工程壁面ghost粒子框架。
- 不改独立A-phi/Laplace研究分支。
- 不重做不锈钢弹性桨、整圈搅拌和TEAM7全套。
- 不加入μ(T)、σ(T)、导电底板、水冷壁电磁或新的skull模型。
- 不跑百秒级客户工况，不以长时间掩盖启动问题。
- 不把本轮客户几何测试称为法国500 mm自然对流复现。

## 2. B｜无浮力背景流：查清它从哪里产生

### B1 已有数据与待检验假设

234 OFF结果：
- 输出初态速度为0；
- 0.1 s速度RMS约0.0756 m/s；
- 0.5 s约0.00502 m/s；
- 1 s约0.00354 m/s。

同轮1 s的浮力ON−OFF矢量差约5.5e-5 m/s，明显小于背景。

日志还显示，静水预处理结束时Umax约0.115 m/s，随后进行了状态清理。候选原因是预处理不充分或重置后压力、密度、位置、自由面与力通道不一致，但这些只是待检验假设。

不可仅凭“没有粒子穿墙”认定静水平衡，也不可仅凭“启动后会衰减”就忽略持续背景流。

### B2 先检查源码与保存事件快照

逐项列出以下操作的实际调用次序与字段读写：
- 静水初始化与预处理；
- 密度/压力更新；
- 自由面标记和大气压处理；
- 重力、黏性及其他ForcePrior通道；
- 速度、位移、密度变化率、Force等清理；
- cell/邻域刷新；
- 时间归零及首次生产声学子步。

同一条OFF短跑至少保存：
1. 预处理结束、清理之前；
2. 每个会改变力学状态的清理阶段之后；
3. 首个生产声学子步前后；
4. 0.01、0.05、0.1 s；
5. 必要时0.5和1 s。

导出粒子稳定ID、Position、Velocity、Mass、Density、Pressure、自由面标记、压力作用项、Gravity、Viscous及ForcePrior。字段名称用真实实现；缺少分项时添加只读诊断，不制造无人使用的伪通道。

大场仅保存上述时刻；事件顺序和规约量可每步记录。

### B3 力平衡按实际离散更新核对

对内部、近壁、近底、自由面分区统计：
- 速度RMS、Umax、质心速度和动能；
- 密度偏离、压力范围；
- 实际离散压力作用与重力/其他项的合力或等效加速度；
- 每个清理操作前后的字段跳变；
- 自由面高度、质心高度及粒子代表体积变化。

必须先明确各字段是力还是加速度，以及其在声学分步中的位置。若压力通道不能直接独立导出，可在实际更新内加只读记录。不能从不同时间层取压力、质量和ForcePrior后称不平衡。

不要求边界粒子逐个满足连续静水公式；重点查是否存在区域性、持续性或由某个重置步骤突然引入的残差。

### B4 有限对照与修复路线

先复用234的OFF作为已有基线；对当前代码重现只需0.2 s。根据B2/B3证据最多选择两种候选，不遍历全部组合。

| 证据 | 候选处理 |
|---|---|
| 清理前较平衡、清理后明显失衡 | 修正具体字段重置及更新次序，保持与生产方案一致 |
| 清理前已有较大残差且仍衰减 | 有限延长无浮力预处理，按力平衡与速度共同选初态 |
| 自由面/边界更新引入压力跳变 | 定位该调用，做最小范围修复，不立即改全局边界算法 |
| 力链残留或重复加入 | 只修对应通道生命周期，重跑同一OFF |
| 初态安静但生产持续生成背景流 | 检查实际压力/密度/位置/邻域推进一致性，再选一处修复 |

预处理延长是诊断手段，不是保证成功的默认方案。若延长到5 s仍存在非衰减背景，应停止继续加时，转向已记录的力/边界残差。

禁止做法：
- 每步将速度归零；
- 在正式结果阶段增加未说明的强阻尼；
- 给ON和OFF不同的初始化；
- 清零压力或统一密度而不核查EOS、重力及自由面一致性；
- 通过改变液面约束让图“安静”，却称仍是原工程模型。

预处理阶段如采用已有松弛，必须标注并在正式观测前关闭；关闭后需要验证不会立即重新激发较大流动。

### B5 背景治理的判定

以下是本轮内部工程判据，不是文献精度标准：

1. 相同1 s观测下，OFF背景RMS若比234约3.54e-3 m/s降低至少一个数量级，可记为明确改善；仅此不能记浮力完整通过。
2. 必须同时检查启动峰值、分区残差和后段趋势，不能只挑速度最低的一帧。
3. 若未降低一个数量级，但已查明背景来源且可可靠分离浮力信号，可记PARTIAL_WITH_CAUSE；不能因此无限扫参数。
4. 只有在背景和诊断误差足够受控后，才解释ON总流场的环流形态。

本项最终交付必须指出：具体哪个操作/区域产生背景、修复改变了什么、修复后残差还剩多少。只给“延长后好一些”不算完整归因。

## 3. C｜用共同初态确认固定温场的浮力贡献

### C1 温场保持234定义

固定Eulerian温场：
`T(r,z)=1273.15 + 20*((r/R)^2-0.5) K`。

按当前粒子位置在实际浮力读取前施加，记录实际T_force_ref。温度恒定的是空间函数，不是初始粒子的温度标签。

先生成一个完整、可重用的共同初态，再分支ON/OFF。若当前重启只保存Position而不保存压力、密度、速度和相关状态，不能称“相同力学快照”；可改为同进程快照恢复或确定性的完整初始化，并核对差异。

### C2 最小运行矩阵

| 工况 | 浮力 | 其余设置 | 时间 |
|---|---|---|---|
| C_OFF | 仅关闭Boussinesq增量 | 共同初态、固定温场 | 1 s |
| C_ON | beta=1e-5 | 完全相同 | 1 s |

不能关闭总重力代替关闭浮力。ON/OFF正式观测的起点和温场施加时刻相同。

只有1 s结果显示响应正在演化、且背景不再明显污染判断时，允许将同一对照延至5 s。不要直接恢复20 s矩阵。

### C3 修正234后处理

234脚本对ON和OFF分别用r的99分位数以及各自z范围生成网格，存在网格并不完全相同的问题。本轮必须：

- 网格边界从固定几何定义一次，ON/OFF共用；
- 同一权重、同一有效支持掩码；
- 导出每个bin支持量、体积和覆盖率；
- 不将容器外、空bin或缺数据补零后参与差值；
- 轴对称投影使用实际轴心，正确计算u_r、u_z；
- 不把移动后的相同粒子ID等同于相同空间位置；
- 区分三维原始场、轴对称平均场及其差异，报告非轴对称残余。

径向bin体积不同，不用等权bin平均冒充体积加权L2。可附等权结果，但必须另名。

正式时刻读取CSV头部真实t_s，不只看文件名。匹配时间有偏差时，报告偏差并以可接受的共同时间采样，不用模糊文件名glob误选时刻。

### C4 Python输出

至少输出0、0.1、0.5、1 s的：
- 温度剖面；
- ON/OFF速度模剖面，共同色标；
- ON/OFF/差值矢量图；
- 中高处u_z(r)曲线；
- 体积加权矢量差绝对RMS及相对L2；
- 启动与后段的背景、信号演化。

ON和OFF矢量使用相同箭头缩放。差值如为可见性使用另一比例，图中明确标出，不能让独立自动缩放造成同强度错觉。

若画流线，必须保留原始矢量图。不要用经过连续性投影的流函数图替代对原始速度的核查。

### C5 验收分级

- RESPONSE_SIGNED：轴/缘差值方向正确，已有234证据，本轮核实其持续性。
- RESPONSE_RESOLVED：矢量差在共同网格、合理重分箱和必要误差估计下稳定，且与实际浮力读取对应。
- FLOW_BUOYANCY_DOMINATED：总流场主要环流不再由OFF背景解释。
- FRENCH_TOPOLOGY_MATCH：另行文献专用配置验证，本轮不授予。

不得把前三者合并成一个“自然对流全部PASS”。若响应接近测量底，仅做一次有针对性的分箱稳健性或OFF重复检查；不无条件重跑所有工况。

## 4. D｜纠正力链误差规约与过宽自动判据

### D1 复算同一份234 CSV

先使用原始数据做独立double后处理，结果应能复现：
- 相对误差约2.88e-4；
- 最大z向绝对残差约9.31e-10 N；
- Fp_pre与Fp_off在该CSV中一致。

上述值是审核基准，不是需要硬编码达到的答案。全矢量字段若比CSV z分量更完整，要单独说明。

对比运行内规约、CSV规约、summary、日志和退出码，查明为什么运行结果为0。检查类型转换时机、平方与累加精度、编译浮点选项、实际运行版本和同步；目前不预设某一个因素是根因。

特别注意：将Float的squaredNorm结果再转换到double，不等同于先将每个分量转double再平方累加。诊断应明确每一步精度。

### D2 修正验收逻辑

删除或禁用“50倍全局估计舍入尺度即可自动解释”的正式验收用途。保留估计值作为诊断信息，不把它当作已证明的误差地板；实际舍入误差不必恰好达到该尺度。

采用以下独立证据：
1. 专用浮力与独立公式的误差；
2. 同快照加入前后/ON−OFF的实际绝对残差；
3. 按实际操作与力尺度估计的舍入范围；
4. 生命周期与实际动量读取；
5. 同状态重复调用无额外累加。

若设置自动容限，必须基于具体算术操作的误差界或明确验证过的数值容限，报告其推导，不能另换一个任意放大系数。

接近零Fb的粒子用绝对误差。禁止把小分母加大常数后自动得到漂亮相对误差。

### D3 实际动量读取的最小验证

不再把“完整声学Δv与ForcePrior*dt不同但不为零”当作读取已证实。优先在真实更新位置记录实际读取的ForcePrior及其对更新式的贡献，并核对实际单位、质量换算和子步系数。

可用同一冻结状态的一个声学子步ON/OFF对照辅助判断，但必须恢复所有会被子步修改的状态。包含压力反馈后的速度差不能强行要求逐点等于a_b*dt。

分项笔记、顶层报告、验收表使用同一作用域和状态。仍只完成通道验证时写CHANNEL_VERIFIED，动量读取保持PARTIAL。

本项只做快照与短子步，不重开整轮长时浮力测试。

## 5. E｜电磁：修正积分域，保留已有场精度成果

### E1 固定已验证的场比较

继续使用234解析柱：
R=0.10 m，总高0.60 m，B0=0.01 T，sigma=16 S/m，f=300 kHz；
dp=0.025和0.0125 m，严格迭代残差约1e-3或更低。

保留同坐标、同权重、未经幅值拟合的复数场误差。现有结果分别约0.462%、0.131%，不因几何功率问题否定这项成果。

本轮优先离线核对积分，不继续加密，不继续主要调Picard容差，不增加柱长。

### E2 从已有粒子数据检查几何矩

对全柱与中部切片分别比较：

`V_discrete=sum(V_i)`
`V_geometry=pi*R^2*H_region`
`M_r2_discrete=sum(V_i*r_i^2)`
`M_r2_geometry=pi*R^4*H_region/2`

常数场检查域体积，r²场检查与弱趋肤焦耳热分布相关的径向积分矩。另用实际解析Q积分，不能用r²近似替代本工况完整参考。

234中部示例：
- 解析参考厚度约0.09 m，必须从实际参数确认；
- 粒子代表体积约0.003250 m³；
- 对应连续几何体积约0.002827 m³；
- 体积相对几何约大14.95%，两档几乎相同。

逐项核对：
- 轴心与柱体z上下限；
- 切片上下界及是否包含端点；
- 选中粒子中心层数及每层代表厚度；
- 径向筛选与粒子代表体积是否跨出R；
- 中部切片边缘粒子的体积是否需要截断；
- 总高/半高是否混用；
- V的来源及是否与实际粒子布置一致。

不能仅将参考高度改成“选中层数×dp”后称原0.09 m切片通过。改变区域可作为另一个显式比较，但原定义的误差需保留。

### E3 先独立积分，再决定是否改生产数据

推荐两步：

1. 建立独立的连续几何积分参考，检查解析Q径向积分的数值收敛。
2. 对规则粒子格点，离线估计其代表单元与圆柱/切片的交体积，构造几何截断积分权重；可用可收敛的子单元积分，不要求引入新CAD库。

若实际粒子布局不是规则格点，不能套用dp立方体交体积；先明确代表体积定义，再选与其一致的独立几何求积。

必须区分：
- 原生产体积权重；
- 仅用于几何域验证的截断求积权重；
- 解析Q在粒子中心的近似；
- 解析Q在代表单元内的积分。

截断权重仍可能保留场内采样误差，不能保证粗分辨率立刻机器精度闭合。用求积细化或解析几何矩检查其可信度。

### E4 禁止用重归一掩盖问题

不得：
- 为让总功率正确而整体乘一个经验比例；
- 改目标功率、电流、sigma或B0来抵消积分误差；
- 未检查影响就把所有流体Mass/Vol改成修正值；
- 将后处理截断权重自动用于EM源/自感离散；
- 删除不利中部区域或只报同采样功率。

如果发现生产粒子体积本身与几何域代表关系有系统问题，应单独报告其对EM源积分、焦耳热和能量传递的影响，另做最小修复与回归。不能把诊断权重修复当作生产积分已修复。

### E5 统一功率误差口径

主报告统一使用：
`e_ref=(P_num-P_ref)/P_ref`
同时给有符号差、绝对差和绝对相对误差。接近零参考另给绝对值。

原max分母指标可保留但改名，禁止与上述指标混用比较。

每档分辨率分别列：
1. 原数值离散功率；
2. 同采样解析离散功率；
3. 连续几何解析功率；
4. 几何截断/修正求积后的数值功率；
5. 相同修正权重下解析功率；
6. 原/修正体积与r²矩。

明确peak/RMS及1/2约定保持不变。旧17.8%与新15.1%的差异需用同一组功率数值、两个分母直接解释。

### E6 本轮验收

- GEOMETRY_MOMENTS_EXPLAINED：体积和r²矩差异能由域、选层和权重解释；
- GEOMETRIC_QUADRATURE_VERIFIED：独立几何求积经过自收敛和已知积分检查；
- FIELD_ACCURACY：保留复场指标及分辨率趋势；
- PHYSICAL_POWER_ACCURACY：必须对同一连续几何域、固定参考口径评价，不能由同采样闭合替代。

内部目标：独立几何求积的常数/r²矩误差优先做到≤1%；若用于判断百分之一级功率差，应进一步证明求积不确定度足够小。此目标是求积检查，不是修改场数据使其达标。

若原功率差被修正求积显著消除，可关闭相应后处理积分缺口，并清楚标注生产权重是否改变。若仍有较大残差，直接用已拆分结果确定剩余误差，禁止再泛称“采样问题”。

TEAM7和局部峰值Q仍为独立遗留，本轮不因此升级为EM全模块验证完成。

## 6. F｜热边界只补一个交叉检查，避免喧宾夺主

### F1 认识当前local_dt_dz范围

234实现是整个底部带一个斜率、整个顶部带一个斜率，不是每个表面位置的局部重构。不得在软件或报告中将其描述为一般二维局部边界算法。

保留已通过的uniform、linear_z、radial_quad、parabolic_z结果。

### F2 一个同时随r、z变化的制造场

建议：
`T(r,z)=1200 + 100*(1+0.5*(r/R)^2)*(2*(z-z_bot)/H-1)^2 K`。

这是诊断场，不是法国文献温度。其竖向梯度随半径变化，可检查每面单斜率的局限。表面参考通量独立积分，辐射使用局部T^4。

仅检查已有候选local_dt_dz，必要时附band对照。沿用逐面功率≤2%、温度更新能量≤1%的内部门限，并同时报告连续误差。

若失败，先记录误差与适用范围。本轮不因此全面开发逐点局部边界，除非已明确是极小的实现缺陷且可直接修复。

### F3 修正拟合有效性

当前只以粒子数量判断ok不足。补检查：
- 加权z方差/拟合分母；
- 有效高度层与条件性；
- 温度和斜率有限性；
- 回退次数、原因及模式。

拟合退化后返回0，不能同时标ok=1。回退必须显式记录，不能静默通过。

固定质量焓继续使用约0.4%的真实口径并监测，不将其列为本轮首要精度攻关。

## 7. G｜何时恢复热演化与真实EM耦合

本轮优先完成B/C与E。满足下列条件后，才做一次规定热回归：
- OFF背景来源已经明确，改善或剩余限制有量化；
- 固定温场响应可分辨，力链诊断与字段一致；
- 匹配外层调度和热时间记录正确；
- 热边界模式身份明确。

回归仍采用10 kW规定热，EM明确NOT_RUN，先1 s，必要时5 s。只选一个已固定配置，不展开多个边界/材料组合。

真实EM热源留到积分问题明确之后。若只是离线积分修复而未影响生产热源，必须说明现有生产热源还需要什么检查，不能自动接入后称全链通过。

法国文献几何/温场与客户热电偶定量需求分别记录。当前sigma=16仍是借用值；30 kW电源输入不等于玻璃吸收功率。不得通过设置初始1000℃宣称测点验证。

## 8. 运行顺序、预算与停止条件

### 8.1 推荐顺序

1. 身份核对及D力链离线重算；
2. B初始化事件审查、OFF短测；
3. E电磁已有数据离线积分，可独立推进；
4. B只选证据支持的一至两种修复；
5. C共同初态ON/OFF；
6. F一个边界交叉场；
7. 满足条件才做G规定热短回归；
8. 汇总报告、补丁与软件接收清单。

### 8.2 默认预算

| 工作 | 默认预算 |
|---|---|
| 力链一致性 | 既有CSV离线复算，加一次快照/单子步验证 |
| OFF基线复现 | 0.2 s，已有1 s结果复用 |
| 背景候选 | 最多两种，各先0.2 s；选定者到1 s |
| 延长预处理 | 有证据时一次，最高5 s；非衰减则停止 |
| 修复后固定T | ON/OFF各1 s；必要时延至5 s |
| EM积分 | 复用两档现有场，优先零新增EM求解 |
| 边界交叉场 | 一个静态制造场 |
| 规定热回归 | 有条件的单配置1–5 s |

只有代码修改实际影响EM场时，才重跑相关两档；不能因为报告更新而机械重跑全部。

### 8.3 停止规则

- 无明确证据，不增加第三种背景修复组合。
- 背景不衰减，不继续单纯延长预处理或正式运行。
- EM体积/几何矩尚不一致，不继续主要调求解器容差。
- 只是指标规约错误，不改浮力公式。
- 只完成离线求积，不称生产EM积分已修。
- 边界交叉场失败，记录限制，不让支线阻塞背景主因定位。
- 时间不足时优先交付已查明原因和可复核数据，不以伪PASS收口。

## 9. 目录、原始证据与软件移交

统一产物根：
`docs/ophelie/reference/doc235_runs/`。

建议子目录：
- A_identity；
- B_off_background；
- C_fixed_T_response；
- D_force_gate；
- E_em_geometry；
- F_boundary_crosscheck；
- G_thermal_regression；
- H_handoff。

必须交付：
1. ROUND235_REVIEW.md与ACCEPTANCE_TABLE.csv；
2. 真实执行命令、参数、退出码、版本与文件校验；
3. 初态/清理/首次子步事件表、关键粒子快照、力平衡统计；
4. 修复前后OFF背景与启动峰值曲线；
5. ON/OFF同网格场、支持权重、CSV头真实时刻及Python图；
6. 力链CSV独立重算脚本、规约修复、判据推导和一致性检查；
7. EM体积/r²矩/三种功率/修正求积对照，独立积分脚本及误差预算；
8. 边界交叉温场逐面功率与拟合有效性统计；
9. 源码补丁，包含实际修改，不只给本机路径；
10. SOFTWARE_HANDOFF_DELTA.md。

Python脚本与报告JSON必须对应同一实际运行版本。234包中脚本与结果字段并不完全对应，本轮禁止把后来临时生成的摘要与旧脚本混装。

运行输出按工况隔离，避免同一个bin/output被后续运行覆盖。大场可外置，但审核必需的稀疏场、CSV、脚本、补丁和日志片段直接打包。

### 9.1 软件机可以选择接收

- 准确命名的误差、能量与运行身份字段；
- 经验证的匹配调度路径及配置说明；
- 已完成回归的具体初态/重置修复；
- 明确区分原权重与几何验证求积的报告功能；
- 边界有效性与回退信息。

### 9.2 不得直接移交为默认

- 未完成实际动量读取验证的“完整浮力PASS”；
- 高达约10%的浮点解释自动门；
- 任意几何均适用的local_dt_dz声明；
- 为凑功率而重归一生产质量/体积；
- 法国双环流、客户1000℃匹配、TEAM7或峰值Q已通过的标签。

## 10. 本轮收口必须回答的五个问题

1. OFF流动主要在预处理、状态清理还是生产推进的哪个环节产生？证据是什么？
2. 修复后背景下降多少，浮力响应是否稳定可辨，是否已能解释总环流？
3. 力链CSV与运行规约为什么不一致，新的判据是否避免误放行？
4. 电磁几何功率缺口中多少被域/体积修正解释，生产积分是否仍受影响？
5. 下一轮是否能进入有明确温场与热源身份的自然对流推进？若不能，只剩哪个具体阻塞？

本轮完成标准是查明并处理上述因果关系，使后续开发进入明确的热源—温度—流动验证；不以增加测试数量或扩大PASS标签作为进展。


<!-- ===== END 235_OPHELIE_OFF_FLOW_EM_GEOMETRY_AND_DIAGNOSTIC_CLOSURE_CURSOR_PLAN_2026-11-22.md ===== -->


<!-- ===== BEGIN 236_OPHELIE_HYDROSTATIC_SUPPORT_AND_EM_GEOMETRIC_QUADRATURE_CURSOR_PLAN_2026-11-23.md ===== -->

# 236｜OPHELIE 静水压力支撑修复与电磁几何功率验证执行计划

- 日期：2026-11-23。
- 执行对象：物理开发机 Cursor。
- 输入基线：第235轮交接包。报告提交为 `e4c4ad963`，同时包含未提交WIP；必须先核实实际源码与二进制，不能仅凭HEAD代表全部结果。
- 当前分支：尊重用户现有 `feature/EM-Biot-Savart` 工作安排，不自动切换或清理工作树。
- 本轮编号：236。本文替代235中尚未完成事项的执行安排，保留235历史记录，不改写旧报告。
- 本轮主目标：完成静水启动的具体归因/最小修复，以及电磁连续几何域的独立功率验证。
- 本轮暂缓：新增规定热演化短跑、真实EM—热—流耦合、长时计算。

## 0. 必须先理解的任务范围

本轮不是再增加一轮泛化审计，也不是换求解器。已有证据表明浮力公式具有较好精度、固定温场有方向合理的响应、电磁解析场随加密改善。

当前两个具体障碍：
1. 无浮力时存在强的整体上下运动，影响总环流解释；
2. 电磁同采样场误差小，但粒子积分与连续几何域未匹配，功率验证尚未闭合。

不能以“已经解释了原因”“程序passed=1”代替完成修复与验证。不能把未实施的截断求积再次仅列为下一轮。

### 0.1 本轮审核基准

| 事实 | 已有数值/证据 | 不能推出的结论 |
|---|---|---|
| 235 fix2预处理结束 | 平均u_z约−0.0446 m/s，RMS约0.0550 m/s | 预处理已达静水平衡 |
| 235 fix2首声学步 | 平均u_z约−0.00329 m/s | 清零速度后的状态已平衡 |
| 约0.01 s | 平均u_z约−0.0739 m/s，RMS约0.0790 m/s | 只是少量局部数值涡 |
| 约0.2 s | 平均u_z约+0.0264 m/s | 流动已经是温差驱动 |
| 0.01 s平均运动动能占比 | 约87.5% | 背景主要是局部环流 |
| 力链double复算 | 相对误差约2.88e-4 | 浮力完全没有加入 |
| 电磁中部体积 | 比连续几何域大约14.95% | 将功率差改名即可通过 |
| 边界交叉温场 | local_dt_dz逐面误差约≤1.2% | 已成为任意二维局部边界算法 |

“清理后速度为零”不证明清理未破坏平衡；它只表示没有立即写入非零速度。压力、密度、位置和力通道仍可能不相容。

### 0.2 本轮不可变约束

- 不调大beta、不降低mu，不用经验电流/功率比例掩盖误差。
- 不修改独立A-phi/Laplace研究分支。
- 不换整个黎曼流体框架，不重新混合多篇论文算子。
- 不新建工程壁面ghost粒子布置框架。
- 不在正式计算阶段每步归零速度或添加未说明的强阻尼。
- 不将离线几何求积权重自动用于生产Mass、Vol或EM源积分。
- 不重做搅拌桨、TEAM7全套、温变物性、导电底板等支线。
- 不要求软件机整体合并物理分支；移交必须择收并有回归。

## 1. A｜冻结版本与完成最小的初始化解耦

### A1 冻结真实基线

1. 核对235已提交代码、WIP补丁和实际运行二进制。
2. 保留用户文件，不执行清理或丢弃改动。
3. 将本轮实际使用的物理代码形成可追溯提交或明确的HEAD+补丁校验；优先正式提交功能完整的基线。
4. 记录构建选项、Real类型、设备、二进制SHA256、Reload、有效配置。
5. 分开提交：初始化解耦、诊断、必要修复、EM求积、报告。不要把未经验证的默认行为与诊断混在一起。

### A2 消除热源开关控制力学启动的隐含差异

235源码在 `paper_flat_slip_top || zero_q` 条件下执行特定清理，又依据zero_q决定是否刷新密度。固定温场与非零热源可能因此进入不同机械初始化路径。

要求：
- 将初始化/清理行为整理为明确的力学策略；
- 热源模式只决定热源，浮力开关只决定Boussinesq增量；
- 保留旧行为为有名的兼容策略，便于复现235；
- 不在解耦时悄悄更换密度、压力、自由面或黏性方案。

可增加一个case局部参数，例如 `mechanics_initialization_mode`；具体名称由现有CLI风格确定。本文建议名称不是已经存在的命令，必须实现、帮助说明和配置落盘后才能在RUN_COMMANDS中使用。

### A3 无需推进的初始化一致性检查

使用相同初温、Reload、力学策略、预处理和浮力关闭状态，分别选择zero_q与规定热模式，在正式热更新之前比较：
- Position、Velocity、Density、Pressure；
- 重力、黏性及ForcePrior状态；
- 自由面标记、邻域更新计数；
- 物理时间、外层预算和首步dt。

此检查只到t0，不执行热演化。非零源的Q可不同，热更新之前的力学状态不应因热源标签而改变。

记录逐字段差异与来源。GPU规约小量差异按实际精度解释，不能只给初始化相同的布尔值。

## 2. B｜静水支撑：查真正的不平衡力，不再只查速度

### B1 固定工况

使用235客户NC几何：
- dp=0.010 m，H=0.180 m，R取实际effective，约0.1415 m；
- rho=2700 kg/m³，mu约14.532 Pa·s；
- beta=0，仅关闭浮力增量，保留真实重力；
- 温度均匀固定，热源/扩散/换热温度推进关闭；
- 匹配外层调度，壁面、自由面、压力路径与235相同；
- 不因无热源而隐式改变力学初始化。

本项是静水保持测试，不是自然对流测试。

### B2 先用已有数据确认整体运动

离线复算235三个候选的：
- 质量加权平均速度向量；
- 总动能；
- 平均平移动能 `0.5*M*|U_mean|^2`；
- 扣除平均平移后的剩余动能；
- 质心位置与速度随时间；
- 液面高度和体积。

“扣除平均速度”仅用于诊断分解，不得把这样处理后的速度作为正式流场，或用来抹掉整体动量问题。

预计能复现fix2在0.01 s约87.5%的平均运动动能占比。不能把这一结果硬编码为答案。

### B3 在真实动量更新位置增加只读诊断

先读清实际分步更新式，列出：
- 压力力/加速度从何处产生；
- Gravity与黏性等ForcePrior的单位和时间层；
- 壁面压力作用包含在哪一通道；
- 速度与位置分别在哪个子步更新；
- 黎曼耗散/边界修正是否另有贡献；
- 实际力到加速度的质量换算及dt系数。

在实际读取位置记录压力、重力、黏性及其他有效贡献，特别是声学第一/第二半步。不要把第二半步结束后名为Force的字段直接当作第一半步压力力。

如果不能以字段直接拆分，可在实际更新函数加诊断记录；通过case选项启用，正式计算关闭时不改变结果。不得新建一个无人读取的力变量后用它证明生产链正确。

### B4 必须做全局竖向动量预算

由真实离散更新构造：
`DeltaPz_observed = sum(m_i * (u_z_after-u_z_before))`。

独立累计本步各通道实际积分的竖向冲量，核对：
`DeltaPz_observed = I_gravity + I_pressure_inner + I_pressure_wall + I_viscous + I_other + residual`。

说明分步时间系数、固定质量假设和其他直接速度修正。若某个算子直接写速度，也必须进入预算，不能当作无来源残差。

内部压力作用不预设严格成对守恒；先按当前离散实测其全局贡献。壁面作用不能遗漏。若压力内外无法完全拆分，先闭合总压力贡献，再局部分解。

静水情况下重力约为M*g，当前M约30.3264 kg，对应约297.5 N。该数只用于数量级核查。正式参考由实际Mass和gravity求得。

目标：解释首步和0.01 s整体向下动量从何而来，而不是仅给相对连续静水公式的误差。

### B5 空间分区与阶段

分区：深内部、侧壁、底部、自由面、底角/顶角。使用几何带或现有标记，明确带宽和粒子重叠关系。

至少记录：
1. 预处理末尾、清理前；
2. 速度归零后；
3. 其他机械字段清理后；
4. 正式t0；
5. 首步真实读入力与子步结束；
6. 0.01、0.05、0.1、0.2 s。

保存Pressure、Density、Position、Mass、Velocity、自由面标记和力贡献。常规压力剖面可与连续静水尺度比较，但离散动量预算是主要证据。

若某次重置使合力明显改变，指出具体字段；若t0已不平衡但重置没有改善，不能说“清理未引入速度，所以清理无关”。

## 3. C｜只选择一个有证据的修复路径

### C1 先完成B，再选择

不得同时尝试一堆边界/密度/压力组合。根据预算选择最强证据对应的一条路径：

| 诊断结果 | 最小处理 |
|---|---|
| 已有完整状态较平衡，特定清理后失衡 | 修正该清理的字段生命周期，保持生产方案一致 |
| 预处理阶段整体仍振荡、支撑逐渐趋稳 | 一次有限延长预处理；保持完整状态进入观测，避免重新破坏平衡 |
| 离散壁面压力支撑明显不足/方向错误 | 定位具体壁面调用、法向/状态/邻域或系数，只修实证缺陷 |
| 内部压力/密度更新引入明显非物理净力 | 检查实际时间层与邻域一致性，不先更换离散框架 |
| ForcePrior残留、漏读或重复 | 只修通道生命周期与读取 |
| 有直接速度修正产生未计冲量 | 显式纳入预算并检查作用范围 |

保持完整状态不等于把较大残余速度称为静止。必须量化残余与后续趋势，避免仅保留一个“比较好看的时刻”。

### C2 一次有限的预处理延长对照

只有预算和已有轨迹支持“未充分松弛”时，授权一次延长：
- 先到1 s；
- 若动量残差、动能包络均明确衰减，可延至最高5 s；
- 若持续振荡或残差不下降，停止加时，转向B定位到的支撑缺口。

松弛阶段如使用既有阻尼，必须记录，并在正式保持前关闭；关闭后观察是否反弹。不得在正式结果阶段靠新增阻尼压住流动。

不要每隔若干步反复归零速度来伪造平衡。预处理和正式阶段的物理时间分别记录。

### C3 修复后同配置静水保持

先0.2 s，若预算和启动明显改善，延至1 s。比较：
- 早期峰值；
- 0.5–1 s后段的平均速度、RMS、剩余动能；
- 全局动量预算及分区支撑；
- 密度、液面与体积趋势。

内部改善目标：1 s OFF RMS相对235约3.54e-3 m/s至少下降一个数量级，同时整体竖向运动不再反复主导。这是改善目标，不是精度认证。

若未达到，但根因和有限修复已明确，记录实际差距；不要调低判据后写完全通过。若仍无法定位，列出预算中哪个未解释通道阻塞，而不是泛称自由面噪声。

## 4. D｜补齐真实浮力动量读取，并有限确认响应

### D1 用B的只读诊断完成此前未完成的读取验证

对固定状态记录专用Fb、ForcePrior加入前后，以及真实动量更新读取的该项贡献。核对只加入一次、没有被覆盖、质量换算和子步系数正确。

保留235已修正的double分量提升后规约。禁止恢复50倍全局舍入尺度放行。

现有0.1%通道相对门限应说明用途与依据；它不能独立证明动量读取。接近零浮力粒子继续使用绝对误差。

验收字段分开：
- FORCE_FORMULA_VERIFIED；
- FORCE_CHANNEL_VERIFIED；
- FORCE_MOMENTUM_READ_VERIFIED。

实际读取证据未完成时，第三项继续PARTIAL，不以常量true或程序退出0替代。

### D2 是否做固定温场ON/OFF

只有C的静水保持改善、或B已建立可可靠分离背景的证据后，再跑一对1 s固定温场：
`T=1273.15+20*((r/R)^2-0.5) K`。

使用相同完整初态/机械初始化，温场施加时刻相同；OFF仅关闭Boussinesq。初态比较包括压力、密度及相关状态，不只比较位置和速度。

继续使用共同几何网格、有效支持、体积权重和真实输出时刻；Python出ON、OFF、差值矢量及u_z(r)。区分方向正确、响应可分辨、总流场由浮力主导三个状态。

若静水背景仍未受控，本轮不必重复这组试验；保留235方向证据，不把它降级成浮力公式错误。

## 5. E｜电磁截断求积：本轮必须实施，不能只重复体积差

### E1 固定物理输入与已知场结果

沿用234/235解析柱：
R=0.10 m，总高0.60 m，B0=0.01 T，sigma=16 S/m，f=300 kHz。
两档dp=0.025和0.0125 m，真实迭代残差相近。

保留：
- 同坐标复场误差约0.462%与0.131%；
- 同采样功率差远小于几何域功率差；
- 中部全粒子体积相对0.09 m参考切片约大14.95%。

不再主要调迭代容差，不增加第三档分辨率，不加柱长。

### E2 首先固定实际几何域

由effective/源码确认：
- 柱轴心、上下界；
- 中部切片上下界，参考厚度是否确为0.09 m；
- 粒子中心筛选条件与端点；
- 每个粒子的代表单元与体积来源。

主指标固定：
`e_ref=(P_num-P_ref)/P_ref`。
同时给绝对W、有符号差和绝对相对误差。max分母指标可附加，但必须另名。

不能把参考高度改成选中层数乘dp，然后称原切片通过。

### E3 独立连续几何参考

用独立脚本计算：
- 体积：pi*R²*H_region；
- r²矩：pi*R⁴*H_region/2；
- 完整趋肤解析Q的连续积分：沿r积分2*pi*r*Q_ref(r)，再乘对应高度。

明确解析复幅值、sigma与1/2约定。用求积细化检查独立参考收敛；参考积分相邻精度档变化优先控制在0.1%以内。该数值是本轮求积预算，不是修改物理结果的理由。

用常数、r²及完整解析Q分别检查，不能只用最容易闭合的常数场。

### E4 实施粒子代表单元与几何域的交体积

对规则格点粒子，可按实际生成方式重建代表单元，计算其与圆柱/切片的交体积。可用可细化的子单元采样或解析截断，不要求引入CAD库。

必须先检查一个容易漏掉的问题：

> 只保留中心在圆柱内的粒子时，这些代表单元的并集未必完整覆盖真实圆柱。仅截掉伸出域外的部分，可能同时留下边界缺口。

因此要求：
- 导出单元并集覆盖情况；
- 量化未覆盖、重叠和超出体积；
- 不把“所有截断权重之和”天然当成完整几何体积；
- 不能为补缺口而整体归一权重。

若只有中部CSV，不足以计算跨切片单元或完整几何求积，先查既有完整产物。确实缺数据时，允许仅重跑原两档EM目标以导出Position、V、复J向量、Q和输入元数据；这属于补齐必要数据，不是新增参数扫描。

### E5 完整域验证求积：覆盖不足时的明确处理

优先保留E4的原粒子权重对照。如其并集不覆盖域，则额外建立覆盖完整圆柱的独立积分点集：
- 在完整几何域生成可细化求积点；
- 从数值场重构复J或Q，明确方法、支持和边界外推；
- 不新增求解器粒子，不赋予虚构物理质量；
- 不用解析解给缺失的数值场直接补值。

重构可采用局部最小二乘或既有可靠插值，但必须检查条件性与支持，并用独立制造场核对其积分误差。J向量重构与J_theta标量重构要明确轴心处理；相位不能丢失。

如果重构不可靠或边界支持不足，报告该区域与求积不确定度，保留几何功率验证OPEN。不能因时间限制用解析值填空后宣布数值解通过。

本轮必须至少完成E3/E4的实际求积和覆盖检查；若需E5，实施一个有限、可复核的完整域方案，或明确给出无法完成的具体数据/支持阻塞。不能再只交原14.95%复算。

### E6 对照表必须包含

针对全柱和中部切片、两档分辨率，分别列：
1. 原生产粒子权重的数值功率；
2. 同原权重的解析采样功率；
3. 连续几何解析功率；
4. 截断交体积下数值功率；
5. 相同截断权重下解析功率；
6. 如实施，完整域独立求积的数值与解析功率；
7. 原/截断/完整域的体积与r²矩、覆盖缺口；
8. 场误差、插值误差与求积不确定度。

单元内Q在中心近似与单元内积分不是同一件事；报告各自误差来源。粗档没有立刻达到百分之一级，不等于求积方案无效，要看误差拆分和细档趋势。

### E7 验收及对生产积分的影响

分开判定：
- FIELD_ACCURACY：既有复场精度；
- GEOMETRY_COVERAGE：域覆盖是否完整；
- GEOMETRIC_QUADRATURE：独立求积是否可信；
- GEOMETRIC_POWER_ACCURACY：同一连续域总功率差；
- PRODUCTION_INTEGRAL_STATUS：生产权重是否改变、是否仍有偏差。

几何求积的常数与r²矩优先做到≤1%；涉及百分之一级功率判断时，进一步记录其不确定度。物理功率误差如未达到原项目期望，原样保留，不修改容限凑PASS。

离线求积验证通过不自动等于生产源积分修复。若原Mass/Vol影响自感源或热输入，单独列出影响链与下一步所需修复；不擅自统一修改生产质量/体积。

TEAM7与局部峰值Q仍未被解析柱替代验证。

## 6. F｜本轮暂停与保留事项

### F1 暂停新增热耦合运行

不重复235的规定热1 s smoke，也不新增真实EM—热—流短跑。它们不能解决本轮的静水支撑与积分缺口。

这不表示停止软件开发。软件机继续既定GUI与封装；物理标签保持真实状态。

### F2 保留已有热边界结果

local_dt_dz的uniform、linear_z、radial_quad、parabolic_z、rz_quad结果及拟合退化检查保留。本轮不追加更多制造场，不改成通用逐点二维边界。

固定质量焓约0.4%的短跑误差继续记录，不作为本轮首要攻关。

### F3 恢复耦合的条件

本轮完成后再决定：
- 静水整体动量预算明确；
- 初始化与热源标签解耦；
- 浮力实际读取有证据；
- EM几何功率误差与生产积分影响有明确结论。

满足后下一轮先恢复规定热的温度—流动，再按独立门接入真实EM。客户几何/测点定量与法国文献环流拓扑分开评价。

## 7. 执行顺序、预算与停止条件

### 7.1 执行顺序

1. A：冻结WIP身份，完成初始化解耦及t0对照。
2. B：用已有事件数据复算整体动量；添加实际读入力与预算。
3. E：复用EM现有数据推进独立几何求积，可与流体诊断分目录执行。
4. C：只选择一个有证据的静水修复路线。
5. D：补齐浮力读取；有条件做一对1 s固定温场。
6. 汇总证据与软件移交建议，不重跑热耦合。

### 7.2 运行预算

| 任务 | 本轮预算 |
|---|---|
| 初始化解耦检查 | t0，不做热推进 |
| 静水力预算 | 一个基线0.2 s |
| 最小修复 | 一个主要候选，先0.2 s，必要时到1 s |
| 延长预处理 | 有证据的一次，先1 s，明确衰减才到5 s |
| 浮力读取 | 快照/单子步 |
| 固定温场响应 | 条件满足才ON/OFF各1 s |
| EM | 优先零新增求解；缺完整场时补导出原两档 |
| 热耦合/其他case | 本轮不新增 |

### 7.3 停止条件

- 预算未闭合，不用速度图替代因果解释。
- 一条修复无效，回到已记录的通道预算，不自动扫描第三、第四种组合。
- 预处理不衰减，停止加时。
- 求积覆盖不完整，停止用截断体积直接代表连续域。
- 只是离线权重改进，不宣称生产EM积分已修。
- 一项阻塞不影响另一项独立工作时，继续完成另一项；不把整轮停在缺少甲方电流上。

## 8. 产物与交接要求

产物根：`docs/ophelie/reference/doc236_runs/`。

建议目录：
- A_identity_and_init；
- B_hydro_momentum；
- C_hydro_fix；
- D_buoyancy_read；
- E_em_quadrature；
- F_handoff。

必须包含：
1. ROUND236_REVIEW.md、ACCEPTANCE_TABLE.csv、RUN_COMMANDS.md；
2. 初始化策略与热源解耦的源码补丁及t0逐字段对照；
3. 真实动量更新式说明、调用位置与每通道冲量CSV；
4. 预处理/清理/首步事件数据，以及整体/分区动量预算；
5. 修复前后OFF质心速度、动能、液面及预算图；
6. 浮力实际读取证据，未完成时明确PARTIAL；
7. EM几何参考、截断权重、覆盖检查、求积收敛和功率对照脚本；
8. 必要的原始EM场与稳定索引，不能只给已有库外路径；
9. 所有新增脚本与摘要JSON对应同一代码版本；
10. 软件机可择收的补丁及不可宣称的物理结果。

原始稀疏场、关键CSV、日志片段和补丁直接入包。大场可外置并给校验，但不能只提供审核者无法访问的本机路径。

### 8.1 验收表建议项目

- INIT_INDEPENDENT_OF_HEAT_MODE；
- HYDRO_MOMENTUM_BUDGET_CLOSED；
- OFF_BACKGROUND_CAUSE_IDENTIFIED；
- OFF_BACKGROUND_REDUCED；
- FORCE_MOMENTUM_READ_VERIFIED；
- FIXED_T_RESPONSE_STATUS；
- EM_GEOMETRY_COVERAGE；
- EM_GEOMETRIC_QUADRATURE_VERIFIED；
- EM_GEOMETRIC_POWER_ACCURACY；
- PRODUCTION_EM_INTEGRAL_STATUS。

标签可用PASS、PARTIAL、FAIL、OPEN、NOT_RUN，必须附作用域和证据。程序smoke_passed不能覆盖这些物理项。

## 9. 最终收口报告必须回答

1. 启动时整体向下的净冲量主要由哪个离散通道/状态不平衡造成？
2. 清理、预处理与正式推进各自贡献是什么？修改后改善多少？
3. zero_q与非零热源是否已经共用同一力学初始化？
4. 专用浮力是否被真实动量更新读取，证据在哪里？
5. EM截断/完整域求积实际实现了什么，功率缺口减少多少，未覆盖区域还有多少？
6. 生产EM积分是否仍需修复？下一轮应继续哪个具体问题，还是可以恢复热耦合？

本轮成功标准是完成具体因果链和实际求积，不是再写一份“基本解释清楚，下一步继续”的总结。不能完成时要指出明确的未完成步骤、原因和已有证据，避免误判项目已收口。


<!-- ===== END 236_OPHELIE_HYDROSTATIC_SUPPORT_AND_EM_GEOMETRIC_QUADRATURE_CURSOR_PLAN_2026-11-23.md ===== -->


<!-- ===== BEGIN 237_OPHELIE_LOCAL_BACKGROUND_AND_EM_QUADRATURE_CLOSURE_CURSOR_PLAN_2026-11-23.md ===== -->

# 237｜OPHELIE 局部背景流定位与电磁几何求积收口开发任务书

- 日期：2026-11-23。
- 执行对象：物理开发机 Cursor。
- 本轮输入：第236轮交接包；报告基线 `845e580de` 加未提交物理WIP。先冻结真实版本，不能将该HEAD直接等同于本轮全部运行代码。
- 当前分支沿用用户安排，不自动切换有未提交改动的工作树。
- 本轮目标：定位剩余局部背景流、完成可信的全柱与中部几何功率验证，并给出是否恢复热演化的明确决策。
- 本文是237轮执行计划；236历史记录保留。

## 0. 开工先读：哪些已解决，哪些不能再重复

### 0.1 保留236的有效成果

1. 延长静水预处理至1 s后，启动阶段整体下沉/回弹显著改善。
2. 当前几何中，首步Force与ForcePrior的反向抵消比例从约0.295改善到约1.036。
3. 0.01 s时整体平均运动动能占比从约87.5%降到约2.1%。
4. 浮力公式、专用通道、真实声学第一半步动量读取已有证据；读取相对残差约5.05e-4。
5. 热边界local_dt_dz的有限制造场验证已通过，暂不扩展支线。
6. 电磁同坐标复场误差从粗档约0.462%降到细档约0.131%，保留该成果。

本轮不再重开浮力公式和力链单元测试全套。只有修改实际影响该链路时，才做针对性回归。

### 0.2 当前两个主缺口

| 缺口 | 现状 | 本轮必须完成 |
|---|---|---|
| 局部无浮力背景流 | 1 s OFF RMS仍约3.41e-3 m/s；整体平移已弱，但局部运动仍强于浮力差 | 空间定位、平均/波动分解、选一个有证据的处理 |
| EM几何功率 | 中部截断已有初步结果，但对子采样数敏感，完整域/全柱未完成 | 统一参考、独立收敛求积、域覆盖与实际生产Q验证 |

不将“背景仍存在”重新解释为浮力公式错误；不将“同采样场很准”解释为连续几何总功率已经通过。

### 0.3 本轮禁止事项

- 不盲目继续延长预处理，不新增多个预处理时长扫描。
- 不同时修改beta、mu、声速、压力、连续性、自由面与黏性。
- 不将总速度扣除平均值后作为正式流场。
- 不通过投影、滤波或强阻尼制造漂亮环流再当作生产结果。
- 不通过电流/功率/体积整体归一掩盖标准例误差。
- 不将离线求积权重直接改写生产Mass/Vol。
- 不改A-phi/Laplace独立研究分支。
- 不换流体框架、不加新的工程壁面ghost布置。
- 不重做TEAM7全套、弹性桨、温变物性和导电底板支线。
- 不把客户几何固定温场测试称为法国自然对流复现。

## 1. A｜冻结版本、配置与证据身份

### A1 实际源码冻结

1. 读取236补丁与当前工作树，明确WIP进入了哪些运行。
2. 保存用户文件，不清理、不丢弃改动。
3. 将本轮使用的源码形成可追溯提交；如确有不能提交项，记录HEAD+完整补丁SHA256。
4. 记录编译器、Real类型、浮点选项、设备、二进制SHA256、Reload、有效参数。
5. 功能提交拆分：证据修正、流场后处理、必要物理修复、EM求积及报告。

### A2 正式对比显式固定初始化

236的auto仍兼容旧行为，会按zero_q等热标签选择不同初始化。

本轮正式流体比较必须显式使用同一已测机械策略，例如已有的：
`--mechanics-initialization-mode=skip_dens_fp`。

同时固定：
- hydro_settle_time=1.0 s；
- 匹配外层实际推进路径；
- 同一Reload和有效物性；
- 同一热源/浮力开关作用域。

不要只把默认预处理改成1 s，而让auto再次引入不同机械清理。软件机示例配置也应在接收时显式指定策略，不要求修改所有历史默认。

初始化不同入口的t0字段比较已做过，本轮无需完整重跑；如修改初始化实现，再补相关对照。

### A3 修正矛盾JSON

236的 `D_fixed_T/D2_DIFF_SUMMARY.json` 与正式报告冲突：
- direction_correct=false；
- buoyancy_dominated=true；
而正式结论为方向正确、尚未由浮力主导。

追溯该文件使用的源数据、脚本和比较空间。若来自移动后的粒子ID直接配对：
- 标注为Lagrangian诊断，禁止作为Eulerian空间流场验收；
- 放入明确的历史/无效证据目录并附原因；
- 不删除历史原始数据；
- 正式结果由唯一的共享网格后处理生成。

顶层报告、验收表、JSON、图和脚本必须一致。不能仅手工改JSON布尔值。

## 2. B｜优先后处理已有场，定位局部背景

### B1 数据来源与时窗

优先使用236预处理1 s后的OFF/ON原始场，保留0.1、0.5、1 s输出以及动量、能量日志。

先确认是否有足够时间采样：
- 两帧只能作瞬时比较，不能据此声称稳定时间均值或随机波动；
- 若缺足够场，补一次OFF保持到1 s的运行；
- 后段0.5–1 s建议保存至少6个共同时间场；
- 轻量U、动能、质心、分区量按真实推进dt累计，使用时间权重。

不因做平均而默认延长到10/20 s。稀疏场不能用于判断声学高频谱；本轮不开展频谱研究。

### B2 原始流场必须保留的分解

在固定Eulerian几何网格上做：

1. 质量加权全域平均速度与平均平移动能；
2. 相对平均平移的剩余运动，仅作诊断；
3. 轴对称平均u_r、u_z、u_theta；
4. 三维场偏离轴对称平均的部分；
5. 后段时间平均流与时间变化量。

必须说明每个平均的顺序、权重与有效掩码。不要把周向平均后的弱流速称作原始三维背景已经降低。

平均流也可能是非物理持续回流；波动也不天然等于无影响噪声。需要与OFF条件、力学预算和浮力信号共同解释。

### B3 几何分区

建议用与dp相联系的固定带：
- 侧壁带；
- 底部带；
- 自由面附近；
- 顶/底角；
- 深内部。

使用互斥分区，或明确重叠统计不能相加。导出每区粒子数/代表体积、质量、动能占比、速度RMS、平均u_z、u_r和Umax。

除总动能占比外，也报告该区占总质量/体积的比例，避免小区域高速度被平均掩盖。

同时给近壁距离和高度上的剖面，找出是某一层局部运动，还是穿过内部的有组织回流。

### B4 与实际力学贡献对应

复用236声学第一半步的实际读入力诊断。

注意：
- Force/ForcePrior是聚合通道，不能直接称为已拆清的内部压力/壁面压力；
- 粒子位于底部，不代表它的Force只来自底壁；
- 总动量闭合不代表局部离散准确。

先在已定位的高背景区记录压力/黏性/重力及可分离的壁面作用。只对该区域添加必要诊断，避免重建一整套新审计系统。

检查背景是否在以下事件后系统变化：
- 密度重整；
- 自由面标记/大气压处理；
- 位置与邻域刷新；
- 压力梯度或核修正；
- 黏性力刷新。

本轮必须给出一张“区域—事件—有效力—速度响应”的对应表。不能只把全部剩余问题叫作自由面噪声。

### B5 Python图与数据

至少交付：
- OFF原始速度模、矢量图及分区能量；
- OFF后段平均矢量与时间变化量；
- ON/OFF同网格、同色标、同箭头缩放对照；
- ON−OFF矢量图与u_z(r)；
- 支持量、覆盖率和真实输出时刻。

差值图可独立放大，但必须标注放大比例。流线是辅助，原始矢量与数值不可省略。

容器外与缺支持区域不补零。网格范围来自固定几何，ON/OFF共用。相同粒子ID不等于相同空间位置。

## 3. C｜只做一个有证据的局部修复或诊断对照

### C1 根据B结果选择，不预设答案

| B的结果 | 本轮允许的下一步 |
|---|---|
| 背景主要是初始化后仍衰减的短暂波动 | 选择明确的观测起点，记录初态与衰减；不无限加预处理 |
| 某个近壁/自由面事件持续生成速度 | 定位对应调用或状态一致性，做最小范围修复 |
| 深内部也持续存在有组织OFF回流 | 检查局部离散力平衡，不能仅靠平均滤波处理 |
| 原始三维背景较强但主要为非轴对称变化 | 量化其对轴对称平均与温度输运的影响，不直接宣布可忽略 |
| 后处理插值/配对制造了较大差值 | 只修后处理，保持生产求解器不变 |

最多选择一个主要候选。其他候选列依据和优先级，不全部并行扫描。

若证据指向某个已有开关，可做单因素诊断，但必须注明改变了哪项数值/边界处理，不能把模型变化当作同一物理工况的直接修复。

### C2 修复后最小运行

同几何、物性、初始化、预处理和调度：
- OFF先0.2 s；
- 无异常且背景指标改善，再到1 s；
- 若需对浮力响应判断，ON配对到同1 s。

现有固定温场保持：
`T=1273.15+20*((r/R)^2-0.5) K`，
每次生产浮力读取前按当前位置施加。OFF仅关闭Boussinesq，保留重力。

如果没有物理修复，仅补齐时间采样，不机械再跑第二组ON/OFF。

### C3 不能只用“RMS下降多少”验收

分别回答：
1. 背景空间来源是否查明？
2. 总三维背景、轴对称背景、时间平均回流、波动各是多少？
3. ON−OFF信号对共同网格与合理重分箱是否稳健？
4. 背景是否显著改变规定温场下的环流判读？
5. 修复是否保持质量、动量与初始化一致性？

可设置“OFF RMS下降一个数量级”为工程改善目标，但不能达不到就不断调参，也不能达到就直接称法国复现通过。

若信号比后处理差异或重复运行变化还小，保持RESPONSE_SIGNED，不授RESPONSE_RESOLVED。误差底只在必要时增加一次有针对性的重分箱/重复测试，不展开新的矩阵。

“总流场由浮力主导”必须有空间和时间证据，不凭ON/OFF RMS接近或不接近单独判定。

## 4. D｜EM参考统一：不再混用弱极限与完整趋肤解

### D1 固定案例

R=0.10 m，总高0.60 m，B0=0.01 T，sigma=16 S/m，f=300 kHz；
dp=0.025和0.0125 m，使用已达到相近严格迭代残差的配置。

本轮不加第三档、不增加柱长、不继续主要调Picard容差。现有同坐标复场精度成果保留。

### D2 主参考必须一致

正式功率对照使用与求解场比较一致的完整趋肤解析Q，计算：
`P_ref = H_region * integral_0^R(2*pi*r*Q_full_skin(r) dr)`。

中部参考区域明确为实际配置的上下界；236使用z∈[0.255,0.345]、厚0.09 m，需核对源码/effective。

弱反馈极限只作辅助校验，不与full_skin混用。约0.4%的参考差对百分之一级精度判断不可忽略。

主相对误差统一：
`e_ref=(P_num-P_ref)/P_ref`；
同时给绝对W和有符号差。旧max分母指标另名，不用来表示前后改善。

### D3 独立参考积分验证

用独立脚本、double精度积分完整解析Q：
- 比较至少两档积分精度；
- 相邻精度变化优先≤0.1%；
- 核查peak/RMS和1/2约定；
- 用弱极限作趋向或量级校验，不替换本工况主参考。

独立脚本不要调用生产求积后又拿相同结果做参考。

## 5. E｜EM截断权重：先消除采样假象，再谈功率精度

### E1 236的结果为何尚不能PASS

细档截断功率误差（当时对弱极限参考）：
- nsub=6：约−0.07%；
- nsub=8：约−0.98%；
- nsub=10：约−1.97%。

该变化说明求积尚未稳定；不得挑nsub=8后标“约1%”。

此外，截断单元之和只是现有粒子单元并集，可能不完整覆盖圆柱；求和接近几何体积也不能排除缺口与重叠相抵。

### E2 轴向交叠必须改为精确计算

规则立方单元的轴向交叠：
`Lz=max(0,min(z_i+dp/2,z_hi)-max(z_i-dp/2,z_lo))`。

不要再在z方向用6/8/10个点做粗布尔采样。0.09 m切片与粒子层不整除时，该采样会造成明显跳变。

先确认代表单元的方向、大小和坐标来自真实生成方式；非规则粒子不能无依据套用dp立方体。

### E3 横向圆—矩形交面积独立求积

可用解析方法、稳定一维积分或可细化二维采样。若采用采样，建议32/128/256等明确递增精度，利用相同xy列缓存，避免对每个z层重复算。

每档输出：
- 单元交面积与交体积；
- 总体积；
- 数值/解析中心Q乘权重的功率；
- 相邻精度变化；
- 覆盖缺口与重叠检查。

求积变化应小于要讨论的物理误差。百分之一级功率结论，优先要求求积细化变化≤0.1%；若未达，报告不确定度，不授精度PASS。

### E4 审核复算的对照值

用236所依赖的原中部粒子CSV，轴向精确交叠、横向细化，并仍以中心Q乘交体积，离线复核大致得到：

| 档位 | 现有单元并集体积误差 | 功率相对full_skin参考 |
|---|---:|---:|
| dp=0.025 | 约−2.55% | 约−4.59% |
| dp=0.0125 | 约−0.90% | 约−1.55% |

这里只是审核对照，不是硬编码目标，也不是完整圆柱功率的最终验证。真实值可能受源CSV截断精度影响，应按本机完整数据重算。

如果不能接近上述量级，优先检查参考、轴心、切片界限、实际xy坐标和权重；禁止悄悄改区域匹配答案。

### E5 中心取样误差与覆盖缺口分开

至少用三种场检查：
- 常数1：体积；
- r²：径向矩；
- 完整解析Q：实际工况积分。

分别比较：
1. 粒子中心值乘交体积；
2. 对单元内部场函数真正积分；
3. 连续完整几何积分。

这能区分中心取样误差、单元覆盖缺口与求积误差。数值场未知时不能直接用解析函数代替；解析单元积分只作为参考拆分。

不能整体归一权重，也不能因中心Q误差存在而宣布EM场错了同样多。

## 6. F｜补齐完整域和全柱数据，完成实际生产Q检查

### F1 不再以只有中部CSV为由延后

先检查现有完整输出。若缺少，授权仅重跑原两档解析柱，导出：
- 所有粒子Position、稳定ID、V；
- 复J各向量分量；
- 生产实际Q；
- sigma、频率、B0、轴心、上下界、源和收敛身份。

这属于必要数据导出，不是新增参数扫描。若字段读取接口缺失，做最小导出修改。

特别注意：此前J_sph可能是J_theta幅值，`J_sph^2/(2*sigma)`不自动等于生产实际Q。必须并列导出与核查，不得只靠重构功率替代生产焦耳热验证。

### F2 三条功率链分开

1. 生产实际Q的积分；
2. 完整数值复J向量重构的 `|J|²/(2*sigma)`；
3. 仅J_theta重构的功率（如需保留旧指标）。

说明它们的离散定义及差异；若edge-flux JouleHeat与重构J的定义不同，报告兼容性和影响，不为对上而强行设置相同。

中部同坐标场精度好，不自动证明全柱、边缘或生产Q都同样准确。

### F3 完整域覆盖方案

现有粒子单元不完全覆盖时，实施一个覆盖真实圆柱的独立积分点集：
- 求积点覆盖完整横截面及实际高度；
- 使用数值场局部重构得到复J或实际Q；
- 明确边界支持、外推比例及条件性；
- 不新增物理求解粒子，不赋予新Mass；
- 不用解析J/Q补填缺失的数值值。

可以优先采用已有可靠插值/局部最小二乘，但须用已知常数、r²及独立光滑场验证其误差。轴心方向与复数相位正确处理。不要让拟合天然套用完整解析解形式而失去独立性。

如果无法可靠重构，允许该项OPEN，但必须实际导出完整场、完成支持/覆盖诊断，并说明哪个区域或接口阻塞。不得再只重复“需要全柱数据”。

### F4 完整域求积验收

对全柱和中部、两档分辨率分别报告：
- 原生产权重功率；
- 相同采样解析功率；
- 截断权重功率；
- 完整域独立求积功率；
- 连续full_skin解析参考；
- 体积、r²矩、覆盖缺口/重叠；
- 重构与求积不确定度；
- 实际生产Q与两种J重构功率的差。

独立求积的常数/r²矩内部目标≤1%；若要声称约1%的物理误差，需要更小且明确的数值不确定度，不能只满足1%体积门就直接通过。

全柱对无限长解析解可能仍有模型差异；保持原几何，不通过加长柱体来回避。中部与全柱分别给结论。

### F5 生产积分影响说明

本轮默认不修改Mass/Vol。必须给出一张影响表：
- 原V用于哪些EM源积分；
- 用于哪些焦耳热积分；
- 用于热源传递/归一的哪个环节；
- 本轮修正仅改变报告，还是实际生产积分；
- 若生产权重仍偏，剩余误差预计影响什么，证据是什么。

离线功率改善不能直接标“生产积分已修”。TEAM7、局部峰值Q继续独立OPEN，不扩大解析柱结论。

## 7. G｜明确本轮完成后如何推进，避免长期停在诊断

### G1 三种可执行结论

| 本轮结果 | 下一步 |
|---|---|
| 局部背景已定位/明显改善，固定T响应可信，EM功率证据完整 | 下一轮恢复一个规定热温度—流动工况，再门控接入真实EM |
| EM已完成，局部背景仍主导但具体来源清楚 | 只修一个明确的局部力学缺陷，不重开浮力公式/EM扫描 |
| 局部背景影响已可量化，EM完整域仍受支持不足阻塞 | 可独立准备规定热主线，但真实EM热源定量仍保持限制 |

本轮不重复规定热1 s smoke。是否恢复热演化以B/C/F结果作决定，不机械等待所有报告标签变成PASS。

### G2 不能无限增加验证

如果局部背景暂时无法完全压低，报告其对目标观测的可量化影响与可用范围。不能既没有定位又用更多长跑拖延，也不能因未达到“背景绝对零”无限阻塞所有独立工作。

客户自然对流温度定量、法国文献环流拓扑、软件审查分别列状态，不将其中一个目标代替另外两个。

## 8. 执行顺序和预算

1. A：冻结版本，清理矛盾证据身份。
2. B：复用236场定位局部背景；不足时补一条OFF时间采样。
3. D/E：独立参考与稳定截断求积，可与流体后处理分目录执行。
4. F：必要时导出原两档全场，完成实际Q与完整域验证。
5. C：依据B只选一个主要修复；需要才做一对1 s固定T。
6. 汇总结果并选择G1路径。

默认预算：
- 流体已有数据优先复用；
- 一条OFF补采样到1 s；
- 一个主要候选，先0.2 s，必要时1 s；
- ON/OFF配对最多一组1 s；
- EM仅原两档必要导出，不增加物理参数；
- 无新增长时热耦合、预处理扫描及其他标准例全套。

停止条件：
- 没有空间/事件证据，不继续增加修复组合；
- 求积未收敛，不选择最有利采样数；
- 完整域缺支持，不用解析值补洞；
- 后处理权重变化，不自动修改生产粒子质量；
- 数据身份矛盾，不授正式验收标签。

## 9. 目录与交付标准

产物根：
`docs/ophelie/reference/doc237_runs/`。

建议目录：
- A_identity；
- B_off_localization；
- C_local_fix；
- D_em_reference；
- E_clip_convergence；
- F_full_domain_power；
- G_handoff。

必须交付：
1. ROUND237_REVIEW.md：具体发现、证据、修改、误差与限制；
2. ACCEPTANCE_TABLE.csv：作用域、状态、指标、参考口径、证据；
3. RUN_COMMANDS.md：真实执行命令、退出码、目录、参数与版本；
4. 源码补丁、Python脚本、二进制/输入身份；
5. 局部背景分区能量、平均/波动、原始矢量图及共同网格数据；
6. 唯一有效的固定T摘要，历史错误JSON附失效说明；
7. 完整趋肤参考积分及收敛表；
8. 轴向精确交叠、横向求积细化、覆盖/重叠与功率表；
9. 全柱/中部的实际Q、复J及对应积分，数据缺失时明确缺哪项；
10. SOFTWARE_HANDOFF_DELTA.md：可以择收什么，哪些物理声明不能升级。

大场可外置，但关键稀疏场、CSV、脚本、日志片段和补丁直接入包；不能仅给无法访问的本机路径。

正式摘要由实际脚本生成，不手工改布尔值凑报告。顶层、分项笔记、JSON、图和CLI输出的状态必须一致。

## 10. 本轮收口必须回答的六个问题

1. 剩余局部背景集中在哪里，是持续平均回流还是时变运动？证据足够吗？
2. 对应哪个实际事件/力学通道？本轮修改是否降低了它？
3. 固定温场信号在修正后处理与背景评估下，究竟验证到哪一级？
4. 电磁截断求积是否收敛，现有粒子单元还缺多少几何覆盖？
5. 全柱/中部的实际生产Q功率误差分别是多少，是否与J重构功率一致？
6. 下一轮具体进入规定热推进，还是只修哪个明确阻塞？不能再只写“继续优化”。

本轮成功标准：不重复已完成的浮力公式验证，用有限工作量完成局部背景归因和可信几何功率证据，使下一轮能够作出具体推进决策。


<!-- ===== END 237_OPHELIE_LOCAL_BACKGROUND_AND_EM_QUADRATURE_CLOSURE_CURSOR_PLAN_2026-11-23.md ===== -->


<!-- ===== BEGIN 238_OPHELIE_OFF_DECAY_DISCRETE_WORK_AND_VALIDATION_CLOSEOUT_CURSOR_PLAN_2026-11-23.md ===== -->

# 238｜OPHELIE 无浮力回流衰减、离散做功归因与电磁验证归档执行计划

- 日期：2026-11-23。
- 执行对象：物理开发机 Cursor。
- 基线：第237轮交接；报告身份为 `845e580de` 加物理WIP，开工必须核实实际源码与二进制。
- 本轮主线：一次3–5 s的无浮力保持测试，判断剩余回流是初态残余衰减还是由数值环节维持；据证据只做一个具体修复。
- 本轮电磁工作：整理现有解析柱验证，不再反复扫描、继续追硬1%或擅自修改生产积分。
- 本文为授权执行的开发与测试范围，不要求每个常规步骤再征求确认。

## 0. 必须理解的现状与本轮边界

### 0.1 不再重开的内容

1. 浮力公式、专用通道和实际声学第一半步读取已在236得到证据。本轮不重做全套浮力单元验证。
2. 1 s预处理显著改善了整体下沉/回弹。本轮固定该已测设置，不扫多个预处理时长。
3. 热边界local_dt_dz有限制造场结果保留，不再扩展一般二维局部边界算法。
4. 解析柱电磁场随加密改善，实际Q与复J重构在本例一致。不要因为流体背景未解决而重新怀疑整个电磁路线。

### 0.2 当前必须回答的问题

237发现OFF已有“芯部上升、近壁下降”：
- OFF RMS约3.41e-3 m/s；
- 侧壁带约20%质量，约45%动能；
- t=1 s的轴对称分量约占动能53%；
- 固定T的ON−OFF矢量差约5.5e-5 m/s。

这些说明背景会干扰自然对流拓扑判读，但并未证明侧壁就是供能源。

237所谓后段平均实际仅有0.5、1 s两帧。两帧不能可靠判定稳定平均、衰减包络或随机波动。

本轮必须区分：
- 尚未衰减的初始动能/压缩储能/重力势能；
- 压力或其他离散环节维持的回流；
- 特定事件产生的直接速度/位置修正；
- 后处理造成的误判。

### 0.3 禁止错误归因

- 净动量/净冲量接近零，不等于没有做功或供能。
- 动能集中在侧壁，不等于侧壁算子就是根因。
- 压力做正功，不自动等于非物理供能；可能是压缩储能向动能的可逆交换。
- 总动能不增长，不等于背景对温度输运可忽略。
- 减去OFF流场或滤掉平均回流，仅是诊断，不是修复后的实际解。
- 不能以程序passed、未穿墙或图形闭合代替物理归因。

## 1. A｜冻结版本与明确机械配置

### A1 版本

核实237代码、WIP、二进制、Reload和effective配置。保存用户文件，不清理或丢弃改动。

本轮实际运行源码应形成可追溯提交。确不能提交的项，记录HEAD+完整补丁SHA256。按诊断、修复、后处理、报告拆分提交。

不自动合并软件分支，不改独立A-phi/Laplace分支。当前分支沿用用户安排。

### A2 固定配置

以237已测客户NC工况为主：
- dp=0.010 m；
- H=0.180 m；
- R按实际effective，约0.1415 m；
- rho=2700 kg/m³，mu约14.532 Pa·s；
- beta=0，仅关闭Boussinesq增量，真实重力保留；
- 均匀固定温度，关闭热源、热扩散及温度换热推进；
- hydro_settle_time=1.0 s；
- 显式使用已有 `--mechanics-initialization-mode=skip_dens_fp`；
- 匹配外层实际步长路径；
- 壁面、自由面、核修正、声速等保持237基线。

正式比较不使用auto初始化策略。热源开关不能再次改变机械初始化。

### A3 ID与数据身份修正

检查 `particle_i` 是稳定ID还是排序后下标。此前“1 s配对位移0.27 m”与毫米每秒级流速不相称，必须核对是否存在排序错配。

要求：
- 跟踪使用真实稳定ID（例如实际实现的OriginalID）；
- 逐帧检查ID唯一性、集合和排序映射；
- 若没有稳定ID，不做逐粒子轨迹结论；
- Eulerian场仍用共同网格，不以稳定ID替代空间插值。

更正历史无效JSON的说明，区分“真实位移导致空间不对应”和“数组下标排序错配”。保留原文件与失效记录，不手工改布尔值掩盖问题。

## 2. B｜实现与生产推进一致的离散做功测量

### B1 先列出实际更新序列

复用236的动量预算代码，明确：
- 声学第一半步；
- 声学第二半步；
- 黏性力计算/保持；
- 重力与其他ForcePrior；
- 位置更新；
- 密度重整、自由面处理；
- transport/shifting（如实际启用）；
- 排序与邻域刷新；
- 任何直接改速度的限制或清理。

不能只测第一半步，遗漏第二半步或外层直接写速度的操作。

逐项列出字段单位、时间层与实际更新式。若ForcePrior混有黏性/重力，读取实际已有分项通道；没有时做最小只读诊断。

### B2 动能基准

对固定质量粒子：
`K=sum(0.5*m_i*|v_i|²)`。

分开报告：
- 总K；
- 平均平移动能；
- 轴对称平均分量动能；
- 三维偏离分量；
- 几何分区动能。

轴对称分解、空间粗粒化和时间平均不可混用。报告网格支持及分箱变化影响，不把高分辨率稀疏bin的拟合程度当作真实轴对称精度。

### B3 每个实际速度更新子步的功

对于实际冲量更新：
`Delta(mv)=sum_k I_k`，

用同一子步的平均速度：
`v_mid=(v_before+v_after)/2`，

计算每通道离散功：
`W_k=sum_i(v_mid_i dot I_k_i)`。

核对：
`DeltaK_observed = sum_k W_k + R_K`。

分项冲量必须来自实际更新使用的力、dt和子步系数，不能把事后DeltaK反推成某个分项功后自证闭合。

如更新不是单纯显式kick，按真实离散式处理：
- 隐式黏性可用该算子的实际速度改变量和对应动能变化，明确不能再与同一黏性力重复计数；
- 任何直接速度投影/限幅/清理单独记为事件功；
- 一个算子同时影响多个量时明确归属，避免重复累计。

阶段中存在粒子排序时，前后数据必须用稳定ID对齐，不能按数组下标相减。

### B4 最少分项

- 重力功；
- 黏性相关功/动能改变；
- 内部压力相关功；
- 壁面压力相关功（可可靠分离时）；
- 其他实际力；
- 直接速度修改的事件功；
- 未解释残差。

若无法可靠拆分内部/壁面压力，先报告聚合压力功，保留该归因限制；不能根据粒子所在区域把所有力命名为壁面力。

固定壁面的无滑移黏性通常应有耗散作用，但以实际符号约定和离散式核对。若有局部正功，先检查分项/时间层，不立即修改黏性公式。

### B5 重力势能、压缩储能与位置事件

记录：
`E_g=-sum(m_i * gravity_vector dot x_i)`，
使用固定势能参考，并检查重力功与DeltaE_g的关系。

若位置在外层单独更新，速度kick的重力功和位置势能变化未必逐步一致；按完整外层核对，残差单独列出。不能为了闭合而按不相容的位移重算力。

若当前压力模型是已知barotropic EOS，可在确认实际EOS与压力路径后计算相应压缩储能。可参考一般关系：
`e(rho)-e(rho_ref)=integral[p(s)/s² ds]`，
但必须核对实际EOS、压力偏置和背景静水处理，不能直接套给扰动压力路径。

若本轮无法可靠构造压缩储能，则：
- 保留压力功和密度/压力变化；
- 不把压力正功直接判为非物理能量注入；
- 用时间演化和具体事件相关性缩小范围。

热焓账本不参与本项闭合；本测试不推进温度。不要将声学储能误记为热输入。

### B6 分区功与跨区输运

全局动能账本是主依据。侧壁等Eulerian分区还存在粒子跨界和动能输运：
- 分区功可以指出作用集中的区域；
- 分区DeltaK不等于本区做功，除非同时计入输运；
- 若不实现输运预算，明确标为定位统计，不授局部能量闭合PASS。

需要逐粒子归因时可用稳定ID的固定粒子集合辅助，但与固定几何区域分开命名。

### B7 仪器短测

先0.02–0.05 s，只验证做功记录、单位、符号、子步与累计。
- 记录实际dt，不能用候选dt代替；
- 核对诊断开/关是否改变物理轨迹；
- 优先double规约，生产Float可保留；
- 不用1 J这种远大于本例动能的固定归一分母。

建议归一尺度：
`max(K_start,K_end,sum(abs(W_k)),E_floor)`，
E_floor只处理接近零，记录其具体J值和选择理由；同时始终给绝对残差J。

内部目标：可解释的累计动能规约残差≤1%上述尺度。未达时先修测量，不直接动流体模块。该门只检查预算，不证明物理正确。

## 3. C｜一条有限的OFF衰减测试

### C1 运行安排

使用A2固定配置，先预处理1 s，再进入正式保持：
- 保持先到3 s；
- 若仍需区分衰减与平台，且未出现明确错误，允许同一路径继续到5 s；
- 总保持上限5 s，本轮不再延长到10/20/100 s。

优先同一连续运行在3 s检查。若框架不支持中途续跑，安排一次5 s运行并在3 s节点输出；不要把只有Position的Reload当完整状态续算。

禁止每次分段重新归零速度或清理压力/密度来“接着跑”。

### C2 输出频率

轻量监测：
- K、分项累计功、E_g、质心、U_rms、Umax；
- 分区动能与速度；
- 时间累计、实际调度、密度/压力范围。

按实际子步/外层累计；物理监测输出间隔建议0.01–0.02 s。不要因输出间隔较粗漏累计做功。

全粒子场稀疏输出：
0、0.1、0.5、1、1.5、2、2.5、3 s；
若到5 s，再输出3.5、4、4.5、5 s。

同一时间窗至少有足够的共同网格场再做时间平均。对轻量监测按实际时间加权积分，不把非均匀时间点简单平均。

### C3 分析窗口

建议0.5–1、1–2、2–3、3–5 s分段。
报告：
- 窗口时间加权K与速度；
- 轴对称回流强度；
- 各分项功、平均功率；
- 黏性耗散与压力/重力交换；
- 液面、密度、压缩储能（如可靠）变化。

比较窗口与误差，不用两帧相近直接定义稳态。不要假设K必须逐步单调：声学与势能交换可能使其振荡。

### C4 结论分支

| 观测 | 判定及行动 |
|---|---|
| K及轴对称回流包络持续下降，无未解释的维持功 | 初态残余衰减为主要解释；选择有依据的初态/观测起点 |
| K维持平台，正功与耗散平衡且可定位到具体事件/通道 | 有维持机制；进入D的一个最小修复 |
| K振荡，压力/势能储能交换可以解释 | 不称持续数值供能；检查衰减包络和观测时窗 |
| 动能预算未闭合 | 先修测量，物理根因保持OPEN |
| 仅侧壁动能高，但功/事件未定位 | 只能称空间定位，不能称边界根因已查明 |

L²/nu约数秒仅作为本轮测试时长的量级依据，不是已知的精确指数衰减常数。禁止把观察曲线强行拟合到预设答案。

## 4. D｜根据做功证据只选一个修复

### D1 允许的最小处理

- 如果初态残余衰减：保留完整的已松弛力学状态，避免重新破坏平衡；观测起点由真实残差与趋势选定。
- 如果密度/压力重置造成重复激发：修正该状态/调用顺序。
- 如果特定壁面/压力贡献维持背景：核对实际接触状态、法向、邻域、系数或时间层，只修有证据的缺陷。
- 如果实际黏性作用异常：先核查力/加速度单位、刷新和分步，再决定是否修实现。
- 如果直接速度/位置修正导致动能改变：显式记录并检查该事件，不能静默视为物理力。

“流在侧壁较强”不足以授权任意修改壁面算子。

本轮最多一个主要修复方案；无证据不展开多个开关组合。可保留一个有明确目的的诊断对照，注明改变的模型/离散，不能冒充同物理修复。

### D2 禁止的处理

- 正式段清零侧壁或全域速度；
- 新加强阻尼压住流场；
- 把OFF回流直接从ON解中扣除后称真实结果；
- 调大beta、降低mu使浮力压过错误背景；
- 同时更换压力、连续性、黏性和自由面；
- 为通过测试整体归一力、质量或体积；
- 修改独立A-phi分支。

### D3 修复验证预算

修复后先0.2 s查启动与预算，再到1 s查局部回流。
若本轮主要结论是后段维持机制，允许只将该修复候选延至同3 s窗口一次；不得重新跑5 s多参数矩阵。

对比相同物理时窗，不能拿基线1 s与修复3 s直接宣称改善。

记录：
- 原始三维背景；
- 轴对称平均回流；
- 动能/分项功；
- 质量、密度、液面与实际调度；
- 是否引入新的守恒或边界问题。

## 5. E｜固定温场与恢复热演化的进入门

### E1 不重做浮力公式证明

236的公式、通道及动量读取证据继续有效。若本轮修改实际更新路径，做针对性回归；否则不再重复静态力链快照全套。

背景已有明确解释/改善后，最多做一对1 s固定温场：
`T=1273.15+20*((r/R)^2-0.5) K`。
同一完整初态、机械策略和时钟起点，OFF仅关闭浮力。

Eulerian温场按当前坐标施加。ON/OFF使用同一共同网格、权重、有效支持、真实输出时刻与稳定ID记录。

### E2 状态分级

- RESPONSE_SIGNED：方向合理；
- RESPONSE_RESOLVED：差值超出后处理/运行误差底且空间时间结构可信；
- FLOW_BUOYANCY_DOMINATED：总主要环流不能由OFF背景解释；
- FRENCH_TOPOLOGY_MATCH：独立法国配置验证，不在本轮授予。

不能只凭总RMS差小就判没有浮力，也不能只凭方向正确就判整体自然对流通过。

### E3 恢复规定热的条件

至少满足：
1. 背景的时间演化与能量来源已能解释；
2. 初始化策略显式固定，热源标签不改变力学启动；
3. 诊断/后处理没有未解决的身份错误；
4. 浮力响应可辨，背景对目标观测的限制已量化。

本轮可以完成后续规定热运行配置和输出清单，但不机械再跑重复1 s smoke。实际恢复温度演化应依据上述结论在下一轮执行。

若背景仍维持且根因未查清，保持阻塞，不带着未解释回流开始长时热输运。

## 6. F｜电磁本轮以验证归档为主，不继续无休止追精度

### F1 复用237成果

| 指标 | dp=0.025 m | dp=0.0125 m |
|---|---:|---:|
| 原生产权重全柱功率误差 | +6.05% | +6.89% |
| 截断权重全柱功率误差 | −4.55% | −1.50% |
| 截断权重中部功率误差 | −4.59% | −1.55% |
| 截断体积误差 | −2.55% | −0.90% |

参考为相同完整趋肤解析解，e_ref=(P_num−P_ref)/P_ref。
生产Q与复J重构在本例代数一致；这不是独立实验验证。

以上数值使用已有完整CSV和求积收敛结果，无须再跑求解器。

### F2 编写解析柱验证报告

报告至少包含：
- 物理模型、源与边界设定；
- 参数、单位、峰值/RMS约定；
- 解析参考、同点复场误差和分辨率趋势；
- 原生产权重、截断权重与几何覆盖的区别；
- 全柱/中部功率及实际Q三条链；
- 求积收敛证据；
- 适用范围与剩余验证项。

明确写：
- 当前生产权重细档仍约+6.89%；
- 约1.5%是离线截断权重结果；
- 现有单元并集仍有几何覆盖缺口；
- 没有修改生产Mass/Vol；
- 不能宣称电磁所有case或峰值Q已验证。

不要追硬1%而新增第三档、继续调Picard或增加柱长。当前解析柱可以形成一项有限范围的验证证据，不要求为写报告消灭所有小误差。

### F3 遗留项登记，不在本轮扩展

- TEAM7；
- 强趋肤范围；
- 局部峰值Q；
- 生产几何权重系统偏差；
- 温变sigma与真实线圈工程热源。

如原生产误差仍在项目内部容许范围内，按真实数字列入报告；不能把单例指标自动等同于整个合同电磁模块验收。

本轮只保留明确的后续需求，不重跑这些支线。

## 7. G｜可执行的运行与开发顺序

1. A：冻结WIP身份，核查稳定ID和有效配置。
2. B：做功诊断实现，0.02–0.05 s仪器短测。
3. C：一次OFF保持到3 s，必要时同路径到5 s。
4. 分析衰减/储能交换/维持机制，必须实际作出分支判断。
5. D：只选一个有证据的修复，有限短回归。
6. E：条件满足才做固定T一对1 s，准备恢复规定热的配置。
7. F：并行整理电磁现有验证报告，不新增EM扫描。
8. 提交本轮结论、补丁与下一轮唯一主要动作。

### G1 本轮默认预算

| 工作 | 上限/约束 |
|---|---|
| 做功仪器验证 | 一次0.02–0.05 s |
| OFF衰减 | 一条运行3–5 s |
| 修复候选 | 一个，先0.2再1 s；必要时同窗口到3 s |
| 固定T响应 | 条件满足才ON/OFF各1 s |
| 新EM求解 | 默认0次，复用237数据 |
| 热耦合长跑 | 0次 |
| 预处理时长扫描 | 0组 |

本轮不是重新做所有标准例。若修改公共内核，补与改动直接相关的既有回归，范围由实际修改决定。

### G2 停止条件

- 动能预算不闭合，停止物理根因定论，先修测量。
- 5 s仍无法区分趋势，不继续加时；报告未解释能量通道与所缺证据。
- 一个修复无效，不能自动扫更多组合；回到实际预算定位。
- 仅做了后处理平均，不声明背景已治理。
- 电磁只写报告，不升级TEAM7或峰值Q状态。
- 未达到进入门，不开展长时热源耦合。

## 8. 产物、目录与验收

产物根：
`docs/ophelie/reference/doc238_runs/`。

建议目录：
- A_identity；
- B_work_instrument；
- C_off_decay；
- D_single_fix；
- E_fixed_T_gate；
- F_em_validation_report；
- G_handoff。

必须交付：
1. ROUND238_REVIEW.md：衰减还是维持，明确结论与证据；
2. ACCEPTANCE_TABLE.csv，作用域/指标/参考/状态完整；
3. RUN_COMMANDS.md，区分已执行与准备命令；
4. 源码补丁、版本、有效配置、二进制与输入校验；
5. 各子步/外层动能预算CSV和独立累计脚本；
6. K、分项功、势能/压缩量、质心与分区趋势图；
7. 至少两个真正包含多帧的共同网格时间窗，不用两帧假装稳定均值；
8. 修复前后同窗口原始速度/矢量图；
9. 固定T未运行时明确NOT_RUN，不沿用旧PASS冒充本轮；
10. 电磁解析柱验证报告，生产与离线权重结果分开；
11. 软件择收说明，避免把诊断功能默认用于正式生产。

若计算压缩储能，附实际EOS推导与适用条件；未计算则明确限制。不要为了形式完整伪造储能闭合。

大场可外置，但关键稀疏场、CSV、脚本、日志片段和补丁直接入包。正式JSON从对应脚本生成；历史无效JSON保留失效说明。

### 8.1 建议验收项目

- STABLE_PARTICLE_ID_VERIFIED；
- DISCRETE_KINETIC_WORK_BUDGET；
- OFF_DECAY_OR_MAINTENANCE_CLASSIFIED；
- BACKGROUND_SOURCE_IDENTIFIED；
- SINGLE_FIX_EFFECT；
- FIXED_T_RESPONSE_SCOPE；
- THERMAL_RESTART_READY；
- EM_ANALYTIC_CYLINDER_REPORT_READY；
- SOFTWARE_SELECTIVE_DELTA_READY。

状态使用PASS/PARTIAL/FAIL/OPEN/NOT_RUN，不能让程序smoke_passed覆盖物理结论。

## 9. 本轮最终必须回答的五个问题

1. 现有无浮力回流在3–5 s尺度上究竟衰减、振荡交换，还是被持续维持？
2. 哪项实际做功或状态/事件解释了它？净冲量与局部做功是否已经区分？
3. 是否完成了一个具体修复，效果和剩余误差是多少？
4. 当前是否能恢复规定热的温度—流动主线？如果不能，只剩哪个明确阻塞？
5. 电磁哪些证据可以进入验证报告，哪些仍只是离线改进或独立遗留？

本轮收口标准：用一次有限、有解释力的试验完成背景流分支判断，并推进一个具体处理；不再以“侧壁附近有流动”“需要继续优化”作为最终结论。


<!-- ===== END 238_OPHELIE_OFF_DECAY_DISCRETE_WORK_AND_VALIDATION_CLOSEOUT_CURSOR_PLAN_2026-11-23.md ===== -->


<!-- ===== BEGIN 239_FRENCH_THESES_EXTRACTION_AND_CURSOR_PLAN.md ===== -->

# 239 法国博士论文物性提取成果及 Cursor 下一步实施记录

日期：2026-11-24。范围：OPHELIE 法国自然对流复现；与 A-phi/Laplace 独立。

## 0. 本文件的性质

本轮已实际提取两篇博士论文共8条物性曲线，交付CSV、校准记录、脚本与原图叠加图。不是只给提取任务。原文参数、图像读数、拟合推断必须分开。

已完成：Jacoutot 图III.8/III.9/III.10/III.12四条曲线，Sauvage 图2.1四条Num曲线；电导率A/B比较；图中电导率关系的拟合推断；几何、边界与参考结果表。

未完成：完整二维温度场与速度场重建、一维温度/速度参考剖面的数字化、任何新的求解器仿真。reference中收录的结果图不是已数字化场。不得把这些标为完成。

本文件连同数据包替代上一条聊天中“物性曲线全部待提取”的状态说明。Cursor不必重新从零提取这8条曲线。

## 1. 文件与来源

- J：Jacoutot，2006，原文件“法国博士论文(2).pdf”，190页。
- S：Sauvage，2009，原文件“Sauvage博士论文.pdf”，174页。
- metadata.json含两份原PDF的SHA256、校准、方法、有效范围与误差估计。
- data/jac_III8_k_trace.csv：导热系数，W/(m K)，正文56/PDF69。
- data/jac_III9_cp_trace.csv：比热，J/(kg K)，正文56/PDF69。
- data/jac_III10_mu_trace.csv：动力黏度，Pa s，正文56/PDF69。
- data/jac_III12_sigma_trace.csv：电导率，S/m，正文57/PDF70。
- data/sauvage_fig2_1_{k,cp,mu,sigma}_vector.csv：S图2.1，PDF34。
- data/property_checkpoints.csv：从提取曲线插值得到的检查点，不是原文数据表。
- data/sigma_A_B_comparison.csv：印刷公式与图曲线的比较。
- data/sigma_B_fit_WITH_EXTRAPOLATION.csv：由图曲线拟合并延伸的候选输入；外推行显式标注。
- sigma_B_fit_inference.json：拟合系数、区间与残差。
- qa/：原图叠加与比较图。
- scripts/extract_properties.py：完整物性提取与再生脚本；从包含原PDF的项目工作目录执行。

### 1.1 提取方法

J的图是栅格图：直接提取原PDF嵌入图片，按坐标轴标定，人工引导黑像素跟踪，用5列中值抑制实验符号交叠与像素噪声。CSV保留原图像素坐标。mu按对数纵轴换算。虚线缺口不伪装成新增原始观测。

S的图2.1是矢量图：直接读取Num路径折线顶点，排除图例线段，按轴刻度映射；CSV保留PDF坐标。矢量提取精度不等于原始物性测量精度。

J图III.12的图内图例写EXP 1，图注称NUM；因此本包只称“图中曲线”，不声称已确认其就是最终代码输入。

### 1.2 有效范围与限制

| J曲线 | 提取温区约K | 点数 | 数字化误差估计 |
|---|---:|---:|---|
| k | 307–1795 | 270 | 约0.10 W/(m K)，另有约12K温标误差 |
| cp | 311–1792 | 267 | 约30 J/(kg K)，温标约15K；陡变区误差更大 |
| mu | 503–1789 | 378 | 约0.10 dex量级，实验符号交叠区可更大 |
| sigma | 954–1568 | 187 | 约0.2 S/m并叠加约12K横轴误差；低值相对误差大 |

这些是提取者的保守误差估计，不是统计置信区间。不得用CSV小数位推断高精度。表中点数不代表独立物理信息数量。

sigma低温段贴近坐标轴，不能可靠恢复；高温段超过图中20 S/m上限，不能直接读出。没有把不可辨识低温值填成零，也没有把上边界当成电导率饱和。

mu图中低温平台读数约8e4 Pa s，而正文上限写1e5 Pa s。两者分别保留，不擅自把提取平台改成1e5。主配置采用正文上限1e5作为上限约束；图中mu曲线本身仍作为图像提取数据。是否将低温平台替换为正文值，应明确配置并记录，不能静默处理。现阶段不为这一差异增加一轮大规模扫描。

### 1.3 使用规则

raw trace不能未经处理直接当高精度材料表。先检查叠加图，再生成稀疏或平滑的solver表；不得覆盖raw。处理时不得改变可分辨转折位置、平台、mu跨数量级趋势。对mu在log10空间处理图像噪声是数据整理选择，不代表原作者一定使用对数插值。其他参数可分段线性。输出明确插值规则。

所有越界都要报告粒子数、温区与处理方式。禁止静默钳位、零值填充、无来源外推。只有数K端点缺口可作为明确记录的端点延伸；sigma缺失温区必须使用单独标记的拟合模型或保持未就绪。

## 2. 原文几何与边界参数

J正文51/PDF64表III.1：

| 项目 | 值 |
|---|---:|
| 坩埚内径/外径/高度 | 500/530/350 mm |
| 坩埚分瓣数/间隙角 | 60 / 2度 |
| 计算玻璃半径/高度 | 250/185 mm |
| 表中玻璃质量 | 102 kg |
| 线圈直径/匝数/总高度 | 570 mm / 7 / 230 mm |
| 频率 | 282 kHz |

570mm直径的具体定义、每匝高度、导体截面、底部偏置尚未从本表确认。保留现有合理近似并记录；禁止补造“文献精确尺寸”。正文实验描述还有不同高度表述，计算几何按185mm及相应图核对，不混入装料高度。

J第三章正文55：rho=2800 kg/m3，beta=3e-5 1/K，mu上限1e5 Pa s。mu>1e4 Pa s用于冷层分析，不能自动实施为冻结/删除粒子。

J正文62–67：

| 边界 | 数值 |
|---|---:|
| 顶面对流系数 | 20 W/(m2 K) |
| 顶面对流参考温度 | 303 K |
| 顶面辐射率 | 0.9 |
| 顶面辐射参考温度 | 423 K |
| 冷却水参考温度 | 303 K |
| 第三章后续主要结果侧壁h | 300 W/(m2 K) |
| 第三章后续主要结果底部h | 300 W/(m2 K) |

向外为正：qtop=20*(Ts-303)+0.9*sigmaSB*(Ts^4-423^4)；qside/qbottom=300*(Tw-303)。辐射全部用K。不同参考温度分开输入。检查粒子表面积权重，不重复计入显式冷壁与等效h的同一热阻。

期刊参数另保留：已核对的期刊表1rho=2750、beta=1e-5，底部h=35 W/(m2 K)，结果总吸收功率50kW。博士论文第三章55kW的边界与结果不能覆盖成“同一个法国case”。

## 3. 电导率A/B：用户已授权两方案比较

### A：印刷公式，明确以K代入

sigmaA(T)=10**(3.7921-3179.8/T_K)，单位S/m。

配置名：JAC_PRINTED_FORMULA_K。1473K校验值42.990877 S/m。这是忠实实现印刷表达式，不等于证明原计算使用此式。

### B：图中曲线及其单独标记的延伸模型

纯读图表：jac_III12_sigma_trace.csv。1473K插值约13.677 S/m。

为使工程温区能运行，附带候选拟合：

log10(sigmaB)=3.76180-3869.11/T_K。

拟合使用图中1100–1550K区间；最大log10残差约0.00664，RMS约0.00212。系数小数是拟合结果，不代表原始材料精度。图外部分属于模型外推。先使用JSON完整系数再统一显示有效数字。

有价值的交叉核对：该拟合外推至1660K得到约26.98 S/m，与J正文71和图III.41的场最大值约27 S/m相容。它支持“图中关系有内部一致性”，不证明该拟合是作者原式，也不能证明低温外推正确。

配置名：JAC_GRAPH_FIT_EXTRAPOLATED。必须与纯读图数据身份区分。不要再擅自加入摄氏温标版或更换系数凑图；若以后新增解释，应单独列为假设。

| T K | sigmaA S/m | 图曲线B S/m |
|---|---:|---:|
| 1000 | 4.095 | 0.819 |
| 1200 | 13.876 | 3.426 |
| 1400 | 33.175 | 9.947 |
| 1473 | 42.991 | 13.677 |
| 1500 | 47.015 | 15.276 |
| 1550 | 55.032 | 18.506 |

## 4. 结果参考与适用范围

J第三章55kW、轴对称、包含Marangoni：正文70–73、图III.36–43、表III.5。

| 指标 | 原文值 |
|---|---:|
| 吸收功率 | 55 kW |
| Tmax | 1660 K |
| Umax量级 | 1e-3 m/s |
| Qmax | 1.2e7 W/m3 |
| Bmax | 约0.007 T |
| sigma场最大 | 约27 S/m |
| 顶面/侧壁/底部散热 | 30.8/20.5/3.7 kW |

流动：内部外围热区上升、主环流在轴线闭合、最外层冷壁下降形成较小外围环流。不是两个简单上下等大的涡。

Marangoni：d(gamma_surface)/dT=-1e-4 N/(m K)。正文比较55kW gamma=0与-1e-4及30kW。无Marangoni分支必须比较gamma=0曲线，不要求硬凑到含Marangoni的速度峰值。若当前没有可信实现，先做纯浮力，不通过增大beta补偿缺失项。

附录3：固定总吸收功率时，sigma整体±15%的该组结果变化<0.3%；低温k变化可导致局部约100K、轴线附近速度约100%差异。这是特定算例的敏感性，不证明sigma温变形式或热源分布不重要。

S图2.1主要是Uox2，不能覆盖J第三章R7T7/Uox1。S第109–112页EREBUS算例引用J物性：45kW稳定；55kW三维失稳，含热电偶、由已有解初始化，该处电磁100s更新、流体1s。不要与期刊稳态或其他文章更新间隔混用。

## 5. Cursor执行步骤：避免重复工作

### P0 冻结与接入

1. 记录当前commit、dirty diff、物理基线和软件交接版本。
2. 导入本包reference_data目录，保留raw CSV和metadata；不要自动重提取。
3. 读取所有CSV并检查温度递增、正值、单位及有效范围。重复温度、异常值必须报告。
4. 新建数据表适配层，不改公共流体核、不影响原A-phi分支。
5. 保存raw表与solver表两个层次；检查mu对数轴、cp陡变和k平台。
6. 输出mu/k/cp/sigma随T的独立图，对照本包QA。图像噪声修整不能当物性标定。

### P1 建立明确的配置

保留PAPER_NC_BASELINE；新增JAC_CH3_55KW_NO_MARANGONI和JAC_CH3_55KW_WITH_MARANGONI。后者若无实现标NOT_READY。

共同参数rho2800、beta3e-5、R0.25、H0.185、f282kHz、边界按本文件。J第三章mu/k/cp来自本包；S曲线单独存档，仅用于其材料及对照。

mu低温延伸及上限采取正文1e5约束并明确标记；保持图像提取平台原值作证据。cp/k端点小缺口单列处理，禁止将sigma大范围缺口按同样方式钳位。

sigmaA与sigmaB_fit各保存独立effective配置；B必须记录哪些温区为外推。禁止改其他参数帮助某一方案胜出。

### P2 离线与固定温度电磁对照，立即可做

1. 在共同温度网格比较A/B；1473K的A应42.99、图表B约13.7，允许数字化精度范围。
2. 使用同一个固定温度场、同几何/粒径/求解设置分别求EM。
3. 若温度场仍来自Tmax1500假设图，继续保留ASSUMED，不改称论文原始场。
4. 两方案保存归一化前功率、归一化因子、归一化后功率、Qmax位置、径向/竖向功率分布、上中下/近壁区域功率占比、残差及连续性指标。
5. 总吸收功率相同不等于空间分布正确。不能只比较总功率或最大值。
6. 保持现有已验证EM开关、自感设置、底板/屏蔽开关相同；不得偷偷加入新结构。

### P3 热模块与热平衡，可与背景流排查同时进行

1. 以本包k(T)/cp(T)核对热扩散与边界能量；可以先不求流动或用轻量轴对称热模型。
2. 检查q符号、辐射温标、不同参考温度、面积权重和重复热阻。
3. 输出Pin、Ptop、Pside、Pbottom、dEth/dt；变cp使用与模型一致的焓积分。
4. 检查Pin-Pout≈dEth/dt；稳态再比较Pin≈Pout。
5. 纯导热场是初猜和热模块检查，不是法国自然对流复现通过。

### P4 238轮背景流问题单独处理

先完成实际EOS的离线能量核算，不再大扫初始化参数。确认重力势能、压缩能、黏性耗散、密度重求和/裁剪/压力重置的能量跳变。重置U不能清除全部能量储备。仅凭压力正功不宣布压力算子是唯一根因。

根据能量账本选择单项修复，再做OFF对照。不要清零速度、加无来源阻尼或改beta来掩盖背景流。P0–P3不必等待此问题全部结束；完整NC可信度须保留此限制。

### P5 流动与完整耦合

顺序：固定温度流动→给定Q的热流动→EM+sigma(T)更新。按物理特征时间、热量与动量演化决定是否继续，不直接以一次短跑没有双环流否决。

先gamma=0；有可信Marangoni实现后再gamma=-1e-4。比较温度场、轴线速度符号、内部上升带、外围下降带、涡分界与持续性。A/B比较首轮只改变sigma，不同时调mu/k/h。

B若更接近文献，可作为复现配置选择；仍不能称已证实印刷错误。两者无法区分时如实报告。

### P6 后处理

Python生成T/U/Q/mu/k/sigma的r-z图和指定剖面。明确轴方向、周向平均方式、有效覆盖、插值及掩膜。禁止跨墙插值、无数据填零再画流线、用强平滑制造闭合涡。流函数和原始矢量方向必须交叉核查。

完整二维文献场尚未提取，当前不能声称已有逐点二维误差。可先比较原图结构、原文量级与散热；一维参考剖面是下一项可独立补充的数据工作。

## 6. 通过标准与下一轮打包

- 数据接入通过：单位/轴/温区正确，来源不混用，越界有记录；不要求消除原图自相矛盾。
- A/B实施通过：相同固定T输入，唯sigma模型不同，归一化前后功率都在。
- 数值能量误差与文献误差分开。先报告真实偏差及时间/空间收敛证据，不根据图片临时发明统一10% PASS。
- FULL_REFERENCE_REPRODUCTION不能由“出现涡”“功率归一化通过”“物性表已接入”自动触发。

下轮包：commit/diff、effective配置、输入表及hash、运行命令、原始CSV/日志、图像QA、A/B比较、热平衡账本、背景流新增证据、未完成项。当前无需重做已通过的全部方腔/TEAM7；只有公共求解器有实质改变时才扩大回归。

## 7. 仍未解决的资料项

精确线圈截面与逐匝坐标；印刷sigma公式与图的真正原因；原始物性表及保密测量；一维剖面数字化；完整二维场；期刊与博士工况逐项最终映射。本包并未替作者补造这些信息。


<!-- ===== END 239_FRENCH_THESES_EXTRACTION_AND_CURSOR_PLAN.md ===== -->


<!-- ===== BEGIN 239A_CLIENT_AND_FRENCH_CASE_SCOPE_ADDENDUM_2026-11-24.md ===== -->

# 239A：甲方算例与法国参考算例的执行范围补充

日期：2026-11-24
修订：法国自然对流复现优先版；替代此前239A的并列执行顺序。
配套文件：239_FRENCH_THESES_EXTRACTION_AND_CURSOR_PLAN.md及239完整数据包。

## 1. 必须先读的范围说明

239轮所有工作均保留，但执行优先级改为：**法国自然对流复现优先，阻碍该复现的共用数值问题紧随其后；甲方独有问题与后续实验标定保留在任务清单中顺序推进。**

- 近期主交付：在可重复运行的法国参考配置中，取得可解释的温度分布、环流方向/结构、速度量级和热平衡证据。
- 法国期刊自然对流结果仍是主要目标；Jacoutot博士第三章55kW是补充参考，不因数据更全就悄悄替代期刊50kW目标。两者分别报告。
- 甲方线：保护已交接配置、软件和弹性搅拌桨；238已有轨迹优先离线归因，法国独立任务不必等甲方几何全部问题解决。
- 共用修复：根据是否影响法国参考工况决定优先级，不以“所有历史测试全部关闭”为进入法国耦合的前置条件。
- 所有旧任务保留，改为“本轮优先/有条件继续/后续保留”；不删除、不得改成已通过。
- 本轮不是两套case同时无计划长跑，也不是只整理数据。前置检查可用后应实际推进法国温度驱动流动与电磁热源耦合诊断。
- “文献复现占主要验收工作量”是当前组织工作的依据，不自动等同于已满足合同各条款或算法已全部验证。

若239正文中“共同参数”“新增配置”等文字存在歧义，按本补充执行：博士论文参数只对新建的对应法国参考配置生效，不适用于所有算例。

## 2. 本轮任务矩阵

| 对象 | 必须推进 | 本轮不自动启动 |
|---|---|---|
| 甲方自然对流case | 保存原基线；238轨迹离线归因；影响法国工况的共用修复优先，甲方独有修复后续保留 | 将全部法国材料和边界覆盖进去；重新大范围调参；长时实验温度拟合 |
| 甲方搅拌case | 保留已交接弹性桨、材质、驱动、几何和生效配置；若共享代码被修改，做最小必要短测 | 重新开发弹性桨、完整整圈运行、搅拌与EM的长时全耦合 |
| 法国期刊及博士参考case | 最高优先：配置接入、OFF/ON短测、sigma A/B、热收支，具备条件后实际温度驱动流动和耦合推进 | 无条件长时三维耦合；立即重写Marangoni或冷壁电磁模型 |
| 标准验证case | 保留既有证据；共享模块实际改变时进行对应的针对性回归 | 无依据重做全部方腔、TEAM7及其他历史测试 |
| 客户软件 | 本文件不改变软件机当前冻结基线 | 将物理试验分支自动合并进交付版本 |

“保留”不是声称物理已经通过全部验证。尤其不得将TEAM7现有进度自动改成PASS。

## 3. 甲方配置保护

首次运行前记录本机已有甲方NC与搅拌配置、文件路径、commit及dirty diff，不根据聊天摘要重新创建一套猜测配置。

必须保留：

1. 甲方STL、坐标、尺寸、液位和线圈布局。
2. 甲方已经确认的物性、边界和驱动信息。
3. 当前材料借用项及其来源。例如甲方曾允许电导率存疑时借用法国数据，不等于所有物性都应替换。
4. 搅拌的elastic_driven模式及已确认不锈钢和转速参数；不得退回刚体动画。
5. 电源输入功率与玻璃实际吸收功率的区别；不得自动把30kW电源输入当成30kW焦耳热。

239不要求将甲方NC改成直径500mm，不要求改用法国7匝570mm线圈，不要求统一55kW，更不要求将甲方实验目标改成1660K。

未知电流、实验测点、散热等仍按现有未知状态记录，不用法国数据伪装成甲方已确认值。

## 4. 法国参考配置身份

至少在现有配置系统中明确区分以下身份，名字可以按项目风格调整：

- CLIENT_NC_BASELINE：甲方自然对流。
- CLIENT_STIRRING_BASELINE：甲方弹性搅拌桨。
- FRENCH_PAPER_NC_BASELINE：原自然对流期刊参考。
- JAC_CH3_55KW_NO_MARANGONI：博士第三章纯浮力分支。
- JAC_CH3_55KW_WITH_MARANGONI：博士第三章含热毛细分支，未具备可信实现时标NOT_READY。

不能把这些配置写成同一文件反复覆盖。可复用同一可执行程序，不强制新建多个C++ test目录或重写框架。

博士第三章参考采用R=0.25m、H=0.185m、rho=2800kg/m3、beta=3e-5/K、f=282kHz及239列出的对应边界；这些只属于该参考身份。

期刊50kW与博士55kW的结果不交叉作为硬性验收值。Sauvage Uox2曲线独立存档，不作为Jacoutot第三章R7T7/Uox1的默认物性。

## 5. 电导率A/B本轮在哪做

首先在法国参考身份下，使用同一固定温度场对比：

- A：印刷公式，以开尔文代入。
- B：图曲线对应的拟合候选；超出可读温区的部分显式标EXTRAPOLATED。

首轮只改变sigma，保持几何、固定T、总吸收功率目标、离散分辨率、自感/屏蔽/底板开关与收敛设置相同。

输出原始功率、归一化系数、归一化后功率、Q空间分布与峰值位置、区域功率占比、残差及连续性指标。

甲方案例本轮不必再重复完整A/B矩阵。后续需要借用某一关系时，另存CLIENT_*_BORROWED_SIGMA配置，说明借用及外推；比较甲方实验，不能以法国环流图是否接近来决定甲方真实电导率。

## 6. 238轮数值归因如何延续

238轮原配置继续保留用于归因，不能一边改几何/物性，一边宣称旧问题已定位。优先离线核算已有轨迹：动能、重力势能、与实际EOS一致的压缩能、黏性耗散、密度重求和/裁剪/压力重置事件的能量变化。

同时在法国参考几何及明确配置下做同类OFF检查，区分：

1. **法国参考中也存在影响目标流动的背景运动**：优先定位及修复；温变物性接入、固定T电磁和离线热检查继续推进。未解决前不宣称真实环流已经通过。
2. **问题主要局限于甲方几何/边界设置**：保留证据与待办，不把甲方该项彻底关闭作为法国流动试跑前置条件。
3. **法国参考仅有较弱且衰减的初始化扰动**：记录量级和衰减，结合ON响应及时间窗评估污染，不要求速度机器零或无限延长静水预处理。

OFF/ON比较必须注明温度、mu、重力、初始化和边界是否相同；若需要改变热边界以维持固定T，明确说明对照范围。至少输出速度范数/高分位值、轴向有向剖面、动能、质心变化、扰动空间位置和随时间演化。若信号无法从背景区分，结论是“暂不能归因于浮力”，不是通过。

不要因出现流线就忽略背景；也不要在证据尚不能影响法国主线时无限抠甲方OFF细节。根据账本选择单项机制修复，不清零速度、不加无来源阻尼凑图。公共模块确有修改时做相应回归，软件冻结版不自动更新。

## 7. 执行顺序与本轮终点：法国主线优先

### 7.1 本轮优先执行

1. 冻结甲方配置；核对法国期刊50kW与博士55kW两套身份。列出几何、物性、热边界、功率、Marangoni和对比图，不混用。物性未知项继续注明借用来源。
2. 接入239已提取曲线并检查单位、温区和插值/外推。先使用Jacoutot相应物性，Sauvage曲线保留对应材料身份。不重新从零提取。
3. 法国配置做OFF短测及固定T浮力ON响应；同时进行固定温度下电导率A/B、电磁空间热源与热边界检查。238已有轨迹离线归因作为辅助。
4. 前述检查若未暴露明显实现错误，即实际推进法国给定Q的热流动，再进入EM+sigma(T)更新的耦合诊断；不要只交“接入成功”便结束本轮。
5. 先比较纯浮力对应参考；已有可信Marangoni时再比较含热毛细参考。缺失项及适用范围写明，不能靠beta补偿。
6. 根据温度场、热收支、热源位置和流动时间序列安排分段续跑或有依据的单项修复。没有持续异常且尚未形成目标热场时，可顺序延长；禁止默认时间越短越好，也不无检查直接长跑。
7. Python输出同坐标温度、速度矢量、Q及物性图。优先判断热点位置、内部上升带、轴线回流、外围下降带，再判断速度量级与散热分配；标准算例准确性仍独立保留。

运行时长由有效时间尺度及演化决定，不能把“需要物理时间”当作无进展，也不能以无限延时替代原因分析。每一段必须说明热场是否仍在变化、背景是否污染、为何继续或停止。

### 7.2 有条件继续

- 法国工况暴露共用数值缺陷时，优先修复并做必要回归。
- 共用修复影响甲方交接算例时，补最小必要短测；不把全部甲方矩阵重新跑一遍。
- 电导率A/B先用固定T区分输入影响；后续耦合若差异明显再逐项比较。不同时调mu、k、h和功率帮助某方案胜出。

### 7.3 后续保留，未取消

- 甲方几何独有的边界/背景流细节，以及实验温度标定。
- 甲方完整整圈与长时电磁-热-搅拌耦合。
- 法国搅拌对应工况的系统复现，优先级在自然对流主线之后。
- 非当前阻塞项的TEAM7等遗留验证；保留真实未通过状态。
- 文献一维剖面和二维场进一步数字化，未交付部分不得猜测补齐。

甲方后续允许对不确定有效换热、吸收功率或存疑物性进行标定；需独立配置、记录依据、比较实验。已确认几何与测点不能为凑图随意改，标准验证算例不能这样标定。

### 7.4 本轮必须给出的结论

本轮重点是“法国复现推进到哪一步，以及下一项真正阻塞是什么”，不是仅交一批输入图。至少报告：

- 法国OFF/ON对照是否能区分浮力响应，若不能，具体污染源与下一项单因素实验。
- sigma A/B实际运行情况及Q空间差异。
- 热边界是否生效、温度场和散热是否具有合理演化。
- 法国流动/耦合是否已经运行；未运行必须说明具体实现阻塞，不能只说甲方待办未清完。
- 环流方向/结构、热点与冷层、速度量级及收支的适用范围。
- 甲方后续工作保留清单及共用代码影响范围。

不要求本轮无论结果如何都宣称双环流复现，但也不允许在前置项可用后自行停止在数据接入阶段。

## 8. 报告要求

每个结果必须携带case_identity、geometry_id、material_source、thermal_bc_id、power_control_mode、Marangoni开关、commit/dirty状态与effective配置。

报告分为：

- CLIENT：238问题定位、修复及回归；搅拌原配置保护状态。
- FRENCH：输入接入、sigma A/B、热收支及图像对照。
- SHARED：实际修改的通用代码及影响范围。
- OPEN：尚未解决项及下一步单项实验。

物性数据、固定T电磁、热平衡、OFF背景流、环流趋势、完整复现分别给状态。某项完成不自动证明其他项通过。

文献一维剖面与完整二维场目前未完成数字化。本轮可用已有原图、方向、量级和热功率表做明确范围的比较；不得声称已有完整逐点二维误差，不要求Cursor用猜测填满数据。

## 9. 给Cursor的一句话指令

按239执行，239A修订版优先约束执行顺序：法国自然对流复现优先，接入数据后实际推进法国OFF/ON、电磁A/B、热平衡与分阶段耦合；共用缺陷按是否阻碍法国主线处理。甲方和其他验证任务全部保留，但甲方独有问题不再成为法国工作的统一前置条件。不得覆盖甲方原配置，不重做已交付提取，不无计划开启双线长跑。


<!-- ===== END 239A_CLIENT_AND_FRENCH_CASE_SCOPE_ADDENDUM_2026-11-24.md ===== -->


<!-- ===== BEGIN 240_OPHELIE_FLOW_THERMAL_CLIENT_HANDOFF_AND_LINUX_BRANCH_CURSOR_PLAN_2026-11-24.md ===== -->

# 240：法国低速自然对流突破、甲方两案补齐与Linux软件分支交接

日期：2026-11-24
执行对象：物理开发机上的Cursor；第12节为Linux软件机交接约定。
前置资料：Doc239 lean交接包、239物性数据包、修订版239A、238本地完整证据。

## 0. 本轮决定与优先级

用户最新要求：法国复现仍是验收主要物理目标；同时补齐甲方自然对流和弹性搅拌两案的软件调用闭环，用户另开分支开发Linux软件。

本轮四个结果：

1. 在法国固定域模型中降低非浮力背景运动，获得可分辨的浮力响应。
2. 让局部k(T)、cp(T)、mu(T)进入真正的计算路径，建立合理冷壁温度初猜。
3. 非均匀T下完成sigma A/B，使用明确的功率控制，检查EM热源传递时序。
4. 交付甲方两个真实可运行算例、配置、输出和重启证据，使Linux软件开发可以立即开始。

这不是重写SPHinXsys，也不是把两个工程case全部调到定量验收精度后才允许开发GUI。软件流程可用、数值模块验证、法国复现和甲方实验温度对齐分别记录。

甲方两案交接准备在本轮开始即开展，不等待法国主线全部结束；GPU运行按依赖顺序排队，避免争抢显存导致失败。

本文件接替239的下一轮执行计划；239数据和全部未完成事项保留。独立A-phi/Laplace研究不改。没有用户取消授权的必要回归/可逆修复，Cursor自主推进，不在每一步完成后反复要求用户重新批准或手工触发下一命令。

## 1. 已审查的239事实，禁止错误归因

Doc239记录的HEAD为845e580de9ad13fc49eda78add5586531b811ab1，分支feature/EM-Biot-Savart；存在WIP，HEAD并不等于所有运行代码。以本机实际状态再次冻结，不假装本轮包包含完整源代码。

| 事实 | 本轮意义 |
|---|---|
| 法国OFF 1s U_rms_end=0.0290414m/s，Umax=0.128176m/s | 背景足以污染约1e-3m/s的参考流动 |
| 对应ON U_rms_end=0.0292271m/s | 浮力信号尚不能由速度范数分辨；不说明浮力力链必然错误 |
| hydro_model=open_free_surface，paper_flat_slip_top=0，total_p_plus_GravityForceCK | 与期刊固定平坦滑移顶面存在模型差异 |
| constant_thermal_rheology=1；P6 mu_min=max=2.39045，cp约1598.77，k约6 | 物性表导入不是完整局部温变物性生效 |
| P2为1473K均匀T下A/B | 接线已查；不能替代非均匀T热源形状比较 |
| 配置描述fixed_absorbed_power，P6实际fixed-current | 必须统一实际生效控制与报告 |
| P6 EM重建功率与流体热量功率差约0.7%，存在交替变化 | 需要查量的时间戳、映射开关与取样路径；不能直接断言漏能 |
| 初始均温1473K、Jac BC下散热约221kW | 初态冷却合理；不是必须新增冷壳厚度的证明 |
| P5 french_energy_budget.csv只有表头；P6有31行数据 | P5仅可引用摘要及日志量，不能声称该CSV支持时间积分 |
| 多个summary passed=1而validation_passed=0 | 技术短跑成功不能写成物理验证通过 |

包未提供本轮全部源码、粒子事件、kinetic_work.csv，故尚未证明具体哪个算子是唯一根因。不要从“声学步骤后速度出现”跳到“已证明压力公式坏了”。

## 2. 配置隔离与本轮执行表

保留以下身份，允许共用可执行文件：

| 身份 | 用途 | 不得混入 |
|---|---|---|
| CLIENT_NC_BASELINE | 甲方几何自然对流，软件可运行基线 | 法国500mm/7匝/55kW未经说明覆盖 |
| CLIENT_STIRRING_ELASTIC_BASELINE | 甲方不锈钢弹性桨及既有驱动 | 刚体动画替代真实弹性求解 |
| FRENCH_PAPER_NC | 期刊50kW目标及其边界 | 博士55kW结果当成同一硬性验收值 |
| JAC_CH3_55KW | 已知物性较完整的博士参考，分gamma=0/含Marangoni | Sauvage Uox2与Jac Uox1物性混用 |
| FRENCH_FIXED_DOMAIN_DIAGNOSTIC | 固定平坦滑移顶面/压力诊断 | 作为甲方真实自由液面模式默认覆盖 |

每次运行写出effective配置和几何/物性表hash；配置中的PENDING_WIRE、声明值与实际生效值必须被程序验证，不得将手写JSON当作执行证据。

## 3. P0：第一时间冻结甲方软件可用基线

### 3.1 只做一次本机盘点

记录git status、HEAD、分支、WIP diff、实际目标与二进制hash、运行所需STL/reload/物性表/相对路径。不得reset --hard、clean -fd、强推或git add整个仓库。

确认每个甲方case最后一次可运行版本及参数；找不到证据则标UNVERIFIED并安排短测，不按聊天记忆重建假配置。

### 3.2 立即形成软件交接候选

- 复用既有validate/prepare/run、JSON和CLI约定。先查本机实现与软件计划，不凭空新设计另一套协议。
- 若只有现有binary，可用薄包装层调用，但参数检查、退出状态、输出目录、日志路径必须真实可用。
- 未有GUI不阻止物理侧输出求解器交接包。
- 产生基线提交时，只提交明确列出的源码、构建文件、配置、表和小型必需资产；大体积结果另打包。
- 脏目录不允许直接共享给第二台机并称其为版本。要提交可审查改动，或交付base+完整补丁及补丁hash；优先可直接构建的提交。

P0完成后软件机即可开始，物理机继续P1–P6。尚未通过真实短跑的功能在交接状态中明确标注。

## 4. P1：将流动诊断收敛为三种配置，而非全开关扫描

### 4.1 先审查现有实现

查清paper_flat_slip_top、reduced_pressure/pprime、既有静水背景场、Riemann边界及TVF相关代码。阅读调用顺序、内外算子和墙接触，不只看开关名称。

不要因为旧开关存在就认定已正确实现；也不要在没有审查前再写一个重复版本。保留原始低耗散Riemann主框架，先在OPHELIE局部配置/适配类隔离。

### 4.2 有限对照矩阵

| 配置 | 顶面 | 压力/重力 | 目的 |
|---|---|---|---|
| F0 | 239开放FS | 原总压+重力 | 已有对照；证据足够时不再重复 |
| F1 | 固定平坦不可穿透、切向滑移 | 保持原压力/重力 | 只改变顶面，识别其影响 |
| F2 | 与F1相同 | 一致的静水分解或扰动压力路径 | F1仍污染时才做；定位基准重力平衡影响 |

除上述变量外，F1与F2分别保持几何、粒子排布、温度、mu、步长策略和输出一致。首先冻结T0常数物性检查动量，不同时加入大范围温变mu，避免把刚度与边界错误混为一谈。

### 4.3 固定平坦滑移顶面的边界含义

- 顶面位置固定；法向速度为0；侧/底仍按目标文献无滑移。
- gamma=0时顶面切向黏性牵引为0，不是切向速度为0。
- 顶面机械条件与顶面对流/辐射分别施加，不能加机械顶面就取消散热。
- 使用现有墙粒子/边界能力，默认不引入需重新布置ghost粒子的工程流程。
- 不能仅每步剪裁z、清零整个速度或施加经验dashpot，冒充正确边界。
- 必要的边界法向投影必须与边界离散一致，并记录反力/能量变化；不将全域非物理耗散藏进边界。
- 核对原自由液面识别、密度修正、压力大气钳位是否仍错误作用在固定顶面；检查角区重复施加。
- 增加顶面切向均匀滑动的针对性检查，确认没有把free-slip误做no-slip。

### 4.4 F2压力分解的实施边界

连续模型：p=p_h+p'，grad(p_h)=rho0*g。参考rho0重力被静水背景平衡，温变浮力保留。

若采用Boussinesq固定域扰动压力形式：明确rho是声学辅助密度，rho0是基准惯性密度，beta(T-Tref)是浮力密度差；EOS返回的是p'还是p_total，必须一致。

必须一起审查：

1. 内部压力梯度与墙压力边界。
2. 连续性、密度重求和及密度/压力裁剪。
3. Riemann左右压力及速度状态。
4. 静水背景是否已加在某处，避免双计或漏计。
5. 顶面/侧壁/底部基准压力/牵引。
6. 诊断压力、真实总压力和施力压力的字段含义。
7. 初始化及每次重启的状态。

不能只关闭GravityForceCK而保留原静水密度场，再称为等价分解；不能将上一步实测残余力永久反向补偿。若完整分解不具备条件，只做清楚标记的隔离诊断，不发布生产正确性结论。

### 4.5 静态残余与OFF/ON运行

先在未移动粒子上输出压力/重力/墙作用与净加速度，按内部、侧壁、底部、顶面、角区分区。输出质量加权范数、有向分量和高分位值；不是再次只看全域合力。

OFF检查：U=0，等温，无热源，无热边界演化，浮力关；t=0、首声学步、0.2s、1s记录；若稳定再延长到5s确认无持续增长。

ON检查：首先使用与239相同的冻结径向热核DT30作为力链回归，与对应OFF同物性/同时间窗；温度非均匀但mu暂固定，其他设置一致。径向热核可导致中心上升，不能要求它复现法国外围热区上升的拓扑。

通过基础响应后，使用带外围内部热区的同一固定T场做OFF/ON。该温度场的来源/假设明确；不把人为初猜称为文献原始场。

不必每个失败配置都继续跑5s；1s已显示持续增长或明显污染时保存证据，转到下一项有依据的实现检查。

### 4.6 本轮工程诊断门槛

以下是本轮内部筛选标准，不是文献或合同规定：参考速度尺度Uref=1e-3m/s。

- 进入可信环流评价前，选定时间窗OFF质量加权U_rms应不超过1e-4m/s，且持续无增长；U的99分位原则上不超过3e-4m/s。若局部角点异常，单独报告，不删除异常点后宣布通过。
- 对应ON应出现与温度梯度相符、持续存在的有向响应；主比较剖面信号建议至少为同位置OFF背景的5倍。OFF接近零时报告绝对差，不用除零制造无限SNR。
- 未达到门槛仍可做探索性热流动，但不得标FLOW_VERIFIED；须说明背景相对目标量级。
- 速度范数ON-OFF不能代替逐位置/同采样算子的矢量与有向剖面对比。
- 任意清零/投影/边界约束能量变化应列入预算。

F1可用则不必强制开发F2。F1、F2都失败时停止新增lid/TVF/密度开关排列，交付具体算子与首步事件证据，进入单项修复，不回到无目的参数扫描。

## 5. P2：真正接入局部温变热物性，先固定粒子

### 5.1 必须看实际算子而不是材料标签

建立property_path_audit.md，列明每个字段生产者、存储位置、消费算子、更新时刻、CPU/GPU路径：

- k_i(T_i)：粒子内导热及边界处理。
- cp_i(T_i)/h_i(T_i)：温度推进、源项、热边界和能量积分。
- mu_i(T_i)：黏性力、时间步及隐式求解。
- sigma_i(T_i)：EM求解、松弛与更新。

日志与CSV同时输出真实最小/均值/最大值；若T场有明显跨温区变化而这些始终等于T0值，则接入未完成。平均物性不是局部变物性。

### 5.2 热离散与能量

固定粒子，保持质量和体积定义一致，先关闭动量更新，隔离热模块；不是仅zero-buoyancy，因为239即使浮力关粒子仍有显著运动。

热方程可写为m_i dh_i/dt=sum_j Q_ij+Qsource_i+Qboundary_i，h(T)=integral cp(T)dT。粒子对通量应反对称，使用与既有验证一致的正导通系数；不要只按邻居算术平均mu/k随意改算子。

边界热流产生的能量和粒子实际焓增必须同账。变cp不能在扩散用cp(T)而边界温升用cp(T0)。若暂采用温度形式，需要以相同积分/迭代精度验证焓收支。

只做必要核对：

1. 均匀T、绝热、无源，应不产生温差与净能量变化。
2. 非均匀T、绝热、无源，内交换守恒，温度趋势合理。
3. 给定Q及Robin/辐射，独立积分热源和热损与焓变化闭合。
4. 对本轮新变系数实现做一次针对性的两层导热/串联热阻核对；若复用已有已验证实现则引用对应证据，避免镜像实现的重复测试。

### 5.3 边界原文值与初态含义

Jac第三章：h_side=h_bottom=300，h_top=20 W/(m2 K)，Twater=Tair=303K，Trad=423K，eps=0.9，R=.25m、H=.185m，吸收55kW。期刊目标单独保存。

均温1473K时理论裸壁散热约222kW是合理瞬态条件。1–2s冷却不是失败，也不是必须引入未知skull厚度的证据。

对应Jac稳态散热20.5kW侧壁、3.7kW底部，以该h反算面积平均表面温度约538K和366K，仅是诊断反算，不能硬施加为Dirichlet条件来制造结果。

全玻璃域含冷玻璃时，由k(T)、mu(T)、热边界及能量演化形成冷层；本轮不增加串联skull热阻以补偿初态。若已有case实际裁去了冷区，应先明确域边界后再定义等效模型。

## 6. P3：快速生成冷壁温度初猜，不等待三维声学跑到热平衡

### 6.1 优先实现顺序

先复用已有固定粒子隐式热求解。若没有且显式受限，则允许隔离的轻量轴对称有限体积/有限差分热计算生成初猜，不能替代SPH最终验证，也不扩展成完整新CFD平台。

轴对称模型应包含(1/r)*d(r*k*dT/dr)/dr及竖直导热；r=0按对称通量，r-z能量积分用2*pi*r体积/面积权重。不能将平面二维方腔当轴对称圆柱。

输入Q依次采用可追溯给定Q和当前EM的轴对称Q；归一化目标功率与几何体积一致。非均匀Q不能在整个温度更新中又均匀化。

非线性k、cp及辐射采用有残差检查的迭代/隐式步进，选择可实施方法即可，不强制一个求解器库。输出稳态残差或真实伪时间，不将伪时间说成实验升温时间。

### 6.2 适用范围与收敛

纯导热结果只作初猜/热账参考。缺少对流时Tmax和热点位置可能不同于文献，不能硬凑1660K或把30.8/20.5/3.7kW设成强制分配。

输出T(r,z)、边界T分布、Ptop/Pside/Pbottom、Pin、焓或稳态残差。映射回SPH后重新核算边界和热源，不假设参考模型与SPH的热账完全一致。

初猜收敛的建议内部标准：连续两个诊断窗口的归一化能量不平衡<1%，Tmax及质量平均T变化各<0.5%；做一次时间步减半或网格加密，确认主要结果差异<2%，否则报告尚未收敛。这些是初猜质量检查，不是法国复现验收。

在冷区特征尺度上至少能辨别温度梯度；dp15mm粗粒子无法声称精确解析毫米冷层。先使用参考热网格生成平滑初猜并保存映射误差，再判断是否确需SPH局部/全局加密，不在物理路径未正确前扫五个分辨率。

禁止无检查把纯导热高温解用于生产：若越过物性有效范围、Boussinesq假设不再适用或热量不平衡，先报告并调整初猜求解，不通过温度硬裁剪隐藏失败。

## 7. P4：局部mu(T)与低速流动重连

1. 优先复用已经验证过的局部变mu及隐式黏性实现；核对当前frozen_q生产入口是否真正调用它，而不是其他入口通过后自动宣布迁移完成。
2. GPU隐式优先；若本轮确实不能可靠完成，允许CPU隐式黏性或已有准确路径作为隔离回退，报告设备、性能及数据传输。准确性优先于全部GPU，禁止用大幅减小mu代替。
3. 先运行常mu与统一粒子mu字段等价性；已有证据适用时引用。再用239/227已验证的变黏度Couette等针对性检查，不要求全部基准重跑。
4. 用P3温度初猜计算mu分布，输出mu_min/max、达上限粒子比例、mu>1e4区域形状与有效流动区域。
5. mu>1e4是文献冷层分析阈值，不自动冻结/删除粒子。若采用冻结等额外模型需独立身份，本轮不默认实施。
6. 真实时间步与黏性更新时间一并审计；隐式并非允许任意大步长，仍检查时间步减半的关键剖面及耗散趋势。
7. 在低背景固定域上先固定T求流动，再解热流动，最后开EM更新。恢复移动粒子时，明确温度/物性更新与cell/邻域重建顺序；不要重复引入已修过的邻域年龄问题。

Jac图mu低温平台约8e4与正文1e5上限的差异保持记录；本轮选定一个有来源配置，不为此扩大扫描。Sauvage Uox2物性不混入。

## 8. P5：非均匀T电导率A/B与恒吸收功率

### 8.1 控制模式必须真实一致

法国复现选fixed_absorbed_power，每次EM重求解后吸收功率目标仍为55kW（博士）或对应期刊50kW。线性、材料固定的谐波EM求解可通过Inew=Iold*sqrt(Ptarget/Praw)调整激励，但复数J、场及Q应一致缩放，记录所有因子。若模型不满足该线性前提，按实际求解处理。

fixed-current保留为独立诊断身份，不用它冒充恒吸收功率。功率控制与更换sigma是两件事，先明确模式，再做A/B。

### 8.2 A/B输入

- A：10**(3.7921-3179.8/T_K)。
- B：使用239附带的图拟合，明确拟合区1100–1550K与外推区域。
- 同一个非均匀温度场；先用P3初猜或已有明确ASSUMED场，不能两方案各自选不同T初态。
- 其余geometry、粒子、底板/屏蔽、自感、容差、功率目标一致。
- 均匀T0下已有接线不重复；本轮关注形状差异与热点位置。
- 输出sigma值和外推粒子比例/体积分数；冷区sigma跨多数量级要检查线性求解残差与预条件，不用改变公式静默削平。

### 8.3 Q传递与交替偏差

给每个量加physical_time、em_generation、temperature_generation、particle_state_generation，至少同时区分：EM原始总功率、EM网格积分、流体实际施加瞬时功率、该热子步实际能量。

检查periodic_refresh时conservative_remap是否从1变0；这可能是模式标志或真实路径变化，须审源码，不能只凭CSV猜测。核查日志交替采样是否一条对应诊断、一条对应实际施加，排除重复施加与阶段错位。

同一状态下固定Q、固定粒子先确认传递。随后移动粒子检查；必要的守恒重映射要保留局部形状并记录修正幅度。

本轮建议固定状态传递相对误差<0.1%，移动状态累计热源能量偏差<1%；这是内部调试指标，达到不代表EM空间场验证完成。记录相对EM功率和相对目标功率的两个偏差。

P6末端恰好触发EM刷新却没有后续热子步时，transfer标NOT_MEASURED_NEXT_THERMAL_STEP，不因为没有下一步数据而判传递已通过；可多走一个必要热子步补测。

### 8.4 后续耦合

先用诊断短间隔确认更新，再根据温度/sigma变化选择更新间隔，不默认0.5s或100s是普遍正确值。初猜之后允许分段续跑；对照时间间隔减半对主要结果的影响，避免无意义地每声学步解EM。

## 9. P6：甲方自然对流与弹性搅拌两案必须完成的交接

### 9.1 甲方自然对流

- 原甲方STL、坐标、液位、线圈配置及材料来源保护。
- 区分电源输入、玻璃吸收、真实EM热源、给定体积热源；用户界面应能显示来源。
- 跑一次有效短窗，覆盖至少一次热更新、一次实际EM更新（如果该模式承诺周期EM）以及结果写出。
- 记录温度、速度、质量/逸出、热收支、求解残差和背景限制。不能把“未崩溃”标为实验温度1000C验证通过。
- 如共用流动修复需用于甲方，独立配置回归；法国固定顶面不能未经说明覆盖甲方原物理边界。

### 9.2 甲方弹性搅拌

- 核对真实elastic_driven，不锈钢参数和既有10rpm等已确认驱动；以本机effective文件为准，不重新猜几何或材质。
- 确认固体内部应变/应力求解实际运行，驱动施加于既定驱动区，非全桨刚体坐标动画。
- 短窗应覆盖固体子步、流固载荷交换及写出，报告驱动角度和扭矩/力/位移/应力的适用量。
- 若短窗不足分辨弹性形变，标“弹性路径运行”，不为展示夸大形变或降低E；由独立固体基准提供正确性证据。
- 若热源用替代模型就标替代模型；真实EM/完整FSI热耦合未运行则保留未验证状态。
- 本轮不强制一整圈，不因此停止软件交接。

### 9.3 重启与可调用性

每案至少一次短跑+结果加载；若软件承诺restart，则做保存-重启的短窗核对，保存粒子ID/排序映射、温度、密度/压力、固体形变/驱动相位、时间和必要EM状态。

若当前还没有可信restart，不伪造支持；标功能缺口，软件先禁用该按钮，提供从初态重跑。不能为了“完整”发布虚假的继续计算。

正常结束有非空summary、结果帧及CSV；错误配置有明确非零退出码或结构化失败，不能输出一个passed=1。每案都提供可直接复制运行的命令，无交互相对路径依赖。

本轮短测时长按子步/刷新覆盖确定并写明，不随意指定1s作为所有case充分证据。可复用仍匹配当前binary/config hash的既有有效结果。

## 10. Python后处理与判据

统一输出T/U/Q/mu/k/cp/sigma的r-z图；甲方几何不轴对称时保留真实剖面与三维结果，不强做轴对称。

- 标注原始粒子与插值网格、掩膜、周向平均、有效覆盖率、坐标轴及重力方向。
- 边界附近不跨墙插值，无数据不填零制造流线。
- 速度矢量与流函数重建交叉核对；若散度/覆盖不满足，流函数标近似，不能据闭合线宣称真双涡。
- 输出轴线、内部上升带、外壁下降带的有向速度及时间序列，固定探针和空间采样定义。
- 法国热核初猜与文献温度分布不同，需展示，不用相同颜色范围掩盖数值差异。
- 期刊与博士、含Marangoni与纯浮力分开比较。Marangoni未实现不是本轮重写大模块的授权，也不掩盖差异。

文献一维剖面与二维场未在239数据包完成数字化。当前原图、方向、量级、散热可用，不能声称逐点二维误差或强行填满图像数据。本轮参考提取若有独立可用资源可以补，但不阻塞其他工作。

## 11. 本轮最小运行清单与停止规则

| ID | 必须/条件 | 核心问题 |
|---|---|---|
| B0 | 必须 | 甲方配置与软件基线冻结、资产完整 |
| F1_OFF | 必须 | 固定滑移顶面是否减小背景 |
| F2_OFF | F1失败且分解实现具备条件时 | 静水分解能否改善平衡 |
| F_ON | 选定可用配置后 | 浮力响应能否与OFF分离 |
| T_FIXED | 必须 | 局部k/cp热路径与能量一致 |
| T_INIT | 必须推进 | 冷壁温度初猜形成，报告有效温区 |
| MU_LOCAL | 必须推进 | 局部mu进入实际入口，刚度处理明确 |
| EM_AB_NONUNIFORM | 必须 | sigma温变形式对空间热源影响 |
| CLIENT_NC_SMOKE | 必须 | 甲方自然对流真实调用/结果闭环 |
| CLIENT_ELASTIC_SMOKE | 必须 | 甲方弹性桨真实调用/结果闭环 |
| F_COUPLED | 前置通过后 | 实际向法国耦合流动推进 |

相同编译/输入错误先修复再跑，不靠降低门槛继续。失败需保存一个足以定位的事件包，不保存十套重复大输出。

不进行dp×CFL×TVF×lid×mu×sigma的笛卡尔扫描。先固定dp15mm诊断；只有热层分辨/收敛证据需要时加一个更细分辨率，不把它变成无穷研究。

任何阻塞不会自动阻止甲方软件交接：可运行路径冻结、未达物理标准的标签保留；但不可将有数值崩溃的路径伪装为默认可用。

## 12. Git分支与Linux软件工作约定

### 12.1 原则

物理机继续当前物理工作；软件机从明确的可运行提交创建独立Linux软件分支。分支名按仓库现状确定，建议feature/ophelie-linux-desktop与物理后续分支区分。

不要为开分支丢弃未提交修改。已有物理WIP必须明确整理成提交或补丁；不能仅给845e580de这个旧HEAD并声称包含239改动。

软件开发可以立即从已有可运行基线开始，不等240全部通过。240后产生v240候选，再选择性迁移已验证物理修复。

### 12.2 提交边界

建议独立提交：基线/配置交接、物性适配、固定顶面、静水分解、热模块、局部黏性、功率控制及Q映射、甲方回归、软件封装。实际按依赖顺序组织，不强行为每个文件一个提交。

不要让GUI改动和试验性力学改动混在同一个提交。后续合并/挑选必须包含必要依赖、表和资产，不能只拣一行修复。先审差异再合并，不强推。

两台机器独立clone/工作目录下，一台commit/push不会自动改另一台文件；另一台pull/merge/rebase/checkout等才引入改变。不要共享同一个网络工作目录同时编辑。

### 12.3 软件机首先接什么

- 求解器可执行入口或既有薄包装；validate/prepare/run/replay实际支持项。
- 两个甲方项目及可编辑参数schema/示例；默认值来源。
- STL/reload/物性表及路径规则和许可证/使用限制（沿用项目已有说明，不另造）。
- 构建环境说明、device选择、二进制与动态库要求；本轮只承诺已有验证的Ubuntu GPU路径，CPU切换支持以实际证据为准，不把Windows打包带进来。
- 成功/失败状态、进度、日志、输出目录、取消行为与结果格式。复用已有协议，不重做不兼容API。
- 真实计算和回放清晰区分；已有结果可用于开发结果视图，但不能把回放说成刚完成求解。

状态建议逐项：workflow_completed、em_ran、thermal_ran、elastic_ran、coupling_exercised、background_flow_status、reference_validation_status、experimental_validation_status。旧passed字段如需兼容，文案只能表示明确范围。

软件内提供原始数值单位和材料/热源来源，不把试验调试开关全部暴露给最终用户。后台effective配置和详细报告供审查追溯。

## 13. 证据与交接打包：精简不能删掉关键证据

交付docs/ophelie/reference/doc240_runs/（可按项目已有结构调整），至少包含：

1. ROUND240_REVIEW.md：结论、变化、阻塞、下一步。
2. CASE_IDENTITY_TABLE.csv及每次effective配置。
3. SOURCE：base提交、完整WIP补丁或提交列表、关键C++/hpp、对应CMake、改动函数索引；禁止只交摘要无源码。
4. 首步静态力残余CSV、少量含稳定粒子ID的关键事件；kinetic_work完整CSV；无需所有时间帧，但决定根因的帧不能删。
5. OFF/ON同窗有向剖面及原始时间序列；修正239中固定0–5s窗口却只跑1s的问题，未覆盖窗标NOT_EVALUATED。
6. 热CSV必须包含数据行与实际末时间；表头文件标EMPTY，不能当证据。输出文件copy之前确认本次run_id，防止沿用共享output目录旧文件。
7. EM A/B和Q传递时序CSV，软件所需两案非空结果帧。
8. property_path_audit.md、局部物性分布/范围及隐式求解残差。
9. SOFTWARE_HANDOFF.md：稳定提交、依赖、命令、功能状态、已知限制。
10. Linux软件机最小输入资产和两案配置；大输出不入Git，数据包保存hash。

程序运行目录为每案每run独立路径，禁止多个作业共用output并靠覆盖backup管理结果。

摘要各项用PASS/FAIL/NOT_RUN/NOT_EVALUATED/EXPLORATORY，写出原因。不能再出现passed=1却让用户误以为validation_passed=1；末端无传递测量、能量参考原点、冷层未解析等分别解释。

## 14. 本轮交付成功的定义

至少得到：

- 两个甲方算例的软件调用交接闭环，或具体到功能/文件的未就绪项；可用版本足以让软件机开工。
- 固定域OFF/ON的一组有归因结果；若未成功，则明确到首步离散作用与下一项修复，而非只说压力噪声。
- 温变物性是否真实生效的证据，固定粒子热账及冷壁初猜的实际结果。
- 非均匀T的电导率A/B，以及恒吸收功率实际生效/尚未生效的真实报告。
- 在前置具备时，确实进入法国温度驱动或EM耦合计算，而不是所有工作停留在1s接线检查。

不承诺本轮一定完成法国所有图和实验温度误差；也不允许因为未达到全部目标就只提供泛泛总结。法国复现、甲方适配、软件功能和标准验证是四种状态，分别给证据。


<!-- ===== END 240_OPHELIE_FLOW_THERMAL_CLIENT_HANDOFF_AND_LINUX_BRANCH_CURSOR_PLAN_2026-11-24.md ===== -->

<!-- ===== BEGIN 241_OPHELIE_FRENCH_TEMPERATURE_FLOW_AND_CLIENT_ELASTIC_RELEASE_CURSOR_PLAN_2026-11-24.md ===== -->

# 241_OPHELIE_FRENCH_TEMPERATURE_FLOW_AND_CLIENT_ELASTIC_RELEASE_CURSOR_PLAN_2026-11-24.md

# 241：法国温度场—环流联合对比、甲方弹性搅拌修复与Linux软件可运行交接

日期：2026-11-24。执行对象：物理开发机上的Cursor。

本任务接替240的下一轮物理开发安排；239物性数据包、239A身份隔离要求继续有效。不要修改独立A-phi/Laplace研究分支。软件开发可同时进行，但本轮必须交付甲方两案的真实运行能力，不能仅交付禁用按钮或历史回放。

## 1. 本轮决定：先修甲方阻塞，同时以温度场组织法国复现

最终目标有两个，分别验收：

1. **甲方软件基线**：甲方自然对流可真实执行电磁—热短跑；甲方搅拌采用不锈钢弹性桨，实际运行到规定短窗并输出可展示结果。软件开发不再被粗reload失败卡住。
2. **法国论文主线**：用同次计算的温度剖面和速度场联合对比，主看期刊Fig.8温度场，辅以Fig.7环流趋势，再追踪Fig.6热源。不再以“画出两个涡”作为唯一推进门槛。

本轮顺序：冻结代码与证据 → 修甲方reload/弹性短测 → 甲方两案交接候选 → 法国历史F2回归 → 论文目标温度数据与当前场联合对比 → 修一个有证据的主问题 → 更新交接。CPU侧数据整理可穿插，单GPU上的正式计算串行排队。

若时间紧，先完成甲方两案交接和法国首组联合对比；不因此伪造法国稳态复现。不要启动新GUI、重写求解器或开展全参数扫描。

## 2. 审核过的事实与不能据此推出的结论

### 2.1 Doc240实际状态

- 法国固定域OFF：1s内速度为零；ON固定热核温差30K，1s内产生浮力响应。两组温度控制不完全相同，不是严格同温场的β对照。
- 当前ON主要采用`--prescribe-t-radial=30`，不是第223—224轮的环带高温温场；不能以当前单回流证明历史双回流能力丢失。
- 法国短耦合计算尚不能证明准稳态温度/环流已经接近论文。
- 局部μ有变化。k/cp表字段存在，不等于局部k/cp已实际进入热更新；cp/k在高温平台不变也不能据此断言接线失败。
- 甲方NC当前只完成0.2s、zero_q、EM关闭短烟，不能称真实电磁加热工作流已过。开放自由面仍有背景运动，须单列状态。
- 甲方弹性搅拌粗reload：玻璃约5255粒子、桨483粒子，短跑发散；历史较密reload桨约4332粒子有短窗PASS，不能替代当前粗reload通过。
- 搅拌入口存在强制dp=0.005的逻辑。它与粗reload的实际尺度可能不一致；这是优先核查假设，不是已经证明的唯一根因。
- 固体子步冻结solid_dt0的改动不能自动视为稳定性修复；当前仍有流体时间步塌缩。
- 摘要中极端负哨兵值可能被isfinite误判为正常；仅`finite_ok=1`不足以证明结果有效。

### 2.2 历史双环流的准确定位

第223轮F2引入人工环带温场；第224轮C3d配置延长到30s，记录为`DIAGNOSTIC_TWO_CELL_OBSERVED`。15–20、20–25、25–30s的U_rms约7.09e-4、6.81e-4、6.87e-4 m/s；末两窗轴心Uz变化约12.2%，未授稳定标签。

历史F2：R=0.25m、H=0.185m、dp=0.015m、β=1e-5，Q=0，关闭TVF，空间冻结温度：

```text
T = T0 + 200 * exp(-((r-0.15)/0.05)^2
                  -((z-z_bottom-0.13)/0.035)^2)
```

对应实现：

```text
tests/extra_source_and_tests/3d_examples/
  test_3d_ophelie_french_natural_convection_frozen_q/
  test_3d_ophelie_french_natural_convection_frozen_q.cpp

tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/
  diagnostics/electromagnetic_ophelie_boussinesq.h
函数：hostPrescribeDoc223AnnularHotRingTemperature
```

温度在粒子移动后按空间坐标重新赋值；不是随粒子携带的单次初温。该函数写温度，不直接写目标速度。

第224轮证据相对包根目录：

```text
doc224_runs/C_f2_late_windows/C_F2_ANNULAR_30s_RUN.md
doc224_runs/C_f2_late_windows/C_RESULT_30s_C3d.md
doc224_runs/C_f2_late_windows/C_F2_annular_hot_beta1e-5_30s_C3d.json
doc224_runs/C_f2_late_windows/C_F2_annular_hot_beta1e-5_30s_C3d_meridional_ur_uz.csv
src_snapshot_N224/
```

当前240入口仍保留环带函数调用。历史包找不到时先定位本地留存；确实没有再列出具体缺件，不按聊天描述重新制造“历史原始数据”。

## 3. 论文对比目标与身份，启动前一次性写清

### 3.1 图号对应的物理含义

| 目标 | 内容 | 本轮用途 |
|---|---|---|
| 期刊Fig.6 | 耦合收敛后的焦耳热源竖直剖面 | 判断高热源区域的位置、范围和热量分配 |
| 期刊Fig.7 | 2D轴对称耦合流动剖面 | 轴降—内部环带升—冷壁附近降，主次回流趋势 |
| 期刊Fig.8 | 2D轴对称耦合的无量纲温度竖直剖面 | 本轮法国复现首要场对比 |
| 期刊Fig.9/10 | 含真实热电偶装置的3D速度/温度剖面 | 当前无热电偶模型不得混用为同一目标 |
| 期刊Fig.12–14 | 热电偶温度数值与实验比较 | 后续独立定量验证，不能用Fig.8代替 |

Fig.8是穿过轴线的竖直截面，不是中间高度的水平截面。期刊结果为稳态耦合解；约20次耦合迭代不是20秒物理时间，不能据此设定结束时间。

### 3.2 配置至少保留四个独立身份

- `CLIENT_NC_BASELINE`：甲方几何、甲方已确认参数和独立标注的借用物性。
- `CLIENT_STIRRING_ELASTIC_BASELINE`：甲方弹性桨，不锈钢物性、10rpm既有驱动。
- `FRENCH_F2_HISTORICAL_DIAGNOSTIC`：第224轮人工环带温场回归。
- `FRENCH_PAPER_NC`：期刊50kW及其边界；博士资料只补未明参数，标记来源。

Jac博士55kW无Marangoni/含Marangoni另列身份，不与期刊50kW共用“论文PASS”。期刊与博士的ρ、β、底部h等不得不加说明地拼装。所有未明项保留ASSUMED及来源。

## 4. P0：冻结当前基线，修复“假正常”诊断

1. 记录HEAD、分支、git status、WIP diff、可执行文件hash、构建命令、设备与驱动；代码状态以实际运行版本为准。
2. 每组独立目录保存请求配置和程序输出的effective配置。旧输出不得覆盖。
3. 修复有效性检测：检查所有粒子位置、速度、密度、温度、固体变形/应力及关键规约量；同时检查计数>0、min<=max、哨兵值、非正密度/体积/时间步、NaN/Inf。
4. `Umax`必须非负，未初始化不得用极端有限值冒充有效。无样本写null+原因，JSON不写NaN/Infinity。
5. 运行失败写非零退出码和fail_reason，不能只看日志最后一行或进程正常退出。长任务使用已有进度/超时机制，保存最后有效状态。

完成条件：现有失败样本会被识别为FAIL，正常样本仍正常；这些检查是防止软件展示坏数据，不是要求重做所有模块验证。

## 5. P1：甲方弹性搅拌——先解决reload与数值尺度一致性

### 5.1 先测量后改代码

对当前粗reload和历史PASS密reload，各输出一份`reload_audit.json`：

- 粒子数、边界框、坐标单位、转轴与中心、STL/reload文件hash。
- 读取的Vol/Mass/ρ，按材料与内部/表层分组；立方根Vol仅作间距参考。
- 内部粒子最近邻距离分布及h、核支持半径；实际dp、系统参考dp、固体参考dp、reload生成dp。
- 初始内邻域/流固接触邻域数分布、零邻居比例、最小间隙和明显重叠。
- 变形梯度修正矩阵奇异/病态计数；若当前实现没有该量，明确是否使用及等效诊断。
- 实际ρ_s、E、ν、固体声速与步长公式输入；驱动约束覆盖的粒子数及其空间区域。

必须追踪`applyClientStirringDevelopmentDefaults`前后dp，以及强制0.005覆盖逻辑。不要用总粒子数反推精确dp，不要仅修改JSON中的dp标签。

### 5.2 正式修复选择

优先选实际reload尺度与求解参数完全一致的路径：

A. 当前粗reload确为合法较粗排布，则根据其生成元数据设置对应系统dp/h/体积/邻域，移除无条件强制0.005。

B. 粗reload不适合桨厚度或已缺失关键结构，则重新生成一致的粒子集；允许沿用历史较密可运行排布作为交付基线，但必须是甲方当前STL/坐标/材质，并重新真实运行。

禁止缩放粒子体积掩盖空间间距不匹配，禁止只调h覆盖稀疏/断开的桨。粗分辨率难以描述薄桨时直接说明，交付采用能可靠运行的密度，不必坚持483粒子版本。

加载时显式检查元数据；缺失元数据可做一次审计后生成manifest。矛盾时prepare修复或validate报错，不静默覆盖。

### 5.3 三步定位，避免弹性与接触同时混淆

| 运行 | 条件 | 目的 |
|---|---|---|
| S0 | 相同弹性桨，静止驱动，单固体/已有最小独立入口 | 排除初始化、材料、固体内部积分错误 |
| S1 | 相同桨，10rpm平滑驱动，暂不施加流体接触载荷 | 检查真实弹性体及驱动区域是否稳定 |
| S2 | 同一桨及驱动，恢复等温流固接触 | 定位流体、接触或耦合时间推进 |

复用现有开关/入口；确无入口时只增加本case局部诊断选项。不再创建重复的大型test。S0/S1可以短测；S2先0.002s，过后0.02s，再0.05s。每阶段保持输出与时间步审计。

正常驱动应限制轴/轮毂等指定区域，桨叶其余区域参与弹性求解。若全部固体粒子每步被规定坐标覆盖，必须标明不是有效弹性响应实现，先修约束范围。

### 5.4 步长与执行设备

记录每个流体步实际dt_f、固体稳定dt_s、子步数、子步和、末子步余量、实际物理时间。步长由当前状态和既有稳定准则约束。

上一轮冻结solid_dt0可保留为明确对照，但正式路径不能在变形恶化后无限忽略更小的稳定步长。采用`min(启动参考上限, 当前有效稳定上限, 剩余时间)`或与现有求解器一致的合理实现；非法/过小步长应停机并定位，不能继续循环假推进。

允许减小流体步长/耦合步长、增大固体子步数，优先准确性。不能靠降低不锈钢E、增大黏度或全域清零速度获取PASS。GPU优先；若GPU固体路径仍阻塞，可用已有可信CPU固体积分配合明确的数据同步，准确报告混合执行及耗时，不假称全GPU。

### 5.5 S2通过标准（工程短跑，非材料/FSI论文验证）

- 实际时间达到0.05s，所有场有效，无时间步塌缩、超时或丢粒子；0.002s通过不能直接冒充正式短窗通过。
- 使用真实甲方几何、不锈钢弹性材料、10rpm平滑驱动。角度按实际转速时间积分，不能用未完成加速阶段的时间乘恒定ω冒充真实角度。
- 报告最大总位移，以及扣除驱动刚体转动后的弹性变形；应力、形变梯度Jacobian有效。
- 速度最大值与ωr、惯性启动、流固相互作用作数量级检查；异常尖峰必须定位，不能只满足“有限”。
- 0.05s只是短窗，不是一整圈。10rpm完成一圈的匀速周期为6s，含启动过程另算；整圈不阻塞首轮软件集成。
- 输出至少初态、中间、末态3帧真实粒子/场结果和CSV；重启若已有接口则验证，没有则明确标不支持，不临时承诺。

S2失败时按S0/S1/S2定位结果只修相应一层；至多比较原步长与减半两档，不做大范围参数扫描。

## 6. P2：甲方NC——补齐真实加热短跑，不把zero_q烟测当交付完成

### 6.1 锁定已有甲方条件

复用甲方几何和effective预设；Doc240基线为R=0.1415m、H=0.18m、dp=0.01m、4匝、100kHz、T0=1273.15K。若本机最新甲方确认数据不同，以可追溯版本为准，说明差异，禁止被法国500mm/7匝自动覆盖。

σ=16等借用值、cp/k来源、未知电流/电压与电源功率必须在配置显示。30kW若为电源输入，不能自动解释为玻璃吸收功率；缺失效率时先使用明确命名的“假设玻璃吸收功率”演示配置，且不宣称实验工况已定。

### 6.2 最小真实流程

1. validate读取真实几何、材料和运行模式。
2. prepare生成或检查粒子及线圈资产；如有缓存，核对hash并报告复用。
3. 在当前玻璃几何/σ上真实执行EM，输出原始吸收功率、功率控制方式、电流幅值约定及归一化因子。
4. 将实际Q映射到热模块，真实推进热量和流体；先0.2s再到1s，无需等待升到实验稳定温度。
5. 输出玻璃温度、速度、热源、真实线圈几何、摘要与进度。

阶段证据分别记录`em_executed`、`thermal_advanced`、`flow_advanced`、`solid_advanced`；常温热场一步未进不能写thermal PASS。

### 6.3 短测判定

- 输入/输出Q功率采用同一时刻、同一体积积分与单位，原始EM功率和控制后功率分别报告；固定吸收功率归一只能验证控制链，不能证明原始场准确。
- 热收支输出储能变化、入热、侧/底/顶对流和辐射、能量残差；同一时间区间积分。短时净冷却可以合理，不能强求温度一定上升。
- 暂用工程筛查：功率映射相对差<=1%，累计能量残差/累计能量收支绝对值之和<=5%；近零分母时报告绝对值及不评估原因。这是本轮接线标准，不是合同10%物理验证。
- 不再将开放FS背景直接称自然对流。保留原客户模式，并做同温场β=0对照；若背景明显淹没信号，采用明确标注的固定平坦滑移顶面工程预设供演示，说明与自由面模型的差异。
- 两种顶面可共存；不把固定域设置全局覆盖客户自由液面case。显著体积/粒子逃逸或异常背景仍须标未解决。

## 7. P3：法国历史F2回归——只回答能力是否仍在

### 7.1 静态核查和运行约束

查当前环带函数、冻结T赋值时序、β/T_ref、力累加与重置、pressure/background组合、C3d顶面实际算子。`paper-flat-slip-top`是组合配置，不能只按名称认为与旧版相同。

先用当前代码跑历史配置，不整体退回旧求解器；使用同一reload及材质/常物性。记录历史缺失项。历史C3d：声学不启用AcousticWallFreeSlip，黏性OpenCup保留，p′+静水背景，不重复完整GravityForce。

必须先validate，排除其他温度初始化分支抢占环带分支；输出实际温度模式和热点坐标。

### 7.2 最小对照

- F2_OFF：相同空间冻结环带T，β=0，Q=0。
- F2_ON：完全相同配置，β=1e-5。先10s，稳定推进后20s；仅为直接比较第224轮后期结果且仍有必要时到30s。

两组温场保持完全相同，不能OFF均温、ON环带热，也不能只一组解热扩散。若OFF已有足以污染目标的流动先停止延长，定位背景，不为凑图继续跑。

输出瞬时与末完整时间窗平均速度；若确实计算到20/30s，再使用与历史相同的5s窗口，不用全程累计平均替代后期窗口。

### 7.3 判断方式

- 对比相同r-z位置的Ur/Uz、轴/环/缘方向和U_rms数量级，检查热点位置与实际温差。
- 同时画原始速度矢量与插值流线；不以流函数符号区域数代替涡数。
- β=0背景U_rms低于ON的10%可视为信号清晰的工程筛查；背景更大标污染，不能从ON减OFF的范数声称得到物理流场。
- 观察到两个回流区域可写`DIAGNOSTIC_TWO_CELL_OBSERVED`，稳定性未过则NONSTEADY；不要强行扩大温差维持两个涡。
- 本测试成功只证明指定温场下的水动力能力，不等于EM—热耦合复现。

若失败，先从温场赋值、顶面、压力、黏性与步长的有效差异中找首个不一致；一次仅变一项。不得新增shifting/δ-SPH/ghost路线来代替归因。

## 8. P4：Fig.8温度场与Fig.7速度场联合对比

### 8.1 先验收参考数据，不重复虚构数字化

查现有Fig.7/8提取数据、图像、坐标标定和239资料。Fig.7速度矢量/流线与Fig.8等温线不是同一种数据，文件不能混用。

若现有数据不完整，本轮补齐：

- 保存原始图页/裁剪和带标定点的覆盖图。
- 记录r/z轴、单位、轴线/壁面/自由面坐标、温度无量纲定义及图例。
- 逐条提取可辨识等温线的标签与点坐标；看不清的标uncertain，不编造数值。
- 无量纲定义必须核对原文；无法确认则只做归一化形状比较，不用默认T/Tmax换算绝对温度。
- 分清数字化误差与文献缺失参数；表中列坐标误差、温度标签误差和数据来源。
- 缺少图像时只能进行文字趋势检查，不能宣称已完成定量图像对齐。

不要将人为补全的二维温场当“论文完整数值解”。稀疏等温线插值可作有不确定性的初猜，不能当高精度全域真值。

### 8.2 每组法国计算必须输出同一时间窗的三联图

使用现有Python管线完善，不要求手动ParaView作为唯一依据：

1. 温度剖面：有单位绝对T、论文定义的无量纲T分别画，叠加参考等温线；所有对比用相同坐标与色标。
2. 速度剖面：Ur/Uz矢量、真实速度幅值色标、插值流线；轴/环/缘分区及数据不足区域明确显示。
3. 热源剖面：原始Q和控制后Q、累计功率径向分布，标热点与高热源区域。

对3D场同时保留真实中心薄片与周向平均r-z场；轴对称参考优先对周向平均，真实切片用于检查非轴对称。不得将整直径镜像后的左右大回流算作两个径向胞。

r-z轴对称流函数若使用，应满足`ur=-(1/r)*∂ψ/∂z`、`uz=(1/r)*∂ψ/∂r`的约定，明确轴线处理。回算速度误差>20%时不用ψ单独证明闭合；不要只修ψ图而忽略原始速度。

### 8.3 比较量与本轮判定边界

| 类别 | 必须报告 |
|---|---|
| 温度 | Tmin/Tmax/体积平均T；热点(r/R,z/H)；中心/中环/外缘竖向剖面；上/中/下水平剖面 |
| 等温线 | 同标签轮廓位置、向轴心弯曲方向、热区面积及匹配误差；参考不可辨识区剔除并标注 |
| 速度 | 轴/内部环带/外缘Uz、U_rms/Umax、主要回流旋向；小涡是否仅单bin宽 |
| 热源 | 体积积分功率、高热源位置和范围、径向/轴向功率占比 |
| 收支 | 热输入与各面散热、储能变化、残差、运行物理时间与统计窗口 |

轮廓位置误差可用H或R归一的距离并列图像提取误差；在拿到可靠参考后再填写数值。未预先约定的指标不临时设置“刚好通过”的阈值。

本轮只区分：参考数据可用/不足、方向一致/不一致、温场接近/明显偏离、准稳态/未达准稳态、论文验证未完成。不得把自定形状误差当合同10%通过。

## 9. P5：法国耦合推进——修温度场的主因，不盲目追两个涡

### 9.1 先完成真实热更新审计

沿调用链确认k(T)、cp(T)、μ(T)、σ(T)在哪一步读取，作用于哪个算子。重点检查`ApplyLaplaceDiffusionDtCK`使用标量cp还是逐粒子cp，接触扩散是否用局部导热率，以及边界热量到粒子温升的换算。

用少量覆盖物性曲线变化区的温度点做函数/单步能量检查；不把不合适的大温差制造到正式论文工况中。热物性跨平台测试属于局部测试，不替代论文运行。

若局部cp未实际生效，优先在OPHELIE热模块局部修复。温变cp储能应采用与数值更新一致的焓差积分，不直接用某时刻cp乘全温差充当精确储能。侧/底/顶热流的面积、法向、单位和重复计数逐项核对。

### 9.2 首轮只选一个明确的法国主配置

选择`FRENCH_PAPER_NC`作为主比较；从期刊锁定已知几何、功率、热边界、顶面机械条件与β/ρ等，未明物性使用239提取的Jac资料并标BORROWED。博士55kW用于辅助，不能替换论文身份。

σ公式A/B先使用同一非均匀冻结T，同一几何和功率约定对比真实EM热源；不让两组热演化不同导致无法归因。上一轮thaw A/B只保留为耦合示例。

A/B不要求两条都跑长时间。根据原始功率、固定功率下Q形状与Fig.6位置趋势选一个主配置，另一个列敏感性候选；不要按温度一个点最接近就宣称公式已证实。

EM更新间隔是本实现数值参数，论文未明说则标ASSUMED；用相同状态检查更新前后σ与Q变化，必要时仅做一个减半间隔短测。不开自感时明确A_ind=off，固定功率归一不能掩盖此模型差异。

### 9.3 逐级热流动耦合，避免起步状态误导

1. 合理初猜可用当前热模块或有边界的轻量轴对称热计算；这是加速收敛，不是预先指定最终速度。
2. 保留固定功率及目标散热边界；先看热收支和热点是否向参考区域发展，再决定延长时间。
3. 短窗后画联合三联图：明显中心热核与参考上部环带热区不符时，先查Q/散热/μ冷壁区，不能只改浮力。
4. T已接近、方向仍错时，回看β/T_ref、浮力符号、局部μ和顶面机械条件。
5. 方向一致、小外围涡不清时，先检查热区位置和冷壁/桨等细节分辨率，再考虑一档加密；不先全域扫描分辨率。
6. 全场仍在变但趋势合理时，才安排延长；不得把1s结果与论文稳态直接判失败，也不得仅靠跑长时间期待错误边界自动修复。

工程准稳态筛查：三个连续等长完整窗口中，平均T变化小于代表温差的1%/窗，主要速度统计变化<5%/窗，储能速率相对输入/散热功率规模<5%，无异常漂移。低信号用绝对量说明。以上是本项目筛查标准，非论文原始标准，未达到就标TRANSIENT。

### 9.4 参数调整允许范围

用户允许对未知条件校准，但只能调整已登记的不确定参数。先给一张参数表：已知固定项、博士借用项、未知待校准项、范围依据、预期影响。

优先校准散热与黏度曲线来源，再判断电磁热源形状；每次只动一类参数。所有校准结果明确标CALIBRATED。禁止改动已知物理单位、强写速度、任意增β、硬设热点来冒充预测结果。

本轮最多两组耦合候选：主配置＋针对已识别最大偏差的一组修复。仍无改善时提交调用链证据和对比，不继续无边界扫参。

## 10. P6：Linux软件交接，必须能在新目录真实启动

### 10.1 软件最小接入契约

复用既有validate/prepare/run、case.json/effective JSON和进度接口。没有统一二进制时用薄包装调用已有target，不要求先完成大型框架迁移。

每个case交付：

- 可构建提交或明确base+完整patch，必要CMake与依赖说明。
- 几何、坐标与单位说明、材料表、reload及其manifest；大资产不必进git，但要有校验和与获取位置。
- 默认项目JSON、生成的effective JSON、完整可复制命令、所需环境变量。
- 初态/过程/末态结果、进度、summary、日志、CSV、字段单位与关联文件清单。
- 模式和能力标签：真实EM/热/流/弹性、GPU或CPU/混合、短窗/稳态、校准/借用参数、验证未完成项。

求解器端应输出真实状态。软件不得根据按钮名称推断“电磁已运行”或把回放记为新计算。

### 10.2 清洁目录复跑门控

使用新的运行目录、相对资产路径，在不依赖开发目录隐含文件的情况下执行两个预设。依次完成validate→prepare或有效缓存复用→真实短跑→加载结果。失败时报告准确缺件，不能以历史结果覆盖输出目录。

甲方NC：真实EM、Q映射和热/流推进达到1s；甲方弹性：真实S2达到0.05s。两组均通过才授`LINUX_SOLVER_HANDOFF_READY`。

该标签表示软件可调用的真实短跑基线，不代表法国文献、甲方1000℃测点或合同模块验证已完成。后续若扩展更长时间仍须独立检查。

若未通过，交接状态仍写BLOCKED并列具体原因；GUI界面开发可继续，但不能把按钮禁用/回放当本轮物理任务完成。

### 10.3 两机分支与版本

物理机与软件机使用各自分支和工作目录。物理改动形成可审查提交后由软件分支显式cherry-pick/合并；git push不会自动改另一台机的工作区，pull/merge才引入变化。不要让两机长期同时向同名共享开发分支混推未协调改动。

交接写明commit及全部必需资产hash；WIP补丁必须记录基线、应用顺序和hash。禁止提交整个临时结果仓库、强推或清理用户未提交文件。

## 11. 执行顺序、停止条件与本轮边界

| 顺序 | 工作 | 完成/停止条件 |
|---|---|---|
| 1 | P0 + reload审计 | 诊断有效、尺度与失败证据清楚 |
| 2 | S0/S1/S2 | 优先解除甲方弹性阻塞；失败停在首个异常算子 |
| 3 | 甲方NC真实EM—热短跑 | 两案形成真实交接候选 |
| 4 | 清洁目录复跑并交接 | 软件机可直接接入；不等法国稳态 |
| 5 | F2历史回归＋参考数据检查 | 确认指定温场能力和图像目标 |
| 6 | 法国三联图＋一次主因修复 | 温度/流动/热源可联合解释 |

CPU整理、源代码审计和参考标定可穿插，不在单GPU上同时启动几个正式求解器。每次失败保留第一个坏步及前一有效状态；不要只上传最后一屏日志。

本轮不重开TEAM7、不重跑已充分验证的标准方腔、不新增ghost粒子、不重写Riemann核心、不追求整圈搅拌、不做高分辨率长时全耦合盲跑。若发现本轮改动会影响已有方腔/EM核心，则仅补对应必要回归，并说明触发原因。

## 12. 必须交付的文件和结论格式

建议根目录`doc241_runs/`，包括：

```text
00_freeze/                 # 代码、二进制、设备、资产hash
01_client_reload_audit/    # 粗/密reload对比与尺度修复
02_client_elastic/         # S0 S1 S2、状态场、逐步dt与失败首帧
03_client_nc/              # 真实EM/热/流、能量与功率映射
04_f2_regression/          # OFF ON、历史对比、原始速度与温度
05_paper_reference/        # Fig6/7/8来源、标定、数字化与不确定性
06_french_joint_compare/   # T/U/Q三联图、剖面、热收支、有效配置
07_linux_handoff/          # 两案预设、资产清单、清洁目录复跑
RUN_COMMANDS.md
ROUND241_DECISION.md
WORK_SUMMARY_ZH.md
NEXT_TASKS.md
```

每一组至少包含effective配置、summary、日志、原始CSV及用于展示的真实帧；Python绘图脚本必须纳入可追溯源码。已有原始数据本机完整保留，精简包提供清单和关键证据即可，不要求把全部历史重新打包。

summary至少区分：

```text
workflow_completed
reached_end_time
requested_end_time_s / actual_end_time_s
em_executed / thermal_advanced / flow_advanced / solid_advanced
finite_fields_ok / geometry_resolution_consistent
short_run_passed / steady_state_reached
paper_comparison_status / experimental_validation_status
execution_backend / source_commit / binary_hash / asset_manifest_hash
fail_reason
```

若某字段尚未实现，不手工填PASS；补实现或标NOT_IMPLEMENTED。

最终总结直接回答七个问题：

1. 弹性失败由什么证据定位，reload/dp不一致是否成立？修了哪一层？
2. 当前甲方弹性是否真实到0.05s，采用什么粒子集、材质与设备？
3. 甲方NC是否真实计算了EM并推进热量，能否从清洁目录启动？
4. 当前代码能否复现历史人工环带双回流；不能时首个差异是什么？
5. Fig.8温度场与当前场最大偏差在哪里，与Fig.6热源是否对应？
6. Fig.7方向不一致究竟是温场驱动不同，还是在相同温场下水动力仍不同？证据不足则写未定。
7. Linux软件机现在可以启用哪些真实计算功能，哪些仍未通过物理/合同验证？

不要再只交“补了几个开关、短烟通过”的总结。本轮应交付可运行的甲方基线，以及能定位法国温度场差异的一套联合证据。

<!-- ===== END 241_OPHELIE_FRENCH_TEMPERATURE_FLOW_AND_CLIENT_ELASTIC_RELEASE_CURSOR_PLAN_2026-11-24.md ===== -->


<!-- ===== BEGIN 242_OPHELIE_FINAL_PHYSICS_PREPARATION_AND_LINUX_ENTRY_CURSOR_PLAN_2026-11-24.md ===== -->

# 242_OPHELIE_FINAL_PHYSICS_PREPARATION_AND_LINUX_ENTRY_CURSOR_PLAN_2026-11-24.md

# 242：Linux软件开发前最后一轮物理整改、场图对比与调用闭环

日期：2026-11-24
执行对象：当前物理开发机上的Cursor。
前置资料：Doc241完整交接包、239物性数据包与239A、240/241任务书、现有Linux软件源码。

## 0. 用户决定与本轮结束方式

用户决定：再做一轮物理与调用准备，随后转入Linux软件开发。**本轮不再扩大研究路线，不以论文稳态完全复现或合同验证全部完成作为软件开工条件。**

本轮要交付四项：

1. 甲方弹性搅拌异常速度/密度的具体定位与针对性修复，保留真实不锈钢弹性计算。
2. 甲方自然对流的低背景工程预设与真实EM—热短跑。
3. 法国F2实际速度图，以及至少一组期刊目标工况的温度—速度—热源对比，结束只报能力标签、不交场图的状态。
4. 两个甲方算例通过现有CLI的validate→prepare→run→结果索引闭环，而非仅手工运行C++ binary。

完成本轮后冻结交接版本，转Linux桌面软件。未完成的物理验证保留明确问题单继续研究，不能伪装PASS，也不能让任务无限追加，拖住软件界面开发。

本文件接替241的下一轮执行安排；未被修改的物性来源、案例身份隔离规则继续有效。独立A-phi/Laplace分支不改。不开发Windows，不新增GUI框架，不要求本轮安装Qt/VTK或完成桌面发行包。

## 1. Doc241事实：必须沿用，不能重新解释成全部成功或全部失败

| 项目 | 实际证据 | 本轮处理 |
|---|---|---|
| 粗reload与参考dp | 粗Vol约1e-6 m³，对应10mm；旧入口强制5mm | 不一致已确认；不重复审计整套历史 |
| 当前弹性桨 | 重生5mm，4329桨粒子、43928玻璃粒子，GPU真实跑0.05s | 保留可运行基线 |
| 弹性流体异常 | 末Umax=1.381m/s；CSV峰值3.511m/s；桨最大速度约0.071m/s | 必须定位，不能以有限值为通过标准 |
| 弹性密度 | rho_ref=2700，末rho_max约3641kg/m³ | 约+35%；需查位置/分布/初始即有或推进产生 |
| 固体步长 | 仍固定启动solid_dt0 | 补当前稳定步长与完整时间覆盖检查 |
| 甲方EM—热 | 假设玻璃吸收15kW，P_joule=15000，P_heat=14991.8W，1s | 接线已通，不必重新研发EM |
| 甲方浮力ON/OFF | U_rms约0.00341157/0.00341174m/s | 原开放FS短窗中背景主导；不宣称自然对流已验证 |
| 法国F2 | OFF 10s零流；ON 10s末7.30e-4、20s末5.19e-4m/s | 清晰浮力响应；尚无本轮双胞场图证据 |
| Fig.7/8 | Doc204数字化已在241包中 | 直接使用并标不确定性，不重新声称没有数据 |
| 论文主线 | 241 P5暂缓，未做温度场叠图 | 本轮必须给实际图，不再因Fig.6缺数字化而全推迟 |
| 软件层 | NC默认heating=none，run.py固定zero-q；搅拌prepare未传dp | 必须同步，底层PASS不等于软件已调用成功 |
| 清洁复跑 | 241是新输出目录直接运行binary，且关闭场输出 | 本轮验证统一CLI与可加载场文件 |

重要归因边界：Doc241同时纠正尺度并加密粒子，不能断言“此前所有失败唯一由dp错误造成”。无需为证明唯一根因继续穷举粗粒子。

## 2. 本轮执行顺序与资源限制

按依赖执行：

1. A：保留241版本、补关键诊断与错误判定。
2. B：甲方弹性最小对照，抓首个异常并只修对应路径。
3. C：甲方NC固定平坦滑移顶面与低背景对照，接回真实EM—热。
4. D：法国F2绘图与Fig.8温度对比；图像整理可与前面穿插。
5. E：现有CLI参数传递与结果索引；调用层修改可穿插，正式闭环测试在物理预设冻结后执行。
6. F：冻结版本、交接Linux软件开发任务。

单GPU正式运行串行；不要同时跑多个solver造成显存/耗时异常。尽量复用241已有结果，缺空间数据才重跑。默认只使用当前分辨率，不开始网格扫描。

诊断矩阵最多：弹性3组基础对照＋1组针对性修复（必要时再加一档dt/2）；NC固定域OFF/ON＋1组真实热耦合；法国F2有数据则零重跑、无数据最多1组ON；期刊目标1组初始联合对比＋1组明确主因修复。超出矩阵必须在总结中说明新证据与必要性，不能自动转为扫参。

## 3. A：一次性冻结与关键诊断修复

### 3.1 冻结已有成功版本

保存实际HEAD、WIP diff、binary hash、设备、构建命令、STL/reload/材料hash。Doc241记录base为16d67c9d0799aeee899a914059fc51d67740f738且另有WIP；不要只检出base就假定含全部修复。

保留原5mm密reload，不覆盖；新实验另建目录。物理修复只限OPHELIE局部适配、case和必要调用层，避免影响SPHinXsys其他case与A-phi研究。

### 3.2 修正会使坏运行变成PASS的源码逻辑

在`test_3d_ophelie_client_stirring_f0_elastic_fsi.cpp`重点处理：

1. 现有`time_ok = t_final >= 0.8 * end_time`不能用于`reached_end_time`。以实际达到目标为准；允许最后一步的小量超调，报告绝对/相对差，不允许少跑20%还标完成。
2. 现有`coupled_ok`把“有转速”与“有流固载荷”用OR相连。分开记录drive_executed、solid_advanced、fluid_contact_evaluated、load_sample_valid；驱动本身不能证明FSI载荷传递。
3. 固体子循环到kSolidSubMax后当前仅WARN并返回，可能尚未覆盖dt_fluid。时间未覆盖必须FAIL/停止或按严格规则继续补完，不得流体继续前进而固体欠时间。
4. 检查所有向量分量和张量分量；hostRotorSnap只检查pos[0]/force[0]不足以审计固体。新增固体Position、Velocity、变形梯度/应力（若有）、Jacobian的有效性；场名以实际代码为准。
5. 区分finite_fields_ok、numerical_quality_ok、workflow_completed、validation_passed。不存在的验证不能用passed=1代替。

不要新增庞大测试框架。用一个缩短运行、一个非法配置和已有坏摘要确认判定有效，再开展正式矩阵。

### 3.3 输出最小诊断

逐流体声学步记录time、dt_f、dt_s_candidate、dt_s_used、solid_covered_time、solid_substeps、Umax、U_rms、密度分位数/极值、max弹性偏差、力/转矩；异常附近频率加密。

对速度与密度前10个异常粒子记录稳定ParticleID、位置、到桨/壁/顶面的距离或接触分区、邻域数、压力。若排序改变索引，用稳定ID追踪，不能把排序后的同一数组下标当同一粒子。

## 4. B：甲方弹性搅拌异常的最小归因与修复

### 4.1 固定条件

采用Doc241可运行5mm dense reload、真实甲方桨与容器、不锈钢E约2e11Pa、rho_s约7900、玻璃rho约2700、当前黏度14.532Pa·s、10rpm与0.05s加速。驱动轴、轴心、坐标单位按当前有效配置打印并核对STL变换，不凭名称假设Z/Y。

保留顶端5cm驱动带但标provisional，画出278个约束粒子分布，确认没有覆盖整桨或桨叶。不能把所有桨粒子强写规定轨迹冒充弹性响应。

本轮不再试483粒子粗桨，不降低E凑稳定，不通过增加黏度或清零流场压制异常。

### 4.2 三组基础对照

| 编号 | 条件 | 时长/输出 | 要回答的问题 |
|---|---|---|---|
| B0 | 241原配置复跑，开启空间场和高频诊断 | 先0.02s，捕获已知约0.01s速度尖峰；必要时到0.05s | 异常从哪些粒子/哪一步出现？ |
| B1 | 同一粒子/物性/流固接触，桨驱动转速=0 | 同B0窗口 | 不旋转是否也有尖峰和压缩？ |
| B2 | 同一弹性桨和10rpm驱动，关闭流体载荷对固体的传递，作为固体诊断 | 0.02s或覆盖首个异常时刻 | 固体内部/驱动是否独立正常？ |

B2是诊断，不能授FSI验证；优先复用既有开关，没有则添加case局部开关。若B0尚未出现异常，先核对是否真正用了同一二进制/粒子，不随意加长到数秒。

输出初态、首个尖峰前后、末态VTP；全程只存5–8帧即可。额外保存速度/密度直方图和异常位置切片。

### 4.3 根据证据选择唯一优先修复

- B1也明显异常：优先检查粒子初始重叠/间隙、初始压力/密度、重力与壁面平衡、流固界面初始化；先查异常位置再决定是否需要短静置。
- B1正常、B2固体异常：检查驱动速度/位置时刻一致性、约束范围、参考构型、固体步长与刚体旋转误差。
- B1/B2正常、B0异常：检查流固压力/黏性载荷方向、界面速度时间层、AverageVelocity的实际定义与赋值、移动后的cell/邻域更新。沿已有成熟FSI路径比较，不按猜测新增经验阻尼。
- 异常只在自由面或容器壁：定位边界贡献与背景，不能直接归咎弹性E或电磁。

只实施证据指向的修复；保留修前修后同一诊断。若多个问题并存，先修首个异常，不把全部模块同时替换。

### 4.4 固体子步修复

恢复当前状态的稳定步长评估，打印其输入；可采用`min(启动参考上限, 当前稳定候选, 剩余流体步时间)`的保守方式，或与本项目已验证固体求解器一致的实现。

候选步长非正、非有限或持续塌缩时停机保存首个坏状态，不能回退到更大的固定dt强行走完。检查子步累计时间与流体dt一致，记录所有截断/子步上限。

若GPU候选规约有问题，用明确的小型CPU复核或可信CPU固体路径定位；准确性优先，混合执行需标明。无需为本轮重新实现GPU隐式固体。

### 4.5 B组判定与延长条件

必须：实际达到0.05s；无丢粒子、非法场、非法Jacobian、子步欠时或超时；输出真实弹性变形与载荷，所有异常定位有结论。

质量筛查同时报告：

- Umax、U_rms、U95/U99及对应粒子比例，不能只看一个极值也不能丢弃极值。
- 速度与桨速之比、U/c0、峰值持续时间和位置；若继续出现1–3.5m/s量级速度而无合理解释，则numerical_quality不通过，不称低速搅拌流场可信。
- 密度相对参考值及相对合理静水/静置基线的偏差。WCSPH静水密度可随压力变化，不能将全部偏差认定为错误；但约35%压缩不能无解释接受。
- 暂用工程目标：内部流体相对静置平衡基线的密度偏差P99<=2%、最大<=5%，超出标WARN/FAIL并定位；这是本轮筛查阈值，不是论文或合同标准。边界层单列但不从全场统计中隐藏。
- 动能、外驱动功、重力/压力释放等来源；如高速度由初始化产生，需证明静置后衰减且之后驱动信号可分辨，不仅解释一句“启动瞬态”。

仅当0.05s质量合理且耗时可接受，允许一次延长到0.1–0.2s检查趋势；不要求整圈6s，本轮不靠长跑掩盖初始异常。

## 5. C：甲方NC低背景工程预设与真实电磁加热

### 5.1 保留原开放FS，新增明确的固定域预设

Doc241开关浮力结果几乎相同，说明该短窗不能识别浮力贡献，不说明实际玻璃不存在自然对流。

保留`CLIENT_NC_OPEN_FS_DIAGNOSTIC`，新增`CLIENT_NC_FLAT_SLIP_ENGINEERING`。复用法国固定平坦不可穿透、切向自由的已存在局部配置，核对组合开关和实际压力/重力路径。

侧/底无滑移；顶部机械条件与散热分开，不取消辐射/对流。程序中明确hydro_model及top_boundary，软件以后展示这一建模假设，不把固定顶面说成真实可变形自由面。

### 5.2 一对严格可比较的固定温场诊断

C_OFF/C_ON采用甲方相同几何、相同粒子、相同空间冻结温场、相同物性与顶面，只有β不同；Q=0。可复用ΔT=20K径向诊断场，明确其为人工诊断，不是甲方实验温场。

先0.2s检查，再到1s；报告U_rms、各分区Uz、最大速度与背景。OFF/ON的温度赋值规则必须一致，不能OFF均温、ON有温差。

工程信号目标：OFF U_rms<=ON的10%，ON方向符合指定温场驱动；若ON还过弱，先报实际绝对量/预期加速度和时长，不增β硬凑。固定域仍异常时查本机配置与法国零背景配置差异，只修首个不一致。

### 5.3 恢复真实EM—热—流短跑

使用C组已确认的工程预设，真实EM、σ=16借用值、T0=1273.15K、当前甲方几何和线圈参数，先推进1s。

继续明确：15kW是ASSUMED玻璃吸收功率；30kW是电源输入，不能默认效率=50%。当前电流/电压未确认，定电流模式不可用；固定吸收功率模式可根据目标求解等效电流，并标该电流为计算结果。

输出原始EM功率、控制后功率、实际入热、σ场、温度与速度场，以及每个散热边界累计能量。初态散热大于15kW而净冷却可以合理，不能要求短时一定升温。

复核同时间区间收支：功率映射差工程目标<=1%，累计能量残差按241统一定义<=5%；近零分母用绝对量，不能用随意大floor压低误差。固定功率吻合仅证明控制链，不证明EM空间精度。

C组通过后，这个明确标注假设的固定域EM预设作为Linux首轮默认候选。原开放FS仍保留研究入口，禁止无提示更换建模条件。

## 6. D：法国必须交实际场图，不能再次只有参考盘点

### 6.1 F2先利用已有20s结果

检查原机输出是否存在F2_ON_20s的r-z速度/温度CSV或VTP；有则直接绘图，不重跑。若已被后续运行覆盖且不可恢复，使用242当前代码重跑一组ON至20s并开启输出；OFF已有清洁结果可复用，除非水动力修改影响该路径。

输出：原始Ur/Uz矢量、带速度幅值的流线、温度等高线；同时展示瞬时末态及一个完整后期时间窗的平均。不要把全程平均当后期稳态。

实际回答：轴心是否下降、内部环带是否上升、最外缘是否下降、是否有两个回流区。没有可靠闭合证据就写direction_observed；双胞不明显也要交图，不继续盲目延长。

历史200K冻结环带是诊断温场，不作为Fig.8温度复现。左右镜像的两个大回流不算径向主次两个胞。

### 6.2 期刊目标工况至少跑一组，并完成温度对比

目标为`FRENCH_PAPER_NC`：期刊50kW、无热电偶装置、平坦滑移顶面，已知参数来自期刊，未明项按239提取博士资料补充并标BORROWED。禁止用Jac55kW含Marangoni结果替代期刊Fig.7/8。

先给出effective身份表：几何、ρ、β、μ/k/cp/σ曲线、频率与线圈、功率控制、侧/底/顶散热、辐射参数、是否Marangoni、自感开关、初始温场、更新间隔。缺失项标ASSUMED，不擅自称全部与论文相同。

温变物性审计只补241尚未证明的关键链：局部cp/k是否进入真实热更新；不能只打印物性数组。若仍标量cp，先明确当前近似及误差风险，优先修本case局部路径；温变cp能量按一致焓差检查。

初猜可来自现有热模块/轻量轴对称热计算，或标明不确定性的参考温场插值。初猜不得被持续冻结成论文结果，正式联合工况必须允许温度演化。

先短窗输出三联图与热收支，检查趋势后再决定一个有限延长窗口。不要把1s瞬态与论文稳态直接作最终误差判定。只有初始/短窗结果时标TRANSIENT，仍需完成叠图说明差异。

### 6.3 参考与图像要求

直接使用241包内`reference_fig7_fig8_Doc204/`：保留原图、标定、等温线标签及不确定性。Fig.8是穿轴线的竖直温度剖面，Fig.10是带热电偶三维计算，不能混用。

核对无量纲温度定义；无法确认时只对可确定的形状/坐标比较，不擅自换算绝对T。Fig.6未数字化不阻塞Fig.8叠图，可以先用Fig.6原图做位置定性比较，标明无定量误差。

每组正式法国结果交：

1. T：绝对温度与论文定义无量纲温度、参考等温线叠加、同坐标同色标。
2. U：周向平均r-z矢量/流线及真实中心薄片，显示非轴对称差异与缺采样区域。
3. Q：原始/控制后焦耳热源与热点位置、径向/轴向功率占比。

至少报告Tmin/Tmax/Tmean、热点位置、中心/中环/近壁竖向剖面、上/中/下横向剖面；轴/环/缘Uz与U_rms；输入、散热、储能、物理时间。

### 6.4 只允许一次有依据的主修复

先根据联合图归类：

- 热区位置错：查Q形状、σ分布、功率映射和边界散热。
- 冷壁层/温度梯度错：查局部k/cp、热边界面积与μ(T)冷区响应。
- T趋势接近而速度方向错：查浮力/压力/黏性及边界路径。
- 主流向对、小外围胞不清：先判当前分辨率与采样是否足以描述，不为凑双胞任意调β。

选择证据最明确的一项改动，给出修前修后同图。未知h/μ曲线等可校准，但预先列范围与来源，并标CALIBRATED；不改已知单位或强写速度。

不要求本轮达论文稳态PASS；要求第一次形成可解释的真实差距图和下一主因，不允许再次以“参考已盘点，比较留后续”结束。

## 7. E：完成软件开工所需的CLI调用链，不在本轮开发桌面

### 7.1 只改已有Python适配层，不另起第二套接口

路径为仓库内`docs/ophelie/software/linux_src/`，241包中对应`software/linux_src/`。以本机实际实现为准，以下问题已在241快照看到：

- `coldcrucible_solver/run.py`：NC固定加入zero-q，仅支持isothermal_hydro/dT20_boussinesq。
- `prepare.py`：搅拌生成仅传`--bodies=all`，dp未传到生成器。
- `examples/natural-convection/case.json`：heating=none，与真实15kW短跑不一致。
- `run.py`：搅拌omega-ramp、clamp等硬编码；优先读取case并写effective，不让GUI参数看似可调却不生效。

不要求统一为新C++ binary；现有薄包装可继续。

### 7.2 新增明确的NC真实加热模式

保留旧诊断模式，新模式明确fixed_absorbed_power、target_absorbed_power_W=15000、power_status=ASSUMED、generator_input_power_W=30000、current_status=unknown、固定域顶面身份及1s结束时间。

字段名可沿现有schema统一，但JSON、validate、命令构建、effective、summary必须一致。真实EM模式不允许混入zero-q。未知电流不阻塞固定吸收功率模式，但必须阻止错误的given_current调用。

有冲突参数validate必须报错，不静默降级等温模式。非法mode不允许fallback成功。

### 7.3 prepare与资产一致性

搅拌prepare必须实际向生成器传dp=0.005；若生成器暂不接受该参数，先补该case的参数入口。读取生成后的体积/坐标近邻尺度和manifest验证，不只是写prepared.json中的标签。

reload检查分玻璃、桨、壁独立统计。Vol^(1/3)是均匀体积粒子的参考，不能独自代替空间最近邻、局部h与拓扑检查。缓存fingerprint加入几何文件内容hash与影响生成的配置，避免替换同名STL仍复用旧粒子。

修复相对资产路径/软链，提供可复制资产清单。运行时不得静默优先使用binary/data里旧桨而忽略项目所选桨；若当前只支持固定几何，validate明确限制并核对hash，不能假称任意STL已支持。

### 7.4 结果与进度

正式演示预设必须开启初态/中间/末态至少3帧VTP，并写结果索引（名称可复用现有机制）：文件相对路径、body、时间、可用字段与单位、帧数、实际结束时间。

建议玻璃Velocity/Pressure/Temperature/Q；桨Position/Velocity/弹性位移/应力或当前可用弹性字段。尚未实现的字段不伪造，列表注明缺失。空间冻结诊断与正式热计算区别显示。

每次运行独立工作/输出目录；如底层仍固定写binary/output，适配层需隔离或安全收集并禁止并发混写，不能索引上一运行残留文件。

汇总CLI不得仅看returncode或passed，应带出实际能力与质量标签。进度、取消、错误退出、结果加载路径沿现有协议，避免另建协议。

### 7.5 一次真正端到端测试

在不含旧结果的新项目目录，用软件正式入口依次运行：

```text
validate → prepare（或hash验证通过的缓存） → run → 读取result index
```

两个case都必须通过该入口，不再以手工调用binary替代。输出完整命令及退出码。确认：

- NC summary显示EM实际执行、热与流实际推进；目标15kW及ASSUMED状态一致。
- 搅拌为5mm真实弹性桨、实际到0.05s，数值质量标签与B组一致。
- 至少3帧存在、可用现有VTK/Python reader读开，时间单调、点数合理、要求字段存在且有效。
- 验证一个非法功率模式、一个dp/缓存不一致负例，会清晰失败，不悄悄切换模式。
- 不要求在本轮打开Qt GUI；reader能读取与GUI完整验收分别记录。

## 8. 明确通过标准：不要再只有一个总PASS

| 标签 | 条件 |
|---|---|
| WORKFLOW_PASS | CLI全流程完成、确实到终止时间、结果可读取 |
| NUMERICAL_QUALITY_PASS | B/C对应速度、密度、时间覆盖及场有效性检查通过 |
| EM_THERMAL_WIRING_PASS | 真实EM/Q/热更新与能量映射检查通过 |
| BUOYANCY_SIGNAL_RESOLVED | 同温场β对照可分辨，不由背景主导 |
| PAPER_COMPARISON_AVAILABLE | 实际计算图与Fig.7/8参考完成对比并列差距 |
| PAPER_VALIDATED | 本轮默认不授；需之后充分稳态、参数和误差证据 |
| LINUX_INTEGRATION_BASELINE_READY | 两案CLI和场文件闭环完成，假设/能力/质量状态完整 |

`LINUX_INTEGRATION_BASELINE_READY`是可集成状态，不自动覆盖NUMERICAL_QUALITY_PASS。若某case质量仍FAIL，可以开发软件、导入/回放已知结果，但默认运行预设须清晰标诊断用途；不得称物理可信或第三方验收已通过。

当异常已定位但本轮修不完：给出首个坏状态、相关源码、复现命令、影响范围、推荐下一刀；不再开启无限新矩阵。用户按本轮结束转软件的决定执行，物理未完成项独立跟踪。

## 9. 本轮禁止的扩张与允许的回归

- 不重开TEAM7/全套电磁验证，不重跑整套方腔，不新增δ-SPH/ghost/shifting路线。
- 不在本轮要求整圈弹性搅拌、千秒加热或高分辨率全耦合长跑。
- 不以固定吸收功率吻合证明EM原始热源正确，不以人工温场双胞证明法国耦合复现。
- 不改不锈钢E为软材料、不任意增β、不清零流体速度凑稳定。
- 不将未知电流换成无来源常数、不将15kW假设写成甲方确认。
- 如果局部修改触及已有核心算子，补一个与修改直接相关的必要回归；不因此重启整个验证项目。

## 10. 本轮交付目录与最终问题

建议`doc242_runs/`：

```text
00_freeze/
01_elastic_diagnostics/        # B0/B1/B2、极值位置、初始/尖峰/末态
02_elastic_fix/                # 针对性修复及时间步/子步覆盖审计
03_client_nc_fixed_domain/    # OFF/ON和真实EM热短跑
04_french_f2_fields/           # 现有数据绘图或唯一必要重跑
05_french_paper_joint/        # Fig8叠图、T/U/Q三联图、误差/差异表
06_cli_e2e/                   # 两案validate prepare run与结果可读测试
07_linux_entry/               # 冻结提交/补丁、资产、示例、能力清单
RUN_COMMANDS.md
ROUND242_DECISION.md
WORK_SUMMARY_ZH.md
KNOWN_LIMITATIONS.md
```

物理曲线和结果文件本机完整保留，精简交接包包含关键CSV、5–8帧必要场结果或可获取资产路径/hash、脚本、有效配置和日志。不能本轮再只打包summary，却丢失空间场。

冻结交接优先形成可审查提交；不要因“用户未说git commit”自动停住本地版本整理。遵守本机仓库规则，提交只包含明确相关源码、配置、脚本和小资产；不强推、不把临时大包和全部结果纳入git。若用base+patch，给应用顺序、hash和构建步骤，不把HEAD当全部运行代码。

本轮最终总结必须直接回答：

1. 1–3.5m/s异常在哪些粒子/哪一步产生，桨静止是否也存在？有无针对性改善？
2. 约35%密度偏差的分布与来源，修后是否仍存在？
3. 固体是否覆盖完整流体时间步，当前稳定dt是否真的参与控制？
4. 甲方NC是否获得低背景浮力响应，实际EM—热短跑采用什么假设？
5. 法国F2到底画出了什么流动，是否有真实双胞证据？
6. Fig.8温度图最大差异在哪里，下一主因是什么？必须附图。
7. 两案是否经统一CLI运行并产出可读取场文件，是否仍有静默忽略参数？
8. Linux机需要的代码、粒子、几何、参数与结果是否齐全？还不能声称什么？

本轮之后的工作切换：Linux机开始桌面装配、参数、启动/停止、进度、结果加载和打包；物理问题按KNOWN_LIMITATIONS单独续进。不要把法国完整论文验收和Linux软件开工继续绑定。

<!-- ===== END 242_OPHELIE_FINAL_PHYSICS_PREPARATION_AND_LINUX_ENTRY_CURSOR_PLAN_2026-11-24.md ===== -->


<!-- ===== BEGIN 243_OPHELIE_72H_PHYSICS_CLOSEOUT_AND_FRENCH_REPRODUCTION_CURSOR_PLAN_2026-11-24.md ===== -->

# 243_OPHELIE_72H_PHYSICS_CLOSEOUT_AND_FRENCH_REPRODUCTION_CURSOR_PLAN_2026-11-24.md

# 243：答辩前三天物理推进计划——法国热结构、搅拌异常与可审查证据

日期：2026-11-24。
执行对象：物理开发机上的Cursor。Linux桌面软件已由另一条工作线开始，本任务不开发GUI。
前置：Doc242实际运行代码/补丁、原始结果、239物性提取包、Doc204 Fig.7/8数字化参考。

## 0. 用户目标、优先级与硬边界

用户还有约三天答辩。法国自然对流是物理展示主目标；甲方搅拌需继续修复，但不得拖住法国主线和Linux软件。

本轮不承诺三天内必然实现论文完全定量复现；目标是最大化得到可信结果的概率，交出可解释的场图与量化差距，同时避免临近答辩引入失控改动。

优先级：

1. **P0：纠正归因与锁定法国工况，核实温变物性真实参与更新。**
2. **P1：用合理非均匀初猜推进法国耦合，形成Fig.8温度场＋Fig.7流动＋Fig.6热源联合证据。**
3. **P2：对甲方搅拌做最小静水诊断，定位第一个异常算子，实施一项有证据的修复。**
4. **P3：补齐可重生成图表、数据、能力标签和软件增量交接，最后一天以复跑与报告为主。**

建议物理GPU时间分配：法国约60%，搅拌约25%，必要回归/整理约15%。不是平均分配每个未完成任务。单GPU正式求解串行，CPU源码审查与图表整理可穿插。

禁止重开全套TEAM7/方腔、新造Maxwell或流体框架、引入ghost/δ-SPH新路线、做全参数/全分辨率扫描。已通过基础验证保留，只有相关核心代码变动才补对应必要回归。独立A-phi/Laplace研究不修改。

## 1. 已核查事实与必须纠正的结论

### 1.1 已经可以直接利用的基础

- 法国固定域：同一环带温场OFF零流、ON有清晰响应。当前能够产生轴降—环升—缘降的趋势，不应继续重复论证浮力公式能不能驱动流动。
- 甲方NC固定平坦滑移顶：同温场OFF静止，ON U_rms约3.39e-5m/s，浮力信号可分辨。开放自由面路径仍有背景问题，保留为独立限制。
- 两案CLI已经能prepare/run并收集场文件。软件可以继续集成，不等待法国论文完全验证。
- 甲方电磁—热短跑使用ASSUMED玻璃吸收15kW；电源30kW不能自动解释成50%效率，等效电流反算也不是甲方实测电流。

### 1.2 法国242结果不能这样归因

242报告将5s散热约154kW、入热50kW归为`heat_loss_dominant_cold_wall`，并推断“继续延时只会更冷，无法形成论文结构”。**该推断未证实，本轮改标：均温起步的短时净冷却，边界实现正确性与最终平衡未定。**

实际5s末：Tmin≈1440.55K，Tmean≈1468.3K，Tmax≈1477.99K。初始1473K，最低下降而最高上升，温度梯度正在形成。

质量约97.21kg、cp约1150J/(kg·K)，热容量约112kJ/K。净散热约104kW对应平均降温约0.93K/s，与5s平均下降约4.7K相符；这不能单独证明h设错。

242法国日志同时写明：`constant_thermal_rheology=1`，μ=4Pa·s、k=4W/(m·K)、cp≈1150J/(kg·K)。这不是完整温变物性论文复现。σ有温变路径也不代表热/黏度全部温变。

5s温度差距不等于最终稳态误差；约20次论文耦合迭代不等于20秒物理时间。不能为了让5s图好看就修改已知边界或增加β。

### 1.3 搅拌仍需承认的问题

- 242静止桨仍出现异常流动，关闭流体向固体载荷也未消除尖峰，因此首先查流体平衡、几何/边界与界面更新，不优先软化弹性材料。
- 加静水IC只降低早期波动；约0.012–0.015s仍可有3.5–4.3m/s尖峰，不能认定“缺初始静水压”解释全部问题。
- CLI末Umax≈1.24m/s，密度最大相对偏差约54%，数值质量未过。
- 当前固体候选dt有打印，真实子步仍取启动dt_s0；不能称当前稳定步长已经控制推进。
- `hold-hydrostatic-ic`每次重置密度/压力仅是诊断，不能作为默认修复。

## 2. 启动前最多一次资料定位与版本冻结

先读当前源码与有效配置，不把旧包HEAD当现在唯一版本。用户已在其他机器推进软件，不覆盖对方分支。

保存：本机HEAD、分支、WIP、binary hash、工具链/设备、关键STL/reload/物性表hash。新实验在独立目录，保留242的可运行版本。

本机需要的资料按下列入口寻找；不依赖ChatGPT页面、sandbox链接或会话内脚本：

| 资料 | 应用位置/实际入口 |
|---|---|
| 当前法国主case | `tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/` |
| 当前甲方弹性case | `tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_f0_elastic_fsi/` |
| OPHELIE共享物性 | `tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/`下`electromagnetic_ophelie_french_material_laws.h`及相关局部算子 |
| 242法国基线 | `doc242_runs/05_french_paper_joint/PAPER_NC_50kW_5s/`，本机通常位于docs/ophelie/reference下 |
| 242搅拌 | `doc242_runs/01_elastic_diagnostics/`、`02_elastic_fix/`、`06_cli_e2e/` |
| Fig.7/8参考 | 242包内`reference_fig7_fig8_Doc204/reference/`及provenance/说明 |
| 物性提取 | 已交付239数据包：Jac/Sauv曲线与公式来源、拟合及QA |
| 软件适配层 | `docs/ophelie/software/linux_src/`，已有validate/prepare/run与result-index |

目录不同时只做一次定位并写入`INPUTS.md`，缺件列具体文件，不遍历所有历史包反复读计划。已有图像和数据不得假称没有；当前模型不支持的材料项明确标未接入。

## 3. P0-A：锁定法国主配置，隔离诊断与论文目标

### 3.1 本轮唯一论文主身份

`FRENCH_PAPER_NC_50KW`：期刊无热电偶装置、无机械搅拌自然对流，平坦固定不可穿透/切向滑移顶面，侧/底无滑移。优先对齐Fig.8温度、Fig.7速度趋势、Fig.6热源。

人工F2环带另保留`FRENCH_F2_DIAGNOSTIC`，不是正式初值来源真值，不把冻结温场图当耦合复现。

Jac博士55kW、Marangoni组别、Sauv另一玻璃组别不得混成同一个验证身份。每次输出物性dataset_id、property_model_id和论文/借用/假设来源。

### 3.2 本轮必须打印的有效参数表

- R/H、坐标轴与重力、粒子dp、初始质量/体积。
- ρ、β、μ/k/cp/σ公式或表、上下界与温标；℃转换仅在入口发生，求解内部用K。
- 线圈半径/匝数/高度/频率、输入电流幅值约定、固定吸收功率50kW及归一因子。
- EM形式、自感开关、σ更新方式和间隔。
- 侧/底/顶h、环境温度、辐射ε与辐射环境温度、边界面积。
- 顶面机械条件、压力/静水背景分解、完整重力是否另施加、Marangoni是否启用。
- 初猜来源、温度是否随时间推进、实际控制分支。

期刊已知参数固定；未知项才用博士补充或假设。参数表是后续所有run共用基线，不每次临时拼装。

## 4. P0-B：局部温变物性接线，先查后补

### 4.1 查真实调用链，不看数组就授PASS

输出`PROPERTY_PATH_AUDIT.md`，每项列：数据来源 → 更新时机 → 字段/函数 → 实际消费算子 → 步长限制 → 能量或动量检查。

| 项 | 必须核实 | 本轮顺序 |
|---|---|---|
| μ(T) | 粒子局部值进入黏性相互作用/隐式系统，不是体积平均一个μ | 第一优先 |
| k(T) | 流体内部与接触导热采用一致局部系数 | 第二优先 |
| cp(T) | 实际温升换算使用局部cp/焓，不是构造函数标量覆盖 | 与k一起审计 |
| σ(T) | 当前温度更新到真实EM解，非仅可视化字段；时间戳明确 | 接回正式耦合 |

特别审查`ApplyLaplaceDiffusionDtCK`及热边界消耗cp的路径。复用已有局部实现，避免重写重复算子。

### 4.2 最小功能验证，不重跑一套论文case

1. 选3–5个覆盖实际温区与曲线转折的温度点，记录μ/k/cp/σ；与已提取表插值结果核对，检查单位、单调区段和截断。
2. 无热源绝热的小型热交换：实际温变cp下按一致焓差检查总能量，目标漂移<=0.5%（工程筛查）；失败先查公式/双计数，不继续正式耦合。
3. 局部μ算子若修改，用已有可信黏性验证的一个必要短回归。恒μ时恢复原结果；局部μ用同一离散的残差/耗散与CPU小规模复核，不新造大型验证项目。
4. 在非均匀正式初猜上打印物性min/mean/max及实际消费统计；高温平台不变可以合理，不能为了min!=max故意使用无物理根据温度。

温变cp应以`h(T)=积分cp(T)dT`定义储能；若使用显式近似，报告离散能量残差和步长敏感性，不能以当前cp乘总温差充当严格焓。

### 4.3 冷区大黏度的处理

μ(T)可能跨多个数量级，必须报告最大黏度比、显/隐式路径、线性残差、实际dt与迭代次数。高μ导致慢，不等于浮力公式错误。

优先沿已有可信隐式黏性实现做局部系数一致化；不在三天内另写全GPU隐式框架。已有CPU/混合可信路径可以作为备选，准确标设备与耗时。

不得静默截断μ、直接冻结冷区速度或把高μ区域删掉。如果需黏度上限作数值延续，只能在独立DIAGNOSTIC身份下逐级接近目标上限；最终展示注明限制，并用目标物性作回查。续进阶段结果不是论文定量验证。

若cp局部化短时间无法完成，不无限卡住主线：先用已验证常cp、局部μ(T)/k(T)做明确标注的中间模型，输出限制和对应误差风险；不得称完整温变物性复现。

## 5. P0-C：热边界审计，不先调h

### 5.1 几何与热量核对

- 先用圆柱几何解析面积：侧2πRH、顶πR²、底πR²，与离散有效面积积分比较；工程目标<=5%，边角归属规则明确。
- 检查侧/底/顶分类重叠及遗漏；机械顶墙不应额外引入第二份热损失。
- 对流与辐射各自报告面积、采用的边界温度、环境温度和符号；净流入与净流出约定统一。
- 顶面辐射用绝对K；冷水/空气/辐射环境不能用同一个默认值静默覆盖。
- 用均匀T的单步检查与解析`hA(T-Tenv)`、`εσSB A(T^4-Trad^4)`核对；此时可检验面积与单位，正式非均匀场需按各面温度积分。
- 边界扣热与总储能采用同一时间区间、质量/体积和cp约定；固定T诊断的维持温场功单列，不混入正式能量守恒。

单步输入功率/边界计算目标差<=1%（排除已量化面积差）；累计正式热收支残差工程目标<=2%，>5%先停下来排查。不用大floor掩盖残差。

### 5.2 审核结果决定是否修h

面积/符号/温度层有实现错误：先修实现，保持文献h。

实现一致但参考不明：保留ASSUMED身份，先运行真实热演化；只有出现可比热平衡且偏差仍与散热相关时，再选择一个有范围依据的未知边界参数敏感性。

不允许由5s散热/入热>1直接推导必须减h；也不允许仅为T*接近0.5–0.9而调低绝对温度。

## 6. P1-A：构造可自由演化的合理非均匀初猜

### 6.1 优先顺序

首选已有热模块或现有轻量轴对称热计算生成冷壁梯度；如果已有成熟初猜文件则直接复用并查身份。不要从零开发通用2D求解器。

热初猜过程可暂冻结速度并使用一份固定Q，只为加快建立热梯度；若用有限体积轴对称辅助工具，须包含r权重、轴线对称、实际顶部/底部/侧壁条件和积分功率。它是初猜生成器，不是SPH验证真值。

Doc204稀疏等温线可作有不确定性的备选初猜，保留插值范围/掩码与假设Tmax，不把未识别区域补成“论文完整场”。

### 6.2 导入与正式推进要求

- 正确处理r/z原点、m单位、粒子采样位置、轴线和越界点。
- 导入后实际统计Tmin/Tmax/均值、剖面、总焓与冷/热区体积分数；输出初态场图。
- 初态速度不指定目标环流，压力/背景初始化与当前固定域路径一致。
- 正式run不得开启`freeze-prescribed-t`；初猜导入一次后，温度由热方程、对流和Q共同演化。
- 不得每步覆写Fig.8温度，或每步重置流体速度，以伪造稳态。
- 若初猜温度超物性有效范围，标记并采用有记录的外推策略，不静默截断后称文献精确值。

### 6.3 初猜独立性

主配置取得可比热结构后，再做第二个合理初猜的短续算（例如适度扰动初猜梯度/均值，保持物性和边界不变）。若二者呈趋同则增强稳态可信度；时长不够只标未验证。不要为完成此项提前铺开多组长跑。

## 7. P1-B：有限矩阵推进法国耦合，不固定要求超长物理时间

### 7.1 正式矩阵

| 组 | 内容 | 是否必须 |
|---|---|---|
| F_BASE | 242常物性5s原始结果，作为历史对照 | 复用，不重跑 |
| F_VAR | 真实局部物性＋合理初猜＋50kW＋当前可信EM路径 | 主任务 |
| F_SEED | 主模型第二初猜短续进 | 主结果有希望后执行 |
| F_SENS | 一个有证据的σ/边界/EM数值敏感性 | 最多选一项，不全开 |

先不改变dp=15mm，不扩展全部线圈/频率/功率。不得同时改变初猜、h、β、功率再声称找到唯一根因。F_VAR是更完整模型，不用于声称单项因果；如果要归因局部μ，用短的相同初猜常μ/变μ对照即可。

### 7.2 运行时长由演化与耗时共同决定

先2s仪器化短跑，测每物理秒耗时、步长、黏性/EM迭代、温度变化与稳定性；只做必要输出，不能每步全量VTP。

健康后逐段推进到10s、30s；有热结构发展且预算允许时才继续到60/120s或现有重启的更长窗口。这些是检查点，不是强制总时长或稳态判定。

如果没有可靠checkpoint，不临时构造只存T的“完整重启”。可以单次run定期保存统计和有限场帧；热初猜导入与完整状态重启分开命名。

每段判断：

- 数值有效、热收支健康、T结构向目标发展：允许继续。
- 高μ令dt/隐式迭代失控：先查局部黏性实现/预处理/容差，不靠把温变物性关掉宣称完成。
- 热结构几乎不变且离目标很远：先判是否已接近平衡、Q/边界/初猜是否有问题；不无限跑。
- 场仍明显发展：标TRANSIENT，不以5s差距定最终误差。

第一天必须交2s/首段图；长任务每30–60分钟写进度/有效物理时间，不允许一天后才发现零推进。

### 7.3 EM方面只做必要的实场核查

当前固定功率50kW控制可以保留；同时报告归一前功率、等效电流及归一后功率，避免50kW吻合掩盖形状问题。

σ(T)先选明确定义的主曲线（239提取的图拟合/印刷公式不混用），写模型ID及有效区间。若需A/B，用同一冻结非均匀T比较原始与固定功率Q分布；最多两次EM解，不同时热演化A/B以免失去可比性。

A_ind=off必须标模型差异，不因本次“皮肤状Q”就声称已完整考虑自感趋肤。只有同温场下Q趋势与Fig.6明显不一致、且现有自感路径已有可信验证时，才将其选为唯一F_SENS；不临时启用未验证自感强行改善图像。

EM更新间隔用固定主值，记录每次采用的T/σ/Q时间戳。需要数值敏感性时仅主值与减半两档；不要与h/σ曲线同时扫描。

## 8. P1-C：统一场图、稳态判定与参考误差

### 8.1 每次同窗输出T/U/Q三联图

- Fig.8：绝对T与原文定义的T*=T/Tmax，分别保留；Tmax绝对来源不明不反推唯一温标。
- Fig.7：r-z周向平均速度矢量/流线与中心薄片；主次胞必须在同一半径剖面判断，左右镜像不能当两个径向胞。
- Fig.6：Q位置/范围、径向与轴向累计功率，定性图像或数字化结果的精度分开说明。

所有图含坐标单位、刻度、色标、图例、实际时间/平均窗口、模型身份、数据/脚本来源。参考只在已识别轮廓/点区域比较，不将稀疏插值场作为全域精确真值。

### 8.2 必交数值

Tmin/Tmax/Tmean、热点r/R和z/H、三条竖向/三条水平温度剖面、主要等温线形状；轴/内部环带/最外缘Uz、Umax/U_rms和回流旋向；Q积分/热点与高热源区域；每面散热和总焓变化。

两组图使用统一色标，另有自适应色标时明确说明。不能仅用各自色标让接近均温的场看起来梯度很强。

### 8.3 工程准稳态筛查，不能伪造论文标准

用三个连续、等长、完整后期窗口，窗口长度根据实际热/流动演化确定并保持不变：

- 储能速率/输入与散热功率规模<=5%；
- Tmean每窗变化<=代表温差的1%；
- 热点位置及主要等温线无明显持续漂移；
- U_rms和主要分区Uz邻窗变化<=5%，低信号同时报告绝对量；
- 热收支残差符合第5节要求，密度/体积/粒子数无异常漂移。

满足才标ENGINEERING_QUASI_STEADY，不自动变成PAPER_VALIDATED。不满足则保留TRANSIENT并显示仍在发展的量。

没有可靠参考误差定义时不给合同10%或自造“论文PASS”。形状接近、方向一致、量级接近分别报告，不能用一个百分比概括全部。

## 9. P2：甲方搅拌只定位第一个异常，避免继续堆初始化补丁

### 9.1 复用已知对照，不重做242矩阵

B1静止桨异常、B2去载荷尖峰仍在是已有证据。保留5mm密粒子、不锈钢材料；不软化E，不用hold-hydrostatic作为默认。

先运行零转速、等温、无EM/无浮力静水诊断。新增局部诊断开关必须明确作用范围，不改其他case默认。

### 9.2 从三种状态区分流体与弹性运动

- S_FIXED：桨完全固定，作为静止墙接触，流体有重力。诊断静水/壁面平衡。
- S_ELASTIC_REST：桨真实弹性、驱动为零，其余相同。判断固体重力/约束/接触运动是否使问题放大。
- S_DRIVEN：只有前两者有清楚结论后，恢复10rpm驱动。

固定桨仅为诊断，不替代最终弹性模型。若S_FIXED仍异常，不优先重写固体积分器。

### 9.3 抓首个坏步而非最终Umax

输出声学步内以下阶段前后粒子级增量：密度求和/重初始化、压力更新、完整重力、内压力力、墙/桨接触压力、黏性、位置更新、cell/邻域更新。以实际调用顺序为准，不凭任务书重排成熟算法。

每阶段记录Umax/U_rms、rho极值及P99、极值ParticleID与位置、压力、邻域数；追踪前10异常粒子。用稳定ID处理排序。

核对：

1. 静水IC后的ρ/p是否在首个密度求和时被重写；若被重写，应修相容的密度/压力/边界路径，不能每步强刷静水解。
2. 重力是否重复加、ForcePrior是否正确清理/累加、静水背景与完整重力是否混用。
3. 几何是否单位或轴变换错误、流体初始是否侵入桨/壁、邻域是否在实际位置更新后刷新。
4. 法向与接触体积是否一致，接触力符号是否正确，桨平均速度与位置时间层是否匹配。
5. c0=8下合理静水压缩应按EOS和深度估计；不能将54%偏差归为正常弱可压缩效应，也不能只按rho/rho_ref无视静水基线。

必要时做唯一一次`dt/2`对照；波动仍同量级则不继续缩dt扫稳定。若c0需敏感性，只在静水配置中比较原值与合理一档，不能用极大c0掩盖不平衡且不报告耗时。

### 9.4 固体dt真实参与控制

将live稳定候选用于实际dt_use约束，非法候选应FAIL而非继续使用更大固定dt。子步和必须覆盖dt_fluid，不满足即停；记录dt_candidate/dt_used/covered_time。

该项是代码欠账，不能预先归为所有高速尖峰根因。S_FIXED异常仍应先修流体路径。

### 9.5 本轮只修一个有证据的主问题

选择第一个明显不平衡阶段，最小局部修复，保持相同初态作修前/修后比较。若仍有第二问题，先保留证据，不把全部离散算子同时替换。

静水目标：相对合理静水基线的密度异常明显降低，无持续自激流动；速度背景足够低于旋转驱动信号。P99密度偏差<=2%、最大<=5%可作工程筛查；边界层单列但不从总统计隐藏。失败不能通过调阈值转PASS。

通过静水后再跑真实弹性0.05s，若质量健康可到0.2s；本轮不强求整圈。报告粒子逃逸、固体Jacobian/应力有效性、桨速与流速量级。仍异常则软件保持诊断标签，不渲染成可信工程结果。

## 10. 三天安排与停止规则

按相对开始时间计，不受机器时区影响。以答辩前至少12小时冻结演示版本为原则。

| 时间段 | 法国主线 | 搅拌/整理 | 当段必须留下 |
|---|---|---|---|
| 0–6小时 | 有效身份、物性调用链、热边界单步审计 | 版本冻结、旧图脚本补全 | 参数表/缺口表/首个明确修复点 |
| 6–24小时 | 局部物性与初猜，F_VAR首段实际推进 | 静水首坏步诊断，限制投入 | 首张真实T/U/Q三联图、进度与耗时 |
| 24–48小时 | 条件健康则分段续算；最多一个敏感性 | 一项搅拌修复和必要回归 | 后期窗口趋势、修前修后、可复跑版本 |
| 48–60小时 | 主模型稳定性/第二初猜短检查，优先复用已有长结果 | 增量交软件机、补报告 | 冻结候选及完整图表 |
| 最后12小时 | 不再改核心物理；只复跑/导出/检查 | 软件集成与演示备份 | 答辩可用包、明确限制和运行说明 |

物性局部化某一项超过约4小时仍无清楚路径时，提交调用链与阻塞点，采用有标签中间模型继续有效计算，不整天盲改。边界审计只做必要解析核对，不扩展到未知结构全建模。

任何正式run在首个非有限量、非正密度/体积、非法固体Jacobian、子步欠时间、明显不守恒时停止并保留首坏步。出现这些问题不允许继续长跑“看看会不会自己稳定”。

若48小时仍未达论文热结构：选择最可信模型与瞬态趋势诚实展示，附差距原因与下一步，不临时强写目标温场/速度。同时保留人工F2作为独立功能证据，但清晰区分。

## 11. 与Linux软件机的接口和增量交接

软件机继续GUI、参数、结果加载、打包；物理机不重构界面和重复CLI。不因软件已开工就冻结必要物理修复，也不让软件实时跟随未经检查的物理WIP。

两机独立分支。物理修复以小提交或base+完整patch交付，写接口是否变化、配置是否兼容、资产是否重生、最小复跑命令。软件机显式cherry-pick/合并；不强推、不覆盖对方未提交代码。

每次最多交一个已检查候选，稳定配置与实验配置分开。状态字段至少包括：workflow_completed、numerical_quality_ok、em_executed、thermal_advanced、flow_advanced、solid_advanced、steady_state_status、paper_comparison_status、assumptions。

第三方展示中必须区分真实计算、历史回放、人工温场诊断、校准结果、文献预测。搅拌质量FAIL不能因为GUI能加载就改名“已验证”。

## 12. 交付物：少读文档，多看可复现证据

目录建议`doc243_runs/`：

```text
00_freeze/                    # 版本、设备、二进制、资产hash
01_model_and_properties/      # 有效身份表、调用链、函数/单步检查
02_boundary_audit/            # 面积、解析热流、能量审计
03_thermal_initial_guess/     # 初猜来源、导入检查、非冻结证明
04_french_coupled/            # F_VAR主线、后期窗、有限敏感性
05_client_hydro_fix/          # 首坏步、固定桨/弹性静止诊断、局部修复
06_figures_and_tables/        # 可直接答辩的PNG/PDF/SVG和CSV
07_software_increment/       # 最小代码/配置/资产增量与回归命令
INPUTS.md
RUN_COMMANDS.md
ROUND243_DECISION.md
KNOWN_LIMITATIONS.md
DEFENSE_EVIDENCE_INDEX.md
```

必须有真正可执行绘图脚本，指定输入目录和输出目录，禁止再交只提示“使用会话内脚本”的占位文件。普通Python绘图即可，不新增大型可视化依赖。

每组保留effective配置、日志、原始CSV、必要VTP、统计窗口和源代码hash。精简包至少包含每个关键结论对应的场帧/CSV，其他完整场可本机保留并列清单/hash，不用重复上传全部历史。

### 答辩证据至少六张图/表

1. **模型对齐表**：期刊已知、博士借用、假设/校准及当前差异。
2. **物性生效图/表**：真实μ/k/cp/σ路径与实际温区，说明仍常量的项。
3. **法国热收支随时间**：输入、三面散热、储能变化及残差。
4. **法国同窗T/U/Q图**：参考与计算，同坐标/色标，时间状态明确。
5. **温度/速度关键剖面与后期变化**：说明接近程度与是否仍在发展。
6. **甲方搅拌异常修复证据**：首坏步位置、修前修后与质量状态，不只展示动画。

若最终只达到部分目标，DEFENSE_EVIDENCE_INDEX必须逐项说明支持什么、不支持什么。单模块验证由已有验证材料另行支撑，本轮不以工程短跑代替标准case/实验误差证明。

## 13. Cursor最终汇报只回答这些问题

1. 法国当前实际使用哪些温变物性，哪些仍为近似？给真实调用证据。
2. 热边界有没有实现错误？还是只是起步净冷却？给面积/单步/累计证据。
3. 初猜从哪里来，是否导入一次后自由演化？运行到了多少物理时间？
4. Fig.8温度结构是否接近、最大偏差在哪里，是否已有准稳态证据？
5. Fig.7轴降—环升—缘降是否在正式耦合中出现？小外围胞证据是否足够？
6. Q位置与Fig.6的差异是什么，当前自感/σ模型有哪些限制？
7. 搅拌首个异常算子定位了吗，修复实际改善多少，质量是否仍FAIL？
8. 交给Linux机的准确版本、兼容配置和不可夸大的结论是什么？

不要用新增开关数量、运行次数、授牌数量概括进展。以实际温度结构、收支、速度场与异常修复证据作为本轮完成标准。

<!-- ===== END 243_OPHELIE_72H_PHYSICS_CLOSEOUT_AND_FRENCH_REPRODUCTION_CURSOR_PLAN_2026-11-24.md ===== -->


<!-- ===== BEGIN 244_OPHELIE_TARGETED_PHYSICS_DEVELOPMENT_AND_TEST_PLAN_2026-11-24.md ===== -->

# 244_OPHELIE_TARGETED_PHYSICS_DEVELOPMENT_AND_TEST_PLAN_2026-11-24.md

# 第244轮：法国温度场整改与搅拌首步异常定位——开发测试任务书

日期：2026-11-24。
执行对象：物理开发机上的Cursor。Linux软件机继续独立开发。
版本：修订版；替代此前第244轮任务书。只执行本文件，旧版未保留的任务不执行。

## 0. 本轮范围：只做两条物理任务，不扩展路线

1. 法国自然对流：纠正空间热平衡判断，恢复符合原文归一关系的温度初猜，检查真实冷区物性和数值刚性，推进自由演化的热—流动计算。
2. 搅拌：用真正固定的桨做静水诊断，分析首个声学步各算子的贡献，定位并修复第一个明确错误。

本轮不开发GUI，不重跑整套方腔/TEAM7，不新建流体框架，不新增ghost/δ-SPH路线，不全域扫描粒子分辨率，不拟合文献参数，不改写目标图或求解输出。

目标是两项具体证据：**法国正确温度定义下的冷区处理和实际场图；搅拌具体哪个算子、哪些粒子首先失衡。** 不再以增加开关或阶段PASS标签代替解决问题。

## 1. Doc243已审查事实：保留成果，纠正结论

### 1.1 保留现有有效改动

- 逐粒子cp热扩散更新`ApplyLaplaceDiffusionDtParticleCpCK`。
- 已接入的局部μ(T)、k(T)和σ(T)路径；仍需补足必要消费/守恒检查。
- 无附加壳阻的30s耦合结果：T约1291–1553K，轴心下降、内部外围区域上升，U_rms约4.42e-4m/s。作为有效研究基线保存，不重跑作废。
- 搅拌跳过运行时dens-sum的尝试已失败，恢复正常默认；不重复同一个开关试验。
- 两案现有CLI和结果索引，不重构软件接口。

### 1.2 必须纠正的报告表述

| 原结论 | 本轮处理 |
|---|---|
| 全域均温1473K散热157kW，故50kW与论文温区不兼容 | 撤回推断；只说明该均温初态净冷却 |
| 0-D平衡温度778/1260K决定空间模型最终温度 | 不采用；集总均温模型不能替代非均温空间解 |
| 缺冷壳资料已被证明是根因 | 未证明；附加热阻降为独立敏感性历史结果，本轮不继续调它 |
| Fig.8初猜1400–1550K失败，故参考温场不可用 | 原初猜改动了Tstar比例，不能用来否定论文温场 |
| 搅拌首坏步定位完成 | 改为仅缩小到流体/边界侧；具体算子仍未定位 |
| hold重刷密度使偏差降低，证明dens-sum是唯一根因 | 不成立；必须直接测该算子前后增量 |

先更新本轮结论，不删除历史原始记录，避免后续继续引用错误推断。

## 2. 启动前冻结与资料定位：最多做一次

保存实际HEAD、分支、WIP、binary hash、设备、STL/reload/物性表hash。Doc243包的旧HEAD不等于现在软件机或物理机代码，不擅自回退。

关键源码：

```text
tests/extra_source_and_tests/3d_examples/
  test_3d_ophelie_french_natural_convection_frozen_q/
  test_3d_ophelie_client_stirring_f0_elastic_fsi/

tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/
  electromagnetic_ophelie_french_material_laws.h
  electromagnetic_ophelie_french_skull_resistance.h
  diagnostics/electromagnetic_ophelie_thermal_diffusion_one_way.h
```

输入：Doc243本机完整运行记录、239物性提取资料、Doc204 Fig.7/8标定/等温线/部分重建数据。一次定位后写入`INPUTS.md`，不要反复翻阅全部历史任务书。

Doc243精简包未装完整VTP与全日志；需要的数据优先在原机找，确实缺失再最短重跑，不从summary反推不存在的粒子证据。

本轮法国主身份：期刊50kW、无机械搅拌、无热电偶装置、固定平坦滑移顶面；侧/底无滑移。博士资料仅补未明确项并标来源，不混入55kW/Marangoni组作为同一验证目标。

## 3. 法国任务F1：用非均匀边界温度判断热平衡

### 3.1 保持当前已知物理参数，先检查实现

正式主配置关闭额外`skull-resistance/free-crust-resistance`，保留其源码和历史结果，不在本轮继续调壳厚。

实际收支应为：

```text
P_in = ∫V Q dV
P_side = ∫side h_side*(T_face-T_cool) dA
P_bottom = ∫bottom h_bottom*(T_face-T_cool) dA
P_top = ∫top [h_top*(T_face-T_air)
              + epsilon*sigma_SB*(T_face^4-T_rad^4)] dA
储能变化率 = P_in - P_side - P_bottom - P_top + 数值残差
```

T_face必须是边界局部温度，不能用内部代表温度、Tmean或Tmax覆盖所有边界。均温0-D只用于公式/初态数量级核对，不决定空间模型是否能形成热平衡。

### 3.2 只补Doc243未证明的关键项

Doc243已完成面积归一和均温解析散热核对，不重做。补输出：

1. 侧/底/顶每个边界粒子的温度、对应面温、面积权重及热流。
2. 各面面积加权温度与温度范围，各面实际功率。
3. 粒子中心温度与几何面温的关系；是否有外推、其距离/方向/梯度来源。
4. 热流积分与实际扣除储能是否对应同一时间段、同一cp/焓定义。
5. 不同面的归属是否重复；机械顶墙是否额外增加一份热损失。

特别检查侧面是否仍把热粒子中心温度直接当冷壁面温；若存在此近似，量化其范围，不直接宣布错误或添加经验热阻。

若要修改面温处理，必须使局部传导和Robin边界一致；不把已知换热条件改成固定温度，不用任意低温边界凑结果。

### 3.3 完成标准

交`BOUNDARY_SPATIAL_AUDIT.md`和原始CSV，明确：已经证实的实现问题、修复内容、仍属分辨率近似的项。没有新错误则保持边界参数不变，继续F2，不因0-D估算再次停止主线。

## 4. 法国任务F2：正确构造Fig.8初猜，正面检查高黏度

### 4.1 保持原文归一关系

现有参考定义：`Tstar = T/Tmax`。正确重建形式：

```text
T_seed(r,z) = Tscale * Tstar_reference(r,z)
```

Tscale精确值未确定时标ASSUMED；原文名义“约1500K”不等于精确Tmax。先登记一个有来源说明的诊断值；不要启动Tscale扫描。

禁止再次把Tstar∈[0.5,0.9]仿射映射到1400–1550K并称保留Fig.8温度结构。旧热温带文件保留为人工诊断，不用于本轮正式参考判断。

在已识别参考点检查`T_seed/Tscale`与原Tstar是否一致。未知区域插值规则、有效区域掩码、外推方法和误差明确记录，不把部分重建当原作者完整数值解。

### 4.2 先查物性与dt，再启动长计算

仅导入正确初猜，输出：

- Tmin/Tmax/Tmean、冷区体积分数、总焓、初态温度剖面。
- μ/k/cp/σ的min/max/P50/P95/P99及空间分布。
- 每条物性曲线有效范围、外推/截断数量和原因。
- 真实声学、热扩散、黏性步长限制；显/隐式黏性实际执行路径。

若没有准确Tscale，先用名义1500K作明确假设的物性/步长诊断，不把结果写为论文绝对温度已确定。

### 4.3 处理数值刚性的顺序

1. 查现有可信隐式黏性是否支持局部μ及当前边界；支持则复用，记录残差、迭代、耗散和实际dt。
2. GPU路径不适配时，可用已有CPU/混合可信实现做小规模验证；准确性优先，不临时编写大型GPU线性求解器。
3. 显式dt合理时先极短自由演化，不因μ大就认定模型错误。
4. 真实目标μ暂不能稳定推进，报告最严粒子/算子/dt与耗时瓶颈；不要静默压缩温差、截断μ或冻结冷区速度来授PASS。

本轮不新增黏度上限扫描。已有上限仅能作为明确标注的数值诊断，不能替代目标物性结果。

## 5. 法国任务F3：最小矩阵与真实热演化

| 编号 | 内容 | 运行要求 |
|---|---|---|
| F0 | 243无附加壳阻30s基线 | 复用，不重跑 |
| F1 | 正确Fig.8初猜的导入/物性/dt检查 | 必做，先不长跑 |
| F2 | 原边界、正确初猜、实际局部物性，极短自由演化 | 初始检查健康后做 |
| F3 | F2健康后的续进；或复用已有可信非均匀热初猜自由演化 | 最多一条主计算线 |

F2/F3不冻结温度、不指定目标速度，不在每步重新刷Fig.8温场。初猜导入与完整状态重启分开命名：仅保存T不能冒充完整checkpoint。

先极短窗口检查，再根据健康度和实测耗时到2/10/30/60s检查点。这些是观察点，不是稳态时间承诺。已有30s数据不重复起步计算以凑运行次数。

仅因10s外围Uz未向上，不能判最终失败；仅因净散热不能判边界错误。可以因时间预算停止，但结论写“尚未达到/无法在预算内完成”，不写“已经证明不可能”。

当前dp15mm可能不足以解析数毫米冷壁层。先查边界面温/热流近似，必要时利用已有轻量轴对称热模块生成更合理初猜或作独立热诊断；不在本轮新建大型辅助求解器，不全域加密扫参。

## 6. 法国输出和判断标准

每个正式结果同窗输出：

- 温度绝对值与Tstar、参考等温线叠加、热点位置；中心/中环/近壁竖向剖面与上/中/下横向剖面。
- r-z周向平均Ur/Uz矢量和流线；保留真实中心薄片识别非轴对称。
- 原始Q与固定功率后的Q分布、功率积分及径向/轴向分配。
- 输入功率、各面散热、总焓变化与残差时间曲线。

报告轴/内部环带/最外缘Uz，不能把宽泛“外围上升”当成最外侧冷壁小胞已经形成；左右镜像也不算两个径向胞。

T/U/Q图必须有坐标单位、刻度、色标、图例、实际时间/平均窗、模型身份。参考和计算同色标同坐标，不用各自色标隐藏差距。

物性接线至少补：3–5温度点函数核对、局部cp的绝热热交换/焓守恒短检查。已有结果充分则复用，不新建完整验证项目。

正式累计能量残差工程目标<=2%，>5%先查实现；低分母用绝对值，不用大floor。准稳态需热收支和连续后期窗口一起判断；未达到标TRANSIENT，不授PAPER_VALIDATED。

本轮完成标准是：正确温度定义、真实冷区处理证据、可读实际场图与明确差距，不强行承诺论文完全复现。

## 7. 搅拌任务S1：真正固定桨的首步静水诊断

### 7.1 固定混杂因素

采用当前5mm密reload，等温、无EM、无浮力，保留实际重力；桨完全固定为静止接触壁。

停转弹性桨仍可能运动；关流体→固体载荷也没有关闭固体→流体接触。已有B1/B2只缩小问题范围，不能替代本诊断。

先检查STL/reload坐标、轴变换、初始玻璃与桨/壁重叠、法向和邻域。无需重新审计全部历史粒子集，只查当前实际资产。

固定桨只作诊断，不把该结果标成弹性FSI验证，不替换甲方最终真实弹性任务。

### 7.2 输出粒子级阶段增量

必须从同一初态获得以下位置的真实量：

```text
S0：邻域、法向、初始rho/p/u、ForcePrior各来源
S1：密度求和/正则前后（若该步发生）
S2：声学第一半步内部压力相互作用
S3：该半步墙/桨接触贡献
S4：速度更新中的完整重力及其他力
S5：连续性/第二半步更新
S6：位置更新、下一次cell与邻域重建
```

这是观测位置，不要求按此表重新排列成熟算法。融合kernel内部可加诊断字段或小规模host复算，不把两个已演化成不同状态的关开算子run直接相减当同一状态贡献。

每粒子记录稳定ID、位置、ρ、p、u、内压力/壁压力/桨压力/重力/黏性加速度、邻域数、法向与接触体积分量。

追踪前10异常粒子，处理排序后的稳定ID，不能用数组下标冒充同一粒子。首步全量CSV应随精简包交付，数据量小，不再只给最终Umax。

### 7.3 静水平衡检查

```text
a_res = a_pressure_inner + a_pressure_contact + g + a_other
```

按内部、容器壁邻域、桨邻域、自由面分组，报告绝对残差、最大值与P99及位置。先看局部是否相消，不能用全局平均抵消证明正确。

直接量出密度求和前后Δρ/Δp；不能因为hold强刷后ρ偏差降低就宣布密度根因已证实。

## 8. 搅拌任务S2：只修第一个明确错误

按证据分支：

- S1破坏状态：检查静水ρ与密度求和/自由面正则相容性，不每步覆写解析解。
- S2/S3失衡：检查压力离散、壁面虚拟状态/加速度、法向和体积、界面几何。
- S4重力或旧ForcePrior重复：修累加/清理位置，不全域删重力。
- S5异常：查连续性与压力/速度时间层及邻域状态。
- S6后变坏：查位置更新、cell/邻域重建、变量同步和排序。

实施最小OPHELIE局部修复，先证明同一首步状态的残差降低，再到0.02/0.05s；首步未解决前不继续一排长窗口。

若证据指出需要一致的静水分解，参照已有法国路径做局部适配；不能只移一个background_hydro_p开关而留下重复重力或不一致边界，也不在没有首步证据前整体移植求解器。

允许一次dt/2核对时间离散影响，不重复扫dt。当前c0=8下静水密度变化按EOS/深度核查，不能把50%以上压缩解释为正常弱可压缩性。

固体live稳定dt应真正约束dt_use，非法候选停机；子步和覆盖完整流体dt。该代码欠账单独修，但它不能解释固定桨流体首步失衡。

## 9. 搅拌任务S3：修复后才恢复真实弹性旋转

固定桨静水健康后，按顺序恢复：弹性桨静止 → 10rpm平滑驱动 → 0.05s短测。保留真实不锈钢E约2e11Pa，禁止软化E、增黏度、清零全场速度来压异常。

输出密度/速度分位数与极值位置、ParticleID、桨速、力/转矩、弹性位移、固体Jacobian/应力有效性、粒子逃逸及每步时间覆盖。

工程筛查可用：相对合理静水基线的密度偏差P99<=2%、最大<=5%；这是本轮质量门，不是合同误差。边界层单列但不能从总统计隐藏。

流速与桨速比较仅作数量级检查，仍有1–4m/s异常时必须解释具体来源，不只以有限数值判断成功。质量未过就保持FAIL，不调阈值。

## 10. 时间约束与停止规则

本轮不再展开新的多天大计划。建议前2小时必须产出：正确初猜映射/冷区物性与dt表，以及固定桨首步增量CSV。

后续优先做：一个有证据的局部修复、法国一条健康主计算线、真实场图。单GPU串行；CPU脚本/报告整理穿插。

停止条件：非有限场、非正ρ/体积、非法固体Jacobian、子步欠时间、严重能量不闭合。保存首个坏状态和此前有效状态，不继续长跑期待自行稳定。

如果实际代码条件不支持某项：交明确调用链、复现命令、残差和具体阻塞；不自行改变研究路线，不因为“诊断已完成”跳过应做的粒子级证据。

答辩前至少留6小时冻结稳定版本、导出和软件复跑；不要在最后时段继续改核心算法。未完成项真实列明，不扩展为额外任务矩阵。

## 11. 交付与软件机关系

建议目录`doc244_runs/`：

```text
00_freeze/
01_spatial_boundary_audit/
02_correct_reference_seed/
03_french_coupled/
04_first_acoustic_step/
05_targeted_fix/
06_elastic_regression/
07_figures/
08_software_increment/
INPUTS.md
RUN_COMMANDS.md
DECISION.md
KNOWN_LIMITATIONS.md
```

必交：

1. 纠正后的热平衡说明与局部面温/通量CSV。
2. 正确Fig.8初猜、掩码、归一化核对与物性/dt表。
3. 法国实际T/U/Q图及热收支，明确瞬态/准稳态状态。
4. 固定桨首步粒子级各算子增量CSV、残差位置图。
5. 首个局部修复的修前/修后比较；无修复则具体阻塞证据。
6. 真实弹性复测结果与质量状态；未进入该阶段就写未执行。
7. 源码/配置/资产hash、必要场帧、真实可执行绘图脚本。

禁止再交仅提示“使用会话内生成脚本”的占位文件。原始求解器CSV不得被绘图脚本覆盖。精简包包含首步数据与关键结论对应场帧；其他完整数据本机保存并列清单。

软件机只接已检查的小提交或base+完整patch，列参数兼容性与最小回归命令。物理失败不改软件状态为通过，软件可运行也不等于物理验证通过。不要覆盖软件机分支或重构GUI。

## 12. 最终总结只回答这些问题

1. 是否已停止用0-D均温结论否定空间热平衡？实际各面温度与热流是多少？
2. Fig.8初猜是否保持T/Tmax，绝对温标和缺失区如何处理？
3. 冷区真实μ范围、最严dt和实际推进路径是什么？
4. 法国温度结构和轴/环/缘流动是否改善，仍有什么差距？
5. 搅拌首个异常究竟来自哪个实际算子、哪些粒子？给增量证据。
6. 修复减少了多少残差/异常，固定桨与真实弹性是否分别通过？
7. 哪些结论已证实，哪些只是合理怀疑或仍未执行？
8. 交给软件机的准确版本、能力与限制是什么？

不以新增开关、试跑次数或标签数量总结进展。本轮以实际热结构和首步异常的可核查证据为准。

<!-- ===== END 244_OPHELIE_TARGETED_PHYSICS_DEVELOPMENT_AND_TEST_PLAN_2026-11-24.md ===== -->


<!-- ===== BEGIN 245_OPHELIE_THESIS_PROPERTIES_AND_PAPER_COMPARISON_PLAN.md ===== -->

# 245_OPHELIE_THESIS_PROPERTIES_AND_PAPER_COMPARISON_PLAN.md

# 第245轮：博士论文物性实际接入、数值修复回归与法国文献温度—环流联合对比

日期：2026-11-25  
执行对象：物理开发电脑上的 Cursor  
适用范围：OPHELIE 法国自然对流分支；不修改独立 A-phi/Laplace 方法分支。  
性质：下一轮开发与测试任务书。本文中的建议接口、目录和阈值均为待落实要求，不代表已有代码已经实现。

## 0. 本轮必须交付什么

本轮不能再仅提交“短测可运行”和状态摘要，必须同时交付以下三类可复核成果：

1. **博士论文物性数据真正进入求解器。** 从原始 CSV、加载路径、插值结果到实际热扩散、热源升温、边界散热、局部黏性和电磁模块的使用位置，形成完整证据。
2. **修复后的计算可信度回归。** 重点检查粒子排序造成的边界量错位、局部物性排序与更新、隐式黏性时间步准确性、能量收支。
3. **每组主要计算自动输出与法国期刊论文相对应的图和原始数据。** 至少包含温度等高线、温度剖面曲线、速度矢量和环流图；不能再让用户仅凭一个速度箭头截图判断。

本轮优先解决法国自然对流主线。软件分支保持已交接的稳定版本。搅拌问题只安排有明确证据的有限回归，不得挤占物性接入与论文对比任务。本文不包含后备调参拟合方案。

## 1. 当前问题与本轮判断原则

### 1.1 已有实质进展

上一轮发现并修复了热边界辅助变量未随 ParticleSort 正确重排的问题。修复前，自由液面权重可能分布到整个高度，且出现负绝对温度的表面外推；修复后的短测中，顶部、底部权重位置恢复合理，负表面温度消失。这是实质性错误，必须保留回归。

当前两秒自然对流测试中，开启浮力有速度信号、关闭浮力速度为零，说明该测试的流动可以归因到浮力驱动。但这不等于温度、环流形态、稳态和论文数值已经验证通过。

### 1.2 尚未落实的关键问题

此前提取了博士论文数据，不等于最新自然对流 case 已经使用。最新审查中，运行所标注的 μ、cp、k 仍指向旧的稀疏近似表，不能用“已采用局部温变物性”替代“已采用博士论文提取数据”的证明。

已有整套 Doc239 case identity 可能同时改变功率、密度、浮力系数和边界条件。不得为了加载博士数据，顺带把法国期刊 50 kW 工况换成博士论文其他工况。

两秒属于短时诊断，不能直接判定法国自然对流无法形成。全域平均黏度对应的 Rayleigh 数也不能单独否定强非均匀黏度温度场中的局部热区流动。

### 1.3 证据等级必须区分

| 状态 | 可以声称 | 不可以声称 |
| --- | --- | --- |
| 原始数据存在 | 已提取/已收到 | 求解器已使用 |
| 表格加载成功 | 文件已读取、插值可查询 | 所有实际算子已消费 |
| 算子消费通过 | 对应模块使用了该数据 | 耦合结果已对齐论文 |
| 短测稳定 | 在指定步长与时间窗可运行 | 长时准确、达到稳态 |
| 图形类似 | 有相近流动特征 | 定量复现完成 |

## 2. 开始前冻结基线，避免再次失去可比性

创建本轮独立工作目录，保留上一轮可执行文件、代码提交号、完整命令、effective 配置、初始温度 CSV、热源文件、物性文件、关键结果和日志。

至少记录：几何半径/高度、粒子间距、初始温度生成方法、总热功率、浮力系数、密度、各面热边界、环境温度、辐射参数、黏性模式、时间步设置、粒子排序周期及代码提交号。

禁止覆盖旧结果后再比较摘要。禁止同一运行编号下替换输入但保留旧图。每次运行都要有唯一 run_id，并在 CSV、图注和 summary 中一致使用。

## 3. 任务 A：整理博士论文物性数据并独立接入

### A1. 先找到真实数据，不能再重画近似曲线替代

查找第239轮提取包、对应 solver tables、原始图和提取脚本。优先使用已经提取、与目标玻璃材料对应的 Jacoutot 数据；其他博士论文材料只能作为明确标注的补充，不得自动混用。

对 μ(T)、k(T)、cp(T)、σ(T) 分别建立清单，包含：

- 文献名称、页码、图号或公式号、玻璃材料标识。
- 原始 CSV 的真实路径、SHA-256、节点数、有效温度范围。
- 横纵坐标单位、对数坐标转换、摄氏度到开尔文转换是否已完成。
- 数字化数据的不确定性，曲线是否存在截断、平台、人工补点。
- 插值方式、超范围处理、上下限、是否做锚点缩放。

**缺少文件必须报错并列明缺失项，不得悄悄回退旧表。** 不得在找不到原数据时重新手填几组数，并继续标记为“博士论文数字化数据”。

### A2. 数据集必须分开命名

建议至少保留以下明确身份；可以沿用项目命名，但含义必须等价：

| 数据集 | 用途 | 标注要求 |
| --- | --- | --- |
| `LEGACY_APPROX_TABLE` | 冻结旧基线 | 旧近似表，不标记为新博士数据 |
| `JAC_THESIS_RAW_BORROWED` | 本轮真实数据接入主线 | 博士原始提取曲线，借用于期刊复现 |
| `JAC_THESIS_SHAPE_CEP_ANCHORED` | 可选、需单独运行 | 若保留曲线形状并按期刊锚点缩放，必须明确写出变换 |

原始博士曲线与期刊表格在某温度点存在差异时，先保留差异。不得把缩放后的曲线再称为博士原始曲线。

σ 的“论文印刷公式”与“图线数字化/拟合”必须保留为不同模型。两者可以做独立测试，但禁止混成一条曲线，禁止默认哪条更接近最终图就更真实。拟合适用区间之外的结果要标注外推。

### A3. 物性选择与工况选择解耦

实现独立的物性数据集选择接口。例如：`material_dataset`、`material_manifest`，具体参数名称由现有代码结构决定。

加载物性只允许改变明确选择的物性函数。不得隐式改变：50 kW/55 kW、β、ρ、h、辐射率、几何、初始温度、是否 Marangoni、是否施加 Lorentz 力。

保留原有整套 case identity 接口也可以，但新的单独物性入口必须可用。对同条件新旧物性对照，自动比较 effective 配置；除物性及其必要派生时间步外出现其他差异，应停止该对照并说明。

### A4. 输出可核查的加载证据

每次运行输出 `material_effective.json`、`material_probe.csv` 和四条曲线图。

probe 至少包含 750、1000、1200、1400、1473、1500 K，以及每张表的最低/最高节点；超范围点仍可用于测试，但必须标记插值、外推或限幅。

CSV 至少包含：dataset_id、物性名称、T_K、原始表插值值、实际返回值、单位、处理方式、是否命中限幅。曲线图同时画旧表与新表，并标出模拟实际覆盖温度范围。

启动日志打印真实解析后的路径、哈希、节点数、范围、插值和外推规则。只打印一个模型名字，不足以验收。

### A5. 通过条件

原始文件可定位、数据清单完整、独立加载成功、探针数值可复算、缺文件测试会明确失败、物性切换不改变其他工况参数。全部满足后，才进入耦合短测。

## 4. 任务 B：确认实际算子使用了这些物性

逐项列出当前实际消费位置，不得仅检查初始化函数：

| 物性/量 | 必查消费位置 | 关键风险 |
| --- | --- | --- |
| μ(T) | 隐式或显式黏性相互作用及壁面黏性项 | 只更新输出数组，算子仍用常数；显隐式重复施加 |
| k(T) | 热扩散离散及接触传热 | 更新后未同步；仍使用旧均值 |
| cp(T) | 热扩散、体积热源升温、边界热损失升温 | 不同热项使用不同 cp，破坏能量闭合 |
| σ(T) | 电磁组装、Joule 热及更新触发 | 只在开始读一次；功率控制掩盖分布差异 |
| ρ、体积、质量 | 热量与温升转换、黏性权重 | 配置与实际粒子字段不一致 |

输出同一组稳定粒子 ID 在物性更新前后和算子执行前的 T、μ、k、cp、σ。选择冷区、过渡区、热区，验证其对应探针函数。

核查 CPU/GPU 同步方向和时机。若黏性暂时走 CPU 隐式，允许准确优先，但要明确时间成本。不得为了宣称“全 GPU”保留有问题的路径。

变 cp 的能量统计采用与实际 cp 模型一致的比焓：`h(T)-h(Tref)=积分 cp(T)dT`。不能继续用固定 cp×温差作为唯一能量依据。至少完成短时绝热热源升温检查或等价的非重复核查。

若本轮自然对流仍使用冻结热源：明确写“σ 尚未反馈该热源”，不能把物性接口接好等同于完成温度—电磁耦合。先完成热流主线，随后安排两种 σ 模型的同一温度快照电磁比较。

## 5. 任务 C：粒子排序与边界量回归

### C1. 对不同类型变量采用正确处理

随粒子输运的状态量应随排序保持 ID 对应；由当前位置定义的边界几何掩码、面法向、面积权重等，必须明确是排序携带还是排序/移动后重新计算。不能机械地把所有变量都注册为 evolving 后就视为结束。

检查热边界掩码、面编号、面积权重、外推距离、T、μ、k、cp、σ、Q 及相关缓存。固定空间热源在粒子移动后也需要正确重新映射，不能误当成永远随粒子移动的材料属性。

### C2. 必须保留的回归数据

1. 按稳定 ID 对齐的排序前后字段比较。
2. 顶面、侧面、底面权重在空间中的分布图与高度范围。
3. 各面有效面积与功率；既查总量，也查位置。
4. 表面温度的最小值、最大值、非有限值和非正绝对温度计数。
5. 明显粒子运动后至少再取一个时刻复查，而非只看初始状态。

顶部外推、底部外推与侧壁直接取粒子温度的处理差异要记录。当前边界壳厚度与粒子间距一起输出；不能把壳内粒子温度直接称作真实外表面温度而不解释。

若发现侧壁热损失偏高，先检查表面温度取值和边界面积，再考虑热阻模型。不得用新增壳层热阻掩盖坐标错位或重复计算冷玻璃层的问题。

## 6. 任务 D：最小测试矩阵，先隔离原因，再延长时间

| 编号 | 内容 | 目的/限制 |
| --- | --- | --- |
| R0 | 冻结上一轮旧物性结果 | 能复用则复用，不重复耗时计算 |
| R1 | 同几何、同初值、同边界、同功率，仅切换新物性 | 测量真实物性改变的影响 |
| R2 | R1 基础上，将实际控制隐式推进的时间步减半 | 检查时间离散敏感性，不能只改未生效的参数 |
| R3 | 通过健康检查后推进到 10 s、30 s，分段输出 | 观察流动和温度演化，不预设已稳态 |
| R4 | 仅在趋势和成本支持时继续一个窗口 | 需写清继续能回答什么问题 |

建议短测先用与既有结果可比的时间窗，例如 0.2 s 和 2 s；时间步减半比较使用相同物理终止时刻。

R1/R2 至少比较 T 场、速度场、Umax、Urms、各面热损失、焓变化和能量残差。报告相对差值及绝对差值；近零量不得只报相对误差。可用主要非零积分量差异 5% 作为本轮初筛提示线，但该值不是论文验收标准；若超过，应定位后再决定是否延长。

核查隐式子步实际累计时间、壁面项、局部 μ 配对方式，以及显式黏性 ForcePrior 是否重复保留。稳定不代表准确。

上一轮约 4.4% 的能量残差要重新按一致焓定义审查，列出贡献与时间积分口径。不得仅提高允许阈值使其通过。若只是统计口径错误，要用同一数据重算并说明；若为更新错误，先修复。

### 初始温度的身份必须写清

Fig.8 数字化重建是部分参考初猜，不是作者原始完整温度场。若仍用 `T=Tscale*Tstar` 且 `Tscale=1500 K` 为假设，要明确写出。若最大 Tstar 只有 0.95，则实际种子 Tmax 是 1425 K，而不是 1500 K。

禁止把该场偷偷做成仿射变换、整体加热或高温截断后，继续称为原 Fig.8 场。缺失区域、插值区域、外推区域必须有 mask 和覆盖率。主线耦合计算中不得冻结温度以维持目标环流；冻结温度仅可单独标注为诊断测试。

## 7. 任务 E：从下一轮开始强制生成论文同类图

不要求等待“算好”才出图。每个主要运行的初始、短时末帧及后续 10/30 s 关键帧，都要自动生成下述结果。保留原始数值和可重复脚本，不能只交截图。

### E1. Fig.8：温度等高线和温度曲线

至少生成：

1. 期刊 Fig.8 原图，保留图号和坐标定义。
2. 同一轴向剖面的计算绝对温度图 T/K。
3. 按论文定义归一化的 T* 图，使用匹配的等值线等级。
4. 参考等值线与计算等值线叠加图。
5. 至少三条轴向线和三条径向线的温度曲线；线的位置在剖面上标出。

曲线 CSV 至少记录 run_id、时间、r、z、归一化坐标、T、T*、参考值是否真实存在、参考来源和不确定性。

**参考只有离散等值线交点时，只画真实参考点。** 不得把自己补出来的连续场当成作者完整温度曲线。可以画数字化重建曲线，但必须与原始交点分别标注。

绝对温度与形状比较分开：Tmax 未完全确定时，仍可比较 T* 等值线位置；不能因此宣称绝对温度已对齐。原文按最大数值归一化的定义要核实并记录，不得默认使用 `(T-Tmin)/(Tmax-Tmin)` 替代 `T/Tmax`。

### E2. Fig.7：速度矢量、流线与环流结构

至少生成轴对称 r-z 视图下：速度大小、稀疏矢量、流线，以及关键径向/高度位置的速度剖面。使用一致的视向、轴方向、容器边界和尺度。

检查轴区、内部上升带、最外侧下降带及回流闭合。不能只用“轴区平均为负、壁区平均为负”就判定没有内部上升带，也不能把左右镜像显示的同一个环流算成两个独立环流。

如果由速度重建流函数，必须写清轴对称关系及采用的是物理速度而非 transport velocity。常用约定为 `u_r=-(1/r)∂ψ/∂z`、`u_z=(1/r)∂ψ/∂r`；其他符号约定可以使用，但须统一。检查轴线处理和离散一致性。流函数重建不能覆盖原始速度的散度或壁面法向残差。

对流线同时输出原始采样点矢量图；注明插值、平滑、掩码和种子规则。没有原始数据支持的区域不画流线。统计剖面重建的连续性误差和壁面法向速度，用于识别后处理误差；这不等同于直接断言粒子穿墙。

至少同时展示一个瞬时场和一个明确时间窗口的平均场。只有单帧出现的次级涡，标记为瞬时特征。已有 Fig.7 方向数字化点可作定性约束，不能伪造原文未给出的速度幅值。

### E3. Fig.6：Joule 热分布

保存原始电磁热源与功率归一化后的热源，分别记录体积分功率、Qmax、位置、轴向/径向分布。明确当前热源是否冻结、σ 是否更新、导体自身感应响应是否启用。

图上标出参考与计算使用的归一化方式。固定总功率后分布仍可能错误，不能用总功率正确替代空间分布验证。

如果本轮没有完成 Fig.6 数字化，先做清楚标注的定性并排比较，不能阻塞 Fig.7/8，也不能虚构量化参考。

### E4. 热收支与演化图

自动生成 Tmean/Tmin/Tmax、Umax/Urms、输入功率、侧壁/底部/顶部对流/辐射功率、累计输入能量、累计损失能量、焓变化及残差曲线。

稳态需用温度变化、热收支与流场变化共同判断。文献耦合迭代次数不能换算为同样的物理秒数；不使用“两秒没有目标图，所以路线不行”作为结论。

### E5. 比较报告必须有具体差距

建立 `paper_comparison_metrics.csv`，每行包含：比较对象、参考来源、参考数值/形态、计算值、误差定义、误差值、数字化不确定性、通过状态及原因。

可量化项目优先包括：真实参考等值线的归一化位置误差、真实剖面交点温度差、热区位置、上升/下降区域位置、已知速度幅值和 Q 分布指标。无参考数值的项目标记“不具备定量参考”，不得填 0 或 PASS。

不同运行使用固定色标；若为展示局部差异另用自动色标，需另存一张并明确标注。所有图注明 run_id、时间、物性集、温度单位、归一化定义和假设。

## 8. 电导率模型的有限比较

完成热流主线数据接入后，选择同一个有完整记录的温度快照，分别用 σ 印刷公式和 σ 图线模型求电磁场。先固定其他条件，报告原始功率、功率控制后的分布、等效驱动幅值和热源位置变化。

若两种模型都归一化到 50 kW，必须保留归一化前结果，否则可能掩盖两模型实际耦合强度差异。没有完成真实反馈更新时，不称为全耦合收敛。

这一比较不能与同时更换热边界、初温、几何等混在一起。不要本轮一次铺开全部组合。依据材料出处、实现正确性和对结果的影响决定后续，而不是仅看哪张图更漂亮。

## 9. 搅拌支线的时间限制

保留静止桨 p-prime 静水平衡回归。已知驱动搅拌仍存在显著速度峰值问题，终点速度正常不能掩盖过程超调。

若安排修复，只针对首次异常出现的粒子、位置、时间和受力项做一个有证据的修复；不再铺开无归因的密度处理组合。刚体规定运动与真弹性桨分别命名，强制重设密度的测试不能作为连续性准确的证据。

本轮法国自然对流任务未完成前，不扩展新的弹性算法和 GUI 功能。

## 10. 执行顺序及停止条件

按下列顺序执行，不能先跑长时计算再补物性证明：

1. 冻结旧配置与结果；找到并核对第239轮实际数据。
2. 独立物性接口、manifest、探针和新旧曲线图。
3. 确认热、黏性等实际算子消费，以及排序/更新正确。
4. 新物性短测，与旧结果同条件比较；立即生成第一批论文对照图。
5. 实际时间步减半检查、能量审查。
6. 健康后推进 10/30 s，持续生成对照图和差距表。
7. 有预算再做同温度快照两种 σ 的电磁比较及一个搅拌定点回归。

出现非有限量、非正绝对温度、边界权重错位、缺失数据静默回退、显隐式重复黏性或明显不合理能量增长，应停止延长并定位。仅仅图形没有双环流不是立即停止的条件；应结合温度、速度、热收支和演化判断。

不得通过关闭物性温变、冻结温度、重置密度、提高错误阈值来把主线标记为通过。必要诊断可以做，但须独立命名并说明不属于最终预测。

## 11. 交付目录与证据要求

建议目录结构如下，允许按项目惯例调整，但内容不得缺失：

```text
doc245_runs/
  00_baseline_freeze/
  01_material_dataset/
  02_operator_consumption/
  03_sort_boundary_regression/
  04_same_condition_comparison/
  05_timestep_check/
  06_main_evolution/
  07_paper_figures/
  08_stirring_limited_regression/
  09_handoff/
```

交付 `INPUTS.md`、`RUN_COMMANDS.md`、`DECISION.md`、`KNOWN_LIMITS.md`；附真实数据集、哈希、effective 配置、原始日志、关键帧字段、对比 CSV、PNG/PDF 和可执行绘图脚本。

绘图脚本要能从已落盘字段复现全部主要图，不接受只有空函数或依赖未交接文件的占位脚本。若精简打包排除大文件，必须给出其实际保存位置、哈希和重建方式；不能只保留 summary。

## 12. 最终总结必须逐条回答

1. 本轮究竟用了哪几个真实 CSV，路径、哈希、节点数和温度范围是什么？
2. 哪些仍是旧近似、借用、缩放或外推？是否存在静默回退？
3. μ、k、cp、σ 各自已进入哪些实际算子，哪些尚未耦合？
4. 更换物性时，有没有同时改变其他工况条件？
5. 粒子排序与热边界错位是否完成回归，运动后是否仍正确？
6. 隐式时间步减半后结果变化多少，能量残差如何解释？
7. 新旧物性使温度、散热和流动发生了哪些可量化变化？
8. Fig.6/7/8 分别已能比较什么，具体差距是什么，参考缺失是什么？
9. 当前只是瞬时演化，还是具备稳态证据？下一段时间推进能回答什么？
10. 哪个版本可继续交给软件分支，必须保留哪些准确性限制？

**本轮完成的标志是：物性数据实际生效、关键数值错误有回归证据、与论文同类的图和差距表已经生成。不是把任务书中的模型名称再写进 summary 一遍。**

<!-- ===== END 245_OPHELIE_THESIS_PROPERTIES_AND_PAPER_COMPARISON_PLAN.md ===== -->


<!-- ===== BEGIN 246_OPHELIE_ORDERED_CALIBRATION_AND_STIRRING_METRICS_PLAN.md ===== -->

# 246_OPHELIE_ORDERED_CALIBRATION_AND_STIRRING_METRICS_PLAN.md

# 第246轮：法国自然对流有序参数组合测试与搅拌指标准备

日期：2026-11-25。执行对象：物理开发电脑上的Cursor。

## 0. 本轮目标与边界

用有限、可比较的试验找到更接近法国自然对流Fig.7/8的温度分层和主次回流。允许校准吸收功率、有效热阻及已明确标注的热源模型。目标是得到真实计算产生、可以追溯的校准结果；改变文献条件的结果标记为“参考法国文献的校准工况”，保留实际参数，不沿用“原论文50 kW工况验证通过”的标签。

本轮不改浮力公式、不引入新的SPH框架、不重做方腔、不全面重写流体求解器、不把目标速度或目标温度持续施加到最终演化中。A-phi/Laplace独立分支和软件稳定交付分支不受影响。所有扫描限于OPHELIE研究case。

必做：修正环流判断和焓基线，完成共同基线与3个核心边界组合。选做：一个侧壁组合、一个热源位置组合、两个功率组合。筛选最多8个唯一配置；不展开全因子乘积。最终只确认最优候选和基线。

## 1. 当前证据与必须纠正的判断

第245轮R10c已出现轴侧下降、内部上升、最外侧下降及上外侧小回流候选；但温度范围约1396～1546 K，Tmin/Tmax≈0.903，与论文大范围0.5～0.8等温区差别明显。R11b均匀1473 K起算的20～30 s平均完整速度场也有内部上升和外围下降。不能因宽外缘平均Uz为负就认定没有环流。

论文主涡视觉区域在上部，当前R10c主涡偏右偏低。当前结果尚为瞬态。论文Fig.7方向样本和涡心ROI均是人工近似，没有可靠速度幅值，不得据此生成速度误差小于10%的结论。

现有skull/crust降低总散热，但未解析的冷层可能在流体图之外；绘制对比时必须标清流体域与等效冷层。不能通过补造边界冷温度填进流体图来得到论文等值线。

## 2. 先做两项短修复，不能再次延伸成新研究

### 2.1 环流诊断

保留旧Uz_axis/Uz_rim，但取消其一票否决。自动输出z/H=0.50、0.75、0.85处Uz(r)，报告内部上升峰位置、峰值、最外带下降速度、轴侧方向。用二维矢量图和流线确认回流，流线不得通过强制无散度修正原速度后冒充原解。

### 2.2 独立初始焓

第245轮R10c：独立焓残差约0.614 MJ，而实际热初猜均温比默认1473 K高约3.95 K，按97.21 kg与cp≈1599换算约0.614 MJ。优先检查诊断是否错误地以均匀1473 K为初始基线。

在物性加载、初温加载和静水准备完成、开始计入热源与散热之前，捕获实际H_initial=sum(m_i*h(T_i))，与后续同一h(T)定义和积分口径比较。若静水准备本身推进热方程，明确它计入哪段时间，不遗漏能量。

用一次短测确认。若只是统计基线问题，只修统计，不改热算子；若仍不闭合，列出剩余绝对残差及具体项，停止新增组合扩张，保留已生成结果。禁止提高阈值强行通过。

## 3. 公共配置：所有组合只改矩阵中明确列出的项

沿用第245轮R10c公共配置：R=0.25 m、H=0.185 m、dp=0.015 m；rho=2750 kg/m³、beta=1e-5 K^-1；JAC_THESIS_RAW_BORROWED的局部mu/k/cp；当前sigma解析模型；implicit_split；A_ind=off；顶部固定平坦滑移工程模型；周期EM更新0.5 s；默认固定玻璃吸收功率50 kW。

快速筛选统一使用同一份第245轮coldwall_skin ASSUMED初温CSV，文件哈希相同，初始速度一致。它是热启动加速初猜，不是参考真值，不改变其幅值、热团位置或填充方法以优待某个组合。freeze_prescribed_t必须为0。不能一组从均匀初温、另一组从热初猜直接排名。

确认阶段必须加入均匀1473 K初温对照，判断候选效果是否只来自预设热团。均匀1473 K是高温均匀启动，不叫室温冷启动。

先核对当前源码中的有效参数。以下5 mm和k=1来自已审查历史实现，若本机实际R10c不同，先报告并把基线改成实际值；不得静默把错误基线继续称为原R10c。

### 3.1 显式基线热参数

| 参数 | 基线值 |
| --- | ---: |
| 冷却侧h_side | 300 W/(m² K) |
| 冷却侧h_bottom | 35 W/(m² K) |
| 顶面对流h_free | 20 W/(m² K) |
| 顶面辐射率epsilon | 0.8 |
| Tcool、Tamb、Trad | 各300 K |
| 侧壁附加层delta_side | 0.005 m |
| 底部附加层delta_bottom | 0.005 m |
| 顶部附加层delta_top | 0.005 m |
| 附加层k_side/bottom/top | 各1 W/(m K) |

### 3.2 必须看有效热阻，不只看厚度

侧面/底部的串联模型为R_total=delta/k+1/h_cool；h_eff=1/R_total。基线h_side_eff=120，h_bottom_eff≈29.787 W/(m² K)。底部厚度5 mm减至2.5 mm，在h_bottom=35下仅使h_eff增至约32.184，增强约8%，不能误称“冷却翻倍”。

因此本轮底部主敏感性方案同时设h_cool_bottom=70和delta_bottom=2.5 mm，有效h_bottom≈59.574，约为基线2倍。它是一个“底部有效冷却组合”，不是分辨h与delta独立贡献的试验。70在论文讨论过的底部敏感性范围内，但不能声称就是Fig.8的确定参数。

顶部包含非线性辐射，不能用相同h_eff公式替代实际表面温度平衡。必须输出Ti、Ts、对流/辐射/导热热流及其闭合。

## 4. 阶段A：按顺序测试边界组合

全部保留50 kW和原始EM热源。不改sigma、不移动线圈、不改mu/beta。具体值见同包CSV。

| ID | 侧壁delta / mm | 底部delta / mm | 底部h_cool | 顶部delta / mm | 意图 |
| --- | ---: | ---: | ---: | ---: | --- |
| C0 | 5 | 5 | 35 | 5 | 显式基线 |
| C1_BOTTOM | 5 | 2.5 | 70 | 5 | 单独增强底部有效冷却 |
| C2_TOP | 5 | 5 | 35 | 2.5 | 单独增强表面冷却 |
| C3_BOTTOM_TOP | 5 | 2.5 | 70 | 2.5 | 同时发展底部冷区和表面冷玻璃 |
| C4_SIDE_BOTTOM_TOP，条件执行 | 2.5 | 2.5 | 70 | 2.5 | 对照C3再增强侧壁冷却，检验外围下降带 |

C4侧壁h_cool保持300，h_side_eff由120变至约171.429。仅当C3仍稳定、最外侧下降带不足或近壁冷区仍明显不足时执行；如果C3已整体快速过冷而上部热区消失，不执行C4。

不能预先指定C3必胜。它可能使整锅冷却，也可能使流速降低；都应如实保留。C1、C2用于区分底部与顶部对温度分层及回流的贡献。

### 4.1 已存在接口与最小适配

历史源码已有：`--skull-resistance`、`--skull-delta-side=`、`--skull-delta-bottom=`、`--skull-k=`、`--skull-h-cool-side=`、`--skull-h-cool-bottom=`、`--free-crust-resistance`、`--free-crust-delta=`、`--free-crust-k=`。必须检查本机--help和实际解析后再使用，不能把未知参数忽略后继续跑。

统一参数顺序：公共参数先、组合覆盖后；禁止脚本后部重复添加默认值覆盖组合。每组启动时打印原始h、delta、k和最终h_eff。明确thermal BC覆盖与skull转换的调用顺序，不重复串联两次。不能把h_eff再次当作h_cool送入同一转换。

## 5. 时间与计算预算

先用C0进行5 s健康检查和墙钟测量。通过后，阶段A统一推进到60 s，保留10、30、60 s的字段和热收支。不得以60 s为准稳态验收。

可以对所有组合使用同一份完整的共同检查点，随后在同一物理时刻切换参数做阶跃敏感性；但这一批必须全部同样处理，并与从初始条件起算的结果分开编号。不能只给某一组继承已发展好的速度场。

没有可靠重启动功能时，直接从同一初温重新跑；不要为扫描重写checkpoint系统。如果已有运行无法保存逐帧字段，至少输出关键帧原始粒子或网格场，而不是只留末态summary。

阶段A先完成C0～C3，四组足以决定主要方向。选择1～2个候选和C0继续到120 s；不对全部组合无条件长跑。若墙钟明显超预算，先停止C4、热源诊断和功率扩展，保留核心四组及图表，不缩短到彼此不同物理时间后强行比较。

## 6. 候选评价：先健康，再温度，最后回流

### 6.1 健康条件

无NaN、负绝对温度、明显粒子逃逸；步长实际生效；能量计算使用实际初始焓；热源转移、边界面积和排序量空间归属合理。健康不等于文献验证通过。

### 6.2 共同评价量

1. T/Tmax等高线、真实参考交点的温度差、Tmin/Tmax及冷区体积比例。
2. 上部热区是否保留，底部/轴侧/近壁的温度是否按目标方向分化；不能只优化Tmin一个粒子。
3. 主回流位置、内部上升带位置与最外侧下降带；通过平均速度场与原始矢量相互核对。
4. P输入、各面散热、焓变化率及温度变化率。
5. 10～30、30～60、60～120 s的变化趋势，是否只在某个瞬间相似。

温度指标采用真实参考交点和同一T/Tmax定义。不对缺失区域最近点外推填0/PASS；参考位置不在数据支持域内应标NA。参考图数字化的坐标与角度有不确定性，不能机械地把10%设为所有形态指标的合格线。

不定义一个把温度、方向和热平衡随意加权的总分掩盖失败项。先淘汰不健康结果，再按温度差距、回流位置、变化趋势列出并排表，保留折中候选。如果没有候选改善，明确报告，不自动增加十几个新参数。

### 6.3 必须统一的绘图口径

每组同一时间输出Fig.8式温度图、绝对温度图、两条参考剖面、Fig.7式平均速度图、至少三条Uz(r)及热收支。统一坐标、参考等值线、归一化和色标。不得采用(T−Tmin)/(Tmax−Tmin)替代T/Tmax制造冷热对比。

## 7. 阶段B：一个上移热源的诊断组合

仅在阶段A有健康候选且主热区/主涡仍偏低时执行。继承该候选边界，设Q1_UPPER：对EM原热源施加非负空间权重

`Q_controlled(x)=C*Q_EM(x)*exp[1.5*(z/H−0.5)]`。

z是相对于玻璃底面的高度，不能直接用世界坐标。C由实际体积积分求得，使本轮规定吸收功率保持不变。顶部相对底部的权重比约4.48，属于有明显区分度的一次诊断，不是已知论文热源。

这是“EM形状基础上的受控热源”，不是原始Maxwell预测；单独编号、单独图注。它只能回答热源上移能否改善温度和涡心，不能用来声称EM验证通过。

实现时保留Q_EM原字段，新增Q_controlled字段，不能覆盖原始EM输出后仍标原场。每次EM刷新后仅加权一次，保存原/新Q分布、体积分、功率重心z_Q、径向重心和上部功率占比。推荐在固定EM粒子/欧拉场上按体积加权归一，再传给流体；粒子采样误差另记录，不把目标功率当作已施加能量。

接口若不存在只新增局部受控热源入口，默认关闭；不得修改原EM算子、线圈几何或其他case。不得通过直接上移初始温度团代替本试验。

若只改变Q就显著改善涡心位置，后续再通过sigma模型或电磁配置寻求相近原始Q；若无明显改善，停止加大加权系数。

## 8. 阶段C：有限功率比较

选一个最好的边界/热源身份，固定其余参数。50 kW结果已存在，只新增：P35=35 kW、P65=65 kW。两者为校准搜索点，不是论文参数范围或实验数据。

不把“当前散热62 kW”当成必须把功率改到62 kW的依据。增加功率可能让整锅更热而分层更差；降低功率也可能只导致整体降温。根据温度形状与热收支共同选择。

采用原始EM模式时通过固定吸收功率控制重新求EM；采用Q1诊断时按目标功率归一Q_controlled，必须保留其受控身份。电源输入功率与玻璃吸收功率分开，不强行设效率。

本轮不再同时扫描beta、黏度倍数、cp、k、sigma和功率。若最佳候选仍只欠流速，下一轮另做一个明确的流变对照，不能在本轮暗改已选物性。

## 9. 最终确认与展示结论

最优候选至少延长一个共同观察窗口，展示后期平均而非仅挑一帧。补均匀1473 K初温对照；若仅热初猜短窗能相似，明确标注瞬态热启动结果，不能称独立稳态复现。

最终表同时列：论文原参数、计算实际参数、哪些被校准、哪个量用于校准、哪些量作额外检查。未达到稳态就报告时刻和演化，不把有限终止时间自动视为平衡。

保留未经修改的R10c/R11b基线。PPT可以展示形态改善，但功率、边界、初态与受控热源身份应可查；不覆盖不理想结果、不把改后的功率写成论文50 kW、不修图改数据。

## 10. 搅拌：先明确要比什么，不扩大到弹性桨研发

法国机械搅拌论文以功率数Np和Nu等宏观指标评价系统，相关Re、Pr及等效黏度定义需从期刊正文核对；博士论文用于补充。先建一张reference table：页码/公式、几何、N、吸收功率、物性定义、特征长度、参考值或关联式及数字化误差。

不同几何的甲方case可计算这些指标，但不能直接拿它的值当法国原几何验证。不要凭记忆确定桨直径或用坩埚直径代替；先核对论文及实际STL、轴向和长度单位。

### 10.1 Re与Pr：计算容易，口径优先

常见搅拌定义Re=rho*N*D²/mu，Pr=mu*cp/k，最终以原文为准。N用转每秒：10 rpm=1/6 s^-1，omega=2*pi*N约1.0472 rad/s。动态黏度与运动黏度不可混用。

输出物性所取温度、局部/平均/等效黏度的定义。论文等效黏度不能默认换成全域均值。若原文定义暂未取得，只输出“工程定义指标，未完成论文同口径比较”。

### 10.2 Np：公式容易，可靠力矩是前提

常见Np=P_mech/(rho*N³*D⁵)，P_mech=omega*M，最终核对原文。M由流体对桨的压力和黏性载荷，绕实际搅拌轴积分；轴原点、单位轴向量、力作用位置、作用反作用符号要明确。

单独输出M_pressure、M_viscous、M_total、P_mech与时间序列，不把驱动惯性力矩、固体内力、离心载荷或电磁吸收热功率当搅拌阻力功率。双边相互作用不得重复计数。

法国论文对照可以采用规定刚性转动桨；不需要为得到Np先解决甲方弹性桨全部FSI问题。若已有旋转流体仍有速度尖峰，先定位该运行质量，不用计算出一个数就宣布Np已验证。扭矩应在桨相位/多窗或完整周期范围报告，短窗不足时明确覆盖比例。

若等效黏度由同一扭矩反推，再用它验证Np-Re关联，可能只是代数恒等关系；这种一致性不是独立验证，必须说明。

### 10.3 Nu：最可能需要新增物理/输出的部分

先核实论文Nu评价的是哪一个换热表面、特征长度与温差。若其对象是冷却搅拌桨，而当前桨绝热，仅有流体温度和力矩就不能得到同口径Nu。

常见工程形式Nu=h*L/k，h=Qdot/(A*DeltaT)，但L、A、DeltaT及参考k需按原文。所需Qdot是流体向指定换热面的热量，不是全域EM热功率。输出面热流、面积、参考温度、壁温及其时间平均；避免温差近零导致虚大Nu。

若流固接触传热已存在，只补面热流积分和指标输出；若不存在，必须新增真实热边界/传热并做能量核对，不能仅补一个Nu后处理公式。数据不足标NA，不用规定h反算Nu再当作独立预测验证。

### 10.4 本轮搅拌的有限交付

先交reference table、Re/Pr计算表、现有力矩和热流输出清单。没有可靠力矩则列Np缺口；没有论文对应热边界则列Nu缺口。自然对流核心四组未完成前，不扩展搅拌参数扫描或新弹性算法。能后处理完成的工作可穿插，但同一GPU重任务顺序运行，避免抢占使所有结果延迟。

## 11. 自动化、记录与交接

用本包CSV生成运行清单。每组独立目录，保存config_requested、config_effective、命令、代码hash、初温hash、物性表hash和原始日志。写一个顺序执行脚本即可，不要搭建复杂优化平台。失败组只记录失败，不把它替换为默认参数重跑而继续沿用原编号。

配置覆盖后先dry-run打印并核对h_eff、Q模型、功率；运行后summary再次写出实际值。建立总表：ID、父配置、变化项、终止时间、健康状态、温度参考误差、上升带、涡心、热收支和是否继续。

打包至少包含：新增/修改源码、可运行命令、CSV矩阵、effective、关键帧T/U/Q和边界面温度、物性证据、能量初始基线、全部同尺度对比图、评价表和搅拌指标缺口表。只交summary不算完成。

最终回答：哪个组合改善了哪个量、代价是什么、是否依赖热初猜、是否仍净冷却、哪些参数偏离原论文、下一步只需补哪个关键点。不得把所有组合统一写成PASS。

<!-- ===== END 246_OPHELIE_ORDERED_CALIBRATION_AND_STIRRING_METRICS_PLAN.md ===== -->


<!-- ===== BEGIN 247_CURSOR_EXECUTION_PLAN.md ===== -->

# 247_CURSOR_EXECUTION_PLAN.md

# 247：法国搅拌论文数字化、参考场初始化与报告图输出

## 0. 目标和范围

本轮直接使用本包的数据，完成法国自然对流与法国搅拌两条参考工况的结果输出。甲方case、Linux软件及A-phi/Laplace研究分支不在本轮修改范围。保留当前已能出现环流的稳定设置，禁止重新叠加未经归因的压力、连续性、shifting或黏性方案。

本包提供：实际人工读图描点、插值得到的二维初猜、无量纲关系、Python计算辅助程序与原图叠加检查图。初猜可用于缩短启动过程；它不能替代实际运行的结果文件。报告图注可简洁写“以文献温度场重建值初始化，演化Δt=…，采用…热源”，无需长篇描述启动过程，但不能标为独立冷启动复现或稳态验证。

优先级：先让已有稳定代码正确接收物性和参考场，输出真实VTP/CSV；然后比较温度轮廓及搅拌指标。不要先重构求解器，不要为凑图修改导出的数值或色标数据。

## 1. 两篇论文的图号必须隔离

| 唯一标识 | 内容 | 本轮用途 |
|---|---|---|
| NC_Fig6 | 自然对流论文，归一化焦耳热分布 | 规定空间热源的形状参考；与实际EM预测热源分开 |
| NC_Fig8 | 自然对流论文，轴对称温度等值线 | 已有参考数据继续使用；不能把搅拌1620 K套过来 |
| NC_Fig10 | 自然对流论文，含测温探针的三维温度分布 | 若未建探针，不声称相同工况 |
| STIR_Fig6 | 搅拌论文，A/B监测点随时间的温度 | 不是焦耳热分布图 |
| STIR_Fig8a/b | 搅拌论文，两个竖直温度截面 | a主对比；b本包未数字化，不伪造数据 |
| STIR_Fig8c | 搅拌论文，穿过桨叶的水平温度截面 | 辅助对比，局部非轴对称明显 |
| STIR_Fig8d | 搅拌论文，自由液面温度 | 第二张推荐对比图，勿误设为桨叶高度 |
| STIR_Fig9 | 搅拌论文，对应四个面的归一化速度/流动图 | Vmax=0.22 m/s；不能把等值线直接当作流函数或数字化矢量场 |
| STIR_Fig10 | 搅拌论文，收敛后的焦耳热分布 | Qmax=6.64e6 W/m³；本轮仅核对，不包含此图数字化网格 |
| STIR_Fig11 | 等温搅拌Np–Re | 独立常物性验证的首选 |
| STIR_Fig15/16 | 冷却桨传热关系 | Fig16优先，Fig15完整Nu口径仍有不确定性 |

搅拌文献：Jacoutot et al., Chemical Engineering Science 63 (2008), 2391–2401，DOI 10.1016/j.ces.2008.01.026。以下页码为印刷页。

## 2. 本包实际提取了什么

`data/*_contour_points.csv`：图上可辨识等值线的人工稀疏描点，字段有曲线编号、点顺序、原图像素、归一化坐标、等值线数值。STIR温度提供 `T_K = 1620 × value_normalized`。这是视觉数字化数据，不是作者原始场数据。

- STIR_Fig8a：159点，9条折线片段，涵盖0.6/0.7/0.8/0.9。
- STIR_Fig8c：91点，9条片段，涵盖0.7/0.8/0.9。
- STIR_Fig8d：78点，5条片段，涵盖0.40/0.45/0.50/0.55，0.55存在内外两条轮廓。
- NC_Fig6：110点，5条片段，仅保留核查后的0.1/0.3/0.9；中间曲线未完整确认，未补成假观测。

`figures/*_trace_QA.png`：描点与原图叠加；Cursor必须先打开核对，再导入。

`data/*_INITIAL_GUESS_grid.csv`：161×121二维网格，散点线性插值、外区补全。它是**初始化重建**，不是可逐点打分的参考真值。

字段解释：

- `x_unit`：竖直图横向坐标由左至右0…1；NC图专指r/R。
- `z_or_y_unit`：竖直图底部0顶部1；水平图为第二平面坐标底部0顶部1。
- `reference_visible_mask`：原图可见区域掩码，黑色搅拌桨示意/遮挡区域已排除。**不是实际流体/固体几何掩码**，特别是自由液面图可能含投影示意；真实掩码必须由STL与粒子类型决定。
- `nearest_extrapolation`：原描点凸包之外，不能当作观测。
- `distance_to_trace_unit`：距最近描点的归一化距离，可用于标识低可信区域。
- `value_uncertainty_estimate`：STIR暂取±0.025、NC暂取±0.05，仅为读图工作估计，不是统计置信区间。坐标估计约±4像素。

二维网格各区域均含重建假设；边界极薄的低温层没有被完整解析。不要把重建网格上数万个点误报为数万个文献观测。插值值不会自动达到Tmax=1620 K：原图没有提供最高温点坐标，禁止为强行达峰对整个场再次归一化。

### 2.1 温度初猜局限

a/c/d分别有约16%/25%/28%的可见网格位于描点凸包外，采用最近值补全。若直接映射，冷壁附近薄层会偏暖；真实Robin边界应随后自行演化，不能把这个薄层称为已复现的凝固层。

Fig8a只有一个竖直剖面，不能确定唯一三维场。默认方案是以a构造轴对称近似初猜，后续三维桨运动产生非轴对称变化；应记录这个假设。c/d先作为对比图，不强行同时拼成连续三维场；尤其不能在不同高度的三个面之间简单硬覆盖。

## 3. 搅拌论文已确认工况

| 参数 | 数值/定义 | 来源 |
|---|---|---|
| rho0 | 2800 kg/m³ | p2392 §2.1 |
| beta | 3e-5 K^-1 | p2392 |
| Boussinesq T0 | 293 K，参考温度；不是启动温度 | p2392 |
| 物性 | μ、cp、k、σ随温度，见Fig2 | p2393 |
| 图8/9/10转速 | 10 rpm = 1/6 rev/s | p2397–2398图注 |
| 图8/9/10焦耳吸收功率 | 60 kW，玻璃内积分功率 | 同上 |
| 图8最高温度 | 1620 K，即1346.85℃ | Fig8 |
| 图9最大速度 | 0.22 m/s | Fig9 |
| 图10最大焦耳热 | 6.64e6 W/m³ | Fig10 |
| 图8/9/10状态 | 作者耦合收敛，t=1000 s | 图注及正文 |
| Fig12/13/14 | Ptot=56 kW，勿用60 kW标签 | 各图注 |
| Fig15/16 | 56与70 kW两组 | 各图注 |
| 自由液面hs | 20 W/(m² K) | p2394 Table1 |
| 水冷参考Tc | 333 K | Table1 |
| 对流环境Ta | 303 K | Table1 |
| 辐射环境Tar | 423 K | Table1 |
| 发射率 | 0.9 | Table1 |
| 非运动冷壁hw | 300 W/(m² K) | p2396 |
| 冷却搅拌桨hw | 2500 W/(m² K) | p2396 |
| 固壁速度 | 无滑移 | §2.1 |
| 自由液面速度 | 滑移 | §2.1 |
| 玻璃“凝固”判据 | μ>1e4 Pa·s | p2397–2398 |
| 流体参考系 | 作者采用旋转参考系 | §2.1 |
| EM更新 | 热步0.2 s，选10热步一次，即约2 s | p2396；SPH不得照搬其热/声学时间步 |

自由液面向外热损失：

`q_out = 0.9*5.67e-8*(T_surface^4 - 423^4) + 20*(T_surface - 303)`，T全部为K。

非运动冷壁：`q_out=300*(Twall-333)`；冷却桨：`q_out=2500*(Twall-333)`。论文只建玻璃域，hw合并多种热阻。若本程序显式建立实体导热和冷却，不要重复叠加同一层热阻。

注意：作者用高黏度形成的低流动层解释热阻。只设置同一个hw，尚不能保证我们得到相同冷壁传热。不得直接把论文阈值当成已实现的相变潜热/固化模型。

论文介绍提到650 mm坩埚、400 kW/300 kHz电源；不等于玻璃直径、电磁吸收功率或所有图的精确实验输入。当前自然对流几何中已确认的500 mm玻璃直径继续保留；搅拌几何及桨直径须从已交接模型/明确文献尺寸复核，不能从650 mm推算。**Re、Np中d为桨直径。**

### 3.1 温变物性落实顺序

1. 使用已交接239/245的曲线文件及配置，记录每条曲线来源，核对是否为本搅拌论文Fig2对应玻璃；不得仅因为都叫法国玻璃就自动认定完全一致。
2. 本包未重新数字化Fig2，不提供新拟合系数。已有曲线能接则先接；来源不一致要记差异，不用未知多项式替代。
3. 温度初始化之后，先计算实际μ(T)、cp(T)、k(T)、σ(T)，再第一次EM求解。不得先按默认σ=16求一次EM后才导入温度。
4. 定期从粒子温度投影到EM网格，记录σ目标、松弛后的实际σ、Q映射以及更新时刻。仅导出了σ(T)图不等于EM已经使用该值。
5. 每次钳位/外推输出数量、体积分数、阈值及来源。禁止无日志冻结μ或σ。

## 4. 自然对流Fig6如何用作热源

焦耳热是体积源Q，单位W/m³，不是温度初值，也没有“只在t=0加一次Q”的物理等价方案。

本轮可并列保留两种模式：

- `q_mode=em_coupled`：实际EM求解、随温度更新σ、映射Q；用于检查真实电磁耦合。
- `q_mode=prescribed_nc_fig6_shape`：用本包图6重建形状，在整个短窗内持续施加；用于检查给定非均匀热源时的温度与流动响应。禁止同一时刻叠加同功率EM热源。

图6只数字化到若干轮廓，低于0.1的大范围区域没有直接数值。当前初猜在凸包外用最近轮廓值随距离指数衰减（归一化衰减尺度0.12）；这是明确的低热源区补全先验，**不是论文公式**。约71%的二维域涉及这类外推，适合做受控热源测试，不适合宣称精确重建全域EM热源。

实施：

1. `r_unit=r/R`, `z_unit=(z-z_bottom)/H`，轴对称复制到真实玻璃粒子。
2. 用真实粒子体积 `V_i=m_i/rho_i`，计算 `S=sum(Qshape_i*V_i)`。
3. `Q_i=Qshape_i*Pabs/S`。自然对流论文该组总焦耳功率50 kW；若采用其他功率，明确记录。
4. 本包 `normalize_source.py` 提供此步骤。采用二维轴对称积分时应使用 `2*pi*r*dr*dz`，不能用普通二维面积求和。
5. 记录有效Pabs、峰值Q、热源重心、上半区功率占比。Fig6归一化图本身不能同时唯一确定原文Qmax与未知的低热源区。
6. 若只把规定源作为EM迭代的初始猜测，首个实际EM收敛解会覆盖它；需明确采用“持续规定源”还是“EM初猜”，不要混称。

可以配合文献温度场初猜运行，但报告名称应为“规定热源下热流动计算”，不能把规定Q图当作电磁模块独立预测结果。

## 5. 搅拌无量纲量：落实到程序

### 5.1 三种功率分开命名

- `P_Joule_W`：玻璃内电磁焦耳热积分，即Ptot，56/60/70 kW。
- `P_viscous_W`：机械耗散功率P，Fig11/12和Np使用。
- `P_heat_removed_stirrer_W`：冷却搅拌桨从玻璃吸收的热功率，即Fig16的Pstir。

禁止把Ptot代入Np，也禁止把机械搅拌功率当成Fig16的Pstir。

### 5.2 机械功率和Re/Np

`N = rpm / 60`，单位rev/s；`omega = 2*pi*N`，单位rad/s。

`D = 0.5*(grad(u) + grad(u).T)`。

`P_viscous = sum(2*mu_i*(D_i:D_i)*V_i)`。

这是论文Eq14/15。稳态时才可作为整体搅拌耗功的近似；短窗内还存在动能变化，不能强行令瞬时耗散等于桨输入功率。

扭矩交叉检查：`P_torque=2*pi*N*C_resisting`。需分清桨表面力、驱动反力及弹性储能；不能直接用电机设定功率。本文因冷却层附近扭矩梯度误差，后续优先采用耗散积分。

`Re = rho*N*d²/mu_ref`

`Np = P_viscous/(rho*N³*d⁵)`

常物性等温验证参考：`Np = 47/Re`（Fig11）。只有同一种桨几何、浸没深度及容器比例才能用47做定量门槛；甲方不同桨只能计算自身Np/Re，不能照抄47。

变量物性时应同时输出：体积平均μ、有效μ与范围，不能挑一种平均值恰好使NpRe=47后称为验证。

论文的等效定义：

`Re_e=47/Np`，`mu_e=rho*N*d²/Re_e = P_viscous/(47*N²*d³)`。

**按该定义必有Np*Re_e=47，所以这个乘积不是独立验证。** Fig11应另跑常物性等温工况：预先指定μ，计算P，再对比Np/Re。

低成本首轮：保持同一桨/网格，给定常物性，试5/10/20 rpm或等价Re点。先检查P和扭矩是否趋稳，再计算误差；受时间限制不能稳态时输出趋势和未完成标签，禁止直接PASS。

### 5.3 传热指标

`Pr_e=mu_e*cp_ref/k_ref`。cp_ref、k_ref的评价温度与平均方法必须记录。

Fig16直接给出：

`P_heat_removed_stirrer/P_Joule = 0.005*Re_e^(1/3)*Pr_e^(1/3)`。

先完成单个10 rpm工况的取量，不要求立刻跑出整条56/70 kW转速曲线。若时间允许，再逐步增加点。

Fig15：`Nu/Pr^(1/3)`与`Re^(1/3)`相关，系数随总焦耳功率变化，不能假设唯一K。正文给 `Vi=mu_wall/mu_bulk`；采用等效黏度时不再额外乘Vi修正。论文没有充分展开Nu的特征长度、面积和温差平均细节，本包辅助程序只输出明确标注的工程定义Nu，不把它冒充已经统一口径的Fig15数值。

图11/16 CSV由论文公式生成，是解析参考线，不是伪装成读出的散点。请保留CSV的source字段。

## 6. Cursor执行顺序与停机条件

### P0：封存、接入与零时间检查

1. 建立隔离case配置，复用当前稳定求解路径；列出本轮git差异。
2. 打开3张STIR叠加图和NC叠加图。不得再次让大模型从头猜读整图。
3. 为每个平面输出：物理坐标原点、方向、宽高、物理单位、旋转角度、与桨几何的对应。Fig8a/b不同切面不能互换。
4. 导入Fig8a：左半与右半分别转换到r/R，轴对称初猜默认取两半平均；仅为初始化，不要求计算结果保持左右镜像。插值到真实粒子；由STL决定实体排除区域。
5. d仅作为自由液面二维校验，不把它设为整个体积的温度。c为桨高度切面校验。
6. t=0立即导出T、μ、σ、k、cp与初始化元数据；先看是否坐标翻转、摄氏/开尔文混用、固体内误填流体。
7. 使用与物性模型一致的内能 `h(T)=integral(cp(T)dT)`做能量差；不能简单用变化cp的 `cp(T)*T`当内能。

### P1：短窗可运行性

1. 先运行0.02 s和0.1 s；通过后到0.5 s，再视稳定性到2 s。均以实际达到时间为准。
2. 使用已有稳定时间步与黏性处理。若低温区大黏度导致刚性，优先查隐式残差、壁面接触和时间步，不先改温度/物性曲线。
3. 温度/σ导入必须早于第一次EM求解；电磁场初始化与更新要有单独日志。
4. 检查NaN、负物性、粒子越界、真实达到时间、EM迭代残差、热源积分及能量项。
5. 每个检查点保存VTP和CSV，终止也必须保存最后有效状态。不得只落图或summary。
6. 粘性耗散进入热方程与否按当前模型明确记录，功率诊断积分本身不自动意味着热源重复施加。

### P2：同口径比较

温度图使用固定1620 K作为归一化分母；同时保留实际K与℃色标。不得用本次计算自身Tmax再归一化来掩盖绝对温差。

对原始轮廓点插值采样，逐线统计预测T/1620与参考值差；避开原图黑色遮挡，单独列出±读图不确定度。不能用初始化重建网格直接作为独立验证的几万点样本。

交付3列图：文献原图、计算同面温度等值线、差异/轮廓偏移图。保留t=0与后续时刻；若图形仅在t=0相近随后迅速离开，应报告漂移，不算稳定复现。

环流图用实际速度矢量和流线另画；不得用温度等值线替代流线。三维速度作平面投影后未必二维无散；不要强行把任意截面都解释为严格轴对称流函数。若需要参考Fig9，先核对旋转参考系：实验室与旋转系满足 `u_abs=u_rel+omega×r`；两种速度都可导出，图中明确使用哪一种。

定量门槛分层：

- 接入通过：单位/坐标正确，数据与实际使用物性一致，文件齐全。
- 运行通过：达到声明时长，无数值中止，热源功率一致。
- 参考初值保持性：报告温度轮廓误差及随时间漂移，不等于冷启动复现。
- 独立搅拌验证：常物性稳定工况的Np/Re及误差；10%只能用于事先约定的指标，不能事后删掉未通过项。
- 耦合复现：实际EM–热–流动演化得到的温度、速度和功率；不得用规定热源那一支的Q来签EM合格。

## 7. 报告展示字段：不要求与文献逐场比较

按用户要求，以下是实际计算字段的分布展示：

- 温度T_K与T_C。
- 实际用于动量方程的动态黏度mu_Pa_s（可用对数色标，注明范围）。
- 实际用于最近一次EM的sigma_used_S_m；另有sigma_target_from_T供内部诊断。
- 电磁B幅值、J幅值与Q_Joule，明确峰值/有效值及复数幅值约定。若J是峰值复幅值，时间平均Q=|J|²/(2σ)；若RMS则无1/2，必须全链一致。
- 速度m/s、桨形状、线圈实体外观。

每帧附带 `physical_time, em_last_update_time, q_mode, property_profile_id, N_rpm, Pabs_W`。规定Q场文件名不得叫EM预测结果。插值得到的背景网格也需保留覆盖掩码，避免容器外无效值被画成物性。

展示σ/μ分布不需要文献对比，也不能据此直接判定两个模块已独立验证通过。

## 8. 最终文件与本轮总结格式

至少产出：

1. `effective_config.json`：所有单位、物性来源、边界、几何、初始化和热源模式。
2. `plane_mapping.json`：归一化数据到物理坐标的映射。
3. `initialization_audit.csv`：T、物性、掩码、钳位范围与统计。
4. `history.csv`：实际时间、功率、温度范围、动能、EM更新时间和残差。
5. 原始结果VTP、同面CSV、温度/速度/物性PNG。
6. `stirring_metrics.json`：三个功率、d、N、Re/Np/Pr和定义。
7. `comparison.csv`：按参考轮廓逐线统计，不混入重建区。
8. `round_summary.md`：完成项、未完成项、明确失败原因、下一步只列最关键的1–2项。

验收展示优先：自然对流Fig8温度、其实际速度/环流；搅拌Fig8a与8d；温变物性与EM展示；搅拌Np/Re单点或少量点。时间不足则明确哪些是参考场初始化的短窗结果，哪些是独立验证。不要再次花整轮只做任务拆分而没有实际结果文件。

<!-- ===== END 247_CURSOR_EXECUTION_PLAN.md ===== -->


<!-- ===== BEGIN OPHELIE_Doc241_round_handoff_20261124_WORK_SUMMARY_ZH.md ===== -->

# OPHELIE_Doc241_round_handoff_20261124_WORK_SUMMARY_ZH.md

# Doc241 本轮做了什么 / 没做什么（打包说明）

日期：2026-11-24。分支：`feature/EM-Biot-Savart` @ `16d67c9d0` + **未提交 WIP**（见 `code/patches/`）。  
正式门控标签：**`LINUX_SOLVER_HANDOFF_READY`**

---

## 一、做了什么

### P0 冻结与诊断
- 冻结 HEAD / status / 二进制 hash / 设备到 `00_freeze/`
- 加强有限性 / 哨兵 / `Umax` 等假正常防护（F0 侧）

### P1 甲方弹性搅拌
- 根因：**粗 Reload** Vol→dp≈0.01，旧逻辑强制系统 dp=0.005 → 尺度不一致（负对照 exit 3）
- 修复：去掉无条件强制；从 Reload `Vol^(1/3)` 推断系统 dp；相对误差>15% → exit 3；写出 `reload_manifest.json`
- Path B：`--dp=0.005` 重生 dense Reload（n_rotor=4329）
- S2：0.002 / 0.02 / **0.05s PASS**（不锈钢弹性，GPU）

### P2 甲方自然对流（真实 EM+热）
- 修 `--client-nc-material` 清零吸收功率 → 保留 CLI `--target-absorbed-power`
- ASSUMED 玻璃吸收 **15 kW**（≠ 电源标称 30 kW）
- EM 1s：`P_joule=15000`，热推进，`em=RUN`
- β=0 对照：U_rms 比≈1 → 标签 **`OPEN_FS_RESIDUAL_DOMINANT`**（开放自由面背景淹没浮力）

### P3 法国历史 F2 环带回归
- OFF 10s：U≡0（背景干净）
- ON 10s：`U_rms_end≈7.3e-4`（与历史 ~7e-4 同量级）
- ON 20s：后期窗 `U_rms_avg≈5.9e-4`，末≈5.2e-4
- **能力回归 PASS**；未授双回流 / 稳态

### P4 论文参考盘点（仅盘点）
- Fig.7/8：复用 Doc204 数字化（受限可用）
- Fig.6：本树无独立数字化 → **不足**
- **未做**与计算结果的叠图对比（属 P5）

### P6 Linux 交接
- 空目录清洁复跑：弹性 0.05s + NC EM 1s **均 PASS**
- 授 **`LINUX_SOLVER_HANDOFF_READY`**
- 软件树：`docs/ophelie/software/`（此前已并入 linux handoff）

### 仓库卫生（穿插）
- 228–240 计划文进 `archive/DOCS_201_250_MERGED.md`；`doc*_runs` gitignore；代码勿堆 docs

---

## 二、没做什么，以及原因

| 未做项 | 原因 |
|---|---|
| **P5** `FRENCH_PAPER_NC` 三联图 + 主因修复 | 本轮明确**暂缓**：先收甲方两案 + Linux 交接；法国论文主线留下一轮 |
| Fig.8/7 **定量叠图对比** | 依赖 P5 正式耦合/或对标跑；P4 只完成参考可用性盘点 |
| Fig.6 热源数字化补齐 | 本树缺独立提取件；不伪造数值场 |
| F2 `DIAGNOSTIC_TWO_CELL` / 准稳态 | 20s 只证明浮力信号与历史量级；未做双胞分区诊断与稳定性窗 |
| 甲方 NC **自然对流环流验证** | β ON/OFF 速度几乎相同，开放 FS 残差主导；短窗不能宣称环流已验证 |
| 弹性 **整圈** / 长时 | 门控仅 S2 0.05s 工程短窗 |
| 合同验证 / 实验 1000°C / 论文稳态复现 | 计划边界外；handoff 标签明确不覆盖 |
| WIP **git commit / push** | 用户未要求提交；包内以 patch + 源码快照交付 |
| 软件机本机 GUI 联调 | 物理机交付可调用基线；对面按 software 文档自行 build/prepare |

---

## 三、软件可声明 vs 不可声明

**可声明**
- 甲方弹性：dense Reload + 真实弹性短窗 0.05s 可从清洁目录启动
- 甲方 NC：真实 EM + 热 1s（ASSUMED 15 kW）可从清洁目录启动
- 法国 F2：指定环带冻结温场下浮力能力仍在

**不可声明**
- “自然对流环流已验证”
- “论文 Fig.7/8 已复现”
- “合同模块 / 实验测点已通过”
- 粗 Reload 483 粒子桨可当交付基线

---

## 四、包内路径速查

```text
README.md                          # 入口
docs/241_*.md                      # 本轮计划
docs/WORK_SUMMARY_ZH.md            # 本文
results/doc241_runs/               # 全部运行记录与判定
code/patches/                      # 相对 HEAD 的 WIP patch
code/snapshot/                     # 改动/新增源码全文
code/BASE_COMMIT.txt
reference_fig7_fig8_Doc204/        # P4 用数字化参考（非作者 CFD）
software/                          # Linux 软件交接树摘要/全文
MANIFEST.md
```

<!-- ===== END OPHELIE_Doc241_round_handoff_20261124_WORK_SUMMARY_ZH.md ===== -->


<!-- ===== BEGIN OPHELIE_Doc242_round_handoff_lean_20261124_WORK_SUMMARY_ZH.md ===== -->

# OPHELIE_Doc242_round_handoff_lean_20261124_WORK_SUMMARY_ZH.md

# Doc242 工作摘要（简体中文）

日期：2026-11-24。分支：`feature/EM-Biot-Savart` @ `16d67c9d0` + WIP。  
门控：**`LINUX_INTEGRATION_BASELINE_READY`**

## 本轮目标

Linux 软件开工前最后一轮物理整改与 CLI 闭环：弹性异常定位/默认、甲方 NC 固定滑移+EM、法国场图与 Fig.8 对比、两案 CLI e2e。**不以论文稳态/合同验收作为软件开工条件。**

## 按阶段做了什么

| 阶段 | 内容 | 结果标签 |
|---|---|---|
| A | 冻结 HEAD/WIP/binary；F0 到达时间/固体覆盖/场审计 | 基线可审 |
| B | B0/B1/B2 对照；静水 IC / hold-hydro；生产默认 IC-only | workflow 可完成；**质量 FAIL** |
| C | 固定平坦滑移顶 OFF/ON 0.2s→1s；真实 EM 15 kW 1s | `BUOYANCY_SIGNAL_RESOLVED`；`EM_THERMAL_WIRING_PASS` |
| D | F2 20s 带场；论文 50 kW 1s+5s 叠 Fig.8 | `direction_observed`；`PAPER_COMPARISON_AVAILABLE`；主因 `heat_loss_dominant_cold_wall` |
| E | CLI：`fixed_absorbed_power_em`、禁 zero-q、prepare 传 dp；负例；E1/E2 GPU | **两案 WORKFLOW_PASS** |
| F | 能力清单与判定冻结 | **`LINUX_INTEGRATION_BASELINE_READY`** |

## 关键数值（便于扫一眼）

- 弹性 CLI E2：`reached_end_time=1`，`u_max≈1.24 m/s`，`rho_rel_max≈0.54`，`numerical_quality_ok=0`
- 甲方 NC EM：`P_joule=15000`，`zero_q=false`，`em_status=RUN`
- 法国 5s：dT≈37 K，T\*_mid≈0.99，散热/P≈3.1

## 最终八问（答）

1. **尖峰**：约 0.012–0.015s 仍可 ~3.5–4.3 m/s，末态 ~1.2 m/s；IC-only 改善早期，未消中段；静止桨见 B1。
2. **密度**：偏差仍在（B ~35–40%；E2 ~54%）；未改 E / 未清零速度。
3. **固体时间**：`solid_time_uncovered=0`，覆盖完整。
4. **甲方 NC**：固定域浮力可分辨；EM 短跑 **ASSUMED 15 kW**（≠ 电源 30 kW）。
5. **F2**：轴↓/环↑/缘↓；不授双胞/稳态。
6. **Fig.8**：T\* 仍近 1；主因冷壁散热过大；附图在 `05_french_paper_joint/*/figures`。
7. **CLI**：validate→prepare→run 两案通过；场 VTP 可读；非法模式拒绝。
8. **齐套**：代码/示例/结果/限制清单齐；**不能**声称质量门、论文验收、合同验证。

## 未授 / 未做

- 未授：`PAPER_VALIDATED`、`NUMERICAL_QUALITY_PASS`
- 未开 Qt GUI；未再延法国时窗；未自动 git commit

## 转交

- 软件：桌面装配 / 进度取消 / 结果加载（物理与论文验收解绑）
- 物理续：弹性尖峰与密度；法国 h/边界审计 — 见 `KNOWN_LIMITATIONS.md`

<!-- ===== END OPHELIE_Doc242_round_handoff_lean_20261124_WORK_SUMMARY_ZH.md ===== -->


<!-- ===== BEGIN OPHELIE_Doc243_round_handoff_lean_20261124_WORK_SUMMARY_ZH.md ===== -->

# OPHELIE_Doc243_round_handoff_lean_20261124_WORK_SUMMARY_ZH.md

# Doc243 工作摘要（简体中文）

日期：2026-11-24。分支：`feature/EM-Biot-Savart`。  
本包：**精简交接**（无 VTP、无完整大 run.log）。

## 目标

72h 物理收口 + 法国优先复现：物性真消费、边界审计、非均匀初猜、耦合矩阵；搅拌只定位首坏。

## 法国 P0/P1

| 项 | 结果 |
|---|---|
| 物性接线 | PASS（`--temperature-dependent-properties` + ParticleCpCK + local μ） |
| P0-C 面积/0-D | 实现 PASS；`bare_walls`@1473K 与 50kW **热力学不兼容**（loss/P≈3.14） |
| 主 IC | coldwall_skin（轴↓缘↑） |
| F_SEED Fig.8 带 | 缘永不升 → 弃 |
| F_SENS skull+crust | loss/P≈**1.25**，主栈；30s 仍 TRANSIENT |
| 0-D Teq（默认壳） | ≈1260 K → **不盲延 60s** |
| Fig.8 / PAPER | **未达成**（壳厚 ASSUMED/OPEN） |

主栈：`coldwall + 变物性 + --skull-resistance --free-crust-resistance`。

## P2 搅拌

| 步 | 结果 |
|---|---|
| 粗 Reload | INVALID（dp=0.01） |
| dp=5mm ω=0 | 跑完；rho_rel≈0.56；u≈1.9 |
| B2 去载荷 | 尖峰不变 → 首坏=**流体** |
| HOLD | rho→0.048；u 仍坏 |
| skip dens-sum | rho **更差** → **失败并撤销默认** |
| 恢复冒烟 | dens_sum 默认已恢复 |

未授：`NUMERICAL_QUALITY_PASS`。中期债：`background_hydro_p`。

## 包内有什么

- 计划书 `243_*.md`、`doc243_runs/` 精简结果（CSV/JSON/MD/SVG/种子）
- 本轮关键关键源码副本 + `doc243_code.diff`
- 本摘要与 README

## 故意不装

- 全部 `*.vtp` / `out/` 粒子场
- 完整大 `run.log`（仅保留各跑 `run_excerpt.txt` 头尾）
- 二进制、build、Reload 粒子库

## 未做 / 未授

- 未授 PAPER_VALIDATED / ENGINEERING_QUASI_STEADY / NUMERICAL_QUALITY_PASS  
- 未自动 git commit / push

<!-- ===== END OPHELIE_Doc243_round_handoff_lean_20261124_WORK_SUMMARY_ZH.md ===== -->

