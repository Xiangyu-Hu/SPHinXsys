# A3 SYCL/CK 移植说明（2026-11-05）

## 状态
**GPU SYCL 构建已通过**（`test_2d_gpu_paper_a3_windows_frame_D4`）。

## 库改动
- `LocalIsotropicDiffusion::InterParticleCoeff`（设备读 `ThermalConductivity`）
- `Robin<DiffusionType>` CK（`diffusion_dynamics_ck.*`）
- Contact 换热变量命名改为 `PhiTransferFrom{ContactBody}`（与经典一致）

## API
`DiffusionRelaxationCK` + RK2 + `Robin<LocalIsotropicDiffusion>` × 内外对流边界。

## 后端 meta
`sycl_ck_gpu`（非 classic）。
