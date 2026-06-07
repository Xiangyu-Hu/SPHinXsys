# CURSOR_APHI_STAGE10_14_TEAM7_MULTIRESOLUTION_PLAN

Stage 10.14 P8 (document only).

## Issue

Standard TEAM7 air domains are large; uniform fine SPH resolution is expensive.

## Strategy

1. Single-resolution simplified TEAM7 scaffold (`test_3d_aphi_ck_simplified_team7_source_driven`) — **done in P6**.
2. Air-box sensitivity (`test_3d_aphi_ck_boundary_support_policy_diagnostic`) — **done in P5**.
3. Two-level MR: fine conductor/coil, coarse far air, transition buffer ≥ 2× dp.
4. Before MR production: Contact source-driven baselines + probe CSV infrastructure (P7 follow-up).

## Diagnostics before trusting MR

- GMRES iteration/residual vs uniform baseline
- Conductor Joule integral drift < 2% vs finest uniform reference
- `B=curlA` dp sweep at resolution jump
- Contact operator consistency if bodies use different dp

## Particle-count estimate (order of magnitude)

Uniform dp=0.1 on 1.2×1.0×0.3 m³ box ≈ few×10³ particles (see P1 @ 360 particles for TEAM7 fractions on same scale). Full TEAM7 air at dp=0.05 may exceed 10⁶ — MR required.
