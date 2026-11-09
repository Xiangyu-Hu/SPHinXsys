# GPU-paper overlays (`extra_src/shared/gpu_paper`)

These files **replace earlier ad-hoc edits under `src/`**. Library tree `src/` stays clean.

| File | Replaces former `src/` patch | Used by |
|---|---|---|
| `gpu_paper_local_isotropic_diffusion.*` | `LocalIsotropicDiffusion` + `InterParticleCoeff` / `dvThermalConductivity` | A3 |
| `gpu_paper_diffusion_robin.*` | CK `Robin<>` contact BC | A3 |
| `gpu_paper_sycl_group.h` | *(new)* `GPU_PAPER_GROUP` → `setWorkGroupSize` | B2 / B4 group scan |
| `gpu_paper_recording.*` | run CSV/JSON recorder | all paper cases |

A5 needs `FSI::InitializeDisplacementCK` from existing `structure_surface_motion.h` in `src/`; the paper case includes that header directly instead of patching `all_fluid_structure_interactions.h`.
