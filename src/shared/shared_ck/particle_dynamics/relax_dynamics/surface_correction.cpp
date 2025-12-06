#include "surface_correction.h"

namespace SPH
{
//=================================================================================================//
LevelsetBounding::LevelsetBounding(NearShapeSurface &body_part)
    : BaseLocalDynamics<BodyPartByCell>(body_part),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      level_set_(body_part.getLevelSetShape().getLevelSet()),
      constrained_distance_(0.5 * getSPHAdaptation().MinimumSpacing()) {}
//=================================================================================================//
LevelsetKernelGradientIntegral::LevelsetKernelGradientIntegral(SPHBody &sph_body, LevelSetShape &level_set_shape)
    : LocalDynamics(sph_body),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_residual_(particles_->registerStateVariable<Vecd>("KernelGradientIntegral")),
      level_set_(level_set_shape.getLevelSet()) {}
//=================================================================================================//
} // namespace SPH
