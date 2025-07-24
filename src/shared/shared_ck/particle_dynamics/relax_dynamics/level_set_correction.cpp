#include "level_set_correction.h"

namespace SPH
{
//=================================================================================================//
LevelsetBounding::LevelsetBounding(NearShapeSurface &body_part)
    : BaseLocalDynamics<BodyPartByCell>(body_part),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      level_set_(body_part.getLevelSetShape().getLevelSet()),
      constrained_distance_(0.5 * sph_body_.getSPHAdaptation().MinimumSpacing()) {}
//=================================================================================================//
LevelsetKernelGradientIntegral::LevelsetKernelGradientIntegral(NearShapeSurface &body_part)
    : BaseLocalDynamics<BodyPartByCell>(body_part),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_residue_(particles_->registerStateVariableOnly<Vecd>("ZeroGradientResidue")),
      level_set_(body_part.getLevelSetShape().getLevelSet()) {}
//=================================================================================================//
} // namespace SPH
