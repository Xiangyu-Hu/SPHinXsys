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
} // namespace SPH
