#include "derived_solid_state.h"

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
DisplacementAndPosition::DisplacementAndPosition(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_pos0_(particles_->registerStateVariableOnlyFrom<Vecd>("InitialPosition", "Position")),
      dv_displacement_(particles_->registerStateVariableOnly<Vecd>("Displacement")) {}
//=============================================================================================//
UpdateDisplacementFromPosition::UpdateDisplacementFromPosition(SPHBody &sph_body)
    : DisplacementAndPosition(sph_body) {}
//=============================================================================================//
UpdatePositionFromDisplacement::UpdatePositionFromDisplacement(SPHBody &sph_body)
    : DisplacementAndPosition(sph_body) {}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
