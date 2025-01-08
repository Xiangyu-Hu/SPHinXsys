#include "active_muscle_dynamics.h"

namespace SPH
{
namespace active_muscle_dynamics
{
//=================================================================================================//
MuscleActivation::MuscleActivation(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      pos0_(particles_->registerStateVariableFrom<Vecd>("InitialPosition", "Position")),
      active_contraction_stress_(particles_->getVariableDataByName<Real>("ActiveContractionStress")){};
//=================================================================================================//
} // namespace active_muscle_dynamics
} // namespace SPH
