#include "emitter_boundary_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
WithinDisposerIndication::WithinDisposerIndication(OrientedBoxByCell &oriented_box_part)
    : BaseLocalDynamics<OrientedBoxByCell>(oriented_box_part),
      sv_oriented_box_(oriented_box_part.svOrientedBox()),
      sv_total_real_particles_(particles_->svTotalRealParticles()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_life_status_(particles_->registerStateVariable<int>("LifeStatus", 0)) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
