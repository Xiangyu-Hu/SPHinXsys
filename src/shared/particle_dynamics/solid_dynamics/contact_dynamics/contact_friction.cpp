#include "contact_friction.h"

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
PairwiseFrictionFromWall::
    PairwiseFrictionFromWall(BaseContactRelation &contact_relation, Real eta)
    : LocalDynamics(contact_relation.getSPHBody()), ContactWithWallData(contact_relation),
      eta_(eta), Vol_(*particles_->getVariableByName<Real>("VolumetricMeasure")),
      mass_(*particles_->getVariableByName<Real>("Mass")),
      vel_(*particles_->getVariableByName<Vecd>("Velocity"))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        wall_vel_n_.push_back(contact_particles_[k]->registerSharedVariable<Vecd>("Velocity"));
        wall_n_.push_back(contact_particles_[k]->template getVariableByName<Vecd>("NormalDirection"));
        wall_Vol_n_.push_back(contact_particles_[k]->getVariableByName<Real>("VolumetricMeasure"));
    }
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
