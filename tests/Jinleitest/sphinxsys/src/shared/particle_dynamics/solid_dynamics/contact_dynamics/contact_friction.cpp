#include "contact_friction.h"

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
PairwiseFrictionFromWall::
    PairwiseFrictionFromWall(BaseContactRelation &contact_relation, Real eta)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      eta_(eta), Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity"))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        wall_vel_n_.push_back(contact_particles_[k]->registerStateVariable<Vecd>("Velocity"));
        wall_n_.push_back(contact_particles_[k]->template getVariableDataByName<Vecd>("NormalDirection"));
        wall_Vol_n_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
    }
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
