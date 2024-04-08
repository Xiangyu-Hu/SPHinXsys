#include "contact_friction.h"

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
PairwiseFrictionFromWall::
    PairwiseFrictionFromWall(BaseContactRelation &contact_relation, Real eta)
    : LocalDynamics(contact_relation.getSPHBody()), ContactWithWallData(contact_relation),
      eta_(eta), Vol_(particles_->Vol_), mass_(particles_->mass_),
      vel_(particles_->vel_)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        wall_vel_n_.push_back(&contact_particles_[k]->vel_);
        wall_n_.push_back(&contact_particles_[k]->n_);
        wall_Vol_n_.push_back(&contact_particles_[k]->Vol_);
    }
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
