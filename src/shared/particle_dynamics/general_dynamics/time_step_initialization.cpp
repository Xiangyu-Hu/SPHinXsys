#include "time_step_initialization.h"

namespace SPH
{
//=================================================================================================//
TimeStepInitialization::TimeStepInitialization(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr)
    : BaseTimeStepInitialization(sph_body, gravity_ptr), GeneralDataDelegateSimple(sph_body),
      pos_(particles_->pos_), force_prior_(particles_->force_prior_), mass_(particles_->mass_) {}
//=================================================================================================//
void TimeStepInitialization::update(size_t index_i, Real dt)
{
    force_prior_[index_i] = mass_[index_i] * gravity_->InducedAcceleration(pos_[index_i]);
}
//=================================================================================================//
} // namespace SPH
