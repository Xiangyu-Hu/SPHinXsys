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
RandomizeParticlePosition::RandomizeParticlePosition(SPHBody &sph_body)
    : LocalDynamics(sph_body), GeneralDataDelegateSimple(sph_body),
      pos_(particles_->pos_), randomize_scale_(sph_body.sph_adaptation_->MinimumSpacing()) {}
//=================================================================================================//
void RandomizeParticlePosition::update(size_t index_i, Real dt)
{
    Vecd &pos_n_i = pos_[index_i];
    for (int k = 0; k < pos_n_i.size(); ++k)
    {
        pos_n_i[k] += dt * (Real(rand() / RAND_MAX) - 0.5) * 2.0 * randomize_scale_;
    }
}
//=================================================================================================//
} // namespace SPH
