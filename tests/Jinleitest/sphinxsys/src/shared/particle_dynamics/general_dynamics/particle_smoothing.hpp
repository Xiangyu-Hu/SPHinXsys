#ifndef PARTICLE_SMOOTHING_HPP
#define PARTICLE_SMOOTHING_HPP

#include "particle_smoothing.h"

namespace SPH
{
//=================================================================================================//
template <typename VariableType>
ParticleSmoothing<VariableType>::
    ParticleSmoothing(BaseInnerRelation &inner_relation, const std::string &variable_name)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      W0_(sph_body_.getSPHAdaptation().getKernel()->W0(ZeroVecd)),
      smoothed_(particles_->template getVariableDataByName<VariableType>(variable_name)),
      temp_(particles_->registerStateVariable<VariableType>(variable_name + "_temp")) {}
//=================================================================================================//
template <typename VariableType>
void ParticleSmoothing<VariableType>::interaction(size_t index_i, Real dt)
{
    Real weight = W0_;
    VariableType summation = W0_ * smoothed_[index_i];
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        summation += inner_neighborhood.W_ij_[n] * smoothed_[index_j];
        weight += inner_neighborhood.W_ij_[n];
    }
    temp_[index_i] = summation / (weight + TinyReal);
}
//=================================================================================================//
template <typename VariableType>
void ParticleSmoothing<VariableType>::update(size_t index_i, Real dt)
{
    smoothed_[index_i] = temp_[index_i];
}
//=================================================================================================//
template <typename VariableType>
ParticleSnapshotAverage<VariableType>::
    ParticleSnapshotAverage(SPHBody &sph_body, const std::string &variable_name)
    : LocalDynamics(sph_body),
      target_variable_(particles_->template getVariableDataByName<VariableType>(variable_name)),
      averaged_variable_(particles_->template registerStateVariable<VariableType>("Averaged" + variable_name))
{
    particles_->addVariableToWrite<VariableType>("Averaged" + variable_name);
}
//=================================================================================================//
template <typename VariableType>
void ParticleSnapshotAverage<VariableType>::setupDynamics(Real dt)
{
    number_of_snapshot_++;
}
//=================================================================================================//
template <typename VariableType>
void ParticleSnapshotAverage<VariableType>::update(size_t index_i, Real dt)
{
    averaged_variable_[index_i] += (target_variable_[index_i] - averaged_variable_[index_i]) / Real(number_of_snapshot_);
}
//=================================================================================================//
} // namespace SPH
#endif // PARTICLE_SMOOTHING_HPP
