#ifndef PARTICLE_SMOOTHING_HPP
#define PARTICLE_SMOOTHING_HPP

#include "particle_smoothing.h"

namespace SPH
{
//=================================================================================================//
template <typename VariableType>
ParticleSmoothing<VariableType>::
    ParticleSmoothing(BaseInnerRelation &inner_relation, const std::string &variable_name)
    : LocalDynamics(inner_relation.getSPHBody()), GeneralDataDelegateInner(inner_relation),
      W0_(sph_body_.sph_adaptation_->getKernel()->W0(ZeroVecd)),
      smoothed_(*particles_->template getVariableByName<VariableType>(variable_name))
{
    particles_->registerVariable(temp_, variable_name + "_temp");
}
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
} // namespace SPH
#endif // PARTICLE_SMOOTHING_HPP
