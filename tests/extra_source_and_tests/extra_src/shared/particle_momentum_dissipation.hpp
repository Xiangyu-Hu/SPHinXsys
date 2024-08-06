#ifndef PARTICLE_MOMENTUM_DISSIPATION_HPP
#define PARTICLE_MOMENTUM_DISSIPATION_HPP
#include "particle_momentum_dissipation.h"

namespace SPH
{
namespace multi_species_continuum
{
//=================================================================================================//
template <typename VariableType>
PorousMediaDampingPairwiseInner<VariableType>::
    PorousMediaDampingPairwiseInner(BaseInnerRelation &inner_relation,
                                    const std::string &variable_name, Real eta)
    : LocalDynamics(inner_relation.getSPHBody()),
      DataDelegateInner(inner_relation),
      Vol_(*particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      mass_(*particles_->getVariableDataByName<Real>("Mass")),
      variable_(*particles_->getVariableDataByName<VariableType>(variable_name)),
      eta_(eta) {}
//=================================================================================================//
template <typename VariableType>
void PorousMediaDampingPairwiseInner<VariableType>::interaction(size_t index_i, Real dt)
{
    Real Vol_i = Vol_[index_i];
    Real mass_i = mass_[index_i];
    VariableType &variable_i = variable_[index_i];

    std::array<Real, MaximumNeighborhoodSize> parameter_b;
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    // forward sweep
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real mass_j = mass_[index_j];

        VariableType variable_derivative = (variable_i - variable_[index_j]);
        parameter_b[n] = eta_ * inner_neighborhood.dW_ij_[n] * Vol_i * Vol_[index_j] * dt / inner_neighborhood.r_ij_[n];

        VariableType increment = parameter_b[n] * variable_derivative / (mass_i * mass_j - parameter_b[n] * (mass_i + mass_j));
        variable_[index_i] += increment * mass_j;
        variable_[index_j] -= increment * mass_i;
    }

    // backward sweep
    for (size_t n = inner_neighborhood.current_size_; n != 0; --n)
    {
        size_t index_j = inner_neighborhood.j_[n - 1];
        Real mass_j = mass_[index_j];

        VariableType variable_derivative = (variable_i - variable_[index_j]);
        VariableType increment = parameter_b[n - 1] * variable_derivative / (mass_i * mass_j - parameter_b[n - 1] * (mass_i + mass_j));

        variable_[index_i] += increment * mass_j;
        variable_[index_j] -= increment * mass_i;
    }
}
//=================================================================================================//
} // namespace multi_species_continuum
} // namespace SPH
#endif // PARTICLE_MOMENTUM_DISSIPATION_HPP