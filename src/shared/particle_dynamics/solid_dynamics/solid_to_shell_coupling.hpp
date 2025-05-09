#ifndef SOLID_TO_SHELL_CONSTRAINT_HPP
#define SOLID_TO_SHELL_CONSTRAINT_HPP

#include "solid_to_shell_coupling.h"

namespace SPH
{
namespace solid_dynamics
{
template <class DynamicsIdentifier>
TotalWeightComputation<DynamicsIdentifier>::
    TotalWeightComputation(DynamicsIdentifier &identifier, BaseContactRelation &contact_relation)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      DataDelegateContact(contact_relation),
      total_weight_(this->particles_->template registerStateVariable<Real>("TotalWeight"))
{
    for (auto *contact_particle : contact_particles_)
        contact_Vol_.emplace_back(contact_particle->template getVariableDataByName<Real>("VolumetricMeasure"));
}

template <class DynamicsIdentifier>
void TotalWeightComputation<DynamicsIdentifier>::update(size_t index_i, Real dt)
{
    Real weight_ttl = 0;
    for (size_t k = 0; k != contact_configuration_.size(); ++k)
    {
        const Real *Vol_k = contact_Vol_[k];
        const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];
            weight_ttl += weight_j;
        }
    }
    total_weight_[index_i] = weight_ttl;
}
//=================================================================================================//
template <class DynamicsIdentifier, typename DataType>
ConsistentMapping<DynamicsIdentifier, DataType>::
    ConsistentMapping(DynamicsIdentifier &identifier,
                      BaseContactRelation &contact_relation,
                      const std::string &variable_name,
                      const std::string &contact_variable_name)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      DataDelegateContact(contact_relation),
      interpolated_quantities_(this->particles_->template registerStateVariable<DataType>(variable_name)),
      total_weight_(this->particles_->template getVariableDataByName<Real>("TotalWeight"))
{
    for (auto *contact_particle : contact_particles_)
    {
        contact_Vol_.emplace_back(contact_particle->template getVariableDataByName<Real>("VolumetricMeasure"));
        contact_data_.emplace_back(contact_particle->template getVariableDataByName<DataType>(contact_variable_name));
    }
};

template <class DynamicsIdentifier, typename DataType>
void ConsistentMapping<DynamicsIdentifier, DataType>::update(size_t index_i, Real)
{
    // only consider particles with contact neighbors
    if (total_weight_[index_i] < TinyReal)
        return;

    DataType interpolated_quantities = ZeroData<DataType>::value;
    for (size_t k = 0; k != contact_configuration_.size(); ++k)
    {
        const Real *Vol_k = contact_Vol_[k];
        const DataType *data_k = contact_data_[k];
        const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];
            interpolated_quantities += weight_j * data_k[index_j];
        }
    }
    interpolated_quantities_[index_i] = interpolated_quantities / total_weight_[index_i];
}
//=================================================================================================//
template <class DynamicsIdentifier, typename DataType>
ConservativeMapping<DynamicsIdentifier, DataType>::
    ConservativeMapping(DynamicsIdentifier &identifier,
                        BaseContactRelation &contact_relation,
                        const std::string &variable_name,
                        const std::string &contact_variable_name)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      DataDelegateContact(contact_relation),
      interpolated_quantities_(this->particles_->template registerStateVariable<DataType>(variable_name)),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure"))
{
    for (auto *contact_particle : contact_particles_)
    {
        contact_total_weight_.emplace_back(contact_particle->template getVariableDataByName<Real>("TotalWeight"));
        contact_data_.emplace_back(contact_particle->template getVariableDataByName<DataType>(contact_variable_name));
    }
};

template <class DynamicsIdentifier, typename DataType>
void ConservativeMapping<DynamicsIdentifier, DataType>::update(size_t index_i, Real)
{
    DataType interpolated_quantities = ZeroData<DataType>::value;
    for (size_t k = 0; k != contact_configuration_.size(); ++k)
    {
        const DataType *data_k = contact_data_[k];
        const Real *total_weight_k = contact_total_weight_[k];
        const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            if (total_weight_k[index_j] < TinyReal)
                continue;
            Real weight_j = contact_neighborhood.W_ij_[n] * Vol_[index_i] / total_weight_k[index_j];
            interpolated_quantities += weight_j * data_k[index_j];
        }
    }
    interpolated_quantities_[index_i] = interpolated_quantities;
}
} // namespace solid_dynamics
} // namespace SPH

#endif // SOLID_TO_SHELL_CONSTRAINT_HPP