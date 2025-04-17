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
template <class DynamicsIdentifier>
InterpolationVelocityConstraint<DynamicsIdentifier>::
    InterpolationVelocityConstraint(DynamicsIdentifier &identifier, BaseContactRelation &contact_relation)
    : MotionConstraint<DynamicsIdentifier>(identifier),
      DataDelegateContact(contact_relation),
      total_weight_(this->particles_->template getVariableDataByName<Real>("TotalWeight"))
{
    for (auto *contact_particle : contact_particles_)
    {
        contact_Vol_.emplace_back(contact_particle->template getVariableDataByName<Real>("VolumetricMeasure"));
        contact_vel_.emplace_back(contact_particle->template getVariableDataByName<Vecd>("Velocity"));
    }
};

template <class DynamicsIdentifier>
void InterpolationVelocityConstraint<DynamicsIdentifier>::update(size_t index_i, Real)
{
    // only consider particles with contact neighbors
    if (total_weight_[index_i] < TinyReal)
        return;

    Vecd vel = Vecd::Zero();
    for (size_t k = 0; k != contact_configuration_.size(); ++k)
    {
        const Real *Vol_k = contact_Vol_[k];
        const Vecd *vel_k = contact_vel_[k];
        const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];
            vel += weight_j * vel_k[index_j];
        }
    }
    this->vel_[index_i] = vel / total_weight_[index_i];
}
//=================================================================================================//
template <class DynamicsIdentifier>
InterpolationForceConstraint<DynamicsIdentifier>::
    InterpolationForceConstraint(DynamicsIdentifier &identifier, BaseContactRelation &contact_relation)
    : BaseForcePrior<DynamicsIdentifier>(identifier, "SolidToShellCouplingForce"),
      DataDelegateContact(contact_relation),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure"))
{
    for (auto *contact_particle : contact_particles_)
    {
        contact_total_weight_.emplace_back(contact_particle->template getVariableDataByName<Real>("TotalWeight"));
        contact_force_.emplace_back(contact_particle->template getVariableDataByName<Vecd>("Force"));
    }
};

template <class DynamicsIdentifier>
void InterpolationForceConstraint<DynamicsIdentifier>::interaction(size_t index_i, Real)
{
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k != contact_configuration_.size(); ++k)
    {
        const Vecd *force_k = contact_force_[k];
        const Real *total_weight_k = contact_total_weight_[k];
        const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            if (total_weight_k[index_j] < TinyReal)
                continue;
            Real weight_j = contact_neighborhood.W_ij_[n] * Vol_[index_i] / total_weight_k[index_j];
            force += weight_j * force_k[index_j];
        }
    }
    this->current_force_[index_i] = force;
}
} // namespace solid_dynamics
} // namespace SPH

#endif // SOLID_TO_SHELL_CONSTRAINT_HPP