#include "solid_to_shell_coupling.h"

namespace SPH
{
namespace solid_dynamics
{
TotalWeightComputation::TotalWeightComputation(BodyPartByParticle &body_part, BaseContactRelation &contact_relation)
    : BaseLocalDynamics<BodyPartByParticle>(body_part),
      DataDelegateContact(contact_relation),
      total_weight_(particles_->registerStateVariable<Real>("TotalWeight"))
{
    for (auto *contact_particle : contact_particles_)
    {
        contact_is_coupled_.emplace_back(contact_particle->registerStateVariable<int>("IsCoupled"));
        contact_Vol_.emplace_back(contact_particle->getVariableDataByName<Real>("VolumetricMeasure"));
    }
};

void TotalWeightComputation::update(size_t index_i, Real dt)
{
    Real weight_ttl = TinyReal;
    for (size_t k = 0; k != contact_configuration_.size(); ++k)
    {
        const int *is_coupled_k = contact_is_coupled_[k];
        const Real *Vol_k = contact_Vol_[k];
        const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            if (is_coupled_k[index_j] == 0)
                continue;

            Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];
            weight_ttl += weight_j;
        }
    }
    total_weight_[index_i] = weight_ttl;
}

InterpolationSolidVelocityConstraint::InterpolationSolidVelocityConstraint(BodyPartByParticle &body_part, BaseContactRelation &contact_relation)
    : MotionConstraint<BodyPartByParticle>(body_part),
      DataDelegateContact(contact_relation),
      total_weight_(particles_->getVariableDataByName<Real>("TotalWeight"))
{
    for (auto *contact_particle : contact_particles_)
    {
        contact_is_coupled_.emplace_back(contact_particle->registerStateVariable<int>("IsCoupled"));
        contact_Vol_.emplace_back(contact_particle->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_vel_.emplace_back(contact_particle->getVariableDataByName<Vecd>("Velocity"));
    }
};

void InterpolationSolidVelocityConstraint::update(size_t index_i, Real)
{
    Vecd vel = Vecd::Zero();
    for (size_t k = 0; k != contact_configuration_.size(); ++k)
    {
        const int *is_coupled_k = contact_is_coupled_[k];
        const Real *Vol_k = contact_Vol_[k];
        const Vecd *vel_k = contact_vel_[k];
        const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            if (is_coupled_k[index_j] == 0)
                continue;

            Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];
            vel += weight_j * vel_k[index_j];
        }
    }
    vel_[index_i] = vel / total_weight_[index_i];
}

InterpolationShellForceConstraint::InterpolationShellForceConstraint(BodyPartByParticle &body_part, BaseContactRelation &contact_relation)
    : BaseForcePrior<BodyPartByParticle>(body_part, "CouplingForce"),
      DataDelegateContact(contact_relation),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure"))
{
    for (auto *contact_particle : contact_particles_)
    {
        contact_total_weight_.emplace_back(contact_particle->getVariableDataByName<Real>("TotalWeight"));
        contact_is_coupled_.emplace_back(contact_particle->registerStateVariable<int>("IsCoupled"));
        contact_force_.emplace_back(contact_particle->getVariableDataByName<Vecd>("Force"));
    }
};

void InterpolationShellForceConstraint::interaction(size_t index_i, Real)
{
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k != contact_configuration_.size(); ++k)
    {
        const int *is_coupled_k = contact_is_coupled_[k];
        const Vecd *force_k = contact_force_[k];
        const Real *total_weight_k = contact_total_weight_[k];
        const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            if (is_coupled_k[index_j] == 0)
                continue;

            Real weight_j = contact_neighborhood.W_ij_[n] * Vol_[index_i] / total_weight_k[index_j];
            force += weight_j * force_k[index_j];
        }
    }
    current_force_[index_i] = force;
}
} // namespace solid_dynamics
} // namespace SPH