#include "solid_to_shell_coupling.h"

namespace SPH
{
namespace solid_dynamics
{
CouplingPart::CouplingPart(BaseContactRelation &contact_relation, const std::string &body_part_name, Real distance_factor)
    : BodyPartByParticle(contact_relation.getSPHBody(), body_part_name),
      DataDelegateContact(contact_relation),
      distance_factor_(distance_factor),
      is_coupled_(base_particles_.registerStateVariable<int>("IsCoupled", 0))
{
    // update contact relation before tagging particles
    dynamic_cast<RealBody *>(&sph_body_)->updateCellLinkedList();
    for (auto *contact_body : contact_bodies_)
        dynamic_cast<RealBody *>(contact_body)->updateCellLinkedList();
    contact_relation.updateConfiguration();

    // tag the particles
    TaggingParticleMethod tagging_particle_method = std::bind(&CouplingPart::tagManually, this, _1);
    tagParticles(tagging_particle_method);
};

void CouplingPart::tagManually(size_t index_i)
{
    Real dp_i = sph_body_.getSPHBodyResolutionRef();
    for (size_t k = 0; k < contact_configuration_.size(); k++)
    {
        Real dp_k = contact_bodies_[k]->getSPHBodyResolutionRef();
        Real dp = 0.5 * (dp_i + dp_k);
        const auto &contact_neighbors = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n < contact_neighbors.current_size_; n++)
        {
            if (contact_neighbors.r_ij_[n] < distance_factor_ * dp)
            {
                is_coupled_[index_i] = 1;
                body_part_particles_.push_back(index_i);
                // stop searching once the particle is found to be coupled
                return;
            }
        }
    }
};

NearestNeighborSolidVelocityConstraint::SearchNearestParticle::SearchNearestParticle(BodyPartByParticle &body_part, BaseContactRelation &contact_relation)
    : BaseLocalDynamics<BodyPartByParticle>(body_part),
      DataDelegateContact(contact_relation)
{
    const auto max = std::numeric_limits<size_t>::max();
    coupling_ids_.resize(particles_->TotalRealParticles(),
                         std::make_pair(max, max));
    for (auto *contact_particle : contact_particles_)
        contact_is_coupled_.emplace_back(contact_particle->registerStateVariable<int>("IsCoupled"));
};

void NearestNeighborSolidVelocityConstraint::SearchNearestParticle::update(size_t index_i, Real)
{
    const auto max = std::numeric_limits<size_t>::max();
    size_t k_min = max;
    size_t j_min = max;
    Real r_min = std::numeric_limits<Real>::max();
    for (size_t k = 0; k < contact_configuration_.size(); k++)
    {
        const auto &contact_neighbors = (*contact_configuration_[k])[index_i];
        const auto *is_coupled_k = contact_is_coupled_[k];
        for (size_t n = 0; n < contact_neighbors.current_size_; n++)
        {
            size_t index_j = contact_neighbors.j_[n];
            if (is_coupled_k[index_j] != 1)
                continue;
            if (contact_neighbors.r_ij_[n] < r_min)
            {
                r_min = contact_neighbors.r_ij_[n];
                k_min = k;
                j_min = index_j;
            }
        }
    }
    coupling_ids_[index_i] = std::make_pair(k_min, j_min);
};

NearestNeighborSolidVelocityConstraint::NearestNeighborSolidVelocityConstraint(BodyPartByParticle &body_part, BaseContactRelation &contact_relation)
    : MotionConstraint<BodyPartByParticle>(body_part),
      DataDelegateContact(contact_relation),
      search_nearest_particle_(body_part, contact_relation)
{
    for (auto *contact_particle : contact_particles_)
        contact_vel_.emplace_back(contact_particle->getVariableDataByName<Vecd>("Velocity"));
};

void NearestNeighborSolidVelocityConstraint::update(size_t index_i, Real)
{
    auto [k, index_j] = search_nearest_particle_.get_nearest_id(index_i);
    // Only update the velocity when the nearest neighbor is found.
    if (const auto max = std::numeric_limits<size_t>::max(); k != max && index_j != max)
        vel_[index_i] = contact_vel_[k][index_j];
}

NearestNeighborShellForceConstraint::SearchNearestParticle::SearchNearestParticle(BodyPartByParticle &body_part, BaseContactRelation &contact_relation)
    : BaseLocalDynamics<BodyPartByParticle>(body_part),
      DataDelegateContact(contact_relation)
{
    const auto max = std::numeric_limits<size_t>::max();
    for (auto *contact_particle : contact_particles_)
    {
        coupling_ids_.emplace_back(particles_->TotalRealParticles(), max);
        contact_is_coupled_.emplace_back(contact_particle->registerStateVariable<int>("IsCoupled"));
    }
};

void NearestNeighborShellForceConstraint::SearchNearestParticle::update(size_t index_i, Real)
{
    const auto max = std::numeric_limits<size_t>::max();
    for (size_t k = 0; k < contact_configuration_.size(); k++)
    {
        const auto &contact_neighbors = (*contact_configuration_[k])[index_i];
        const auto *is_coupled_k = contact_is_coupled_[k];

        size_t j_min = max;
        Real r_min = std::numeric_limits<Real>::max();

        for (size_t n = 0; n < contact_neighbors.current_size_; n++)
        {
            size_t index_j = contact_neighbors.j_[n];
            if (is_coupled_k[index_j] != 1)
                continue;
            if (contact_neighbors.r_ij_[n] < r_min)
            {
                r_min = contact_neighbors.r_ij_[n];
                j_min = index_j;
            }
        }

        coupling_ids_[k][index_i] = j_min;
    }
};

NearestNeighborShellForceConstraint::NearestNeighborShellForceConstraint(BodyPartByParticle &body_part, BaseContactRelation &contact_relation)
    : BaseForcePrior<BodyPartByParticle>(body_part, "CouplingForce"),
      DataDelegateContact(contact_relation),
      search_nearest_particle_(body_part, contact_relation)
{
    for (auto *contact_particle : contact_particles_)
        contact_force_.emplace_back(contact_particle->getVariableDataByName<Vecd>("Force"));
};

void NearestNeighborShellForceConstraint::interaction(size_t index_i, Real)
{
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        size_t index_j = search_nearest_particle_.get_nearest_id(index_i, k);
        // Only update the force when the nearest neighbor is found.
        if (index_j != std::numeric_limits<size_t>::max())
            force += contact_force_[k][index_j];
    }
    current_force_[index_i] = force;
}

InterpolationSolidVelocityConstraint::InterpolationSolidVelocityConstraint(BodyPartByParticle &body_part, BaseContactRelation &contact_relation)
    : MotionConstraint<BodyPartByParticle>(body_part),
      DataDelegateContact(contact_relation)
{
    for (auto *contact_particle : contact_particles_)
    {
        contact_is_coupled_.emplace_back(contact_particle->getVariableDataByName<int>("IsCoupled"));
        contact_Vol_.emplace_back(contact_particle->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_vel_.emplace_back(contact_particle->getVariableDataByName<Vecd>("Velocity"));
    }
};

void InterpolationSolidVelocityConstraint::update(size_t index_i, Real)
{
    Vecd vel = Vecd::Zero();
    Real weight_ttl = 0;
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
            weight_ttl += weight_j;
            vel += weight_j * vel_k[index_j];
        }
    }
    vel_[index_i] = vel / weight_ttl;
}

InterpolationShellForceConstraint::InterpolationShellForceConstraint(BodyPartByParticle &body_part, BaseContactRelation &contact_relation)
    : BaseForcePrior<BodyPartByParticle>(body_part, "CouplingForce"),
      DataDelegateContact(contact_relation)
{
    for (auto *contact_particle : contact_particles_)
    {
        contact_is_coupled_.emplace_back(contact_particle->getVariableDataByName<int>("IsCoupled"));
        contact_Vol_.emplace_back(contact_particle->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_force_.emplace_back(contact_particle->getVariableDataByName<Vecd>("Force"));
    }
};

void InterpolationShellForceConstraint::interaction(size_t index_i, Real)
{
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k != contact_configuration_.size(); ++k)
    {
        Real weight_ttl_k = 0;
        Vecd force_k_ttl = Vecd::Zero();

        const int *is_coupled_k = contact_is_coupled_[k];
        const Real *Vol_k = contact_Vol_[k];
        const Vecd *force_k = contact_force_[k];
        const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            if (is_coupled_k[index_j] == 0)
                continue;

            Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];
            weight_ttl_k += weight_j;
            force_k_ttl += weight_j * force_k[index_j];
        }

        force += force_k_ttl / weight_ttl_k;
    }
    current_force_[index_i] = force;
}
} // namespace solid_dynamics
} // namespace SPH