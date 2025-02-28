#include "solid_to_shell_coupling.h"

namespace SPH
{
namespace solid_dynamics
{
NearestNeighborSolidVelocityConstraint::SearchNearestParticle::SearchNearestParticle(BodyPartByParticle &body_part, BaseContactRelation &contact_relation)
    : BaseLocalDynamics<BodyPartByParticle>(body_part),
      DataDelegateContact(contact_relation)
{
    const auto inf = std::numeric_limits<size_t>::infinity();
    coupling_ids_.resize(particles_->TotalRealParticles(),
                         std::make_pair(inf, inf));
    for (auto *contact_particle : contact_particles_)
        contact_is_coupled_.emplace_back(contact_particle->registerStateVariable<int>("IsCoupled"));
};

void NearestNeighborSolidVelocityConstraint::SearchNearestParticle::update(size_t index_i, Real)
{
    const auto inf = std::numeric_limits<size_t>::infinity();
    size_t k_min = inf;
    size_t j_min = inf;
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
    if (const auto inf = std::numeric_limits<size_t>::infinity(); k != inf && index_j != inf)
        vel_[index_i] = contact_vel_[k][index_j];
}

NearestNeighborShellForceConstraint::SearchNearestParticle::SearchNearestParticle(BodyPartByParticle &body_part, BaseContactRelation &contact_relation)
    : BaseLocalDynamics<BodyPartByParticle>(body_part),
      DataDelegateContact(contact_relation)
{
    const auto inf = std::numeric_limits<size_t>::infinity();
    for (auto *contact_particle : contact_particles_)
    {
        coupling_ids_.emplace_back(StdLargeVec<size_t>(particles_->TotalRealParticles(), inf));
        contact_is_coupled_.emplace_back(contact_particle->registerStateVariable<int>("IsCoupled"));
    }
};

void NearestNeighborShellForceConstraint::SearchNearestParticle::update(size_t index_i, Real)
{
    const auto inf = std::numeric_limits<size_t>::infinity();
    for (size_t k = 0; k < contact_configuration_.size(); k++)
    {
        const auto &contact_neighbors = (*contact_configuration_[k])[index_i];
        const auto *is_coupled_k = contact_is_coupled_[k];

        size_t j_min = inf;
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
    for (size_t k = 0; k != contact_particles_.size(); ++k)
        contact_force_.emplace_back(contact_particles_[k]->getVariableDataByName<Vecd>("Force"));
};

void NearestNeighborShellForceConstraint::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        size_t index_j = search_nearest_particle_.get_nearest_id(index_i, k);
        // Only update the force when the nearest neighbor is found.
        if (index_j != std::numeric_limits<size_t>::infinity())
            force += contact_force_[k][index_j];
    }
    current_force_[index_i] = force;
}
} // namespace solid_dynamics
} // namespace SPH