#include "velocity_gradient.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
VelocityGradient<Contact<Wall>>::VelocityGradient(BaseContactRelation &wall_contact_relation)
    : InteractionWithWall<VelocityGradient>(wall_contact_relation),
      distance_from_wall_(*particles_->getVariableByName<Vecd>("DistanceFromWall")) {}
//=================================================================================================//
void VelocityGradient<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    Matd vel_grad = Matd::Zero();
    const Vecd &distance_from_wall = distance_from_wall_[index_i];
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &vel_ave_k = *(wall_vel_ave_[k]);
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            const Vecd &e_ij = contact_neighborhood.e_ij_[n];

            Vecd distance_diff = distance_from_wall - contact_neighborhood.r_ij_[n] * e_ij;
            Real factor = 1.0 - distance_from_wall.dot(distance_diff) / distance_from_wall.squaredNorm();
            Vecd nablaW_ijV_j = contact_neighborhood.dW_ijV_j_[n] * e_ij;
            vel_grad -= factor * (vel_[index_i] - vel_ave_k[index_j]) * nablaW_ijV_j.transpose();
        }
    }

    vel_grad_[index_i] += vel_grad;
}
//=================================================================================================//
DistanceFromWall::DistanceFromWall(BaseContactRelation &wall_contact_relation)
    : LocalDynamics(wall_contact_relation.getSPHBody()), FSIContactData(wall_contact_relation),
      spacing_ref_(sph_body_.sph_adaptation_->ReferenceSpacing()),
      distance_from_wall_(*particles_->registerSharedVariable<Vecd>("DistanceFromWall"))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        wall_n_.push_back(&(contact_particles_[k]->n_));
        wall_phi_.push_back(contact_particles_[k]->getVariableByName<Real>("SignedDistance"));
    }
}
//=================================================================================================//
void DistanceFromWall::interaction(size_t index_i, Real dt)
{
    Vecd distance = 100.0 * spacing_ref_ * Vecd::Ones();
    Vecd normal = Vecd::Ones();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &n_k = *(wall_n_[k]);
        StdLargeVec<Real> &phi_k = *(wall_phi_[k]);
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            const Vecd &e_ij = contact_neighborhood.e_ij_[n];

            Vecd temp = contact_neighborhood.r_ij_[n] * e_ij + phi_k[index_j] * n_k[index_j];
            if (temp.squaredNorm() < distance.squaredNorm())
            {
                distance = temp;       // more reliable distance
                normal = n_k[index_j]; // more reliable normal
            }
        }
    }

    Vecd normal_distance = distance.dot(normal) * normal;
    Real limiter = SMIN(3.0 * (distance - normal_distance).norm() / spacing_ref_, 1.0);
    distance_from_wall_[index_i] = (1.0 - limiter) * normal_distance + limiter * distance;
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
