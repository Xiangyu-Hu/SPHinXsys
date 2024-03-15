#include "velocity_gradient.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
VelocityGradient<Contact<Wall>>::VelocityGradient(BaseContactRelation &wall_contact_relation)
    : InteractionWithWall<VelocityGradient>(wall_contact_relation) {}
//=================================================================================================//
void VelocityGradient<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    Matd vel_grad = Matd::Zero();
    Vecd distance_to_fluid = Vecd::Zero();
    Vecd distance_from_surface = Vecd::Zero();
    Real total_weights = 1.0e-6;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &vel_ave_k = *(wall_vel_ave_[k]);
        StdLargeVec<Vecd> &n_k = *(wall_n_[k]);
        StdLargeVec<Real> &phi_k = *(wall_phi_[k]);
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            const Vecd &e_ij = contact_neighborhood.e_ij_[n];
            Real W_ij = contact_neighborhood.W_ij_[n];

            distance_to_fluid += W_ij * contact_neighborhood.r_ij_[n] * e_ij;
            distance_from_surface += W_ij * phi_k[index_j] * n_k[index_j];
            total_weights += W_ij;

            Vecd nablaW_ijV_j = contact_neighborhood.dW_ijV_j_[n] * e_ij;
            vel_grad -=  (vel_[index_i] - vel_ave_k[index_j]) * nablaW_ijV_j.transpose();
        }
    }

    distance_to_fluid /= total_weights;
    distance_from_surface /= total_weights;

    if (index_i == 0)
    {
        Vecd distance_diff = distance_to_fluid + distance_from_surface;
        Vecd dummy = vel_[index_i] * distance_diff.dot(distance_from_surface) / distance_diff.squaredNorm();
        dummy += Vecd::Zero();
    }

    vel_grad_[index_i] += ReflectiveFactor(distance_to_fluid, distance_from_surface) * vel_grad;
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
