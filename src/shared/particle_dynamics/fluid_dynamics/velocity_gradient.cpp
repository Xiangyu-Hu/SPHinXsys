#include "velocity_gradient.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
void VelocityGradient<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    Matd vel_grad = Matd::Zero();
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

            Vecd distance_to_fluid = contact_neighborhood.r_ij_[n] * e_ij;
            Vecd distance_from_surface = phi_k[index_j] * n_k[index_j];
            Vecd nablaW_ijV_j = contact_neighborhood.dW_ijV_j_[n] * e_ij;
            vel_grad -= ReflectiveFactor(distance_to_fluid, distance_from_surface) *
                        (vel_[index_i] - vel_ave_k[index_j]) * nablaW_ijV_j.transpose();
        }
    }

    vel_grad_[index_i] += vel_grad;
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
