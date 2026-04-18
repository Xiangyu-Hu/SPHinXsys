#include "density_gradient.hpp"
#include "data_type.h"

template class SPH::fluid_dynamics::DensityGradient<
    SPH::Inner<SPH::LinearGradientCorrection>>;

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
DensityGradient<Contact<Wall>>::DensityGradient(BaseContactRelation &wall_contact_relation)
    : InteractionWithWall<DensityGradient>(wall_contact_relation),
      distance_from_wall_(particles_->getVariableDataByName<Vecd>("DistanceFromWall")),
      rho_(this->particles_->template getVariableDataByName<Real>("Density")),
      rho_grad_(this->particles_->template getVariableDataByName<Vecd>("DensityGradient"))
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        wall_rho_.push_back(
            this->contact_particles_[k]->template getVariableDataByName<Real>("Density"));
    }
}
//=================================================================================================//
void DensityGradient<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    Vecd density_grad = Vecd::Zero();
    const Vecd &distance_from_wall = distance_from_wall_[index_i];
    Real rho_i = rho_[index_i];

    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real *Vol_k = wall_Vol_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];

        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            const Vecd &e_ij = contact_neighborhood.e_ij_[n];
            Real rho_j = wall_rho_[k][index_j];
            Vecd nablaW_ijV_j = contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * e_ij;
            Vecd distance_diff = distance_from_wall - contact_neighborhood.r_ij_[n] * e_ij;
            Real factor = 1.0 - distance_from_wall.dot(distance_diff) / distance_from_wall.squaredNorm();
            density_grad -= factor * (rho_i - rho_j) * nablaW_ijV_j;
        }
    }

    rho_grad_[index_i] += density_grad;
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
