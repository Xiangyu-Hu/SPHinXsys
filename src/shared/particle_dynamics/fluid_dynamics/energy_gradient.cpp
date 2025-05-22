#include "energy_gradient.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
EnergyGradient<Contact<Wall>>::EnergyGradient(BaseContactRelation &wall_contact_relation)
    : InteractionWithWall<EnergyGradient>(wall_contact_relation),
      distance_from_wall_(particles_->getVariableDataByName<Vecd>("DistanceFromWall")),
      energy_(this->particles_->template getVariableDataByName<Real>("Energy")),
      energy_grad_(this->particles_->template getVariableDataByName<Vecd>("EnergyGradient")) {}
//=================================================================================================//
void EnergyGradient<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    Vecd energy_grad = Vecd::Zero();
    const Vecd &distance_from_wall = distance_from_wall_[index_i];
    Real energy_i = energy_[index_i];

    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real *Vol_k = wall_Vol_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];

        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            const Vecd &e_ij = contact_neighborhood.e_ij_[n];

            Vecd distance_diff = distance_from_wall - contact_neighborhood.r_ij_[n] * e_ij;
            Real factor = 1.0 - distance_from_wall.dot(distance_diff) / distance_from_wall.squaredNorm();

            Real energy_j = energy_[index_j]; // e_j değerini aldık

            Vecd nablaW_ijV_j = contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * e_ij;

            energy_grad -= factor * (energy_i - energy_j) * nablaW_ijV_j;
        }
    }

    energy_grad_[index_i] += energy_grad;
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
