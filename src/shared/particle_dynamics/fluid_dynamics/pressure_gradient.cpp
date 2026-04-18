// pressure_gradient.cpp
#include "pressure_gradient.hpp"
#include "data_type.h"
#pragma message(".cpp is compiled")

template class SPH::fluid_dynamics::PressureGradient<
    SPH::Inner<SPH::LinearGradientCorrection>>;

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
PressureGradient<Contact<Wall>>::PressureGradient(
    BaseContactRelation &wall_contact_relation)
    : InteractionWithWall<PressureGradient>(wall_contact_relation),
      distance_from_wall_(particles_->getVariableDataByName<Vecd>("DistanceFromWall")),
      p_(this->particles_->template getVariableDataByName<Real>("Pressure")),
      p_grad_(this->particles_->template getVariableDataByName<Vecd>("PressureGradient"))
{
    // obtain pressure data from each wall particle set
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        wall_p_.push_back(this->contact_particles_[k]->template getVariableDataByName<Real>("Pressure"));
    }
}
//=================================================================================================//
void PressureGradient<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    Vecd pressure_grad = Vecd::Zero();
    const Vecd &distance_from_wall = distance_from_wall_[index_i];
    Real p_i = p_[index_i];

    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real *Vol_k = wall_Vol_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];

        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            const Vecd &e_ij = contact_neighborhood.e_ij_[n];

            Real p_j = wall_p_[k][index_j];
            Vecd nablaW_ijV_j = contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * e_ij;

            Vecd distance_diff = distance_from_wall - contact_neighborhood.r_ij_[n] * e_ij;
            Real factor = 1.0 - distance_from_wall.dot(distance_diff) / distance_from_wall.squaredNorm();

            pressure_grad -= factor * (p_i - p_j) * nablaW_ijV_j;
        }
    }

    p_grad_[index_i] += pressure_grad;
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
