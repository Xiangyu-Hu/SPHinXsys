#ifndef FLUID_STRUCTURE_INTERACTION_HPP
#define FLUID_STRUCTURE_INTERACTION_HPP

#include "fluid_structure_interaction.h"

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
template <class FluidIntegration2ndHalfType>
PressureForceFromFluid<FluidIntegration2ndHalfType>::
    PressureForceFromFluid(BaseContactRelation &contact_relation)
    : BaseForceFromFluid(contact_relation, "PressureForceFromFluid"),
      vel_ave_(solid_.AverageVelocity(particles_)),
      acc_ave_(solid_.AverageAcceleration(particles_)),
      n_(particles_->getVariableDataByName<Vecd>("NormalDirection"))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_rho_.push_back(contact_particles_[k]->template getVariableDataByName<Real>("Density"));
        contact_mass_.push_back(contact_particles_[k]->template getVariableDataByName<Real>("Mass"));
        contact_vel_.push_back(contact_particles_[k]->template getVariableDataByName<Vecd>("Velocity"));
        contact_Vol_.push_back(contact_particles_[k]->template getVariableDataByName<Real>("VolumetricMeasure"));
        contact_p_.push_back(contact_particles_[k]->template getVariableDataByName<Real>("Pressure"));
        contact_force_prior_.push_back(contact_particles_[k]->template getVariableDataByName<Vecd>("ForcePrior"));
        riemann_solvers_.push_back(RiemannSolverType(*contact_fluids_[k], *contact_fluids_[k]));
    }
}
//=================================================================================================//
template <class FluidIntegration2ndHalfType>
void PressureForceFromFluid<FluidIntegration2ndHalfType>::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real *Vol_k = contact_Vol_[k];
        Real *rho_k = contact_rho_[k];
        Real *mass_k = contact_mass_[k];
        Real *p_k = contact_p_[k];
        Vecd *vel_k = contact_vel_[k];
        Vecd *force_prior_k = contact_force_prior_[k];
        RiemannSolverType &riemann_solvers_k = riemann_solvers_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];
            Real r_ij = contact_neighborhood.r_ij_[n];
            Real face_wall_external_acceleration =
                (force_prior_k[index_j] / mass_k[index_j] - acc_ave_[index_i]).dot(e_ij);
            Real p_j_in_wall = p_k[index_j] + rho_k[index_j] * r_ij * SMAX(Real(0), face_wall_external_acceleration);
            Real u_jump = 2.0 * (vel_k[index_j] - vel_ave_[index_i]).dot(n_[index_i]);
            force -= (riemann_solvers_k.DissipativePJump(u_jump) * n_[index_i] + (p_j_in_wall + p_k[index_j]) * e_ij) *
                     contact_neighborhood.dW_ij_[n] * Vol_k[index_j];
        }
    }
    force_from_fluid_[index_i] = force * Vol_[index_i];
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
#endif // FLUID_STRUCTURE_INTERACTION_HPP