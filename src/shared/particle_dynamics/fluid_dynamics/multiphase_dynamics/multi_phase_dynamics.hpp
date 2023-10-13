#pragma once

#include "multi_phase_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType>
Integration1stHalf<Contact<>, RiemannSolverType, KernelCorrectionType>::
    Integration1stHalf(BaseContactRelation &contact_relation)
    : BaseIntegration<FluidContactData>(contact_relation),
      correction_(particles_)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_corrections_.push_back(KernelCorrectionType(contact_particles_[k]));
        Fluid &contact_fluid = DynamicCast<Fluid>(this, contact_particles_[k]->getBaseMaterial());
        riemann_solvers_.push_back(RiemannSolverType(fluid_, contact_fluid));
        contact_p_.push_back(contact_particles_[k]->template getVariableByName<Real>("Pressure"));
    }
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType>
void MultiPhaseMomentumInterface<RiemannSolverType, KernelCorrectionType>::interaction(size_t index_i, Real dt)
{
    Vecd acceleration = Vecd::Zero();
    Real rho_dissipation(0);
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        StdLargeVec<Real> &p_k = *(this->contact_p_[k]);
        KernelCorrectionType &correction_k = contact_corrections_[k];
        RiemannSolverType &riemann_solver_k = riemann_solvers_[k];
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd &e_ij = contact_neighborhood.e_ij_[n];
            Real dW_ijV_j = contact_neighborhood.dW_ijV_j_[n];

            acceleration -= riemann_solver_k.AverageP(this->p_[index_i] * correction_k(index_j), p_k[index_j] * correction_(index_i)) *
                            2.0 * e_ij * dW_ijV_j;
            rho_dissipation += riemann_solver_k.DissipativeUJump(this->p_[index_i] - p_k[index_j]) * dW_ijV_j;
        }
    }
    this->acc_[index_i] += acceleration / this->rho_[index_i];
    this->drho_dt_[index_i] += rho_dissipation * this->rho_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType>
MultiPhaseContinuityInterface<RiemannSolverType>::
    MultiPhaseContinuityInterface(BaseContactRelation &contact_relation)
    : BaseIntegration<MultiPhaseContactData>(contact_relation)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        Fluid &contact_fluid = DynamicCast<Fluid>(this, contact_particles_[k]->getBaseMaterial());
        riemann_solvers_.push_back(RiemannSolverType(fluid_, contact_fluid));
        contact_vel_.push_back(contact_particles_[k]->template getVariableByName<Vecd>("Velocity"));
    }
}
//=================================================================================================//
template <class RiemannSolverType>
void MultiPhaseContinuityInterface<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Real density_change_rate = 0.0;
    Vecd p_dissipation = Vecd::Zero();
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &vel_k = *(this->contact_vel_[k]);
        RiemannSolverType &riemann_solver_k = riemann_solvers_[k];
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd &e_ij = contact_neighborhood.e_ij_[n];
            Real dW_ijV_j = contact_neighborhood.dW_ijV_j_[n];

            Vecd vel_ave = riemann_solver_k.AverageV(this->vel_[index_i], vel_k[index_j]);
            density_change_rate += 2.0 * (this->vel_[index_i] - vel_ave).dot(e_ij) * dW_ijV_j;
            Real u_jump = (this->vel_[index_i] - vel_k[index_j]).dot(e_ij);
            p_dissipation += riemann_solver_k.DissipativePJump(u_jump) * dW_ijV_j * e_ij;
        }
    }
    this->drho_dt_[index_i] += density_change_rate * this->rho_[index_i];
    this->acc_[index_i] += p_dissipation / this->rho_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
