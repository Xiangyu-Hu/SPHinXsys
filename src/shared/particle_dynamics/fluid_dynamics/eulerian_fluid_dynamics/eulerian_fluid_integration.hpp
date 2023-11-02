#ifndef EULERIAN_FLUID_INTEGRATION_HPP
#define EULERIAN_FLUID_INTEGRATION_HPP

#include "eulerian_fluid_integration.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
EulerianIntegration<DataDelegationType>::EulerianIntegration(BaseRelationType &base_relation)
    : BaseIntegration<DataDelegationType>(base_relation),
      mom_(*this->particles_->template registerSharedVariable<Vecd>("Momentum")),
      dmom_dt_(*this->particles_->template registerSharedVariable<Vecd>("MomentumChangeRate")) {}
//=================================================================================================//
template <class RiemannSolverType>
EulerianIntegration1stHalf<Inner<>, RiemannSolverType>::
    EulerianIntegration1stHalf(BaseInnerRelation &inner_relation, Real limiter_parameter)
    : EulerianIntegration<FluidDataInner>(inner_relation),
      riemann_solver_(this->fluid_, this->fluid_, limiter_parameter) {}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration1stHalf<Inner<>, RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
    Vecd momentum_change_rate = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];

        FluidState state_j(rho_[index_j], vel_[index_j], p_[index_j]);
        FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
        Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);
        Matd convect_flux = rho_star * interface_state.vel_ * interface_state.vel_.transpose();
        momentum_change_rate -= 2.0 * (convect_flux + interface_state.p_ * Matd::Identity()) * e_ij * dW_ijV_j;
    }
    dmom_dt_[index_i] = momentum_change_rate;
}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration1stHalf<Inner<>, RiemannSolverType>::update(size_t index_i, Real dt)
{
    mom_[index_i] += (dmom_dt_[index_i] + rho_[index_i] * acc_prior_[index_i]) * dt;
    vel_[index_i] = mom_[index_i] / rho_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType>
EulerianIntegration1stHalf<ContactWall<>, RiemannSolverType>::
    EulerianIntegration1stHalf(BaseContactRelation &wall_contact_relation, Real limiter_parameter)
    : InteractionWithWall<EulerianIntegration>(wall_contact_relation),
      riemann_solver_(fluid_, fluid_, limiter_parameter) {}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration1stHalf<ContactWall<>, RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
    Vecd momentum_change_rate = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &n_k = *(wall_n_[k]);
        Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];

            Vecd vel_in_wall = -state_i.vel_;
            Real p_in_wall = state_i.p_;
            Real rho_in_wall = state_i.rho_;
            FluidState state_j(rho_in_wall, vel_in_wall, p_in_wall);
            FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, n_k[index_j]);
            Real rho_star = fluid_.DensityFromPressure(interface_state.p_);
            Matd convect_flux = rho_star * interface_state.vel_ * interface_state.vel_.transpose();
            momentum_change_rate -= 2.0 * (convect_flux + interface_state.p_ * Matd::Identity()) * e_ij * dW_ijV_j;
        }
    }
    dmom_dt_[index_i] += momentum_change_rate;
}
//=================================================================================================//
template <class RiemannSolverType>
EulerianIntegration2ndHalf<Inner<>, RiemannSolverType>::
    EulerianIntegration2ndHalf(BaseInnerRelation &inner_relation, Real limiter_parameter)
    : EulerianIntegration<FluidDataInner>(inner_relation),
      riemann_solver_(this->fluid_, this->fluid_, limiter_parameter) {}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration2ndHalf<Inner<>, RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
    Real density_change_rate = 0.0;
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];

        FluidState state_j(rho_[index_j], vel_[index_j], p_[index_j]);
        FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);

        Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);
        density_change_rate -= 2.0 * (rho_star * interface_state.vel_).dot(e_ij) * dW_ijV_j;
    }
    drho_dt_[index_i] = density_change_rate;
}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration2ndHalf<Inner<>, RiemannSolverType>::update(size_t index_i, Real dt)
{
    rho_[index_i] += drho_dt_[index_i] * dt;
    p_[index_i] = fluid_.getPressure(rho_[index_i]);
}
//=================================================================================================//
template <class RiemannSolverType>
EulerianIntegration2ndHalf<ContactWall<>, RiemannSolverType>::
    EulerianIntegration2ndHalf(BaseContactRelation &wall_contact_relation, Real limiter_parameter)
    : InteractionWithWall<EulerianIntegration>(wall_contact_relation),
      riemann_solver_(this->fluid_, this->fluid_, limiter_parameter){};
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration2ndHalf<ContactWall<>, RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    FluidState state_i(this->rho_[index_i], this->vel_[index_i], this->p_[index_i]);
    Real density_change_rate = 0.0;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &n_k = *(this->wall_n_[k]);
        Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];

            Vecd vel_in_wall = -state_i.vel_;
            Real p_in_wall = state_i.p_;
            Real rho_in_wall = state_i.rho_;

            FluidState state_j(rho_in_wall, vel_in_wall, p_in_wall);
            FluidStarState interface_state = this->riemann_solver_.getInterfaceState(state_i, state_j, n_k[index_j]);
            Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);

            density_change_rate -= 2.0 * (rho_star * interface_state.vel_).dot(e_ij) * dW_ijV_j;
        }
    }
    this->drho_dt_[index_i] += density_change_rate;
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // EULERIAN_FLUID_INTEGRATION_HPP