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
      mom_(this->particles_->template registerStateVariable<Vecd>("Momentum")),
      dmom_dt_(this->particles_->template registerStateVariable<Vecd>("MomentumChangeRate")),
      dmass_dt_(this->particles_->template registerStateVariable<Real>("MassChangeRate")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")) {}
//=================================================================================================//
template <class RiemannSolverType>
EulerianIntegration1stHalf<Inner<>, RiemannSolverType>::
    EulerianIntegration1stHalf(BaseInnerRelation &inner_relation, Real limiter_parameter)
    : EulerianIntegration<DataDelegateInner>(inner_relation),
      riemann_solver_(this->fluid_, this->fluid_, limiter_parameter) {}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration1stHalf<Inner<>, RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    FluidStateIn state_i(rho_[index_i], vel_[index_i], p_[index_i]);
    Vecd momentum_change_rate = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];

        FluidStateIn state_j(rho_[index_j], vel_[index_j], p_[index_j]);
        FluidStateOut interface_state = riemann_solver_.InterfaceState(state_i, state_j, e_ij);
        Matd convect_flux = interface_state.rho_ * interface_state.vel_ * interface_state.vel_.transpose();
        
        momentum_change_rate -= 2.0 * Vol_[index_i] * (convect_flux + interface_state.p_ * Matd::Identity()) * e_ij * dW_ijV_j;
    }
    dmom_dt_[index_i] = momentum_change_rate;
}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration1stHalf<Inner<>, RiemannSolverType>::update(size_t index_i, Real dt)
{
    mom_[index_i] += (dmom_dt_[index_i] + force_prior_[index_i]) * dt;
    vel_[index_i] = mom_[index_i] / mass_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType>
EulerianIntegration1stHalf<Contact<Wall>, RiemannSolverType>::
    EulerianIntegration1stHalf(BaseContactRelation &wall_contact_relation, Real limiter_parameter)
    : BaseEulerianIntegrationWithWall(wall_contact_relation),
      riemann_solver_(fluid_, fluid_, limiter_parameter) {}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration1stHalf<Contact<Wall>, RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    FluidStateIn state_i(rho_[index_i], vel_[index_i], p_[index_i]);
    Vecd momentum_change_rate = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Vecd *n_k = wall_n_[k];
        Real *Vol_k = wall_Vol_[k];
        Vecd *vel_ave_k = wall_vel_ave_[k];
        Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ij_[n] * Vol_k[index_j];

            Vecd vel_j_in_wall = 2.0 * vel_ave_k[index_j] - state_i.vel_;
            Real p_j_in_wall = state_i.p_;
            Real rho_in_wall = state_i.rho_;
            FluidStateIn state_j(rho_in_wall, vel_j_in_wall, p_j_in_wall);
            FluidStateOut interface_state = riemann_solver_.InterfaceState(state_i, state_j, n_k[index_j]);
            Matd convect_flux = interface_state.rho_ * interface_state.vel_ * interface_state.vel_.transpose();
            momentum_change_rate -= 2.0 * Vol_[index_i] * (convect_flux + interface_state.p_ * Matd::Identity()) * e_ij * dW_ijV_j;
        }
    }
    dmom_dt_[index_i] += momentum_change_rate;
}
//=================================================================================================//
template <class RiemannSolverType>
EulerianIntegration2ndHalf<Inner<>, RiemannSolverType>::
    EulerianIntegration2ndHalf(BaseInnerRelation &inner_relation, Real limiter_parameter)
    : EulerianIntegration<DataDelegateInner>(inner_relation),
      riemann_solver_(this->fluid_, this->fluid_, limiter_parameter) {}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration2ndHalf<Inner<>, RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    FluidStateIn state_i(rho_[index_i], vel_[index_i], p_[index_i]);
    Real mass_change_rate = 0.0;
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];

        FluidStateIn state_j(rho_[index_j], vel_[index_j], p_[index_j]);
        FluidStateOut interface_state = riemann_solver_.InterfaceState(state_i, state_j, e_ij);
        
        mass_change_rate -= 2.0 * Vol_[index_i] * (interface_state.rho_ * interface_state.vel_).dot(e_ij) * dW_ijV_j;
    }
    dmass_dt_[index_i] = mass_change_rate;
}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration2ndHalf<Inner<>, RiemannSolverType>::update(size_t index_i, Real dt)
{
    mass_[index_i] += dmass_dt_[index_i] * dt;
    rho_[index_i] = mass_[index_i] / Vol_[index_i];
    p_[index_i] = fluid_.getPressure(rho_[index_i]);
}
//=================================================================================================//
template <class RiemannSolverType>
EulerianIntegration2ndHalf<Contact<Wall>, RiemannSolverType>::
    EulerianIntegration2ndHalf(BaseContactRelation &wall_contact_relation, Real limiter_parameter)
    : BaseEulerianIntegrationWithWall(wall_contact_relation),
      riemann_solver_(this->fluid_, this->fluid_, limiter_parameter){};
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration2ndHalf<Contact<Wall>, RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    FluidStateIn state_i(this->rho_[index_i], this->vel_[index_i], this->p_[index_i]);
    Real mass_change_rate = 0.0;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Vecd *n_k = this->wall_n_[k];
        Real *Vol_k = this->wall_Vol_[k];
        Vecd *vel_ave_k = wall_vel_ave_[k];
        Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ij_[n] * Vol_k[index_j];

            Vecd vel_j_in_wall = 2.0 * vel_ave_k[index_j] - state_i.vel_;
            Real p_j_in_wall = state_i.p_;
            Real rho_in_wall = state_i.rho_;

            FluidStateIn state_j(rho_in_wall, vel_j_in_wall, p_j_in_wall);
            FluidStateOut interface_state = this->riemann_solver_.InterfaceState(state_i, state_j, n_k[index_j]);
            mass_change_rate -= 2.0 * this->Vol_[index_i] * (interface_state.rho_ * interface_state.vel_).dot(e_ij) * dW_ijV_j;
        }
    }
    this->dmass_dt_[index_i] += mass_change_rate;
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // EULERIAN_FLUID_INTEGRATION_HPP
