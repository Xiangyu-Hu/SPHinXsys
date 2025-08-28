#ifndef EULERIAN_COMPRESSIBLE_FLUID_INTEGRATION_HPP
#define EULERIAN_COMPRESSIBLE_FLUID_INTEGRATION_HPP

#include "eulerian_compressible_fluid_integration.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class RiemannSolverType>
EulerianCompressibleIntegration1stHalf<RiemannSolverType>::
    EulerianCompressibleIntegration1stHalf(BaseInnerRelation &inner_relation, Real limiter_parameter)
    : BaseIntegrationInCompressible(inner_relation),
      riemann_solver_(compressible_fluid_, compressible_fluid_, limiter_parameter) {}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianCompressibleIntegration1stHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Real energy_per_volume_i = E_[index_i] / Vol_[index_i];
    CompressibleFluidState state_i(rho_[index_i], vel_[index_i], p_[index_i], energy_per_volume_i);
    Vecd momentum_change_rate = force_prior_[index_i];
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];

        Real energy_per_volume_j = E_[index_j] / Vol_[index_j];
        CompressibleFluidState state_j(rho_[index_j], vel_[index_j], p_[index_j], energy_per_volume_j);
        CompressibleFluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
        Matd convect_flux = interface_state.rho_ * interface_state.vel_ * interface_state.vel_.transpose();
        momentum_change_rate -= 2.0 * Vol_[index_i] * dW_ijV_j * (convect_flux + interface_state.p_ * Matd::Identity()) * e_ij;
    }
    force_[index_i] = momentum_change_rate;
}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianCompressibleIntegration1stHalf<RiemannSolverType>::update(size_t index_i, Real dt)
{
    mom_[index_i] += force_[index_i] * dt;
    vel_[index_i] = mom_[index_i] / mass_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType>
EulerianCompressibleIntegration2ndHalf<RiemannSolverType>::
    EulerianCompressibleIntegration2ndHalf(BaseInnerRelation &inner_relation, Real limiter_parameter)
    : BaseIntegrationInCompressible(inner_relation),
      riemann_solver_(compressible_fluid_, compressible_fluid_, limiter_parameter) {}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianCompressibleIntegration2ndHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Real energy_per_volume_i = E_[index_i] / Vol_[index_i];
    CompressibleFluidState state_i(rho_[index_i], vel_[index_i], p_[index_i], energy_per_volume_i);
    Real mass_change_rate = 0.0;
    Real energy_change_rate = force_prior_[index_i].dot(vel_[index_i]); // TODO: not conservative formulation
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];

        Real energy_per_volume_j = E_[index_j] / Vol_[index_j];
        CompressibleFluidState state_j(rho_[index_j], vel_[index_j], p_[index_j], energy_per_volume_j);
        CompressibleFluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);

        mass_change_rate -= 2.0 * Vol_[index_i] * dW_ijV_j * (interface_state.rho_ * interface_state.vel_).dot(e_ij);
        energy_change_rate -= 2.0 * Vol_[index_i] * dW_ijV_j * ((interface_state.E_ + interface_state.p_) * interface_state.vel_).dot(e_ij);
    }
    dmass_dt_[index_i] = mass_change_rate;
    dE_dt_[index_i] = energy_change_rate;
}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianCompressibleIntegration2ndHalf<RiemannSolverType>::update(size_t index_i, Real dt)
{
    E_[index_i] += dE_dt_[index_i] * dt;
    mass_[index_i] += dmass_dt_[index_i] * dt;
    rho_[index_i] = mass_[index_i] / Vol_[index_i];
    Real rho_e = E_[index_i] / Vol_[index_i] - 0.5 * (mom_[index_i] / mass_[index_i]).squaredNorm() * rho_[index_i];
    p_[index_i] = compressible_fluid_.getPressure(rho_[index_i], rho_e);
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // EULERIAN_COMPRESSIBLE_FLUID_INTEGRATION_HPP