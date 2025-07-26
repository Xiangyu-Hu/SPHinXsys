#ifndef TURBULENCEMODEL_HPP
#define TURBULENCEMODEL_HPP
#include "turbulence_model.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class RiemannSolverType>
KEpsilonStd1stHalf<RiemannSolverType>::KEpsilonStd1stHalf(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator, Real limiter_parameter)
    : StdWallFunctionFVM(inner_relation, ghost_creator),
      dK_dt_(this->particles_->template registerStateVariableData<Real>("TKEChangeRate")),
      wall_adjacent_cell_flag_(this->particles_->template getVariableDataByName<Real>("FlagForWallAdjacentCells")),
      strain_rate_(this->particles_->template registerStateVariableData<Real>("StrainRate")),
      riemann_solver_(this->fluid_, this->fluid_, limiter_parameter)
{
}
//=================================================================================================//
template <class RiemannSolverType>
void KEpsilonStd1stHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    ExendedFluidState state_i(rho_[index_i], vel_[index_i], p_[index_i], K_[index_i], Eps_[index_i]);
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    Matd K_prod = Matd::Zero(), vel_matrix = Matd::Zero();
    Matd strain_tensor = Matd::Zero(), strain_rate_modulus = Matd::Zero();
    K_prod_[index_i] = 0.0, strain_rate_[index_i] = 0.0;
    Real K_adv = 0.0, K_lap = 0.0;
    vel_gradient_mat_[index_i] = Matd::Zero();
    mu_t_[index_i] = rho_[index_i] * C_mu_ * ((K_[index_i] * K_[index_i]) / (Eps_[index_i]));

    if (wall_adjacent_cell_flag_[index_i] == 1.0)
    {
        nearwallquantities(index_i);
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Real dW_ij = inner_neighborhood.dW_ij_[n];
            Real r_ij = inner_neighborhood.r_ij_[n];
            Vecd &e_ij = inner_neighborhood.e_ij_[n];
            ExendedFluidState state_j(rho_[index_j], vel_[index_j], p_[index_j], K_[index_j], Eps_[index_j]);
            ExtendedFluidStarState interface_state = riemann_solver_.getExtendedInterfaceState(state_i, state_j, e_ij);

            Real mu_t_avg = (2.0 * mu_t_[index_i] * mu_t_[index_j]) / (mu_t_[index_i] + mu_t_[index_j]);

            K_adv += -2.0 * (dW_ij * Vol_[index_j] * interface_state.rho_) * (interface_state.K_) * (interface_state.vel_.dot(e_ij));
            K_lap += 2.0 * dW_ij * Vol_[index_j] * ((viscosity_.ReferenceViscosity() + mu_t_avg / sigma_k_) * (K_[index_i] - K_[index_j]) / (r_ij));
        }
        strain_tensor = 0.5 * (vel_gradient_mat_[index_i] + vel_gradient_mat_[index_i].transpose());
        strain_rate_modulus = 2.0 * strain_tensor.array() * strain_tensor.array();
        strain_rate_[index_i] = std::sqrt(strain_rate_modulus.sum());

        K_prod_[index_i] = K_prod_p_[index_i];
        Eps_[index_i] = Eps_p_[index_i];

        dK_dt_[index_i] = K_adv + K_prod_[index_i] - rho_[index_i] * Eps_[index_i] + K_lap;
    }
    else
    {
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Real dW_ij = inner_neighborhood.dW_ij_[n];
            Real r_ij = inner_neighborhood.r_ij_[n];
            Vecd &e_ij = inner_neighborhood.e_ij_[n];
            ExendedFluidState state_j(rho_[index_j], vel_[index_j], p_[index_j], K_[index_j], Eps_[index_j]);
            ExtendedFluidStarState interface_state = riemann_solver_.getExtendedInterfaceState(state_i, state_j, e_ij);

            Real mu_t_avg = (2.0 * mu_t_[index_i] * mu_t_[index_j]) / (mu_t_[index_i] + mu_t_[index_j]);

            K_adv += -2.0 * (dW_ij * Vol_[index_j] * interface_state.rho_) * (interface_state.K_) * (interface_state.vel_.dot(e_ij));
            K_lap += 2.0 * dW_ij * Vol_[index_j] * ((viscosity_.ReferenceViscosity() + mu_t_avg / sigma_k_) * (K_[index_i] - K_[index_j]) / (r_ij));

            vel_matrix = (vel_[index_i] - vel_[index_j]) * e_ij.transpose();
            vel_gradient_mat_[index_i] += dW_ij * Vol_[index_j] * vel_matrix;
        }
        strain_tensor = 0.5 * (vel_gradient_mat_[index_i] + vel_gradient_mat_[index_i].transpose());
        strain_rate_modulus = 2.0 * strain_tensor.array() * strain_tensor.array();
        strain_rate_[index_i] = std::sqrt(strain_rate_modulus.sum());

        K_prod = (mu_t_[index_i] * strain_rate_modulus);
        K_prod_[index_i] = K_prod.sum();

        dK_dt_[index_i] = K_adv + K_prod_[index_i] - rho_[index_i] * Eps_[index_i] + K_lap;
    }
}
//=================================================================================================//
template <class RiemannSolverType>
void KEpsilonStd1stHalf<RiemannSolverType>::update(size_t index_i, Real dt)
{
    K_[index_i] += (dK_dt_[index_i] / rho_[index_i]) * dt;
}
//=================================================================================================//
template <class RiemannSolverType>
KEpsilonStd2ndHalf<RiemannSolverType>::KEpsilonStd2ndHalf(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator, Real limiter_parameter)
    : BaseTurbulence(inner_relation, ghost_creator),
      dEps_dt_(this->particles_->template registerStateVariableData<Real>("DissipationChangeRate")),
      wall_adjacent_cell_flag_(this->particles_->template getVariableDataByName<Real>("FlagForWallAdjacentCells")),
      riemann_solver_(this->fluid_, this->fluid_),
      vel_gradient_mat_(this->particles_->template getVariableDataByName<Matd>("VelocityGradient"))
{
}
//=================================================================================================//
template <class RiemannSolverType>
void KEpsilonStd2ndHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    ExendedFluidState state_i(rho_[index_i], vel_[index_i], p_[index_i], K_[index_i], Eps_[index_i]);
    Real Eps_changerate = 0.0;
    Real Eps_adv = 0.0, Eps_lap = 0.0, Eps_prod = 0.0, Eps_destruction = 0.0;
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    if (wall_adjacent_cell_flag_[index_i] != 1)
    {
        mu_t_[index_i] = rho_[index_i] * C_mu_ * ((K_[index_i] * K_[index_i]) / (Eps_[index_i]));

        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Real dW_ij = inner_neighborhood.dW_ij_[n];
            Real r_ij = inner_neighborhood.r_ij_[n];
            Vecd &e_ij = inner_neighborhood.e_ij_[n];
            ExendedFluidState state_j(rho_[index_j], vel_[index_j], p_[index_j], K_[index_j], Eps_[index_j]);
            ExtendedFluidStarState interface_state = riemann_solver_.getExtendedInterfaceState(state_i, state_j, e_ij);

            Real mu_t_avg = (2.0 * mu_t_[index_i] * mu_t_[index_j]) / (mu_t_[index_i] + mu_t_[index_j]);

            Eps_adv += -2.0 * (dW_ij * Vol_[index_j] * interface_state.rho_) * (interface_state.Eps_) * (interface_state.vel_.dot(e_ij));
            Eps_lap += 2.0 * dW_ij * Vol_[index_j] * (viscosity_.ReferenceViscosity() + mu_t_avg / sigma_eps_) * ((Eps_[index_i] - Eps_[index_j]) / (r_ij));

            Eps_changerate = Eps_adv + Eps_lap;
        }
        Eps_prod = C1_eps_ * (Eps_[index_i] / (K_[index_i])) * K_prod_[index_i];
        Eps_destruction = -C2_eps_ * rho_[index_i] * (Eps_[index_i] * Eps_[index_i]) / (K_[index_i]);
        dEps_dt_[index_i] = Eps_changerate + Eps_prod + Eps_destruction;
    }
}
//=================================================================================================//
template <class RiemannSolverType>
void KEpsilonStd2ndHalf<RiemannSolverType>::update(size_t index_i, Real dt)
{
    if (wall_adjacent_cell_flag_[index_i] != 1)
    {
        Eps_[index_i] += (dEps_dt_[index_i] / rho_[index_i]) * dt;
    }
}
} // namespace fluid_dynamics

} // namespace SPH
#endif // TURBULENCEMODEL_HPP