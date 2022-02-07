/**
 * @file 	eulerian_compressible_fluid_dynamics_inner.hpp
 * @author	Zhentong Wang,Chi Zhang and Xiangyu Hu
 */

#include "eulerian_compressible_fluid_dynamics_inner.h"

//=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	namespace eulerian_compressible_fluid_dynamics
	{
		//=================================================================================================//
		template <class RiemannSolverType>
		BasePressureRelaxationInner<RiemannSolverType>::
			BasePressureRelaxationInner(BaseBodyRelationInner &inner_relation)
			: BasePressureRelaxation(inner_relation),
			  riemann_solver_(*material_, *material_) {}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BasePressureRelaxationInner<RiemannSolverType>::Interaction(size_t index_i, Real dt)
		{
			CompressibleFluidState state_i(rho_n_[index_i], vel_n_[index_i], p_[index_i], E_[index_i]);
			Vecd momentum_change_rate = dmom_dt_prior_[index_i];
			Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ij = inner_neighborhood.dW_ij_[n];
				Vecd &e_ij = inner_neighborhood.e_ij_[n];

				CompressibleFluidState state_j(rho_n_[index_j], vel_n_[index_j], p_[index_j], E_[index_j]);
				CompressibleFluidState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
				Vecd vel_star = interface_state.vel_;
				Real p_star = interface_state.p_;
				Real rho_star = interface_state.rho_;

				momentum_change_rate -= 2.0 * Vol_[index_j] *
										(SimTK::outer(rho_star * vel_star, vel_star) + p_star * Matd(1.0)) * e_ij * dW_ij;
			}
			dmom_dt_[index_i] = momentum_change_rate;
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		BaseDensityAndEnergyRelaxationInner<RiemannSolverType>::
			BaseDensityAndEnergyRelaxationInner(BaseBodyRelationInner &inner_relation)
			: BaseDensityAndEnergyRelaxation(inner_relation),
			  riemann_solver_(*material_, *material_) {}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseDensityAndEnergyRelaxationInner<RiemannSolverType>::Interaction(size_t index_i, Real dt)
		{
			CompressibleFluidState state_i(rho_n_[index_i], vel_n_[index_i], p_[index_i], E_[index_i]);
			Real density_change_rate = 0.0;
			Real energy_change_rate = dE_dt_prior_[index_i];
			Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd &e_ij = inner_neighborhood.e_ij_[n];
				Real dW_ij = inner_neighborhood.dW_ij_[n];

				CompressibleFluidState state_j(rho_n_[index_j], vel_n_[index_j], p_[index_j], E_[index_j]);
				CompressibleFluidState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
				//Vecd vel_star = interface_state.get_state_vel();
				Vecd vel_star = interface_state.vel_;
				Real p_star = interface_state.p_;
				Real rho_star = interface_state.rho_;
				Real E_star = interface_state.E_;

				density_change_rate -= 2.0 * Vol_[index_j] * dot(rho_star * vel_star, e_ij) * dW_ij;
				energy_change_rate -= 2.0 * Vol_[index_j] * dot(E_star * vel_star + p_star * vel_star, e_ij) * dW_ij;
			}
			drho_dt_[index_i] = density_change_rate;
			dE_dt_[index_i] = energy_change_rate;
		};
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//