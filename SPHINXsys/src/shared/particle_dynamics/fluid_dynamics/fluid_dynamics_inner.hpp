/**
 * @file 	fluid_dynamics_inner.hpp
 * @author	Chi ZHang and Xiangyu Hu
 */

#pragma once

#include "fluid_dynamics_inner.h"

//=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		template <class RiemannSolverType>
		BasePressureRelaxationInner<RiemannSolverType>::
			BasePressureRelaxationInner(BaseBodyRelationInner &inner_relation)
			: BasePressureRelaxation(inner_relation),
			  riemann_solver_(fluid_, fluid_) {}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BasePressureRelaxationInner<RiemannSolverType>::interaction(size_t index_i, Real dt)
		{
			FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
			Vecd acceleration = acc_prior_[index_i];
			Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ij = inner_neighborhood.dW_ij_[n];
				Vecd &e_ij = inner_neighborhood.e_ij_[n];

				FluidState state_j(rho_[index_j], vel_[index_j], p_[index_j]);
				Real p_star = riemann_solver_.getPStar(state_i, state_j, e_ij);
				acceleration -= 2.0 * p_star * Vol_[index_j] * dW_ij * e_ij / state_i.rho_;
			}
			acc_[index_i] = acceleration;
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		BaseDensityRelaxationInner<RiemannSolverType>::
			BaseDensityRelaxationInner(BaseBodyRelationInner &inner_relation)
			: BaseDensityRelaxation(inner_relation),
			  riemann_solver_(fluid_, fluid_) {}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseDensityRelaxationInner<RiemannSolverType>::interaction(size_t index_i, Real dt)
		{
			FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
			Real density_change_rate = 0.0;
			Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd &e_ij = inner_neighborhood.e_ij_[n];
				Real dW_ij = inner_neighborhood.dW_ij_[n];

				FluidState state_j(rho_[index_j], vel_[index_j], p_[index_j]);
				Vecd vel_star = riemann_solver_.getVStar(state_i, state_j, e_ij);
				density_change_rate += 2.0 * state_i.rho_ * Vol_[index_j] * dot(state_i.vel_ - vel_star, e_ij) * dW_ij;
			}
			drho_dt_[index_i] = density_change_rate;
		};
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//