/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	eulerian_weakly_compressible_fluid_dynamics_inner.hpp
 * @brief 	Here, we define the algorithm classes for weakly compressible fluid dynamics within the body.
 * @details We consider here weakly compressible fluids.
 *			TODO: It seems that the eulerian and Lagrangian formulation can be merged together
 * @author	Zhentong Wang, Chi ZHang and Xiangyu Hu
 */

#ifndef EULERIAN_WEAKLY_COMPRESSIBLE_FLUID_DYNAMICS_INNER_HPP
#define EULERIAN_WEAKLY_COMPRESSIBLE_FLUID_DYNAMICS_INNER_HPP

#include "eulerian_weakly_compressible_fluid_dynamics_inner.h"

//=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	namespace eulerian_weakly_compressible_fluid_dynamics
	{
		//=================================================================================================//
		template <class RiemannSolverType>
		BaseIntegration1stHalf<RiemannSolverType>::BaseIntegration1stHalf(BaseInnerRelation &inner_relation)
			: BaseIntegration(inner_relation), riemann_solver_(this->fluid_, this->fluid_) {}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration1stHalf<RiemannSolverType>::initialization(size_t index_i, Real dt)
		{
			rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
			p_[index_i] = fluid_.getPressure(rho_[index_i]);
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration1stHalf<RiemannSolverType>::update(size_t index_i, Real dt)
		{
			mom_[index_i] += dmom_dt_[index_i] * dt;
			vel_[index_i] = mom_[index_i] / rho_[index_i];
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration1stHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
		{
			FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
			Vecd momentum_change_rate = dmom_dt_prior_[index_i];
			Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
				Vecd &e_ij = inner_neighborhood.e_ij_[n];

				FluidState state_j(rho_[index_j], vel_[index_j], p_[index_j]);
				FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
				Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);

				momentum_change_rate -= 2.0 *
					((rho_star * interface_state.vel_) * interface_state.vel_.transpose() + interface_state.p_ * Matd::Identity()) * e_ij * dW_ijV_j;
			}
			dmom_dt_[index_i] = momentum_change_rate;
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		BaseIntegration2ndHalf<RiemannSolverType>::BaseIntegration2ndHalf(BaseInnerRelation &inner_relation)
			: BaseIntegration(inner_relation), riemann_solver_(fluid_, fluid_) {}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration2ndHalf<RiemannSolverType>::update(size_t index_i, Real dt)
		{
			rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration2ndHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
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
		};
		//=================================================================================================//
	}
	//=================================================================================================//
}
#endif // EULERIAN_WEAKLY_COMPRESSIBLE_FLUID_DYNAMICS_INNER_HPP
