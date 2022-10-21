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
 *  HU1527/12-1 and Hu1527/12-4												*
 *                                                                          *
 * Portions copyright (c) 2017-2020 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file    eulerian_compressible_fluid_dynamics_inner.h
 * @brief 	Here, we define the algorithm classes for eulerian fluid dynamics within the body.
 * @details We consider here compressible fluids.
 * @author	Zhentong Wang, Chi ZHang and Xiangyu Hu
 * @version	1.0
 *			Try to implement EIGEN libaary for base vector, matrix and 
 *			linear algebra operation.  
 *			-- Chi ZHANG
 */

#ifndef EULERIAN_COMPRESSIBLE_FLUID_DYNAMICS_INNER_HPP
#define EULERIAN_COMPRESSIBLE_FLUID_DYNAMICS_INNER_HPP

#include "eulerian_compressible_fluid_dynamics_inner.h"

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
		void BasePressureRelaxationInner<RiemannSolverType>::interaction(size_t index_i, Real dt)
		{
			CompressibleFluidState state_i(rho_[index_i], vel_[index_i], p_[index_i], E_[index_i]);
			Vecd momentum_change_rate = dmom_dt_prior_[index_i];
			Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ij = inner_neighborhood.dW_ij_[n];
				Vecd &e_ij = inner_neighborhood.e_ij_[n];

				CompressibleFluidState state_j(rho_[index_j], vel_[index_j], p_[index_j], E_[index_j]);
				CompressibleFluidState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
				Vecd vel_star = interface_state.vel_;
				Real p_star = interface_state.p_;
				Real rho_star = interface_state.rho_;

				momentum_change_rate -= 2.0 * Vol_[index_j] * dW_ij *
										((rho_star * vel_star) * vel_star.transpose() + p_star * Matd::Identity()) * e_ij;
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
		void BaseDensityAndEnergyRelaxationInner<RiemannSolverType>::interaction(size_t index_i, Real dt)
		{
			CompressibleFluidState state_i(rho_[index_i], vel_[index_i], p_[index_i], E_[index_i]);
			Real density_change_rate = 0.0;
			Real energy_change_rate = dE_dt_prior_[index_i];
			Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd &e_ij = inner_neighborhood.e_ij_[n];
				Real dW_ij = inner_neighborhood.dW_ij_[n];

				CompressibleFluidState state_j(rho_[index_j], vel_[index_j], p_[index_j], E_[index_j]);
				CompressibleFluidState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
				//Vecd vel_star = interface_state.get_state_vel();
				Vecd vel_star = interface_state.vel_;
				Real p_star = interface_state.p_;
				Real rho_star = interface_state.rho_;
				Real E_star = interface_state.E_;

				density_change_rate -= 2.0 * Vol_[index_j] * dW_ij * (rho_star * vel_star).dot(e_ij);
				energy_change_rate -= 2.0 * Vol_[index_j] * dW_ij * (E_star * vel_star + p_star * vel_star).dot(e_ij);
			}
			drho_dt_[index_i] = density_change_rate;
			dE_dt_[index_i] = energy_change_rate;
		};
		//=================================================================================================//
	}
	//=================================================================================================//
}
#endif //EULERIAN_COMPRESSIBLE_FLUID_DYNAMICS_INNER_HPP
//=================================================================================================//