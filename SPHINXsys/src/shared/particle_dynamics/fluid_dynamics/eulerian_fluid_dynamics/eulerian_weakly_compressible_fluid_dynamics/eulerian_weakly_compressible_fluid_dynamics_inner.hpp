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
 * @file 	eulerian_weakly_compressible_fluid_dynamics_inner.hpp
 * @brief 	Here, we define the algorithm classes for weakly compressible fluid dynamics within the body.
 * @details We consider here weakly compressible fluids.
 *			TODO: It seems that the eulerian and Lagrangian formulation can be merged together
 * @author	Zhentong Wang, Chi ZHang and Xiangyu Hu
 * @version	1.0
 *			Try to implement EIGEN libaary for base vector, matrix and 
 *			linear algebra operation.  
 *			-- Chi ZHANG
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
		template<class RiemannSolverType>
		BasePressureRelaxationInner<RiemannSolverType>::
			BasePressureRelaxationInner(BaseBodyRelationInner &inner_relation) :
			BasePressureRelaxation(inner_relation),
			riemann_solver_(*material_, *material_) {}
		//=================================================================================================//
		template<class RiemannSolverType>
		void BasePressureRelaxationInner<RiemannSolverType>::interaction(size_t index_i, Real dt)
		{
			FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
			Vecd momentum_change_rate = dmom_dt_prior_[index_i];
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ij = inner_neighborhood.dW_ij_[n];
				Vecd& e_ij = inner_neighborhood.e_ij_[n];

				FluidState state_j(rho_[index_j], vel_[index_j], p_[index_j]);
				FluidState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
				Real p_star = interface_state.p_;
				Vecd vel_star = interface_state.vel_;
				Real rho_star = material_->DensityFromPressure(p_star);

				momentum_change_rate -= 2.0 * Vol_[index_j] * dW_ij *
					(rho_star * vel_star * vel_star.transpose() + p_star * Matd::Identity()) * e_ij;
			}
			dmom_dt_[index_i] = momentum_change_rate;
		}
		//=================================================================================================//
		template<class RiemannSolverType>
		BaseDensityAndEnergyRelaxationInner<RiemannSolverType>::
			BaseDensityAndEnergyRelaxationInner(BaseBodyRelationInner &inner_relation) :
			BaseDensityAndEnergyRelaxation(inner_relation),
			riemann_solver_(*material_, *material_) {}
		//=================================================================================================//
		template<class RiemannSolverType>
		void BaseDensityAndEnergyRelaxationInner<RiemannSolverType>::interaction(size_t index_i, Real dt)
		{
			FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
			Real density_change_rate = 0.0;
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd& e_ij = inner_neighborhood.e_ij_[n];
				Real dW_ij = inner_neighborhood.dW_ij_[n];

				FluidState state_j(rho_[index_j], vel_[index_j], p_[index_j]);
				FluidState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
				Real p_star = interface_state.p_;
				Vecd vel_star = interface_state.vel_;
				Real rho_star = material_->DensityFromPressure(p_star);

				density_change_rate -= 2.0 * Vol_[index_j] * dW_ij * rho_star * vel_star.dot(e_ij);
			}
			drho_dt_[index_i] = density_change_rate;
		};
		//=================================================================================================//
	}
	//=================================================================================================//
}
#endif //EULERIAN_WEAKLY_COMPRESSIBLE_FLUID_DYNAMICS_INNER_HPP
//=================================================================================================//