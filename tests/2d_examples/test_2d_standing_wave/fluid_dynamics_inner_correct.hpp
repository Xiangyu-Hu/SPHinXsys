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
#pragma once

#include "fluid_dynamics_inner_correct.h"


namespace SPH
{
	//=====================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		template <class RiemannSolverType>
		BaseIntegration1stHalfCorrect<RiemannSolverType>::BaseIntegration1stHalfCorrect(BaseInnerRelation& inner_relation)
			: BaseIntegration1stHalf<RiemannSolverType>(inner_relation),
			B_(*this->particles_->getVariableByName<Matd>("WeightedCorrectionMatrix"))
		{
			this->particles_->registerVariable(p_B_, "CorrectedPressure");
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration1stHalfCorrect<RiemannSolverType>::initialization(size_t index_i, Real dt)
		{
			BaseIntegration1stHalf<RiemannSolverType>::initialization(index_i, dt);
			p_B_[index_i] = this->p_[index_i] * this->B_[index_i];
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration1stHalfCorrect<RiemannSolverType>::update(size_t index_i, Real dt)
		{
			BaseIntegration1stHalf<RiemannSolverType>::update(index_i, dt);
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		Vecd BaseIntegration1stHalfCorrect<RiemannSolverType>::computeNonConservativeAcceleration(size_t index_i)
		{
			Vecd acceleration = this->acc_prior_[index_i] * this->rho_[index_i];
			const Neighborhood& inner_neighborhood = this->inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
				const Vecd& e_ij = inner_neighborhood.e_ij_[n];

				acceleration += (p_B_[index_i] - p_B_[index_j]) * dW_ijV_j * e_ij;
			}
			return acceleration / this->rho_[index_i];
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration1stHalfCorrect<RiemannSolverType>::interaction(size_t index_i, Real dt)
		{
			Vecd acceleration = Vecd::Zero();
			Real rho_dissipation(0);
			const Neighborhood& inner_neighborhood = this->inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
				const Vecd& e_ij = inner_neighborhood.e_ij_[n];

				acceleration -= (p_B_[index_i] + p_B_[index_j]) * dW_ijV_j * e_ij;
				rho_dissipation += this->riemann_solver_.DissipativeUJump(this->p_[index_i] - this->p_[index_j]) * dW_ijV_j;
			}
			this->acc_[index_i] += acceleration / this->rho_[index_i];
			this->drho_dt_[index_i] = rho_dissipation * this->rho_[index_i];
		}
	}
	//=================================================================================================//
}
//=================================================================================================//