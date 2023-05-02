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

#include "fluid_dynamics_complex_wkgc.h"

namespace SPH
{
	//=====================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		template <class BaseIntegration1stHalfCorrectType>
		void BaseIntegration1stHalfCorrectWithWall<BaseIntegration1stHalfCorrectType>::interaction(size_t index_i, Real dt)
		{
			BaseIntegration1stHalfCorrectType::interaction(index_i, dt);

			Vecd acc_prior_i = this->acc_prior_[index_i];

			Vecd acceleration = Vecd::Zero();
			Real rho_dissipation(0);
			for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd>& acc_ave_k = *(this->wall_acc_ave_[k]);
				Neighborhood& wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
				{
					size_t index_j = wall_neighborhood.j_[n];
					Vecd& e_ij = wall_neighborhood.e_ij_[n];
					Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];
					Real r_ij = wall_neighborhood.r_ij_[n];

					Real face_wall_external_acceleration = (acc_prior_i - acc_ave_k[index_j]).dot(-e_ij);
					Real p_in_wall = this->p_[index_i] + this->rho_[index_i] * r_ij * SMAX(0.0, face_wall_external_acceleration);
					acceleration -= (this->p_[index_i] + p_in_wall) * this->B_[index_i] * e_ij * dW_ijV_j;
					rho_dissipation += this->riemann_solver_.DissipativeUJump(this->p_[index_i] - p_in_wall) * dW_ijV_j;
				}
			}
			this->acc_[index_i] += acceleration / this->rho_[index_i];
			this->drho_dt_[index_i] += rho_dissipation * this->rho_[index_i];
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//