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
 * @file 	k-epsilon_turbulent_model.hpp
 * @brief 	Here.
 * @author	 Xiangyu Hu
 */
#pragma once

#include "k-epsilon_turbulent_model_inner.h"

namespace SPH
{
	//=====================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		void K_TurtbulentModelInner::interaction(size_t index_i, Real dt)
		{
				Vecd vel_i = vel_[index_i];
				Real rho_i = rho_[index_i];
				Real turbu_mu_i = turbu_mu_[index_i];
				Real turbu_k_i = turbu_k_[index_i];

				dk_dt_[index_i] = 0.0;;
				Real k_derivative(0.0);
				Real k_lap(0.0);
				Matd strain_rate = Matd::Zero();
				Matd Re_stress = Matd::Zero();
				velocity_gradient[index_i] = Matd::Zero();

				const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
				for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				{
					size_t index_j = inner_neighborhood.j_[n];
					//Vecd nablaW_ijV_j = 0.5* inner_neighborhood.dW_ijV_j_[n] * (B_[index_i] + B_[index_j]) * inner_neighborhood.e_ij_[n];
					//Vecd nablaW_ijV_j =  inner_neighborhood.dW_ijV_j_[n] * B_[index_i]  * inner_neighborhood.e_ij_[n];
					Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] *  inner_neighborhood.e_ij_[n];
					velocity_gradient[index_i] += -(vel_i - vel_[index_j]) * nablaW_ijV_j.transpose();

					k_derivative = (turbu_k_i - turbu_k_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
					k_lap += 2.0 * (mu_+turbu_mu_i/sigma_k)*k_derivative * inner_neighborhood.dW_ijV_j_[n]/ rho_i;
				}

				strain_rate = 0.5 * (velocity_gradient[index_i].transpose() + velocity_gradient[index_i]);
				Re_stress = 2.0 * strain_rate * turbu_mu_i / rho_i - (2.0 / 3.0) * turbu_k_i * Matd::Identity();
				Matd k_production_matrix = Re_stress.array() * velocity_gradient[index_i].array();
				k_production_[index_i] = k_production_matrix.sum();

				dk_dt_[index_i] = k_production_[index_i] - turbu_epsilon_[index_i] + k_lap;

				//** for test */
				k_diffusion_[index_i] = k_lap;
				vel_x_[index_i] = vel_[index_i][0];
		}
		//=================================================================================================//
		void E_TurtbulentModelInner::
			interaction(size_t index_i, Real dt)
		{
			Real rho_i = rho_[index_i];
			Real turbu_mu_i = turbu_mu_[index_i];
			Real turbu_k_i = turbu_k_[index_i];
			Real turbu_epsilon_i = turbu_epsilon_[index_i];

			dE_dt_[index_i] = 0.0;
			Real epsilon_production(0.0);
			Real epsilon_derivative(0.0);
			Real epsilon_lap(0.0);
			Real epsilon_dissipation(0.0);
			const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				epsilon_derivative = (turbu_epsilon_i - turbu_epsilon_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
				epsilon_lap += 2.0 * (mu_ + turbu_mu_i / sigma_E) * epsilon_derivative * inner_neighborhood.dW_ijV_j_[n] / rho_i;
			}

			epsilon_production = C_l * turbu_epsilon_i * k_production_[index_i] / turbu_k_i;
			epsilon_dissipation = C_2 * turbu_epsilon_i * turbu_epsilon_i / turbu_k_i;

			dE_dt_[index_i] = epsilon_production - epsilon_dissipation + epsilon_lap;
			
			//** for test */
			ep_production[index_i] = epsilon_production;
			ep_dissipation_[index_i] = epsilon_dissipation;
			ep_diffusion_[index_i] = epsilon_lap;

		}
		//=================================================================================================//
		void TKEnergyAccInner::interaction(size_t index_i, Real dt)
		{
			Real turbu_k_i = turbu_k_[index_i];
			Vecd acceleration = Vecd::Zero();
			Vecd k_gradient = Vecd::Zero();
			
			const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
				//strong form 
				//k_gradient += -1.0*(turbu_k_i - turbu_k_[index_j]) * nablaW_ijV_j;
				//weak form 
				k_gradient += 1.0 * (turbu_k_i + turbu_k_[index_j]) * nablaW_ijV_j;
			}
			acceleration = -1.0 * (2.0 / 3.0) * k_gradient;
			acc_prior_[index_i] += acceleration;
			
			//for test
			tke_acc_inner_[index_i] = acceleration;
			test_k_grad_rslt_[index_i] = k_gradient;
		}
		//=================================================================================================//
		void TurbuViscousAccInner::
			interaction(size_t index_i, Real dt)
		{
			Real turbu_mu_i = turbu_mu_[index_i];
			Vecd acceleration = Vecd::Zero();
			Vecd vel_derivative = Vecd::Zero();
			const Neighborhood& inner_neighborhood = inner_configuration_[index_i];

			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				vel_derivative = (vel_[index_i] - vel_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
				acceleration += 2.0 * (mu_+ turbu_mu_i) * vel_derivative * inner_neighborhood.dW_ijV_j_[n];
				//acceleration += 2.0 * (mu_ ) * vel_derivative * inner_neighborhood.dW_ijV_j_[n];

			}
			//for test
			visc_acc_inner_[index_i] = acceleration / rho_[index_i];

			acc_prior_[index_i] += acceleration / rho_[index_i];
		}
		//=================================================================================================//




	}
	//=================================================================================================//
}
//=================================================================================================//