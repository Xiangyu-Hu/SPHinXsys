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

				Real mu_eff_i = turbu_mu_[index_i]/sigma_k + mu_;

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
					//** Strong form *
					velocity_gradient[index_i] += -(vel_i - vel_[index_j]) * nablaW_ijV_j.transpose();
					//** Weak form *
					//velocity_gradient[index_i] += (vel_i + vel_[index_j]) * nablaW_ijV_j.transpose();

					Real mu_eff_j = turbu_mu_[index_j] / sigma_k + mu_;
					Real mu_harmo = 2 * mu_eff_i * mu_eff_j / (mu_eff_i + mu_eff_j);
					k_derivative = (turbu_k_i - turbu_k_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
					k_lap += 2.0 * mu_harmo * k_derivative * inner_neighborhood.dW_ijV_j_[n]/ rho_i;
				}

				strain_rate = 0.5 * (velocity_gradient[index_i].transpose() + velocity_gradient[index_i]);
				Re_stress = 2.0 * strain_rate * turbu_mu_i / rho_i - (2.0 / 3.0) * turbu_k_i * Matd::Identity();
				Matd k_production_matrix = Re_stress.array() * velocity_gradient[index_i].array();
				k_production_[index_i] = k_production_matrix.sum();

				dk_dt_[index_i] = k_production_[index_i] - turbu_epsilon_[index_i] + k_lap;

				//**Temp treatment to see the influence of S region * 
				//Real correct_factor_PS = 0.0011;
				//if (vel_[index_i][0] < 0.9)// ** Particle P and S *
				//{
				//	dk_dt_[index_i] = dk_dt_[index_i] + correct_factor_PS;
				//}


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

			Real mu_eff_i = turbu_mu_[index_i] / sigma_E + mu_;

			dE_dt_[index_i] = 0.0;
			Real epsilon_production(0.0);
			Real epsilon_derivative(0.0);
			Real epsilon_lap(0.0);
			Real epsilon_dissipation(0.0);
			const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real mu_eff_j = turbu_mu_[index_j] / sigma_E + mu_;
				Real mu_harmo = 2 * mu_eff_i * mu_eff_j / (mu_eff_i + mu_eff_j);
				epsilon_derivative = (turbu_epsilon_i - turbu_epsilon_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
				epsilon_lap += 2.0 * mu_harmo * epsilon_derivative * inner_neighborhood.dW_ijV_j_[n] / rho_i;
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
			
			//std::cout << "hahaha" << std::endl;

			const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
				//** strong form * 
				//k_gradient += -1.0*(turbu_k_i - turbu_k_[index_j]) * nablaW_ijV_j;
				//** weak form * 
				k_gradient += -1.0 * (- 1.0) * (turbu_k_i + turbu_k_[index_j]) * nablaW_ijV_j;
			}
			
			acceleration = -1.0 * (2.0 / 3.0) * k_gradient;

			//**Temp treatment to see the influence of TKE * 
			//Real correct_factor_P = 0.06;
			//Real correct_factor_S = 0.01;
			//if (is_near_wall_P1_[index_i]==1)// ** Particle P *
			//{
			//	if (pos_[index_i][1]<1.0)
			//	{
			//		acceleration[1] = acceleration[1] + correct_factor_P;
			//	}
			//	else if (pos_[index_i][1] >= 1.0)
			//	{
			//		acceleration[1] = acceleration[1] - correct_factor_P;
			//	}
			//}
			//else if (is_near_wall_P1_[index_i] == 0 && is_near_wall_P2_[index_i] == 10) // ** Particle S *
			//{
			//	if (pos_[index_i][1] < 1.0)
			//	{
			//		acceleration[1] = acceleration[1] + correct_factor_S;
			//	}
			//	else if (pos_[index_i][1] >= 1.0)
			//	{
			//		acceleration[1] = acceleration[1] - correct_factor_S;
			//	}
			//}
			
			acc_prior_[index_i] += acceleration;
			
			//for test
			tke_acc_inner_[index_i] = acceleration;
			test_k_grad_rslt_[index_i] = k_gradient;
		}
		//=================================================================================================//
		void TurbuViscousAccInner::
			interaction(size_t index_i, Real dt)
		{
			Real mu_eff_i = turbu_mu_[index_i] + mu_;

			Matd shear_stress_i = Matd::Zero();
			Matd shear_stress_j = Matd::Zero();

			Vecd acceleration = Vecd::Zero();
			Vecd vel_derivative = Vecd::Zero();
			const Neighborhood& inner_neighborhood = inner_configuration_[index_i];

			//** for test *
			Real acc_j_wall_side = 0.0;
			Real acc_j_inner_side = 0.0;
			Real visc_force_j_wall_side = 0.0;
			Real visc_force_j_inner_side = 0.0;

			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				const Vecd& e_ij = inner_neighborhood.e_ij_[n];

				Real mu_eff_j = turbu_mu_[index_j] + mu_;
				Real mu_harmo = 2 * mu_eff_i * mu_eff_j / (mu_eff_i + mu_eff_j);
				vel_derivative = (vel_[index_i] - vel_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
				
				shear_stress_j = mu_harmo * (vel_derivative * e_ij.transpose());
				shear_stress_i += shear_stress_j;
				Vecd acc_j = 2.0 * mu_harmo * vel_derivative * inner_neighborhood.dW_ijV_j_[n];
				acceleration += acc_j;
				//if (index_i == 490)
				//{
				//	std::cout << "------" << std::endl;
				//	std::cout << "index_j=" << n << std::endl;
				//	std::cout << "dW_ijV_j_=" << inner_neighborhood.dW_ijV_j_[n] << std::endl;
				//	std::cout << "s_s_j=" << shear_stress_j << std::endl;
				//	std::cout << "shear_stress_i=" << shear_stress_i << std::endl;
				//	std::cout << "e_ij=" << e_ij << std::endl;
				//	std::cout << "r_ij/l0=" << inner_neighborhood.r_ij_[n] / 0.1 << std::endl;
				//	std::cout << "Acc_j=" << acc_j << std::endl;
				//	std::cout << "------" << std::endl;
				//}

				//** for test
				//if (index_i > 200 && (turbu_mu_[index_i] - 0.001025458) < TinyReal)
				//{
				//	std::cout << "---INNER---ParticleNo.=" << index_i<< std::endl;
				//	std::cout << "Acc_j=" << acc_j << std::endl;
				//	std::cout << "acceleration=" << acceleration << std::endl;
				//	if (vel_[index_i][0] < vel_[index_j][0])
				//	{
				//		std::cout << "J belongs to inner side" << std::endl;
				//		acc_j_inner_side = acc_j_inner_side+ acc_j[0];
				//		visc_force_j_inner_side = visc_force_j_inner_side + mu_harmo * vel_derivative[0];
				//		std::cout << "acc_j_inner_side=" << acc_j_inner_side<< std::endl;
				//		std::cout << "visc_force_j_inner_side=" << visc_force_j_inner_side << std::endl;
				//	}
				//	else if(vel_[index_i][0] > vel_[index_j][0])
				//	{
				//		std::cout << "J belongs to wall side" << std::endl;
				//		acc_j_wall_side = acc_j_wall_side + acc_j[0];
				//		visc_force_j_wall_side = visc_force_j_wall_side + mu_harmo * vel_derivative[0];
				//		std::cout << "acc_j_wall_side=" << acc_j_wall_side << std::endl;
				//		std::cout << "visc_force_j_wall_side=" << visc_force_j_wall_side << std::endl;
				//	}
				//	std::cout << "-----" << std::endl;
				//}

				//acceleration += 2.0 * mu_harmo * vel_derivative * inner_neighborhood.dW_ijV_j_[n];
				//acceleration += 2.0 * (mu_ ) * vel_derivative * inner_neighborhood.dW_ijV_j_[n];


			}
			//for test
			visc_acc_inner_[index_i] = acceleration / rho_[index_i];
			shear_stress_[index_i] = shear_stress_i;

			acc_prior_[index_i] += acceleration / rho_[index_i];
		}
		//=================================================================================================//




	}
	//=================================================================================================//
}
//=================================================================================================//