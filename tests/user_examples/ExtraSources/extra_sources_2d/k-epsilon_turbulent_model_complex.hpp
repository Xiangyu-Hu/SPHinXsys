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
 * @file 	k-epsilon_turbulent_model_complex.hpp
 * @brief 	Here.
 * @author	 Xiangyu Hu
 */
#pragma once

#include "k-epsilon_turbulent_model_complex.h"

namespace SPH
{
	//=====================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		void K_TurtbulentModelComplex::
			interaction(size_t index_i, Real dt)
		{
			K_TurtbulentModelInner::interaction(index_i, dt);
			Vecd vel_i = vel_[index_i];
			Real rho_i = rho_[index_i];
			Real turbu_mu_i = turbu_mu_[index_i];
			Real turbu_k_i = turbu_k_[index_i];
			
			//dk_dt_[index_i] = 0.0;
			Real k_production(0.0);
			Real k_derivative(0.0);
			Real k_lap(0.0);
			Matd strain_rate = Matd::Zero();
			Matd Re_stress = Matd::Zero();
			//Matd velocity_gradient = Matd::Zero();

			for (size_t k = 0; k < FluidContactData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd>& vel_ave_k = *(contact_vel_ave_[k]);
				Neighborhood& contact_neighborhood = (*FluidContactData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd nablaW_ijV_j = contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
					velocity_gradient[index_i] += -2.3*(vel_i - vel_ave_k[index_j]) * nablaW_ijV_j.transpose();

					/** With standard wall function, diffusion of k to wall is zero */
					k_derivative = 0.0;
					k_lap += 0.0;
				}
			}
			strain_rate = 0.5 * (velocity_gradient[index_i].transpose() + velocity_gradient[index_i]);
			Re_stress = 2.0 * strain_rate * turbu_mu_i / rho_i - (2.0 / 3.0) * turbu_k_i * Matd::Identity();
			Matd k_production_matrix = Re_stress.array() * velocity_gradient[index_i].array();
			k_production = k_production_matrix.sum();

			/** With standard wall function, epilson on wall is zero */
			dk_dt_[index_i] += k_production - 0.0 + k_lap;

			/** Store the production of K for the use of Epsilon euqation and output*/
			k_production_[index_i] = k_production;
			
			//** for test */
			lap_k_[index_i] += k_lap;
			lap_k_term_[index_i] = 0.0;
			lap_k_term_[index_i] = (mu_ + turbu_mu_i / sigma_k) * lap_k_[index_i];
		}
		//=================================================================================================//
		void E_TurtbulentModelComplex::
			interaction(size_t index_i, Real dt)
		{
			E_TurtbulentModelInner::interaction(index_i, dt);
			Real rho_i = rho_[index_i];
			Real turbu_mu_i = turbu_mu_[index_i];
			Real turbu_k_i = turbu_k_[index_i];
			Real turbu_epsilon_i = turbu_epsilon_[index_i];

			Real epsilon_production(0.0);
			Real epsilon_derivative(0.0);
			Real epsilon_lap(0.0);
			for (size_t k = 0; k < FluidContactData::contact_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = (*FluidContactData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					/** With standard wall function, diffusion of epsilon to wall is zero */
					epsilon_derivative = 0.0;
					epsilon_lap += 0.0;
				}
			}
			epsilon_production = C_l * turbu_epsilon_i * k_production_[index_i] / turbu_k_i;
			/** With standard wall function, epilson on wall is zero */
			dE_dt_[index_i] += epsilon_production - 0.0 + epsilon_lap;
		}
		//=================================================================================================//
		void TKEnergyAccComplex::
			interaction(size_t index_i, Real dt)
		{
			TKEnergyAccInner::interaction(index_i, dt);

			Real turbu_k_i = turbu_k_[index_i];
			Vecd acceleration = Vecd::Zero();
			Vecd k_gradient = Vecd::Zero();
			for (size_t k = 0; k < FluidContactData::contact_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = (*FluidContactData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd nablaW_ijV_j = contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
					//the value of k on wall in this part needs discusstion!
					k_gradient += -1.0 * 2.3 *(turbu_k_i - 0.0) * nablaW_ijV_j;
					//k_gradient +=  (turbu_k_i - 0.0) * nablaW_ijV_j;
					acceleration -= (2.0 / 3.0) * k_gradient;
				}
			}
			//if (surface_indicator_[index_i] == 0 && pos_[index_i][0] <= 5.95)//To prevent kernel truncation near outlet
			
			//According to NHT book, dkdy for wall part should be zero
			//acc_prior_[index_i] += acceleration;
			
			//for test
			tke_acc_wall_[index_i] = acceleration;
			test_k_grad_rslt_[index_i] += k_gradient;
		}
		//=================================================================================================//
		template <class TurbuViscousAccInnerType>
		void BaseTurbuViscousAccWithWall<TurbuViscousAccInnerType>::
			interaction(size_t index_i, Real dt)
		{
			TurbuViscousAccInnerType::interaction(index_i, dt);
			Real turbu_mu_i = this->turbu_mu_[index_i];
			Real rho_i = this->rho_[index_i];
			const Vecd& vel_i = this->vel_[index_i];
			const Vecd& vel_fric_i = this->velo_friction_[index_i];
			Vecd direction_vel_fric = vel_fric_i.normalized();
			
			Real y_plus_i = this->wall_Y_plus_[index_i];

			Vecd acceleration = Vecd::Zero();
			Vecd vel_derivative = Vecd::Zero();
			for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd>& vel_ave_k = *(this->wall_vel_ave_[k]);
				Neighborhood& contact_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Real r_ij = contact_neighborhood.r_ij_[n];

					//vel_derivative = 2.0 * (vel_i - vel_ave_k[index_j]) / (r_ij + 0.01 * this->smoothing_length_);
					//acceleration += 2.0 * (this->mu_+ turbu_mu_i) * vel_derivative * contact_neighborhood.dW_ijV_j_[n] / rho_i;
					//This is to check whether the wall-sub-nearest fluid particles fric, velo. is zero or not
					if (index_i > 2000 && GlobalStaticVariables::physical_time_ > 5. && vel_fric_i.dot(vel_fric_i) == 0.0+TinyReal)
					{
						system("pause");
						std::cout << index_j << std::endl;
					}
					vel_derivative = 2.0 * vel_fric_i.dot(vel_fric_i)* direction_vel_fric;
					acceleration +=  vel_derivative * contact_neighborhood.dW_ijV_j_[n] ;
				}
			}
			//for test

			this->visc_acc_wall_[index_i] = acceleration;

			this->acc_prior_[index_i] += acceleration;
		}
		//=================================================================================================//
		void StandardWallFunctionCorrection::interaction(size_t index_i, Real dt)
		{
			if (is_migrate_[index_i] == 1) //If this particle has started migrating 
			{
				if (is_near_wall_P2_[index_i] == 1) //if it is in P2 region
				{
					is_migrate_[index_i] = 1;  //keep migration status
				}
				else if (is_near_wall_P2_[index_i] == 0) //if it is out of P2
				{
					is_migrate_[index_i] = 0; //ends migration status
				}
			}
			is_near_wall_P2_[index_i] = 0;

			index_nearest[index_i] = 0;
			velo_tan_[index_i] = 0.0;
			velo_friction_[index_i] = Vecd::Zero();
			wall_Y_plus_[index_i] = 0.0;
			wall_Y_star_[index_i] = 0.0;
			distance_to_wall[index_i] = 0.0;
			is_near_wall_P1_[index_i] = 0;

			Real r_wall_normal = 0.0;
			Real r_wall_normal_temp = 0.0;
			Real r_min = 1.0e3;
			Real velo_fric(0.0);
			const Vecd& vel_i = vel_[index_i];
			Real rho_i = rho_[index_i];

			Vecd e_ij_t = Vecd::Zero();
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd>& n_k = *(contact_n_[k]);

				Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Real r_ij = contact_neighborhood.r_ij_[n];
					Vecd& e_ij = contact_neighborhood.e_ij_[n];
					Vecd& n_k_j = n_k[index_j];

					//The distance to wall is 0.5L0 smaller than the rij_normal 
					r_wall_normal_temp = abs(n_k_j.dot(r_ij * e_ij)) - 0.5 * particle_spacing_;
					if (r_wall_normal_temp <= 0.0 + TinyReal)
					{
						std::cout << "r_wall_normal_temp <= 0.0" << std::endl;
						system("pause");
					}
					if (r_wall_normal_temp <= intial_distance_to_wall && r_ij < r_min)
					{
						r_min = r_ij; //Find the nearest wall particle
						r_wall_normal = r_wall_normal_temp;
						distance_to_wall[index_i] = r_wall_normal;
						index_nearest[index_i] = index_j;
					}
					Vecd n_k_j_nearest = n_k[index_nearest[index_i]];
					if (dimension_ == 2)
					{
						e_ij_t[0] = n_k_j_nearest[1];
						e_ij_t[1] = n_k_j_nearest[0] * (-1.0);
						//std::cout << e_ij_t << std::endl;
					}
				}
			}
			if (r_wall_normal < 1.0 * particle_spacing_ &&
				r_wall_normal > 0.0 * particle_spacing_ + TinyReal)
			{
				is_near_wall_P1_[index_i] = 1;
			}
			if (r_wall_normal < (cutoff_radius_ - 0.5 * particle_spacing_) &&
				r_wall_normal > 0.0 * particle_spacing_ + TinyReal)
			{
				Real velo_tan = 0.0; //tangible velo for fluid particle i
				velo_tan = abs(e_ij_t.dot(vel_i));
				velo_tan_[index_i] = velo_tan;
				coefficientA = 0.0;
				coefficientB = 0.0;
				coefficientA = velo_tan * Karman + TinyReal;
				coefficientB = (turbu_const_E * rho_i * r_wall_normal) / mu_;
				velo_fric = getFrictionVelo(0.0, 3.0, 1e-6);
				checkFrictionVelo(velo_fric, 1e-2);

				velo_friction_[index_i] = velo_fric* e_ij_t;
				//friction velocity have the same direction of vel_i, if not, change its direction
				if (vel_i.dot(velo_friction_[index_i]) < 0.0)
					velo_friction_[index_i] = -1.0*velo_friction_[index_i];

			}
			if (r_wall_normal <= 1.5 * particle_spacing_ &&
				r_wall_normal > 1.0 * particle_spacing_)
			{
				is_near_wall_P2_[index_i] = 1;
			}

			if (is_near_wall_P1_[index_i] == 0 && is_near_wall_P1_pre_[index_i] == 1)
			{
				is_migrate_[index_i] = 1;
				//std::cout <<  index_i << "particle starts migrating" << std::endl;
			}

			is_near_wall_P1_pre_[index_i] = 0;
			is_near_wall_P1_pre_[index_i] = is_near_wall_P1_[index_i];

			if (is_near_wall_P1_[index_i] == 1&& pos_[index_i][0]>=0.0) // this is a temporal treamtment.
				//if (r_wall_normal < (cutoff_radius_ - 0.5 * particle_spacing_) &&
				//r_wall_normal > 0.0 * particle_spacing_ + TinyReal && pos_[index_i][0] >= 0.0)
			{
				turbu_k_[index_i] = velo_fric * velo_fric / sqrt(C_mu);
				turbu_epsilon_[index_i] = pow(C_mu, 0.75) * pow(turbu_k_[index_i], 1.5) / (Karman * r_wall_normal);
				wall_Y_plus_[index_i] = r_wall_normal * velo_fric * rho_i / mu_;
				wall_Y_star_[index_i] = r_wall_normal * pow(C_mu, 0.25) * pow(turbu_k_[index_i], 0.5) * rho_i / mu_;
			}
		}
	}
	//=================================================================================================//
}
//=================================================================================================//