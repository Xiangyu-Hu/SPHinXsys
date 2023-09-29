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
		void GetVelocityGradientComplex::
			interaction(size_t index_i, Real dt)
		{
			GetVelocityGradientInner::interaction(index_i, dt);
			Vecd vel_i = vel_[index_i];
			velocity_gradient_wall[index_i] = Matd::Zero();
			for (size_t k = 0; k < FluidContactData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd>& vel_ave_k = *(contact_vel_ave_[k]);
				Neighborhood& contact_neighborhood = (*FluidContactData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd nablaW_ijV_j = contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
					velocity_gradient_[index_i] += -2.0*(vel_i - vel_ave_k[index_j]) * nablaW_ijV_j.transpose();
					//velocity_gradient_[index_i] += -0.0 * (vel_i - vel_ave_k[index_j]) * nablaW_ijV_j.transpose();

					//for test
					velocity_gradient_wall[index_i] += -2.0 * (vel_i - vel_ave_k[index_j]) * nablaW_ijV_j.transpose();
				}
			}
		}
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
			//velocity_gradient_wall[index_i] = Matd::Zero();
			/*
			for (size_t k = 0; k < FluidContactData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd>& vel_ave_k = *(contact_vel_ave_[k]);
				Neighborhood& contact_neighborhood = (*FluidContactData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd nablaW_ijV_j = contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
					//velocity_gradient_[index_i] += -2.3*(vel_i - vel_ave_k[index_j]) * nablaW_ijV_j.transpose();
					velocity_gradient_[index_i] += -0.0 * (vel_i - vel_ave_k[index_j]) * nablaW_ijV_j.transpose();

					//for test 
					velocity_gradient_wall[index_i] += -0.0 * (vel_i - vel_ave_k[index_j]) * nablaW_ijV_j.transpose();


					// With standard wall function, diffusion of k to wall is zero 
					k_derivative = 0.0;
					k_lap += 0.0;
				}
			}
	        */
			strain_rate = 0.5 * (velocity_gradient_[index_i].transpose() + velocity_gradient_[index_i]);
			Re_stress = 2.0 * strain_rate * turbu_mu_i / rho_i - (2.0 / 3.0) * turbu_k_i * Matd::Identity();
			Matd k_production_matrix = Re_stress.array() * velocity_gradient_[index_i].array();
			k_production = k_production_matrix.sum();

			std::cout << "K complex should not be executed, please check" << std::endl;
			system("pause");
			/** With standard wall function, epilson on wall is zero */
			dk_dt_[index_i] += k_production - 0.0 + 0.0;

			/** Store the production of K for the use of Epsilon euqation and output*/
			k_production_[index_i] = k_production;
			
			//** for test */
			k_diffusion_[index_i] += k_lap;
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
			std::cout << "Eplison complex should not be executed, please check" << std::endl;
			system("pause");
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
					//** weak form * 
					k_gradient += -1.0 * (-1.0) * (turbu_k_i + turbu_k_i) * nablaW_ijV_j;
				}
			}
			acceleration = -1.0 * (2.0 / 3.0) * k_gradient;

			acc_prior_[index_i] += acceleration;
			
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
			Vecd e_tau = vel_fric_i.normalized();
			//std::cout << "e_tau=" << e_tau << std::endl;
			Real y_plus_i = this->wall_Y_plus_[index_i];
			Real y_p = this->distance_to_wall_[index_i];

			Matd shear_stress_i_wall = Matd::Zero();

			Real u_plus_i =0.0;
			Real mu_w = 0.0;
			Real mu_p = 0.0;
			Real theta = 0.0;

			Vecd acceleration = Vecd::Zero();
			Vecd vel_derivative = Vecd::Zero();
			for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd>& vel_ave_k = *(this->wall_vel_ave_[k]);
				Neighborhood& contact_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
				StdLargeVec<Vecd>& n_k = *(this->wall_n_[k]);

				//** for test
				if (contact_neighborhood.current_size_ != 0)
				{
					u_plus_i = log(this->turbu_const_E * y_plus_i) / this->Karman;
					mu_w = y_plus_i * this->mu_ / u_plus_i;
					mu_p = this->mu_ + turbu_mu_i;
					theta = y_p* vel_fric_i.norm()/ vel_i.norm();
					if (index_i > 200 && 0)
					{
						std::cout << "******" << std::endl;
						std::cout << "y_plus_i=" << y_plus_i << std::endl << "y_p=" << y_p << std::endl;
						std::cout << "u_plus_i=" << u_plus_i << std::endl;
						std::cout << "u_p=" << vel_i << std::endl;
						std::cout << "mu_w=" << mu_w << std::endl << "mu_p=" << mu_p << std::endl;
						std::cout << "theta=" << theta << std::endl;
						std::cout << "-----" << std::endl;
					}
				}

				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					
					size_t index_j = contact_neighborhood.j_[n];
					Real r_ij = contact_neighborhood.r_ij_[n];
					Vecd& e_ij = contact_neighborhood.e_ij_[n];
					//**Calculate the direction matrix of wall shear stress
					Vecd e_n = n_k[index_j];
					Matd direc_matrix =Matd::Zero();
					direc_matrix = e_tau * e_n.transpose() + (e_tau * e_n.transpose()).transpose();
					//std::cout << "direc_matrix=" << direc_matrix << std::endl;

					//This is to check whether the wall-sub-nearest fluid particles fric, velo. is zero or not
					if (index_i > 2000 && GlobalStaticVariables::physical_time_ > 5. && vel_fric_i.dot(vel_fric_i) <= 0.0+TinyReal&& contact_neighborhood.current_size_>2)
					{
						system("pause");
						std::cout << index_j << std::endl;
						std::cout << vel_fric_i << std::endl;
						std::cout << y_p << std::endl;
						std::cout << contact_neighborhood.current_size_ << std::endl;
					}

					shear_stress_i_wall += rho_i* vel_fric_i.dot(vel_fric_i) * direc_matrix;

					//** Note I think there is something wrong with the direction of dW_ijV_j_[n]
					Vecd acc_j = -1.0 * -1.0 * 2.0 * vel_fric_i.dot(vel_fric_i) * direc_matrix * e_ij * contact_neighborhood.dW_ijV_j_[n];
					acceleration += acc_j;

					//** for test
					if (index_i > 200&&0)
					{
						//std::cout << "mu_eff=" << mu_eff << std::endl ;
						//std::cout << "2*mu_p/rij=" << 2*mu_p / r_ij << std::endl;
						//std::cout << "0.5 * r_ij - theta=" << 0.5 * r_ij - theta << std::endl;
						std::cout << "Acc_j=" << acc_j << std::endl;
						std::cout << "acceleration=" << acceleration << std::endl;
						std::cout << "-----" << std::endl;
					}
				}
			}
			//for test
			Real wall_viscous_factor = 1.0;
			this->visc_acc_wall_[index_i] = wall_viscous_factor * acceleration;
			this->shear_stress_[index_i] += shear_stress_i_wall;
			this->shear_stress_wall_[index_i] = shear_stress_i_wall;

			this->acc_prior_[index_i] += wall_viscous_factor  * acceleration;
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
			distance_to_wall_[index_i] = 0.0;
			is_near_wall_P1_[index_i] = 0;

			Real r_wall_normal = 0.0;
			Real r_wall_normal_temp = 0.0;
			Real r_min = 1.0e3;
			Real velo_fric(0.0);
			const Vecd& vel_i = vel_[index_i];
			Real rho_i = rho_[index_i];

			Real coefficientA = 0.0;
			Real coefficientB = 0.0;

			Vecd e_ij_t = Vecd::Zero();
			Vecd e_n = Vecd::Zero();
			Vecd n_k_j_nearest = Vecd::Zero();
			Matd direc_matrix = Matd::Zero();
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
					//if (r_wall_normal_temp <= intial_distance_to_wall && r_ij < r_min)
					if (r_ij < r_min)
					{
						r_min = r_ij; //Find the nearest wall particle
						r_wall_normal = r_wall_normal_temp;
						distance_to_wall_[index_i] = r_wall_normal;
						index_nearest[index_i] = index_j;
					}
					if (distance_to_wall_[index_i] <= 0.0 + TinyReal)
					{
						std::cout << "strange" << std::endl;
						system("pause");
					}
					n_k_j_nearest = n_k[index_nearest[index_i]];
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
			if (r_wall_normal < (cutoff_radius_ - 0.5 * particle_spacing_)+TinyReal &&
				r_wall_normal > 0.0 * particle_spacing_ + TinyReal)
			{
				is_near_wall_P2_[index_i] = 10;
				Real velo_tan = 0.0; //tangible velo for fluid particle i
				velo_tan = abs(e_ij_t.dot(vel_i));
				velo_tan_[index_i] = velo_tan;
				
				
				//** Wilcox method *
				/*
				coefficientA = velo_tan * Karman + TinyReal;
				coefficientB = (turbu_const_E * rho_i * r_wall_normal) / mu_;				
				velo_fric = getFrictionVelo(0.0, 3.0, 1e-6, coefficientA, coefficientB);
				checkFrictionVelo(velo_fric, 1e-2, coefficientA, coefficientB);
				*/
				//** Fluent method *
				velo_fric = sqrt(abs(Karman * velo_tan * pow(C_mu, 0.25) * pow(turbu_k_[index_i], 0.5) /
					log(turbu_const_E * pow(C_mu, 0.25) * pow(turbu_k_[index_i], 0.5) * r_wall_normal * rho_i / mu_)));

				if (velo_fric != static_cast<Real>(velo_fric)) 
				{
					std::cout << "不是实数" << std::endl;
					std::cout << "velo_fric=" << velo_fric << "velo_tan=" << velo_tan << std::endl;
					std::cout << "turbu_k_=" << pow(turbu_k_[index_i], 0.5) << std::endl;
					std::cout << "sum=" << (Karman * velo_tan * pow(C_mu, 0.25) * pow(turbu_k_[index_i], 0.5) /
						log(turbu_const_E * pow(C_mu, 0.25) * pow(turbu_k_[index_i], 0.5) * r_wall_normal * rho_i / mu_)) << std::endl;
					std::cout << "numerator=" << Karman * velo_tan * pow(C_mu, 0.25) * pow(turbu_k_[index_i], 0.5) << std::endl;
					std::cout << "denominator=" << log(turbu_const_E * pow(C_mu, 0.25) * pow(turbu_k_[index_i], 0.5) * r_wall_normal * rho_i / mu_) << std::endl;

					system("pause");
				}
				velo_friction_[index_i] = velo_fric* e_ij_t;
				//friction velocity have the same direction of vel_i, if not, change its direction
				if (vel_i.dot(velo_friction_[index_i]) < 0.0)
					velo_friction_[index_i] = -1.0*velo_friction_[index_i];
				
				//**Calcualte Y+, including P layer and SUB layer
				wall_Y_plus_[index_i] = r_wall_normal * velo_fric * rho_i / mu_;

			}

			if (is_near_wall_P1_[index_i] == 0 && is_near_wall_P1_pre_[index_i] == 1)
			{
				is_migrate_[index_i] = 1;
				//std::cout <<  index_i << "particle starts migrating" << std::endl;
			}

			is_near_wall_P1_pre_[index_i] = 0;
			is_near_wall_P1_pre_[index_i] = is_near_wall_P1_[index_i];

			direc_matrix = velo_friction_[index_i].normalized() * n_k_j_nearest.transpose();
			//&& pos_[index_i][0]>=0.0
			if (is_near_wall_P1_[index_i] == 1) // this is a temporal treamtment.
			{
				//** Wilcox method *
				//turbu_k_[index_i] = velo_fric * velo_fric / sqrt(C_mu);

				turbu_epsilon_[index_i] = pow(C_mu, 0.75) * pow(turbu_k_[index_i], 1.5) / (Karman * r_wall_normal);
				//wall_Y_plus_[index_i] = r_wall_normal * velo_fric * rho_i / mu_;
				wall_Y_star_[index_i] = r_wall_normal * pow(C_mu, 0.25) * pow(turbu_k_[index_i], 0.5) * rho_i / mu_;
				
				//** Fluent and OpenFoam method *
				Real denominator = pow(C_mu, 0.25) * pow(turbu_k_[index_i], 0.5) * Karman * r_wall_normal;
				//if (GlobalStaticVariables::physical_time_ >= 0.065 && index_i == 499)
				//{
				//	system("pause");
				//}
				velocity_gradient_[index_i] = direc_matrix * velo_fric * velo_fric / denominator;
				k_production_[index_i] = rho_i* pow(velo_fric,4)/ denominator;
			}
		}
	}
	//=================================================================================================//
}
//=================================================================================================//