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
 * @brief 	
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

					//for test
					velocity_gradient_wall[index_i] += -2.0 * (vel_i - vel_ave_k[index_j]) * nablaW_ijV_j.transpose();
				}
			}
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
			Real y_plus_i = this->wall_Y_plus_[index_i];
			Real y_p = this->distance_to_wall_[index_i];

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

				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Real r_ij = contact_neighborhood.r_ij_[n];
					Vecd& e_ij = contact_neighborhood.e_ij_[n];
					//** Calculate the direction matrix of wall shear stress *
					Vecd e_n = n_k[index_j];
					Matd direc_matrix =Matd::Zero();
					direc_matrix = e_tau * e_n.transpose() + (e_tau * e_n.transpose()).transpose();

					//** This is to check whether the wall-sub-nearest fluid particles fric, velo. is zero *
					if (index_i > 2000 && GlobalStaticVariables::physical_time_ > 5. && vel_fric_i.dot(vel_fric_i) <= 0.0+TinyReal&& contact_neighborhood.current_size_>2)
					{
						system("pause");
						std::cout << index_j << std::endl;
						std::cout << vel_fric_i << std::endl;
						std::cout << y_p << std::endl;
						std::cout << contact_neighborhood.current_size_ << std::endl;
					}

					Vecd acc_j = -1.0 * -1.0 * 2.0 * vel_fric_i.dot(vel_fric_i) * direc_matrix * e_ij * contact_neighborhood.dW_ijV_j_[n];
					acceleration += acc_j;
				}
			}
			this->acc_prior_[index_i] +=  acceleration;
			//** For test *
			this->visc_acc_wall_[index_i] = acceleration;
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
			y_p_[index_i] = 0.0;
			Real r_dummy_normal = 0.0;
			Real r_dummy_normal_temp = 0.0;
			Real r_min = 1.0e3;
			Real velo_fric(0.0);
			const Vecd& vel_i = vel_[index_i];
			Real rho_i = rho_[index_i];
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

					//** The distance to dummy interface is 0.5 dp smaller than the r_ij_normal *  
					r_dummy_normal_temp = abs(n_k_j.dot(r_ij * e_ij)) - 0.5 * particle_spacing_;
					if (r_dummy_normal_temp <= 0.0 + TinyReal)
					{
						std::cout << "r_dummy_normal_temp <= 0.0" << std::endl;
						system("pause");
					}
					if (r_ij < r_min)
					{
						r_min = r_ij; //** Find the nearest wall particle *
						r_dummy_normal = r_dummy_normal_temp;
						distance_to_wall_[index_i] = r_dummy_normal;
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
					}
				}
			}
			if (r_dummy_normal < 1.0 * particle_spacing_ &&
				r_dummy_normal > 0.0 * particle_spacing_ + TinyReal)
			{
				is_near_wall_P1_[index_i] = 1;
			}
			if (r_dummy_normal < (cutoff_radius_ - 0.5 * particle_spacing_)+TinyReal &&
				r_dummy_normal > 0.0 * particle_spacing_ + TinyReal)
			{
				is_near_wall_P2_[index_i] = 10;
				Real velo_tan = 0.0; //** tangitial velo for fluid particle i *
				velo_tan = abs(e_ij_t.dot(vel_i));
				velo_tan_[index_i] = velo_tan;

				y_p_[index_i] = r_dummy_normal + offset_dist_;

				velo_fric = sqrt(abs(Karman * velo_tan * pow(C_mu, 0.25) * pow(turbu_k_[index_i], 0.5) /
					log(turbu_const_E * pow(C_mu, 0.25) * pow(turbu_k_[index_i], 0.5) * y_p_[index_i] * rho_i / mu_)));

				if (velo_fric != static_cast<Real>(velo_fric)) 
				{
					std::cout << "friction velocity is not a real, please check" << std::endl;
					std::cout << "velo_fric=" << velo_fric << "velo_tan=" << velo_tan << std::endl;
					std::cout << "turbu_k_=" << pow(turbu_k_[index_i], 0.5) << std::endl;
					std::cout << "sum=" << (Karman * velo_tan * pow(C_mu, 0.25) * pow(turbu_k_[index_i], 0.5) /
						log(turbu_const_E * pow(C_mu, 0.25) * pow(turbu_k_[index_i], 0.5) * r_dummy_normal * rho_i / mu_)) << std::endl;
					std::cout << "numerator=" << Karman * velo_tan * pow(C_mu, 0.25) * pow(turbu_k_[index_i], 0.5) << std::endl;
					std::cout << "denominator=" << log(turbu_const_E * pow(C_mu, 0.25) * pow(turbu_k_[index_i], 0.5) * r_dummy_normal * rho_i / mu_) << std::endl;

					system("pause");
				}
				velo_friction_[index_i] = velo_fric* e_ij_t;
				//** friction velocity have the same direction of vel_i, if not, change its direction *
				if (vel_i.dot(velo_friction_[index_i]) < 0.0)
					velo_friction_[index_i] = -1.0*velo_friction_[index_i];
				
				//** Calcualte Y+, including P layer and SUB layer *
				wall_Y_plus_[index_i] = y_p_[index_i] * velo_fric * rho_i / mu_;
			}
			if (is_near_wall_P1_[index_i] == 0 && is_near_wall_P1_pre_[index_i] == 1)
			{
				is_migrate_[index_i] = 1;
			}
			is_near_wall_P1_pre_[index_i] = 0;
			is_near_wall_P1_pre_[index_i] = is_near_wall_P1_[index_i];

			direc_matrix = velo_friction_[index_i].normalized() * n_k_j_nearest.transpose();
			
			if (is_near_wall_P1_[index_i] == 1) // ** Correct the near wall values *
			{
				turbu_epsilon_[index_i] = pow(C_mu, 0.75) * pow(turbu_k_[index_i], 1.5) / (Karman * y_p_[index_i]);
				wall_Y_star_[index_i] = y_p_[index_i] * pow(C_mu, 0.25) * pow(turbu_k_[index_i], 0.5) * rho_i / mu_;
				Real denominator = pow(C_mu, 0.25) * pow(turbu_k_[index_i], 0.5) * Karman * y_p_[index_i];
				velocity_gradient_[index_i] = direc_matrix * velo_fric * velo_fric / denominator;
				k_production_[index_i] = rho_i* pow(velo_fric,4)/ denominator;
			}
		}
	}
	//=================================================================================================//
}
//=================================================================================================//