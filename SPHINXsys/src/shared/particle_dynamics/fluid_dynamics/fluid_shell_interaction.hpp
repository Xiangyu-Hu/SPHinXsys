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
 *  HU1527/12-1 and HU1527/12-4												*
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
* @file     fluid_shell_interaction.h
* @brief    Here, we define the algorithm classes for the interaction between fluid and
*			thin structure, plate and shell.
* @author	Chi ZHang and Xiangyu Hu
*/
#pragma once

#include "fluid_shell_interaction.h"

//=================================================================================================//
namespace SPH
{
//=================================================================================================//
	namespace fluid_dynamics
	{
    	//=================================================================================================//
		template <class BaseIntegrationType>
		template <class BaseBodyRelationType, typename... Args>
		InteractionWithShell<BaseIntegrationType>::InteractionWithShell(BaseContactRelation &contact_relation, BaseBodyRelationType &base_body_relation, Args &&...args)
			: BaseIntegrationType(base_body_relation, std::forward<Args>(args)...)
			, FluidShellData(contact_relation)
		{
			if (&base_body_relation.sph_body_ != &contact_relation.sph_body_)
			{
				std::cout << "\n Error: the two body_realtions do not have the same source body!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}

			shell_thickness_.reserve(FluidShellData::contact_particles_.size());
			shell_n_.reserve(FluidShellData::contact_particles_.size());
			shell_vel_ave_.reserve(FluidShellData::contact_particles_.size());
			shell_acc_ave_.reserve(FluidShellData::contact_particles_.size());

			for (const auto& cp: FluidShellData::contact_particles_)
			{
				shell_thickness_.push_back(&(cp->thickness_));
				shell_n_.push_back(&(cp->n_));
				shell_vel_ave_.push_back(cp->AverageVelocity());
				shell_acc_ave_.push_back(cp->AverageAcceleration());
			}
		}
		//=================================================================================================//
		template <class Integration1stHalfType>
		void BaseFluidShellIntegration1stHalf<Integration1stHalfType>::interaction(size_t index_i, Real dt)
		{
			Integration1stHalfType::interaction(index_i, dt);
			Vecd acc_prior_i = computeNonConservativeAcceleration(index_i);

			Vecd acceleration = Vecd::Zero();
			Real rho_dissipation = 0.0;
			for (size_t k = 0; k < FluidShellData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd> &vel_ave_k = *(this->shell_vel_ave_[k]);
				StdLargeVec<Vecd> &acc_ave_k = *(this->shell_acc_ave_[k]);
				StdLargeVec<Vecd> &n_k = *(this->shell_n_[k]);
				StdLargeVec<Real> &thickness_k = *(this->shell_thickness_[k]);

				Neighborhood &shell_neighborhood = (*FluidShellData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != shell_neighborhood.current_size_; ++n)
				{
					size_t index_j = shell_neighborhood.j_[n];
					Vecd &e_ij = shell_neighborhood.e_ij_[n];
					Real dW_ijV_j = shell_neighborhood.dW_ijV_j_[n];
					Real r_ij = shell_neighborhood.r_ij_[n];

					Real face_shell_external_acceleration = (acc_prior_i - acc_ave_k[index_j]).dot(-e_ij);
					Real p_in_shell = this->p_[index_i] + this->rho_[index_i] * r_ij * SMAX(0.0, face_shell_external_acceleration);
					acceleration -= (this->p_[index_i] + p_in_shell) * e_ij * dW_ijV_j * thickness_k[index_j];
					rho_dissipation += this->riemann_solver_.DissipativeUJump(this->p_[index_i] - p_in_shell) * dW_ijV_j * thickness_k[index_j];
				}
			}
			this->acc_[index_i] += acceleration / this->rho_[index_i];
			this->drho_dt_[index_i] += rho_dissipation * this->rho_[index_i];
		}
		//=================================================================================================//
		template <class Integration1stHalfType>
		Vecd BaseFluidShellIntegration1stHalf<Integration1stHalfType>::computeNonConservativeAcceleration(size_t index_i)
		{
			return this->acc_prior_[index_i];
		}
		//=================================================================================================//
		template <class BaseIntegration2ndHalfType>
		void BaseFluidShellIntegration2ndHalf<BaseIntegration2ndHalfType>::interaction(size_t index_i, Real dt)
		{
			BaseIntegration2ndHalfType::interaction(index_i, dt);

			Real density_change_rate = 0.0;
			Vecd p_dissipation = Vecd::Zero();

			for (size_t k = 0; k < FluidShellData::contact_configuration_.size(); ++k)
			{
				Vecd &acc_prior_i = this->acc_prior_[index_i];

				StdLargeVec<Vecd> &vel_ave_k = *(this->shell_vel_ave_[k]);
				StdLargeVec<Vecd> &acc_ave_k = *(this->shell_acc_ave_[k]);
				StdLargeVec<Vecd> &n_k = *(this->shell_n_[k]);
				StdLargeVec<Real> &thickness_k = *(this->shell_thickness_[k]);
				Neighborhood &shell_neighborhood = (*FluidShellData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != shell_neighborhood.current_size_; ++n)
				{
					size_t index_j = shell_neighborhood.j_[n];
					Vecd &e_ij = shell_neighborhood.e_ij_[n];
					Real r_ij = shell_neighborhood.r_ij_[n];
					Real dW_ijV_j = shell_neighborhood.dW_ijV_j_[n];
					Vecd correct_n = SGN( e_ij.dot(n_k[index_j]) ) * n_k[index_j];

					Vecd vel_in_shell = 2.0 * vel_ave_k[index_j] - this->vel_[index_i];
					Real u_jump = (this->vel_[index_i] - vel_in_shell).dot(correct_n);
					density_change_rate += u_jump * dW_ijV_j * thickness_k[index_j];
					p_dissipation += this->riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * thickness_k[index_j] * e_ij;
				}
			}
			this->drho_dt_[index_i] += density_change_rate * this->rho_[index_i];
			this->acc_[index_i] += p_dissipation / this->rho_[index_i];
		}
		//=================================================================================================//
		template <class ViscousAccelerationInnerType>
		void BaseViscousAccelerationWithShell<ViscousAccelerationInnerType>::interaction(size_t index_i, Real dt)
		{
			ViscousAccelerationInnerType::interaction(index_i, dt);

			Real rho_i = this->rho_[index_i];
			const Vecd &vel_i = this->vel_[index_i];

			Vecd acceleration = Vecd::Zero();
			Vecd vel_derivative = Vecd::Zero();
			for (size_t k = 0; k < FluidShellData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd> &vel_ave_k = *(this->shell_vel_ave_[k]);
				StdLargeVec<Real> &thickness_k = *(this->shell_thickness_[k]);
				Neighborhood &contact_neighborhood = (*FluidShellData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Real r_ij = contact_neighborhood.r_ij_[n];

					vel_derivative = 2.0 * (vel_i - vel_ave_k[index_j]) / (r_ij + 0.01 * this->smoothing_length_);
					acceleration += 2.0 * this->mu_ * vel_derivative * contact_neighborhood.dW_ijV_j_[n] * thickness_k[index_j] / rho_i;
				}
			}

			this->acc_prior_[index_i] += acceleration;
		}
    	//=================================================================================================//
    }
//=================================================================================================//
}