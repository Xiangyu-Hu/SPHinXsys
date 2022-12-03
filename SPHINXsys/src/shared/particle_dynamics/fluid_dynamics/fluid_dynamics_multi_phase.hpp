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
 * @file fluid_dynamics_multi_phase.hpp
 * @brief Here, we define the algorithm classes for the dynamics involving multiple fluids.
 * @author	Chi ZHang and Xiangyu Hu
 */

#pragma once

#include "fluid_dynamics_multi_phase.h"

namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		template <class RelaxationInnerType>
		RelaxationMultiPhase<RelaxationInnerType>::
			RelaxationMultiPhase(BaseInnerRelation &inner_relation,
								 BaseContactRelation &contact_relation)
			: RelaxationInnerType(inner_relation), MultiPhaseContactData(contact_relation)
		{
			if (&inner_relation.sph_body_ != &contact_relation.sph_body_)
			{
				std::cout << "\n Error: the two body_relations do not have the same source body!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}

			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_fluids_.push_back(&contact_particles_[k]->fluid_);
				contact_p_.push_back(&(contact_particles_[k]->p_));
				contact_rho_n_.push_back(&(contact_particles_[k]->rho_));
				contact_vel_n_.push_back(&(contact_particles_[k]->vel_));
			}
		}
		//=================================================================================================//
		template <class Integration1stHalfType>
		BaseMultiPhaseIntegration1stHalf<Integration1stHalfType>::
			BaseMultiPhaseIntegration1stHalf(BaseInnerRelation &inner_relation,
											 BaseContactRelation &contact_relation)
			: RelaxationMultiPhase<Integration1stHalfType>(inner_relation, contact_relation)
		{
			for (size_t k = 0; k != this->contact_particles_.size(); ++k)
			{
				riemann_solvers_.push_back(CurrentRiemannSolver(this->fluid_, *this->contact_fluids_[k]));
			}
		}
		//=================================================================================================//
		template <class Integration1stHalfType>
		BaseMultiPhaseIntegration1stHalf<Integration1stHalfType>::
			BaseMultiPhaseIntegration1stHalf(ComplexRelation &complex_relation)
			: BaseMultiPhaseIntegration1stHalf(complex_relation.inner_relation_, complex_relation.contact_relation_) {}
		//=================================================================================================//
		template <class Integration1stHalfType>
		void BaseMultiPhaseIntegration1stHalf<Integration1stHalfType>::interaction(size_t index_i, Real dt)
		{
			Integration1stHalfType::interaction(index_i, dt);

			Vecd acceleration = Vecd::Zero();
			Real rho_dissipation(0);
			for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &p_k = *(this->contact_p_[k]);
				CurrentRiemannSolver &riemann_solver_k = riemann_solvers_[k];
				Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd &e_ij = contact_neighborhood.e_ij_[n];
					Real dW_ijV_j = contact_neighborhood.dW_ijV_j_[n];

					acceleration -= 2.0 * riemann_solver_k.AverageP(this->p_[index_i], p_k[index_j]) * e_ij * dW_ijV_j;
					rho_dissipation += riemann_solver_k.DissipativeUJump(this->p_[index_i] - p_k[index_j]) * dW_ijV_j;
				}
			}
			this->acc_[index_i] += acceleration / this->rho_[index_i];
			this->drho_dt_[index_i] += rho_dissipation * this->rho_[index_i];
		}
		//=================================================================================================//
		template <class Integration1stHalfType>
		Vecd BaseMultiPhaseIntegration1stHalf<Integration1stHalfType>::
			computeNonConservativeAcceleration(size_t index_i)
		{
			Vecd acceleration = this->rho_[index_i] *
								Integration1stHalfType::computeNonConservativeAcceleration(index_i);

			for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &p_k = *(this->contact_p_[k]);
				CurrentRiemannSolver &riemann_solver_k = riemann_solvers_[k];
				Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd &e_ij = contact_neighborhood.e_ij_[n];
					Real dW_ijV_j = contact_neighborhood.dW_ijV_j_[n];

					Real p_ave_k = riemann_solver_k.AverageP(this->p_[index_i], p_k[index_j]);
					acceleration += 2.0 * (this->p_[index_i] - p_ave_k) * dW_ijV_j * e_ij;
				}
			}
			return acceleration / this->rho_[index_i];
		}
		//=================================================================================================//
		template <class Integration2ndHalfType>
		BaseMultiPhaseIntegration2ndHalf<Integration2ndHalfType>::
			BaseMultiPhaseIntegration2ndHalf(BaseInnerRelation &inner_relation,
											 BaseContactRelation &contact_relation)
			: RelaxationMultiPhase<Integration2ndHalfType>(inner_relation, contact_relation)
		{
			for (size_t k = 0; k != this->contact_particles_.size(); ++k)
			{
				riemann_solvers_.push_back(CurrentRiemannSolver(this->fluid_, *this->contact_fluids_[k]));
			}
		}
		//=================================================================================================//
		template <class Integration2ndHalfType>
		BaseMultiPhaseIntegration2ndHalf<Integration2ndHalfType>::
			BaseMultiPhaseIntegration2ndHalf(ComplexRelation &complex_relation)
			: BaseMultiPhaseIntegration2ndHalf(complex_relation.inner_relation_, complex_relation.contact_relation_) {}
		//=================================================================================================//
		template <class Integration2ndHalfType>
		void BaseMultiPhaseIntegration2ndHalf<Integration2ndHalfType>::interaction(size_t index_i, Real dt)
		{
			Integration2ndHalfType::interaction(index_i, dt);

			Real density_change_rate = 0.0;
			Vecd p_dissipation = Vecd::Zero();
			for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &rho_k = *(this->contact_rho_n_[k]);
				StdLargeVec<Real> &p_k = *(this->contact_p_[k]);
				StdLargeVec<Vecd> &vel_k = *(this->contact_vel_n_[k]);
				CurrentRiemannSolver &riemann_solver_k = riemann_solvers_[k];
				Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd &e_ij = contact_neighborhood.e_ij_[n];
					Real dW_ijV_j = contact_neighborhood.dW_ijV_j_[n];

					Vecd vel_ave = riemann_solver_k.AverageV(this->vel_[index_i], vel_k[index_j]);
					density_change_rate += 2.0 * (this->vel_[index_i] - vel_ave).dot(e_ij) * dW_ijV_j;
					Real u_jump = u_jump = (this->vel_[index_i] - vel_k[index_j]).dot(e_ij);
					p_dissipation += riemann_solver_k.DissipativePJump(u_jump) * dW_ijV_j * e_ij;
				}
			}
			this->drho_dt_[index_i] += density_change_rate * this->rho_[index_i];
			this->acc_[index_i] += p_dissipation / this->rho_[index_i];
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//