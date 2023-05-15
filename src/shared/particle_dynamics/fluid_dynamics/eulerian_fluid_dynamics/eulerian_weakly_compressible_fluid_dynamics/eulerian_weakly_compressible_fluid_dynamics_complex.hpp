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
 * @file 	eulerian_weakly_compressible_fluid_dynamics_complex.hpp
 * @brief 	Here, we define the algorithm classes for complex weakly compressible fluid dynamics,
 * 			which is involving with either solid walls (with suffix WithWall)
 * 			or/and other bodies treated as wall for the fluid (with suffix Complex).
 * @author	Zhentong Wang, Chi ZHang and Xiangyu Hu
 */

#ifndef EULERIAN_WEAKLY_COMPRESSIBLE_FLUID_DYNAMICS_COMPLEX_HPP
#define EULERIAN_WEAKLY_COMPRESSIBLE_FLUID_DYNAMICS_COMPLEX_HPP

#include "eulerian_weakly_compressible_fluid_dynamics_complex.h"

//=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	namespace eulerian_weakly_compressible_fluid_dynamics
	{
		//=================================================================================================//
		template <class BaseIntegrationType>
		template <class BaseBodyRelationType>
		InteractionWithWall<BaseIntegrationType>::
			InteractionWithWall(BaseBodyRelationType &base_body_relation, BaseContactRelation &wall_contact_relation)
			: BaseIntegrationType(base_body_relation), WCFluidWallData(wall_contact_relation)
		{
			if (&base_body_relation.getSPHBody() != &wall_contact_relation.getSPHBody())
			{
				std::cout << "\n Error: the two body_relations do not have the same source body!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}

			for (size_t k = 0; k != WCFluidWallData::contact_particles_.size(); ++k)
			{
				Real rho_0_k = WCFluidWallData::contact_bodies_[k]->base_material_->ReferenceDensity();
				wall_inv_rho0_.push_back(1.0 / rho_0_k);
				wall_vel_ave_.push_back(WCFluidWallData::contact_particles_[k]->AverageVelocity());
				wall_acc_ave_.push_back(WCFluidWallData::contact_particles_[k]->AverageAcceleration());
				wall_n_.push_back(&(WCFluidWallData::contact_particles_[k]->n_));
			}
		}
		//=================================================================================================//
		template <class BaseViscousAccelerationType>
		template <class BaseBodyRelationType>
		ViscousWithWall<BaseViscousAccelerationType>::
			ViscousWithWall(BaseBodyRelationType &base_body_relation,
							BaseContactRelation &wall_contact_relation)
			: InteractionWithWall<BaseViscousAccelerationType>(base_body_relation, wall_contact_relation) {}
		//=================================================================================================//
		template <class BaseViscousAccelerationType>
		void ViscousWithWall<BaseViscousAccelerationType>::
			interaction(size_t index_i, Real dt)
		{
			BaseViscousAccelerationType::interaction(index_i, dt);

			Real rho_i = this->rho_[index_i];
			const Vecd &vel_i = this->vel_[index_i];

			Vecd acceleration = Vecd::Zero();
			Vecd vel_derivative = Vecd::Zero();
			for (size_t k = 0; k < WCFluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd> &vel_ave_k = *(this->wall_vel_ave_[k]);
				Neighborhood &contact_neighborhood = (*WCFluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Real r_ij = contact_neighborhood.r_ij_[n];

					vel_derivative = 2.0 * (vel_i - vel_ave_k[index_j]) / (r_ij + 0.01 * this->smoothing_length_);
					acceleration += 2.0 * this->mu_ * vel_derivative * contact_neighborhood.dW_ijV_j_[n] / rho_i;
				}
			}

			this->dmom_dt_prior_[index_i] += acceleration * rho_i;
		}
		//=================================================================================================//
		template <class BaseViscousAccelerationType>
		BaseViscousAccelerationWithWall<BaseViscousAccelerationType>::
			BaseViscousAccelerationWithWall(ComplexRelation &fluid_wall_relation)
			: BaseViscousAccelerationType(fluid_wall_relation.getInnerRelation(),
										  fluid_wall_relation.getContactRelation()) {}
		//=================================================================================================//
		template <class BaseViscousAccelerationType>
		BaseViscousAccelerationWithWall<BaseViscousAccelerationType>::
			BaseViscousAccelerationWithWall(BaseInnerRelation &fluid_inner_relation,
											BaseContactRelation &wall_contact_relation)
			: BaseViscousAccelerationType(fluid_inner_relation, wall_contact_relation) {}
		//=================================================================================================//
		template <class BaseViscousAccelerationType>
		BaseViscousAccelerationWithWall<BaseViscousAccelerationType>::
			BaseViscousAccelerationWithWall(ComplexRelation &fluid_complex_relation,
											BaseContactRelation &wall_contact_relation)
			: BaseViscousAccelerationType(fluid_complex_relation, wall_contact_relation) {}
		//=================================================================================================//
		template <class BaseIntegration1stHalfType>
		template <class BaseBodyRelationType>
		BaseIntegration1stHalfWithWall<BaseIntegration1stHalfType>::
			BaseIntegration1stHalfWithWall(BaseBodyRelationType &base_body_relation,
										   BaseContactRelation &wall_contact_relation)
			: InteractionWithWall<BaseIntegration1stHalfType>(base_body_relation, wall_contact_relation) {}
		//=================================================================================================//
		template <class BaseIntegration1stHalfType>
		void BaseIntegration1stHalfWithWall<BaseIntegration1stHalfType>::
			interaction(size_t index_i, Real dt)
		{
			BaseIntegration1stHalfType::interaction(index_i, dt);

			FluidState state_i(this->rho_[index_i], this->vel_[index_i], this->p_[index_i]);

			Vecd momentum_change_rate = Vecd::Zero();
			for (size_t k = 0; k < WCFluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd> &n_k = *(this->wall_n_[k]);
				Neighborhood &wall_neighborhood = (*WCFluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
				{
					size_t index_j = wall_neighborhood.j_[n];
					Vecd &e_ij = wall_neighborhood.e_ij_[n];
					Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];

					Vecd vel_in_wall = -state_i.vel_;
					Real p_in_wall = state_i.p_;
					Real rho_in_wall = state_i.rho_;
					FluidState state_j(rho_in_wall, vel_in_wall, p_in_wall);
					FluidStarState interface_state = this->riemann_solver_.getInterfaceState(state_i, state_j, n_k[index_j]);
					Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);

					momentum_change_rate -= 2.0 * ((rho_star * interface_state.vel_) * interface_state.vel_.transpose() + interface_state.p_ * Matd::Identity()) * e_ij * dW_ijV_j;
				}
			}
			this->dmom_dt_[index_i] += momentum_change_rate;
		}
		//=================================================================================================//
		template <class BaseIntegration2ndHalfType>
		template <class BaseBodyRelationType>
		BaseIntegration2ndHalfWithWall<BaseIntegration2ndHalfType>::
			BaseIntegration2ndHalfWithWall(BaseBodyRelationType &base_body_relation,
										   BaseContactRelation &wall_contact_relation)
			: InteractionWithWall<BaseIntegration2ndHalfType>(base_body_relation, wall_contact_relation) {}
		//=================================================================================================//
		template <class BaseIntegration2ndHalfType>
		void BaseIntegration2ndHalfWithWall<BaseIntegration2ndHalfType>::
			interaction(size_t index_i, Real dt)
		{
			BaseIntegration2ndHalfType::interaction(index_i, dt);

			FluidState state_i(this->rho_[index_i], this->vel_[index_i], this->p_[index_i]);
			Real density_change_rate = 0.0;
			for (size_t k = 0; k < WCFluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd> &n_k = *(this->wall_n_[k]);
				Neighborhood &wall_neighborhood = (*WCFluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
				{
					size_t index_j = wall_neighborhood.j_[n];
					Vecd &e_ij = wall_neighborhood.e_ij_[n];
					Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];

					Vecd vel_in_wall = -state_i.vel_;
					Real p_in_wall = state_i.p_;
					Real rho_in_wall = state_i.rho_;

					FluidState state_j(rho_in_wall, vel_in_wall, p_in_wall);
					FluidStarState interface_state = this->riemann_solver_.getInterfaceState(state_i, state_j, n_k[index_j]);
					Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);

					density_change_rate -= 2.0 * (rho_star * interface_state.vel_).dot(e_ij) * dW_ijV_j;
				}
			}
			this->drho_dt_[index_i] += density_change_rate;
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
#endif // EULERIAN_WEAKLY_COMPRESSIBLE_FLUID_DYNAMICS_COMPLEX_HPP