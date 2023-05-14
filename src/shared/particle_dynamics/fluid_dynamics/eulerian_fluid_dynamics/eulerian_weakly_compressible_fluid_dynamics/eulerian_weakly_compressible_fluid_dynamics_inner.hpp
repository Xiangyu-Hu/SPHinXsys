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
 * @file 	eulerian_weakly_compressible_fluid_dynamics_inner.hpp
 * @brief 	Here, we define the algorithm classes for weakly compressible fluid dynamics within the body.
 * @details We consider here weakly compressible fluids.
 *			TODO: It seems that the eulerian and Lagrangian formulation can be merged together
 * @author	Zhentong Wang, Chi ZHang and Xiangyu Hu
 */

#ifndef EULERIAN_WEAKLY_COMPRESSIBLE_FLUID_DYNAMICS_INNER_HPP
#define EULERIAN_WEAKLY_COMPRESSIBLE_FLUID_DYNAMICS_INNER_HPP

#include "eulerian_weakly_compressible_fluid_dynamics_inner.h"

//=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	namespace eulerian_weakly_compressible_fluid_dynamics
	{
		//=================================================================================================//
		void ViscousAccelerationInner::
			interaction(size_t index_i, Real dt)
		{
			Real rho_i = rho_[index_i];
			const Vecd &vel_i = vel_[index_i];

			Vecd acceleration = Vecd::Zero();
			Vecd vel_derivative = Vecd::Zero();
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				// viscous force
				vel_derivative = (vel_i - vel_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
				acceleration += 2.0 * mu_ * vel_derivative * inner_neighborhood.dW_ijV_j_[n] / rho_i;
			}

			dmom_dt_prior_[index_i] += rho_i * acceleration;
		}
		//=================================================================================================//
		void NonReflectiveBoundaryVariableCorrection::
			interaction(size_t index_i, Real dt)
		{
			Shape &body_shape = *sph_body_.body_shape_;
			if (surface_indicator_[index_i] == 1)
			{
				Vecd normal_direction = body_shape.findNormalDirection(pos_[index_i]);
				n_[index_i] = normal_direction;
				Real velocity_farfield_normal = vel_farfield_.dot(n_[index_i]);
				Real velocity_boundary_normal = vel_[index_i].dot(n_[index_i]);

				// judge it is the inflow condition
				if (n_[index_i][0] <= 0.0 || fabs(n_[index_i][1]) > fabs(n_[index_i][0]))
				{
					// supersonic inflow condition
					if (fabs(velocity_boundary_normal) >= sound_speed_)
					{
						vel_[index_i] = vel_farfield_;
						p_[index_i] = p_farfield_;
						rho_[index_i] = rho_farfield_;
						mom_[index_i] = rho_[index_i] * vel_[index_i];
					}

					// subsonic inflow condition
					if (fabs(velocity_boundary_normal) < sound_speed_)
					{
						Real inner_weight_summation = 0.0;
						Real rho_summation = 0.0;
						Real p_summation = 0.0;
						Real vel_normal_summation(0.0);
						size_t total_inner_neighbor_particles = 0;
						Neighborhood &inner_neighborhood = inner_configuration_[index_i];
						for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
						{
							size_t index_j = inner_neighborhood.j_[n];
							if (surface_indicator_[index_j] != 1)
							{
								Real W_ij = inner_neighborhood.W_ij_[n];
								inner_weight_summation += W_ij * Vol_[index_j];
								rho_summation += rho_[index_j];
								vel_normal_summation += vel_[index_j].dot(n_[index_i]);
								p_summation += p_[index_j];
								total_inner_neighbor_particles += 1;
							}
						}
						Real rho_average = rho_summation / (total_inner_neighbor_particles + TinyReal);
						Real vel_normal_average = vel_normal_summation / (total_inner_neighbor_particles + TinyReal);
						Real p_average = p_summation / (total_inner_neighbor_particles + TinyReal);

						p_[index_i] = p_average * inner_weight_summation + p_farfield_ * (1.0 - inner_weight_summation);
						rho_[index_i] = rho_average * inner_weight_summation + rho_farfield_ * (1.0 - inner_weight_summation);
						Real vel_normal = vel_normal_average * inner_weight_summation + velocity_farfield_normal * (1.0 - inner_weight_summation);
						vel_[index_i] = vel_normal * n_[index_i] + (vel_farfield_ - velocity_farfield_normal * n_[index_i]);
						mom_[index_i] = rho_[index_i] * vel_[index_i];
					}
				}
				// judge it is the outflow condition
				else
				{
					// supersonic outflow condition
					if (fabs(velocity_boundary_normal) >= sound_speed_)
					{
						Real rho_summation = 0.0;
						Real p_summation = 0.0;
						Vecd vel_summation = Vecd::Zero();
						size_t total_inner_neighbor_particles = 0;
						Neighborhood &inner_neighborhood = inner_configuration_[index_i];
						for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
						{
							size_t index_j = inner_neighborhood.j_[n];
							if (surface_indicator_[index_j] != 1)
							{
								rho_summation += rho_[index_j];
								vel_summation += vel_[index_j];
								p_summation += p_[index_j];
								total_inner_neighbor_particles += 1;
							}
						}
						Real rho_average = rho_summation / (total_inner_neighbor_particles + TinyReal);
						Vecd vel_average = vel_summation / (total_inner_neighbor_particles + TinyReal);
						Real p_average = p_summation / (total_inner_neighbor_particles + TinyReal);

						p_[index_i] = p_average;
						rho_[index_i] = rho_average;
						vel_[index_i] = vel_average;
						mom_[index_i] = rho_[index_i] * vel_[index_i];
					}

					// subsonic outflow condition
					if (fabs(velocity_boundary_normal) < sound_speed_)
					{
						Real inner_weight_summation = 0.0;
						Real rho_summation = 0.0;
						Real p_summation = 0.0;
						Real vel_normal_summation(0.0);
						Vecd vel_tangential_summation = Vecd::Zero();
						size_t total_inner_neighbor_particles = 0;
						Neighborhood &inner_neighborhood = inner_configuration_[index_i];
						for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
						{
							size_t index_j = inner_neighborhood.j_[n];
							if (surface_indicator_[index_j] != 1)
							{
								Real W_ij = inner_neighborhood.W_ij_[n];
								inner_weight_summation += W_ij * Vol_[index_j];
								rho_summation += rho_[index_j];
								vel_normal_summation += vel_[index_j].dot(n_[index_i]);
								vel_tangential_summation += vel_[index_j] - (vel_[index_j].dot(n_[index_i])) * n_[index_i];
								p_summation += p_[index_j];
								total_inner_neighbor_particles += 1;
							}
						}
						Real rho_average = rho_summation / (total_inner_neighbor_particles + TinyReal);
						Real vel_normal_average = vel_normal_summation / (total_inner_neighbor_particles + TinyReal);
						Vecd vel_tangential_average = vel_tangential_summation / (total_inner_neighbor_particles + TinyReal);
						Real p_average = p_summation / (total_inner_neighbor_particles + TinyReal);

						p_[index_i] = p_average * inner_weight_summation + p_farfield_ * (1.0 - inner_weight_summation);
						rho_[index_i] = rho_average * inner_weight_summation + rho_farfield_ * (1.0 - inner_weight_summation);
						Real vel_normal = vel_normal_average * inner_weight_summation + velocity_farfield_normal * (1.0 - inner_weight_summation);
						vel_[index_i] = vel_normal * n_[index_i] + vel_tangential_average;
						mom_[index_i] = rho_[index_i] * vel_[index_i];
					}
				}
			}
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		BaseIntegration1stHalf<RiemannSolverType>::BaseIntegration1stHalf(BaseInnerRelation &inner_relation)
			: BaseIntegration(inner_relation), riemann_solver_(this->fluid_, this->fluid_) {}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration1stHalf<RiemannSolverType>::initialization(size_t index_i, Real dt)
		{
			rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
			p_[index_i] = fluid_.getPressure(rho_[index_i]);
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration1stHalf<RiemannSolverType>::update(size_t index_i, Real dt)
		{
			mom_[index_i] += dmom_dt_[index_i] * dt;
			vel_[index_i] = mom_[index_i] / rho_[index_i];
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration1stHalf<RiemannSolverType>::
			interaction(size_t index_i, Real dt)
		{
			FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
			Vecd momentum_change_rate = dmom_dt_prior_[index_i];
			Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
				Vecd &e_ij = inner_neighborhood.e_ij_[n];

				FluidState state_j(rho_[index_j], vel_[index_j], p_[index_j]);
				FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
				Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);

				momentum_change_rate -= 2.0 *
										((rho_star * interface_state.vel_) * interface_state.vel_.transpose() + interface_state.p_ * Matd::Identity()) * e_ij * dW_ijV_j;
			}
			dmom_dt_[index_i] = momentum_change_rate;
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		BaseIntegration2ndHalf<RiemannSolverType>::BaseIntegration2ndHalf(BaseInnerRelation &inner_relation)
			: BaseIntegration(inner_relation), riemann_solver_(fluid_, fluid_) {}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration2ndHalf<RiemannSolverType>::update(size_t index_i, Real dt)
		{
			rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration2ndHalf<RiemannSolverType>::
			interaction(size_t index_i, Real dt)
		{
			FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
			Real density_change_rate = 0.0;
			Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd &e_ij = inner_neighborhood.e_ij_[n];
				Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];

				FluidState state_j(rho_[index_j], vel_[index_j], p_[index_j]);
				FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);

				Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);
				density_change_rate -= 2.0 * (rho_star * interface_state.vel_).dot(e_ij) * dW_ijV_j;
			}
			drho_dt_[index_i] = density_change_rate;
		};
		//=================================================================================================//
	}
	//=================================================================================================//
}
#endif // EULERIAN_WEAKLY_COMPRESSIBLE_FLUID_DYNAMICS_INNER_HPP
