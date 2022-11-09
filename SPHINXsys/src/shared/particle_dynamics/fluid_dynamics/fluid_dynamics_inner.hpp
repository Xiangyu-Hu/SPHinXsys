/**
 * @file 	fluid_dynamics_inner.hpp
 * @author	Chi Zhang and Xiangyu Hu
 */

#pragma once

#include "fluid_dynamics_inner.h"

//=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		template <class RiemannSolverType>
		BaseIntegration1stHalf<RiemannSolverType>::BaseIntegration1stHalf(BaseInnerRelation &inner_relation)
			: BaseIntegration(inner_relation), riemann_solver_(fluid_, fluid_) {}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration1stHalf<RiemannSolverType>::initialization(size_t index_i, Real dt)
		{
			rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
			p_[index_i] = fluid_.getPressure(rho_[index_i]);
			pos_[index_i] += vel_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration1stHalf<RiemannSolverType>::update(size_t index_i, Real dt)
		{
			vel_[index_i] += (acc_prior_[index_i] + acc_[index_i]) * dt;
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		Vecd BaseIntegration1stHalf<RiemannSolverType>::computeNonConservativeAcceleration(size_t index_i)
		{
			Vecd acceleration = acc_prior_[index_i] * rho_[index_i];
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
				const Vecd &e_ij = inner_neighborhood.e_ij_[n];

				acceleration += (p_[index_i] - p_[index_j]) * dW_ijV_j * e_ij;
			}
			return acceleration / rho_[index_i];
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration1stHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
		{
			Vecd acceleration(0);
			Real rho_dissipation(0);
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
				const Vecd &e_ij = inner_neighborhood.e_ij_[n];

				acceleration -= (p_[index_i] + p_[index_j]) * dW_ijV_j * e_ij;
				rho_dissipation += riemann_solver_.DissipativeUJump(p_[index_i] - p_[index_j]) * dW_ijV_j;
			}
			acc_[index_i] += acceleration / rho_[index_i];
			drho_dt_[index_i] = rho_dissipation * rho_[index_i];
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		BaseIntegration2ndHalf<RiemannSolverType>::BaseIntegration2ndHalf(BaseInnerRelation &inner_relation)
			: BaseIntegration(inner_relation), riemann_solver_(fluid_, fluid_),
			  Vol_(particles_->Vol_), mass_(particles_->mass_) {}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration2ndHalf<RiemannSolverType>::initialization(size_t index_i, Real dt)
		{
			pos_[index_i] += vel_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration2ndHalf<RiemannSolverType>::update(size_t index_i, Real dt)
		{
			rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
			Vol_[index_i] = mass_[index_i] / rho_[index_i];
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseIntegration2ndHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
		{
			Real density_change_rate(0);
			Vecd p_dissipation(0);
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				const Vecd &e_ij = inner_neighborhood.e_ij_[n];
				Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];

				Real u_jump = dot(vel_[index_i] - vel_[index_j], e_ij);
				density_change_rate += u_jump * dW_ijV_j;
				p_dissipation += riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * e_ij;
			}
			drho_dt_[index_i] += density_change_rate * rho_[index_i];
			acc_[index_i] = p_dissipation / rho_[index_i];
		};
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//