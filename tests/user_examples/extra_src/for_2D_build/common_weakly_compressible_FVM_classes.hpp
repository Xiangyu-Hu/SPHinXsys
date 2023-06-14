/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*                                                                         *
 * ------------------------------------------------------------------------*/
#pragma once

#include "common_weakly_compressible_FVM_classes.h"


namespace SPH
{
	//=================================================================================================//
	template<class RiemannSolverType>
	BaseViscousAccelerationInnerInFVM<RiemannSolverType>::
    BaseViscousAccelerationInnerInFVM(BaseInnerRelationInFVM &inner_relation, Real limiter_parameter)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInnerInFVM<BaseParticles>(inner_relation),
      fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->base_material_)), riemann_solver_(fluid_, fluid_, limiter_parameter), rho_(particles_->rho_),
      p_(*particles_->getVariableByName<Real>("Pressure")), vel_(particles_->vel_), acc_prior_(particles_->acc_prior_), pos_(particles_->pos_),
      dmom_dt_prior_(*particles_->getVariableByName<Vecd>("OtherMomentumChangeRate")),
      mu_(fluid_.ReferenceViscosity()){};
	//=================================================================================================//
	template<class RiemannSolverType>
	void BaseViscousAccelerationInnerInFVM<RiemannSolverType>::interaction(size_t index_i, Real dt)
	{
		Real rho_i = rho_[index_i];
		const Vecd& vel_i = vel_[index_i];
		FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
		Vecd acceleration = Vecd::Zero();
		Vecd vel_derivative = Vecd::Zero();
		const NeighborhoodInFVM& inner_neighborhood = inner_configuration_in_FVM_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			//viscous force in inner particles
			if (inner_neighborhood.boundary_type_[n] == 2)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd vel_j = vel_[index_j];
				vel_derivative = (vel_i - vel_j) / (inner_neighborhood.r_ij_[n] + TinyReal);
				acceleration += 2.0 * mu_ * vel_derivative * inner_neighborhood.dW_ijV_j_[n] / rho_i;
			}
			//viscous force between fluid partile and wall boundary
			if (inner_neighborhood.boundary_type_[n] == 3)
			{
				Vecd vel_j = -vel_i;
				vel_derivative = (vel_i - vel_j) / (inner_neighborhood.r_ij_[n] + TinyReal);
				acceleration += 2.0 * mu_ * vel_derivative * inner_neighborhood.dW_ijV_j_[n] / rho_i;
			}
			//viscous force between fluid partile and far-field boundary condition
			if (inner_neighborhood.boundary_type_[n] == 9)
			{
				Vecd e_ij = inner_neighborhood.e_ij_[n];
				Vecd far_field_velocity(1.0, 0.0);
				Real far_field_density = 1.0;
				Real far_field_pressure = fluid_.getPressure(far_field_density);
				FluidState state_farfield(far_field_density, far_field_velocity, far_field_pressure);
				FluidStarState state_farfield_star = riemann_solver_.getInterfaceState(state_i, state_farfield, e_ij);
				Vecd vel_star = state_farfield_star.vel_;
				Real vel_normal = vel_star.dot(-e_ij);
				Vecd vel_boundary_condition = vel_normal * (-e_ij) + far_field_velocity - far_field_velocity.dot(-e_ij) * (-e_ij);

				vel_derivative = (vel_i - vel_boundary_condition) / (inner_neighborhood.r_ij_[n] + TinyReal);
				acceleration += 2.0 * mu_ * vel_derivative * inner_neighborhood.dW_ijV_j_[n] / rho_i;
			}
            //viscous force between fluid partile and velocity inlet condition, Here we set parabolic fucntion velocity inlet
            if (inner_neighborhood.boundary_type_[n] == 10)
            {
                Real U_f = 1.0; //characteristic velocity is set as 1
                Real h = 1.0; // the height of the inflow domain size
                Vecd parabolic_velocity_inlet = Vecd::Zero();
                parabolic_velocity_inlet[0] = 1.5 * U_f * (1.0 - pos_[index_i][1] * pos_[index_i][1] / pow(0.5 * h, 2));
                Vecd vel_j = parabolic_velocity_inlet;
                vel_derivative = (vel_i - vel_j) / (inner_neighborhood.r_ij_[n] + TinyReal);
                acceleration += 2.0 * mu_ * vel_derivative * inner_neighborhood.dW_ijV_j_[n] / rho_i;
            }
		}
		dmom_dt_prior_[index_i] += rho_i * acceleration;
	}
	//=================================================================================================//
	template <class RiemannSolverType>
	BaseIntegration1stHalfInFVM<RiemannSolverType>::BaseIntegration1stHalfInFVM(BaseInnerRelationInFVM &inner_relation, Real limiter_parameter)
    : BaseRelaxationInFVM(inner_relation), riemann_solver_(fluid_, fluid_, limiter_parameter) {}
	//=================================================================================================//
	template <class RiemannSolverType>
	void BaseIntegration1stHalfInFVM<RiemannSolverType>::initialization(size_t index_i, Real dt)
	{
		rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
		p_[index_i] = fluid_.getPressure(rho_[index_i]);
	}
	//=================================================================================================//
	template <class RiemannSolverType>
	void BaseIntegration1stHalfInFVM<RiemannSolverType>::interaction(size_t index_i, Real dt)
	{
		FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
		Vecd momentum_change_rate = dmom_dt_prior_[index_i];
		NeighborhoodInFVM& inner_neighborhood = inner_configuration_in_FVM_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			Vecd& e_ij = inner_neighborhood.e_ij_[n];
			Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
			if (inner_neighborhood.boundary_type_[n] == 2)
			{
				FluidState state_j(rho_[index_j], vel_[index_j], p_[index_j]);
				FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
				Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);

				momentum_change_rate -= 2.0 * ((rho_star * interface_state.vel_) * interface_state.vel_.transpose() + interface_state.p_ * Matd::Identity()) * e_ij * dW_ijV_j;
			}

			if (inner_neighborhood.boundary_type_[n] == 3)
			{
				Vecd vel_in_wall = -state_i.vel_;
				Real p_in_wall = state_i.p_;
				Real rho_in_wall = state_i.rho_;
				FluidState state_j(rho_in_wall, vel_in_wall, p_in_wall);
				FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
				Real rho_star = fluid_.DensityFromPressure(interface_state.p_);

				momentum_change_rate -= 2.0 * ((rho_star * interface_state.vel_) * interface_state.vel_.transpose() + interface_state.p_ * Matd::Identity()) * e_ij * dW_ijV_j;
			}

			if (inner_neighborhood.boundary_type_[n] == 9)
			{
				Vecd far_field_velocity(1.0, 0.0);
				Real far_field_density = 1.0;
				Real far_field_pressure = fluid_.getPressure(far_field_density);
				FluidState state_boundary(far_field_density, far_field_velocity, far_field_pressure);
				FluidStarState boundary_star_state = riemann_solver_.getInterfaceState(state_i, state_boundary, e_ij);
				Vecd vel_star = boundary_star_state.vel_;
				Real p_j = boundary_star_state.p_;
				Real vel_normal = vel_star.dot(-e_ij);
				Vecd vel_j = vel_normal * (-e_ij) + far_field_velocity - far_field_velocity.dot(-e_ij) * (-e_ij);
				Real rho_j = fluid_.DensityFromPressure(p_j);
				FluidState state_j(rho_j, vel_j, p_j);
				FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
				Real rho_star = fluid_.DensityFromPressure(interface_state.p_);

				momentum_change_rate -= 2.0 * ((rho_star * interface_state.vel_) * interface_state.vel_.transpose() + interface_state.p_ * Matd::Identity()) * e_ij * dW_ijV_j;
			}
            if (inner_neighborhood.boundary_type_[n] == 10)
            {
                Real U_f = 1.0; //characteristic velocity is set as 1
                Real h = 1.0;   // the height of the inflow domain size
                Vecd parabolic_velocity_inlet = Vecd::Zero();
                parabolic_velocity_inlet[0] = 1.5 * U_f * (1.0 - pos_[index_i][1] * pos_[index_i][1] / pow(0.5 * h, 2));
                Vecd vel_inlet = parabolic_velocity_inlet;
                Real p_inlet = state_i.p_;
                Real rho_inlet = state_i.rho_;
                FluidState state_j(rho_inlet, vel_inlet, p_inlet);
                FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
                Real rho_star = fluid_.DensityFromPressure(interface_state.p_);

                momentum_change_rate -= 2.0 * ((rho_star * interface_state.vel_) * interface_state.vel_.transpose() + interface_state.p_ * Matd::Identity()) * e_ij * dW_ijV_j;
            }
		}
		dmom_dt_[index_i] = momentum_change_rate;
	}
	//=================================================================================================//
	template <class RiemannSolverType>
	void BaseIntegration1stHalfInFVM<RiemannSolverType>::update(size_t index_i, Real dt)
	{
		mom_[index_i] += dmom_dt_[index_i] * dt;
		vel_[index_i] = mom_[index_i] / rho_[index_i];
	}
	//=================================================================================================//
	template <class RiemannSolverType>
	BaseIntegration2ndHalfInFVM<RiemannSolverType>::BaseIntegration2ndHalfInFVM(BaseInnerRelationInFVM &inner_relation, Real limiter_parameter)
    : BaseRelaxationInFVM(inner_relation), riemann_solver_(fluid_, fluid_, limiter_parameter) {}
	//=================================================================================================//
	template <class RiemannSolverType>
	void BaseIntegration2ndHalfInFVM<RiemannSolverType>::interaction(size_t index_i, Real dt)
	{
		FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
		Real density_change_rate = 0.0;
		NeighborhoodInFVM& inner_neighborhood = inner_configuration_in_FVM_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			Vecd& e_ij = inner_neighborhood.e_ij_[n];
			Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
			if (inner_neighborhood.boundary_type_[n] == 2)
			{
				FluidState state_j(rho_[index_j], vel_[index_j], p_[index_j]);
				FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
				Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);

				density_change_rate -= 2.0 * (rho_star * interface_state.vel_).dot(e_ij) * dW_ijV_j;
			}

			if (inner_neighborhood.boundary_type_[n] == 3)
			{
				Vecd vel_in_wall = -state_i.vel_;
				Real p_in_wall = state_i.p_;
				Real rho_in_wall = state_i.rho_;
				FluidState state_j(rho_in_wall, vel_in_wall, p_in_wall);
				FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
				Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);

				density_change_rate -= 2.0 * (rho_star * interface_state.vel_).dot(e_ij) * dW_ijV_j;
			}

			if (inner_neighborhood.boundary_type_[n] == 9)
			{
				Vecd far_field_velocity(1.0, 0.0);
				Real far_field_density = 1.0;
				Real far_field_pressure = fluid_.getPressure(far_field_density);
				FluidState state_boundary(far_field_density, far_field_velocity, far_field_pressure);
				FluidStarState boundary_star_state = riemann_solver_.getInterfaceState(state_i, state_boundary, e_ij);
				Vecd vel_star = boundary_star_state.vel_;
				Real vel_normal = vel_star.dot(-e_ij);
				Real p_j = boundary_star_state.p_;
				Vecd vel_j = vel_normal * (-e_ij) + far_field_velocity - far_field_velocity.dot(-e_ij) * (-e_ij);
				Real rho_j = fluid_.DensityFromPressure(p_j);
				FluidState state_j(rho_j, vel_j, p_j);
				FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
				Real rho_star = fluid_.DensityFromPressure(interface_state.p_);

				density_change_rate -= 2.0 * (rho_star * interface_state.vel_).dot(e_ij) * dW_ijV_j;
			}
            if (inner_neighborhood.boundary_type_[n] == 10)
            {
                Real U_f = 1.0; //characteristic velocity is set as 1
                Real h = 1.0;   // the height of the inflow domain size
                Vecd parabolic_velocity_inlet = Vecd::Zero();
                parabolic_velocity_inlet[0] = 1.5 * U_f * (1.0 - pos_[index_i][1] * pos_[index_i][1] / pow(0.5 * h, 2));
                Vecd vel_inlet = parabolic_velocity_inlet;
                Real p_inlet = state_i.p_;
                Real rho_inlet = state_i.rho_;
                FluidState state_j(rho_inlet, vel_inlet, p_inlet);
                FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
                Real rho_star = fluid_.DensityFromPressure(interface_state.p_);

                density_change_rate -= 2.0 * (rho_star * interface_state.vel_).dot(e_ij) * dW_ijV_j;
            }
		}
		drho_dt_[index_i] = density_change_rate;
	}
	//=================================================================================================//
	template <class RiemannSolverType>
	void BaseIntegration2ndHalfInFVM<RiemannSolverType>::update(size_t index_i, Real dt)
	{
		rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
	}
	//=================================================================================================//
}
//=================================================================================================//