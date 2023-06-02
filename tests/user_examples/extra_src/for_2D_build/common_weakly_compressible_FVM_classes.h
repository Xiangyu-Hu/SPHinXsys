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
  * @file 	common_weakly_compressible_FVM_classes.h
  * @brief 	Here, we define the common weakly compressible classes for fluid dynamics in FVM.
  * @author	Zhentong Wang and Xiangyu Hu
  */

#ifndef COMMON_WEAKLY_COMPRESSIBLE_FVM_CLASSES_H
#define COMMON_WEAKLY_COMPRESSIBLE_FVM_CLASSES_H
#include "common_weakly_compressible_Eulerian_classes.h"
#include "common_shared_FVM_classes.h"
namespace SPH
{
    /**
    * @class WCAcousticTimeStepSizeInFVM
    * @brief Computing the acoustic time step size
    */
    class WCAcousticTimeStepSizeInFVM : public fluid_dynamics::AcousticTimeStepSize
    {
      protected:
        StdLargeVec<Real> &rho_, &p_;
        StdLargeVec<Vecd> &vel_;
        Fluid &fluid_;

      public:
        explicit WCAcousticTimeStepSizeInFVM(SPHBody &sph_body);
        virtual ~WCAcousticTimeStepSizeInFVM(){};
        virtual Real outputResult(Real reduced_value) override;
    };

	/**
	* @class BaseViscousAccelerationInnerInFVM
	* @brief the viscosity force induced acceleration
	*/
	template<class RiemannSolverType>
	class BaseViscousAccelerationInnerInFVM : public LocalDynamics, public DataDelegateInnerInFVM<FluidParticles>
	{
	public:
		explicit BaseViscousAccelerationInnerInFVM(BaseInnerRelationInFVM& inner_relation, Real limiter_parameter=30.0)
		  : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInnerInFVM<FluidParticles>(inner_relation),
			fluid_(particles_->fluid_), riemann_solver_(fluid_, fluid_, limiter_parameter), rho_(particles_->rho_), 
			p_(particles_->p_), vel_(particles_->vel_), acc_prior_(particles_->acc_prior_), pos_(particles_->pos_),
			dmom_dt_prior_(*particles_->getVariableByName<Vecd>("OtherMomentumChangeRate")),
			mu_(particles_->fluid_.ReferenceViscosity()){};
		virtual ~BaseViscousAccelerationInnerInFVM() {};
		void interaction(size_t index_i, Real dt = 0.0)
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
		};
	protected:
		Fluid& fluid_;
		RiemannSolverType riemann_solver_;
		StdLargeVec<Real>& rho_,& p_;
		StdLargeVec<Vecd>& vel_, & acc_prior_, & pos_,& dmom_dt_prior_;
		Real mu_;
	};
	using ViscousAccelerationRiemannInnerInFVM = BaseViscousAccelerationInnerInFVM<AcousticRiemannSolverInEulerianMethod>;

	/**
	* @class BaseRelaxationInFVM
	* @brief Pure abstract base class for all fluid relaxation schemes
	*/
	class BaseRelaxationInFVM : public LocalDynamics, public DataDelegateInnerInFVM<FluidParticles>
	{
	public:
		explicit BaseRelaxationInFVM(BaseInnerRelationInFVM &inner_relation);
		virtual ~BaseRelaxationInFVM() {};
	protected:
		Fluid& fluid_;
		StdLargeVec<Real>& p_, & rho_, & drho_dt_;
		StdLargeVec<Vecd>& vel_, & mom_, & dmom_dt_, & dmom_dt_prior_, & pos_;
		Real mu_;
	};

	/**
	* @class BaseIntegration1stHalfInFVM
	* @brief Template class for pressure relaxation scheme with the Riemann solver
	* as template variable
	*/
	template <class RiemannSolverType>
	class BaseIntegration1stHalfInFVM : public BaseRelaxationInFVM
	{
	public:
		explicit BaseIntegration1stHalfInFVM(BaseInnerRelationInFVM &inner_relation, Real limiter_parameter = 15.0)
			: BaseRelaxationInFVM(inner_relation), riemann_solver_(fluid_, fluid_, limiter_parameter) {};
		virtual ~BaseIntegration1stHalfInFVM() {};
		RiemannSolverType riemann_solver_;
		void initialization(size_t index_i, Real dt = 0.0)
		{
			rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
			p_[index_i] = fluid_.getPressure(rho_[index_i]);
		};
		void interaction(size_t index_i, Real dt = 0.0)
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
		};
		void update(size_t index_i, Real dt = 0.0)
		{
			mom_[index_i] += dmom_dt_[index_i] * dt;
			vel_[index_i] = mom_[index_i] / rho_[index_i];
		};
	};
	using Integration1stHalfAcousticRiemannInFVM = BaseIntegration1stHalfInFVM<AcousticRiemannSolverInEulerianMethod>;

	/**
	* @class BaseIntegration2ndHalfInFVM
	* @brief Template class for pressure relaxation scheme with the Riemann solver
	* as template variable
	*/
	template <class RiemannSolverType>
	class BaseIntegration2ndHalfInFVM : public BaseRelaxationInFVM
	{
	public:
		explicit BaseIntegration2ndHalfInFVM(BaseInnerRelationInFVM &inner_relation, Real limiter_parameter = 15.0)
			: BaseRelaxationInFVM(inner_relation), riemann_solver_(fluid_, fluid_, limiter_parameter) {};
		virtual ~BaseIntegration2ndHalfInFVM() {};
		RiemannSolverType riemann_solver_;
		void interaction(size_t index_i, Real dt = 0.0)
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
		};
		void update(size_t index_i, Real dt = 0.0)
		{
			rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
		};
	};
	using Integration2ndHalfAcousticRiemannInFVM = BaseIntegration2ndHalfInFVM<AcousticRiemannSolverInEulerianMethod>;

	/**
	* @class BaseFluidForceOnSolidInFVM
	* @brief Base class for computing the forces from the fluid
	*/
	class BaseForceFromFluidInFVM : public LocalDynamics, public DataDelegateInnerInFVM<FluidParticles>
	{
	public:
		explicit BaseForceFromFluidInFVM(BaseInnerRelationInFVM& inner_relation);
		virtual ~BaseForceFromFluidInFVM() {};
		StdLargeVec<Vecd>& getForceFromFluid() { return force_from_fluid_; };

	protected:
		StdLargeVec<Real>& Vol_;
		StdLargeVec<Vecd> force_from_fluid_;
	};

	/**
	* @class ViscousForceFromFluidInFVM
	* @brief Computing the viscous force from the fluid
	*/
	class ViscousForceFromFluidInFVM : public BaseForceFromFluidInFVM
	{
	public:
		explicit ViscousForceFromFluidInFVM(BaseInnerRelationInFVM& inner_relation);
		virtual ~ViscousForceFromFluidInFVM() {};
		void interaction(size_t index_i, Real dt = 0.0);
	protected:
		Fluid& fluid_;
		StdLargeVec<Vecd>& vel_;
		Real mu_;
	};

	/**
	 * @class BasePressureForceAccelerationFromFluidInFVM
	* @brief Template class fro computing the pressure force from the fluid with different Riemann solvers in FVM.
	* The pressure force is added on the viscous force of the latter is computed.
	* time step size compared to the fluid dynamics
	*/
	template <class RiemannSolverType>
	class BasePressureForceAccelerationFromFluidInFVM : public BaseForceFromFluidInFVM
	{
	public:
		explicit BasePressureForceAccelerationFromFluidInFVM(BaseInnerRelationInFVM& inner_relation)
			: BaseForceFromFluidInFVM(inner_relation), vel_(particles_->vel_), p_(particles_->p_),
			fluid_(particles_->fluid_), rho_(particles_->rho_), riemann_solver_(fluid_, fluid_)
		{
			particles_->registerVariable(force_from_fluid_, "PressureForceFromFluid");
		};
		Fluid& fluid_;
		StdLargeVec<Vecd>& vel_;
		StdLargeVec<Real>& p_, & rho_;
		RiemannSolverType riemann_solver_;
		virtual ~BasePressureForceAccelerationFromFluidInFVM() {};

		void interaction(size_t index_i, Real dt = 0.0)
		{
			Real Vol_i = Vol_[index_i];
			Vecd& vel_i = vel_[index_i];

			Vecd force = Vecd::Zero();
			const NeighborhoodInFVM& inner_neighborhood = inner_configuration_in_FVM_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
				if (inner_neighborhood.boundary_type_[n] == 3)
				{
					Vecd e_ij = inner_neighborhood.e_ij_[n];
					Vecd vel_in_wall = -vel_i;
					Real p_in_wall = p_[index_i];
					Real rho_in_wall = fluid_.DensityFromPressure(p_in_wall);
					FluidState state_j(rho_in_wall, vel_in_wall, p_in_wall);
					FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
					Real p_star = interface_state.p_;
					force -= 2.0 * (-e_ij) * p_star * Vol_i * inner_neighborhood.dW_ijV_j_[n];
				}
			}
			force_from_fluid_[index_i] = force;
		};
	};
	using PressureForceAccelerationFromFluidInFVM = BasePressureForceAccelerationFromFluidInFVM<NoRiemannSolverInEulerianMethod>;
	using PressureForceAccelerationFromFluidRiemannInFVM = BasePressureForceAccelerationFromFluidInFVM<AcousticRiemannSolverInEulerianMethod>;

	/**
	* @class BaseAllForceAccelerationFromFluidInFVM
	* @brief template class for computing force from fluid with updated viscous force in FVM
	*/
	template <class PressureForceType>
	class BaseAllForceAccelerationFromFluidInFVM : public PressureForceType
	{
	public:
		template <class ViscousForceFromFluidType>
		BaseAllForceAccelerationFromFluidInFVM(BaseInnerRelationInFVM& inner_relation,
			ViscousForceFromFluidType& viscous_force_from_fluid)
			: PressureForceType(inner_relation),
			viscous_force_from_fluid_(viscous_force_from_fluid.getForceFromFluid())
		{
			this->particles_->registerVariable(this->force_from_fluid_, "AllForceFromFluid");
		};
		virtual ~BaseAllForceAccelerationFromFluidInFVM() {};

		void interaction(size_t index_i, Real dt = 0.0)
		{
			PressureForceType::interaction(index_i, dt);
			this->force_from_fluid_[index_i] += viscous_force_from_fluid_[index_i];
		};

	protected:
		StdLargeVec<Vecd>& viscous_force_from_fluid_;
	};
	using AllForceAccelerationFromFluid = BaseAllForceAccelerationFromFluidInFVM<PressureForceAccelerationFromFluidInFVM>;
	using AllForceAccelerationFromFluidRiemann = BaseAllForceAccelerationFromFluidInFVM<PressureForceAccelerationFromFluidRiemannInFVM>;
}
#endif // COMMON_WEAKLY_COMPRESSIBLE_FVM_CLASSES_H