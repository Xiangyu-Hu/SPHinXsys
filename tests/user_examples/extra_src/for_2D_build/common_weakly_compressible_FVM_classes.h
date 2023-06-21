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
#include "common_weakly_compressible_eulerian_classes.h"
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
	class BaseViscousAccelerationInnerInFVM : public LocalDynamics, public DataDelegateInnerInFVM<BaseParticles>
	{
	public:
		explicit BaseViscousAccelerationInnerInFVM(BaseInnerRelationInFVM& inner_relation, Real limiter_parameter=30.0);
		virtual ~BaseViscousAccelerationInnerInFVM() {};
		void interaction(size_t index_i, Real dt = 0.0);
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
	class BaseRelaxationInFVM : public LocalDynamics, public DataDelegateInnerInFVM<BaseParticles>
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
		explicit BaseIntegration1stHalfInFVM(BaseInnerRelationInFVM &inner_relation, Real limiter_parameter = 15.0);
		virtual ~BaseIntegration1stHalfInFVM() {};
		RiemannSolverType riemann_solver_;
		void initialization(size_t index_i, Real dt = 0.0);
		void interaction(size_t index_i, Real dt = 0.0);
		void update(size_t index_i, Real dt = 0.0);
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
		explicit BaseIntegration2ndHalfInFVM(BaseInnerRelationInFVM &inner_relation, Real limiter_parameter = 15.0);
		virtual ~BaseIntegration2ndHalfInFVM() {};
		RiemannSolverType riemann_solver_;
		void interaction(size_t index_i, Real dt = 0.0);
		void update(size_t index_i, Real dt = 0.0);
	};
	using Integration2ndHalfAcousticRiemannInFVM = BaseIntegration2ndHalfInFVM<AcousticRiemannSolverInEulerianMethod>;

	/**
	* @class BaseFluidForceOnSolidInFVM
	* @brief Base class for computing the forces from the fluid
	*/
	class BaseForceFromFluidInFVM : public LocalDynamics, public DataDelegateInnerInFVM<BaseParticles>
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
			: BaseForceFromFluidInFVM(inner_relation),  fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())), vel_(particles_->vel_),
			p_(*particles_->getVariableByName<Real>("Pressure")), rho_(particles_->rho_), riemann_solver_(fluid_, fluid_)
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
	using PressureForceAccelerationFromFluidInFVM = BasePressureForceAccelerationFromFluidInFVM<NoRiemannSolverInWCEulerianMethod>;
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