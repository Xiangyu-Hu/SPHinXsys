/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
* @file 	eulerian_fluid_dynamics_inner.h
* @brief 	Here, we define the algorithm classes for fluid dynamics within the body. 
* @details 	We consider here compressible fluids. 
* @author	Chi ZHang, Xiangyu Hu and Zhentong Wang
*/

#pragma once

#include "fluid_dynamics_inner.h"

#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "riemann_solver.h"

namespace SPH
{
	namespace eulerian_fluid_dynamics
	{
		typedef DataDelegateSimple<EulerianFluidBody, CompressibleFluidParticles, CompressibleFluid> CompressibleFluidDataSimple;
		typedef DataDelegateInner<EulerianFluidBody, CompressibleFluidParticles, CompressibleFluid> CompressibleFluidDataInner;

		class CompressibleFlowTimeStepInitialization
			: public ParticleDynamicsSimple,
			  public CompressibleFluidDataSimple
		{
		private:
			UniquePtrKeeper<Gravity> gravity_ptr_keeper_;

		public:
			explicit CompressibleFlowTimeStepInitialization(SPHBody &sph_body);
			CompressibleFlowTimeStepInitialization(SPHBody &sph_body, Gravity &gravity);
			virtual ~CompressibleFlowTimeStepInitialization(){};

		protected:
			StdLargeVec<Real> &rho_n_, &dE_dt_prior_;
			StdLargeVec<Vecd> &pos_n_, &vel_n_, &dmom_dt_prior_;
			Gravity *gravity_;
			virtual void setupDynamics(Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class CompressibleFluidInitialCondition
		 * @brief  Set initial condition for a fluid body.
		 * This is a abstract class to be override for case specific initial conditions
		 */
		class CompressibleFluidInitialCondition
			: public ParticleDynamicsSimple,
			  public CompressibleFluidDataSimple
		{
		public:
			explicit CompressibleFluidInitialCondition(EulerianFluidBody &body);
			virtual ~CompressibleFluidInitialCondition(){};

		protected:
			StdLargeVec<Vecd> &pos_n_, &vel_n_, &mom_;
			StdLargeVec<Real> &rho_n_, &E_, &p_;
			Real gamma_;
		};

		/**
		 * @class ViscousAccelerationInner
		 * @brief  the viscosity force induced acceleration
		 */
		class ViscousAccelerationInner
			: public InteractionDynamics,
			  public CompressibleFluidDataInner
		{
		public:
			explicit ViscousAccelerationInner(BaseBodyRelationInner &inner_relation);
			virtual ~ViscousAccelerationInner(){};

		protected:
			Real mu_;
			Real smoothing_length_;
			StdLargeVec<Real> &Vol_, &rho_n_, &p_, &mass_, &dE_dt_prior_;
			StdLargeVec<Vecd> &vel_n_, &dmom_dt_prior_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class AcousticTimeStepSize
		* @brief Computing the acoustic time step size
		*/
		class AcousticTimeStepSize : public ParticleDynamicsReduce<Real, ReduceMax>, public CompressibleFluidDataSimple
		{
		public:
			explicit AcousticTimeStepSize(EulerianFluidBody &body);
			virtual ~AcousticTimeStepSize(){};

		protected:
			StdLargeVec<Real> &rho_n_, &p_;
			StdLargeVec<Vecd> &vel_n_;
			Real smoothing_length_;
			Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
			Real OutputResult(Real reduced_value) override;
		};

		/**
		 * @class BaseRelaxation
		 * @brief Pure abstract base class for all fluid relaxation schemes
		 */
		class BaseRelaxation : public ParticleDynamics1Level, public CompressibleFluidDataInner
		{
		public:
			explicit BaseRelaxation(BaseBodyRelationInner &inner_relation);
			virtual ~BaseRelaxation(){};

		protected:
			StdLargeVec<Real> &Vol_, &rho_n_, &p_, &drho_dt_, &E_, &dE_dt_, &dE_dt_prior_;
			StdLargeVec<Vecd> &vel_n_, &mom_, &dmom_dt_, &dmom_dt_prior_;
		};

		/**
		 * @class BasePressureRelaxation
		 * @brief Abstract base class for all pressure relaxation schemes
		 */
		class BasePressureRelaxation : public BaseRelaxation
		{
		public:
			explicit BasePressureRelaxation(BaseBodyRelationInner &inner_relation);
			virtual ~BasePressureRelaxation(){};

		protected:
			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class BasePressureRelaxationInner
		 * @brief Template class for pressure relaxation scheme with the Riemann solver
		 * as template variable
		 */
		template <class RiemannSolverType>
		class BasePressureRelaxationInner : public BasePressureRelaxation
		{
		public:
			explicit BasePressureRelaxationInner(BaseBodyRelationInner &inner_relation);
			virtual ~BasePressureRelaxationInner(){};
			RiemannSolverType riemann_solver_;

		protected:
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};
		using PressureRelaxationHLLCRiemannInner = BasePressureRelaxationInner<HLLCRiemannSolver>;
		using PressureRelaxationHLLCWithLimiterRiemannInner = BasePressureRelaxationInner<HLLCWithLimiterRiemannSolver>;

		/**
		 * @class BaseDensityAndEnergyRelaxation
		 * @brief Abstract base class for all density relaxation schemes
		 */
		class BaseDensityAndEnergyRelaxation : public BaseRelaxation
		{
		public:
			explicit BaseDensityAndEnergyRelaxation(BaseBodyRelationInner &inner_relation);
			virtual ~BaseDensityAndEnergyRelaxation(){};

		protected:
			virtual void Initialization(size_t index_i, Real dt = 0.0) override{};
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class BaseDensityAndEnergyRelaxationInner
		 * @brief  Template density relaxation scheme in HLLC Riemann solver with and without limiter
		 */
		template <class RiemannSolverType>
		class BaseDensityAndEnergyRelaxationInner : public BaseDensityAndEnergyRelaxation
		{
		public:
			explicit BaseDensityAndEnergyRelaxationInner(BaseBodyRelationInner &inner_relation);
			virtual ~BaseDensityAndEnergyRelaxationInner(){};
			RiemannSolverType riemann_solver_;

		protected:
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};
		using DensityAndEnergyRelaxationHLLCRiemannInner = BaseDensityAndEnergyRelaxationInner<HLLCRiemannSolver>;
		using DensityAndEnergyRelaxationHLLCWithLimiterRiemannInner = BaseDensityAndEnergyRelaxationInner<HLLCWithLimiterRiemannSolver>;
	}
}
