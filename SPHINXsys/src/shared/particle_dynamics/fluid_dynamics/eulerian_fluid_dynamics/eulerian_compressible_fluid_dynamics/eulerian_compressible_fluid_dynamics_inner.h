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
 * @file    eulerian_compressible_fluid_dynamics_inner.h
 * @brief 	Here, we define the algorithm classes for eulerian fluid dynamics within the body.
 * @details We consider here compressible fluids.
 * @author	Zhentong Wang, Chi ZHang and Xiangyu Hu
 */

#pragma once

#include "fluid_dynamics_inner.h"

#include "all_particle_dynamics.h"
#include "all_general_dynamics.h"
#include "base_kernel.h"
#include "external_force.h"
#include "riemann_solver.h"
#include "compressible_fluid.h"

namespace SPH
{
	namespace eulerian_compressible_fluid_dynamics
	{
		typedef DataDelegateSimple<CompressibleFluidParticles> CompressibleFluidDataSimple;
		typedef DataDelegateInner<CompressibleFluidParticles> CompressibleFluidDataInner;

		class CompressibleFlowTimeStepInitialization : public BaseTimeStepInitialization, public CompressibleFluidDataSimple
		{
		protected:
			StdLargeVec<Real> &rho_, &dE_dt_prior_;
			StdLargeVec<Vecd> &pos_, &vel_, &dmom_dt_prior_;

		public:
			CompressibleFlowTimeStepInitialization(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd::Zero()));
			virtual ~CompressibleFlowTimeStepInitialization(){};

			void update(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class CompressibleFluidInitialCondition
		 * @brief  Set initial condition for a fluid body.
		 * This is a abstract class to be override for case specific initial conditions
		 */
		class CompressibleFluidInitialCondition : public LocalDynamics, public CompressibleFluidDataSimple
		{
		public:
			explicit CompressibleFluidInitialCondition(SPHBody &sph_body);
			virtual ~CompressibleFluidInitialCondition(){};

		protected:
			StdLargeVec<Vecd> &pos_, &vel_, &mom_;
			StdLargeVec<Real> &rho_, &E_, &p_;
			Real gamma_;
		};

		/**
		 * @class ViscousAccelerationInner
		 * @brief  the viscosity force induced acceleration
		 */
		class ViscousAccelerationInner : public LocalDynamics, public CompressibleFluidDataInner
		{
		public:
			explicit ViscousAccelerationInner(BaseInnerRelation &inner_relation);
			virtual ~ViscousAccelerationInner(){};
			
			inline void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Real> &Vol_, &rho_, &p_, &mass_, &dE_dt_prior_;
			StdLargeVec<Vecd> &vel_, &dmom_dt_prior_;
			Real mu_;
			Real smoothing_length_;
		};

		/**
		 * @class AcousticTimeStepSize
		 * @brief Computing the acoustic time step size
		 */
		class AcousticTimeStepSize : public LocalDynamicsReduce<Real, ReduceMax>, public CompressibleFluidDataSimple
		{
		protected:
			CompressibleFluid &compressible_fluid_;
			StdLargeVec<Real> &rho_, &p_;
			StdLargeVec<Vecd> &vel_;
			Real smoothing_length_;

		public:
			explicit AcousticTimeStepSize(SPHBody &sph_body);
			virtual ~AcousticTimeStepSize(){};

			Real reduce(size_t index_i, Real dt = 0.0);
			virtual Real outputResult(Real reduced_value) override;
		};

		/**
		 * @class BaseIntegration
		 * @brief Pure abstract base class for all fluid relaxation schemes
		 */
		class BaseIntegration : public LocalDynamics, public CompressibleFluidDataInner
		{
		public:
			explicit BaseIntegration(BaseInnerRelation &inner_relation);
			virtual ~BaseIntegration(){};

		protected:
			CompressibleFluid &compressible_fluid_;
			StdLargeVec<Real> &Vol_, &rho_, &p_, &drho_dt_, &E_, &dE_dt_, &dE_dt_prior_;
			StdLargeVec<Vecd> &vel_, &mom_, &dmom_dt_, &dmom_dt_prior_;
		};

		/**
		 * @class BaseIntegration1stHalf
		 * @brief Template class for pressure relaxation scheme with the Riemann solver
		 * as template variable
		 */
		template <class RiemannSolverType>
		class BaseIntegration1stHalf : public BaseIntegration
		{
		public:
			explicit BaseIntegration1stHalf(BaseInnerRelation &inner_relation);
			virtual ~BaseIntegration1stHalf(){};
			RiemannSolverType riemann_solver_;
			void initialization(size_t index_i, Real dt = 0.0);
			
			inline void interaction(size_t index_i, Real dt = 0.0);
			
			void update(size_t index_i, Real dt = 0.0);
		};
		using Integration1stHalfHLLCRiemann = BaseIntegration1stHalf<HLLCRiemannSolver>;
		using Integration1stHalfHLLCWithLimiterRiemann = BaseIntegration1stHalf<HLLCWithLimiterRiemannSolver>;

		/**
		 * @class BaseIntegration2ndHalf
		 * @brief  Template density relaxation scheme in HLLC Riemann solver with and without limiter
		 */
		template <class RiemannSolverType>
		class BaseIntegration2ndHalf : public BaseIntegration
		{
		public:
			explicit BaseIntegration2ndHalf(BaseInnerRelation &inner_relation);
			virtual ~BaseIntegration2ndHalf(){};
			RiemannSolverType riemann_solver_;
			
			inline void interaction(size_t index_i, Real dt = 0.0);
			
			void update(size_t index_i, Real dt = 0.0);
		};
		using Integration2ndHalfHLLCRiemann = BaseIntegration2ndHalf<HLLCRiemannSolver>;
		using Integration2ndHalfHLLCWithLimiterRiemann = BaseIntegration2ndHalf<HLLCWithLimiterRiemannSolver>;
	}
}
