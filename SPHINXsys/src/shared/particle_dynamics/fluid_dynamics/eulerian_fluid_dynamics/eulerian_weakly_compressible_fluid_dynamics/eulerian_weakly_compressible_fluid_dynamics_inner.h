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
 * @file 	eulerian_weakly_compressible_fluid_dynamics_inner.h
 * @brief 	Here, we define the algorithm classes for weakly compressible fluid dynamics within the body.
 * @details We consider here weakly compressible fluids.
 *			TODO: It seems that the eulerian and Lagrangian formulation can be merged together
 * @author	Zhentong Wang, Chi ZHang and Xiangyu Hu
 */

#pragma once

#include "fluid_dynamics_inner.h"

#include "all_particle_dynamics.h"
#include "all_general_dynamics.h"
#include "base_kernel.h"
#include "external_force.h"
#include "riemann_solver.h"

namespace SPH
{
	namespace eulerian_weakly_compressible_fluid_dynamics
	{
		typedef DataDelegateSimple<WeaklyCompressibleFluidParticles> EulerianWeaklyCompressibleFluidDataSimple;
		typedef DataDelegateInner<WeaklyCompressibleFluidParticles> EulerianWeaklyCompressibleFluidDataInner;

		class EulerianFlowTimeStepInitialization
			: public BaseTimeStepInitialization,
			  public EulerianWeaklyCompressibleFluidDataSimple
		{
		protected:
			StdLargeVec<Real> &rho_;
			StdLargeVec<Vecd> &pos_, &dmom_dt_prior_;

		public:
			EulerianFlowTimeStepInitialization(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd::Zero()));
			virtual ~EulerianFlowTimeStepInitialization(){};

			void update(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class CompressibleFluidInitialCondition
		 * @brief  Set initial condition for a fluid body.
		 * This is a abstract class to be override for case specific initial conditions
		 */
		class WeaklyCompressibleFluidInitialCondition : public LocalDynamics, public EulerianWeaklyCompressibleFluidDataSimple
		{
		public:
			explicit WeaklyCompressibleFluidInitialCondition(SPHBody &sph_body);
			virtual ~WeaklyCompressibleFluidInitialCondition(){};

		protected:
			StdLargeVec<Vecd> &pos_, &vel_, &mom_;
			StdLargeVec<Real> &rho_, &p_;
		};

		/**
		 * @class ViscousAccelerationInner
		 * @brief  the viscosity force induced acceleration
		 */
		class ViscousAccelerationInner
			: public LocalDynamics,
			  public EulerianWeaklyCompressibleFluidDataInner
		{
		public:
			explicit ViscousAccelerationInner(BaseInnerRelation &inner_relation);
			virtual ~ViscousAccelerationInner(){};
			
			inline void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Real> &Vol_, &rho_, &p_;
			StdLargeVec<Vecd> &vel_, &dmom_dt_prior_;
			Real mu_;
			Real smoothing_length_;
		};

		/**
		 * @class AcousticTimeStepSize
		 * @brief Computing the acoustic time step size
		 */
		class AcousticTimeStepSize : public LocalDynamicsReduce<Real, ReduceMax>,
									 public EulerianWeaklyCompressibleFluidDataSimple
		{
		protected:
			Fluid &fluid_;
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
		class BaseIntegration : public LocalDynamics, public EulerianWeaklyCompressibleFluidDataInner
		{
		public:
			explicit BaseIntegration(BaseInnerRelation &inner_relation);
			virtual ~BaseIntegration(){};

		protected:
			Fluid &fluid_;
			StdLargeVec<Real> &Vol_, &mass_, &rho_, &p_, &drho_dt_;
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
		using Integration1stHalf = BaseIntegration1stHalf<NoRiemannSolver>;
		/** define the mostly used pressure relaxation scheme using Riemann solver */
		using Integration1stHalfAcousticRiemann = BaseIntegration1stHalf<AcousticRiemannSolver>;
		using Integration1stHalfHLLCRiemann = BaseIntegration1stHalf<HLLCRiemannSolverInWeaklyCompressibleFluid>;
		using Integration1stHalfHLLCWithLimiterRiemann = BaseIntegration1stHalf<HLLCRiemannSolverWithLimiterInWeaklyCompressibleFluid>;

		using Integration1stHalfHLLCWithLimiterRiemann = BaseIntegration1stHalf<HLLCRiemannSolverWithLimiterInWeaklyCompressibleFluid>;

		/**
		 * @class Integration2ndHalf
		 * @brief  Template density relaxation scheme with different Riemann solver
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
		using Integration2ndHalfHLLCWithLimiterRiemann = BaseIntegration2ndHalf<HLLCRiemannSolverWithLimiterInWeaklyCompressibleFluid>;

		/**
		 * @class NonReflectiveBoundaryVariableCorrection
		 * @brief this function is applied to non_reflective flows
		 * @brief modify the velocity of particles with far-field velocity under non_reflective boundary condition
		 */
		class NonReflectiveBoundaryVariableCorrection : public LocalDynamics, public EulerianWeaklyCompressibleFluidDataInner
		{
		public:
			NonReflectiveBoundaryVariableCorrection(BaseInnerRelation &inner_relation)
				: LocalDynamics(inner_relation.getSPHBody()), EulerianWeaklyCompressibleFluidDataInner(inner_relation),
				  rho_(particles_->rho_), p_(particles_->p_), mass_(particles_->mass_), Vol_(particles_->Vol_),
				  vel_(particles_->vel_), mom_(particles_->mom_), pos_(particles_->pos_),
				  surface_indicator_(particles_->surface_indicator_)
			{
				particles_->registerVariable(n_, "NormalDirection");
			};
			virtual ~NonReflectiveBoundaryVariableCorrection(){};
			
			inline void interaction(size_t index_i, Real dt = 0.0);

		protected:
			Real p_farfield_, rho_farfield_, gamma_, sound_speed_;
			Vecd vel_farfield_;
			StdLargeVec<Real> &rho_, &p_, &mass_, &Vol_;
			StdLargeVec<Vecd> &vel_, &mom_, &pos_;
			StdLargeVec<Vecd> n_;
			StdLargeVec<int> &surface_indicator_;
		};
	}
}
