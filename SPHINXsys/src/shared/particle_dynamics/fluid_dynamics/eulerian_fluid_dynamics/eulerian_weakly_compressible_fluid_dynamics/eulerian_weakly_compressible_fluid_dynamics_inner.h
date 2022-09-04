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
 * @file 	eulerian_weakly_compressible_fluid_dynamics_inner.h
 * @brief 	Here, we define the algorithm classes for weakly compressible fluid dynamics within the body.
 * @details 	We consider here weakly compressible fluids.
 *			TODO: It seems that the eulerian and Lagrangian formulation can be merged together
 * @author	Zhentong Wang,Chi Zhang and Xiangyu Hu
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
		typedef DataDelegateSimple<EulerianFluidBody, WeaklyCompressibleFluidParticles, Fluid> EulerianWeaklyCompressibleFluidDataSimple;
		typedef DataDelegateInner<EulerianFluidBody, WeaklyCompressibleFluidParticles, Fluid> EulerianWeaklyCompressibleFluidDataInner;

		class EulerianFlowTimeStepInitialization
			: public BaseTimeStepInitialization,
			  public EulerianWeaklyCompressibleFluidDataSimple
		{
		protected:
			StdLargeVec<Real> &rho_;
			StdLargeVec<Vecd> &pos_, &dmom_dt_prior_;

		public:
			EulerianFlowTimeStepInitialization(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd(0)));
			virtual ~EulerianFlowTimeStepInitialization(){};

			void update(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class FreeSurfaceIndicationInner
		 * @brief  indicate the particles near the free surface of a fluid body.
		 * Note that, SPHinXsys does not require this function for simulating general free surface flow problems.
		 * However, some other applications may use this function, such as transport velocity formulation,
		 * for masking some function which is only applicable for the bulk of the fluid body.
		 */
		class FreeSurfaceIndicationInner
			: public InteractionDynamicsWithUpdate,
			  public EulerianWeaklyCompressibleFluidDataInner
		{
		public:
			explicit FreeSurfaceIndicationInner(BaseBodyRelationInner &inner_relation, Real threshold = 0.75);
			virtual ~FreeSurfaceIndicationInner(){};

		protected:
			Real threshold_by_dimensions_;
			StdLargeVec<Real> &Vol_;
			StdLargeVec<int> &surface_indicator_;
			StdLargeVec<Real> pos_div_;
			Real smoothing_length_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
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
			explicit ViscousAccelerationInner(BaseBodyRelationInner &inner_relation);
			virtual ~ViscousAccelerationInner(){};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			Real mu_;
			Real smoothing_length_;
			StdLargeVec<Real> &Vol_, &rho_, &p_;
			StdLargeVec<Vecd> &vel_, &dmom_dt_prior_;
		};

		/**
		 * @class AcousticTimeStepSize
		 * @brief Computing the acoustic time step size
		 */
		class AcousticTimeStepSize : public LocalDynamicsReduce<Real, ReduceMax>,
									 public EulerianWeaklyCompressibleFluidDataSimple
		{
		protected:
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
		 * @class BaseRelaxation
		 * @brief Pure abstract base class for all fluid relaxation schemes
		 */
		class BaseRelaxation : public ParticleDynamics1Level, public EulerianWeaklyCompressibleFluidDataInner
		{
		public:
			explicit BaseRelaxation(BaseBodyRelationInner &inner_relation);
			virtual ~BaseRelaxation(){};

		protected:
			StdLargeVec<Real> &Vol_, &mass_, &rho_, &p_, &drho_dt_;
			StdLargeVec<Vecd> &vel_, &mom_, &dmom_dt_, &dmom_dt_prior_;
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
		using PressureRelaxationInner = BasePressureRelaxationInner<NoRiemannSolver>;
		/** define the mostly used pressure relaxation scheme using Riemann solver */
		using PressureRelaxationAcousticRiemannInner = BasePressureRelaxationInner<AcousticRiemannSolver>;
		using PressureRelaxationHLLCRiemannInner = BasePressureRelaxationInner<HLLCRiemannSolverInWeaklyCompressibleFluid>;
		using PressureRelaxationHLLCWithLimiterRiemannInner = BasePressureRelaxationInner<HLLCRiemannSolverWithLimiterInWeaklyCompressibleFluid>;

		/**
		 * @class BaseDensityRelaxation
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
		 * @class DensityRelaxationInner
		 * @brief  Template density relaxation scheme with different Riemann solver
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
		using DensityAndEnergyRelaxationInner = BaseDensityAndEnergyRelaxationInner<NoRiemannSolver>;
		/** define the mostly used density relaxation scheme using Riemann solver */
		using DensityAndEnergyRelaxationAcousticRiemannInner = BaseDensityAndEnergyRelaxationInner<AcousticRiemannSolver>;
		using DensityAndEnergyRelaxationHLLCRiemannInner = BaseDensityAndEnergyRelaxationInner<HLLCRiemannSolverInWeaklyCompressibleFluid>;
		using DensityAndEnergyRelaxationHLLCWithLimiterRiemannInner = BaseDensityAndEnergyRelaxationInner<HLLCRiemannSolverWithLimiterInWeaklyCompressibleFluid>;

		/**
		 * @class NonReflectiveBoundaryVariableCorrection
		 * @brief this function is applied to non_reflective flows
		 * @brief modify the velocity of particles with far-field velocity under non_reflective boundary condition
		 */
		class NonReflectiveBoundaryVariableCorrection : public LocalDynamics, public EulerianWeaklyCompressibleFluidDataInner
		{
		public:
			NonReflectiveBoundaryVariableCorrection(BaseBodyRelationInner &inner_relation)
				: LocalDynamics(inner_relation.sph_body_), EulerianWeaklyCompressibleFluidDataInner(inner_relation),
				  rho_(particles_->rho_), p_(particles_->p_), vel_(particles_->vel_),
				  mom_(particles_->mom_), pos_(particles_->pos_), mass_(particles_->mass_), Vol_(particles_->Vol_),
				  surface_indicator_(particles_->surface_indicator_)
			{
				particles_->registerVariable(n_, "NormalDirection");
			};
			virtual ~NonReflectiveBoundaryVariableCorrection(){};
			void interaction(size_t index_i, Real dt = 0.0);

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
