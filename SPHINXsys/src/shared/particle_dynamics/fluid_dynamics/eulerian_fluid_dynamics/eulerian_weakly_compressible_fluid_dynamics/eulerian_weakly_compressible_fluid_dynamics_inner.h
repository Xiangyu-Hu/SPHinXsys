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
* @author	Zhentong Wang,Chi Zhang and Xiangyu Hu
*/

#pragma once

#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "riemann_solver.h"

namespace SPH
{
	namespace eulerian_weakly_compressible_fluid_dynamics
	{
		typedef DataDelegateSimple<EulerianFluidBody, WeaklyCompressibleFluidParticles, Fluid> EulerianWeaklyCompressibleFluidDataSimple;
		typedef DataDelegateInner<EulerianFluidBody, WeaklyCompressibleFluidParticles, Fluid> EulerianWeaklyCompressibleFluidDataInner;

		class EulerianFlowTimeStepInitialization
			: public ParticleDynamicsSimple,
			public EulerianWeaklyCompressibleFluidDataSimple
		{
		private:
			UniquePtrKeeper<Gravity> gravity_ptr_keeper_;

		public:
			explicit EulerianFlowTimeStepInitialization(SPHBody &sph_body);
			EulerianFlowTimeStepInitialization(SPHBody &sph_body, Gravity &gravity);
			virtual ~EulerianFlowTimeStepInitialization() {};

		protected:
			StdLargeVec<Real> &rho_n_, &mass_;
			StdLargeVec<Vecd> &pos_n_, &vel_n_, &dmom_dt_prior_;
			Gravity *gravity_;
			virtual void setupDynamics(Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class FreeSurfaceIndicationInner
		* @brief  indicate the particles near the free surface of a fluid body.
		* Note that, SPHinXsys does not require this function for simulating general free surface flow problems.
		* However, some other applications may use this function, such as transport velocity formulation,
		* for masking some function which is only applicable for the bulk of the fluid body.
		*/
		class FreeSurfaceIndicationInner
			: public InteractionDynamicsWithUpdate, public EulerianWeaklyCompressibleFluidDataInner
		{
		public:
			explicit FreeSurfaceIndicationInner(BaseBodyRelationInner &inner_relation, Real thereshold = 0.75);
			virtual ~FreeSurfaceIndicationInner() {};

		protected:
			Real thereshold_by_dimensions_;
			StdLargeVec<Real>& Vol_;
			StdLargeVec<int>& surface_indicator_;
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
			: public InteractionDynamics, public EulerianWeaklyCompressibleFluidDataInner
		{
		public:
			explicit ViscousAccelerationInner(BaseBodyRelationInner &inner_relation);
			virtual ~ViscousAccelerationInner() {};
		protected:
			Real mu_;
			Real smoothing_length_;
			StdLargeVec<Real> &Vol_, &rho_n_, &p_;
			StdLargeVec<Vecd> &vel_n_, &dmom_dt_prior_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class AcousticTimeStepSize
		* @brief Computing the acoustic time step size
		*/
		class AcousticTimeStepSize :
			public ParticleDynamicsReduce<Real, ReduceMax>, public EulerianWeaklyCompressibleFluidDataSimple
		{
		public:
			explicit AcousticTimeStepSize(EulerianFluidBody &fluid_body);
			virtual ~AcousticTimeStepSize() {};
		protected:
			StdLargeVec<Real>& rho_n_, &p_;
			StdLargeVec<Vecd>& vel_n_;
			Real smoothing_length_;
			Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
			Real OutputResult(Real reduced_value) override;
		};

		/**
		* @class VorticityInner
		* @brief  compute vorticity in the fluid field
		*/
		class VorticityInner
			: public InteractionDynamics, public EulerianWeaklyCompressibleFluidDataInner
		{
		public:
			explicit VorticityInner(BaseBodyRelationInner &inner_relation);
			virtual ~VorticityInner() {};
		protected:
			StdLargeVec<Real>& Vol_;
			StdLargeVec<Vecd>& vel_n_;
			StdLargeVec<AngularVecd> vorticity_;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class BaseRelaxation
		 * @brief Pure abstract base class for all fluid relaxation schemes
		 */
		class BaseRelaxation : public ParticleDynamics1Level, public EulerianWeaklyCompressibleFluidDataInner
		{
		public:
			explicit BaseRelaxation(BaseBodyRelationInner &inner_relation);
			virtual ~BaseRelaxation() {};
		protected:
			StdLargeVec<Real>& Vol_, &mass_, &rho_n_, &p_, &drho_dt_;
			StdLargeVec<Vecd>&vel_n_, &mom_, &dmom_dt_, &dmom_dt_prior_;
		};

		/**
		 * @class BasePressureRelaxation
		 * @brief Abstract base class for all pressure relaxation schemes
		 */
		class BasePressureRelaxation : public BaseRelaxation
		{
		public:
			explicit BasePressureRelaxation(BaseBodyRelationInner &inner_relation);
			virtual ~BasePressureRelaxation() {};
		protected:
			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class BasePressureRelaxationInner
		 * @brief Template class for pressure relaxation scheme with the Riemann solver
		 * as template variable
		 */
		template<class RiemannSolverType>
		class BasePressureRelaxationInner : public BasePressureRelaxation
		{
		public:
			explicit BasePressureRelaxationInner(BaseBodyRelationInner &inner_relation);
			virtual ~BasePressureRelaxationInner() {};
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
			virtual ~BaseDensityAndEnergyRelaxation() {};
		protected:
			virtual void Initialization(size_t index_i, Real dt = 0.0) override {};
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class DensityRelaxationInner
		 * @brief  Template density relaxation scheme with different Riemann solver
		 */
		template<class RiemannSolverType>
		class BaseDensityAndEnergyRelaxationInner : public BaseDensityAndEnergyRelaxation
		{
		public:
			explicit BaseDensityAndEnergyRelaxationInner(BaseBodyRelationInner &inner_relation);
			virtual ~BaseDensityAndEnergyRelaxationInner() {};
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
		class NonReflectiveBoundaryVariableCorrection :public InteractionDynamics, public EulerianWeaklyCompressibleFluidDataInner
		{
		public:
			NonReflectiveBoundaryVariableCorrection(BaseBodyRelationInner &inner_relation) :
				InteractionDynamics(*inner_relation.sph_body_), EulerianWeaklyCompressibleFluidDataInner(inner_relation),
				rho_n_(particles_->rho_n_), p_(particles_->p_), vel_n_(particles_->vel_n_),
				mom_(particles_->mom_), pos_n_(particles_->pos_n_), mass_(particles_->mass_), Vol_(particles_->Vol_),
				surface_indicator_(particles_->surface_indicator_)
			{
				particles_->registerAVariable<Vecd>(n_, "NormalDirection");
			};
			virtual ~NonReflectiveBoundaryVariableCorrection() {};
		protected:
			Real p_farfield_, rho_farfield_, gamma_, sound_speed_;
			Vecd vel_farfield_;
			StdLargeVec<Real>& rho_n_, &p_, &mass_, &Vol_;
			StdLargeVec<Vecd>& vel_n_, &mom_, &pos_n_;
			StdLargeVec<Vecd> n_;
			StdLargeVec<int>& surface_indicator_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};
	}
}
