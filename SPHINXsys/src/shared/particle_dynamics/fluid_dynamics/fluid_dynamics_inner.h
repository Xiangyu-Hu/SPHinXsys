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
 * @file 	fluid_dynamics_inner.h
 * @brief 	Here, we define the algorithm classes for fluid dynamics within the body.
 * @details 	We consider here weakly compressible fluids. The algorithms may be
 * 			different for free surface flow and the one without free surface.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef FLUID_DYNAMICS_INNER_H
#define FLUID_DYNAMICS_INNER_H

#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "all_body_relations.h"
#include "fluid_body.h"
#include "fluid_particles.h"
#include "weakly_compressible_fluid.h"
#include "riemann_solver.h"

namespace SPH
{
	namespace fluid_dynamics
	{
		typedef DataDelegateSimple<FluidBody, FluidParticles, Fluid> FluidDataSimple;
		typedef DataDelegateInner<FluidBody, FluidParticles, Fluid> FluidDataInner;

		/**
		 * @class FluidInitialCondition
		 * @brief  Set initial condition for a fluid body.
		 * This is a abstract class to be override for case specific initial conditions
		 */
		class FluidInitialCondition : public LocalDynamics, public FluidDataSimple
		{
		public:
			explicit FluidInitialCondition(SPHBody &sph_body);
			virtual ~FluidInitialCondition(){};

		protected:
			StdLargeVec<Vecd> &pos_, &vel_;
		};

		/**
		 * @class DensitySummationInner
		 * @brief  computing density by summation
		 */
		class DensitySummationInner : public LocalDynamics, public FluidDataInner
		{
		public:
			explicit DensitySummationInner(BaseBodyRelationInner &inner_relation);
			virtual ~DensitySummationInner(){};
			void interaction(size_t index_i, Real dt = 0.0);
			void update(size_t index_i, Real dt = 0.0);

		protected:
			Real W0_, rho0_, inv_sigma0_;
			StdLargeVec<Real> &Vol_, &rho_, &mass_, &rho_sum_;
			virtual Real ReinitializedDensity(Real rho_sum, Real rho_0, Real rho_n) { return rho_sum; };
		};

		/**
		* @class DensitySummationInnerVariableSmoothingLength
		* @brief  computing density by summation with variable smoothing length
		*/
		class DensitySummationInnerVariableSmoothingLength: public DensitySummationInner
		{
		public:
			explicit DensitySummationInnerVariableSmoothingLength(BaseBodyRelationInner& inner_relation);
			virtual ~DensitySummationInnerVariableSmoothingLength() {};

			void interaction(size_t index_i, Real dt = 0.0);
		protected:
			StdLargeVec<Real>&h_ratio_;
			StdLargeVec<Real> inv_sigma0_;
			SPHBody* sph_body_;
			virtual void setupDynamics(Real dt = 0.0) override;
		};
		
		/**
		 * @class BaseViscousAccelerationInner
		 * @brief Base class for the viscosity force induced acceleration
		 */
		class BaseViscousAccelerationInner : public LocalDynamics, public FluidDataInner
		{
		public:
			explicit BaseViscousAccelerationInner(BaseBodyRelationInner &inner_relation);
			virtual ~BaseViscousAccelerationInner(){};

		protected:
			Real mu_;
			Real smoothing_length_;
			StdLargeVec<Real> &Vol_, &rho_, &p_;
			StdLargeVec<Vecd> &vel_, &acc_prior_;
		};

		/**
		 * @class ViscousAccelerationInner
		 * @brief  the viscosity force induced acceleration
		 */
		class ViscousAccelerationInner : public BaseViscousAccelerationInner
		{
		public:
			explicit ViscousAccelerationInner(BaseBodyRelationInner &inner_relation)
				: BaseViscousAccelerationInner(inner_relation){};
			virtual ~ViscousAccelerationInner(){};
			void interaction(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class AngularConservativeViscousAccelerationInner
		 * @brief the viscosity force induced acceleration, a formulation for conserving
		 * angular momentum, to be tested for its practical applications.
		 */
		class AngularConservativeViscousAccelerationInner : public BaseViscousAccelerationInner
		{
		public:
			explicit AngularConservativeViscousAccelerationInner(BaseBodyRelationInner &inner_relation)
				: BaseViscousAccelerationInner(inner_relation){};
			virtual ~AngularConservativeViscousAccelerationInner(){};
			void interaction(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class TransportVelocityCorrectionInner
		 * @brief transport velocity correction
		 */
		class TransportVelocityCorrectionInner : public LocalDynamics, public FluidDataInner
		{
		public:
			explicit TransportVelocityCorrectionInner(BaseBodyRelationInner &inner_relation, Real coefficient = 7.0);
			virtual ~TransportVelocityCorrectionInner(){};
			virtual void setupDynamics(Real dt = 0.0) override;
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Real> &Vol_, &rho_;
			StdLargeVec<Vecd> &pos_;
			StdLargeVec<int> &surface_indicator_;
			Real p_background_;
			const Real coefficient_;
		};

		/**
		 * @class AcousticTimeStepSize
		 * @brief Computing the acoustic time step size
		 */
		class AcousticTimeStepSize : public LocalDynamicsReduce<Real, ReduceMax>, public FluidDataSimple
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
		* @class AcousticTimeStepSizeVariableSmoothingLength
		* @brief Computing the acoustic time step size with variable smoothing length
		*/
		class AcousticTimeStepSizeVariableSmoothingLength : public AcousticTimeStepSize
		{
		public:
			explicit AcousticTimeStepSizeVariableSmoothingLength(SPHBody &sph_body);
			virtual ~AcousticTimeStepSizeVariableSmoothingLength() {};
		};

		/**
		 * @class AdvectionTimeStepSizeForImplicitViscosity
		 * @brief Computing the advection time step size when viscosity is handled implicitly
		 */
		class AdvectionTimeStepSizeForImplicitViscosity : public LocalDynamicsReduce<Real, ReduceMax>, public FluidDataSimple
		{
		protected:
			Real smoothing_length_;
			StdLargeVec<Vecd> &vel_;

		public:
			explicit AdvectionTimeStepSizeForImplicitViscosity(SPHBody &sph_body, Real U_max);
			virtual ~AdvectionTimeStepSizeForImplicitViscosity(){};

			Real reduce(size_t index_i, Real dt = 0.0);
			virtual Real outputResult(Real reduced_value) override;
		};

		/**
		 * @class AdvectionTimeStepSize
		 * @brief Computing the advection time step size
		 */
		class AdvectionTimeStepSize : public AdvectionTimeStepSizeForImplicitViscosity
		{
		public:
			explicit AdvectionTimeStepSize(SPHBody &sph_body, Real U_max);
			virtual ~AdvectionTimeStepSize(){};

			Real reduce(size_t index_i, Real dt = 0.0);
		};

		/**
		* @class AdvectionTimeStepSizeVariableSmoothingLength
		* @brief Computing the advection time step size with variable smoothing length
		*/
		class AdvectionTimeStepSizeVariableSmoothingLength : public AdvectionTimeStepSize
		{
		public:
			explicit AdvectionTimeStepSizeVariableSmoothingLength(SPHBody &sph_body, Real U_max);
			virtual ~AdvectionTimeStepSizeVariableSmoothingLength() {};
		};

		/**
		 * @class VorticityInner
		 * @brief  compute vorticity in the fluid field
		 */
		class VorticityInner : public LocalDynamics, public FluidDataInner
		{
		public:
			explicit VorticityInner(BaseBodyRelationInner &inner_relation);
			virtual ~VorticityInner(){};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Real> &Vol_;
			StdLargeVec<Vecd> &vel_;
			StdLargeVec<AngularVecd> vorticity_;
		};

		/**
		 * @class BaseRelaxation
		 * @brief Pure abstract base class for all fluid relaxation schemes
		 */
		class BaseRelaxation : public LocalDynamics, public FluidDataInner
		{
		public:
			explicit BaseRelaxation(BaseBodyRelationInner &inner_relation);
			virtual ~BaseRelaxation(){};

		protected:
			StdLargeVec<Real> &Vol_, &mass_, &rho_, &p_, &drho_dt_;
			StdLargeVec<Vecd> &pos_, &vel_, &acc_, &acc_prior_;
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
			void initialization(size_t index_i, Real dt = 0.0);
			void update(size_t index_i, Real dt = 0.0);

		protected:
			virtual Vecd computeNonConservativeAcceleration(size_t index_i);
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
			void interaction(size_t index_i, Real dt = 0.0);
		};
		using PressureRelaxationInner = BasePressureRelaxationInner<NoRiemannSolver>;
		/** define the mostly used pressure relaxation scheme using Riemann solver */
		using PressureRelaxationRiemannInner = BasePressureRelaxationInner<AcousticRiemannSolver>;
		using PressureRelaxationDissipativeRiemannInner = BasePressureRelaxationInner<DissipativeRiemannSolver>;

		/**
		 * @class BaseDensityRelaxation
		 * @brief Abstract base class for all density relaxation schemes
		 */
		class BaseDensityRelaxation : public BaseRelaxation
		{
		public:
			explicit BaseDensityRelaxation(BaseBodyRelationInner &inner_relation);
			virtual ~BaseDensityRelaxation(){};
			virtual void initialization(size_t index_i, Real dt = 0.0);
			virtual void update(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class DensityRelaxationInner
		 * @brief  Template density relaxation scheme with different Riemann solver
		 */
		template <class RiemannSolverType>
		class BaseDensityRelaxationInner : public BaseDensityRelaxation
		{
		public:
			explicit BaseDensityRelaxationInner(BaseBodyRelationInner &inner_relation);
			virtual ~BaseDensityRelaxationInner(){};
			RiemannSolverType riemann_solver_;
			void interaction(size_t index_i, Real dt = 0.0);
		};
		using DensityRelaxationInner = BaseDensityRelaxationInner<NoRiemannSolver>;
		/** define the mostly used density relaxation scheme using Riemann solver */
		using DensityRelaxationRiemannInner = BaseDensityRelaxationInner<AcousticRiemannSolver>;
		using DensityRelaxationDissipativeRiemannInner = BaseDensityRelaxationInner<DissipativeRiemannSolver>;

		/**
		 * @class PressureRelaxationInnerOldroyd_B
		 * @brief Pressure relaxation scheme with the mostly used Riemann solver.
		 */
		class PressureRelaxationInnerOldroyd_B : public PressureRelaxationDissipativeRiemannInner
		{
		public:
			explicit PressureRelaxationInnerOldroyd_B(BaseBodyRelationInner &inner_relation);
			virtual ~PressureRelaxationInnerOldroyd_B(){};
			void initialization(size_t index_i, Real dt = 0.0);
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Matd> &tau_, &dtau_dt_;
		};

		/**
		 * @class DensityRelaxationInnerOldroyd_B
		 * @brief Density relaxation scheme with the mostly used Riemann solver.
		 */
		class DensityRelaxationInnerOldroyd_B : public DensityRelaxationDissipativeRiemannInner
		{
		public:
			explicit DensityRelaxationInnerOldroyd_B(BaseBodyRelationInner &inner_relation);
			virtual ~DensityRelaxationInnerOldroyd_B(){};
			void interaction(size_t index_i, Real dt = 0.0);
			void update(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Matd> &tau_, &dtau_dt_;
			Real mu_p_, lambda_;
		};
	}
}
#endif // FLUID_DYNAMICS_INNER_H