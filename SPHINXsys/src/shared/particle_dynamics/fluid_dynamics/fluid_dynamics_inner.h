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
		class FluidInitialCondition : public ParticleDynamicsSimple, public FluidDataSimple
		{
		public:
			explicit FluidInitialCondition(FluidBody &fluid_body);
			virtual ~FluidInitialCondition(){};

		protected:
			StdLargeVec<Vecd> &pos_n_, &vel_n_;
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
			  public FluidDataInner
		{
		public:
			explicit FreeSurfaceIndicationInner(BaseBodyRelationInner &inner_relation, Real thereshold = 0.75);
			virtual ~FreeSurfaceIndicationInner(){};

		protected:
			Real thereshold_by_dimensions_;
			StdLargeVec<Real> &Vol_;
			StdLargeVec<int> &surface_indicator_;
			StdLargeVec<Real> pos_div_;
			Real smoothing_length_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
        * @class SpatialTemporalFreeSurfaceIdentification
        * @brief using the spatial-temporal method to indicate the surface particles to avoid mis-judgement.
        */
		template <class FreeSurfaceIdentification>
		class SpatialTemporalFreeSurfaceIdentification : public FreeSurfaceIdentification
		{
		public:
			template <typename... ConstructorArgs>
			explicit SpatialTemporalFreeSurfaceIdentification(ConstructorArgs &&...args);
			virtual ~SpatialTemporalFreeSurfaceIdentification(){};

		protected:
			StdLargeVec<int> previous_surface_indicator_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;

			void checkNearPreviousFreeSurface(size_t index_i);
		};
		using SpatialTemporalFreeSurfaceIdentificationInner =
			SpatialTemporalFreeSurfaceIdentification<FreeSurfaceIndicationInner>;

		/**
		* @class DensitySummationInner
		* @brief  computing density by summation
		*/
		class DensitySummationInner : public InteractionDynamicsWithUpdate, public FluidDataInner
		{
		public:
			explicit DensitySummationInner(BaseBodyRelationInner &inner_relation);
			virtual ~DensitySummationInner(){};

		protected:
			Real W0_, rho0_, inv_sigma0_;
			StdLargeVec<Real> &Vol_, &rho_n_, &mass_, &rho_sum_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
			virtual Real ReinitializedDensity(Real rho_sum, Real rho_0, Real rho_n) { return rho_sum; };
		};

		/**
		 * @class DensitySummationFreeSurfaceInner
		 * @brief computing density by summation with a re-normalization for free surface flows 
		 */
		class DensitySummationFreeSurfaceInner : public DensitySummationInner
		{
		public:
			explicit DensitySummationFreeSurfaceInner(BaseBodyRelationInner &inner_relation)
				: DensitySummationInner(inner_relation){};
			virtual ~DensitySummationFreeSurfaceInner(){};

		protected:
			virtual Real ReinitializedDensity(Real rho_sum, Real rho_0, Real rho_n) override
			{
				return rho_sum + SMAX(0.0, (rho_n - rho_sum)) * rho_0 / rho_n;
			};
		};
		/**
		 * @class DensitySummationFreeStreamInner
		 * @brief the density of three-layer surface particles is calculated by DensitySummationFreeSurface,
		 * @brief and the density of other internal particles is obtained by DensitySummation.
		 * @brief applied in free stream flow without wall boundary.
		 */
		class DensitySummationFreeStreamInner : public DensitySummationFreeSurfaceInner
		{
		public:
			explicit DensitySummationFreeStreamInner(BaseBodyRelationInner &inner_relation)
				: DensitySummationFreeSurfaceInner(inner_relation),
				  surface_indicator_(*particles_->getVariableByName<indexInteger, int>("SurfaceIndicator")){};
			virtual ~DensitySummationFreeStreamInner(){};

		protected:
			StdLargeVec<int> &surface_indicator_;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class ViscousAccelerationInner
		 * @brief  the viscosity force induced acceleration
		 */
		class ViscousAccelerationInner : public InteractionDynamics, public FluidDataInner
		{
		public:
			explicit ViscousAccelerationInner(BaseBodyRelationInner &inner_relation);
			virtual ~ViscousAccelerationInner(){};

		protected:
			Real mu_;
			Real smoothing_length_;
			StdLargeVec<Real> &Vol_, &rho_n_, &p_;
			StdLargeVec<Vecd> &vel_n_, &dvel_dt_prior_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class AngularConservativeViscousAccelerationInner
		 * @brief the viscosity force induced acceleration, a formulation for conserving
		 * angular momentum, to be tested for its practical applications.
		 */
		class AngularConservativeViscousAccelerationInner : public ViscousAccelerationInner
		{
		public:
			explicit AngularConservativeViscousAccelerationInner(BaseBodyRelationInner &inner_relation)
				: ViscousAccelerationInner(inner_relation){};
			virtual ~AngularConservativeViscousAccelerationInner(){};

		protected:
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class TransportVelocityCorrectionInner
		 * @brief transport velocity correction
		 */
		class TransportVelocityCorrectionInner : public InteractionDynamics, public FluidDataInner
		{
		public:
			explicit TransportVelocityCorrectionInner(BaseBodyRelationInner &inner_relation);
			virtual ~TransportVelocityCorrectionInner(){};

		protected:
			StdLargeVec<Real> &Vol_, &rho_n_;
			StdLargeVec<Vecd> &pos_n_;
			StdLargeVec<int> &surface_indicator_;
			Real p_background_;

			virtual void setupDynamics(Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class AcousticTimeStepSize
		* @brief Computing the acoustic time step size
		*/
		class AcousticTimeStepSize : public ParticleDynamicsReduce<Real, ReduceMax>, public FluidDataSimple
		{
		public:
			explicit AcousticTimeStepSize(FluidBody &fluid_body);
			virtual ~AcousticTimeStepSize(){};

		protected:
			StdLargeVec<Real> &rho_n_, &p_;
			StdLargeVec<Vecd> &vel_n_;
			Real smoothing_length_;
			Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
			Real OutputResult(Real reduced_value) override;
		};

		/**
		* @class AdvectionTimeStepSize
		* @brief Computing the advection time step size
		*/
		class AdvectionTimeStepSize : public ParticleDynamicsReduce<Real, ReduceMax>, public FluidDataSimple
		{
		public:
			explicit AdvectionTimeStepSize(FluidBody &fluid_body, Real U_max);
			virtual ~AdvectionTimeStepSize(){};

		protected:
			Real smoothing_length_;
			StdLargeVec<Vecd> &vel_n_;
			Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
			Real OutputResult(Real reduced_value) override;
		};

		/**
		* @class AdvectionTimeStepSizeForImplicitViscosity
		* @brief Computing the advection time step size when viscosity is handled implicitly
		*/
		class AdvectionTimeStepSizeForImplicitViscosity : public AdvectionTimeStepSize
		{
		public:
			explicit AdvectionTimeStepSizeForImplicitViscosity(FluidBody &fluid_body, Real U_max);
			virtual ~AdvectionTimeStepSizeForImplicitViscosity(){};
		};

		/**
		* @class VorticityInner
		* @brief  compute vorticity in the fluid field
		*/
		class VorticityInner : public InteractionDynamics, public FluidDataInner
		{
		public:
			explicit VorticityInner(BaseBodyRelationInner &inner_relation);
			virtual ~VorticityInner(){};

		protected:
			StdLargeVec<Real> &Vol_;
			StdLargeVec<Vecd> &vel_n_;
			StdLargeVec<AngularVecd> vorticity_;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class BaseRelaxation
		 * @brief Pure abstract base class for all fluid relaxation schemes
		 */
		class BaseRelaxation : public ParticleDynamics1Level, public FluidDataInner
		{
		public:
			explicit BaseRelaxation(BaseBodyRelationInner &inner_relation);
			virtual ~BaseRelaxation(){};

		protected:
			StdLargeVec<Real> &Vol_, &mass_, &rho_n_, &p_, &drho_dt_;
			StdLargeVec<Vecd> &pos_n_, &vel_n_, &dvel_dt_, &dvel_dt_prior_;
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

		protected:
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
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

		protected:
			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
        * @class FreeStreamBoundaryVelocityCorrection
        * @brief this function is applied to freestream flows
		* @brief modify the velocity of free surface particles with far-field velocity
        */
		class FreeStreamBoundaryVelocityCorrection : public ParticleDynamicsSimple, public FluidDataInner
		{
		public:
			explicit FreeStreamBoundaryVelocityCorrection(BaseBodyRelationInner &inner_relation)
				: ParticleDynamicsSimple(*inner_relation.sph_body_),
				  FluidDataInner(inner_relation), u_ref_(1.0), t_ref_(2.0),
				  rho_ref_(material_->ReferenceDensity()), rho_sum(particles_->rho_sum_),
				  vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
				  surface_indicator_(*particles_->getVariableByName<indexInteger, int>("SurfaceIndicator")){};
			virtual ~FreeStreamBoundaryVelocityCorrection(){};

		protected:
			Real u_ref_, t_ref_, rho_ref_;
			StdLargeVec<Real> &rho_sum;
			StdLargeVec<Vecd> &vel_n_, &dvel_dt_;
			StdLargeVec<int> &surface_indicator_;

			virtual void Update(size_t index_i, Real dt = 0.0) override;
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

		protected:
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};
		using DensityRelaxationInner = BaseDensityRelaxationInner<NoRiemannSolver>;
		/** define the mostly used density relaxation scheme using Riemann solver */
		using DensityRelaxationRiemannInner = BaseDensityRelaxationInner<AcousticRiemannSolver>;

		/**
		 * @class Oldroyd_B_FluidInitialCondition
		 * @brief  set initial condition for Oldroyd_B_Fluid dynamics
		 * This is a abstract class to be override for case specific initial conditions
		 */
		class Oldroyd_B_FluidInitialCondition : public ParticleDynamicsSimple, public FluidDataSimple
		{
		public:
			explicit Oldroyd_B_FluidInitialCondition(FluidBody &fluid_body)
				: ParticleDynamicsSimple(fluid_body), FluidDataSimple(fluid_body){};
			virtual ~Oldroyd_B_FluidInitialCondition(){};
		};

		/**
		* @class PressureRelaxationRiemannInnerOldroyd_B
		* @brief Pressure relaxation scheme with the mostly used Riemann solver.
		*/
		class PressureRelaxationRiemannInnerOldroyd_B : public PressureRelaxationRiemannInner
		{
		public:
			explicit PressureRelaxationRiemannInnerOldroyd_B(BaseBodyRelationInner &inner_relation);
			virtual ~PressureRelaxationRiemannInnerOldroyd_B(){};

		protected:
			StdLargeVec<Matd> &tau_, &dtau_dt_;
			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class DensityRelaxationRiemannInnerOldroyd_B
		* @brief Density relaxation scheme with the mostly used Riemann solver.
		*/
		class DensityRelaxationRiemannInnerOldroyd_B : public DensityRelaxationRiemannInner
		{
		public:
			explicit DensityRelaxationRiemannInnerOldroyd_B(BaseBodyRelationInner &inner_relation);
			virtual ~DensityRelaxationRiemannInnerOldroyd_B(){};

		protected:
			StdLargeVec<Matd> &tau_, &dtau_dt_;
			Real mu_p_, lambda_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class FlowRelaxationBuffer
		 * @brief Flow buffer in which the particles relaxes to a given target velocity profile.
		 * This technique will be used for applying several boundary conditions,
		 * such as freestream, inflow, damping boundary conditions.
		 */
		class FlowRelaxationBuffer : public PartDynamicsByCell, public FluidDataSimple
		{
		public:
			FlowRelaxationBuffer(FluidBody &fluid_body, BodyPartByCell &body_part);
			virtual ~FlowRelaxationBuffer(){};

		protected:
			StdLargeVec<Vecd> &pos_n_, &vel_n_;
			/** default value is 0.1 suggests reaching  target inflow velocity in about 10 time steps */
			Real relaxation_rate_;

			/** inflow profile to be defined in applications */
			virtual Vecd getTargetVelocity(Vecd &position, Vecd &velocity) = 0;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class InflowBoundaryCondition
		 * @brief inflow boundary condition which relaxes
		 * the particles to a given velocity profile.
		 */
		class InflowBoundaryCondition : public FlowRelaxationBuffer
		{
		public:
			InflowBoundaryCondition(FluidBody &fluid_body, BodyPartByCell &body_part)
				: FlowRelaxationBuffer(fluid_body, body_part){};
			;
			virtual ~InflowBoundaryCondition(){};
		};

		/**
		 * @class DampingBoundaryCondition
		 * @brief damping boundary condition which relaxes 
		 * the particles to zero velocity profile.
		 */
		class DampingBoundaryCondition : public PartDynamicsByCell, public FluidDataSimple
		{
		public:
			DampingBoundaryCondition(FluidBody &fluid_body, BodyRegionByCell &body_part);
			virtual ~DampingBoundaryCondition(){};

		protected:
			StdLargeVec<Vecd> &pos_n_, &vel_n_;
			/** default value is 0.1 suggests reaching  target inflow velocity in about 10 time steps */
			Real strength_;
			BoundingBox damping_zone_bounds_;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		};

		/**
		 * @class StaticConfinementDensity
		 * @brief static confinement condition for density summation
		 */
		class StaticConfinementDensity : public PartDynamicsByCell, public FluidDataSimple
		{
		public:
			StaticConfinementDensity(FluidBody &fluid_body, NearShapeSurface &near_surface);
			virtual ~StaticConfinementDensity(){};

		protected:
			Real rho0_, inv_sigma0_;
			StdLargeVec<Real> &mass_, &rho_sum_;
			StdLargeVec<Vecd> &pos_n_;
			LevelSetShape *level_set_shape_;

			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class StaticConfinementPressureRelaxation
		 * @brief static confinement condition for pressure relaxation
		 */
		class StaticConfinementPressureRelaxation : public PartDynamicsByCell, public FluidDataSimple
		{
		public:
			StaticConfinementPressureRelaxation(FluidBody &fluid_body, NearShapeSurface &near_surface);
			virtual ~StaticConfinementPressureRelaxation(){};

		protected:
			StdLargeVec<Real> &rho_n_, &p_;
			StdLargeVec<Vecd> &pos_n_, &vel_n_, &dvel_dt_;
			LevelSetShape *level_set_shape_;
			AcousticRiemannSolver riemann_solver_;

			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class StaticConfinementDensityRelaxation
		 * @brief static confinement condition for density relaxation
		 */
		class StaticConfinementDensityRelaxation : public PartDynamicsByCell, public FluidDataSimple
		{
		public:
			StaticConfinementDensityRelaxation(FluidBody &fluid_body, NearShapeSurface &near_surface);
			virtual ~StaticConfinementDensityRelaxation(){};

		protected:
			StdLargeVec<Real> &rho_n_, &p_, &drho_dt_;
			StdLargeVec<Vecd> &pos_n_, &vel_n_;
			LevelSetShape *level_set_shape_;
			AcousticRiemannSolver riemann_solver_;

			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class StaticConfinement
		 * @brief Static confined boundary condition for complex structures.
		 */
		class StaticConfinement
		{
		public:
			StaticConfinementDensity density_summation_;
			StaticConfinementPressureRelaxation pressure_relaxation_;
			StaticConfinementDensityRelaxation density_relaxation_;

			StaticConfinement(FluidBody &fluid_body, NearShapeSurface &near_surface);
			virtual ~StaticConfinement(){};
		};

		/**
		 * @class InflowCondition
		 * @brief Inflow boundary condition which impose target velocity profile.
		 * The body part region is required to have parallel lower- and upper-bound surfaces.
		 */
		class InflowCondition : public PartSimpleDynamicsByParticle, public FluidDataSimple
		{
		public:
			explicit InflowCondition(FluidBody &fluid_body, BodyPartByParticle &body_part);
			virtual ~InflowCondition(){};

		protected:
			StdLargeVec<Vecd> &pos_n_, &vel_n_;

			/** inflow velocity profile to be defined in applications */
			virtual Vecd getTargetVelocity(Vecd &position, Vecd &velocity) = 0;
			virtual void Update(size_t unsorted_index_i, Real dt = 0.0) override;
		};

		/**
		 * @class EmitterInflowCondition
		 * @brief Inflow boundary condition imposed on an emitter, in which pressure and density profile are imposed too.
		 * The body part region is required to have parallel lower- and upper-bound surfaces.
		 */
		class EmitterInflowCondition : public InflowCondition
		{
		public:
			explicit EmitterInflowCondition(FluidBody &fluid_body, BodyPartByParticle &body_part);
			virtual ~EmitterInflowCondition(){};

		protected:
			StdLargeVec<Real> &rho_n_, &p_;
			/** inflow pressure condition */
			Real inflow_pressure_;
			Real rho0_;

			virtual void Update(size_t unsorted_index_i, Real dt = 0.0) override;
		};

		/**
		 * @class EmitterInflowInjecting
		 * @brief Inject particles into the computational domain.
		 */
		class EmitterInflowInjecting : public PartSimpleDynamicsByParticle, public FluidDataSimple
		{
		public:
			explicit EmitterInflowInjecting(FluidBody &fluid_body, BodyRegionByParticle &body_part,
											size_t body_buffer_width, int axis_direction, bool positive);
			virtual ~EmitterInflowInjecting(){};

			/** This class is only implemented in sequential due to memory conflicts. */
			virtual void parallel_exec(Real dt = 0.0) override { exec(); };

		protected:
			StdLargeVec<Vecd> &pos_n_;
			StdLargeVec<Real> &rho_n_, &p_;
			const int axis_; /**< the axis direction for bounding*/
			Vecd periodic_translation_;
			size_t body_buffer_width_;
			BoundingBox body_part_bounds_;

			virtual void checkLowerBound(size_t unsorted_index_i, Real dt = 0.0);
			virtual void checkUpperBound(size_t unsorted_index_i, Real dt = 0.0);
			ParticleFunctor checking_bound_;

			virtual void Update(size_t unsorted_index_i, Real dt = 0.0) override
			{
				checking_bound_(unsorted_index_i, dt);
			};
		};

		/**
		 * @class FreeSurfaceProbeOnFluidBody
		 * @brief Probe the free surface profile for a fluid body part by reduced operation.
		 */
		class FreeSurfaceProbeOnFluidBody : public PartDynamicsByCellReduce<Real, ReduceMax>,
											public FluidDataSimple
		{
		public:
			FreeSurfaceProbeOnFluidBody(FluidBody &fluid_body, BodyPartByCell &body_part)
				: PartDynamicsByCellReduce<Real, ReduceMax>(fluid_body, body_part), FluidDataSimple(fluid_body),
				  pos_n_(particles_->pos_n_)
			{
				quantity_name_ = "FreeSurfaceProbeOnFluidBody";
				initial_reference_ = 0.0;
			}
			virtual ~FreeSurfaceProbeOnFluidBody(){};

		protected:
			StdLargeVec<Vecd> &pos_n_;
			virtual void SetupReduce() override{};
			virtual Real ReduceFunction(size_t index_i, Real dt = 0.0) override { return pos_n_[index_i][1]; };
		};
		/**
		 * @class ColorFunctionGradientInner
		 * @brief  indicate the particles near the interface of a fluid-fluid interaction and computing norm
		 */
		class ColorFunctionGradientInner : public InteractionDynamics, public FluidDataInner
		{
		public:
			explicit ColorFunctionGradientInner(BaseBodyRelationInner &inner_relation);
			virtual ~ColorFunctionGradientInner(){};

		protected:
			Real thereshold_by_dimensions_;
			StdLargeVec<Real> &Vol_;
			StdLargeVec<int> &surface_indicator_;
			StdLargeVec<Vecd> color_grad_;
			StdLargeVec<Vecd> surface_norm_;
			StdLargeVec<Real> &pos_div_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class ColorFunctionGradientInterplationInner
		 * @brief  the viscous force induced acceleration
		 */
		class ColorFunctionGradientInterplationInner : public InteractionDynamics, public FluidDataInner
		{
		public:
			explicit ColorFunctionGradientInterplationInner(BaseBodyRelationInner &inner_relation);
			virtual ~ColorFunctionGradientInterplationInner(){};

		protected:
			Real thereshold_by_dimensions_;
			StdLargeVec<Real> &Vol_;
			StdLargeVec<int> &surface_indicator_;
			StdLargeVec<Vecd> &color_grad_;
			StdLargeVec<Vecd> &surface_norm_;
			StdLargeVec<Real> &pos_div_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class SurfaceTensionAccelerationInner
		 * @brief  the viscous force induced acceleration
		 */
		class SurfaceTensionAccelerationInner : public InteractionDynamics, public FluidDataInner
		{
		public:
			SurfaceTensionAccelerationInner(BaseBodyRelationInner &inner_relation, Real gamma);
			explicit SurfaceTensionAccelerationInner(BaseBodyRelationInner &inner_relation);
			virtual ~SurfaceTensionAccelerationInner(){};

		protected:
			Real gamma_;
			StdLargeVec<Real> &Vol_, &mass_;
			StdLargeVec<Vecd> &dvel_dt_prior_;
			StdLargeVec<int> &surface_indicator_;
			StdLargeVec<Vecd> &color_grad_;
			StdLargeVec<Vecd> &surface_norm_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};
	}
}
#endif //FLUID_DYNAMICS_INNER_H