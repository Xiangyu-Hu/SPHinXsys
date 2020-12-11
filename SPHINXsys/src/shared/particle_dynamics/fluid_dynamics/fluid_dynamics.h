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
* @file fluid_dynamics.h
* @brief Here, we define the algorithm classes for fluid dynamics. 
* @details We consider here weakly compressible fluids. The algorithms may be
* different for free surface flow and flow without free surface.   
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "all_particle_dynamics.h"
#include "weakly_compressible_fluid.h"
#include "base_kernel.h"

namespace SPH
{
	namespace fluid_dynamics
	{
		typedef DataDelegateSimple<FluidBody, FluidParticles, Fluid> FluidDataDelegateSimple;

		typedef DataDelegateInner<FluidBody, FluidParticles, Fluid> FluidDataDelegateInner;

		typedef DataDelegateComplex<FluidBody, FluidParticles, Fluid, SolidBody, SolidParticles> FluidDataDelegateComplex;

		/**
		 * @class FluidInitialCondition
		 * @brief  Set initial condition for a fluid body.
		 * This is a abstract class to be override for case specific initial conditions
		 */
		class FluidInitialCondition
			: public ParticleDynamicsSimple, public FluidDataDelegateSimple
		{
		public:
			FluidInitialCondition(FluidBody* body);
			virtual ~FluidInitialCondition() {};
		protected:
			StdLargeVec<Vecd>& pos_n_, & vel_n_;
		};

		/**
		* @class RegularFreeSurfaceIndication
		* @brief  indciate the regular free surface of a fluid body.
		* Should not be used for violent free surface with wave breaking.
		*/
		class RegularFreeSurfaceIndication
			: public InteractionDynamics, public FluidDataDelegateComplex
		{
		public:
			RegularFreeSurfaceIndication(SPHBodyComplexRelation* body_complex_relation);
			virtual ~RegularFreeSurfaceIndication() {};

		protected:
			StdVec<Real> contact_inv_rho_0_;
			StdLargeVec<Real>& Vol_, & pos_div_;
			StdVec<StdLargeVec<Real>*> contact_mass_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class DensityBySummation
		* @brief  computing density by summation
		*/
		class DensityBySummation
			: public InteractionDynamics, public FluidDataDelegateComplex
		{
		public:
			DensityBySummation(SPHBodyComplexRelation* body_complex_relation);
			virtual ~DensityBySummation() {};

		protected:
			Real W0_, rho_0_, inv_sigma_0_;
			StdVec<Real> contact_inv_rho_0_;
			StdLargeVec<Real>& Vol_, & rho_n_, & mass_;
			StdVec<StdLargeVec<Real>*> contact_mass_;
			
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual Real ReinitializedDensity(Real rho_sum, Real rho_0, Real rho_n) { return rho_sum; };
		};

		/**
		* @class DensityBySummationFreeSurface
		* @brief  computing density by summation for free surface flows
		*/
		class DensityBySummationFreeSurface : public DensityBySummation
		{
		public:
			DensityBySummationFreeSurface(SPHBodyComplexRelation* body_complex_relation)
				: DensityBySummation(body_complex_relation) {};
			virtual ~DensityBySummationFreeSurface() {};
		protected:
			virtual Real ReinitializedDensity(Real rho_sum, Real rho_0, Real rho_n) override {
				return rho_sum + SMAX(0.0, (rho_n - rho_sum)) * rho_0 / rho_n;
			};
		};

		/**
		 * @class ViscousAcceleration
		 * @brief  the viscosity force induced acceleration
		 */
		class ViscousAcceleration
			: public InteractionDynamics, public FluidDataDelegateComplex
		{
		public:
			ViscousAcceleration(SPHBodyComplexRelation* body_complex_relation);
			virtual ~ViscousAcceleration() {};
		protected:
			//viscosity
			Real mu_;
			Real smoothing_length_;
			StdLargeVec<Real> &Vol_, &rho_n_, &p_;
			StdLargeVec<Vecd> &vel_n_, &dvel_dt_others_;
			StdVec<StdLargeVec<Real>*> contact_Vol_;
			StdVec<StdLargeVec<Vecd>*> contact_vel_ave_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class AngularConservativeViscousAcceleration
		 * @brief  the viscosity force induced acceleration
		 */
		class AngularConservativeViscousAcceleration : public ViscousAcceleration
		{
		public:
			AngularConservativeViscousAcceleration(SPHBodyComplexRelation* body_complex_relation)
				: ViscousAcceleration(body_complex_relation) {};
			virtual ~AngularConservativeViscousAcceleration() {};
		protected:
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class TransportVelocityCorrection
		 * @brief  transport velocity correction
		 */
		class TransportVelocityCorrection
			: public InteractionDynamics, public FluidDataDelegateComplex
		{
		public:
			TransportVelocityCorrection(SPHBodyComplexRelation* body_complex_relation, StdLargeVec<Vecd>& dvel_dt_trans);
			virtual ~TransportVelocityCorrection() {};
		protected:
			StdLargeVec<Real>& Vol_, & rho_n_;
			StdLargeVec<Vecd>& pos_n_, & dvel_dt_trans_;
			StdVec<StdLargeVec<Real>*> contact_Vol_;
			Real p_background_;

			virtual void setupDynamics(Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0)  override;
		};

		/**
		 * @class TransportVelocityStress
		 * @brief  Transport velocity induced stress
		 */
		class TransportVelocityStress : public TransportVelocityCorrection
		{
		public:
			TransportVelocityStress(SPHBodyComplexRelation* body_complex_relation, StdLargeVec<Vecd>& dvel_dt_trans);
			virtual ~TransportVelocityStress() {};

		protected:
			StdLargeVec<Vecd>& vel_n_, & dvel_dt_others_;
			virtual void Interaction(size_t index_i, Real dt)  override;
		};

		/**
		 * @class TransportVelocityFormulation
		 * @brief  transport velocity formulation including velocity correction
		 * and transport velocity stress.
		 */
		class TransportVelocityFormulation
		{
		public:
			TransportVelocityFormulation(SPHBodyComplexRelation* body_complex_relation);
			virtual ~TransportVelocityFormulation() {};

			TransportVelocityCorrection correction_;
			TransportVelocityStress stress_;
		protected:
			StdLargeVec<Vecd> dvel_dt_trans_;
		};
		/**
		 * @class TotalMechanicalEnergy
		 * @brief  Compute the total mechanical energy
		 */
		class TotalMechanicalEnergy
			: public ParticleDynamicsReduce<Real, ReduceSum<Real>>, public FluidDataDelegateSimple
		{
		public:
			explicit TotalMechanicalEnergy(FluidBody* body, Gravity* gravity);
			virtual ~TotalMechanicalEnergy() {};
		protected:
			StdLargeVec<Real>& mass_;
			StdLargeVec<Vecd>& vel_n_, & pos_n_;
			Gravity* gravity_;
			Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class AcousticTimeStepSize
		* @brief Computing the acoustic time step size
		*/
		class AcousticTimeStepSize
			: public ParticleDynamicsReduce<Real, ReduceMax>, public FluidDataDelegateSimple
		{
		public:
			explicit AcousticTimeStepSize(FluidBody* body);
			virtual ~AcousticTimeStepSize() {};
		protected:
			StdLargeVec<Real>& rho_n_, & p_;
			StdLargeVec<Vecd>& vel_n_;
			Real smoothing_length_;
			Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
			Real OutputResult(Real reduced_value) override;
		};

		/**
		* @class AdvectionTimeStepSize
		* @brief Computing the advection time step size
		*/
		class AdvectionTimeStepSize : public AcousticTimeStepSize
		{
		public:
			explicit AdvectionTimeStepSize(FluidBody* body, Real U_max);
			virtual ~AdvectionTimeStepSize() {};
		protected:
			Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
			Real OutputResult(Real reduced_value) override;
		};

		/**
		* @class VorticityInFluidField
		* @brief  compute vorticity in fluid field (without consider wall boundary effect)
		*/
		class VorticityInFluidField
			: public InteractionDynamics, public FluidDataDelegateInner
		{
		public:
			VorticityInFluidField(SPHBodyInnerRelation* body_inner_relation);
			virtual ~VorticityInFluidField() {};
		protected:
			StdLargeVec<Real>& Vol_;
			StdLargeVec<Vecd>& vel_n_;
			StdLargeVec<Vecd>	vorticity_;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class PressureRelaxationFirstHalfRiemann
		 * @brief  first half of the pressure relaxation scheme with Riemann solver
		 * computing first half step displacement, density increment and full step velocity
		 */
		class PressureRelaxationFirstHalfRiemann
			: public ParticleDynamics1Level, public FluidDataDelegateComplex
		{
		public:
			PressureRelaxationFirstHalfRiemann(SPHBodyComplexRelation* body_complex_relation);
			virtual ~PressureRelaxationFirstHalfRiemann() {};
		protected:
			StdLargeVec<Real>& Vol_, & mass_, & rho_n_, & p_, & drho_dt_;
			StdLargeVec<Vecd>& pos_n_, & vel_n_, & dvel_dt_, & dvel_dt_others_;
			StdVec<StdLargeVec<Real>*> contact_Vol_;
			StdVec<StdLargeVec<Vecd>*> contact_vel_ave_, contact_dvel_dt_ave_, contact_n_;

			virtual Real getPStar(Vecd& e_ij, Vecd& vel_i, Real p_i, Real rho_i,
				Vecd& vel_j, Real p_j, Real rho_j);

			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class PressureRelaxationFirstHalf
		* @brief  first half of the pressure relaxation scheme without using Riemann solver.
		*/
		class PressureRelaxationFirstHalf : public PressureRelaxationFirstHalfRiemann
		{
		public:
			PressureRelaxationFirstHalf(SPHBodyComplexRelation* body_complex_relation)
				: PressureRelaxationFirstHalfRiemann(body_complex_relation) {};
			virtual ~PressureRelaxationFirstHalf() {};
		protected:
			virtual Real getPStar(Vecd& e_ij, Vecd& vel_i, Real p_i, Real rho_i,
				Vecd& vel_j, Real p_j, Real rho_j) override;
		};

		/**
		 * @class PressureRelaxationSecondHalfRiemann
		 * @brief  second half of the pressure relaxation scheme with Riemann solver
		 * computing second half step displacement, density increment
		 */
		class PressureRelaxationSecondHalfRiemann
			: public PressureRelaxationFirstHalfRiemann
		{
		public:
			PressureRelaxationSecondHalfRiemann(SPHBodyComplexRelation* body_complex_relation)
				: PressureRelaxationFirstHalfRiemann(body_complex_relation) {};
			virtual ~PressureRelaxationSecondHalfRiemann() {};
		protected:
			virtual Vecd getVStar(Vecd& e_ij, Vecd& vel_i, Real p_i, Real rho_i,
				Vecd& vel_j, Real p_j, Real rho_j);
			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class PressureRelaxationSecondHalf
		* @brief  second half of the pressure relaxation scheme without using Riemann solver.
		* The difference from the free surface version is that no Riemann problem is applied
		*/
		class PressureRelaxationSecondHalf : public PressureRelaxationSecondHalfRiemann
		{
		public:
			PressureRelaxationSecondHalf(SPHBodyComplexRelation* body_complex_relation)
				: PressureRelaxationSecondHalfRiemann(body_complex_relation) {};
			virtual ~PressureRelaxationSecondHalf() {};
		protected:
			virtual Vecd getVStar(Vecd& e_ij, Vecd& vel_i, Real p_i, Real rho_i,
				Vecd& vel_j, Real p_j, Real rho_j) override;
		};

		/**
		 * @class FluidInitialCondition
		 * @brief  set initial condition for Oldroyd_B_Fluid dynamics
		 * This is a abstract class to be override for case specific initial conditions
		 */
		class Oldroyd_B_FluidInitialCondition
			: public ParticleDynamicsSimple, public FluidDataDelegateSimple
		{
		public:
			Oldroyd_B_FluidInitialCondition(FluidBody* body)
				: ParticleDynamicsSimple(body), FluidDataDelegateSimple(body) {};
			virtual ~Oldroyd_B_FluidInitialCondition() {};
		};

		/**
		* @class PressureRelaxationFirstHalfOldroyd_B
		* @brief  first half of the pressure relaxation scheme using Riemann solver.
		*/
		class PressureRelaxationFirstHalfOldroyd_B : public PressureRelaxationFirstHalfRiemann
		{
		public:
			PressureRelaxationFirstHalfOldroyd_B(SPHBodyComplexRelation* body_complex_relation);
			virtual ~PressureRelaxationFirstHalfOldroyd_B() {};
		protected:
			StdLargeVec<Matd>& tau_, & dtau_dt_;
			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class PressureRelaxationSecondHalfOldroyd_B
		* @brief  second half of the pressure relaxation scheme using Riemann solver.
		*/
		class PressureRelaxationSecondHalfOldroyd_B : public PressureRelaxationSecondHalfRiemann
		{
		public:
			PressureRelaxationSecondHalfOldroyd_B(SPHBodyComplexRelation* body_complex_relation);
			virtual ~PressureRelaxationSecondHalfOldroyd_B() {};
		protected:
			StdLargeVec<Matd>& tau_, & dtau_dt_;
			Real mu_p_, lambda_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class FlowRelaxationBuffer
		 * @brief Flow buffer in which 
		 * the particles relaxes to a given target velocity profile.
		 * This technique will be used for applying several boundary conditions,
		 * such as freestream, inflow, damping boundary conditions.
		 */
		class FlowRelaxationBuffer
			: public PartDynamicsByCell, public FluidDataDelegateSimple
		{
		public:
			FlowRelaxationBuffer(FluidBody* body, BodyPartByCell* body_part);
			virtual ~FlowRelaxationBuffer() {};
		protected:
			StdLargeVec<Vecd>& pos_n_, & vel_n_;
			/** default value is 0.1 suggests reaching  target inflow velocity in about 10 time steps */
			Real relaxation_rate_;

			/** inflow profile to be defined in applications */
			virtual Vecd getTargetVelocity(Vecd& position, Vecd& velocity) = 0;
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
			InflowBoundaryCondition(FluidBody* body, BodyPartByCell* body_part) :
				FlowRelaxationBuffer(body, body_part) {};;
			virtual ~InflowBoundaryCondition() {};
		};

		/**
		 * @class DampingBoundaryCondition
		 * @brief damping boundary condition which relaxes 
		 * the particles to zero velocity profile.
		 */
		class DampingBoundaryCondition
			: public PartDynamicsByCell, public FluidDataDelegateSimple
		{
		public:
			DampingBoundaryCondition(FluidBody* body, BodyPartByCell* body_part);
			virtual ~DampingBoundaryCondition() {};
		protected:
			StdLargeVec<Vecd>& pos_n_, & vel_n_;
			/** default value is 0.1 suggests reaching  target inflow velocity in about 10 time steps */
			Real strength_;
			Vecd damping_zone_lower_bound_;
			Vecd damping_zone_upper_bound_;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		};

		/**
		 * @class EmitterInflowCondition
		 * @brief Inflow boundary condition.
		 * The body part region is required to
		 * have parallel lower- and upper-bound surfaces.
		 */
		class EmitterInflowCondition
			: public PartDynamicsByParticle, public FluidDataDelegateSimple
		{
		public:
			explicit EmitterInflowCondition(FluidBody* body, BodyPartByParticle* body_part);
			virtual ~EmitterInflowCondition() {};
		protected:
			StdLargeVec<Real>& rho_n_, & p_;
			StdLargeVec<Vecd>& pos_n_, & vel_n_;
			/** inflow pressure condition */
			Real inflow_pressure_;
			Real rho_0_;

			/** inflow velocity profile to be defined in applications */
			virtual Vecd getTargetVelocity(Vecd& position, Vecd& velocity) = 0;
			/** inflow parameters to be defined in applications */
			virtual void SetInflowParameters() = 0;

			virtual void Update(size_t unsorted_index_i, Real dt = 0.0) override;
		};

		/**
		 * @class EmitterInflowInjecting
		 * @brief Inject particles into the computational domain.
		 */
		class EmitterInflowInjecting
			: public PartDynamicsByParticle, public FluidDataDelegateSimple
		{
		public:
			/**
			 * @brief Constructor.
			 * @param[in] fluid body.
			 * @param[in] body part by particles.
			 * @param[in] axis direction of in flow: 0, 1, 2 for x-, y- and z-axis.
			 * @param[in] direction sign of the inflow: true for positive direction.
			 */
			explicit EmitterInflowInjecting(FluidBody* body, BodyPartByParticle* body_part,
				size_t body_buffer_width, int axis_direction, bool positive);
			virtual ~EmitterInflowInjecting() {};

			/** This class is only implemented in sequential due to memory conflicts. */
			virtual void parallel_exec(Real dt = 0.0) override { exec(); };
		protected:
			StdLargeVec<Vecd>& pos_n_;
			/** the axis direction for bounding*/
			const int axis_;
			/** lower and upper bound for checking */
			Vecd body_part_lower_bound_, body_part_upper_bound_;
			/** periodic translation*/
			Vecd periodic_translation_;
			size_t body_buffer_width_;

			virtual void checkLowerBound(size_t unsorted_index_i, Real dt = 0.0);
			virtual void checkUpperBound(size_t unsorted_index_i, Real dt = 0.0);
			ParticleFunctor checking_bound_;

			virtual void Update(size_t unsorted_index_i, Real dt = 0.0) override {
				checking_bound_(unsorted_index_i, dt);
			};
		};

		/**
        * @class ViscousAccelerationWallModel
        * @brief  the viscosity force induced acceleration with wall modeling
        */
		class ViscousAccelerationWallModel : public ViscousAcceleration
		{
		public:
			StdLargeVec<Vecd> gradient_p_; /**< outer region pressure gradient */
			StdLargeVec<Matd> gradient_vel_; /**< outer region velocity gradient */

			ViscousAccelerationWallModel(SPHBodyComplexRelation* body_complex_relation);
			virtual ~ViscousAccelerationWallModel() {};
		protected:
			StdVec<StdLargeVec<Vecd>*> contact_n_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};
		/**
		 * @class FreeSurfaceProbeOnFluidBody
		 * @brief Probe the free surface profile for a fluid body part by reduced operation.
		 */
		class FreeSurfaceProbeOnFluidBody : public PartDynamicsByCellReduce<Real, ReduceMax>,
			public FluidDataDelegateSimple
		{
		public:
			FreeSurfaceProbeOnFluidBody(FluidBody* body, BodyPartByCell* body_part)
				: PartDynamicsByCellReduce<Real, ReduceMax>(body, body_part), FluidDataDelegateSimple(body),
				pos_n_(particles_->pos_n_)
			{
				initial_reference_ = 0.0;
			}
			virtual ~FreeSurfaceProbeOnFluidBody() {};
		protected:
			StdLargeVec<Vecd>& pos_n_;
			virtual void SetupReduce() override {};
			virtual Real ReduceFunction(size_t index_i, Real dt = 0.0) override { return pos_n_[index_i][1]; };
		};
	}
}
