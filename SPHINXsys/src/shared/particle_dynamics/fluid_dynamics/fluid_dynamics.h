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
		typedef ParticleDynamics<void, FluidBody, FluidParticles, WeaklyCompressibleFluid>
			WeaklyCompressibleFluidDynamics;

		typedef ParticleDynamicsSimple<FluidBody, FluidParticles, WeaklyCompressibleFluid>
			WeaklyCompressibleFluidDynamicsSimple;

		typedef ParticleDynamicsInner<FluidBody, FluidParticles, WeaklyCompressibleFluid>  
			WeaklyCompressibleFluidDynamicsInner;

		typedef ParticleDynamicsComplex<FluidBody, FluidParticles,
			WeaklyCompressibleFluid, SolidBody, SolidParticles> WeaklyCompressibleFluidDynamicsComplex;

		typedef ParticleDynamicsComplexWithUpdate<FluidBody, FluidParticles,
			WeaklyCompressibleFluid, SolidBody, SolidParticles> WeaklyCompressibleFluidDynamicsComplexWithUpdate;

		typedef ParticleDynamicsComplex1Level<FluidBody, FluidParticles,
			WeaklyCompressibleFluid, SolidBody, SolidParticles> WeaklyCompressibleFluidDynamicsComplex1Level;

		template <class ReturnType>
		using WeaklyCompressibleFluidDynamicsSum = ParticleDynamicsReduce<ReturnType, ReduceSum<ReturnType>, FluidBody,
			FluidParticles, WeaklyCompressibleFluid>;

		typedef ParticleDynamicsReduce<Real, ReduceMax, FluidBody,
			FluidParticles, WeaklyCompressibleFluid> WeaklyCompressibleFluidDynamicsMaximum;

		typedef PartDynamicsByCell<FluidBody,
			FluidParticles, BodyPartByCell> WeaklyCompressibleFluidPartDynamicsByCell;

		typedef PartDynamicsByParticle<FluidBody,	FluidParticles, BodyPartByParticle, WeaklyCompressibleFluid> 
			WeaklyCompressibleFluidPartDynamicsByParticle;

		typedef ParticleDynamicsSimple<FluidBody, ViscoelasticFluidParticles, Oldroyd_B_Fluid>
			Oldroyd_B_FluidDynamicsSimple;
			
		typedef ParticleDynamicsComplexSplitting<FluidBody, FluidParticles,
			WeaklyCompressibleFluid, SolidBody, SolidParticles> SplittingFluidDynamicsComplex;
		/**
		 * @class WeaklyCompressibleFluidInitialCondition
		 * @brief  Set default initial condition for a fluid body with weakly compressible fluid.
		 * This is a abstract class to be override for case specific initial conditions
		 */
		class WeaklyCompressibleFluidInitialCondition : public WeaklyCompressibleFluidDynamicsSimple
		{
		public:
			WeaklyCompressibleFluidInitialCondition(FluidBody *body)
				: WeaklyCompressibleFluidDynamicsSimple(body) {};
			virtual ~WeaklyCompressibleFluidInitialCondition() {};
		};

		/**
		* @class DensityBySummation
		* @brief  computing density by summation
		*/
		class DensityBySummation : public WeaklyCompressibleFluidDynamicsComplex
		{
		protected:
			Real W0_;
			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual Real reinitializeDensity(Real rho_sum, Real rho_0, Real rho_n) { return rho_sum; };
		public:
			DensityBySummation(SPHBodyComplexRelation* body_complex_relation)
				: WeaklyCompressibleFluidDynamicsComplex(body_complex_relation) {
				W0_ = body_->kernel_->W(Vecd(0));
			};
			virtual ~DensityBySummation() {};
		};

		/**
		* @class DensityBySummationFreeSurface
		* @brief  computing density by summation for free surface flows
		*/
		class DensityBySummationFreeSurface : public DensityBySummation
		{
		protected:
			virtual Real reinitializeDensity(Real rho_sum, Real rho_0, Real rho_n) override {	
				return rho_sum + SMAX(0.0, (rho_n - rho_sum)) * rho_0 / rho_n;
			};
		public:
			DensityBySummationFreeSurface(SPHBodyComplexRelation* body_complex_relation)
				: DensityBySummation(body_complex_relation) {};
			virtual ~DensityBySummationFreeSurface() {};
		};

		/**
		 * @class ComputingViscousAcceleration
		 * @brief  the viscosity force induced acceleration
		 */
		class ComputingViscousAcceleration : public WeaklyCompressibleFluidDynamicsComplex
		{
		protected:
			//viscosity
			Real mu_;
			Real smoothing_length_;

			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			ComputingViscousAcceleration(SPHBodyComplexRelation* body_complex_relation)
				: WeaklyCompressibleFluidDynamicsComplex(body_complex_relation) {
				mu_ = material_->getReferenceViscosity();
				smoothing_length_ = body_->kernel_->GetSmoothingLength();
			};
			virtual ~ComputingViscousAcceleration() {};
		};

		/**
		 * @class ComputingAngularConservativeViscousAcceleration
		 * @brief  the viscosity force induced acceleration
		 */
		class ComputingAngularConservativeViscousAcceleration : public WeaklyCompressibleFluidDynamicsComplex
		{
		protected:
			//viscosity
			Real mu_;
			Real smoothing_length_;

			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			ComputingAngularConservativeViscousAcceleration(SPHBodyComplexRelation* body_complex_relation)
				: WeaklyCompressibleFluidDynamicsComplex(body_complex_relation) {
				mu_ = material_->getReferenceViscosity();
				smoothing_length_ = body_->kernel_->GetSmoothingLength();
			};
			virtual ~ComputingAngularConservativeViscousAcceleration() {};
		};
		
		/**
		 * @class TransportVelocityCorrection
		 * @brief  transport velocity correction
		 */
		class TransportVelocityCorrection : public WeaklyCompressibleFluidDynamicsComplex
		{
		public:
			TransportVelocityCorrection(SPHBodyComplexRelation* body_complex_relation, StdLargeVec<Vecd>& dvel_dt_trans)
				: WeaklyCompressibleFluidDynamicsComplex(body_complex_relation), dvel_dt_trans_(dvel_dt_trans),
				p_background_(0) {};
			virtual ~TransportVelocityCorrection() {};

		protected:
			StdLargeVec<Vecd>& dvel_dt_trans_;
			Real p_background_;

			virtual void setupDynamics(Real dt = 0.0) override;
			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0)  override;
		};

		/**
		 * @class TransportVelocityStress
		 * @brief  Transport velocity induced stress
		 */
		class TransportVelocityStress : public TransportVelocityCorrection
		{
		public:
			TransportVelocityStress(SPHBodyComplexRelation* body_complex_relation, StdLargeVec<Vecd>& dvel_dt_trans)
				: TransportVelocityCorrection(body_complex_relation, dvel_dt_trans) {};
			virtual ~TransportVelocityStress() {};

		protected:
			virtual void ComplexInteraction(size_t index_particle_i, Real dt)  override;
		};

		/**
		 * @class TransportVelocityFormulation
		 * @brief  transport velocity formulation including velocity correction
		 * and transport veocity stress.
		 */
		class TransportVelocityFormulation : public WeaklyCompressibleFluidDynamics
		{
		public:
			TransportVelocityFormulation(SPHBodyComplexRelation* body_complex_relation);
			virtual ~TransportVelocityFormulation() {};
	
			TransportVelocityCorrection correction_;
			TransportVelocityStress stress_;

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		protected:
			StdLargeVec<Vecd> dvel_dt_trans_;
		};
		/**
		 * @class TotalMechanicalEnergy
		 * @brief  Compute the total mechanical energy
		 */
		class TotalMechanicalEnergy  : public WeaklyCompressibleFluidDynamicsSum<Real>
		{
		protected:
			Gravity* gravity_;
			Real ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit TotalMechanicalEnergy(FluidBody* body, Gravity* gravity);
			virtual ~TotalMechanicalEnergy() {};
		};

		/**
		* @class GetAcousticTimeStepSize
		* @brief Computing the acoustic time step size
		* computing time step size
		*/
		class GetAcousticTimeStepSize : public WeaklyCompressibleFluidDynamicsMaximum
		{
		protected:
			Real smoothing_length_;
			Real ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;
			Real OutputResult(Real reduced_value) override;
		public:
			explicit GetAcousticTimeStepSize(FluidBody* body);
			virtual ~GetAcousticTimeStepSize() {};
		};

		/**
		* @class GetAdvectionTimeStepSize
		* @brief Computing the advection time step size
		* computing time step size
		*/
		class GetAdvectionTimeStepSize 	: public GetAcousticTimeStepSize
		{
		protected:
			Real ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;
			Real OutputResult(Real reduced_value) override;
		public:
			explicit GetAdvectionTimeStepSize(FluidBody* body, Real U_max);
			virtual ~GetAdvectionTimeStepSize() {};
		};

		/**
		* @class ComputingVorticityInFluidField
		* @brief  compute vorticity in fluid field (without consider wall boundary effect)
		*/
		class ComputingVorticityInFluidField : public WeaklyCompressibleFluidDynamicsInner
		{
		protected:
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			ComputingVorticityInFluidField(SPHBodyInnerRelation* body_inner_relation)
				: WeaklyCompressibleFluidDynamicsInner(body_inner_relation) {};
			virtual ~ComputingVorticityInFluidField() {};
		};

		/**
		 * @class PressureRelaxationFirstHalfRiemann
		 * @brief  first half of the pressure relaxation scheme with Riemann solver
		 * computing first half step displacement, density increment and full step velocity
		 */
		class PressureRelaxationFirstHalfRiemann : public WeaklyCompressibleFluidDynamicsComplex1Level
		{
		protected:
			virtual Real getPStar(Vecd& e_ij, Vecd& vel_i, Real p_i, Real rho_i,
				Vecd& vel_j, Real p_j, Real rho_j);

			virtual void Initialization(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			PressureRelaxationFirstHalfRiemann(SPHBodyComplexRelation* body_complex_relation)
				: WeaklyCompressibleFluidDynamicsComplex1Level(body_complex_relation) {};
			virtual ~PressureRelaxationFirstHalfRiemann() {};
		};

		/**
		* @class PressureRelaxationFirstHalf
		* @brief  first half of the pressure relaxation scheme without using Riemann solver.
		*/
		class PressureRelaxationFirstHalf : public PressureRelaxationFirstHalfRiemann
		{
		protected:
			virtual Real getPStar(Vecd& e_ij, Vecd& vel_i, Real p_i, Real rho_i,
				Vecd& vel_j, Real p_j, Real rho_j) override;
		public:
			PressureRelaxationFirstHalf(SPHBodyComplexRelation* body_complex_relation)
				: PressureRelaxationFirstHalfRiemann(body_complex_relation) {};
			virtual ~PressureRelaxationFirstHalf() {};
		};

		/**
		 * @class PressureRelaxationSecondHalfRiemann
		 * @brief  second half of the pressure relaxation scheme with Riemann solver
		 * computing second half step displacement, density increment
		 */
		class PressureRelaxationSecondHalfRiemann : public WeaklyCompressibleFluidDynamicsComplex1Level
		{
		protected:
			virtual Vecd getVStar(Vecd& e_ij, Vecd& vel_i, Real p_i, Real rho_i,
				Vecd& vel_j, Real p_j, Real rho_j);
			virtual void Initialization(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			PressureRelaxationSecondHalfRiemann(SPHBodyComplexRelation* body_complex_relation)
				: WeaklyCompressibleFluidDynamicsComplex1Level(body_complex_relation) {};
			virtual ~PressureRelaxationSecondHalfRiemann() {};
		};

		/**
		* @class PressureRelaxationSecondHalf
		* @brief  second half of the pressure relaxation scheme without using Riemann solver.
		* The difference from the free surface version is that no Riemann problem is applied
		*/
		class PressureRelaxationSecondHalf : public PressureRelaxationSecondHalfRiemann
		{
		protected:
			virtual Vecd getVStar(Vecd& e_ij, Vecd& vel_i, Real p_i, Real rho_i,
				Vecd& vel_j, Real p_j, Real rho_j) override;
		public:
			PressureRelaxationSecondHalf(SPHBodyComplexRelation* body_complex_relation)
				: PressureRelaxationSecondHalfRiemann(body_complex_relation) {};
			virtual ~PressureRelaxationSecondHalf() {};
		};

		/**
		 * @class WeaklyCompressibleFluidInitialCondition
		 * @brief  set initial condition for Oldroyd_B_Fluid dynamics
		 * This is a abstract class to be override for case specific initial conditions
		 */
		class Oldroyd_B_FluidInitialCondition : public Oldroyd_B_FluidDynamicsSimple
		{
		public:
			Oldroyd_B_FluidInitialCondition(FluidBody *body)
				: Oldroyd_B_FluidDynamicsSimple(body) {};
			virtual ~Oldroyd_B_FluidInitialCondition() {};
		};

		/**
		* @class PressureRelaxationFirstHalfOldroyd_B
		* @brief  first half of the pressure relaxation scheme without using Riemann solver.
		*/
		class PressureRelaxationFirstHalfOldroyd_B : public PressureRelaxationFirstHalfRiemann
		{
		protected:
			ViscoelasticFluidParticles *viscoelastic_fluid_particles_;
			virtual void Initialization(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			PressureRelaxationFirstHalfOldroyd_B(SPHBodyComplexRelation* body_complex_relation);
			virtual ~PressureRelaxationFirstHalfOldroyd_B() {};
		};

		/**
		* @class PressureRelaxationSecondHalfOldroyd_B
		* @brief  second half of the pressure relaxation scheme without using Riemann solver.
		* The difference from the free surface version is that no Riemann problem is applied
		*/
		class PressureRelaxationSecondHalfOldroyd_B : public PressureRelaxationSecondHalfRiemann
		{
		protected:
			ViscoelasticFluidParticles* viscoelastic_fluid_particles_;
			Real mu_p_, lambda_;

			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			PressureRelaxationSecondHalfOldroyd_B(SPHBodyComplexRelation* body_complex_relation);
			virtual ~PressureRelaxationSecondHalfOldroyd_B() {};
		};

		/**
		 * @class InflowBoundaryCondition
		 * @brief inflow boundary condition which relaxes 
		 * the particles to a given velocity profile.
		 */
		class InflowBoundaryCondition : public WeaklyCompressibleFluidPartDynamicsByCell
		{
		protected:
			/** default value is 0.1 suggests reaching  target inflow velocity in about 10 time steps */
			Real constrain_strength_;

			/** inflow profile to be defined in applications */
			virtual Vecd GetInflowVelocity(Vecd& position, Vecd& velocity) = 0;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			InflowBoundaryCondition(FluidBody* body, BodyPartByCell* body_part)
				: WeaklyCompressibleFluidPartDynamicsByCell(body, body_part),
				constrain_strength_(0.1) {};
			virtual ~InflowBoundaryCondition() {};
		};

		/**
		 * @class DampingBoundaryCondition
		 * @brief damping boundary condition which relaxes 
		 * the particles to zero velocity profile.
		 */
		class DampingBoundaryCondition : public WeaklyCompressibleFluidPartDynamicsByCell
		{
		protected:
			/** dedault value is 0.1 suggests reaching  target inflow velocity in about 10 time steps */
			Real strength_;
			Vecd damping_zone_lower_bound_;
			Vecd damping_zone_upper_bound_;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			DampingBoundaryCondition(FluidBody* body, BodyPartByCell* body_part)
				: WeaklyCompressibleFluidPartDynamicsByCell(body, body_part),
				strength_(5.0) 
			{
				body_part->BodyPartBounds(damping_zone_lower_bound_, damping_zone_upper_bound_);
			};
			virtual ~DampingBoundaryCondition() {};
		};

		/**
		 * @class EmitterInflowCondition
		 * @brief Inflow boundary condition.
		 * The body part region is required to 
		 * have parallel lower- and upper-bound surfaces. 
		 */
		class EmitterInflowCondition : public WeaklyCompressibleFluidPartDynamicsByParticle
		{
		protected:
			/** inflow pressure condition */
			Real inflow_pressure_;

			/** inflow velocity profile to be defined in applications */
			virtual Vecd GetInflowVelocity(Vecd& position, Vecd& velocity) = 0;
			/** inflow parameters to be defined in applications */
			virtual void SetInflowParameters() = 0;

			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit EmitterInflowCondition(FluidBody* body, BodyPartByParticle* body_part)
				: WeaklyCompressibleFluidPartDynamicsByParticle(body, body_part),
				inflow_pressure_(0) {};
			virtual ~EmitterInflowCondition() {};
		};

		/**
		 * @class EmitterInflowInjecting
		 * @brief Inject particles into the computational domain.
		 */
		class EmitterInflowInjecting : public WeaklyCompressibleFluidPartDynamicsByParticle
		{
		protected:
			/** the axis direction for bounding*/
			const int axis_;
			/** lower and upper bound for checking */
			Vecd body_part_lower_bound_, body_part_upper_bound_;
			/** periodic translation*/
			Vecd periodic_translation_;
			size_t body_buffer_size_;

			virtual void CheckLowerBound(size_t index_particle_i, Real dt = 0.0);
			virtual void CheckUpperBound(size_t index_particle_i, Real dt = 0.0);
			InnerFunctor checking_bound_;

			virtual void Update(size_t index_particle_i, Real dt = 0.0) override {
				checking_bound_(index_particle_i, dt); };
		public:
			/**
			 * @brief Constructor.
			 * @param[in] fluid body.
			 * @param[in] body part by particles.
			 * @param[in] axis direction of in flow: 0, 1, 2 for x-, y- and z-axis.
			 * @param[in] direction sign of the inlfow: ture for positive direction.
			 */
			explicit EmitterInflowInjecting(FluidBody* body, BodyPartByParticle* body_part,
				size_t body_buffer_size, int axis_direction, bool positive);
			virtual ~EmitterInflowInjecting() {};

			/** This class is only implemented in sequential due to memory conflicts. */
			virtual void parallel_exec(Real dt = 0.0) override { exec(); };
		};
		/**
		 * @class ImplicitComputingViscousAcceleration
		 * @brief  compute the viscous acceleration with implicit algorithm with splitting cell method.
		 */
		class ImplicitComputingViscousAcceleration : public WeaklyCompressibleFluidDynamicsComplex1Level
		{
		protected:
			//viscosity
			Real mu_;
			Real smoothing_length_;

			virtual void Initialization(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			ImplicitComputingViscousAcceleration(SPHBodyComplexRelation* body_complex_relation)
				: WeaklyCompressibleFluidDynamicsComplex1Level(body_complex_relation)
			{
				mu_ = material_->getReferenceViscosity();
				smoothing_length_ = body_->kernel_->GetSmoothingLength();
			};
			virtual ~ImplicitComputingViscousAcceleration() {};
		};
	}
}