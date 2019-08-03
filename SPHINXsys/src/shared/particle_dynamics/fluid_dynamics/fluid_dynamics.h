/**
* @file fluid_dynamics.h
* @brief Here, we define the algorithm classes for fluid dynamics. 
* @details We consider here a weakly compresible fluids. The algorithms may be
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
		typedef ParticleDynamicsSimple<FluidBody, FluidParticles, WeaklyCompressibleFluid> 
			WeaklyCompressibleFluidDynamicsSimple;
		typedef ParticleDynamicsComplex<FluidBody, FluidParticles,
			WeaklyCompressibleFluid, SolidBody, SolidParticles> WeaklyCompressibleFluidDynamicsComplex;
		typedef ParticleDynamicsComplexWithUpdate<FluidBody, FluidParticles,
			WeaklyCompressibleFluid, SolidBody, SolidParticles> WeaklyCompressibleFluidDynamicsComplexWithUpdate;
		template <class ReturnType>
		using WeaklyCompressibleFluidDynamicsSum = ParticleDynamicsReduce<ReturnType, ReduceSum<ReturnType>, FluidBody,
			FluidParticles, WeaklyCompressibleFluid>;
		typedef ParticleDynamicsReduce<Real, ReduceMax, FluidBody,
			FluidParticles, WeaklyCompressibleFluid> WeaklyCompressibleFluidDynamicsMaximum;
		typedef ParticleDynamicsInner<FluidBody, FluidParticles, WeaklyCompressibleFluid>  WeaklyCompressibleFluidDynamicsInner;
		typedef ParticleDynamicsComplex2Levels<FluidBody, FluidParticles,
			WeaklyCompressibleFluid, SolidBody, SolidParticles> WeaklyCompressibleFluidDynamicsComplex2Levels;
		typedef EulerianConstraint<FluidBody,
			FluidParticles, FluidBodyPart> WeaklyCompressibleFluidConstraint;

		typedef ParticleDynamicsSimple<FluidBody, ViscoelasticFluidParticles, Oldroyd_B_Fluid>
			Oldroyd_B_FluidDynamicsSimple;
		typedef ParticleDynamicsComplex2Levels<FluidBody, ViscoelasticFluidParticles,
			Oldroyd_B_Fluid, SolidBody, SolidParticles> Oldroyd_B_FluidDynamicsComplex2Levels;

		/**
		 * @class WeaklyCompressibleFluidInitialCondition
		 * @brief  Set default initial condition for a fluid body with weakly compressible fluid.
		 * The default is to set all particle at rest
		 */
		class WeaklyCompressibleFluidInitialCondition : public WeaklyCompressibleFluidDynamicsSimple
		{
		protected:
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			WeaklyCompressibleFluidInitialCondition(FluidBody *body)
				: WeaklyCompressibleFluidDynamicsSimple(body) {};
			virtual ~WeaklyCompressibleFluidInitialCondition() {};
		};

		/**
		 * @class ResetFluidCondition
		 * @brief  Set user defined condition for a fluid body with weakly compressible fluid.
		 */
		class ResetFluidCondition : public WeaklyCompressibleFluidDynamicsSimple
		{
		protected:
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			ResetFluidCondition(FluidBody *body)
				: WeaklyCompressibleFluidDynamicsSimple(body) {};
			virtual ~ResetFluidCondition() {};
		};

		/**
		 * @class InitialNumberDensity
		 * @brief  Initialize particle number density
		 */
		class InitialNumberDensity : public WeaklyCompressibleFluidDynamicsComplex
		{
		protected:
			Real W0_;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;
		public:
			InitialNumberDensity(FluidBody *body, StdVec<SolidBody*> interacting_bodies)
				: WeaklyCompressibleFluidDynamicsComplex(body, interacting_bodies) {
				W0_ = body->kernel_->W(Vecd(0)); };
			virtual ~InitialNumberDensity() {};
		};

		/**
		 * @class DensityBySummation
		 * @brief  compute fluid particle density by summation, used for internal flows
		 */
		class DensityBySummation : public WeaklyCompressibleFluidDynamicsComplexWithUpdate
		{
		protected:
			Real W0_;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			DensityBySummation(FluidBody *body, StdVec<SolidBody*> interacting_bodies)
			: WeaklyCompressibleFluidDynamicsComplexWithUpdate(body, interacting_bodies) { W0_ = body->kernel_->W(Vecd(0)); };
			virtual ~DensityBySummation() {};
		};

		/**
		 * @class DensityBySummationFreeSurface
		 * @brief  compute fluid particle density by summation, 
		 * used for free surface flows
		 */
		class DensityBySummationFreeSurface : public DensityBySummation
		{
		protected:
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ContactInteraction(size_t index_particle_i, size_t contact_body_index, Real dt = 0.0) override;

		public:
			DensityBySummationFreeSurface(FluidBody *body, StdVec<SolidBody*> interacting_bodies)
				: DensityBySummation(body, interacting_bodies) {};
			virtual ~DensityBySummationFreeSurface() {};

		};

		/**
		 * @class DivergenceCorrection
		 * @brief  obtained divergence correction factor for each fluid particle
		 * so the divegnce of linear velocity field is reproduced
		 * note that this correction does not work for free surface flow
		 */
		class DivergenceCorrection : public WeaklyCompressibleFluidDynamicsComplexWithUpdate
		{
		protected:
			Real dimension_;

			virtual void InnerInteraction(size_t index_particle_i,	Real dt = 0.0) override;
			virtual void ContactInteraction(size_t index_particle_i, size_t contact_body_index, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			DivergenceCorrection(FluidBody *body, StdVec<SolidBody*> interacting_bodies)
				: WeaklyCompressibleFluidDynamicsComplexWithUpdate(body, interacting_bodies) {	
				dimension_ = Real(Vecd(0).size());	};
		};

		/**
		 * @class ComputingViscousAcceleration
		 * @brief  the viscosity force induecd accelerarion
		 */
		class ComputingViscousAcceleration : public WeaklyCompressibleFluidDynamicsComplex
		{
		protected:
			//viscousity
			Real mu_;
			Real smoothing_length_;

			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;
		public:
			ComputingViscousAcceleration(FluidBody *body, StdVec<SolidBody*> interacting_bodies)
				: WeaklyCompressibleFluidDynamicsComplex(body, interacting_bodies) {
				mu_ = material_->mu_;
				smoothing_length_ = body->kernel_->GetSmoothingLength();
			};
			virtual ~ComputingViscousAcceleration() {};
		};

		/**
		 * @class TransportVelocityStress
		 * @brief  Transport velocity induced stress
		 */
		class TransportVelocityStress : public WeaklyCompressibleFluidDynamicsComplex
		{
		protected:
			virtual void InnerInteraction(size_t index_particle_i, Real dt)  override;
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)  override;

		public:
			TransportVelocityStress(FluidBody *body, StdVec<SolidBody*> interacting_bodies)
				: WeaklyCompressibleFluidDynamicsComplex(body, interacting_bodies) {};
			virtual ~TransportVelocityStress() {};
		};

		/**
		 * @class TransportVelocityCorrection
		 * @brief  transport velocty correction
		 */
		class TransportVelocityCorrection : public WeaklyCompressibleFluidDynamicsComplex
		{
		protected:
			Real p_background_;

			virtual void SetupDynamics(Real dt = 0.0) override;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0)  override;
			virtual void ContactInteraction(size_t index_particle_i,
				size_t interacting_body_index, Real dt = 0.0)  override;

		public:
			TransportVelocityCorrection(FluidBody *body, StdVec<SolidBody*> interacting_bodies)
				: WeaklyCompressibleFluidDynamicsComplex(body, interacting_bodies) {};
			virtual ~TransportVelocityCorrection() {};
		};

		/**
		 * @class TotalMechanicalEnergy
		 * @brief  Compute the total mechanical energy
		 */
		class TotalMechanicalEnergy  : public WeaklyCompressibleFluidDynamicsSum<Real>
		{
		protected:
			Real average_farctor_;
			Real potential_;
			Real ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit TotalMechanicalEnergy(FluidBody* body, ExternalForce *external_force);
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
			explicit GetAdvectionTimeStepSize(FluidBody* body, Real U_f);
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
			ComputingVorticityInFluidField(FluidBody *body) 
				: WeaklyCompressibleFluidDynamicsInner(body) {};
			virtual ~ComputingVorticityInFluidField() {};
		};

		/**
		 * @class PressureRelaxationVerletFreeSurface
		 * @brief  pressure relaxation scheme for free surface problems
		 * computing inner interaction using the Verlet time stepping
		 * in which the density updating is separated into two pecices
		 */
		class PressureRelaxationVerletFreeSurface : public WeaklyCompressibleFluidDynamicsComplex2Levels
		{
		protected:
			ExternalForce *external_force_;
			Real mu_;

			virtual void Initialization(size_t index_particle_i, Real dt = 0.0) override;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;
			virtual void Intermediate(size_t index_particle_i, Real dt = 0.0) override;
			virtual void InnerInteraction2nd(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ContactInteraction2nd(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			PressureRelaxationVerletFreeSurface(FluidBody *body, StdVec<SolidBody*> interacting_bodies, ExternalForce *external_force)
				: WeaklyCompressibleFluidDynamicsComplex2Levels(body, interacting_bodies), external_force_(external_force)
			{
				mu_ = material_->mu_;
			};
			virtual ~PressureRelaxationVerletFreeSurface() {};
		};

		/**
		 * @class PressureRelaxationVerlet
		 * @brief computing inner interaction using the Verlet time stepping
		 *  in which the density updating is separated into two pecices
		 */
		class PressureRelaxationVerlet : public PressureRelaxationVerletFreeSurface
		{
		protected:
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;

		public:
			PressureRelaxationVerlet(FluidBody *body, StdVec<SolidBody*> interacting_bodies, ExternalForce *external_force)
				: PressureRelaxationVerletFreeSurface(body, interacting_bodies, external_force) {};
			PressureRelaxationVerlet(FluidBody *body, StdVec<SolidBody*> interacting_bodies)
				: PressureRelaxationVerletFreeSurface(body, interacting_bodies, new ExternalForce) {};
			virtual ~PressureRelaxationVerlet() {};
		};

		/**
		 * @class WeaklyCompressibleFluidInitialCondition
		 * @brief  set initial condition for Oldroyd_B_Fluid dynamics
		 */
		class Oldroyd_B_FluidInitialCondition : public Oldroyd_B_FluidDynamicsSimple
		{
		protected:
			//default for set all particle at rest
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			Oldroyd_B_FluidInitialCondition(FluidBody *body)
				: Oldroyd_B_FluidDynamicsSimple(body) {};
			virtual ~Oldroyd_B_FluidInitialCondition() {};
		};

		/**
		 * @class VerletOldroyd_B_Fluid
		 * @brief computing inner interaction using the Verlet time stepping
		 * in which the density updating is separated into two pecices
		 */
		class VerletOldroyd_B_Fluid : public Oldroyd_B_FluidDynamicsComplex2Levels
		{
		protected:
			ExternalForce *external_force_;
			Real mu_, mu_p_, lambda_, p_background_;

			virtual void Initialization(size_t index_particle_i, Real dt = 0.0) override;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;
			virtual void Intermediate(size_t index_particle_i, Real dt = 0.0) override;
			virtual void InnerInteraction2nd(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ContactInteraction2nd(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;

		public:
			VerletOldroyd_B_Fluid(FluidBody *body,
				StdVec<SolidBody*> interacting_bodies, ExternalForce *external_force);
			virtual ~VerletOldroyd_B_Fluid() {};
		};

		/**
		 * @class InflowBoundaryCondition
		 * @briefinflow boundary condition
		 */
		class InflowBoundaryCondition : public WeaklyCompressibleFluidConstraint
		{
		protected:
			/** dedault value is 0.1 suggests reaching  target inflow velocity in about 10 time steps */
			Real constrain_strength_;

			/** inflow profile to be defined in applications */
			virtual Vecd GetInflowVelocity(Vecd &position, Vecd &velocity) = 0;
			virtual void ConstraintAParticle(size_t index_particle_i, Real dt = 0.0) override;
		public:
			InflowBoundaryCondition(FluidBody* body, FluidBodyPart *body_part)
				: EulerianConstraint<FluidBody, FluidParticles, FluidBodyPart>(body, body_part),
				constrain_strength_(0.1) {};
			virtual ~InflowBoundaryCondition() {};
		};
	}
}
