#pragma once
#include "all_particle_dynamics.h"
#include "elastic_solid.h"
#include "weakly_compressible_fluid.h"
#include "base_kernel.h"

namespace SPH
{
	namespace solid_dynamics
	{
		/**
		 * @class NormalDirectionSummation
		 * @brief Computing surface normal direction of a body
		 * the value are valid for all particles but only 
		 * near the surface particle will used
		 * note that the normal is normalized.
		 */
		class NormalDirectionSummation : 
			public ParticleDynamicsComplex<SolidBody, SolidBodyParticles, SolidBody, SolidBodyParticles>
		{
		protected:
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;

		public:
			NormalDirectionSummation(SolidBody *body, StdVec<SolidBody*> interacting_bodies)
				: ParticleDynamicsComplex(body, interacting_bodies) {};
			virtual ~NormalDirectionSummation() {};
		};

		/**
		* @class NormalDirectionReNormalization
		* @brief Computing surface normal direction of a body
		 * using a second order algorithm
		 */
		class NormalDirectionReNormalization : 
			public ParticleDynamicsComplexWithUpdate<ElasticBody, ElasticBodyParticles, SolidBody, SolidBodyParticles>
		{
		protected:
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ContactInteraction(size_t index_particle_i, size_t contact_body_index, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			NormalDirectionReNormalization(ElasticBody *body, StdVec<SolidBody*> interacting_bodies)
				: ParticleDynamicsComplexWithUpdate(body, interacting_bodies) {};
			virtual ~NormalDirectionReNormalization() {};
		};

		/**
		* @class UpdateElasticNormalDirection
		* @brief update particle normal directions for elastic solid
		*/
		class UpdateElasticNormalDirection : public ParticleDynamicsSimple<ElasticBody, ElasticBodyParticles>
		{
		protected:
			virtual void ParticleUpdate(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit UpdateElasticNormalDirection(ElasticBody *elastic_body)
				: ParticleDynamicsSimple<ElasticBody, ElasticBodyParticles>(elastic_body) {};
			virtual ~UpdateElasticNormalDirection() {};
		};

		/**
		* @class InitializeDisplacement
		* @brief initialize the displacement for computing average velocity.
		* This class is for FSI applications to achieve smaller solid dynamics
		* time step size compared to the fluid dynamics
		*/
		class InitializeDisplacement : public ParticleDynamicsSimple<ElasticBody, ElasticBodyParticles>
		{
		protected:
			virtual void ParticleUpdate(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit InitializeDisplacement(ElasticBody *body) : ParticleDynamicsSimple(body) {};
			virtual ~InitializeDisplacement() {};
		};

		/**
		* @class UpdateAverageVelocity
		* @brief Computing average velocity.
		* This class is for FSI applications to achieve smaller solid dynamics
		* time step size compared to the fluid dynamics
		*/
		class UpdateAverageVelocity : public ParticleDynamicsSimple<ElasticBody, ElasticBodyParticles>
		{
		protected:
			virtual void ParticleUpdate(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit UpdateAverageVelocity(ElasticBody *body) : ParticleDynamicsSimple(body) {};
			virtual ~UpdateAverageVelocity() {};
		};

		/**
		* @class FluidPressureForceOnSolid
		* @brief Computing the pressure force from the fluid
		* This class is for FSI applications to achieve smaller solid dynamics
		* time step size compared to the fluid dynamics
		*/
		class FluidPressureForceOnSolid 
			: public ParticleDynamicsComplex<SolidBody, SolidBodyParticles, WeaklyCompressibleFluidBody, WeaklyCompressibleFluidParticles>
		{
		protected:
			ExternalForce *external_force_;
			WeaklyCompressibleFluid* material_;

			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;

		public:
			FluidPressureForceOnSolid(SolidBody *body, StdVec<WeaklyCompressibleFluidBody*> interacting_bodies,
				WeaklyCompressibleFluid* material, ExternalForce *external_force)
				: ParticleDynamicsComplex(body, interacting_bodies), material_(material), external_force_(external_force) {};
			FluidPressureForceOnSolid(SolidBody *body, StdVec<WeaklyCompressibleFluidBody*> interacting_bodies,
				WeaklyCompressibleFluid* material) : ParticleDynamicsComplex(body, interacting_bodies), 
				material_(material), external_force_(new ExternalForce) {};
			virtual ~FluidPressureForceOnSolid() {};
		};

		/**
		* @class FluidPressureForceOnSolidFreeSurface
		* @brief Computing the pressure force from the free surface flow
		* This class is for FSI applications to achieve smaller solid dynamics
		* time step size compared to the fluid dynamics
		*/
		class FluidPressureForceOnSolidFreeSurface : public FluidPressureForceOnSolid
		{
		protected:
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;

		public:
			FluidPressureForceOnSolidFreeSurface(SolidBody *body, StdVec<WeaklyCompressibleFluidBody*> interacting_bodies,
				WeaklyCompressibleFluid* material, ExternalForce *external_force)
				: FluidPressureForceOnSolid(body, interacting_bodies, material, external_force) {};
			FluidPressureForceOnSolidFreeSurface(SolidBody *body, StdVec<WeaklyCompressibleFluidBody*> interacting_bodies,
				WeaklyCompressibleFluid* material)
				: FluidPressureForceOnSolid(body, interacting_bodies, material) {};
			virtual ~FluidPressureForceOnSolidFreeSurface() {};
		};

		/**
		* @class TotalViscousForceOnSolid
		* @brief Computing the total viscous force from fluid
		*/
		class TotalViscousForceOnSolid : 
			public ParticleDynamicsSum<Vecd, SolidBody, SolidBodyParticles>
		{
		protected:
			Vecd ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit TotalViscousForceOnSolid(SolidBody* body);
			virtual ~TotalViscousForceOnSolid() {};
		};

		/**
		 * @class TotalForceOnSolid
		 * @brief Computing total force from fluid.
		 */
		class TotalForceOnSolid : public TotalViscousForceOnSolid
		{
		protected:
			Vecd ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit TotalForceOnSolid(SolidBody* body) : TotalViscousForceOnSolid(body) {};
			virtual ~TotalForceOnSolid() {};
		};

		/**
		* @class FluidViscousForceOnSolid
		* @brief Computing the viscous force from the fluid
		*/
		class FluidViscousForceOnSolid : 
			public ParticleDynamicsComplex<SolidBody, SolidBodyParticles, WeaklyCompressibleFluidBody, WeaklyCompressibleFluidParticles>
		{
		protected:
			Real mu_;
			Real smoothing_length_;

			//dynamics of a particle
			//to be realized in specific algorithms
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;

		public:
			FluidViscousForceOnSolid(SolidBody *body, StdVec<WeaklyCompressibleFluidBody*> interacting_bodies,
				WeaklyCompressibleFluid* fluid_material)
				: ParticleDynamicsComplex(body, interacting_bodies) {
				mu_ = fluid_material->mu_;
				/** the smmothing lenght should be discuss more. */
				smoothing_length_ = powern(2.0, body->refinement_level_)*body->kernel_->GetSmoothingLength();
			};
			virtual ~FluidViscousForceOnSolid() {};
		};

		/**
		* @class GetAcousticTimeStepSize
		* @brief Computing the acoustic time step size
		* computing time step size
		*/
		class GetAcousticTimeStepSize : public ParticleDynamicsMinimum<ElasticBody, ElasticBodyParticles>
		{
		protected:
			Real smoothing_length_;
			Real ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit GetAcousticTimeStepSize(ElasticBody* body);
			virtual ~GetAcousticTimeStepSize() {};
		};

		/**
		* @class CorrectConfiguration
		* @brief obtain the corrected initial configuration in strong form
		*/
		class CorrectConfiguration 
			: public ParticleDynamicsComplexWithUpdate<ElasticBody, ElasticBodyParticles, SolidBody, SolidBodyParticles>
		{
		protected:
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			CorrectConfiguration(ElasticBody *body, StdVec<SolidBody*> interacting_bodies)
				: ParticleDynamicsComplexWithUpdate(body, interacting_bodies) {};
			virtual ~CorrectConfiguration() {};
		};

		/**
		* @class DeformationGradientTensorBySummation
		* @brief computing deformation gradient tensor by summation
		*/
		class DeformationGradientTensorBySummation 
			: public ParticleDynamicsComplex<ElasticBody, ElasticBodyParticles, SolidBody, SolidBodyParticles>
		{
		protected:
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;

		public:
			DeformationGradientTensorBySummation(ElasticBody *body, StdVec<SolidBody*> interacting_bodies)
				: ParticleDynamicsComplex(body, interacting_bodies) {};
			virtual ~DeformationGradientTensorBySummation() {};
		};

		/**
		* @class StressRelaxation
		* @brief computing stress relaxation process by verlet time stepping
		*/
		class StressRelaxation 
			: public ParticleDynamicsComplex2Levels<ElasticBody, ElasticBodyParticles, ElasticBody, ElasticBodyParticles>
		{
		protected:
			Real smoothing_length_;
			ElasticSolid *material_;

			virtual void Initialization(size_t index_particle_i, Real dt = 0.0) override;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;
			virtual void Intermediate(size_t index_particle_i, Real dt = 0.0) override;
			virtual void InnerInteraction2nd(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ContactInteraction2nd(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			StressRelaxation(ElasticBody *body, StdVec<ElasticBody*> interacting_bodies)
				: ParticleDynamicsComplex2Levels(body, interacting_bodies) {
				material_ = body->material_;
				smoothing_length_ = body->kernel_->GetSmoothingLength();
			};
			virtual ~StressRelaxation() {};
		};

		/**
		* @class StressRelaxationFirstStep
		* @brief computing stress relaxation process by verlet time stepping
		* This is the first step
		*/
		class StressRelaxationFirstStep : public ParticleDynamicsInner1Level<ElasticBody, ElasticBodyParticles>
		{
		protected:
			ElasticSolid *material_;

			virtual void Initialization(size_t index_particle_i, Real dt = 0.0) override;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			StressRelaxationFirstStep(ElasticBody *body) : ParticleDynamicsInner1Level(body) {
				material_ = body->material_;
			};
			virtual ~StressRelaxationFirstStep() {};
		};

		/**
		* @class StressRelaxationSecondStep
		* @brief computing stress relaxation process by verlet time stepping
		* This is the second step
		*/
		class StressRelaxationSecondStep : public ParticleDynamicsInner1Level<ElasticBody, ElasticBodyParticles>
		{
		protected:
			virtual void Initialization(size_t index_particle_i, Real dt = 0.0) override;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			StressRelaxationSecondStep(ElasticBody *body) : ParticleDynamicsInner1Level(body) {};
			virtual ~StressRelaxationSecondStep() {};
		};

		/**
		* @class StressInConstrinedElasticBodyFirstHalf
		* @brief computing stress within constrined elastic body
		* This is the first half step
		*/
		class StressInConstrinedElasticBodyFirstHalf 
			: public ParticleDynamicsSimple<ElasticBody, ElasticBodyParticles>
		{
		protected:
			ElasticSolid *material_;
			virtual void ParticleUpdate(size_t index_particle_i, Real dt = 0.0) override;

		public:
			StressInConstrinedElasticBodyFirstHalf(ElasticBody *body)
				: ParticleDynamicsSimple(body) { material_ = body->material_; };
			virtual ~StressInConstrinedElasticBodyFirstHalf() {};
		};

		/**
		* @class StressInConstrinedElasticBodySecondHalf
		* @brief computing stress within constrined elastic body
		* This is the second half step
		*/
		//computing stress within constrined elastic body
		class StressInConstrinedElasticBodySecondHalf 
			: public ParticleDynamicsContact<ElasticBody, ElasticBodyParticles, ElasticBody, ElasticBodyParticles>
		{
		protected:
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;

		public:
			StressInConstrinedElasticBodySecondHalf(ElasticBody *body, StdVec<ElasticBody*> interacting_bodies)
				: ParticleDynamicsContact(body, interacting_bodies) {};
			virtual ~StressInConstrinedElasticBodySecondHalf() {};
		};

		/**@class ConstrainSolidBodyRegion
		 * @brief Constriant a solid body part with prescribed motion.
		 * Note the average values for FSI are prescirbed also.
		 */
		class ConstrainSolidBodyRegion 
			: public LagrangianConstraint<SolidBody, SolidBodyParticles, SolidBodyPart>
		{
		protected:
			virtual Vecd GetDisplacement(Vecd &pos) { return Vecd(0); };
			virtual Vecd GetVelocity(Vecd &pos) { return Vecd(0); };
			virtual Vecd GetAcceleration(Vecd &pos) { return Vecd(0); };
			virtual void ConstraintAParticle(size_t index_particle_i,
				Real dt = 0.0) override;
		public:
			ConstrainSolidBodyRegion(SolidBody *body, SolidBodyPart *body_part)
				: LagrangianConstraint<SolidBody, SolidBodyParticles, SolidBodyPart>(body, body_part) {};
			virtual ~ConstrainSolidBodyRegion() {};
		};

		/**@class ImposeExternalForce
		 * @brief impose external force on a solid body part
		 * by add extra acceleration
		 */
		class ImposeExternalForce
			: public LagrangianConstraint<SolidBody, SolidBodyParticles, SolidBodyPartForSimbody>
		{
		protected:
			/**
			 * @brief acceleration will be specified by the application
			 */
			virtual Vecd GetAcceleration(Vecd &pos) = 0;
			virtual void ConstraintAParticle(size_t index_particle_i,
				Real dt = 0.0) override;
		public:
			ImposeExternalForce(SolidBody *body, SolidBodyPartForSimbody *body_part)
				: LagrangianConstraint<SolidBody, SolidBodyParticles, SolidBodyPartForSimbody>(body, body_part) {};
			virtual ~ImposeExternalForce() {};
		};

		/**
		 * @class ConstrianSoildBodyPartBySimBody
		 * @brief Constrain a solid body part from the motion
		 * computed from Simbody.
		 */
		class ConstrianSoildBodyPartBySimBody
			: public LagrangianConstraint<SolidBody, SolidBodyParticles, SolidBodyPartForSimbody>
		{
		protected:
			SimTK::MultibodySystem &MBsystem_;
			SimTK::MobilizedBody& mobod_;
			SimTK::RungeKuttaMersonIntegrator &integ_;
			SimTK::Force::DiscreteForces &force_on_bodies_;
			const SimTK::State *simbody_state_;
			Vec3 initial_mobod_origin_location_;

			virtual void  PrepareConstraint() override;
			void virtual ConstraintAParticle(size_t index_particle_i,
				Real dt = 0.0) override;
		public:
			ConstrianSoildBodyPartBySimBody(SolidBody *body,
				SolidBodyPartForSimbody *body_part,
				SimTK::MultibodySystem &MBsystem,
				SimTK::MobilizedBody &mobod,
				SimTK::Force::DiscreteForces &force_on_bodies,
				SimTK::RungeKuttaMersonIntegrator &integ);
			virtual ~ConstrianSoildBodyPartBySimBody() {};
		};

		/**
		 * @class ForceOnSolidBodyPartForSimBody
		 * @brief Compute the force acting on the solid body part
		 * for applying to simbody forces latter
		 */
		class ForceOnSolidBodyPartForSimBody 
			: public LagrangianConstraintSum<SpatialVec, SolidBody, SolidBodyParticles, SolidBodyPartForSimbody>
		{
		protected:
			SimTK::MultibodySystem &MBsystem_;
			SimTK::MobilizedBody& mobod_;
			SimTK::RungeKuttaMersonIntegrator &integ_;
			SimTK::Force::DiscreteForces &force_on_bodies_;
			const SimTK::State *simbody_state_;
			Vec3 current_mobod_origin_location_;

			virtual void SetupReduce() override;
			virtual SpatialVec ReduceFunction(size_t index_particle_i,	Real dt = 0.0) override;
		public:
			ForceOnSolidBodyPartForSimBody(SolidBody *body,
				SolidBodyPartForSimbody *body_part,
				SimTK::MultibodySystem &MBsystem,
				SimTK::MobilizedBody &mobod,
				SimTK::Force::DiscreteForces &force_on_bodies,
				SimTK::RungeKuttaMersonIntegrator &integ);
			virtual ~ForceOnSolidBodyPartForSimBody() {};
		};

		/**
		 * @class ForceOnElasticBodyPartForSimBody
		 * @brief Compute the force acting on the elastic solid body part
		 * for applying to simbody forces latter
		 */
		class ForceOnElasticBodyPartForSimBody 
			: public LagrangianConstraintSum<SpatialVec, ElasticBody, ElasticBodyParticles, SolidBodyPartForSimbody>
		{
		protected:
			SimTK::MultibodySystem &MBsystem_;
			SimTK::MobilizedBody& mobod_;
			SimTK::RungeKuttaMersonIntegrator &integ_;
			SimTK::Force::DiscreteForces &force_on_bodies_;
			const SimTK::State *simbody_state_;
			Vec3 current_mobod_origin_location_;

			virtual void SetupReduce() override;
			virtual SpatialVec ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			ForceOnElasticBodyPartForSimBody(ElasticBody *body,
				SolidBodyPartForSimbody *body_part,
				SimTK::MultibodySystem &MBsystem,
				SimTK::MobilizedBody &mobod,
				SimTK::Force::DiscreteForces &force_on_bodies,
				SimTK::RungeKuttaMersonIntegrator &integ);
			virtual ~ForceOnElasticBodyPartForSimBody() {};
		};

		/**
		* @class DampingBySplittingAlgorithm
		* @brief Velocity damping by splitting scheme
		* this method modify the total acceleration and velocity directly
		*/
		class DampingBySplittingAlgorithm
			: public ParticleDynamicsInnerSplitting<ElasticBody, ElasticBodyParticles>
		{
			//viscousity
			Real eta_;

		protected:
			ElasticSolid *material_;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			DampingBySplittingAlgorithm(ElasticBody *elastic_body, Real eta);
			virtual ~DampingBySplittingAlgorithm() {};
		};
	}
}