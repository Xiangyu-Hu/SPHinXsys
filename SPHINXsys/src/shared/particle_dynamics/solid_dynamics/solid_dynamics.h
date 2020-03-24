/**
* @file 	solid_dynamics.h
* @brief 	Here, we define the algorithm classes for solid dynamics. 
* @details 	We consider here a weakly compresible fluids. .   
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
* @version	0.1
*/
#pragma once
#include "all_particle_dynamics.h"
#include "elastic_solid.h"
#include "weakly_compressible_fluid.h"
#include "base_kernel.h"

namespace SPH
{
	namespace solid_dynamics
	{
		typedef ParticleDynamicsSimple<SolidBody, SolidParticles, Solid> SolidDynamicsSimple;

		template <class ReturnType>
		using SolidDynamicsSum = ParticleDynamicsReduce<ReturnType, ReduceSum<ReturnType>, SolidBody, SolidParticles>;

		template <class ReturnType>
		using SolidDynamicsConstraintForSimbodySum = ConstraintByParticleReduce<ReturnType,
			ReduceSum<ReturnType>, SolidBody, SolidParticles, SolidBodyPartForSimbody>;

		template <class ReturnType>
		using ElasticSolidDynamicsConstraintForSimbodySum = ConstraintByParticleReduce<ReturnType,
			ReduceSum<ReturnType>, SolidBody, ElasticSolidParticles, SolidBodyPartForSimbody>;

		typedef ParticleDynamicsInner<SolidBody, SolidParticles, Solid> SolidDynamicsInner;

		typedef ParticleDynamicsComplex<SolidBody, SolidParticles, Solid,
			SolidBody, SolidParticles> SolidDynamicsComplex;

		typedef ParticleDynamicsComplexWithUpdate<SolidBody, SolidParticles, Solid,
			SolidBody, SolidParticles> SolidDynamicsComplexWithUpdate;

		typedef ParticleDynamicsComplexWithUpdate<SolidBody, SolidParticles, Solid,
			SolidBody, SolidParticles> SolidDynamicsComplexWithUpdate;
		
		typedef ParticleDynamicsComplex<SolidBody, SolidParticles, Solid, 
			FluidBody, FluidParticles, WeaklyCompressibleFluid> FSIDynamicsComplex;

		typedef ParticleDynamicsSimple<SolidBody, ElasticSolidParticles, ElasticSolid> ElasticSolidDynamicsSimple;

		typedef ParticleDynamicsReduce<Real, ReduceMin, SolidBody, ElasticSolidParticles, ElasticSolid> ElasticSolidDynamicsMinimum;

		typedef ParticleDynamicsInner1Level<SolidBody, ElasticSolidParticles, ElasticSolid> ElasticSolidDynamicsInner1Level;

		typedef ParticleDynamicsInnerSplitting<SolidBody, ElasticSolidParticles, ElasticSolid> ElasticSolidDynamicsInnerSplitting;

		typedef ParticleDynamicsContact<SolidBody, ElasticSolidParticles, ElasticSolid, 
			SolidBody, ElasticSolidParticles> ElasticSolidDynamicsContact;
		
		typedef ParticleDynamicsComplex<SolidBody, ElasticSolidParticles, ElasticSolid,
			SolidBody, SolidParticles> ElasticSolidDynamicsComplex;

		typedef ParticleDynamicsInner<SolidBody, ElasticSolidParticles, ElasticSolid> ElasticSolidDynamicsInner;
		typedef ParticleDynamicsComplexWithUpdate<SolidBody, ElasticSolidParticles, ElasticSolid,
			SolidBody, SolidParticles> ElasticSolidDynamicsComplexWithUpdate;

		/**
		 * @class SolidDynamicsInitialCondition
		 * @brief  set initial condition for solid fluid body
		 * This is a abstract class to be override for case specific initial conditions.
		 */
		class SolidDynamicsInitialCondition : public SolidDynamicsSimple
		{
		public:
			SolidDynamicsInitialCondition(SolidBody *body)
				: SolidDynamicsSimple(body) {};
			virtual ~SolidDynamicsInitialCondition() {};
		};
		/**
		 * @class NormalDirectionSummation
		 * @brief Computing surface normal direction of a body
		 * the value are valid for all particles but only 
		 * near the surface particle will used
		 * note that the normal is normalized.
		 */
		class NormalDirectionSummation : public SolidDynamicsComplex
		{
		protected:
			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			NormalDirectionSummation(SolidBody *body, StdVec<SolidBody*> interacting_bodies)
				: SolidDynamicsComplex(body, interacting_bodies) {};
			virtual ~NormalDirectionSummation() {};
		};

		/**
		* @class NormalDirectionReNormalization
		* @brief Computing surface normal direction of a body
		 * using a second order algorithm
		 */
		class NormalDirectionReNormalization : public SolidDynamicsComplex
		{
		protected:
			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			NormalDirectionReNormalization(SolidBody *body, StdVec<SolidBody*> interacting_bodies)
				: SolidDynamicsComplex(body, interacting_bodies) {};
			virtual ~NormalDirectionReNormalization() {};
		};

		/**
		 * @class ElasticSolidDynamicsInitialCondition
		 * @brief  set initial condition for a solid body with different material
		 * This is a abstract class to be override for case specific initial conditions.
		 */
		class ElasticSolidDynamicsInitialCondition : public ElasticSolidDynamicsSimple
		{
		public:
			ElasticSolidDynamicsInitialCondition(SolidBody *body)
				: ElasticSolidDynamicsSimple(body) {};
			virtual ~ElasticSolidDynamicsInitialCondition() {};
		};

		/**
		* @class UpdateElasticNormalDirection
		* @brief update particle normal directions for elastic solid
		*/
		class UpdateElasticNormalDirection : public ElasticSolidDynamicsSimple
		{
		protected:
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit UpdateElasticNormalDirection(SolidBody *elastic_body)
				: ElasticSolidDynamicsSimple(elastic_body) {};
			virtual ~UpdateElasticNormalDirection() {};
		};

		/**
		* @class InitializeDisplacement
		* @brief initialize the displacement for computing average velocity.
		* This class is for FSI applications to achieve smaller solid dynamics
		* time step size compared to the fluid dynamics
		*/
		class InitializeDisplacement : public ElasticSolidDynamicsSimple
		{
		protected:
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit InitializeDisplacement(SolidBody *body) : ElasticSolidDynamicsSimple(body) {};
			virtual ~InitializeDisplacement() {};
		};

		/**
		* @class UpdateAverageVelocity
		* @brief Computing average velocity.
		* This class is for FSI applications to achieve smaller solid dynamics
		* time step size compared to the fluid dynamics
		*/
		class UpdateAverageVelocity : public ElasticSolidDynamicsSimple
		{
		protected:
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit UpdateAverageVelocity(SolidBody *body) : ElasticSolidDynamicsSimple(body) {};
			virtual ~UpdateAverageVelocity() {};
		};

		/**
		* @class FluidViscousForceOnSolid
		* @brief Computing the viscous force from the fluid
		*/
		class FluidViscousForceOnSolid : public FSIDynamicsComplex
		{
		protected:
			Real mu_;
			Real smoothing_length_;

			//dynamics of a particle
			//to be realized in specific algorithms
			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override;

		public:
			FluidViscousForceOnSolid(SolidBody* body, StdVec<FluidBody*> interacting_bodies)
				: FSIDynamicsComplex(body, interacting_bodies) {
				//more work should be done for more general cases with multiple resolutions
				//and for fluids with different viscosities
				mu_ = interacting_material_[0]->getReferenceViscosity();
				/** the smmothing lenght should be discuss more. */
				smoothing_length_ = powern(2.0, body_->refinement_level_) * body->kernel_->GetSmoothingLength();
			};
			virtual ~FluidViscousForceOnSolid() {};
		};

		/**
		* @class FluidAngularConservativeViscousForceOnSolid
		* @brief Computing the viscous force from the fluid
		*/
		class FluidAngularConservativeViscousForceOnSolid : public FluidViscousForceOnSolid
		{
		protected:
			Real mu_;
			Real smoothing_length_;

			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			FluidAngularConservativeViscousForceOnSolid(SolidBody* body, StdVec<FluidBody*> interacting_bodies)
				: FluidViscousForceOnSolid(body, interacting_bodies) {};
			virtual ~FluidAngularConservativeViscousForceOnSolid() {};
		};

		/**
		* @class FluidPressureForceOnSolid
		* @brief Computing the pressure force from the fluid.
		* The pressrue force is added on the viscous force of the latter is computed. 
		* This class is for FSI applications to achieve smaller solid dynamics
		* time step size compared to the fluid dynamics
		*/
		class FluidPressureForceOnSolid : public FSIDynamicsComplex
		{
		protected:
			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			FluidPressureForceOnSolid(SolidBody *body, StdVec<FluidBody*> interacting_bodies)
				: FSIDynamicsComplex(body, interacting_bodies) {};
			virtual ~FluidPressureForceOnSolid() {};
		};

		/**
		* @class TotalViscousForceOnSolid
		* @brief Computing the total viscous force from fluid
		*/
		class TotalViscousForceOnSolid : public SolidDynamicsSum<Vecd>
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
		* @class GetAcousticTimeStepSize
		* @brief Computing the acoustic time step size
		* computing time step size
		*/
		class GetAcousticTimeStepSize : public ElasticSolidDynamicsMinimum
		{
		protected:
			Real smoothing_length_;
			Real ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit GetAcousticTimeStepSize(SolidBody* body);
			virtual ~GetAcousticTimeStepSize() {};
		};

		/**
		* @class CorrectConfiguration
		* @brief obtain the corrected initial configuration in strong form
		*/
		class CorrectConfiguration : public SolidDynamicsInner
		{
		protected:
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			CorrectConfiguration(SolidBody *body) : SolidDynamicsInner(body) {};
			virtual ~CorrectConfiguration() {};
		};

		/**
		* @class DeformationGradientTensorBySummation
		* @brief computing deformation gradient tensor by summation
		*/
		class DeformationGradientTensorBySummation : public ElasticSolidDynamicsInner
		{
		protected:
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;

		public:
			DeformationGradientTensorBySummation(SolidBody *body)
				: ElasticSolidDynamicsInner(body) {};
			virtual ~DeformationGradientTensorBySummation() {};
		};

		/**
		* @class StressRelaxationFirstHalf
		* @brief computing stress relaxation process by verlet time stepping
		* This is the first step
		*/
		class StressRelaxationFirstHalf 
			: public ParticleDynamicsInner1Level<SolidBody, ElasticSolidParticles, ElasticSolid>
		{
		protected:
			Real numerical_viscosity_;

			virtual void Initialization(size_t index_particle_i, Real dt = 0.0) override;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			StressRelaxationFirstHalf(SolidBody *body) : ParticleDynamicsInner1Level(body) {
				numerical_viscosity_ = material_->getNumericalViscosity(body_->kernel_->GetSmoothingLength());
			};
			virtual ~StressRelaxationFirstHalf() {};
			void setupDampingStressFactor(Real alpha = 1.0){numerical_viscosity_ *= alpha;}
		};

		/**
		* @class StressRelaxationSecondHalf
		* @brief computing stress relaxation process by verlet time stepping
		* This is the second step
		*/
		class StressRelaxationSecondHalf : public ElasticSolidDynamicsInner1Level
		{
		protected:
			virtual void Initialization(size_t index_particle_i, Real dt = 0.0) override;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			StressRelaxationSecondHalf(SolidBody *body) : ElasticSolidDynamicsInner1Level(body) {};
			virtual ~StressRelaxationSecondHalf() {};
		};

		/**@class ConstrainSolidBodyRegion
		 * @brief Constriant a solid body part with prescribed motion.
		 * Note the average values for FSI are prescirbed also.
		 */
		class ConstrainSolidBodyRegion 
			: public ConstraintByParticle<SolidBody, SolidParticles, BodyPartByParticle>
		{
		protected:
			virtual Vecd GetDisplacement(Vecd &pos) { return Vecd(0); };
			virtual Vecd GetVelocity(Vecd &pos) { return Vecd(0); };
			virtual Vecd GetAcceleration(Vecd &pos) { return Vecd(0); };
			virtual void ConstraintAParticle(size_t index_particle_i,
				Real dt = 0.0) override;
		public:
			ConstrainSolidBodyRegion(SolidBody *body, BodyPartByParticle*body_part)
				: ConstraintByParticle<SolidBody, SolidParticles, BodyPartByParticle>(body, body_part) {};
			virtual ~ConstrainSolidBodyRegion() {};
		};
		/**@class ConstrainSolidBodyRegionSinusoidalMotion
		 * @brief Constriant a solid body part with prescribed motion.
		 * Note the average values for FSI are prescirbed also.
		 */
		class ConstrainSolidBodyRegionSinusoidalMotion 
			: public ConstraintByParticle<SolidBody, SolidParticles, BodyPartByParticle>
		{
		protected:
			Real h_m_ = 0.1;
			Real f_ = 1.0;
			Real phi_ = 0.0;
			int id_ = 1;
			virtual Vecd GetDisplacement(Vecd &pos);
			virtual Vecd GetVelocity(Vecd &pos);
			virtual Vecd GetAcceleration(Vecd &pos);
			virtual void ConstraintAParticle(size_t index_particle_i,
				Real dt = 0.0) override;
		public:
			ConstrainSolidBodyRegionSinusoidalMotion(SolidBody *body, BodyPartByParticle*body_part,
				Real h_m, Real f, Real phi, int idex)
				: ConstraintByParticle<SolidBody, SolidParticles, BodyPartByParticle>(body, body_part) {
					h_m_ = h_m;
					f_ = f;
					phi_ = phi;
					id_ = idex;
				};
			virtual ~ConstrainSolidBodyRegionSinusoidalMotion() {};
		};
		/**@class ConstrainSolidBodyRegion
		 * @brief Constriant a solid body part with prescribed motion.
		 * Note the average values for FSI are prescirbed also.
		 */
		class constrainNormDirichletBoundary 
			: public ConstraintByParticle<SolidBody, SolidParticles, BodyPartByParticle>
		{
		protected:
			int axis_id_;
			virtual Vecd GetDisplacement(Vecd &pos) { return Vecd(0); };
			virtual Vecd GetVelocity(Vecd &pos) { return Vecd(0); };
			virtual Vecd GetAcceleration(Vecd &pos) { return Vecd(0); };
			virtual void ConstraintAParticle(size_t index_particle_i,
				Real dt = 0.0) override;
		public:
			constrainNormDirichletBoundary(SolidBody *body, BodyPartByParticle *body_part, int axis_id)
				: ConstraintByParticle<SolidBody, SolidParticles, BodyPartByParticle>(body, body_part) {axis_id_ = axis_id;};
			virtual ~constrainNormDirichletBoundary() {};
		};
		/**@class ImposeExternalForce
		 * @brief impose external force on a solid body part
		 * by add extra acceleration
		 */
		class ImposeExternalForce
			: public ConstraintByParticle<SolidBody, SolidParticles, SolidBodyPartForSimbody>
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
				: ConstraintByParticle<SolidBody, SolidParticles, SolidBodyPartForSimbody>(body, body_part) {};
			virtual ~ImposeExternalForce() {};
		};

		/**
		 * @class ConstrianSoildBodyPartBySimBody
		 * @brief Constrain a solid body part from the motion
		 * computed from Simbody.
		 */
		class ConstrianSoildBodyPartBySimBody
			: public ConstraintByParticle<SolidBody, SolidParticles, SolidBodyPartForSimbody>
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
		class ForceOnSolidBodyPartForSimBody : public SolidDynamicsConstraintForSimbodySum<SpatialVec>
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
		class ForceOnElasticBodyPartForSimBody : public ElasticSolidDynamicsConstraintForSimbodySum<SpatialVec>
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
			ForceOnElasticBodyPartForSimBody(SolidBody *body,
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
		class DampingBySplittingAlgorithm : public ElasticSolidDynamicsInnerSplitting
		{
			//viscousity
			Real eta_;

		protected:
			ElasticSolid *material_;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			DampingBySplittingAlgorithm(SolidBody *elastic_body, Real eta);
			virtual ~DampingBySplittingAlgorithm() {};
		};
	}
}