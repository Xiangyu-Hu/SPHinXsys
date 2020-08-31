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
* @file 	solid_dynamics.h
* @brief 	Here, we define the algorithm classes for solid dynamics. 
* @details 	We consider here a weakly compressible solids.   
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
		using SolidDynamicsConstraintForSimbodySum = PartDynamicsByParticleReduce<ReturnType,
			ReduceSum<ReturnType>, SolidBody, SolidParticles, SolidBodyPartForSimbody>;

		template <class ReturnType>
		using ElasticSolidDynamicsConstraintForSimbodySum = PartDynamicsByParticleReduce<ReturnType,
			ReduceSum<ReturnType>, SolidBody, ElasticSolidParticles, SolidBodyPartForSimbody>;

		typedef ParticleDynamicsInner<SolidBody, SolidParticles, Solid> SolidDynamicsInner;

		typedef ParticleDynamicsComplex<SolidBody, SolidParticles, Solid,
			SolidBody, SolidParticles> SolidDynamicsComplex;

		typedef ParticleDynamicsComplexWithUpdate<SolidBody, SolidParticles, Solid,
			SolidBody, SolidParticles> SolidDynamicsComplexWithUpdate;

		typedef ParticleDynamicsComplexWithUpdate<SolidBody, SolidParticles, Solid,
			SolidBody, SolidParticles> SolidDynamicsComplexWithUpdate;
		
		typedef ParticleDynamicsContact<SolidBody, SolidParticles, Solid, 
			FluidBody, FluidParticles, WeaklyCompressibleFluid> FSIDynamics;

		typedef ParticleDynamicsSimple<SolidBody, ElasticSolidParticles, ElasticSolid> ElasticSolidDynamicsSimple;

		typedef ParticleDynamicsReduce<Real, ReduceMin, SolidBody, ElasticSolidParticles, ElasticSolid> ElasticSolidDynamicsMinimum;

		typedef ParticleDynamicsInner1Level<SolidBody, ElasticSolidParticles, ElasticSolid> ElasticSolidDynamicsInner1Level;

		typedef ParticleDynamicsInnerSplitting<SolidBody, ElasticSolidParticles, ElasticSolid> ElasticSolidDynamicsInnerSplitting;

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
			NormalDirectionSummation(SPHBodyComplexRelation* body_complex_relation)
				: SolidDynamicsComplex(body_complex_relation) {};
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
			NormalDirectionReNormalization(SPHBodyComplexRelation* body_complex_relation)
				: SolidDynamicsComplex(body_complex_relation) {};
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
		class FluidViscousForceOnSolid : public FSIDynamics
		{
		protected:
			Real mu_;
			Real smoothing_length_;

			//dynamics of a particle
			//to be realized in specific algorithms
			virtual void ContactInteraction(size_t index_particle_i, Real dt = 0.0) override;

		public:
			FluidViscousForceOnSolid(SPHBodyContactRelation* body_contact_relation)
				: FSIDynamics(body_contact_relation) {
				//more work should be done for more general cases with multiple resolutions
				//and for fluids with different viscosities
				mu_ = contact_material_[0]->getReferenceViscosity();
				/** the smoothing length should be discuss more. */
				smoothing_length_ = powern(2.0, body_->refinement_level_) * body_->kernel_->GetSmoothingLength();
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
			virtual void ContactInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			FluidAngularConservativeViscousForceOnSolid(SPHBodyContactRelation* body_contact_relation)
				: FluidViscousForceOnSolid(body_contact_relation) {};
			virtual ~FluidAngularConservativeViscousForceOnSolid() {};
		};

		/**
		* @class FluidPressureForceOnSolid
		* @brief Computing the pressure force from the fluid.
		* The pressrue force is added on the viscous force of the latter is computed. 
		* This class is for FSI applications to achieve smaller solid dynamics
		* time step size compared to the fluid dynamics
		*/
		class FluidPressureForceOnSolid : public FSIDynamics
		{
		protected:
			virtual void ContactInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			FluidPressureForceOnSolid(SPHBodyContactRelation* body_contact_relation)
				: FSIDynamics(body_contact_relation) {};
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
			CorrectConfiguration(SPHBodyInnerRelation* body_inner_relation) 
				: SolidDynamicsInner(body_inner_relation) {};
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
			DeformationGradientTensorBySummation(SPHBodyInnerRelation* body_inner_relation)
				: ElasticSolidDynamicsInner(body_inner_relation) {};
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
			StressRelaxationFirstHalf(SPHBodyInnerRelation* body_inner_relation) 
				: ParticleDynamicsInner1Level(body_inner_relation) {
				numerical_viscosity_ 
					= material_->getNumericalViscosity(body_->kernel_->GetSmoothingLength());
			};
			virtual ~StressRelaxationFirstHalf() {};
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
			StressRelaxationSecondHalf(SPHBodyInnerRelation* body_inner_relation) 
				: ElasticSolidDynamicsInner1Level(body_inner_relation) {};
			virtual ~StressRelaxationSecondHalf() {};
		};

		/**@class ConstrainSolidBodyRegion
		 * @brief Constrain a solid body part with prescribed motion.
		 * Note the average values for FSI are prescirbed also.
		 */
		class ConstrainSolidBodyRegion 
			: public PartDynamicsByParticle<SolidBody, SolidParticles, BodyPartByParticle>
		{
		protected:
			virtual Vecd GetDisplacement(Vecd &pos) { return Vecd(0); };
			virtual Vecd GetVelocity(Vecd &pos) { return Vecd(0); };
			virtual Vecd GetAcceleration(Vecd &pos) { return Vecd(0); };
			virtual void Update(size_t index_particle_i,
				Real dt = 0.0) override;
		public:
			ConstrainSolidBodyRegion(SolidBody *body, BodyPartByParticle*body_part)
				: PartDynamicsByParticle<SolidBody, 
				SolidParticles, BodyPartByParticle>(body, body_part) {};
			virtual ~ConstrainSolidBodyRegion() {};
		};
		/**@class ConstrainSolidBodyRegionSinusoidalMotion
		 * @brief Constrain a solid body part with prescribed motion.
		 * Note the average values for FSI are prescirbed also.
		 */
		class ConstrainSolidBodyRegionSinusoidalMotion 
			: public PartDynamicsByParticle<SolidBody, SolidParticles, BodyPartByParticle>
		{
		protected:
			Real h_m_ = 0.1;
			Real f_ = 1.0;
			Real phi_ = 0.0;
			int id_ = 1;
			virtual Vecd GetDisplacement(Vecd &pos);
			virtual Vecd GetVelocity(Vecd &pos);
			virtual Vecd GetAcceleration(Vecd &pos);
			virtual void Update(size_t index_particle_i,
				Real dt = 0.0) override;
		public:
			ConstrainSolidBodyRegionSinusoidalMotion(SolidBody *body, BodyPartByParticle*body_part,
				Real h_m, Real f, Real phi, int idex)
				: PartDynamicsByParticle<SolidBody, SolidParticles, BodyPartByParticle>(body, body_part) {
					h_m_ = h_m;
					f_ = f;
					phi_ = phi;
					id_ = idex;
				};
			virtual ~ConstrainSolidBodyRegionSinusoidalMotion() {};
		};
		/**@class ConstrainSolidBodyRegion
		 * @brief Constrain a solid body part with prescribed motion.
		 * Note the average values for FSI are prescirbed also.
		 */
		class constrainNormDirichletBoundary 
			: public PartDynamicsByParticle<SolidBody, SolidParticles, BodyPartByParticle>
		{
		protected:
			int axis_id_;
			virtual Vecd GetDisplacement(Vecd &pos) { return Vecd(0); };
			virtual Vecd GetVelocity(Vecd &pos) { return Vecd(0); };
			virtual Vecd GetAcceleration(Vecd &pos) { return Vecd(0); };
			virtual void Update(size_t index_particle_i,
				Real dt = 0.0) override;
		public:
			constrainNormDirichletBoundary(SolidBody *body, BodyPartByParticle *body_part, int axis_id)
				: PartDynamicsByParticle<SolidBody, SolidParticles, BodyPartByParticle>(body, body_part) {axis_id_ = axis_id;};
			virtual ~constrainNormDirichletBoundary() {};
		};
		/**@class ImposeExternalForce
		 * @brief impose external force on a solid body part
		 * by add extra acceleration
		 */
		class ImposeExternalForce
			: public PartDynamicsByParticle<SolidBody, SolidParticles, SolidBodyPartForSimbody>
		{
		protected:
			/**
			 * @brief acceleration will be specified by the application
			 */
			virtual Vecd GetAcceleration(Vecd &pos) = 0;
			virtual void Update(size_t index_particle_i,
				Real dt = 0.0) override;
		public:
			ImposeExternalForce(SolidBody *body, SolidBodyPartForSimbody *body_part)
				: PartDynamicsByParticle<SolidBody, SolidParticles, SolidBodyPartForSimbody>(body, body_part) {};
			virtual ~ImposeExternalForce() {};
		};

		/**
		 * @class ConstrainSolidBodyPartBySimBody
		 * @brief Constrain a solid body part from the motion
		 * computed from Simbody.
		 */
		class ConstrainSolidBodyPartBySimBody
			: public PartDynamicsByParticle<SolidBody, SolidParticles, SolidBodyPartForSimbody>
		{
		protected:
			SimTK::MultibodySystem &MBsystem_;
			SimTK::MobilizedBody& mobod_;
			SimTK::RungeKuttaMersonIntegrator &integ_;
			SimTK::Force::DiscreteForces &force_on_bodies_;
			const SimTK::State *simbody_state_;
			Vec3d initial_mobod_origin_location_;

			virtual void setupDynamics(Real dt = 0.0) override;
			void virtual Update(size_t index_particle_i,
				Real dt = 0.0) override;
		public:
			ConstrainSolidBodyPartBySimBody(SolidBody *body,
				SolidBodyPartForSimbody *body_part,
				SimTK::MultibodySystem &MBsystem,
				SimTK::MobilizedBody &mobod,
				SimTK::Force::DiscreteForces &force_on_bodies,
				SimTK::RungeKuttaMersonIntegrator &integ);
			virtual ~ConstrainSolidBodyPartBySimBody() {};
		};
		/**
		 * @class ConstrainNormalDirectionforSoildBodyPartBySimBody
		 * @brief Constrain normal directin for a solid body part from the motion
		 * computed from Simbody.
		 */
		class ConstrainNormalDirectionforSoildBodyPartBySimBody
			: public ConstrainSolidBodyPartBySimBody
		{
			void virtual Update(size_t index_particle_i, Real dt = 0.0) override;
			public:
			ConstrainNormalDirectionforSoildBodyPartBySimBody(SolidBody *body,
				SolidBodyPartForSimbody *body_part,
				SimTK::MultibodySystem &MBsystem,
				SimTK::MobilizedBody &mobod,
				SimTK::Force::DiscreteForces &force_on_bodies,
				SimTK::RungeKuttaMersonIntegrator &integ)
			: ConstrainSolidBodyPartBySimBody(body, body_part, MBsystem, mobod, force_on_bodies, integ){};
			virtual ~ConstrainNormalDirectionforSoildBodyPartBySimBody() {};
		};
		/**
		 * @class ForceOnSolidBodyPartForSimBody
		 * @brief Compute the force acting on the solid body part
		 * for applying to simbody forces latter
		 */
		class ForceOnSolidBodyPartForSimBody : public SolidDynamicsConstraintForSimbodySum<SimTK::SpatialVec>
		{
		protected:
			SimTK::MultibodySystem &MBsystem_;
			SimTK::MobilizedBody& mobod_;
			SimTK::RungeKuttaMersonIntegrator &integ_;
			SimTK::Force::DiscreteForces &force_on_bodies_;
			const SimTK::State *simbody_state_;
			Vec3d current_mobod_origin_location_;

			virtual void SetupReduce() override;
			virtual SimTK::SpatialVec ReduceFunction(size_t index_particle_i,	Real dt = 0.0) override;
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
		class ForceOnElasticBodyPartForSimBody : public ElasticSolidDynamicsConstraintForSimbodySum<SimTK::SpatialVec>
		{
		protected:
			SimTK::MultibodySystem &MBsystem_;
			SimTK::MobilizedBody& mobod_;
			SimTK::RungeKuttaMersonIntegrator &integ_;
			SimTK::Force::DiscreteForces &force_on_bodies_;
			const SimTK::State *simbody_state_;
			Vec3d current_mobod_origin_location_;

			virtual void SetupReduce() override;
			virtual SimTK::SpatialVec ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;
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
		protected:
			//viscosity
			Real eta_;
			virtual void setupDynamics(Real dt = 0.0) override;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			DampingBySplittingAlgorithm(SPHBodyInnerRelation* body_inner_relation);
			virtual ~DampingBySplittingAlgorithm() {};
		};


		/**
		* @class DampingBySplittingAlgorithm
		* @brief Velocity damping by splitting scheme
		* this method modify the total acceleration and velocity directly
		*/
		class DampingBySplittingPairwise : public ElasticSolidDynamicsInnerSplitting
		{
		protected:
			//viscosity
			Real eta_;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			DampingBySplittingPairwise(SPHBodyInnerRelation* body_inner_relation);
			virtual ~DampingBySplittingPairwise() {};
		};

		/**
		* @class DampingBySplittingWithRandomChoice
		* @brief Velocity damping by splitting scheme
		* this method modify the total acceleration and velocity directly
		*/
		class DampingBySplittingWithRandomChoice : public DampingBySplittingAlgorithm
		{
		protected:
			Real random_ratio_;
			bool RandomChoice();
		public:
			DampingBySplittingWithRandomChoice(SPHBodyInnerRelation* body_inner_relation, Real random_ratio);
			virtual ~DampingBySplittingWithRandomChoice() {};

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;

		};
	}
}