/**
* @file 	relax_dynamics.h
* @brief 	This is the classes of particle relaxation in order to produce body fitted
* 			initial particle distribution.   
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "math.h"

#include "all_particle_dynamics.h"
#include "neighbor_relation.h"
#include "base_kernel.h"
#include "mesh_cell_linked_list.h"

namespace SPH
{
	namespace relax_dynamics
	{
		template <class ParticlesType = BaseParticles>
		using RelaxDynamicsSimple = ParticleDynamicsSimple<SPHBody, ParticlesType>;

		template <class ParticlesType = BaseParticles>
		using RelaxDynamicsMin = ParticleDynamicsReduce<Real, ReduceMin, SPHBody, ParticlesType>;

		template <class ReturnType, class ParticlesType = BaseParticles>
		using RelaxDynamicsSum = ParticleDynamicsReduce<ReturnType, ReduceSum<ReturnType>, SPHBody,	ParticlesType>;

		template <class ParticlesType = BaseParticles>
		using RelaxDynamicsInner1Level = ParticleDynamicsInner1Level<SPHBody, ParticlesType>;
	
		template <class ParticlesType = BaseParticles>
		using RelaxDynamicsInner =  ParticleDynamicsInner<SPHBody, ParticlesType>;

		template <class ParticlesType = BaseParticles>
		using RelaxSurfaceConstraint = ConstraintByParticle<SPHBody, ParticlesType, BodySurface>;

		template <class ParticlesType = BaseParticles>
		using RelaxConstranitByCell = ConstraintByCell<SPHBody, ParticlesType, NearBodySurface>;

		template <class ParticlesType = BaseParticles, class InteractingParticlesType = BaseParticles>
		using RelaxDynamicsComplex =  
			ParticleDynamicsComplex<SPHBody, ParticlesType, BaseMaterial, SPHBody, ParticlesType>;

		template <class ParticlesType = BaseParticles, class InteractingParticlesType = BaseParticles>
		using RelaxDynamicsComplex1Level 
			= ParticleDynamicsComplex1Level<SPHBody, ParticlesType, BaseMaterial, SPHBody, ParticlesType>;
	
		/**
		* @class GetTimeStepSize
		* @brief relaxation dynamics for particle initialization
		* computing time step size
		*/
		class GetTimeStepSize : public RelaxDynamicsMin<BaseParticles>
		{
		protected:
			Real smoothing_length_;
			Real ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit GetTimeStepSize(SPHBody* body);
			virtual ~GetTimeStepSize() {};
		};

		/**
		* @class PhysicsRelaxationInner
		* @brief simple algorithm for physics relaxation
		* without considering contact interaction.
		* this is usually used for solid like bodies
		*/
		class PhysicsRelaxationInner : public ParticleDynamicsInner<SPHBody>
		{
		protected:
			Real smoothing_length_;
			Real particle_spacing_;
			Real sound_speed_;
			Real p0_;
			Real p_star_;
			Real back_mesh_spacing_;
			Real mass_;
			/** Background mesh.*/
			MeshBackground* mesh_background_;

			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			PhysicsRelaxationInner(SPHBody *body);
			virtual ~PhysicsRelaxationInner() {};
		};
		/**
		* @class PhysicsRelaxationComplex
		* @brief position verlet algorithm for physics relaxation
		* with considering contact interaction
		* this is usually used for fluid like bodies
		*/
		class PhysicsRelaxationComplex : public RelaxDynamicsComplex1Level<BaseParticles>
		{
		protected:
			Real eta_;
			Real p0_;
			Real sound_speed_;
			Real p_star_;
			Real mass_;

			virtual void Initialization(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			PhysicsRelaxationComplex(SPHBody* body, StdVec<SPHBody*> interacting_bodies);
			virtual ~PhysicsRelaxationComplex() {};
		};

		/**
		* @class BodySurfaceBounding
		* @brief constrain surface particles by
		* map contrained particles to geometry face and
		* r = r + phi * norm (vector distance to face)
		*/
		class BodySurfaceBounding : public RelaxConstranitByCell<BaseParticles>
		{
		protected:
			virtual void ConstraintAParticle(size_t index_particle_i, Real dt = 0.0) override;
		public:
			BodySurfaceBounding(SPHBody *body, NearBodySurface* body_part)
				:RelaxConstranitByCell<BaseParticles>(body, body_part) {};
			virtual ~BodySurfaceBounding() {};
		};

		/**
		* @class ConstriantSurfaceParticles
		* @brief constrain surface particles by
		* map contrained particles to geometry face and
		* r = r + phi * norm (vector distance to face)
		*/
		class ConstriantSurfaceParticles : public RelaxSurfaceConstraint<BaseParticles>
		{
		protected:
			virtual void ConstraintAParticle(size_t index_particle_i, Real dt = 0.0) override;
		public:
			ConstriantSurfaceParticles(SPHBody* body, BodySurface* body_part)
				:RelaxSurfaceConstraint<BaseParticles>(body, body_part) {};
			virtual ~ConstriantSurfaceParticles() {};
		};

		/**
		* @class computeNumberDensitybySummation
		* @brief  compute the particle number density by summation.
		*/
		class computeNumberDensitybySummation : public RelaxDynamicsComplex<BaseParticles>
		{
		protected:
			Real W0_;
			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			computeNumberDensitybySummation(SPHBody *body, StdVec<SPHBody*> interacting_bodies)
				: RelaxDynamicsComplex<BaseParticles>(body, interacting_bodies) {
				W0_ = body->kernel_->W(Vecd(0));
			};
			virtual ~computeNumberDensitybySummation() {};
		};

		/**
		 * @class getAveragedParticleNumberDensity
		 * @brief  Compute the Everaged Particle Number Density.
		 */
		class getAveragedParticleNumberDensity  : public RelaxDynamicsSum<Real, BaseParticles>
		{
		protected:
			Real average_farctor_;
			Real ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit getAveragedParticleNumberDensity(SPHBody* body);
			virtual ~getAveragedParticleNumberDensity() {};
		};
		/**
		* @class FinalizingParticleRelaxation
		* @brief update the number density after relaxation.
		*/
		class FinalizingParticleRelaxation : public RelaxDynamicsSimple<BaseParticles>
		{
		protected:
			/** the average particle number desnity. */
			Real sigma_;
			/** the method to compute average particle number desnity. */
			getAveragedParticleNumberDensity* get_average_number_density_;

			/** the function for set global parameters for the particle dynamics */
			virtual void SetupDynamics(Real dt = 0.0) override { 
				sigma_ = get_average_number_density_->parallel_exec(); 
			};
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit FinalizingParticleRelaxation(SPHBody *body)
				: RelaxDynamicsSimple<BaseParticles>(body), sigma_(0.0){
				get_average_number_density_ = new getAveragedParticleNumberDensity(body);
			};
			virtual ~FinalizingParticleRelaxation() {};
		};
	}
}