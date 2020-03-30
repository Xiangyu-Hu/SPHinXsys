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
#include "relax_body.h"
#include "relax_body_particles.h"
#include "neighbor_relation.h"
#include "base_kernel.h"
#include "mesh_cell_linked_list.h"

namespace SPH
{
	namespace relax_dynamics
	{
		typedef ParticleDynamicsSimple<RelaxBody, RelaxBodyParticles, BaseMaterial> RelaxDynamicsSimple;
		typedef ParticleDynamicsReduce<Real, ReduceMin, RelaxBody, RelaxBodyParticles> RelaxDynamicsMin;
		template <class ReturnType>
		using RelaxDynamicsSum = ParticleDynamicsReduce<ReturnType, ReduceSum<ReturnType>, RelaxBody,
			RelaxBodyParticles, BaseMaterial>;
		typedef ParticleDynamicsInner1Level<RelaxBody, RelaxBodyParticles, BaseMaterial>  RelaxDynamicsInner1Level;
		typedef ParticleDynamicsInner<RelaxBody, RelaxBodyParticles, BaseMaterial>  RelaxDynamicsInner;
		typedef ParticleDynamicsComplex1Level<RelaxBody, RelaxBodyParticles, BaseMaterial,
			RelaxBody, RelaxBodyParticles> RelaxDynamicsComplex1Level;
		typedef ConstraintByParticle<RelaxBody, RelaxBodyParticles, RelaxBodySurface> RelaxLagrangianSurfaceConstraint;
		typedef ConstraintByParticle<RelaxBody, RelaxBodyParticles, RelaxBodySingularity> RelaxLagrangianSingularityConstraint;
		typedef ParticleDynamicsComplex<RelaxBody, RelaxBodyParticles,
			BaseMaterial, RelaxBody, RelaxBodyParticles> RelaxDynamicsComplex;
		typedef ParticleDynamicsCellListSplitting<RelaxBody, RelaxBodyParticles, BaseMaterial> RelaxDynamicsSplitCellList;
	
		/**
		* @class GetTimeStepSize
		* @brief relaxation dynamics for particle initialization
		* computing time step size
		*/
		class GetTimeStepSize : public RelaxDynamicsMin
		{
		protected:
			Real smoothing_length_;
			Real ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit GetTimeStepSize(RelaxBody* body);
			virtual ~GetTimeStepSize() {};
		};

		/**
		* @class PhysicsRelaxationInner
		* @brief position verlet algorithm for physics relaxation
		* without considering contact interaction.
		* this is usually used for solid like bodies
		*/
		class PhysicsRelaxationInner : public RelaxDynamicsInner1Level
		{
		protected:
			Real smoothing_length_;
			Real sound_speed_;
			Real eta_;
			Real p0_;
			Real p_star_;
			Real back_mesh_spacing_;
			Real mass_;

			virtual void Initialization(size_t index_particle_i, Real dt = 0.0) override;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			PhysicsRelaxationInner(RelaxBody *body);
			virtual ~PhysicsRelaxationInner() {};
			void setupBackgroundPressure(Real p_b){p_star_ = p_b;}
		};
		/**
		* @class PhysicsRelaxationComplex
		* @brief position verlet algorithm for physics relaxation
		* with considering contact interaction
		* this is usually used for fluid like bodies
		*/
		class PhysicsRelaxationComplex : public RelaxDynamicsComplex1Level
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
			PhysicsRelaxationComplex(RelaxBody *body, StdVec<RelaxBody*> interacting_bodies);
			virtual ~PhysicsRelaxationComplex() {};
			void setupBackgroundPressure(Real p_b){p_star_ = p_b;}
		};

		/**
		* @brief constrain surface particles by
		* map contrained particles to geometry face and
		* r = r + phi * norm (vector distance to face)
		*/
		class ConstriantSurfaceParticles : public RelaxLagrangianSurfaceConstraint
		{
		protected:
			virtual void ConstraintAParticle(size_t index_particle_i,
				Real dt = 0.0) override;
		public:
			ConstriantSurfaceParticles(RelaxBody *body, RelaxBodySurface *body_part)
				:RelaxLagrangianSurfaceConstraint(body, body_part) {};
			virtual ~ConstriantSurfaceParticles() {};
		};
		/**
		* @brief constrain Singularity particles by
		* map contrained particles to geometry face and
		* r = r + phi * norm (vector distance to face)
		*/
		class ConstriantSingularityParticles : public RelaxLagrangianSingularityConstraint
		{
		protected:
			virtual void ConstraintAParticle(size_t index_particle_i,
				Real dt = 0.0) override;
		public:
			ConstriantSingularityParticles(RelaxBody *body, RelaxBodySingularity *body_part)
				:RelaxLagrangianSingularityConstraint(body, body_part) {};
			virtual ~ConstriantSingularityParticles() {};
		};

		/**
		* @class computeNumberDensitybySummation
		* @brief  compute the particle number density by summation.
		*/
		class computeNumberDensitybySummation : public RelaxDynamicsComplex
		{
		protected:
			Real W0_;
			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			computeNumberDensitybySummation(RelaxBody *body, StdVec<RelaxBody*> interacting_bodies)
				: RelaxDynamicsComplex(body, interacting_bodies) {
				W0_ = body->kernel_->W(Vecd(0));
			};
			virtual ~computeNumberDensitybySummation() {};
		};

		/**
		 * @class getAveragedParticleNumberDensity
		 * @brief  Compute the Everaged Particle Number Density.
		 */
		class getAveragedParticleNumberDensity  : public RelaxDynamicsSum<Real>
		{
		protected:
			Real average_farctor_;
			Real ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit getAveragedParticleNumberDensity(RelaxBody* body);
			virtual ~getAveragedParticleNumberDensity() {};
		};
		/**
		* @class updateNumberDensity
		* @brief update the number density after relaxation.
		*/
		class updateNumberDensity : public RelaxDynamicsSimple
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
			explicit updateNumberDensity(RelaxBody *body) : RelaxDynamicsSimple(body),
				sigma_(0.0){
				get_average_number_density_ = new getAveragedParticleNumberDensity(body);
			};
			virtual ~updateNumberDensity() {};
		};
	}
}