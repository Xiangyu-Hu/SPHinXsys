/**
* @file 	relax_dynamics.h
* @brief 	This is the classes of particle relaxation in order to produce body fitted
* initial particle distribution.   
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "all_particle_dynamics.h"

namespace SPH
{
	namespace relax_dynamics
	{
		typedef ParticleDynamicsSimple<RelaxBody, RelaxBodyParticles, Material> RelaxDynamicsSimple;
		typedef ParticleDynamicsReduce<Real, ReduceMin, RelaxBody, RelaxBodyParticles> RelaxDynamicsMin;
		typedef ParticleDynamicsInner1Level<RelaxBody, RelaxBodyParticles, Material>  RelaxDynamicsInner1Level;
		typedef ParticleDynamicsComplex1Level<RelaxBody, RelaxBodyParticles, Material,
			RelaxBody, RelaxBodyParticles> RelaxDynamicsComplex1Level;
		typedef ConstraintByParticle<RelaxBody, RelaxBodyParticles, RelaxBodySurface> RelaxLagrangianConstraint;
		
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
			Real p_star_;
			Real back_mesh_spacing_;

			virtual void Initialization(size_t index_particle_i, Real dt = 0.0) override;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			PhysicsRelaxationInner(RelaxBody *body);
			virtual ~PhysicsRelaxationInner() {};
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
			Real eta_, p_star_;

			virtual void Initialization(size_t index_particle_i, Real dt = 0.0) override;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;

		public:
			PhysicsRelaxationComplex(RelaxBody *body, StdVec<RelaxBody*> interacting_bodies);
			virtual ~PhysicsRelaxationComplex() {};
		};

		/**
		* @brief constrain surface particles by
		* map contrained particles to geometry face and
		* r = r + phi * norm (vector distance to face)
		*/
		class ConstriantSurfaceParticles : public RelaxLagrangianConstraint
		{
		protected:
			virtual void ConstraintAParticle(size_t index_particle_i,
				Real dt = 0.0) override;
		public:
			ConstriantSurfaceParticles(RelaxBody *body, RelaxBodySurface *body_part)
				:RelaxLagrangianConstraint(body, body_part) {};
			virtual ~ConstriantSurfaceParticles() {};
		};
	}
}