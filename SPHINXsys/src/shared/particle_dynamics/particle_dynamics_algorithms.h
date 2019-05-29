/**
* @file 	particle_dynamics_algorithms.h
* @brief 	This is the class for algorithms which are using at least two particle dynamics
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "base_particle_dynamics.h"
#include "base_particle_dynamics.hpp"

namespace SPH {

	/**
	* @class ParticleDynamicsInner1Level
	* @brief This class includes an initialization, an inner interaction and an update steps
	*/
	template <class BodyType, class ParticlesType>
	class ParticleDynamicsInner1Level : 
		public ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType>
	{
		size_t number_of_particles_;
		virtual void SetupDynamics(Real dt = 0.0) {};
		virtual void Initialization(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual void Update(size_t index_particle_i, Real dt = 0.0) = 0;

	public:
		ParticleDynamicsInner1Level(BodyType* body) : ParticleDynamicsWithInnerConfigurations(body) {
			number_of_particles_ = body->number_of_particles_;
		};
		virtual ~ParticleDynamicsInner1Level() {};

		virtual void exec(Real dt = 0.0);
		virtual void parallel_exec(Real dt = 0.0);
	};

	/**
	* @class ParticleDynamicsComplex
	* @brief compex operations involving both inner and contact particle dynamics
	*/
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	class ParticleDynamicsComplex 
		: public ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
	{
	protected:
		size_t number_of_particles_;
		virtual void SetupDynamics(Real dt = 0.0) {};
		virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) = 0;

		virtual void LoopingInnerInteraction_exec(Real dt = 0.0);
		virtual void LoopingContactInteraction_exec(Real dt = 0.0);

		virtual void LoopingInnerInteraction_parallel_exec(Real dt = 0.0);
		virtual void LoopingContactInteraction_parallel_exec(Real dt = 0.0);
	public:
		ParticleDynamicsComplex(BodyType *body, StdVec<InteractingBodytype*> interacting_bodies)
			: ParticleDynamicsWithContactConfigurations(body, interacting_bodies) {
			number_of_particles_ = body->number_of_particles_;
		};
		virtual ~ParticleDynamicsComplex() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class ParticleDynamicsComplexWithUpdate
	* @brief an extra particle updating is inlcuded
	*/
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	class ParticleDynamicsComplexWithUpdate 
		: public ParticleDynamicsComplex<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
	{
	protected:
		virtual void Update(size_t index_particle_i, Real dt = 0.0) = 0;

		virtual void LoopingUpdate_exec(Real dt = 0.0);
		virtual void LoopingUpdate_parallel_exec(Real dt = 0.0);	
	public:
		ParticleDynamicsComplexWithUpdate(BodyType *body, StdVec<InteractingBodytype*> interacting_bodies)
			: ParticleDynamicsComplex(body, interacting_bodies) {};
		virtual ~ParticleDynamicsComplexWithUpdate() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class ParticleDynamicsComplex1Level
	* @brief compex operations involving one level with both inner
	* and contact particle dynamics, initilaize and update
	*/
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	class ParticleDynamicsComplex1Level 
		: public ParticleDynamicsComplexWithUpdate<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
	{
	protected:
		virtual void Initialization(size_t index_particle_i, Real dt = 0.0) = 0;

		virtual void LoopingInitialization_exec(Real dt = 0.0);
		virtual void LoopingInitialization_parallel_exec(Real dt = 0.0);
	public:
		ParticleDynamicsComplex1Level(BodyType *body, StdVec<InteractingBodytype*> interacting_bodies)
			: ParticleDynamicsComplexWithUpdate(body, interacting_bodies) {};
		virtual ~ParticleDynamicsComplex1Level() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class ParticleDynamicsComplex2Levels
	* @brief compex operations involving 2 level with both inner
	* and contact particle dynamics, initilaize and update
	*/
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	class ParticleDynamicsComplex2Levels
		: public ParticleDynamicsComplex1Level<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
	{
	protected:
		virtual void Intermediate(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual void InnerInteraction2nd(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual void ContactInteraction2nd(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) = 0;

		virtual void LoopingIntermediate_exec(Real dt = 0.0);
		virtual void LoopingIntermediate_parallel_exec(Real dt = 0.0);

		virtual void LoopingInnerInteraction2nd_exec(Real dt = 0.0);
		virtual void LoopingContactInteraction2nd_exec(Real dt = 0.0);

		virtual void LoopingInnerInteraction2nd_parallel_exec(Real dt = 0.0);
		virtual void LoopingContactInteraction2nd_parallel_exec(Real dt = 0.0);
	public:
		ParticleDynamicsComplex2Levels(BodyType *body, StdVec<InteractingBodytype*> interacting_bodies)
			: ParticleDynamicsComplex1Level(body, interacting_bodies) {};
		virtual ~ParticleDynamicsComplex2Levels() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class ParticleDynamicsComplexSplitting
	* @brief inner interaction use splitting scheme
	*/
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	class ParticleDynamicsComplexSplitting : 
		public ParticleDynamicsComplex<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
	{
	protected:
		int numer_of_lists_;

		virtual void LoopingInnerInteraction_exec(Real dt = 0.0) override;
		virtual void LoopingInnerInteraction_parallel_exec(Real dt = 0.0) override;
	public:
		explicit ParticleDynamicsComplexSplitting(BodyType *body, StdVec<InteractingBodytype*> interacting_bodies)
			: ParticleDynamicsWithInnerConfigurations(body) {
			numer_of_lists_ = powern(3, Vecd(0).size());
		};
		virtual ~ParticleDynamicsComplexSplitting() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};
}
