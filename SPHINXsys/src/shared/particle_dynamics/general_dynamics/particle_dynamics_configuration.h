/**
 * @file 	particle_dynamics_configuration.h
 * @brief 	This the class for update configurations
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#pragma once

#include "base_particle_dynamics.h"
#include "base_particle_dynamics.hpp"

namespace SPH {

	class SPHSystem;

	/**
	 * @class ParticleDynamicsConfiguration
	 * @brief Update both inner and contact configurations
	 */
	class ParticleDynamicsConfiguration : public ParticleDynamics<void, SPHBody, Particles>
	{
	protected:

	public:
		ParticleDynamicsConfiguration(SPHBody *body);
		virtual ~ParticleDynamicsConfiguration() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	};

	/**
	 * @class ParticleDynamicsConfiguration
	 * @brief Update both contact configurations only
	 */
	class ParticleDynamicsContactConfiguration : public ParticleDynamics<void, SPHBody, Particles>
	{
	protected:

	public:
		ParticleDynamicsContactConfiguration(SPHBody *body);
		virtual ~ParticleDynamicsContactConfiguration() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	};

	/**
	 * @class ParticleDynamicsInteractionConfiguration
	 * @brief Update interaction configuration only
	 */
	class ParticleDynamicsInteractionConfiguration : public ParticleDynamics<void, SPHBody, Particles>
	{
	protected:
		
		SPHBodyVector interacting_bodies_;

	public:
		ParticleDynamicsInteractionConfiguration(SPHBody *body,
			SPHBodyVector interacting_bodies);
		virtual ~ParticleDynamicsInteractionConfiguration() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	};

	/**
	 * @class ParticleDynamicsInnerConfiguration
	 * @brief Update inner configuration only
	 */
	 //update particle inner configuration only
	class ParticleDynamicsInnerConfiguration : public ParticleDynamics<void, SPHBody, Particles>
	{
	protected:

	public:
		ParticleDynamicsInnerConfiguration(SPHBody *body);
		virtual ~ParticleDynamicsInnerConfiguration() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	};

	/**
	 * @class ParticleDynamicsCellLinkedList
	 * @brief Update cell linked list
	 */
	class ParticleDynamicsCellLinkedList : public ParticleDynamics<void, SPHBody, Particles>
	{
	protected:

	public:
		ParticleDynamicsCellLinkedList(SPHBody *body);
		virtual ~ParticleDynamicsCellLinkedList() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	};

	/**
	 * @class ParticleDynamicsByCellParticleLists
	 * @brief Update by cell particle list
	 */
	class ParticleDynamicsByCellParticleLists : public ParticleDynamics<void, SPHBody, Particles>
	{
		SPHSystem *system_;

	protected:

	public:
		ParticleDynamicsByCellParticleLists(SPHSystem *system, SPHBody *body);
		virtual ~ParticleDynamicsByCellParticleLists() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	};
}
