/**
 * @file 	particle_dynamics_configuration.h
 * @brief 	This the class for update configurations
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#pragma once

#include "base_particle_dynamics.h"
#include "base_particle_dynamics.hpp"
#include "base_kernel.h"

namespace SPH {

	class SPHSystem;

	/** Functor for cofiguration operation. */
	typedef std::function<void(CellList*, Real)> ConfigurationFunctor;
	/** Iterators for inner functors with splitting for configuration dynamics. sequential computing. */
	void ConfigurationIteratorSplitting(SplitCellLists& split_cell_lists,
		ConfigurationFunctor& configuration_functor, Real dt = 0.0);
	/** Iterators for inner functors with splitting for configuration dynamics. parallel computing. */
	void ConfigurationIteratorSplitting_parallel(SplitCellLists& split_cell_lists,
		ConfigurationFunctor& configuration_functor, Real dt = 0.0);

	/**
	 * @class ParticleDynamicsConfigrationSplitting
	 * @brief This is for using splitting algorihm to update particle configuration
	 */
	class ConfigurationDynamicsSplitting
		: public ParticleDynamicsWithInnerConfigurations<SPHBody, BaseParticles>
	{
	protected:
		MeshCellLinkedList* mesh_cell_linked_list_;
		matrix_cell cell_linked_lists_;
		Vecu number_of_cells_;
		Kernel* kernel_;
		Real cutoff_radius_;
		Real cell_spacing_;
		Vecd mesh_lower_bound_, mesh_upper_bound_;

		virtual void ConfigurationInteraction(CellList* cell_list, Real dt = 0.0) = 0;
		ConfigurationFunctor functor_configuration_interaction_;
	public:
		explicit ConfigurationDynamicsSplitting(SPHBody* body);
		virtual ~ConfigurationDynamicsSplitting() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	 * @class ParticleDynamicsConfigrationWithUpdateSplitting
	 * @brief This is for using splitting algorihm to update particle configuration
	 */
	class ConfigrationDynamicsWithUpdateSplitting
		: public ConfigurationDynamicsSplitting
	{
	protected:
		virtual void Update(size_t index_particle_i, Real dt = 0.0) = 0;
		InnerFunctor functor_update_;
	public:
		explicit ConfigrationDynamicsWithUpdateSplitting(SPHBody* body) 
			: ConfigurationDynamicsSplitting(body),
			functor_update_(std::bind(&ConfigrationDynamicsWithUpdateSplitting::Update, this, _1, _2)) {};
		virtual ~ConfigrationDynamicsWithUpdateSplitting() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
 * @class ParticleSortingSplitting
 * @brief This is for using splitting algorihm to update particle configuration
 */
	class ParticleSortingSplitting
		: public ConfigurationDynamicsSplitting
	{
	protected:
		Real W0_;
		virtual void ConfigurationInteraction(CellList* cell_list, Real dt = 0.0) override;
	public:
		ParticleSortingSplitting(SPHBody* body)
			: ConfigurationDynamicsSplitting(body), W0_(kernel_->W(Vecd(0))) {};
		virtual ~ParticleSortingSplitting() {};
	};

	/**
	 * @class InnerConfigurationSplitting
	 * @brief This is for using splitting algorihm to update particle configuration
	 */
	class InnerConfigurationSplitting
		: public ConfigrationDynamicsWithUpdateSplitting
	{
	protected:
		size_t real_particles_bound_;
		void buildNeighborList(ListData list_data, 
			Neighborhood& neighborhood_here, ConcurrentListDataVector list_data_vector);
		virtual void SetupDynamics(Real dt = 0.0) override {
			real_particles_bound_ = particles_->real_particles_bound_;
		};
		virtual void ConfigurationInteraction(CellList* cell_list, Real dt = 0.0) override;
		virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
	public:
		InnerConfigurationSplitting(SPHBody* body)
			: ConfigrationDynamicsWithUpdateSplitting(body),
			real_particles_bound_(body->base_particles_->real_particles_bound_){};
		virtual ~InnerConfigurationSplitting() {};
	};

	/**
	 * @class ParticleDynamicsConfiguration
	 * @brief Update both inner and contact configurations
	 */
	class ParticleDynamicsConfiguration : public ParticleDynamics<void, SPHBody>
	{
	protected:
		InnerConfigurationSplitting* update_inner_configuration_;
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
	class ParticleDynamicsContactConfiguration : public ParticleDynamics<void, SPHBody>
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
	class ParticleDynamicsInteractionConfiguration : public ParticleDynamics<void, SPHBody>
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
	class ParticleDynamicsInnerConfiguration : public ParticleDynamics<void, SPHBody>
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
	class ParticleDynamicsCellLinkedList : public ParticleDynamics<void, SPHBody>
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
	class ParticleDynamicsByCellParticleLists : public ParticleDynamics<void, SPHBody>
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
