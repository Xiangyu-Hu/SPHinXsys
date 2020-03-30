/**
 * @file 	particle_dynamics_configuration.h
 * @brief 	This the class for update configurations
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#pragma once

#include "particle_dynamics_algorithms.h"
#include "particle_dynamics_algorithms.hpp"
#include "base_kernel.h"

namespace SPH {

	class SPHSystem;

	/** Functor for cofiguration operation. */
	typedef std::function<void(CellList*, Real)> ConfigurationFunctor;
	/** Iterators for inner functors with splitting for configuration dynamics. sequential computing. */
	void ConfigurationIteratorSplit(SplitCellLists& split_cell_lists,
		ConfigurationFunctor& configuration_functor, Real dt = 0.0);
	/** Iterators for inner functors with splitting for configuration dynamics. parallel computing. */
	void ConfigurationIteratorSplit_parallel(SplitCellLists& split_cell_lists,
		ConfigurationFunctor& configuration_functor, Real dt = 0.0);
	/**
	 * @class ConfigurationDynamicsInner
	 * @brief This is for using splitting algorihm to update particle configuration
	 */
	template<class NeighborRelationType = NeighborRelation>
	class ConfigurationDynamicsInner
		: public ParticleDynamicsInner<SPHBody, BaseParticles>
	{
	protected:
		BaseMeshCellLinkedList* mesh_cell_linked_list_;
		matrix_cell cell_linked_lists_;
		Vecu number_of_cells_;
		Real cell_spacing_;
		Kernel* kernel_;

		virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
	public:
		explicit ConfigurationDynamicsInner(SPHBody* body, BaseMeshCellLinkedList* mesh_cell_linked_list)
			: ParticleDynamicsInner<SPHBody, BaseParticles>(body)
		{
			mesh_cell_linked_list_ = mesh_cell_linked_list;
			cell_linked_lists_ = mesh_cell_linked_list->getCellLinkedLists();
			number_of_cells_ = mesh_cell_linked_list->getNumberOfCells();
			cell_spacing_ = mesh_cell_linked_list->getCellSpacing();
			kernel_ = body_->kernel_;
		};
		virtual ~ConfigurationDynamicsInner() {};
	};

	/**
	 * @class ConfigurationDynamicsSplit
	 * @brief This is for using splitting algorihm to update particle configuration
	 */
	class ConfigurationDynamicsSplit
		: public ParticleDynamicsWithInnerConfigurations<SPHBody, BaseParticles>
	{
	protected:
		BaseMeshCellLinkedList* mesh_cell_linked_list_;
		matrix_cell cell_linked_lists_;
		Vecu number_of_cells_;
		Kernel* kernel_;
		Real cell_spacing_;

		virtual void ConfigurationInteraction(CellList* cell_list, Real dt = 0.0) = 0;
		ConfigurationFunctor functor_configuration_interaction_;
	public:
		explicit ConfigurationDynamicsSplit(SPHBody* body);
		virtual ~ConfigurationDynamicsSplit() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	 * @class ConfigrationDynamicsWithUpdateSplit
	 * @brief This is for using splitting algorihm to update particle configuration
	 */
	class ConfigrationDynamicsWithUpdateSplit
		: public ConfigurationDynamicsSplit
	{
	protected:
		virtual void Update(size_t index_particle_i, Real dt = 0.0) = 0;
		InnerFunctor functor_update_;
	public:
		explicit ConfigrationDynamicsWithUpdateSplit(SPHBody* body) 
			: ConfigurationDynamicsSplit(body),
			functor_update_(std::bind(&ConfigrationDynamicsWithUpdateSplit::Update, this, _1, _2)) {};
		virtual ~ConfigrationDynamicsWithUpdateSplit() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
 * @class ParticleSortingSplit
 * @brief This is for using splitting algorihm to update particle configuration
 */
	class ParticleSortingSplit
		: public ConfigurationDynamicsSplit
	{
	protected:
		Real W0_;
		virtual void ConfigurationInteraction(CellList* cell_list, Real dt = 0.0) override;
	public:
		ParticleSortingSplit(SPHBody* body)
			: ConfigurationDynamicsSplit(body), W0_(kernel_->W(Vecd(0))) {};
		virtual ~ParticleSortingSplit() {};
	};

	/**
	 * @class ParticleDynamicsConfiguration
	 * @brief Update both inner and contact configurations
	 */
	class ParticleDynamicsConfiguration : public ParticleDynamics<void, SPHBody>
	{
	public:
		ParticleDynamicsConfiguration(SPHBody *body);
		virtual ~ParticleDynamicsConfiguration() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	};

	/**
	 * @class ParticleDynamicsContactConfiguration
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
