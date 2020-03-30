/**
 * @file 	particle_dynamics_configuration.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "particle_dynamics_configuration.h"
#include "sph_system.h"
#include "base_body.h"

namespace SPH {
	//=================================================================================================//
	void ConfigurationIteratorSplit(SplitCellLists& split_cell_lists,
		ConfigurationFunctor& configuration_functor, Real dt)
	{
		for (size_t k = 0; k != split_cell_lists.size(); ++k) {
			StdLargeVec<CellList*>& cell_lists = split_cell_lists[k];
			for (size_t l = 0; l != cell_lists.size(); ++l)
				configuration_functor(cell_lists[l], dt);
		}
	}
	//=================================================================================================//
	void ConfigurationIteratorSplit_parallel(SplitCellLists& split_cell_lists,
		ConfigurationFunctor& configuration_functor, Real dt)
	{
		for (size_t k = 0; k != split_cell_lists.size(); ++k) {
			StdLargeVec<CellList*>& cell_lists = split_cell_lists[k];
			parallel_for(blocked_range<size_t>(0, cell_lists.size()),
				[&](const blocked_range<size_t>& r) {
					for (size_t l = r.begin(); l < r.end(); ++l) 
						configuration_functor(cell_lists[l], dt);
				}, ap);
		}
	}
	//=================================================================================================//
	ConfigurationDynamicsSplit::ConfigurationDynamicsSplit(SPHBody* body)
		: ParticleDynamicsWithInnerConfigurations<SPHBody, BaseParticles>(body),
		functor_configuration_interaction_(std::bind(&ConfigurationDynamicsSplit::ConfigurationInteraction, this, _1, _2))
	{
		mesh_cell_linked_list_ = body->base_mesh_cell_linked_list_;
		cell_linked_lists_ = mesh_cell_linked_list_->getCellLinkedLists();
		number_of_cells_ = mesh_cell_linked_list_->getNumberOfCells();
		cell_spacing_ = mesh_cell_linked_list_->getCellSpacing();
		kernel_ = body_->kernel_;
	}
	//=================================================================================================//
	void ConfigurationDynamicsSplit::exec(Real dt)
	{
		SetupDynamics(dt);
		ConfigurationIteratorSplit(split_cell_lists_, functor_configuration_interaction_, dt);
	}
	//=================================================================================================//
	void ConfigurationDynamicsSplit::parallel_exec(Real dt)
	{
		SetupDynamics(dt);
		ConfigurationIteratorSplit_parallel(split_cell_lists_, functor_configuration_interaction_, dt);
	}
	//=================================================================================================//
	void ConfigrationDynamicsWithUpdateSplit::exec(Real dt)
	{
		ConfigurationDynamicsSplit::exec(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void ConfigrationDynamicsWithUpdateSplit::parallel_exec(Real dt)
	{
		ConfigurationDynamicsSplit::parallel_exec(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	ParticleDynamicsConfiguration::ParticleDynamicsConfiguration(SPHBody *body)
		: ParticleDynamics<void, SPHBody>(body)
	{
	}
//=================================================================================================//
	void ParticleDynamicsConfiguration::exec(Real dt)
	{
		body_->UpdateInnerConfiguration();
		body_->UpdateContactConfiguration();
	}
//=================================================================================================//
	void ParticleDynamicsConfiguration::parallel_exec(Real dt)
	{
		body_->UpdateInnerConfiguration();
		body_->UpdateContactConfiguration();
	}
//=================================================================================================//
	ParticleDynamicsContactConfiguration
		::ParticleDynamicsContactConfiguration(SPHBody *body)
		: ParticleDynamics<void, SPHBody>(body)
	{

	}
//=================================================================================================//
	void ParticleDynamicsContactConfiguration::exec(Real dt)
	{
		body_->UpdateContactConfiguration();
	}
//=================================================================================================//
	void ParticleDynamicsContactConfiguration::parallel_exec(Real dt)
	{
		body_->UpdateContactConfiguration();
	}
//=================================================================================================//
	ParticleDynamicsInteractionConfiguration
		::ParticleDynamicsInteractionConfiguration(SPHBody *body, 
			SPHBodyVector interacting_bodies)
		: ParticleDynamics<void, SPHBody>(body), interacting_bodies_(interacting_bodies)
	{
		//building interacting topology
		interacting_bodies_ = interacting_bodies;
	}
//=================================================================================================//
	void ParticleDynamicsInteractionConfiguration::exec(Real dt)
	{
		body_->UpdateInteractionConfiguration(interacting_bodies_);
	}
//=================================================================================================//
	void ParticleDynamicsInteractionConfiguration::parallel_exec(Real dt)
	{
		body_->UpdateInteractionConfiguration(interacting_bodies_);
	}
//=================================================================================================//
	ParticleDynamicsInnerConfiguration
		::ParticleDynamicsInnerConfiguration(SPHBody *body)
		: ParticleDynamics<void, SPHBody>(body)
	{

	}
//=================================================================================================//
	void ParticleDynamicsInnerConfiguration::exec(Real dt)
	{
		body_->UpdateInnerConfiguration();
	}
//=================================================================================================//
	void ParticleDynamicsInnerConfiguration::parallel_exec(Real dt)
	{
		body_->UpdateInnerConfiguration();
	}
//=================================================================================================//
	ParticleDynamicsCellLinkedList
		::ParticleDynamicsCellLinkedList(SPHBody *body)
		: ParticleDynamics<void, SPHBody>(body)
	{

	}
//=================================================================================================//
	void ParticleDynamicsCellLinkedList::exec(Real dt)
	{
		body_->UpdateCellLinkedList();
	}
//=================================================================================================//
	void ParticleDynamicsCellLinkedList::parallel_exec(Real dt)
	{
		body_->UpdateCellLinkedList();
	}
//=================================================================================================//
	ParticleDynamicsByCellParticleLists
		::ParticleDynamicsByCellParticleLists(SPHSystem *system, SPHBody *body)
		: ParticleDynamics<void, SPHBody>(body), system_(system)
	{

	}
//=================================================================================================//
	void ParticleDynamicsByCellParticleLists::exec(Real dt)
	{
		cout << "\n ParticleDynamicsByCellParticleLists not done yet. Exit the program! \n";
		exit(0);
	}
//=================================================================================================//
	void ParticleDynamicsByCellParticleLists::parallel_exec(Real dt)
	{
		cout << "\n ParticleDynamicsByCellParticleLists not done yet. Exit the program! \n";
		exit(0);
	}
//=================================================================================================//
}
//=================================================================================================//