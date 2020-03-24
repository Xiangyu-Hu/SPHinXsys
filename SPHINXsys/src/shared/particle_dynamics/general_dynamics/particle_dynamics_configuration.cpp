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
	void ConfigurationIteratorSplitting(SplitCellLists& split_cell_lists,
		ConfigurationFunctor& configuration_functor, Real dt)
	{
		//one sweeping only
		for (size_t k = 0; k != split_cell_lists.size(); ++k) {
			StdLargeVec<CellList*>& cell_lists = split_cell_lists[k];
			for (size_t l = 0; l != cell_lists.size(); ++l)
				configuration_functor(cell_lists[l], dt);
		}
	}
	//=================================================================================================//
	void ConfigurationIteratorSplitting_parallel(SplitCellLists& split_cell_lists,
		ConfigurationFunctor& configuration_functor, Real dt)
	{
		//one sweeping only
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
	ConfigurationDynamicsSplitting::ConfigurationDynamicsSplitting(SPHBody* body)
		: ParticleDynamicsWithInnerConfigurations<SPHBody, BaseParticles>(body),
		functor_configuration_interaction_(std::bind(&ConfigurationDynamicsSplitting::ConfigurationInteraction, this, _1, _2))
	{
		mesh_cell_linked_list_ = body->mesh_cell_linked_list_;
		cell_linked_lists_ = mesh_cell_linked_list_->cell_linked_lists_;
		number_of_cells_ = mesh_cell_linked_list_->GetNumberOfCells();
		cell_spacing_ = mesh_cell_linked_list_->GetCellSpacing();
		mesh_lower_bound_ = mesh_cell_linked_list_->GetLowerBound();
		mesh_upper_bound_ = mesh_cell_linked_list_->GetUpperBound();
		cutoff_radius_ = mesh_cell_linked_list_->GetCellSpacing();
		kernel_ = body_->kernel_;
	};
	//=================================================================================================//
	void ConfigurationDynamicsSplitting::exec(Real dt)
	{
		SetupDynamics(dt);
		ConfigurationIteratorSplitting(split_cell_lists_, functor_configuration_interaction_, dt);
	}
	//=================================================================================================//
	void ConfigurationDynamicsSplitting::parallel_exec(Real dt)
	{
		SetupDynamics(dt);
		ConfigurationIteratorSplitting_parallel(split_cell_lists_, functor_configuration_interaction_, dt);
	}
	//=================================================================================================//
	void ConfigrationDynamicsWithUpdateSplitting::exec(Real dt)
	{
		ConfigurationDynamicsSplitting::exec(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void ConfigrationDynamicsWithUpdateSplitting::parallel_exec(Real dt)
	{
		ConfigurationDynamicsSplitting::parallel_exec(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void InnerConfigurationSplitting::buildNeighborList(ListData list_data,
		Neighborhood& neighborhood_here, ConcurrentListDataVector list_data_targets)
	{
		size_t particle_index_here = list_data.first;
		NeighborList& neighbor_list_here = std::get<0>(neighborhood_here);
		size_t previous_count_of_neigbors = std::get<2>(neighborhood_here);
		for (size_t n = 0; n != list_data_targets.size(); ++n)
		{
			size_t particle_index_there = list_data_targets[n].first;
			if (particle_index_there > particle_index_here)
			{
				if (particle_index_there == 0) {
					double aa = 0.0;
				}
				//displacement pointing from neighboring particle to origin particle
				Vecd displacement = list_data.second - list_data_targets[n].second;
				if (displacement.norm() <= cutoff_radius_) {
					//neigbor particles for the original particle
					std::get<1>(neighborhood_here) >= previous_count_of_neigbors ?
						neighbor_list_here.push_back(new NeighboringParticle(*kernel_, displacement, particle_index_there))
						: neighbor_list_here[std::get<1>(neighborhood_here)]->Reset(*kernel_, displacement, particle_index_there);
					std::get<1>(neighborhood_here)++;
					//neighbor particles for the target particle
					if (particle_index_there < real_particles_bound_) {
						Neighborhood& neighborhood_there = (*current_inner_configuration_)[particle_index_there];
						NeighborList& neighbor_list_there = std::get<0>(neighborhood_there);
						std::get<1>(neighborhood_there) >= std::get<2>(neighborhood_there) ?
							neighbor_list_there.push_back(new NeighboringParticle(*kernel_, displacement * (-1.0), particle_index_here))
							: neighbor_list_there[std::get<1>(neighborhood_there)]->Reset(*kernel_, displacement * (-1.0), particle_index_here);
						std::get<1>(neighborhood_there)++;
					}
				}
			}
		}
	}
	//=================================================================================================//
	void InnerConfigurationSplitting::Update(size_t index_particle_i, Real dt)
	{
		Neighborhood& neighborhood = (*current_inner_configuration_)[index_particle_i];
		size_t current_count_of_neighbors = std::get<1>(neighborhood);
		std::get<0>(neighborhood).resize(current_count_of_neighbors);
		std::get<2>(neighborhood) = current_count_of_neighbors;
		std::get<1>(neighborhood) = 0;
	}
	//=================================================================================================//
	ParticleDynamicsConfiguration::ParticleDynamicsConfiguration(SPHBody *body)
		: ParticleDynamics<void, SPHBody>(body)
	{
		update_inner_configuration_ = new InnerConfigurationSplitting(body);
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
//		update_inner_configuration_->parallel_exec();
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