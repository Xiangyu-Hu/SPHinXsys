/**
 * @file 	particle_dynamics_configuration.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "particle_dynamics_configuration.h"
#include "base_body.h"

namespace SPH {
	//=================================================================================================//
	template <class NeighborRelationType>
	void ConfigurationDynamicsInner<NeighborRelationType>::InnerInteraction(size_t index_particle_i, Real dt)
	{
		StdLargeVec<BaseParticleData>& base_particle_data = particles_->base_particle_data_;
		BaseParticleData& base_particle_data_i = base_particle_data[index_particle_i];
		Vecu cell_location 
			= mesh_cell_linked_list_->GridIndexesFromPosition(base_particle_data_i.pos_n_);
		int i = (int)cell_location[0];
		int j = (int)cell_location[1];

		Neighborhood& neighborhood = (*inner_configuration_)[index_particle_i];
		NeighborList& neighbor_list = std::get<0>(neighborhood);
		size_t previous_count_of_neigbors = std::get<2>(neighborhood);

		for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells_[0]) - 1); ++l)
			for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells_[1]) - 1); ++m)
			{
				ConcurrentListDataVector& target_particles = cell_linked_lists_[l][m].particle_data_lists_;
				for (size_t n = 0; n != target_particles.size(); ++n)
				{
					//displacement pointing from neighboring particle to origin particle
					Vecd displacement = base_particle_data_i.pos_n_ - target_particles[n].second;
					if (displacement.norm() <= cell_spacing_ && index_particle_i != target_particles[n].first)
					{
						std::get<1>(neighborhood) >= neighbor_list.size() ?
							neighbor_list.emplace_back(new NeighborRelationType(base_particle_data, *kernel_,
								displacement, index_particle_i, target_particles[n].first))
							: neighbor_list[std::get<1>(neighborhood)]->resetRelation(base_particle_data , *kernel_,
								displacement, index_particle_i, target_particles[n].first);
						std::get<1>(neighborhood)++;
					}
				}
			}
		std::get<2>(neighborhood) = std::get<1>(neighborhood);
		std::get<1>(neighborhood) = 0;
	}
	//template definitions should be instantiated here
	template class ConfigurationDynamicsInner<NeighborRelation>;
	template class ConfigurationDynamicsInner<NeighborRelationWithVariableSmoothingLength>;
	//=================================================================================================//
	void ConfigurationDynamicsContact
		::ComplexInteraction(size_t index_particle_i, Real dt)
	{
		StdLargeVec<BaseParticleData>& base_particle_data = particles_->base_particle_data_;
		BaseParticleData& base_particle_data_i = base_particle_data[index_particle_i];

		/** Contact interaction. */
		for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
		{
			StdLargeVec<BaseParticleData>& target_base_particle_data
				= (*interacting_particles_[k]).base_particle_data_;
			BaseMeshCellLinkedList* mesh_cell_linked_list = target_mesh_cell_linked_lists_[k];
			Vecu target_number_of_cells = mesh_cell_linked_list->getNumberOfCells();
			int search_range
				= mesh_cell_linked_list->ComputingSearchRage(body_->refinement_level_,
					interacting_bodies_[k]->refinement_level_);
			Kernel& current_kernel = mesh_cell_linked_list->ChoosingKernel(body_->kernel_,
				interacting_bodies_[k]->kernel_);
			Real cutoff_radius = current_kernel.GetCutOffRadius();
			Vecu target_cell_location
				= mesh_cell_linked_list->GridIndexesFromPosition(base_particle_data_i.pos_n_);
			int i = (int)target_cell_location[0];
			int j = (int)target_cell_location[1];

			matrix_cell target_cell_linked_lists
				= mesh_cell_linked_list->getCellLinkedLists();
			Neighborhood& neighborhood = (*current_interacting_configuration_[k])[index_particle_i];
			NeighborList& neighbor_list = std::get<0>(neighborhood);
			size_t previous_count_of_neigbors = std::get<2>(neighborhood);

			for (int l = SMAX(i - search_range, 0); l <= SMIN(i + search_range, int(target_number_of_cells[0]) - 1); ++l)
				for (int m = SMAX(j - search_range, 0); m <= SMIN(j + search_range, int(target_number_of_cells[1]) - 1); ++m)
				{
					ConcurrentListDataVector& target_particles
						= target_cell_linked_lists[l][m].particle_data_lists_;
					for (size_t n = 0; n < target_particles.size(); n++)
					{
						BaseParticleData& base_particle_data_j = target_base_particle_data[target_particles[n].first];
						//displacement pointing from neighboring particle to origin particle
						Vecd displacement = base_particle_data_i.pos_n_ - target_particles[n].second;
						if (checkNeighbor(displacement.norm(), cutoff_radius, base_particle_data_i, base_particle_data_j))
						{
							std::get<1>(neighborhood) >= neighbor_list.size() ?
								neighbor_list.emplace_back(new NeighborRelation(base_particle_data, current_kernel,
									displacement, index_particle_i, target_particles[n].first))
								: neighbor_list[std::get<1>(neighborhood)]->resetRelation(base_particle_data,
									current_kernel, displacement, index_particle_i, target_particles[n].first);
							std::get<1>(neighborhood)++;
						}
					}
				}
			size_t current_count_of_neighbors = std::get<1>(neighborhood);
			std::get<2>(neighborhood) = current_count_of_neighbors;
			std::get<1>(neighborhood) = 0;
			if (current_count_of_neighbors != 0)
				indexes_interacting_particles_[k]->push_back(index_particle_i);
		}
	}
	//=================================================================================================//
	void ParticleSortingSplit::ConfigurationInteraction(CellList* cell_list, Real dt)
	{
		Vecu& cell_location = cell_list->cell_location_;
		int i = (int)cell_location[0];
		int j = (int)cell_location[1];

		ConcurrentListDataVector& particle_data_lists_i
			= cell_list->particle_data_lists_;
		for (size_t list_index_i = 0; list_index_i != particle_data_lists_i.size(); ++list_index_i)
		{
			size_t index_i = particle_data_lists_i[list_index_i].first;
			Real index_i_in_real = Real(index_i);
			Real sigma = W0_;
			Real index_mean = W0_ * index_i_in_real;
			for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells_[0]) - 1); ++l)
				for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells_[1]) - 1); ++m)
				{
					CellList& target_cell_list = cell_linked_lists_[l][m];
					ConcurrentListDataVector& target_particle_data_lists = target_cell_list.particle_data_lists_;
					for (size_t list_index_j = 0; list_index_j != target_particle_data_lists.size(); ++list_index_j)
					{
						//displacement pointing from neighboring particle to origin particle
						Vecd displacement = particle_data_lists_i[list_index_i].second
							- target_particle_data_lists[list_index_j].second;
						size_t index_j = target_particle_data_lists[list_index_j].first;
						if (displacement.norm() <= cell_spacing_ &&
							particle_data_lists_i[list_index_i].first != target_particle_data_lists[list_index_j].first)
						{
							Real index_j_in_real = Real(index_j);
							Real W_ij = kernel_->W(displacement);
							sigma += W_ij;
							index_mean += W_ij * index_j_in_real;
						}
					}
				}
			index_mean = index_mean / sigma;

			for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells_[0]) - 1); ++l)
				for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells_[1]) - 1); ++m)
				{
					CellList& cell_list_there = cell_linked_lists_[l][m];
					ConcurrentListDataVector& particle_data_lists_there = cell_list_there.particle_data_lists_;
					for (size_t list_index_there = 0; list_index_there != particle_data_lists_there.size(); ++list_index_there)
					{
						//displacement pointing from neighboring particle to origin particle
						Vecd displacement_there = particle_data_lists_i[list_index_i].second
							- particle_data_lists_there[list_index_there].second;
						if (displacement_there.norm() <= cell_spacing_ &&
							particle_data_lists_i[list_index_i].first != particle_data_lists_there[list_index_there].first)
						{
							Real W_ij_there = kernel_->W(displacement_there);

							size_t particle_index_here
								= particle_data_lists_i[list_index_i].first;
							Real index_here_in_real = Real(particle_index_here);

							size_t particle_index_there
								= particle_data_lists_there[list_index_there].first;
							Real index_there_in_real = Real(particle_index_there);

							Real index_mean_difference
								= (W0_ - W_ij_there) * (index_there_in_real - index_here_in_real) / sigma;
							Real new_index_mean = index_mean + index_mean_difference;
							Real index_sqr_difference
								= (W0_ - W_ij_there)
								* (index_there_in_real * index_there_in_real - index_here_in_real * index_here_in_real) / sigma
								+ index_mean * index_mean - new_index_mean * new_index_mean;

							if (index_sqr_difference < 0.0 && particles_->allowSwapping(particle_index_here, particle_index_there)) {
								particles_->swapParticles(particle_index_here, particle_index_there);
								std::swap(particle_data_lists_i[list_index_i].first,
									particle_data_lists_there[list_index_there].first);
								index_mean = new_index_mean;
							}
						}
					}
				}
	
		}
	}
	//=================================================================================================//
}
//=================================================================================================//