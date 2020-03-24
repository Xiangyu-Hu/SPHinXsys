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
	void ParticleSortingSplitting::ConfigurationInteraction(CellList* cell_list, Real dt)
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
						if (displacement.norm() <= cutoff_radius_ &&
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
						if (displacement_there.norm() <= cutoff_radius_ &&
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
	void InnerConfigurationSplitting::ConfigurationInteraction(CellList* cell_list, Real dt)
	{
		int i = (int)cell_list->cell_location_[0];
		int j = (int)cell_list->cell_location_[1];

		ConcurrentListDataVector& particle_data_lists
			= cell_list->particle_data_lists_;
		StdVec<pair<size_t, Vecd>> indexes(0);
		for (size_t num = 0; num != cell_list->real_particles_in_cell_; ++num)
		{
			ListData& list_data = particle_data_lists[num];
			size_t particle_index_here = list_data.first;
			if (list_data.first == 0) {
				double aa = 0.0;
			}
			Neighborhood& neighborhood = (*current_inner_configuration_)[particle_index_here];
			NeighborList& neighbor_list = std::get<0>(neighborhood);
			size_t previous_count_of_neigbors = std::get<2>(neighborhood);

			for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells_[0]) - 1); ++l)
				for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells_[1]) - 1); ++m)
				{
					buildNeighborList(list_data, neighborhood, cell_linked_lists_[l][m].particle_data_lists_);
				}
		}
	}
	//=================================================================================================//
}
//=================================================================================================//