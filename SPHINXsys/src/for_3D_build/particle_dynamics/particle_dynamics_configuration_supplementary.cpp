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
		int k = (int)cell_location[2];

		Neighborhood& neighborhood = (*inner_configuration_)[index_particle_i];
		NeighborList& neighbor_list = std::get<0>(neighborhood);
		size_t previous_count_of_neigbors = std::get<2>(neighborhood);

		for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells_[0]) - 1); ++l)
		{
			for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells_[1]) - 1); ++m)
			{
				for (int q = SMAX(k - 1, 0); q <= SMIN(k + 1, int(number_of_cells_[2]) - 1); ++q)
				{
					ConcurrentListDataVector& target_particles = cell_linked_lists_[l][m][q].particle_data_lists_;
					for (size_t n = 0; n != target_particles.size(); ++n)
					{
						//displacement pointing from neighboring particle to origin particle
						Vecd displacement = base_particle_data_i.pos_n_ - target_particles[n].second;
						if (displacement.norm() <= cell_spacing_ && index_particle_i != target_particles[n].first)
						{
							std::get<1>(neighborhood) >= previous_count_of_neigbors ?
								neighbor_list.push_back(new NeighborRelationType(base_particle_data, *kernel_,
									displacement, index_particle_i, target_particles[n].first))
								: neighbor_list[std::get<1>(neighborhood)]->resetRelation(base_particle_data, *kernel_,
									displacement, index_particle_i, target_particles[n].first);
							std::get<1>(neighborhood)++;
						}
					}
				}
			}
		}
		size_t current_count_of_neighbors = std::get<1>(neighborhood);
		neighbor_list.resize(current_count_of_neighbors);
		std::get<2>(neighborhood) = current_count_of_neighbors;
		std::get<1>(neighborhood) = 0;
	}
	//template definitions should be instantiated here
	template class ConfigurationDynamicsInner<NeighborRelation>;
	template class ConfigurationDynamicsInner<NeighborRelationWithVariableSmoothingLength>;
	//=================================================================================================//
	void ParticleSortingSplit::ConfigurationInteraction(CellList* cell_list_here, Real dt)
	{
		cout << "\n The function "
			<< "ParticleSortingSplitting::ConfigurationInteractions"
			<< " is not done in 3D. Exit the program! \n";
		exit(0);
	}
	//=================================================================================================//
}
//=================================================================================================//