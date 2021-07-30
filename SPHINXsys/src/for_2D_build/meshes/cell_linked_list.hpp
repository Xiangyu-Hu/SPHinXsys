/**
 * @file 	body_relation.hpp
 * @brief 	Here, Functions belong to BaseBody, RealBody and FictitiousBody are given.
 * @author	hi ZHang and Xiangyu Hu
 */

#pragma once

#include "base_particles.h"
#include "cell_linked_list.h"

namespace SPH
{
	//=================================================================================================//
	template<typename GetParticleIndex, typename GetSearchDepth, typename GetNeighborRelation>
	void CellLinkedList::searchNeighborsByParticles(size_t total_real_particles, BaseParticles& source_particles, 
		ParticleConfiguration& particle_configuration, GetParticleIndex& get_particle_index, 
		GetSearchDepth& get_search_depth, GetNeighborRelation& get_neighbor_relation)
	{
		parallel_for(blocked_range<size_t>(0, total_real_particles),
			[&](const blocked_range<size_t>& r) {
				StdLargeVec<Vecd>& pos_n = source_particles.pos_n_;
				for (size_t num = r.begin(); num != r.end(); ++num) {
					size_t index_i = get_particle_index(num);
					Vecd& particle_position = pos_n[index_i];
					int search_depth = get_search_depth(index_i);
					Vecu target_cell_index = CellIndexFromPosition(particle_position);
					int i = (int)target_cell_index[0];
					int j = (int)target_cell_index[1];

					Neighborhood& neighborhood = particle_configuration[index_i];
					for (int l = SMAX(i - search_depth, 0); l <= SMIN(i + search_depth, int(number_of_cells_[0]) - 1); ++l)
						for (int m = SMAX(j - search_depth, 0); m <= SMIN(j + search_depth, int(number_of_cells_[1]) - 1); ++m)
						{
							ListDataVector& target_particles = cell_linked_lists_[l][m].cell_list_data_;
							for (const ListData& list_data : target_particles)
							{
								//displacement pointing from neighboring particle to origin particle
								Vecd displacement = particle_position - list_data.second;
								get_neighbor_relation(neighborhood, displacement, index_i, list_data.first);
							}
						}
				}
			}, ap);
	}
	//=================================================================================================//
}
