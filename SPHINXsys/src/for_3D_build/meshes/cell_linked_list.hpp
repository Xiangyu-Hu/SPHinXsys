/**
 * @file 	cell_linked_list.hpp
 * @brief 	Here, template functions in CellLinkedList are defined.
 * @author	Chi Zhang and Xiangyu Hu
 */

#pragma once

#include "base_particles.h"
#include "base_kernel.h"
#include "cell_linked_list.h"
#include "particle_iterators.h"
#include "mesh_iterators.hpp"

namespace SPH
{
	//=================================================================================================//
	template <class DynamicsRange, typename GetSearchDepth, typename GetNeighborRelation>
	void CellLinkedList::searchNeighborsByParticles(DynamicsRange &dynamics_range,
													ParticleConfiguration &particle_configuration,
													GetSearchDepth &get_search_depth,
													GetNeighborRelation &get_neighbor_relation)
	{
		StdLargeVec<Vecd> &pos = dynamics_range.getBaseParticles().pos_;
		particle_parallel_for(dynamics_range.LoopRange(),
							  [&](size_t index_i)
							  {
								  int search_depth = get_search_depth(index_i);
								  Vecu target_cell_index = CellIndexFromPosition(pos[index_i]);
								  int i = (int)target_cell_index[0];
								  int j = (int)target_cell_index[1];
								  int k = (int)target_cell_index[2];

								  Neighborhood &neighborhood = particle_configuration[index_i];
								  mesh_for_each(Vec3i(SMAX(i - search_depth, 0), SMAX(j - search_depth, 0), SMAX(k - search_depth, 0)),
												Vec3i(SMIN(i + search_depth + 1, (int)number_of_cells_[0]),
													  SMIN(j + search_depth + 1, (int)number_of_cells_[1]),
													  SMIN(k + search_depth + 1, (int)number_of_cells_[2])),
												[&](int l, int m, int n)
												{
													ListDataVector &target_particles = cell_data_lists_[l][m][n];
													for (const ListData &list_data : target_particles)
													{
														get_neighbor_relation(neighborhood, pos[index_i], index_i, list_data);
													}
												});
							  });
	}
	//=================================================================================================//
}
