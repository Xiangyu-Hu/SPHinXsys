/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	cell_linked_list.hpp
 * @brief 	Here gives the classes for managing cell linked lists. This is the basic class
 * 			for building the particle configurations.
 * @details The cell linked list saves for each body a list of particles
 * 			located within the cell.
 * @author	Chi ZHang, Yongchuan and Xiangyu Hu
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
	void CellLinkedList::searchNeighborsByParticles(
		DynamicsRange &dynamics_range, ParticleConfiguration &particle_configuration,
		GetSearchDepth &get_search_depth, GetNeighborRelation &get_neighbor_relation)
	{
		StdLargeVec<Vecd> &pos = dynamics_range.getBaseParticles().pos_;
		particle_for(execution::ParallelPolicy(), dynamics_range.LoopRange(),
					 [&](size_t index_i)
					 {
						 int search_depth = get_search_depth(index_i);
						 Vecu target_cell_index = CellIndexFromPosition(pos[index_i]);
						 int i = (int)target_cell_index[0];
						 int j = (int)target_cell_index[1];
						 int k = (int)target_cell_index[2];

						 Neighborhood &neighborhood = particle_configuration[index_i];
						 mesh_for_each(
							 Vec3i(SMAX(i - search_depth, 0), SMAX(j - search_depth, 0), SMAX(k - search_depth, 0)),
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
