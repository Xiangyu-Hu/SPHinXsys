/* -----------------------------------------------------------------------------*
 *                               SPHinXsys                                      *
 * -----------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle    *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for       *
 * physical accurate simulation and aims to model coupled industrial dynamic    *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH      *
 * (smoothed particle hydrodynamics), a meshless computational method using     *
 * particle discretization.                                                     *
 *                                                                              *
 * SPHinXsys is partially funded by German Research Foundation                  *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,               *
 * HU1527/12-1 and HU1527/12-4.                                                 *
 *                                                                              *
 * Portions copyright (c) 2017-2022 Technical University of Munich and          *
 * the authors' affiliations.                                                   *
 *                                                                              *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may      *
 * not use this file except in compliance with the License. You may obtain a    *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.           *
 *                                                                              *
 * -----------------------------------------------------------------------------*/
/**
 * @file particle_iterators.h
 * @brief This is for the base functions for particle iterator.
 * @author  Xiangyu Hu
 */

#ifndef PARTICLE_ITERATORS_H
#define PARTICLE_ITERATORS_H

#include "base_data_package.h"
#include "sph_data_containers.h"

namespace SPH
{
	//----------------------------------------------------------------------
	//	Body-wise iterators (for sequential and parallel computing).
	//----------------------------------------------------------------------

	template <class LocalDynamicsFunction>
	void particle_for(const size_t &all_real_particles,
					  const LocalDynamicsFunction &local_dynamics_function)
	{
		for (size_t i = 0; i < all_real_particles; ++i)
			local_dynamics_function(i);
	};

	template <class LocalDynamicsFunction>
	void particle_parallel_for(const size_t &all_real_particles,
							   const LocalDynamicsFunction &local_dynamics_function)
	{
		parallel_for(
			IndexRange(0, all_real_particles),
			[&](const IndexRange &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					local_dynamics_function(i);
				}
			},
			ap);
	};

	//----------------------------------------------------------------------
	//	BodypartByParticle-wise iterators (for sequential and parallel computing).
	//----------------------------------------------------------------------

	template <class LocalDynamicsFunction>
	void particle_for(const IndexVector &body_part_particles,
					  const LocalDynamicsFunction &local_dynamics_function)
	{
		for (size_t i = 0; i < body_part_particles.size(); ++i)
			local_dynamics_function(body_part_particles[i]);
	};

	template <class LocalDynamicsFunction>
	void particle_parallel_for(const IndexVector &body_part_particles,
							   const LocalDynamicsFunction &local_dynamics_function)
	{
		parallel_for(
			IndexRange(0, body_part_particles.size()),
			[&](const IndexRange &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					local_dynamics_function(body_part_particles[i]);
				}
			},
			ap);
	};

	//----------------------------------------------------------------------
	//	BodypartByCell-wise iterators (for sequential and parallel computing).
	//----------------------------------------------------------------------

	template <class LocalDynamicsFunction>
	void particle_for(const ConcurrentIndexesInCells &body_part_cells,
					  const LocalDynamicsFunction &local_dynamics_function)
	{
		for (size_t i = 0; i != body_part_cells.size(); ++i)
		{
			ConcurrentIndexVector &particle_indexes = *body_part_cells[i];
			for (size_t num = 0; num < particle_indexes.size(); ++num)
			{
				local_dynamics_function(particle_indexes[num]);
			}
		}
	};

	template <class LocalDynamicsFunction>
	void particle_parallel_for(const ConcurrentIndexesInCells &body_part_cells,
							   const LocalDynamicsFunction &local_dynamics_function)
	{
		parallel_for(
			IndexRange(0, body_part_cells.size()),
			[&](const IndexRange &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					ConcurrentIndexVector &particle_indexes = *body_part_cells[i];
					for (size_t num = 0; num < particle_indexes.size(); ++num)
					{
						local_dynamics_function(particle_indexes[num]);
					}
				}
			},
			ap);
	};

	//----------------------------------------------------------------------
	//	BodypartByCell-wise iterators on cells (for sequential and parallel computing).
	//----------------------------------------------------------------------

	template <class LocalDynamicsFunction>
	void cell_list_for(const DataListsInCells &body_part_cells,
					   const LocalDynamicsFunction &local_dynamics_function)
	{
		for (size_t i = 0; i != body_part_cells.size(); ++i)
			local_dynamics_function(body_part_cells[i]);
	};

	template <class LocalDynamicsFunction>
	void cell_list_parallel_for(const DataListsInCells &body_part_cells,
								const LocalDynamicsFunction &local_dynamics_function)
	{
		parallel_for(
			IndexRange(0, body_part_cells.size()),
			[&](const IndexRange &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					local_dynamics_function(body_part_cells[i]);
				}
			},
			ap);
	};

	template <class LocalDynamicsFunction>
	void particle_for_split(SplitCellLists &split_cell_lists,
							const LocalDynamicsFunction &local_dynamics_function)
	{
		// forward sweeping
		for (size_t k = 0; k != split_cell_lists.size(); ++k)
		{
			ConcurrentCellLists &cell_lists = split_cell_lists[k];
			for (size_t l = 0; l != cell_lists.size(); ++l)
			{
				ConcurrentIndexVector &particle_indexes = *cell_lists[l];
				for (size_t i = 0; i != particle_indexes.size(); ++i)
				{
					local_dynamics_function(particle_indexes[i]);
				}
			}
		}

		// backward sweeping
		for (size_t k = split_cell_lists.size(); k != 0; --k)
		{
			ConcurrentCellLists &cell_lists = split_cell_lists[k - 1];
			for (size_t l = 0; l != cell_lists.size(); ++l)
			{
				ConcurrentIndexVector &particle_indexes = *cell_lists[l];
				for (size_t i = particle_indexes.size(); i != 0; --i)
				{
					local_dynamics_function(particle_indexes[i - 1]);
				}
			}
		}
	};

	//----------------------------------------------------------------------
	//	Splitting algorithm (for sequential and parallel computing).
	//----------------------------------------------------------------------
	template <class LocalDynamicsFunction>
	void particle_parallel_for_split(SplitCellLists &split_cell_lists,
									 const LocalDynamicsFunction &local_dynamics_function)
	{
		// forward sweeping
		for (size_t k = 0; k != split_cell_lists.size(); ++k)
		{
			ConcurrentCellLists &cell_lists = split_cell_lists[k];
			parallel_for(
				blocked_range<size_t>(0, cell_lists.size()),
				[&](const blocked_range<size_t> &r)
				{
					for (size_t l = r.begin(); l < r.end(); ++l)
					{
						ConcurrentIndexVector &particle_indexes = *cell_lists[l];
						for (size_t i = 0; i < particle_indexes.size(); ++i)
						{
							local_dynamics_function(particle_indexes[i]);
						}
					}
				},
				ap);
		}

		// backward sweeping
		for (size_t k = split_cell_lists.size(); k != 0; --k)
		{
			ConcurrentCellLists &cell_lists = split_cell_lists[k - 1];
			parallel_for(
				blocked_range<size_t>(0, cell_lists.size()),
				[&](const blocked_range<size_t> &r)
				{
					for (size_t l = r.begin(); l < r.end(); ++l)
					{
						ConcurrentIndexVector &particle_indexes = *cell_lists[l];
						for (size_t i = particle_indexes.size(); i != 0; --i)
						{
							local_dynamics_function(particle_indexes[i - 1]);
						}
					}
				},
				ap);
		}
	}

	//----------------------------------------------------------------------
	//	Body-wise reduce iterators (for sequential and parallel computing).
	//----------------------------------------------------------------------

	template <class ReturnType, typename Operation, class LocalDynamicsFunction>
	ReturnType particle_reduce(const size_t &all_real_particles, ReturnType temp, Operation &operation,
							   const LocalDynamicsFunction &local_dynamics_function)
	{
		for (size_t i = 0; i < all_real_particles; ++i)
		{
			temp = operation(temp, local_dynamics_function(i));
		}
		return temp;
	};

	template <class ReturnType, typename Operation, class LocalDynamicsFunction>
	ReturnType particle_parallel_reduce(const size_t &all_real_particles, ReturnType temp, Operation &operation,
										const LocalDynamicsFunction &local_dynamics_function)
	{
		return parallel_reduce(
			IndexRange(0, all_real_particles),
			temp, [&](const IndexRange &r, ReturnType temp0) -> ReturnType
			{
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					temp0 = operation(temp0, local_dynamics_function(i));
				}
				return temp0; },
			[&](const ReturnType &x, const ReturnType &y) -> ReturnType
			{
				return operation(x, y);
			});
	};

	//----------------------------------------------------------------------
	//	BodypartByParticle-wise reduce iterators (for sequential and parallel computing).
	//----------------------------------------------------------------------

	template <class ReturnType, typename Operation, class LocalDynamicsFunction>
	ReturnType particle_reduce(const IndexVector &body_part_particles, ReturnType temp, Operation &operation,
							   const LocalDynamicsFunction &local_dynamics_function)
	{
		for (size_t i = 0; i < body_part_particles.size(); ++i)
		{
			temp = operation(temp, local_dynamics_function(body_part_particles[i]));
		}
		return temp;
	};

	template <class ReturnType, typename Operation, class LocalDynamicsFunction>
	ReturnType particle_parallel_reduce(const IndexVector &body_part_particles, ReturnType temp, Operation &operation,
										const LocalDynamicsFunction &local_dynamics_function)
	{
		return parallel_reduce(
			IndexRange(0, body_part_particles.size()),
			temp,
			[&](const IndexRange &r, ReturnType temp0) -> ReturnType
			{
				for (size_t n = r.begin(); n != r.end(); ++n)
				{
					temp0 = operation(temp0, local_dynamics_function(body_part_particles[n]));
				}
				return temp0;
			},
			[&](const ReturnType &x, const ReturnType &y) -> ReturnType
			{
				return operation(x, y);
			});
	};

	//----------------------------------------------------------------------
	//	BodypartByCell-wise reduce iterators (for sequential and parallel computing).
	//----------------------------------------------------------------------

	template <class ReturnType, typename Operation, class LocalDynamicsFunction>
	ReturnType particle_reduce(const ConcurrentIndexesInCells &body_part_cells, ReturnType temp, Operation &operation,
							   const LocalDynamicsFunction &local_dynamics_function)
	{
		for (size_t i = 0; i != body_part_cells.size(); ++i)
		{
			ConcurrentIndexVector &particle_indexes = *body_part_cells[i];
			for (size_t num = 0; num < particle_indexes.size(); ++num)
			{
				temp = operation(temp, local_dynamics_function(particle_indexes[num]));
			}
		}

		return temp;
	};

	template <class ReturnType, typename Operation, class LocalDynamicsFunction>
	ReturnType particle_parallel_reduce(const ConcurrentIndexesInCells &body_part_cells, ReturnType temp, Operation &operation,
										const LocalDynamicsFunction &local_dynamics_function)
	{
		return parallel_reduce(
			IndexRange(0, body_part_cells.size()),
			temp,
			[&](const IndexRange &r, ReturnType temp0) -> ReturnType
			{
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					ConcurrentIndexVector &particle_indexes = *body_part_cells[i];
					for (size_t num = 0; num < particle_indexes.size(); ++num)
					{
						temp0 = operation(temp0, local_dynamics_function(particle_indexes[num]));
					}
				}
				return temp0;
			},
			[&](const ReturnType &x, const ReturnType &y) -> ReturnType
			{ return operation(x, y); });
	};
}
#endif // PARTICLE_ITERATORS_H