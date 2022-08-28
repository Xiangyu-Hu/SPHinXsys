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

	size_t SizeOfLoopRange(const size_t &all_real_particles);
	size_t SizeOfLoopRange(const IndexVector &body_part_particles);
	size_t SizeOfLoopRange(const CellLists &body_part_cells);
	//----------------------------------------------------------------------
	//	Body-wise iterators (for sequential and parallel computing).
	//----------------------------------------------------------------------

	template <class LocalDynamics>
	void particle_for(const size_t &all_real_particles, LocalDynamics &local_dynamics,
					  void (LocalDynamics::*func_ptr)(size_t, Real), Real dt = 0.0)
	{
		for (size_t i = 0; i < all_real_particles; ++i)
			(local_dynamics.*func_ptr)(i, dt);
	};

	template <class LocalDynamics>
	void particle_parallel_for(const size_t &all_real_particles, LocalDynamics &local_dynamics,
							   void (LocalDynamics::*func_ptr)(size_t, Real), Real dt = 0.0)
	{
		parallel_for(
			IndexRange(0, all_real_particles),
			[&](const IndexRange &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					(local_dynamics.*func_ptr)(i, dt);
				}
			},
			ap);
	};

	//----------------------------------------------------------------------
	//	Body-wise reduce iterators (for sequential and parallel computing).
	//----------------------------------------------------------------------

	template <class ReturnType, typename ReduceOperation, class LocalDynamics>
	ReturnType particle_reduce(const size_t &all_real_particles, ReturnType temp, ReduceOperation &reduce_operation,
							   LocalDynamics &local_dynamics, ReturnType (LocalDynamics::*func_ptr)(size_t, Real),
							   Real dt = 0.0)
	{
		for (size_t i = 0; i < all_real_particles; ++i)
		{
			temp = reduce_operation(temp, (local_dynamics.*func_ptr)(i, dt));
		}
		return temp;
	};

	template <class ReturnType, typename ReduceOperation, class LocalDynamics>
	ReturnType particle_parallel_reduce(const size_t &all_real_particles, ReturnType temp, ReduceOperation &reduce_operation,
										LocalDynamics &local_dynamics, ReturnType (LocalDynamics::*func_ptr)(size_t, Real),
										Real dt = 0.0)
	{
		return parallel_reduce(
			IndexRange(0, all_real_particles),
			temp, [&](const IndexRange &r, ReturnType temp0) -> ReturnType
			{
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					temp0 = reduce_operation(temp0, (local_dynamics.*func_ptr)(i, dt));
				}
				return temp0; },
			[&](ReturnType x, ReturnType y) -> ReturnType
			{
				return reduce_operation(x, y);
			});
	};

	//----------------------------------------------------------------------
	//	BodypartByParticle-wise reduce iterators (for sequential and parallel computing).
	//----------------------------------------------------------------------

	template <class ReturnType, typename ReduceOperation, class LocalDynamics>
	ReturnType particle_reduce(const IndexVector &body_part_particles, ReturnType temp, ReduceOperation &reduce_operation,
							   LocalDynamics &local_dynamics, ReturnType (LocalDynamics::*func_ptr)(size_t, Real),
							   Real dt = 0.0)
	{
		for (size_t i = 0; i < body_part_particles.size(); ++i)
		{
			temp = reduce_operation(temp, (local_dynamics.*func_ptr)(body_part_particles[i], dt));
		}
		return temp;
	};

	template <class ReturnType, typename ReduceOperation, class LocalDynamics>
	ReturnType particle_parallel_reduce(const IndexVector &body_part_particles, ReturnType temp, ReduceOperation &reduce_operation,
										LocalDynamics &local_dynamics, ReturnType (LocalDynamics::*func_ptr)(size_t, Real),
										Real dt = 0.0)
	{
		return parallel_reduce(
			IndexRange(0, body_part_particles.size()),
			temp,
			[&](const IndexRange &r, ReturnType temp0) -> ReturnType
			{
				for (size_t n = r.begin(); n != r.end(); ++n)
				{
					temp0 = reduce_operation(temp0, (local_dynamics.*func_ptr)(body_part_particles[n], dt));
				}
				return temp0;
			},
			[&](ReturnType x, ReturnType y) -> ReturnType
			{
				return reduce_operation(x, y);
			});
	};

	//----------------------------------------------------------------------
	//	BodypartByCell-wise reduce iterators (for sequential and parallel computing).
	//----------------------------------------------------------------------

	template <class ReturnType, typename ReduceOperation, class LocalDynamics>
	ReturnType particle_reduce(const CellLists &body_part_cells, ReturnType temp, ReduceOperation &reduce_operation,
							   LocalDynamics &local_dynamics, ReturnType (LocalDynamics::*func_ptr)(size_t, Real),
							   Real dt = 0.0)
	{
		for (size_t i = 0; i != body_part_cells.size(); ++i)
		{
			ListDataVector &list_data = body_part_cells[i]->cell_list_data_;
			for (size_t num = 0; num < list_data.size(); ++num)
			{
				temp = reduce_operation(temp, (local_dynamics.*func_ptr)(list_data[num].first, dt));
			}
		}

		return temp;
	};

	template <class ReturnType, typename ReduceOperation, class LocalDynamics>
	ReturnType particle_parallel_reduce(const CellLists &body_part_cells, ReturnType temp, ReduceOperation &reduce_operation,
										LocalDynamics &local_dynamics, ReturnType (LocalDynamics::*func_ptr)(size_t, Real),
										Real dt = 0.0)
	{
		return parallel_reduce(
			IndexRange(0, body_part_cells.size()),
			temp,
			[&](const IndexRange &r, ReturnType temp0) -> ReturnType
			{
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					ListDataVector &list_data = body_part_cells[i]->cell_list_data_;
					for (size_t num = 0; num < list_data.size(); ++num)
					{
						temp0 = reduce_operation(temp0, (local_dynamics.*func_ptr)(list_data[num].first, dt));
					}
				}
				return temp0;
			},
			[&](ReturnType x, ReturnType y) -> ReturnType
			{ return reduce_operation(x, y); });
	};
}
#endif // PARTICLE_ITERATORS_H