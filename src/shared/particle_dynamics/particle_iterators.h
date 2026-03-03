/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	particle_iterators.h
 * @brief 	This is for the base functions for particle iterator.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef PARTICLE_ITERATORS_H
#define PARTICLE_ITERATORS_H

#include "base_data_type_package.h"
#include "implementation.h"
#include "sphinxsys_containers.h"

namespace SPH
{
using namespace execution;

template <class ExecutionPolicy, typename DynamicsRange, class LocalDynamicsFunction>
void particle_for(const ExecutionPolicy &execution_policy, const DynamicsRange &dynamics_range,
                  const LocalDynamicsFunction &local_dynamics_function)
{
    std::cout << "\n Error: ExecutionPolicy, DynamicsRange or LocalDynamicsFunction not defined for particle_for !" << std::endl;
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    exit(1);
};

/**
 * Range-wise iterators (for sequential and parallel computing).
 */

template <class LocalDynamicsFunction>
inline void particle_for(const SequencedPolicy &seq, const IndexRange &particles_range,
                         const LocalDynamicsFunction &local_dynamics_function)
{
    for (size_t i = particles_range.begin(); i < particles_range.end(); ++i)
        local_dynamics_function(i);
};

template <class LocalDynamicsFunction>
inline void particle_for(const ParallelPolicy &par_host, const IndexRange &particles_range,
                         const LocalDynamicsFunction &local_dynamics_function)
{
    parallel_for(
        particles_range,
        [&](const IndexRange &r)
        {
            for (size_t i = r.begin(); i < r.end(); ++i)
            {
                local_dynamics_function(i);
            }
        },
        ap);
};

/**
 * Bodypart By Particle-wise iterators (for sequential and parallel computing).
 */
template <class LocalDynamicsFunction>
inline void particle_for(const SequencedPolicy &seq, const IndexVector &body_part_particles,
                         const LocalDynamicsFunction &local_dynamics_function)
{
    for (size_t i = 0; i < body_part_particles.size(); ++i)
        local_dynamics_function(body_part_particles[i]);
};

template <class LocalDynamicsFunction>
inline void particle_for(const ParallelPolicy &par_host, const IndexVector &body_part_particles,
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
/**
 * Bodypart By Cell-wise iterators (for sequential and parallel computing).
 */
template <class LocalDynamicsFunction>
inline void particle_for(const SequencedPolicy &seq, const ConcurrentCellLists &body_part_cells,
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
}

template <class LocalDynamicsFunction>
inline void particle_for(const ParallelPolicy &par_host, const ConcurrentCellLists &body_part_cells,
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
/**
 * BodypartByCell-wise iterators on cells (for sequential and parallel computing).
 */
template <class LocalDynamicsFunction>
inline void particle_for(const SequencedPolicy &seq, const DataListsInCells &body_part_cells,
                         const LocalDynamicsFunction &local_dynamics_function)
{
    for (size_t i = 0; i != body_part_cells.size(); ++i)
        local_dynamics_function(body_part_cells[i]);
};

template <class LocalDynamicsFunction>
inline void particle_for(const ParallelPolicy &par_host, const DataListsInCells &body_part_cells,
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

template <class ExecutionPolicy, typename DynamicsRange, class ReturnType,
          typename Operation, class LocalDynamicsFunction>
void particle_reduce(const ExecutionPolicy &execution_policy, const DynamicsRange &dynamics_range,
                     const LocalDynamicsFunction &local_dynamics_function)
{
    std::cout << "\n Error: ExecutionPolicy, DynamicsRange or LocalDynamicsFunction not defined for particle dynamics !" << std::endl;
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    exit(1);
};

/**
 * Body-wise reduce iterators (for sequential and parallel computing).
 */
template <class ReturnType, typename Operation, class LocalDynamicsFunction>
inline ReturnType particle_reduce(const SequencedPolicy &seq, const IndexRange &particles_range,
                                  ReturnType temp, Operation &&operation,
                                  const LocalDynamicsFunction &local_dynamics_function)
{
    for (size_t i = particles_range.begin(); i < particles_range.end(); ++i)
    {
        temp = operation(temp, local_dynamics_function(i));
    }
    return temp;
}

template <class ReturnType, typename Operation, class LocalDynamicsFunction>
inline ReturnType particle_reduce(const ParallelPolicy &par_host, const IndexRange &particles_range,
                                  ReturnType temp, Operation &&operation,
                                  const LocalDynamicsFunction &local_dynamics_function)
{
    return parallel_reduce(
        particles_range,
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
/**
 * BodypartByParticle-wise reduce iterators (for sequential and parallel computing).
 */
template <class ReturnType, typename Operation, class LocalDynamicsFunction>
inline ReturnType particle_reduce(const SequencedPolicy &seq, const IndexVector &body_part_particles,
                                  ReturnType temp, Operation &&operation,
                                  const LocalDynamicsFunction &local_dynamics_function)
{
    for (size_t i = 0; i < body_part_particles.size(); ++i)
    {
        temp = operation(temp, local_dynamics_function(body_part_particles[i]));
    }
    return temp;
}

template <class ReturnType, typename Operation, class LocalDynamicsFunction>
inline ReturnType particle_reduce(const ParallelPolicy &par_host, const IndexVector &body_part_particles,
                                  ReturnType temp, Operation &&operation,
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

/**
 * BodypartByCell-wise reduce iterators (for sequential and parallel computing).
 */
template <class ReturnType, typename Operation, class LocalDynamicsFunction>
inline ReturnType particle_reduce(const SequencedPolicy &seq, const ConcurrentCellLists &body_part_cells,
                                  ReturnType temp, Operation &&operation,
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
}

template <class ReturnType, typename Operation, class LocalDynamicsFunction>
inline ReturnType particle_reduce(const ParallelPolicy &par_host, const ConcurrentCellLists &body_part_cells,
                                  ReturnType temp, Operation &&operation,
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
}
} // namespace SPH
#endif // PARTICLE_ITERATORS_H
