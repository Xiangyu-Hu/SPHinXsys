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
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	particle_iterators_ck.h
 * @brief 	This is for the base functions for particle iterator.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef PARTICLE_ITERATORS_CK_H
#define PARTICLE_ITERATORS_CK_H

#include "implementation.h"
#include "loop_range.h"

#include <numeric>

namespace SPH
{
using namespace execution;

template <class DynamicsIdentifier, class UnaryFunc>
void particle_for(const LoopRangeCK<SequencedPolicy, DynamicsIdentifier> &loop_range,
                  const UnaryFunc &unary_func)
{
    for (size_t i = 0; i < loop_range.LoopBound(); ++i)
        loop_range.computeUnit(unary_func, i);
};

template <class DynamicsIdentifier, class UnaryFunc>
void particle_for(const LoopRangeCK<ParallelPolicy, DynamicsIdentifier> &loop_range,
                  const UnaryFunc &unary_func)
{
    parallel_for(
        IndexRange(0, loop_range.LoopBound()),
        [&](const IndexRange &r)
        {
            for (size_t i = r.begin(); i < r.end(); ++i)
            {
                loop_range.computeUnit(unary_func, i);
            }
        },
        ap);
};

template <class DynamicsIdentifier, class UnaryFunc>
void particle_for(const LoopRangeCK<ParallelPolicy, DynamicsIdentifier, Splitting> &loop_range,
                  const UnaryFunc &unary_func)
{
    const Arrayi array3s = 3 * Arrayi::Ones();
    const BaseMesh mesh(loop_range.getMesh());

    // forward sweeping
    for (UnsignedInt k = 0; k < array3s.prod(); k++)
    {
        // get the corresponding 2D/3D split cell index (m, n)
        // e.g., for k = 0, split_cell_index = (0,0), for k = 3, split_cell_index = (1,0), etc.
        const Arrayi split_cell_index = mesh.transfer1DtoMeshIndex(array3s, k);
        // get the number of cells belonging to the split cell k
        // i_max = (M - m - 1) / 3 + 1, j_max = (N - n - 1) / 3 + 1
        // e.g. all_cells = (M,N) = (6, 9), (m, n) = (1, 1), then i_max = 2, j_max = 3
        const Arrayi all_cells_k = (mesh.AllCells() - split_cell_index - Arrayi::Ones()) / 3 + Arrayi::Ones();
        const UnsignedInt number_of_cells = all_cells_k.prod(); // i_max * j_max

        parallel_for(
            IndexRange(0, number_of_cells),
            [&](const IndexRange &r)
            {
                for (size_t i = r.begin(); i < r.end(); ++i)
                {
                    // get the 2D/3D cell index of the l-th cell in the split cell k
                    // (i , j) = (m + 3 * (l / j_max), n + 3 * l % i_max)
                    // e.g. all_cells = (M,N) = (6, 9), (m, n) = (1, 1), l = 0, then (i, j) = (1, 1)
                    // l = 1, then (i, j) = (1, 4), l = 3, then (i, j) = (4, 1), etc.
                    const Arrayi cell_index = split_cell_index + 3 * mesh.transfer1DtoMeshIndex(all_cells_k, l);
                    UnsignedInt linear_index = mesh_offset + mesh.LinearCellIndexFromCellIndex(cell_index);
                    loop_range.computeUnit(unary_func, linear_index);
                }
            },
            ap);
    }

    // backward sweeping
    for (UnsignedInt k = array3s.prod(); k != 0; --k)
    {
        const Arrayi split_cell_index = mesh.transfer1DtoMeshIndex(array3s, k - 1);
        const Arrayi all_cells_k = (mesh.AllCells() - split_cell_index - Arrayi::Ones()) / 3 + Arrayi::Ones();
        const UnsignedInt number_of_cells = all_cells_k.prod();

        parallel_for(
            IndexRange(0, number_of_cells),
            [&](const IndexRange &r)
            {
                for (size_t i = r.begin(); i < r.end(); ++i)
                {
                    // get the 2D/3D cell index of the l-th cell in the split cell k
                    // (i , j) = (m + 3 * (l / j_max), n + 3 * l % i_max)
                    // e.g. all_cells = (M,N) = (6, 9), (m, n) = (1, 1), l = 0, then (i, j) = (1, 1)
                    // l = 1, then (i, j) = (1, 4), l = 3, then (i, j) = (4, 1), etc.
                    const Arrayi cell_index = split_cell_index + 3 * mesh.transfer1DtoMeshIndex(all_cells_k, l);
                    UnsignedInt linear_index = mesh_offset + mesh.LinearCellIndexFromCellIndex(cell_index);
                    loop_range.computeUnit(unary_func, linear_index);
                }
            },
            ap);
    }
};

template <typename Operation, class DynamicsIdentifier, class ReturnType, class UnaryFunc>
ReturnType particle_reduce(const LoopRangeCK<SequencedPolicy, DynamicsIdentifier> &loop_range,
                           ReturnType temp, const UnaryFunc &unary_func)
{
    Operation operation;
    ReturnType temp0 = temp;
    for (size_t i = 0; i < loop_range.LoopBound(); ++i)
    {
        temp0 = operation(temp0, loop_range.computeUnit(temp, operation, unary_func, i));
    }
    return temp0;
}

template <typename Operation, class DynamicsIdentifier, class ReturnType, class UnaryFunc>
ReturnType particle_reduce(const LoopRangeCK<ParallelPolicy, DynamicsIdentifier> &loop_range,
                           ReturnType temp, const UnaryFunc &unary_func)
{
    Operation operation;
    return parallel_reduce(
        IndexRange(0, loop_range.LoopBound()), temp,
        [&](const IndexRange &r, ReturnType temp0) -> ReturnType
        {
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					temp0 = operation(temp0, loop_range.computeUnit(temp, operation, unary_func, i));
				}
				return temp0; },
        [&](const ReturnType &x, const ReturnType &y) -> ReturnType
        {
            return operation(x, y);
        });
};

template <typename T, typename Op>
T exclusive_scan(const SequencedPolicy &seq_policy, T *first, T *d_first, UnsignedInt d_size, Op op)
{
    UnsignedInt scan_size = d_size - 1;
    std::exclusive_scan(first, first + d_size, d_first, T{0}, op);
    return d_first[scan_size];
}

template <typename T, typename Op>
T exclusive_scan(const ParallelPolicy &par_policy, T *first, T *d_first, UnsignedInt d_size, Op op)
{
    // Exclusive scan is the same as inclusive, but shifted by one
    UnsignedInt scan_size = d_size - 1;
    d_first[0] = T{0};
    using range_type = tbb::blocked_range<UnsignedInt>;
    tbb::parallel_scan(
        range_type(0, scan_size), d_first[0],
        [=](const range_type &r, T sum, bool is_final_scan) -> T
        {
            T tmp = sum;
            for (UnsignedInt i = r.begin(); i < r.end(); ++i)
            {
                tmp = op(tmp, first[i]);
                if (is_final_scan)
                {
                    d_first[i + 1] = tmp;
                }
            }
            return tmp;
        },
        [&](const T &a, const T &b)
        {
            return op(a, b);
        });
    return d_first[scan_size];
}
} // namespace SPH
#endif // PARTICLE_ITERATORS_CK_H
