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
 * @file mesh_iterators.h
 * @brief This is for the base functions for mesh iterator.
 * There are two types of functions: one is for static ranged
 * which are defined by template parameters,
 * the other for dynamics ranges which are input parameters.
 * @author  Xiangyu Hu
 */

#ifndef MESH_ITERATORS_H
#define MESH_ITERATORS_H

#include "base_data_type_package.h"

#include "execution_policy.h"

namespace SPH
{
template <typename FunctionOnEach>
void mesh_for_each(const Arrayi &lower, const Arrayi &upper, const FunctionOnEach &function);
template <typename FunctionOnEach>
void mesh_for_column_major(const Arrayi &lower, const Arrayi &upper, const FunctionOnEach &function);
template <typename FunctionOnEach>
Arrayi mesh_find_if(const Arrayi &lower, const Arrayi &upper, const FunctionOnEach &function);
template <typename FunctionOnEach>
bool mesh_any_of(const Arrayi &lower, const Arrayi &upper, const FunctionOnEach &function)
{
    return mesh_find_if(lower, upper, function).matrix() != upper.matrix();
};

using MeshRange = std::pair<Arrayi, Arrayi>;
/** Iterator on the mesh by looping index. sequential computing. */
template <typename LocalFunction, typename... Args>
void mesh_for(const MeshRange &mesh_range, const LocalFunction &local_function, Args &&...args);
/** Iterator on the mesh by looping index. parallel computing. */
template <typename LocalFunction, typename... Args>
void mesh_parallel_for(const MeshRange &mesh_range, const LocalFunction &local_function, Args &&...args);

/** Iterator on the mesh by looping index. parallel computing. */
template <typename LocalFunction, typename... Args>
void mesh_for(const execution::SequencedPolicy &seq, const MeshRange &mesh_range,
              const LocalFunction &local_function, Args &&...args)
{
    mesh_for(mesh_range, local_function, std::forward<Args>(args)...);
};

template <typename LocalFunction, typename... Args>
void mesh_for(const execution::ParallelPolicy &par_host, const MeshRange &mesh_range,
              const LocalFunction &local_function, Args &&...args)
{
    mesh_parallel_for(mesh_range, local_function, std::forward<Args>(args)...);
};

template <typename FunctionOnData>
void package_for(const execution::SequencedPolicy &seq, UnsignedInt start_index,
                 UnsignedInt end_index, const FunctionOnData &function)
{
    for (size_t i = start_index; i != end_index; ++i)
        function(i);
}

template <typename FunctionOnData>
void package_for(const execution::ParallelPolicy &par_host, UnsignedInt start_index,
                 UnsignedInt end_index, const FunctionOnData &function)
{
    parallel_for(IndexRange(start_index, end_index), [&](const IndexRange &r)
                 {
                    for (size_t i = r.begin(); i != r.end(); ++i)
                    {
                        function(i);
                    } }, ap);
}

template <typename FunctionOnData>
void package_for(const execution::ParallelDevicePolicy &par_device,
                 UnsignedInt start_index, UnsignedInt end_index,
                 const FunctionOnData &function);
} // namespace SPH
#endif // MESH_ITERATORS_H