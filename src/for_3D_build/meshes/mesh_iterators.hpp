/**
 * @file 	mesh_iterators.hpp
 * @brief 	Here, Functions belong to mesh iterators.
 * @author	hi Zhang and Xiangyu Hu
 */

#pragma once

#include "mesh_iterators.h"

namespace SPH
{
//=================================================================================================//
template <typename FunctionOnEach>
void mesh_for_each(const Array3i &lower, const Array3i &upper, const FunctionOnEach &function)
{
    for (int l = lower[0]; l != upper[0]; ++l)
        for (int m = lower[1]; m != upper[1]; ++m)
            for (int n = lower[2]; n != upper[2]; ++n)
            {
                function(Array3i(l, m, n));
            }
}
//=================================================================================================//
template <typename FunctionOnEach>
void mesh_for_column_major(const Array3i &lower, const Array3i &upper, const FunctionOnEach &function)
{
    for (int l = lower[2]; l != upper[2]; ++l)
        for (int m = lower[1]; m != upper[1]; ++m)
            for (int n = lower[0]; n != upper[0]; ++n)
            {
                function(Array3i(n, m, l));
            }
}
//=================================================================================================//
template <typename FunctionOnEach>
Array3i mesh_find_if(const Array3i &lower, const Array3i &upper, const FunctionOnEach &function)
{
    for (int l = lower[0]; l != upper[0]; ++l)
        for (int m = lower[1]; m != upper[1]; ++m)
            for (int n = lower[2]; n != upper[2]; ++n)
            {
                if (function(Array3i(l, m, n)))
                    return Array3i(l, m, n);
            }
    return upper;
}
//=================================================================================================//
template <typename LocalFunction, typename... Args>
void mesh_for(const MeshRange &mesh_range, const LocalFunction &local_function, Args &&...args)
{
    for (int i = (mesh_range.first)[0]; i != (mesh_range.second)[0]; ++i)
        for (int j = (mesh_range.first)[1]; j != (mesh_range.second)[1]; ++j)
            for (int k = (mesh_range.first)[2]; k != (mesh_range.second)[2]; ++k)
            {
                local_function(Array3i(i, j, k));
            }
}
//=================================================================================================//
template <typename LocalFunction, typename... Args>
void mesh_parallel_for(const MeshRange &mesh_range, const LocalFunction &local_function, Args &&...args)
{
    parallel_for(
        IndexRange3d((mesh_range.first)[0], (mesh_range.second)[0],
                     (mesh_range.first)[1], (mesh_range.second)[1],
                     (mesh_range.first)[2], (mesh_range.second)[2]),
        [&](const IndexRange3d &r)
        {
            for (size_t i = r.pages().begin(); i != r.pages().end(); ++i)
                for (size_t j = r.rows().begin(); j != r.rows().end(); ++j)
                    for (size_t k = r.cols().begin(); k != r.cols().end(); ++k)
                    {
                        local_function(Array3i(i, j, k));
                    }
        },
        ap);
}
//=================================================================================================//
} // namespace SPH
