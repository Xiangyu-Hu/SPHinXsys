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
template <int lower0, int upper0,
          int lower1, int upper1, typename FunctionOnEach>
inline void mesh_for_each2d(const FunctionOnEach &function)
{
    for (int l = lower0; l != upper0; ++l)
        for (int m = lower1; m != upper1; ++m)
        {
            function(l, m);
        }
}
//=================================================================================================//
template <int lower0, int upper0,
          int lower1, int upper1, typename CheckOnEach>
inline Array2i mesh_find_if2d(const CheckOnEach &function)
{
    for (int l = lower0; l != upper0; ++l)
        for (int m = lower1; m != upper1; ++m)
        {
            if (function(l, m))
                return Array2i(l, m);
        }
    return Array2i(upper0, upper1);
}
//=================================================================================================//
template <typename FunctionOnEach>
void mesh_for_each(const Array2i &lower, const Array2i &upper, const FunctionOnEach &function)
{
    for (int l = lower[0]; l != upper[0]; ++l)
        for (int m = lower[1]; m != upper[1]; ++m)
        {
            function(Array2i(l, m));
        }
}
//=================================================================================================//
template <typename FunctionOnEach>
Array2i mesh_find_if(const Array2i &lower, const Array2i &upper, const FunctionOnEach &function)
{
    for (int l = lower[0]; l != upper[0]; ++l)
        for (int m = lower[1]; m != upper[1]; ++m)
        {
            if (function(l, m))
                return Array2i(l, m);
        }
    return upper;
}
//=================================================================================================//
template <typename LocalFunction, typename... Args>
void mesh_for(const MeshRange &mesh_range, const LocalFunction &local_function, Args &&...args)
{
    for (int i = (mesh_range.first)[0]; i != (mesh_range.second)[0]; ++i)
        for (int j = (mesh_range.first)[1]; j != (mesh_range.second)[1]; ++j)
        {
            local_function(Array2i(i, j));
        }
}
//=================================================================================================//
template <typename LocalFunction, typename... Args>
void mesh_parallel_for(const MeshRange &mesh_range, const LocalFunction &local_function, Args &&...args)
{
    parallel_for(
        IndexRange2d((mesh_range.first)[0], (mesh_range.second)[0],
                     (mesh_range.first)[1], (mesh_range.second)[1]),
        [&](const IndexRange2d &r)
        {
            for (size_t i = r.rows().begin(); i != r.rows().end(); ++i)
                for (size_t j = r.cols().begin(); j != r.cols().end(); ++j)
                {
                    local_function(Array2i(i, j));
                }
        },
        ap);
}
//=================================================================================================//
} // namespace SPH
