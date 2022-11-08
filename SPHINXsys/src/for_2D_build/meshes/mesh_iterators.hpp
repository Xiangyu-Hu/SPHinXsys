/**
 * @file 	mesh_iterators.hpp
 * @brief 	Here, Functions belong to mesh iterators.
 * @author	hi ZHang and Xiangyu Hu
 */

#pragma once

#include "mesh_iterators.h"

namespace SPH
{
    //=================================================================================================//
    template <typename LocalFunction, typename... Args>
    void mesh_for(const MeshRange &mesh_range, const LocalFunction &local_function, Args &&...args)
    {
        for (size_t i = (mesh_range.first)[0]; i != (mesh_range.second)[0]; ++i)
            for (size_t j = (mesh_range.first)[1]; j != (mesh_range.second)[1]; ++j)
            {
                local_function(i, j);
            }
    }
    //=================================================================================================//
    template <typename LocalFunction, typename... Args>
    void mesh_parallel_for(const MeshRange &mesh_range, const LocalFunction &local_function, Args &&...args)
    {
        parallel_for(
            blocked_range2d<size_t>((mesh_range.first)[0], (mesh_range.second)[0],
                                    (mesh_range.first)[1], (mesh_range.second)[1]),
            [&](const blocked_range2d<size_t> &r)
            {
                for (size_t i = r.rows().begin(); i != r.rows().end(); ++i)
                    for (size_t j = r.cols().begin(); j != r.cols().end(); ++j)
                    {
                        local_function(i, j);
                    }
            },
            ap);
    }
    //=================================================================================================//
}
