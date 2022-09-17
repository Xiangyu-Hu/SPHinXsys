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
                for (size_t k = (mesh_range.first)[2]; k != (mesh_range.second)[2]; ++k)
                {
                    local_function(i, j, k);
                }
    }
    //=================================================================================================//
    template <typename LocalFunction, typename... Args>
    void mesh_parallel_for(const MeshRange &mesh_range, const LocalFunction &local_function, Args &&...args)
    {
        parallel_for(
            blocked_range3d<size_t>((mesh_range.first)[0], (mesh_range.second)[0],
                                    (mesh_range.first)[1], (mesh_range.second)[1],
                                    (mesh_range.first)[2], (mesh_range.second)[2]),
            [&](const blocked_range3d<size_t> &r)
            {
                for (size_t i = r.pages().begin(); i != r.pages().end(); ++i)
                    for (size_t j = r.rows().begin(); j != r.rows().end(); ++j)
                        for (size_t k = r.cols().begin(); k != r.cols().end(); ++k)
                        {
                            local_function(i, j, k);
                        }
            },
            ap);
    }
    //=================================================================================================//
}
