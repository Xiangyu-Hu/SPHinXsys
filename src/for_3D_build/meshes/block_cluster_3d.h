#include "mesh_iterators.h"
#include <unordered_set>

#pragma once

namespace SPH
{
using Coord3D = std::array<int, 3>;
using Block3D = std::pair<Coord3D, Coord3D>; // {min_corner}, {max_corner}

struct Coord3DHash
{
    std::size_t operator()(const Coord3D &coord) const
    {
        return std::hash<int>()(coord[0]) ^ (std::hash<int>()(coord[1]) << 1) ^ (std::hash<int>()(coord[2]) << 2);
    }
};

std::vector<Block3D> clusterActiveCells3D(const std::vector<Coord3D> &activeCells);
} // namespace SPH
