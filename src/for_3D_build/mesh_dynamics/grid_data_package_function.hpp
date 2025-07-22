#ifndef GRID_DATA_PACKAGE_FUNCTIONS_3D_HPP
#define GRID_DATA_PACKAGE_FUNCTIONS_3D_HPP

#include "grid_data_package_function.h"

namespace SPH
{
//=============================================================================================//
template <int PKG_SIZE>
inline std::pair<size_t, Array3i> NeighbourIndexShift(const Array3i shift_index, const CellNeighborhood3d &neighbour)
{
    std::pair<size_t, Array3i> result;
    Array3i neighbour_index = (shift_index + PKG_SIZE * Array3i::Ones()) / PKG_SIZE;
    result.first = neighbour[neighbour_index[0]][neighbour_index[1]][neighbour_index[2]];
    result.second = (shift_index + PKG_SIZE * Array3i::Ones()) - neighbour_index * PKG_SIZE;

    return result;
}
//=============================================================================================//
} // namespace SPH
#endif // GRID_DATA_PACKAGE_FUNCTIONS_3D_HPP