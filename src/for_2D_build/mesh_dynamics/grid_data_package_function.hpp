#ifndef GRID_DATA_PACKAGE_FUNCTIONS_2D_HPP
#define GRID_DATA_PACKAGE_FUNCTIONS_2D_HPP

#include "grid_data_package_function.h"

namespace SPH
{
//=============================================================================================//
template <int PKG_SIZE>
std::pair<size_t, Array2i> NeighbourIndexShift(const Array2i shift_index, const CellNeighborhood2d &neighbour)
{
    std::pair<size_t, Arrayi> result;
    Arrayi neighbour_index = (shift_index + PKG_SIZE * Arrayi::Ones()) / PKG_SIZE;
    result.first = neighbour[neighbour_index[0]][neighbour_index[1]];
    result.second = (shift_index + PKG_SIZE * Arrayi::Ones()) - neighbour_index * PKG_SIZE;

    return result;
}
//=============================================================================================//
} // namespace SPH
#endif // GRID_DATA_PACKAGE_FUNCTIONS_2D_HPP