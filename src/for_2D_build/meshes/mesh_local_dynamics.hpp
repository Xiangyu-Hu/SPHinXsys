#ifndef MESH_LOCAL_DYNAMICS_2D_HPP
#define MESH_LOCAL_DYNAMICS_2D_HPP

#include "mesh_local_dynamics.h"

namespace SPH
{
//=============================================================================================//
template <typename FunctionOnData>
void BaseMeshLocalDynamics::for_each_cell_data(const FunctionOnData &function)
{
    for (int i = 0; i != pkg_size; ++i)
        for (int j = 0; j != pkg_size; ++j)
        {
            function(i, j);
        }
}
//=============================================================================================//
template <typename InDataType, typename OutDataType>
void UpdateLevelSetGradient::computeGradient(MeshVariable<InDataType> &in_variable,
                                             MeshVariable<OutDataType> &out_variable,
                                             const size_t package_index)
{
    auto in_variable_data = in_variable.DataField();
    auto out_variable_data = out_variable.DataField();

    auto &neighborhood = cell_neighborhood_[package_index];
    auto &pkg_data = out_variable_data[package_index];

    for_each_cell_data(
        [&](int i, int j)
        {
            std::pair<size_t, Arrayi> x1 = NeighbourIndexShift(Arrayi(i + 1, j), neighborhood);
            std::pair<size_t, Arrayi> x2 = NeighbourIndexShift(Arrayi(i - 1, j), neighborhood);
            std::pair<size_t, Arrayi> y1 = NeighbourIndexShift(Arrayi(i, j + 1), neighborhood);
            std::pair<size_t, Arrayi> y2 = NeighbourIndexShift(Arrayi(i, j - 1), neighborhood);
            Real dphidx = (in_variable_data[x1.first][x1.second[0]][x1.second[1]] -
                           in_variable_data[x2.first][x2.second[0]][x2.second[1]]);
            Real dphidy = (in_variable_data[y1.first][y1.second[0]][y1.second[1]] -
                           in_variable_data[y2.first][y2.second[0]][y2.second[1]]);

            pkg_data[i][j] = 0.5 * Vecd(dphidx, dphidy) / data_spacing_;
        });
}
//=============================================================================================//
} // namespace SPH
//=============================================================================================//
#endif //MESH_LOCAL_DYNAMICS_3D_HPP