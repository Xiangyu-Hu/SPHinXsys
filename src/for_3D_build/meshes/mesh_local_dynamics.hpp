#ifndef MESH_LOCAL_DYNAMICS_3D_HPP
#define MESH_LOCAL_DYNAMICS_3D_HPP

#include "mesh_local_dynamics.h"

namespace SPH
{
#define VAR_AT(target, package_index, grid_index)  \
  target[package_index][grid_index[0]][grid_index[1]][grid_index[2]]
#define GET_NEIGHBOR_VAL(target, neighbor_index) \
  VAR_AT(target, neighbor_index.first, neighbor_index.second)
//=============================================================================================//
template <typename FunctionOnData>
void BaseMeshLocalDynamics::for_each_cell_data(const FunctionOnData &function)
{
    for (int i = 0; i != pkg_size; ++i)
        for (int j = 0; j != pkg_size; ++j)
            for (int k = 0; k != pkg_size; ++k)
            {
                function(i, j, k);
            }
}
//=============================================================================================//
template <typename DataType>
DataType BaseMeshLocalDynamics::DataValueFromGlobalIndex(MeshVariableData<DataType> *mesh_variable_data,
                                                         const Arrayi &global_grid_index,
                                                         MeshWithGridDataPackagesType *data_mesh,
                                                         size_t *cell_package_index)
{
    Arrayi cell_index_on_mesh_ = Arrayi::Zero();
    Arrayi local_data_index = Arrayi::Zero();
    for (int n = 0; n != 3; n++)
    {
        size_t cell_index_in_this_direction = global_grid_index[n] / pkg_size;
        cell_index_on_mesh_[n] = cell_index_in_this_direction;
        local_data_index[n] = global_grid_index[n] - cell_index_in_this_direction * pkg_size;
    }
    size_t package_index = data_mesh->PackageIndexFromCellIndex(cell_package_index, cell_index_on_mesh_);
    auto &data = mesh_variable_data[package_index];
    return data[local_data_index[0]][local_data_index[1]][local_data_index[2]];
}
//=============================================================================================//
template <typename DataType>
DataType BaseMeshLocalDynamics::CornerAverage(MeshVariableData<DataType> *mesh_variable_data,
                                              Arrayi addrs_index, Arrayi corner_direction,
                                              CellNeighborhood &neighborhood)
{
    DataType average = ZeroData<DataType>::value;
    for (int i = 0; i != 2; ++i)
        for (int j = 0; j != 2; ++j)
            for (int k = 0; k != 2; ++k)
            {
                int x_index = addrs_index[0] + i * corner_direction[0];
                int y_index = addrs_index[1] + j * corner_direction[1];
                int z_index = addrs_index[2] + k * corner_direction[2];
                std::pair<size_t, Arrayi> neighbour_index = NeighbourIndexShift(Arrayi(x_index, y_index, z_index), neighborhood);
                average += GET_NEIGHBOR_VAL(mesh_variable_data, neighbour_index);
            }
    return average * 0.125;
}
//=============================================================================================//
template <class DataType>
DataType ProbeMesh::probeMesh(MeshVariableData<DataType> *mesh_variable_data, const Vecd &position)
{
    Arrayi cell_index = index_handler_->CellIndexFromPosition(position);
    size_t package_index = index_handler_->PackageIndexFromCellIndex(cell_package_index_, cell_index);
    return package_index > 1 ? probeDataPackage(mesh_variable_data, package_index, cell_index, position)
                             : mesh_variable_data[package_index][0][0][0];
}
//=============================================================================================//
template <class DataType>
DataType ProbeMesh::probeDataPackage(MeshVariableData<DataType> *mesh_variable_data,
                                     size_t package_index,
                                     const Arrayi &cell_index,
                                     const Vecd &position)
{
    Arrayi data_index = index_handler_->DataIndexFromPosition(cell_index, position);
    Vecd data_position = index_handler_->DataPositionFromIndex(cell_index, data_index);
    Vecd alpha = (position - data_position) / index_handler_->data_spacing_;
    Vecd beta = Vecd::Ones() - alpha;

    auto &neighborhood = cell_neighborhood_[package_index];
    NeighbourIndex neighbour_index_1 = BaseMeshLocalDynamics::NeighbourIndexShift(data_index + Arrayi(0, 0, 0), neighborhood);
    NeighbourIndex neighbour_index_2 = BaseMeshLocalDynamics::NeighbourIndexShift(data_index + Arrayi(1, 0, 0), neighborhood);
    NeighbourIndex neighbour_index_3 = BaseMeshLocalDynamics::NeighbourIndexShift(data_index + Arrayi(0, 1, 0), neighborhood);
    NeighbourIndex neighbour_index_4 = BaseMeshLocalDynamics::NeighbourIndexShift(data_index + Arrayi(1, 1, 0), neighborhood);

    DataType bilinear_1 = GET_NEIGHBOR_VAL(mesh_variable_data, neighbour_index_1) * beta[0] * beta[1] +
                          GET_NEIGHBOR_VAL(mesh_variable_data, neighbour_index_2) * alpha[0] * beta[1] +
                          GET_NEIGHBOR_VAL(mesh_variable_data, neighbour_index_3) * beta[0] * alpha[1] +
                          GET_NEIGHBOR_VAL(mesh_variable_data, neighbour_index_4) * alpha[0] * alpha[1];

    neighbour_index_1 = BaseMeshLocalDynamics::NeighbourIndexShift(data_index + Arrayi(0, 0, 1), neighborhood);
    neighbour_index_2 = BaseMeshLocalDynamics::NeighbourIndexShift(data_index + Arrayi(1, 0, 1), neighborhood);
    neighbour_index_3 = BaseMeshLocalDynamics::NeighbourIndexShift(data_index + Arrayi(0, 1, 1), neighborhood);
    neighbour_index_4 = BaseMeshLocalDynamics::NeighbourIndexShift(data_index + Arrayi(1, 1, 1), neighborhood);

    DataType bilinear_2 = GET_NEIGHBOR_VAL(mesh_variable_data, neighbour_index_1) * beta[0] * beta[1] +
                          GET_NEIGHBOR_VAL(mesh_variable_data, neighbour_index_2) * alpha[0] * beta[1] +
                          GET_NEIGHBOR_VAL(mesh_variable_data, neighbour_index_3) * beta[0] * alpha[1] +
                          GET_NEIGHBOR_VAL(mesh_variable_data, neighbour_index_4) * alpha[0] * alpha[1];

    return bilinear_1 * beta[2] + bilinear_2 * alpha[2];
}
//=============================================================================================//
} // namespace SPH
//=============================================================================================//
#endif //MESH_LOCAL_DYNAMICS_3D_HPP