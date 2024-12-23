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
template <class DataType>
DataType BaseMeshLocalDynamics::probeMesh(MeshWithGridDataPackagesType &probe_mesh_,
                                          MeshVariableData<DataType> *mesh_variable_data,
                                          const Vecd &position)
{
    Arrayi cell_index = probe_mesh_.CellIndexFromPosition(position);
    size_t package_index = probe_mesh_.PackageIndexFromCellIndex(cell_index);
    return probe_mesh_.isInnerDataPackage(cell_index) ? probeDataPackage(probe_mesh_, mesh_variable_data, package_index, cell_index, position)
                                                      : mesh_variable_data[package_index][0][0];
}
//=============================================================================================//
template <class DataType>
DataType BaseMeshLocalDynamics::probeDataPackage(MeshWithGridDataPackagesType &probe_mesh_,
                                                 MeshVariableData<DataType> *mesh_variable_data,
                                                 size_t package_index,
                                                 const Arrayi &cell_index,
                                                 const Vecd &position)
{
    Arrayi data_index = probe_mesh_.DataIndexFromPosition(cell_index, position);
    Vecd data_position = probe_mesh_.DataPositionFromIndex(cell_index, data_index);
    Vecd alpha = (position - data_position) / probe_mesh_.DataSpacing();
    Vecd beta = Vecd::Ones() - alpha;

    auto &neighborhood = probe_mesh_.cell_neighborhood_[package_index];
    std::pair<size_t, Arrayi> neighbour_index_1 = NeighbourIndexShift(Arrayi(data_index[0], data_index[1]), neighborhood);
    std::pair<size_t, Arrayi> neighbour_index_2 = NeighbourIndexShift(Arrayi(data_index[0] + 1, data_index[1]), neighborhood);
    std::pair<size_t, Arrayi> neighbour_index_3 = NeighbourIndexShift(Arrayi(data_index[0], data_index[1] + 1), neighborhood);
    std::pair<size_t, Arrayi> neighbour_index_4 = NeighbourIndexShift(Arrayi(data_index[0] + 1, data_index[1] + 1), neighborhood);

    DataType bilinear = mesh_variable_data[neighbour_index_1.first][neighbour_index_1.second[0]][neighbour_index_1.second[1]] * beta[0] * beta[1] +
                        mesh_variable_data[neighbour_index_2.first][neighbour_index_2.second[0]][neighbour_index_2.second[1]] * alpha[0] * beta[1] +
                        mesh_variable_data[neighbour_index_3.first][neighbour_index_3.second[0]][neighbour_index_3.second[1]] * beta[0] * alpha[1] +
                        mesh_variable_data[neighbour_index_4.first][neighbour_index_4.second[0]][neighbour_index_4.second[1]] * alpha[0] * alpha[1];

    return bilinear;
}
//=============================================================================================//
template <typename DataType>
DataType BaseMeshLocalDynamics::DataValueFromGlobalIndex(MeshVariableData<DataType> *mesh_variable_data,
                                                         const Arrayi &global_grid_index,
                                                         MeshWithGridDataPackagesType *data_mesh)
{
    Arrayi cell_index_on_mesh_ = Arrayi::Zero();
    Arrayi local_data_index = Arrayi::Zero();
    for (int n = 0; n != 2; n++)
    {
        size_t cell_index_in_this_direction = global_grid_index[n] / pkg_size;
        cell_index_on_mesh_[n] = cell_index_in_this_direction;
        local_data_index[n] = global_grid_index[n] - cell_index_in_this_direction * pkg_size;
    }
    size_t package_index = data_mesh->PackageIndexFromCellIndex(cell_index_on_mesh_);
    auto &data = mesh_variable_data[package_index];
    return data[local_data_index[0]][local_data_index[1]];
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
        {
            int x_index = addrs_index[0] + i * corner_direction[0];
            int y_index = addrs_index[1] + j * corner_direction[1];
            std::pair<size_t, Arrayi> neighbour_index = NeighbourIndexShift(Arrayi(x_index, y_index), neighborhood);
            average += mesh_variable_data[neighbour_index.first][neighbour_index.second[0]][neighbour_index.second[1]];
        }
    return average * 0.25;
}
//=============================================================================================//
template <typename DataType, typename FunctionByPosition>
void UpdateKernelIntegrals::UpdateKernel::assignByPosition(MeshVariableData<DataType> *mesh_variable_data,
                                            const Arrayi &cell_index,
                                            MeshWithGridDataPackagesType *data_mesh,
                                            const FunctionByPosition &function_by_position)
{
    size_t package_index = data_mesh_->PackageIndexFromCellIndex(cell_index);
    auto &pkg_data = mesh_variable_data[package_index];
    for (int i = 0; i != pkg_size; ++i)
        for (int j = 0; j != pkg_size; ++j)
        {
            Vec2d position = data_mesh_->DataPositionFromIndex(cell_index, Arrayi(i, j));
            pkg_data[i][j] = function_by_position(position);
        }
}
//=============================================================================================//
} // namespace SPH
//=============================================================================================//
#endif //MESH_LOCAL_DYNAMICS_3D_HPP