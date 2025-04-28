#ifndef MESH_LOCAL_DYNAMICS_2D_HPP
#define MESH_LOCAL_DYNAMICS_2D_HPP

#include "mesh_local_dynamics.h"

namespace SPH
{
#define VAR_AT(target, package_index, grid_index)  \
  target[package_index][grid_index[0]][grid_index[1]]
#define GET_NEIGHBOR_VAL(target, neighbor_index) \
  VAR_AT(target, neighbor_index.first, neighbor_index.second)
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
inline std::pair<size_t, Arrayi> BaseMeshLocalDynamics::
    NeighbourIndexShift(const Arrayi shift_index, const CellNeighborhood &neighbour)
{
    std::pair<size_t, Arrayi> result;
    Arrayi neighbour_index = (shift_index + pkg_size * Arrayi::Ones()) / pkg_size;
    result.first = neighbour[neighbour_index[0]][neighbour_index[1]];
    result.second = (shift_index + pkg_size * Arrayi::Ones()) - neighbour_index * pkg_size;

    return result;
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
    for (int n = 0; n != 2; n++)
    {
        size_t cell_index_in_this_direction = global_grid_index[n] / pkg_size;
        cell_index_on_mesh_[n] = cell_index_in_this_direction;
        local_data_index[n] = global_grid_index[n] - cell_index_in_this_direction * pkg_size;
    }
    size_t package_index = data_mesh->PackageIndexFromCellIndex(cell_package_index, cell_index_on_mesh_);
    auto &data = mesh_variable_data[package_index];
    return data[local_data_index[0]][local_data_index[1]];
}
//=============================================================================================//
template <typename DataType>
DataType BaseMeshLocalDynamics::CornerAverage(MeshVariableData<DataType> *mesh_variable_data,
                                              Arrayi addrs_index, Arrayi corner_direction,
                                              CellNeighborhood &neighborhood, DataType zero)
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
template <class DataType>
DataType ProbeMesh::probeMesh(MeshVariableData<DataType> *mesh_variable_data, const Vecd &position)
{
    Arrayi cell_index = index_handler_->CellIndexFromPosition(position);
    size_t package_index = index_handler_->PackageIndexFromCellIndex(cell_package_index_, cell_index);
    return package_index > 1 ? probeDataPackage(mesh_variable_data, package_index, cell_index, position)
                             : mesh_variable_data[package_index][0][0];
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
    NeighbourIndex neighbour_index_1 = BaseMeshLocalDynamics::NeighbourIndexShift(
        data_index + Arrayi(0, 0), neighborhood);
    NeighbourIndex neighbour_index_2 = BaseMeshLocalDynamics::NeighbourIndexShift(
        data_index + Arrayi(1, 0), neighborhood);
    NeighbourIndex neighbour_index_3 = BaseMeshLocalDynamics::NeighbourIndexShift(
        data_index + Arrayi(0, 1), neighborhood);
    NeighbourIndex neighbour_index_4 = BaseMeshLocalDynamics::NeighbourIndexShift(
        data_index + Arrayi(1, 1), neighborhood);

    DataType bilinear = GET_NEIGHBOR_VAL(mesh_variable_data, neighbour_index_1) * beta[0] * beta[1] +
                        GET_NEIGHBOR_VAL(mesh_variable_data, neighbour_index_2) * alpha[0] * beta[1] +
                        GET_NEIGHBOR_VAL(mesh_variable_data, neighbour_index_3) * beta[0] * alpha[1] +
                        GET_NEIGHBOR_VAL(mesh_variable_data, neighbour_index_4) * alpha[0] * alpha[1];

    return bilinear;
}
//=============================================================================================//
inline void UpdateLevelSetGradient::UpdateKernel::update(const size_t &package_index)
{
    auto &neighborhood = cell_neighborhood_[package_index];
    auto &pkg_data = phi_gradient_[package_index];

    BaseMeshLocalDynamics::for_each_cell_data(
        [&](int i, int j)
        {
            NeighbourIndex x1 = BaseMeshLocalDynamics::NeighbourIndexShift(
                Arrayi(i + 1, j), neighborhood);
            NeighbourIndex x2 = BaseMeshLocalDynamics::NeighbourIndexShift(
                Arrayi(i - 1, j), neighborhood);
            NeighbourIndex y1 = BaseMeshLocalDynamics::NeighbourIndexShift(
                Arrayi(i, j + 1), neighborhood);
            NeighbourIndex y2 = BaseMeshLocalDynamics::NeighbourIndexShift(
                Arrayi(i, j - 1), neighborhood);
            Real dphidx = GET_NEIGHBOR_VAL(phi_, x1) - GET_NEIGHBOR_VAL(phi_, x2);
            Real dphidy = GET_NEIGHBOR_VAL(phi_, y1) - GET_NEIGHBOR_VAL(phi_, y2);
            pkg_data[i][j] = 0.5 * Vecd(dphidx, dphidy) / data_spacing_;
        });
}
//=============================================================================================//
template <class KernelType>
template <typename DataType, typename FunctionOnData>
void UpdateKernelIntegrals<KernelType>::UpdateKernel::assignByPosition(MeshVariable<DataType> &mesh_variable,
                                                         const Arrayi &cell_index,
                                                         const FunctionByPosition &function_by_position)
{
    size_t package_index = index_handler_->PackageIndexFromCellIndex(cell_index);
    auto &pkg_data = mesh_variable[package_index];
    BaseMeshLocalDynamics::for_each_cell_data(
      [&](int i, int j)
      {
          Vec2d position = index_handler_->DataPositionFromIndex(cell_index, Arrayi(i, j));
          pkg_data[i][j] = function_by_position(position);
      });
}
//=============================================================================================//
template <class KernelType>
Real UpdateKernelIntegrals<KernelType>::UpdateKernel::computeKernelIntegral(const Vecd &position)
{
    Real phi = probe_signed_distance_.update(position);
    Real cutoff_radius = kernel_->CutOffRadius(global_h_ratio_);
    Real threshold = cutoff_radius + data_spacing_; // consider that interface's half width is the data spacing

    Real integral(0);
    if (fabs(phi) < threshold)
    {
        Arrayi global_index_ = mesh_data_.CellIndexFromPositionOnGlobalMesh(position);
        mesh_for_each2d<-3, 4>(
            [&](int i, int j, int k)
            {
                Arrayi neighbor_index = Arrayi(global_index_[0] + i, global_index_[1] + j, global_index_[2] + k);
                Real phi_neighbor = mesh_data_.DataValueFromGlobalIndex(phi_, neighbor_index);
                if (phi_neighbor > -data_spacing_)
                {
                    Vecd phi_gradient = mesh_data_.DataValueFromGlobalIndex(phi_gradient_, neighbor_index);
                    Vecd integral_position = mesh_data_.GridPositionFromIndexOnGlobalMesh(neighbor_index);
                    Vecd displacement = position - integral_position;
                    Real distance = displacement.norm();
                    if (distance < cutoff_radius)
                        integral += kernel_->W(global_h_ratio_, distance, displacement) *
                                    CutCellVolumeFraction(phi_neighbor, phi_gradient, data_spacing_);
                }
            });
    }
    return phi > threshold ? 1.0 : integral * data_spacing_ * data_spacing_ * data_spacing_;
}
//=============================================================================================//
template <class KernelType>
Vecd UpdateKernelIntegrals<KernelType>::UpdateKernel::computeKernelGradientIntegral(const Vecd &position)
{
    Real phi = probe_signed_distance_.update(position);
    Real cutoff_radius = kernel_->CutOffRadius(global_h_ratio_);
    Real threshold = cutoff_radius + data_spacing_;

    Vecd integral = Vecd::Zero();
    if (fabs(phi) < threshold)
    {
        Arrayi global_index_ = mesh_data_.CellIndexFromPositionOnGlobalMesh(position);
        mesh_for_each3d<-3, 4>(
            [&](int i, int j, int k)
            {
                Arrayi neighbor_index = Arrayi(global_index_[0] + i, global_index_[1] + j, global_index_[2] + k);
                Real phi_neighbor = mesh_data_.DataValueFromGlobalIndex(phi_, neighbor_index);
                if (phi_neighbor > -data_spacing_)
                {
                    Vecd phi_gradient = mesh_data_.DataValueFromGlobalIndex(phi_gradient_, neighbor_index);
                    Vecd integral_position = mesh_data_.GridPositionFromIndexOnGlobalMesh(neighbor_index);
                    Vecd displacement = position - integral_position;
                    Real distance = displacement.norm();
                    if (distance < cutoff_radius)
                        integral += kernel_->dW(global_h_ratio_, distance, displacement) *
                                    CutCellVolumeFraction(phi_neighbor, phi_gradient, data_spacing_) *
                                    displacement / (distance + TinyReal);
                }
            });
    }

    return integral * data_spacing_ * data_spacing_ * data_spacing_;
}
//=============================================================================================//
template <class KernelType>
Matd UpdateKernelIntegrals<KernelType>::UpdateKernel::computeKernelSecondGradientIntegral(const Vecd& position)
{
    Real phi = probe_signed_distance_.update(position);
    Real cutoff_radius = kernel_->CutOffRadius(global_h_ratio_);
    Real threshold = cutoff_radius + data_spacing_;

    Matd integral = Matd::Zero();
    if (fabs(phi) < threshold)
    {
        Arrayi global_index_ = mesh_data_.CellIndexFromPositionOnGlobalMesh(position);
        mesh_for_each3d<-3, 4>(
            [&](int i, int j, int k)
            {
                Arrayi neighbor_index = Arrayi(global_index_[0] + i, global_index_[1] + j, global_index_[2] + k);
                Real phi_neighbor = mesh_data_.DataValueFromGlobalIndex(phi_, neighbor_index);
                if (phi_neighbor > -data_spacing_)
                {
                    Vecd phi_gradient = mesh_data_.DataValueFromGlobalIndex(phi_gradient_, neighbor_index);
                    Vecd integral_position = mesh_data_.GridPositionFromIndexOnGlobalMesh(neighbor_index);
                    Vecd displacement = position - integral_position;
                    Real distance = displacement.norm();
                    if (distance < cutoff_radius)
                        integral += kernel_->d2W(global_h_ratio_, distance, displacement) *
                        CutCellVolumeFraction(phi_neighbor, phi_gradient, data_spacing_) *
                        displacement * displacement.transpose() / (distance * distance + TinyReal);
                }
            });
    }

    return integral * data_spacing_ * data_spacing_ * data_spacing_;
}
//=============================================================================================//
inline void ReinitializeLevelSet::UpdateKernel::update(const size_t &package_index)
{
    auto &phi_addrs = phi_[package_index];
    auto &near_interface_id_addrs = near_interface_id_[package_index];
    auto &neighborhood = cell_neighborhood_[package_index];

    BaseMeshLocalDynamics::for_each_cell_data(
        [&](int i, int j)
        {
            // only reinitialize non cut cells
            if (near_interface_id_addrs[i][j] != 0)
            {
                Real phi_0 = phi_addrs[i][j];
                Real sign = phi_0 / sqrt(phi_0 * phi_0 + data_spacing_ * data_spacing_);
                NeighbourIndex x1 = BaseMeshLocalDynamics::NeighbourIndexShift(
                    Arrayi(i + 1, j), neighborhood);
                NeighbourIndex x2 = BaseMeshLocalDynamics::NeighbourIndexShift(
                    Arrayi(i - 1, j), neighborhood);
                NeighbourIndex y1 = BaseMeshLocalDynamics::NeighbourIndexShift(
                    Arrayi(i, j + 1), neighborhood);
                NeighbourIndex y2 = BaseMeshLocalDynamics::NeighbourIndexShift(
                    Arrayi(i, j - 1), neighborhood);
                Real dv_x = upwindDifference(sign, GET_NEIGHBOR_VAL(phi_, x1) - phi_0,
                                             phi_0 - GET_NEIGHBOR_VAL(phi_, x2));
                Real dv_y = upwindDifference(sign, GET_NEIGHBOR_VAL(phi_, y1) - phi_0,
                                             phi_0 - GET_NEIGHBOR_VAL(phi_, y2));
                phi_addrs[i][j] -= 0.5 * sign * (Vec2d(dv_x, dv_y).norm() - data_spacing_);
            }
        });
}
//=============================================================================================//
inline void MarkNearInterface::UpdateKernel::update(const size_t &package_index,
                                                    Real small_shift_factor)
{
    Real small_shift = data_spacing_ * small_shift_factor; 
    auto &phi_addrs = phi_[package_index];
    auto &near_interface_id_addrs = near_interface_id_[package_index];

    // corner averages, note that the first row and first column are not used
    PackageDataMatrix<Real, 5> corner_averages;
    mesh_for_each2d<0, 5>(
        [&](int i, int j)
        {
            corner_averages[i][j] = BaseMeshLocalDynamics::CornerAverage(phi_, Arrayi(i, j), Arrayi(-1, -1),
                                                                         cell_neighborhood_[package_index],
                                                                         (Real)0);
        });

    BaseMeshLocalDynamics::for_each_cell_data(
        [&](int i, int j)
        {
            // first assume far cells
            Real phi_0 = phi_addrs[i][j];
            int near_interface_id = phi_0 > 0.0 ? 2 : -2;
            if (fabs(phi_0) < small_shift)
            {
                near_interface_id = 0;
                Real phi_average_0 = corner_averages[i][j];
                // find inner and outer cut cells
                mesh_for_each2d<0, 2>(
                    [&](int l, int m)
                    {
                        Real phi_average = corner_averages[i + l][j + m];
                        if ((phi_average_0 - small_shift) * (phi_average - small_shift) < 0.0)
                            near_interface_id = 1;
                        if ((phi_average_0 + small_shift) * (phi_average + small_shift) < 0.0)
                            near_interface_id = -1;
                    });
                // find zero cut cells
                mesh_for_each2d<0, 2>(
                    [&](int l, int m)
                    {
                        Real phi_average = corner_averages[i + l][j + m];
                        if (phi_average_0 * phi_average < 0.0)
                            near_interface_id = 0;
                    });
            }
            // assign this is to package
            near_interface_id_addrs[i][j] = near_interface_id;
        });
}
//=============================================================================================//
inline void RedistanceInterface::UpdateKernel::update(const size_t &package_index)
{
    BaseMeshLocalDynamics::for_each_cell_data(
        [&](int i, int j)
        {
            int near_interface_id = near_interface_id_[package_index][i][j];
            if (near_interface_id == 0)
            {
                bool positive_band = false;
                bool negative_band = false;
                mesh_for_each2d<-1, 2>(
                    [&](int r, int s)
                    {
                        NeighbourIndex neighbour_index = BaseMeshLocalDynamics::NeighbourIndexShift(
                            Arrayi(i + r, j + s), cell_neighborhood_[package_index]);
                        int neighbor_near_interface_id =
                            GET_NEIGHBOR_VAL(near_interface_id_, neighbour_index);
                        if (neighbor_near_interface_id >= 1)
                            positive_band = true;
                        if (neighbor_near_interface_id <= -1)
                            negative_band = true;
                    });
                if (positive_band == false)
                {
                    Real min_distance_p = 5.0 * data_spacing_;
                    mesh_for_each2d<-4, 5>(
                        [&](int x, int y)
                        {
                            NeighbourIndex neighbour_index = BaseMeshLocalDynamics::NeighbourIndexShift(
                                Arrayi(i + x, j + y), cell_neighborhood_[package_index]);
                            auto &neighbor_phi = phi_[neighbour_index.first];
                            auto &neighbor_phi_gradient = phi_gradient_[neighbour_index.first];
                            auto &neighbor_near_interface_id = near_interface_id_[neighbour_index.first];
                            if (neighbor_near_interface_id[neighbour_index.second[0]][neighbour_index.second[1]] >= 1)
                            {
                                Real phi_p_ = neighbor_phi[neighbour_index.second[0]][neighbour_index.second[1]];
                                Vecd norm_to_face = neighbor_phi_gradient[neighbour_index.second[0]][neighbour_index.second[1]];
                                norm_to_face /= norm_to_face.norm() + TinyReal;
                                min_distance_p = SMIN(min_distance_p, (Vecd((Real)x, (Real)y) * data_spacing_ + phi_p_ * norm_to_face).norm());
                            }
                        });
                    phi_[package_index][i][j] = -min_distance_p;
                    // this immediate switch of near interface id
                    // does not intervening with the identification of unresolved interface
                    // based on the assumption that positive false_and negative bands are not close to each other
                    near_interface_id_[package_index][i][j] = -1;
                }
                if (negative_band == false)
                {
                    Real min_distance_n = 5.0 * data_spacing_;
                    mesh_for_each2d<-4, 5>(
                        [&](int x, int y)
                        {
                            NeighbourIndex neighbour_index = BaseMeshLocalDynamics::NeighbourIndexShift(
                                Arrayi(i + x, j + y), cell_neighborhood_[package_index]);
                            auto &neighbor_phi = phi_[neighbour_index.first];
                            auto &neighbor_phi_gradient = phi_gradient_[neighbour_index.first];
                            auto &neighbor_near_interface_id = near_interface_id_[neighbour_index.first];
                            if (neighbor_near_interface_id[neighbour_index.second[0]][neighbour_index.second[1]] <= -1)
                            {
                                Real phi_n_ = neighbor_phi[neighbour_index.second[0]][neighbour_index.second[1]];
                                Vecd norm_to_face = neighbor_phi_gradient[neighbour_index.second[0]][neighbour_index.second[1]];
                                norm_to_face /= norm_to_face.norm() + TinyReal;
                                min_distance_n = SMIN(min_distance_n, (Vecd((Real)x, (Real)y) * data_spacing_ - phi_n_ * norm_to_face).norm());
                            }
                        });
                    phi_[package_index][i][j] = min_distance_n;
                    // this immediate switch of near interface id
                    // does not intervening with the identification of unresolved interface
                    // based on the assumption that positive false_and negative bands are not close to each other
                    near_interface_id_[package_index][i][j] = 1;
                }
            }
        });
}
//=============================================================================================//
inline void DiffuseLevelSetSign::UpdateKernel::update(const size_t &package_index)
{
    BaseMeshLocalDynamics::for_each_cell_data(
        [&](int i, int j)
        {
            // near interface cells are not considered
            if (abs(near_interface_id_[package_index][i][j]) > 1)
            {
                mesh_find_if2d<-1, 2>(
                    [&](int l, int m) -> bool
                    {
                        NeighbourIndex neighbour_index =
                            BaseMeshLocalDynamics::NeighbourIndexShift(
                                Arrayi(i + l, j + m), cell_neighborhood_[package_index]);
                        int near_interface_id =
                            GET_NEIGHBOR_VAL(near_interface_id_, neighbour_index);
                        bool is_found = abs(near_interface_id) == 1;
                        if (is_found)
                        {
                            Real phi_0 = phi_[package_index][i][j];
                            near_interface_id_[package_index][i][j] = near_interface_id;
                            phi_[package_index][i][j] =
                                near_interface_id == 1 ? fabs(phi_0) : -fabs(phi_0);
                        }
                        return is_found;
                    });
            }
        });
}
//=============================================================================================//
} // namespace SPH
//=============================================================================================//
#endif //MESH_LOCAL_DYNAMICS_2D_HPP