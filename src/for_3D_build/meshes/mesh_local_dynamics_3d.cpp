#include "mesh_local_dynamics.h"

#include "mesh_iterators.h"

namespace SPH
{
//=============================================================================================//
size_t BaseMeshLocalDynamics::SortIndexFromCellIndex(const Arrayi &cell_index)
{
    return cell_index[0] * all_cells_[1] * all_cells_[2] + cell_index[1] * all_cells_[2] + cell_index[2];
}
//=============================================================================================//
Arrayi BaseMeshLocalDynamics::CellIndexFromSortIndex(const size_t &sort_index)
{
    Array3i cell_index;
    cell_index[0] = sort_index / (all_cells_[1] * all_cells_[2]);
    cell_index[1] = (sort_index / all_cells_[2]) % all_cells_[1];
    cell_index[2] = sort_index % all_cells_[2];

    return cell_index;
}
//=============================================================================================//
void InitializeDataForSingularPackage::update(const size_t package_index, Real far_field_level_set)
{
    auto &phi = phi_.Data()[package_index];
    auto &near_interface_id = near_interface_id_.Data()[package_index];
    auto &phi_gradient = phi_gradient_.Data()[package_index];
    auto &kernel_weight = kernel_weight_.Data()[package_index];
    auto &kernel_gradient = kernel_gradient_.Data()[package_index];

    mesh_data_.for_each_cell_data(
        [&](int i, int j, int k)
        {
            phi[i][j][k] = far_field_level_set;
            near_interface_id[i][j][k] = far_field_level_set < 0.0 ? -2 : 2;
            phi_gradient[i][j][k] = Vecd::Ones();
            kernel_weight[i][j][k] = far_field_level_set < 0.0 ? 0 : 1.0;
            kernel_gradient[i][j][k] = Vec3d::Zero();
        });
}
//=============================================================================================//
bool TagACellIsInnerPackage::isInnerPackage(const Arrayi &cell_index)
{
    return mesh_any_of(
        Array3i::Zero().max(cell_index - Array3i::Ones()),
        all_cells_.min(cell_index + 2 * Array3i::Ones()),
        [&](int l, int m, int n)
        {
            return mesh_data_.isInnerDataPackage(Arrayi(l, m, n)); // actually a core test here, because only core pkgs are assigned
        });
}
//=============================================================================================//
void InitializeCellNeighborhood::update(const size_t &package_index)
{
    size_t sort_index = mesh_data_.occupied_data_pkgs_[package_index - 2].first;
    Arrayi cell_index = CellIndexFromSortIndex(sort_index);
    CellNeighborhood &current = mesh_data_.cell_neighborhood_[package_index];
    std::pair<Arrayi, int> &metadata = mesh_data_.meta_data_cell_[package_index];
    metadata.first = cell_index;
    metadata.second = mesh_data_.occupied_data_pkgs_[package_index - 2].second;
    for (int l = -1; l < 2; l++)
        for (int m = -1; m < 2; m++)
            for (int n = -1; n < 2; n++)
            {
                current[l + 1][m + 1][n + 1] = mesh_data_.PackageIndexFromCellIndex(cell_index + Arrayi(l, m, n));
            }
}
//=============================================================================================//
void InitializeBasicDataForAPackage::update(const size_t &package_index)
{
    auto &phi = phi_.Data()[package_index];
    auto &near_interface_id = near_interface_id_.Data()[package_index];
    Arrayi cell_index = mesh_data_.meta_data_cell_[package_index].first;
    mesh_data_.for_each_cell_data(
        [&](int i, int j, int k)
        {
            Vec3d position = mesh_data_.DataPositionFromIndex(cell_index, Array3i(i, j, k));
            phi[i][j][k] = shape_.findSignedDistance(position);
            near_interface_id[i][j][k] = phi[i][j][k] < 0.0 ? -2 : 2;
        });
}
//=============================================================================================//
Real UpdateKernelIntegrals::computeKernelIntegral(const Vecd &position)
{
    Real phi = probeSignedDistance(position);
    Real cutoff_radius = kernel_.CutOffRadius(global_h_ratio_);
    Real threshold = cutoff_radius + data_spacing_; // consider that interface's half width is the data spacing

    Real integral(0);
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
                        integral += kernel_.W(global_h_ratio_, distance, displacement) *
                                    CutCellVolumeFraction(phi_neighbor, phi_gradient, data_spacing_);
                }
            });
    }
    return phi > threshold ? 1.0 : integral * data_spacing_ * data_spacing_ * data_spacing_;
}
//=============================================================================================//
Vecd UpdateKernelIntegrals::computeKernelGradientIntegral(const Vecd &position)
{
    Real phi = probeSignedDistance(position);
    Real cutoff_radius = kernel_.CutOffRadius(global_h_ratio_);
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
                        integral += kernel_.dW(global_h_ratio_, distance, displacement) *
                                    CutCellVolumeFraction(phi_neighbor, phi_gradient, data_spacing_) *
                                    displacement / (distance + TinyReal);
                }
            });
    }

    return integral * data_spacing_ * data_spacing_ * data_spacing_;
}
//=============================================================================================//
Matd UpdateKernelIntegrals::computeKernelSecondGradientIntegral(const Vecd &position)
{
    Real phi = probeSignedDistance(position);
    Real cutoff_radius = kernel_.CutOffRadius(global_h_ratio_);
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
                        integral += kernel_.d2W(global_h_ratio_, distance, displacement) *
                                    CutCellVolumeFraction(phi_neighbor, phi_gradient, data_spacing_) *
                                    displacement * displacement.transpose() / (distance * distance + TinyReal);
                }
            });
    }

    return integral * data_spacing_ * data_spacing_ * data_spacing_;
}
//=============================================================================================//
void ReinitializeLevelSet::update(const size_t &package_index)
{
    auto phi_data = phi_.Data();
    auto &phi_addrs = phi_data[package_index];
    auto &near_interface_id_addrs = near_interface_id_.Data()[package_index];
    auto &neighborhood = mesh_data_.cell_neighborhood_[package_index];

    mesh_data_.for_each_cell_data(
        [&](int i, int j, int k)
        {
            // only reinitialize non cut cells
            if (near_interface_id_addrs[i][j][k] != 0)
            {
                Real phi_0 = phi_addrs[i][j][k];
                Real sign = phi_0 / sqrt(phi_0 * phi_0 + data_spacing_ * data_spacing_);
                using NeighbourIndex = std::pair<size_t, Arrayi>; /**< stores shifted neighbour info: (size_t)package index, (arrayi)local grid index. */
                NeighbourIndex x1 = mesh_data_.NeighbourIndexShift(Arrayi(i + 1, j, k), neighborhood);
                NeighbourIndex x2 = mesh_data_.NeighbourIndexShift(Arrayi(i - 1, j, k), neighborhood);
                NeighbourIndex y1 = mesh_data_.NeighbourIndexShift(Arrayi(i, j + 1, k), neighborhood);
                NeighbourIndex y2 = mesh_data_.NeighbourIndexShift(Arrayi(i, j - 1, k), neighborhood);
                NeighbourIndex z1 = mesh_data_.NeighbourIndexShift(Arrayi(i, j, k + 1), neighborhood);
                NeighbourIndex z2 = mesh_data_.NeighbourIndexShift(Arrayi(i, j, k - 1), neighborhood);
                Real dv_x = upwindDifference(sign, phi_data[x1.first][x1.second[0]][x1.second[1]][x1.second[2]] - phi_0, phi_0 - phi_data[x2.first][x2.second[0]][x2.second[1]][x2.second[2]]);
                Real dv_y = upwindDifference(sign, phi_data[y1.first][y1.second[0]][y1.second[1]][y1.second[2]] - phi_0, phi_0 - phi_data[y2.first][y2.second[0]][y2.second[1]][y2.second[2]]);
                Real dv_z = upwindDifference(sign, phi_data[z1.first][z1.second[0]][z1.second[1]][z1.second[2]] - phi_0, phi_0 - phi_data[z2.first][z2.second[0]][z2.second[1]][z2.second[2]]);
                phi_addrs[i][j][k] -= 0.3 * sign * (Vec3d(dv_x, dv_y, dv_z).norm() - data_spacing_);
            }
        });
}
//=============================================================================================//
void MarkNearInterface::update(const size_t &package_index)
{
    auto &phi_addrs = phi_.Data()[package_index];
    auto &near_interface_id_addrs = near_interface_id_.Data()[package_index];
    auto neighborhood = mesh_data_.cell_neighborhood_[package_index];

    // corner averages, note that the first row and first column are not used
    PackageDataMatrix<Real, 5> corner_averages;
    mesh_for_each3d<0, 5>(
        [&](int i, int j, int k)
        {
            corner_averages[i][j][k] = mesh_data_.CornerAverage(phi_, Arrayi(i, j, k), Arrayi(-1, -1, -1), neighborhood);
        });

    mesh_data_.for_each_cell_data(
        [&](int i, int j, int k)
        {
            // first assume far cells
            Real phi_0 = phi_addrs[i][j][k];
            int near_interface_id = phi_0 > 0.0 ? 2 : -2;
            if (fabs(phi_0) < small_shift)
            {
                near_interface_id = 0;
                Real phi_average_0 = corner_averages[i][j][k];
                // find inner and outer cut cells
                mesh_for_each3d<0, 2>(
                    [&](int l, int m, int n)
                    {
                        Real phi_average = corner_averages[i + l][j + m][k + n];
                        if ((phi_average_0 - small_shift) * (phi_average - small_shift) < 0.0)
                            near_interface_id = 1;
                        if ((phi_average_0 + small_shift) * (phi_average + small_shift) < 0.0)
                            near_interface_id = -1;
                    });
                // find zero cut cells
                mesh_for_each3d<0, 2>(
                    [&](int l, int m, int n)
                    {
                        Real phi_average = corner_averages[i + l][j + m][k + n];
                        if (phi_average_0 * phi_average < 0.0)
                            near_interface_id = 0;
                    });
            }
            // assign this is to package
            near_interface_id_addrs[i][j][k] = near_interface_id;
        });
}
//=============================================================================================//
void RedistanceInterface::update(const size_t &package_index)
{
    auto phi_data = phi_.Data();
    auto near_interface_id_data = near_interface_id_.Data();
    auto &neighborhood = mesh_data_.cell_neighborhood_[package_index];

    using NeighbourIndex = std::pair<size_t, Arrayi>; /**< stores shifted neighbour info: (size_t)package index, (arrayi)local grid index. */
    mesh_data_.for_each_cell_data(
        [&](int i, int j, int k)
        {
            int near_interface_id = near_interface_id_data[package_index][i][j][k];
            if (near_interface_id == 0)
            {
                bool positive_band = false;
                bool negative_band = false;
                mesh_for_each3d<-1, 2>(
                    [&](int r, int s, int t)
                    {
                        NeighbourIndex neighbour_index = mesh_data_.NeighbourIndexShift(Arrayi(i + r, j + s, k + t), neighborhood);
                        int neighbor_near_interface_id = near_interface_id_data[neighbour_index.first][neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]];
                        if (neighbor_near_interface_id >= 1)
                            positive_band = true;
                        if (neighbor_near_interface_id <= -1)
                            negative_band = true;
                    });
                if (positive_band == false)
                {
                    Real min_distance_p = 5.0 * data_spacing_;
                    mesh_for_each3d<-4, 5>(
                        [&](int x, int y, int z)
                        {
                            NeighbourIndex neighbour_index = mesh_data_.NeighbourIndexShift(Arrayi(i + x, j + y, k + z), neighborhood);
                            auto &neighbor_phi = phi_.Data()[neighbour_index.first];
                            auto &neighbor_phi_gradient = phi_gradient_.Data()[neighbour_index.first];
                            auto &neighbor_near_interface_id = near_interface_id_.Data()[neighbour_index.first];
                            if (neighbor_near_interface_id[neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]] >= 1)
                            {
                                Real phi_p_ = neighbor_phi[neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]];
                                Vecd norm_to_face = neighbor_phi_gradient[neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]];
                                norm_to_face /= norm_to_face.norm() + TinyReal;
                                min_distance_p = SMIN(min_distance_p, (Vecd((Real)x, (Real)y, Real(z)) * data_spacing_ + phi_p_ * norm_to_face).norm());
                            }
                        });
                    phi_data[package_index][i][j][k] = -min_distance_p;
                    // this immediate switch of near interface id
                    // does not intervening with the identification of unresolved interface
                    // based on the assumption that positive false_and negative bands are not close to each other
                    near_interface_id_data[package_index][i][j][k] = -1;
                }
                if (negative_band == false)
                {
                    Real min_distance_n = 5.0 * data_spacing_;
                    mesh_for_each3d<-4, 5>(
                        [&](int x, int y, int z)
                        {
                            NeighbourIndex neighbour_index = mesh_data_.NeighbourIndexShift(Arrayi(i + x, j + y, k + z), neighborhood);
                            auto &neighbor_phi = phi_.Data()[neighbour_index.first];
                            auto &neighbor_phi_gradient = phi_gradient_.Data()[neighbour_index.first];
                            auto &neighbor_near_interface_id = near_interface_id_.Data()[neighbour_index.first];
                            if (neighbor_near_interface_id[neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]] <= -1)
                            {
                                Real phi_n_ = neighbor_phi[neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]];
                                Vecd norm_to_face = neighbor_phi_gradient[neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]];
                                norm_to_face /= norm_to_face.norm() + TinyReal;
                                min_distance_n = SMIN(min_distance_n, (Vecd((Real)x, (Real)y, Real(z)) * data_spacing_ - phi_n_ * norm_to_face).norm());
                            }
                        });
                    phi_data[package_index][i][j][k] = min_distance_n;
                    // this immediate switch of near interface id
                    // does not intervening with the identification of unresolved interface
                    // based on the assumption that positive false_and negative bands are not close to each other
                    near_interface_id_data[package_index][i][j][k] = 1;
                }
            }
        });
}
//=============================================================================================//
void DiffuseLevelSetSign::update(const size_t &package_index)
{
    auto phi_data = phi_.Data();
    auto near_interface_id_data = near_interface_id_.Data();
    auto &neighborhood = mesh_data_.cell_neighborhood_[package_index];

    mesh_data_.for_each_cell_data(
        [&](int i, int j, int k)
        {
            // near interface cells are not considered
            if (abs(near_interface_id_data[package_index][i][j][k]) > 1)
            {
                mesh_find_if3d<-1, 2>(
                    [&](int l, int m, int n) -> bool
                    {
                        using NeighbourIndex = std::pair<size_t, Arrayi>; /**< stores shifted neighbour info: (size_t)package index, (arrayi)local grid index. */
                        NeighbourIndex neighbour_index = mesh_data_.NeighbourIndexShift(Arrayi(i + l, j + m, k + n), neighborhood);
                        int near_interface_id = near_interface_id_data[neighbour_index.first][neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]];
                        bool is_found = abs(near_interface_id) == 1;
                        if (is_found)
                        {
                            Real phi_0 = phi_data[package_index][i][j][k];
                            near_interface_id_data[package_index][i][j][k] = near_interface_id;
                            phi_data[package_index][i][j][k] = near_interface_id == 1 ? fabs(phi_0) : -fabs(phi_0);
                        }
                        return is_found;
                    });
            }
        });
}
//=============================================================================================//
} // namespace SPH
//=============================================================================================//