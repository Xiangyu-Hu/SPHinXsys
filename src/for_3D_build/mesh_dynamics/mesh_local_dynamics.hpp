#ifndef MESH_LOCAL_DYNAMICS_3D_HPP
#define MESH_LOCAL_DYNAMICS_3D_HPP

#include "grid_data_package_function.hpp"
#include "mesh_local_dynamics.h"

namespace SPH
{
//=================================================================================================//
inline void NearInterfaceCellTagging::UpdateKernel::update(const UnsignedInt &package_index)
{
    UnsignedInt sort_index = data_mesh_->getOccupiedDataPackages()[package_index - num_singular_pkgs_].first;
    Arrayi cell_index = base_dynamics->CellIndexFromSortIndex(sort_index);
    UnsignedInt index_1d = data_mesh_->LinearCellIndex(cell_index);

    MeshVariableData<Real> &grid_phi = phi_[package_index];
    Real phi0 = grid_phi[0][0][0];
    cell_contain_id_[index_1d] = phi0 > 0.0 ? 1 : -1;
    bool is_sign_changed = mesh_any_of3d<0, pkg_size>(
        [&](int i, int j, int k) -> bool
        {
            return grid_phi[i][j][k] * phi0 < 0.0;
        });
    if (is_sign_changed)
        cell_contain_id_[index_1d] = 0;
}
//=================================================================================================//
inline void CellContainDiffusion::UpdateKernel::update(const Arrayi &cell_index)
{
    UnsignedInt index_1d = data_mesh_->LinearCellIndex(cell_index);
    if (cell_contain_id_[index_1d] == 2)
    {
        if (mesh_any_of(
                Arrayi::Zero().max(cell_index - Arrayi::Ones()),
                data_mesh_->AllCells().min(cell_index + 2 * Arrayi::Ones()),
                [&](int l, int m, int n)
                {
                    UnsignedInt neighbor_1d = data_mesh_->LinearCellIndex(Arrayi(l, m, n));
                    return cell_contain_id_[neighbor_1d] == -1;
                }))
        {
            cell_contain_id_[index_1d] = -1;
            cell_package_index_[index_1d] = 0; // inside far field package updated
            AtomicRef<UnsignedInt> count_modified_cells(*count_modified_);
            ++count_modified_cells;
        }
        else if (mesh_any_of(
                     Arrayi::Zero().max(cell_index - Arrayi::Ones()),
                     data_mesh_->AllCells().min(cell_index + 2 * Arrayi::Ones()),
                     [&](int l, int m, int n)
                     {
                         UnsignedInt neighbor_1d = data_mesh_->LinearCellIndex(Arrayi(l, m, n));
                         return cell_contain_id_[neighbor_1d] == 1;
                     }))
        {
            cell_contain_id_[index_1d] = 1;
            cell_package_index_[index_1d] = 1; // outside far field package updated
            AtomicRef<UnsignedInt> count_modified_cells(*count_modified_);
            ++count_modified_cells;
        }
    }
}
//=============================================================================================//
inline void UpdateLevelSetGradient::UpdateKernel::update(const UnsignedInt &package_index)
{
    auto &neighborhood = cell_neighborhood_[package_index];
    auto &pkg_data = phi_gradient_[package_index];

    mesh_for_each3d<0, pkg_size>(
        [&](int i, int j, int k)
        {
            PackageGridPair x1 = NeighbourIndexShift<pkg_size>(
                Arrayi(i + 1, j, k), neighborhood);
            PackageGridPair x2 = NeighbourIndexShift<pkg_size>(
                Arrayi(i - 1, j, k), neighborhood);
            PackageGridPair y1 = NeighbourIndexShift<pkg_size>(
                Arrayi(i, j + 1, k), neighborhood);
            PackageGridPair y2 = NeighbourIndexShift<pkg_size>(
                Arrayi(i, j - 1, k), neighborhood);
            PackageGridPair z1 = NeighbourIndexShift<pkg_size>(
                Arrayi(i, j, k + 1), neighborhood);
            PackageGridPair z2 = NeighbourIndexShift<pkg_size>(
                Arrayi(i, j, k - 1), neighborhood);
            Real dphidx = phi_[x1.first][x1.second[0]][x1.second[1]][x1.second[2]] -
                          phi_[x2.first][x2.second[0]][x2.second[1]][x2.second[2]];
            Real dphidy = phi_[y1.first][y1.second[0]][y1.second[1]][y1.second[2]] -
                          phi_[y2.first][y2.second[0]][y2.second[1]][y2.second[2]];
            Real dphidz = phi_[z1.first][z1.second[0]][z1.second[1]][z1.second[2]] -
                          phi_[z2.first][z2.second[0]][z2.second[1]][z2.second[2]];
            pkg_data[i][j][k] = 0.5 * Vecd(dphidx, dphidy, dphidz) / data_spacing_;
        });
}
//=============================================================================================//
inline void UpdateKernelIntegrals::initializeSingularPackages(UnsignedInt package_index, Real far_field_level_set)
{
    auto &kernel_weight = mv_kernel_weight_.Data()[package_index];
    auto &kernel_gradient = mv_kernel_gradient_.Data()[package_index];
    auto &kernel_second_gradient = mv_kernel_second_gradient_.Data()[package_index];

    mesh_for_each3d<0, pkg_size>(
        [&](int i, int j, int k)
        {
            kernel_weight[i][j][k] = far_field_level_set < 0.0 ? 0 : 1.0;
            kernel_gradient[i][j][k] = Vec3d::Zero();
            kernel_second_gradient[i][j][k] = Mat3d::Zero();
        });
}
//=============================================================================================//
template <typename DataType, typename FunctionByGrid>
void UpdateKernelIntegrals::UpdateKernel::assignByGrid(
    MeshVariableData<DataType> *mesh_variable, const Arrayi &cell_index, const FunctionByGrid &function_by_grid)
{
    UnsignedInt package_index = index_handler_->PackageIndexFromCellIndex(cell_pkg_index_, cell_index);
    auto &pkg_data = mesh_variable[package_index];
    mesh_for_each3d<0, pkg_size>(
        [&](int i, int j, int k)
        {
            pkg_data[i][j][k] = function_by_grid(Array3i(i, j, k));
        });
}
//=============================================================================================//
inline Real UpdateKernelIntegrals::UpdateKernel::
    computeKernelIntegral(const UnsignedInt &package_index, const Arrayi &grid_index)
{
    Real phi = phi_[package_index][grid_index[0]][grid_index[1]][grid_index[2]];

    Real integral(0);
    if (fabs(phi) < cutoff_radius_)
    {
        mesh_for_each_neighbor3d(
            depth_,
            [&](int i, int j, int k)
            {
                PackageGridPair neighbor_meta = GeneralNeighbourIndexShift<pkg_size>(
                    package_index, cell_neighborhood_, grid_index + Arrayi(i, j, k));
                Real phi_neighbor = phi_[neighbor_meta.first]
                                        [neighbor_meta.second[0]]
                                        [neighbor_meta.second[1]]
                                        [neighbor_meta.second[2]];
                if (phi_neighbor > -data_spacing_)
                {
                    Vecd phi_gradient = phi_gradient_[neighbor_meta.first]
                                                     [neighbor_meta.second[0]]
                                                     [neighbor_meta.second[1]]
                                                     [neighbor_meta.second[2]];
                    Vecd displacement = -Arrayi(i, j, k).cast<Real>().matrix() * data_spacing_;
                    Real distance = displacement.norm();
                    if (distance < cutoff_radius_)
                        integral += kernel_.W(displacement) *
                                    CutCellVolumeFraction(phi_neighbor, phi_gradient, data_spacing_);
                }
            });
    }
    return phi > cutoff_radius_ ? 1.0 : integral * data_spacing_ * data_spacing_ * data_spacing_;
}
//=============================================================================================//
inline Vecd UpdateKernelIntegrals::UpdateKernel::
    computeKernelGradientIntegral(const UnsignedInt &package_index, const Arrayi &grid_index)
{
    Real phi = phi_[package_index][grid_index[0]][grid_index[1]][grid_index[2]];

    Vecd integral = Vecd::Zero();
    if (fabs(phi) < cutoff_radius_)
    {
        mesh_for_each_neighbor3d(
            depth_,
            [&](int i, int j, int k)
            {
                PackageGridPair neighbor_meta = GeneralNeighbourIndexShift<pkg_size>(
                    package_index, cell_neighborhood_, grid_index + Arrayi(i, j, k));
                Real phi_neighbor = phi_[neighbor_meta.first]
                                        [neighbor_meta.second[0]]
                                        [neighbor_meta.second[1]]
                                        [neighbor_meta.second[2]];
                if (phi_neighbor > -data_spacing_)
                {
                    Vecd phi_gradient = phi_gradient_[neighbor_meta.first]
                                                     [neighbor_meta.second[0]]
                                                     [neighbor_meta.second[1]]
                                                     [neighbor_meta.second[2]];
                    Vecd displacement = -Arrayi(i, j, k).cast<Real>().matrix() * data_spacing_;
                    Real distance = displacement.norm();
                    if (distance < cutoff_radius_)
                        integral += kernel_.dW(displacement) *
                                    CutCellVolumeFraction(phi_neighbor, phi_gradient, data_spacing_) *
                                    displacement / (distance + TinyReal);
                }
            });
    }

    return integral * data_spacing_ * data_spacing_ * data_spacing_;
}
//=============================================================================================//
inline Matd UpdateKernelIntegrals::UpdateKernel::
    computeKernelSecondGradientIntegral(const UnsignedInt &package_index, const Arrayi &grid_index)
{
    Real phi = phi_[package_index][grid_index[0]][grid_index[1]][grid_index[2]];

    Matd integral = Matd::Zero();
    if (fabs(phi) < cutoff_radius_)
    {
        mesh_for_each_neighbor3d(
            depth_,
            [&](int i, int j, int k)
            {
                PackageGridPair neighbor_meta = GeneralNeighbourIndexShift<pkg_size>(
                    package_index, cell_neighborhood_, grid_index + Arrayi(i, j, k));
                Real phi_neighbor = phi_[neighbor_meta.first]
                                        [neighbor_meta.second[0]]
                                        [neighbor_meta.second[1]]
                                        [neighbor_meta.second[2]];
                if (phi_neighbor > -data_spacing_)
                {
                    Vecd phi_gradient = phi_gradient_[neighbor_meta.first]
                                                     [neighbor_meta.second[0]]
                                                     [neighbor_meta.second[1]]
                                                     [neighbor_meta.second[2]];
                    Vecd displacement = -Arrayi(i, j, k).cast<Real>().matrix() * data_spacing_;
                    Real distance = displacement.norm();
                    if (distance < cutoff_radius_)
                        integral += kernel_.d2W(displacement) *
                                    CutCellVolumeFraction(phi_neighbor, phi_gradient, data_spacing_) *
                                    displacement * displacement.transpose() / (distance * distance + TinyReal);
                }
            });
    }
    return integral * data_spacing_ * data_spacing_ * data_spacing_;
}
//=============================================================================================//
inline void ReinitializeLevelSet::UpdateKernel::update(const UnsignedInt &package_index)
{
    auto &phi_addrs = phi_[package_index];
    auto &near_interface_id_addrs = near_interface_id_[package_index];
    auto &neighborhood = cell_neighborhood_[package_index];

    mesh_for_each3d<0, pkg_size>(
        [&](int i, int j, int k)
        {
            // only reinitialize non cut cells
            if (near_interface_id_addrs[i][j][k] != 0)
            {
                Real phi_0 = phi_addrs[i][j][k];
                Real sign = phi_0 / sqrt(phi_0 * phi_0 + data_spacing_ * data_spacing_);
                PackageGridPair x1 = NeighbourIndexShift<pkg_size>(
                    Arrayi(i + 1, j, k), neighborhood);
                PackageGridPair x2 = NeighbourIndexShift<pkg_size>(
                    Arrayi(i - 1, j, k), neighborhood);
                PackageGridPair y1 = NeighbourIndexShift<pkg_size>(
                    Arrayi(i, j + 1, k), neighborhood);
                PackageGridPair y2 = NeighbourIndexShift<pkg_size>(
                    Arrayi(i, j - 1, k), neighborhood);
                PackageGridPair z1 = NeighbourIndexShift<pkg_size>(
                    Arrayi(i, j, k + 1), neighborhood);
                PackageGridPair z2 = NeighbourIndexShift<pkg_size>(
                    Arrayi(i, j, k - 1), neighborhood);
                Real dv_x = upwindDifference(
                    sign,
                    phi_[x1.first][x1.second[0]][x1.second[1]][x1.second[2]] - phi_0,
                    phi_0 - phi_[x2.first][x2.second[0]][x2.second[1]][x2.second[2]]);
                Real dv_y = upwindDifference(
                    sign,
                    phi_[y1.first][y1.second[0]][y1.second[1]][y1.second[2]] - phi_0,
                    phi_0 - phi_[y2.first][y2.second[0]][y2.second[1]][y2.second[2]]);
                Real dv_z = upwindDifference(
                    sign,
                    phi_[z1.first][z1.second[0]][z1.second[1]][z1.second[2]] - phi_0,
                    phi_0 - phi_[z2.first][z2.second[0]][z2.second[1]][z2.second[2]]);
                phi_addrs[i][j][k] -= 0.3 * sign * (Vec3d(dv_x, dv_y, dv_z).norm() - data_spacing_);
            }
        });
}
//=============================================================================================//
inline void MarkCutInterfaces::UpdateKernel::update(const UnsignedInt &package_index, Real dt)
{
    auto &phi_addrs = phi_[package_index];
    auto &near_interface_id_addrs = near_interface_id_[package_index];

    // corner averages, note that the first row and first column are not used
    PackageDataMatrix<Real, 5> corner_averages;
    mesh_for_each3d<0, 5>(
        [&](int i, int j, int k)
        {
            corner_averages[i][j][k] = CornerAverage(phi_, Arrayi(i, j, k), Arrayi(-1, -1, -1),
                                                     cell_neighborhood_[package_index], (Real)0);
        });

    mesh_for_each3d<0, pkg_size>(
        [&](int i, int j, int k)
        {
            // first assume far cells
            Real phi_0 = phi_addrs[i][j][k];
            int near_interface_id = phi_0 > 0.0 ? 2 : -2;
            if (fabs(phi_0) < perturbation_)
            {
                near_interface_id = 0;
                Real phi_average_0 = corner_averages[i][j][k];
                // find inner and outer cut cells
                mesh_for_each3d<0, 2>(
                    [&](int l, int m, int n)
                    {
                        Real phi_average = corner_averages[i + l][j + m][k + n];
                        if ((phi_average_0 - perturbation_) * (phi_average - perturbation_) < 0.0)
                            near_interface_id = 1;
                        if ((phi_average_0 + perturbation_) * (phi_average + perturbation_) < 0.0)
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
inline void MarkNearInterface::UpdateKernel::update(const UnsignedInt &package_index, Real dt)
{
    mesh_for_each3d<0, pkg_size>(
        [&](int i, int j, int k)
        {
            near_interface_id_[package_index][i][j][k] = 3; // undetermined
            Real phi0 = phi_[package_index][i][j][k];
            if (ABS(phi0) < 2.0 * threshold_) // only consider data close to the interface
            {
                bool is_sign_changed = mesh_any_of3d<-1, 2>( // check in the 3x3x3 neighborhood
                    [&](int l, int m, int n) -> bool
                    {
                        PackageGridPair neighbour_index = NeighbourIndexShift<pkg_size>(
                            Arrayi(i + l, j + m, k + n), cell_neighborhood_[package_index]);

                        return phi0 * phi_[neighbour_index.first]
                                          [neighbour_index.second[0]]
                                          [neighbour_index.second[1]]
                                          [neighbour_index.second[2]] <
                               0.0;
                    });

                if (is_sign_changed)
                {
                    if (ABS(phi0) < threshold_)
                    {
                        near_interface_id_[package_index][i][j][k] = 0; // cut cell
                    }
                }
                else
                {
                    near_interface_id_[package_index][i][j][k] = phi0 > 0.0 ? 1 : -1; // in the band
                }
            }
        });
}
//=============================================================================================//
inline void RedistanceInterface::UpdateKernel::update(const UnsignedInt &package_index)
{
    mesh_for_each3d<0, pkg_size>(
        [&](int i, int j, int k)
        {
            int near_interface_id = near_interface_id_[package_index][i][j][k];
            if (near_interface_id == 0)
            {
                bool positive_band = false;
                bool negative_band = false;
                mesh_for_each3d<-1, 2>(
                    [&](int r, int s, int t)
                    {
                        PackageGridPair neighbour_index = NeighbourIndexShift<pkg_size>(
                            Arrayi(i + r, j + s, k + t), cell_neighborhood_[package_index]);
                        int neighbor_near_interface_id = near_interface_id_[neighbour_index.first]
                                                                           [neighbour_index.second[0]]
                                                                           [neighbour_index.second[1]]
                                                                           [neighbour_index.second[2]];
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
                            PackageGridPair neighbour_index = NeighbourIndexShift<pkg_size>(
                                Arrayi(i + x, j + y, k + z), cell_neighborhood_[package_index]);
                            auto &neighbor_phi = phi_[neighbour_index.first];
                            auto &neighbor_phi_gradient = phi_gradient_[neighbour_index.first];
                            auto &neighbor_near_interface_id = near_interface_id_[neighbour_index.first];
                            if (neighbor_near_interface_id[neighbour_index.second[0]]
                                                          [neighbour_index.second[1]]
                                                          [neighbour_index.second[2]] >= 1)
                            {
                                Real phi_p_ = neighbor_phi[neighbour_index.second[0]]
                                                          [neighbour_index.second[1]]
                                                          [neighbour_index.second[2]];
                                Vecd norm_to_face = neighbor_phi_gradient[neighbour_index.second[0]]
                                                                         [neighbour_index.second[1]]
                                                                         [neighbour_index.second[2]];
                                norm_to_face /= norm_to_face.norm() + TinyReal;
                                min_distance_p = SMIN(
                                    min_distance_p,
                                    (Vecd((Real)x, (Real)y, Real(z)) * data_spacing_ + phi_p_ * norm_to_face).norm());
                            }
                        });
                    phi_[package_index][i][j][k] = -min_distance_p;
                    // this immediate switch of near interface id
                    // does not intervening with the identification of unresolved interface
                    // based on the assumption that positive false_and negative bands are not close to each other
                    near_interface_id_[package_index][i][j][k] = -1;
                }
                if (negative_band == false)
                {
                    Real min_distance_n = 5.0 * data_spacing_;
                    mesh_for_each3d<-4, 5>(
                        [&](int x, int y, int z)
                        {
                            PackageGridPair neighbour_index = NeighbourIndexShift<pkg_size>(
                                Arrayi(i + x, j + y, k + z), cell_neighborhood_[package_index]);
                            auto &neighbor_phi = phi_[neighbour_index.first];
                            auto &neighbor_phi_gradient = phi_gradient_[neighbour_index.first];
                            auto &neighbor_near_interface_id = near_interface_id_[neighbour_index.first];
                            if (neighbor_near_interface_id[neighbour_index.second[0]]
                                                          [neighbour_index.second[1]]
                                                          [neighbour_index.second[2]] <= -1)
                            {
                                Real phi_n_ = neighbor_phi[neighbour_index.second[0]]
                                                          [neighbour_index.second[1]]
                                                          [neighbour_index.second[2]];
                                Vecd norm_to_face = neighbor_phi_gradient[neighbour_index.second[0]]
                                                                         [neighbour_index.second[1]]
                                                                         [neighbour_index.second[2]];
                                norm_to_face /= norm_to_face.norm() + TinyReal;
                                min_distance_n = SMIN(
                                    min_distance_n,
                                    (Vecd((Real)x, (Real)y, Real(z)) * data_spacing_ - phi_n_ * norm_to_face).norm());
                            }
                        });
                    phi_[package_index][i][j][k] = min_distance_n;
                    // this immediate switch of near interface id
                    // does not intervening with the identification of unresolved interface
                    // based on the assumption that positive false_and negative bands are not close to each other
                    near_interface_id_[package_index][i][j][k] = 1;
                }
            }
        });
}
//=============================================================================================//
inline void DiffuseLevelSetSign::UpdateKernel::update(const UnsignedInt &package_index)
{
    mesh_for_each3d<0, pkg_size>(
        [&](int i, int j, int k)
        {
            if (near_interface_id_[package_index][i][j][k] == 3) // check undetermined only
            {
                Real phi = phi_[package_index][i][j][k];
                if (mesh_any_of3d<-1, 2>(
                        [&](int l, int m, int n) -> bool
                        {
                            PackageGridPair neighbour_index = NeighbourIndexShift<pkg_size>(
                                Arrayi(i + l, j + m, k + n), cell_neighborhood_[package_index]);
                            return near_interface_id_[neighbour_index.first]
                                                     [neighbour_index.second[0]]
                                                     [neighbour_index.second[1]]
                                                     [neighbour_index.second[2]] == 1;
                        }))
                {
                    phi_[package_index][i][j][k] = ABS(phi);
                    near_interface_id_[package_index][i][j][k] = 1; // mark as positive band
                    AtomicRef<UnsignedInt> count_modified_data(*count_modified_);
                    ++count_modified_data;
                }
                else if (mesh_any_of3d<-1, 2>(
                             [&](int l, int m, int n) -> bool
                             {
                                 PackageGridPair neighbour_index = NeighbourIndexShift<pkg_size>(
                                     Arrayi(i + l, j + m, k + n), cell_neighborhood_[package_index]);
                                 return near_interface_id_[neighbour_index.first]
                                                          [neighbour_index.second[0]]
                                                          [neighbour_index.second[1]]
                                                          [neighbour_index.second[2]] == -1;
                             }))
                {
                    phi_[package_index][i][j][k] = -ABS(phi);
                    near_interface_id_[package_index][i][j][k] = -1; // mark as negative band
                    AtomicRef<UnsignedInt> count_modified_data(*count_modified_);
                    ++count_modified_data;
                }
            }
        });
}
//=============================================================================================//
} // namespace SPH
//=============================================================================================//
#endif // MESH_LOCAL_DYNAMICS_3D_HPP