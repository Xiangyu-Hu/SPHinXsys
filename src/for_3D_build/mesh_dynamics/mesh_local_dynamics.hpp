#ifndef MESH_LOCAL_DYNAMICS_3D_HPP
#define MESH_LOCAL_DYNAMICS_3D_HPP

#include "grid_data_package_function.hpp"
#include "mesh_local_dynamics.h"

namespace SPH
{
//=============================================================================================//
inline void UpdateLevelSetGradient::UpdateKernel::update(const size_t &package_index)
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
template <class KernelType>
template <typename DataType, typename FunctionByGrid>
void UpdateKernelIntegrals<KernelType>::UpdateKernel::
    assignByGrid(MeshVariableData<DataType> *mesh_variable,
                 const Arrayi &cell_index,
                 const FunctionByGrid &function_by_grid)
{
    size_t package_index = index_handler_->PackageIndexFromCellIndex(cell_package_index_, cell_index);
    auto &pkg_data = mesh_variable[package_index];
    mesh_for_each3d<0, pkg_size>(
        [&](int i, int j, int k)
        {
            pkg_data[i][j][k] = function_by_grid(Array3i(i, j, k));
        });
}
//=============================================================================================//
template <class KernelType>
Real UpdateKernelIntegrals<KernelType>::UpdateKernel::
    computeKernelIntegral(const size_t &package_index, const Arrayi &grid_index)
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
                        integral += kernel_->W(global_h_ratio_, distance, displacement) *
                                    CutCellVolumeFraction(phi_neighbor, phi_gradient, data_spacing_);
                }
            });
    }
    return phi > cutoff_radius_ ? 1.0 : integral * data_spacing_ * data_spacing_ * data_spacing_;
}
//=============================================================================================//
template <class KernelType>
Vecd UpdateKernelIntegrals<KernelType>::UpdateKernel::
    computeKernelGradientIntegral(const size_t &package_index, const Arrayi &grid_index)
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
Matd UpdateKernelIntegrals<KernelType>::UpdateKernel::
    computeKernelSecondGradientIntegral(const size_t &package_index, const Arrayi &grid_index)
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
inline void MarkNearInterface::UpdateKernel::update(const size_t &package_index, Real small_shift_factor)
{
    Real small_shift = data_spacing_ * small_shift_factor;
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
inline void RedistanceInterface::UpdateKernel::update(const size_t &package_index)
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
inline void DiffuseLevelSetSign::UpdateKernel::update(const size_t &package_index)
{
    mesh_for_each3d<0, pkg_size>(
        [&](int i, int j, int k)
        {
            // near interface cells are not considered
            if (abs(near_interface_id_[package_index][i][j][k]) > 1)
            {
                mesh_find_if3d<-1, 2>(
                    [&](int l, int m, int n) -> bool
                    {
                        PackageGridPair neighbour_index = NeighbourIndexShift<pkg_size>(
                            Arrayi(i + l, j + m, k + n), cell_neighborhood_[package_index]);
                        int near_interface_id = near_interface_id_[neighbour_index.first]
                                                                  [neighbour_index.second[0]]
                                                                  [neighbour_index.second[1]]
                                                                  [neighbour_index.second[2]];
                        bool is_found = abs(near_interface_id) == 1;
                        if (is_found)
                        {
                            Real phi_0 = phi_[package_index][i][j][k];
                            near_interface_id_[package_index][i][j][k] = near_interface_id;
                            phi_[package_index][i][j][k] = near_interface_id == 1 ? fabs(phi_0) : -fabs(phi_0);
                        }
                        return is_found;
                    });
            }
        });
}
//=============================================================================================//
} // namespace SPH
//=============================================================================================//
#endif // MESH_LOCAL_DYNAMICS_3D_HPP