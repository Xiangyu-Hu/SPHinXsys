#ifndef MESH_LOCAL_DYNAMICS_2D_HPP
#define MESH_LOCAL_DYNAMICS_2D_HPP

#include "grid_data_package_function.hpp"
#include "mesh_local_dynamics.h"

namespace SPH
{
//=================================================================================================//
inline void ConsistencyCorrection::UpdateKernel::update(const size_t &package_index)
{
    MeshVariableData<Real> &grid_phi = phi_[package_index];
    mesh_for_each2d<0, pkg_size_minus1>(
        [&](int i, int j)
        {
            Real phi0 = grid_phi[i][j];
            Real phi0_abs = ABS(phi0);
            mesh_for_each2d<0, 2>(
                [&](int l, int m)
                {
                    Real phi1 = grid_phi[i + l][j + m];
                    Real phi1_abs = ABS(phi1);
                    if ( // inconsistency criterion: level set must be continuous
                        phi0 * phi1 < 0.0 && phi0_abs + phi1_abs > threshold_)
                    {
                        Real sgn = phi0 < phi1 ? SGN(phi0) : SGN(phi1);
                        grid_phi[i][j] = sgn * phi0_abs;
                        grid_phi[i + l][j + m] = sgn * phi1_abs;
                    }
                });
        });
}
//=================================================================================================//
inline void NearSurfaceCellContainTagging::UpdateKernel::update(const size_t &package_index)
{
    size_t sort_index = data_mesh_->occupied_data_pkgs_[package_index - num_singular_pkgs_].first;
    Arrayi cell_index = base_dynamics->CellIndexFromSortIndex(sort_index);
    UnsignedInt index_1d = data_mesh_->transferMeshIndexTo1D(data_mesh_->AllCells(), cell_index);

    MeshVariableData<Real> &grid_phi = phi_[package_index];
    Real phi0 = grid_phi[0][0];
    cell_near_interface_id_[index_1d] = phi0 > 0.0 ? 1 : -1;
    bool is_sign_changed = mesh_any_of2d<0, pkg_size>(
        [&](int i, int j) -> bool
        {
            return grid_phi[i][j] * phi0 < 0.0;
        });
    if (is_sign_changed)
        cell_near_interface_id_[index_1d] = 0;
}
//=================================================================================================//
inline void CellContainDiffusion::UpdateKernel::update(const Arrayi &cell_index)
{
    UnsignedInt index_1d = data_mesh_->transferMeshIndexTo1D(data_mesh_->AllCells(), cell_index);
    if (cell_near_interface_id_[index_1d] == 2)
    {
        if (mesh_any_of(
                Arrayi::Zero().max(cell_index - Arrayi::Ones()),
                data_mesh_->AllCells().min(cell_index + 2 * Arrayi::Ones()),
                [&](int l, int m)
                {
                    UnsignedInt neighbor_1d = data_mesh_->transferMeshIndexTo1D(data_mesh_->AllCells(), Arrayi(l, m));
                    return cell_near_interface_id_[neighbor_1d] == -1;
                }))
        {
            cell_near_interface_id_[index_1d] = -1;
            AtomicRef<UnsignedInt> count_modified_cells(*count_modified_);
            ++count_modified_cells;
        }
        else if (mesh_any_of(
                     Arrayi::Zero().max(cell_index - Arrayi::Ones()),
                     data_mesh_->AllCells().min(cell_index + 2 * Arrayi::Ones()),
                     [&](int l, int m)
                     {
                         UnsignedInt neighbor_1d = data_mesh_->transferMeshIndexTo1D(data_mesh_->AllCells(), Arrayi(l, m));
                         return cell_near_interface_id_[neighbor_1d] == 1;
                     }))
        {
            cell_near_interface_id_[index_1d] = 1;
            cell_package_index_[index_1d] = 1; // outside far field package updated
            AtomicRef<UnsignedInt> count_modified_cells(*count_modified_);
            ++count_modified_cells;
        }
    }
}
//=============================================================================================//
inline void UpdateLevelSetGradient::UpdateKernel::update(const size_t &package_index)
{
    auto &neighborhood = cell_neighborhood_[package_index];
    auto &pkg_data = phi_gradient_[package_index];

    mesh_for_each2d<0, pkg_size>(
        [&](int i, int j)
        {
            PackageGridPair x1 = NeighbourIndexShift<pkg_size>(
                Arrayi(i + 1, j), neighborhood);
            PackageGridPair x2 = NeighbourIndexShift<pkg_size>(
                Arrayi(i - 1, j), neighborhood);
            PackageGridPair y1 = NeighbourIndexShift<pkg_size>(
                Arrayi(i, j + 1), neighborhood);
            PackageGridPair y2 = NeighbourIndexShift<pkg_size>(
                Arrayi(i, j - 1), neighborhood);
            Real dphidx = phi_[x1.first][x1.second[0]][x1.second[1]] - phi_[x2.first][x2.second[0]][x2.second[1]];
            Real dphidy = phi_[y1.first][y1.second[0]][y1.second[1]] - phi_[y2.first][y2.second[0]][y2.second[1]];
            pkg_data[i][j] = 0.5 * Vecd(dphidx, dphidy) / data_spacing_;
        });
}
//=============================================================================================//
inline void UpdateKernelIntegrals::initializeSingularPackages(
    size_t package_index, Real far_field_level_set)
{
    auto &kernel_weight = kernel_weight_.Data()[package_index];
    auto &kernel_gradient = kernel_gradient_.Data()[package_index];
    auto &kernel_second_gradient = kernel_second_gradient_.Data()[package_index];

    mesh_for_each2d<0, pkg_size>(
        [&](int i, int j)
        {
            kernel_weight[i][j] = far_field_level_set < 0.0 ? 0 : 1.0;
            kernel_gradient[i][j] = Vec2d::Zero();
            kernel_second_gradient[i][j] = Mat2d::Zero();
        });
}
//=============================================================================================//
template <typename DataType, typename FunctionByGrid>
void UpdateKernelIntegrals::UpdateKernel::assignByGrid(MeshVariableData<DataType> *mesh_variable,
                                                       const Arrayi &cell_index,
                                                       const FunctionByGrid &function_by_grid)
{
    size_t package_index = index_handler_->PackageIndexFromCellIndex(cell_package_index_, cell_index);
    auto &pkg_data = mesh_variable[package_index];
    mesh_for_each2d<0, pkg_size>(
        [&](int i, int j)
        {
            pkg_data[i][j] = function_by_grid(Array2i(i, j));
        });
}
//=============================================================================================//
inline Real UpdateKernelIntegrals::UpdateKernel::
    computeKernelIntegral(const size_t &package_index, const Arrayi &grid_index)
{
    Real phi = phi_[package_index][grid_index[0]][grid_index[1]];

    Real integral(0);
    if (fabs(phi) < cutoff_radius_)
    {
        mesh_for_each_neighbor2d(
            depth_,
            [&](int i, int j)
            {
                PackageGridPair neighbor_meta = GeneralNeighbourIndexShift<pkg_size>(
                    package_index, cell_neighborhood_, grid_index + Arrayi(i, j));
                Real phi_neighbor = phi_[neighbor_meta.first][neighbor_meta.second[0]][neighbor_meta.second[1]];
                if (phi_neighbor > -data_spacing_)
                {
                    Vecd phi_gradient = phi_gradient_[neighbor_meta.first][neighbor_meta.second[0]][neighbor_meta.second[1]];
                    Vecd displacement = -Arrayi(i, j).cast<Real>().matrix() * data_spacing_;
                    Real distance = displacement.norm();
                    if (distance < cutoff_radius_)
                        integral += kernel_->W(global_h_ratio_, distance, displacement) *
                                    CutCellVolumeFraction(phi_neighbor, phi_gradient, data_spacing_);
                }
            });
    }
    return phi > cutoff_radius_ ? 1.0 : integral * data_spacing_ * data_spacing_;
}
//=============================================================================================//
inline Vecd UpdateKernelIntegrals::UpdateKernel::
    computeKernelGradientIntegral(const size_t &package_index, const Arrayi &grid_index)
{
    Real phi = phi_[package_index][grid_index[0]][grid_index[1]];

    Vecd integral = Vecd::Zero();
    if (fabs(phi) < cutoff_radius_)
    {
        mesh_for_each_neighbor2d(
            depth_,
            [&](int i, int j)
            {
                PackageGridPair neighbor_meta = GeneralNeighbourIndexShift<pkg_size>(
                    package_index, cell_neighborhood_, grid_index + Arrayi(i, j));
                Real phi_neighbor = phi_[neighbor_meta.first][neighbor_meta.second[0]][neighbor_meta.second[1]];
                if (phi_neighbor > -data_spacing_)
                {
                    Vecd phi_gradient = phi_gradient_[neighbor_meta.first][neighbor_meta.second[0]][neighbor_meta.second[1]];
                    Vecd displacement = -Arrayi(i, j).cast<Real>().matrix() * data_spacing_;
                    Real distance = displacement.norm();
                    if (distance < cutoff_radius_)
                        integral += kernel_->dW(global_h_ratio_, distance, displacement) *
                                    CutCellVolumeFraction(phi_neighbor, phi_gradient, data_spacing_) *
                                    displacement / (distance + TinyReal);
                }
            });
    }

    return integral * data_spacing_ * data_spacing_;
}
//=============================================================================================//
inline Matd UpdateKernelIntegrals::UpdateKernel::
    computeKernelSecondGradientIntegral(const size_t &package_index, const Arrayi &grid_index)
{
    Real phi = phi_[package_index][grid_index[0]][grid_index[1]];

    Matd integral = Matd::Zero();
    if (fabs(phi) < cutoff_radius_)
    {
        mesh_for_each_neighbor2d(
            depth_,
            [&](int i, int j)
            {
                PackageGridPair neighbor_meta = GeneralNeighbourIndexShift<pkg_size>(
                    package_index, cell_neighborhood_, grid_index + Arrayi(i, j));
                Real phi_neighbor = phi_[neighbor_meta.first][neighbor_meta.second[0]][neighbor_meta.second[1]];
                if (phi_neighbor > -data_spacing_)
                {
                    Vecd phi_gradient = phi_gradient_[neighbor_meta.first][neighbor_meta.second[0]][neighbor_meta.second[1]];
                    Vecd displacement = -Arrayi(i, j).cast<Real>().matrix() * data_spacing_;
                    Real distance = displacement.norm();
                    if (distance < cutoff_radius_)
                        integral += kernel_->d2W(global_h_ratio_, distance, displacement) *
                                    CutCellVolumeFraction(phi_neighbor, phi_gradient, data_spacing_) *
                                    displacement * displacement.transpose() / (distance * distance + TinyReal);
                }
            });
    }
    return integral * data_spacing_ * data_spacing_;
}
//=============================================================================================//
inline void ReinitializeLevelSet::UpdateKernel::update(const size_t &package_index)
{
    auto &phi_addrs = phi_[package_index];
    auto &near_interface_id_addrs = near_interface_id_[package_index];
    auto &neighborhood = cell_neighborhood_[package_index];

    mesh_for_each2d<0, pkg_size>(
        [&](int i, int j)
        {
            // only reinitialize non cut cells
            if (near_interface_id_addrs[i][j] != 0)
            {
                Real phi_0 = phi_addrs[i][j];
                Real sign = phi_0 / sqrt(phi_0 * phi_0 + data_spacing_ * data_spacing_);
                PackageGridPair x1 = NeighbourIndexShift<pkg_size>(
                    Arrayi(i + 1, j), neighborhood);
                PackageGridPair x2 = NeighbourIndexShift<pkg_size>(
                    Arrayi(i - 1, j), neighborhood);
                PackageGridPair y1 = NeighbourIndexShift<pkg_size>(
                    Arrayi(i, j + 1), neighborhood);
                PackageGridPair y2 = NeighbourIndexShift<pkg_size>(
                    Arrayi(i, j - 1), neighborhood);
                Real dv_x = upwindDifference(sign, phi_[x1.first][x1.second[0]][x1.second[1]] - phi_0,
                                             phi_0 - phi_[x2.first][x2.second[0]][x2.second[1]]);
                Real dv_y = upwindDifference(sign, phi_[y1.first][y1.second[0]][y1.second[1]] - phi_0,
                                             phi_0 - phi_[y2.first][y2.second[0]][y2.second[1]]);
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
            corner_averages[i][j] = CornerAverage(phi_, Arrayi(i, j), Arrayi(-1, -1),
                                                  cell_neighborhood_[package_index], (Real)0);
        });

    mesh_for_each2d<0, pkg_size>(
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
    mesh_for_each2d<0, pkg_size>(
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
                        PackageGridPair neighbour_index = NeighbourIndexShift<pkg_size>(
                            Arrayi(i + r, j + s), cell_neighborhood_[package_index]);
                        int neighbor_near_interface_id =
                            near_interface_id_[neighbour_index.first][neighbour_index.second[0]][neighbour_index.second[1]];
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
                            PackageGridPair neighbour_index = NeighbourIndexShift<pkg_size>(
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
                            PackageGridPair neighbour_index = NeighbourIndexShift<pkg_size>(
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
    mesh_for_each2d<0, pkg_size>(
        [&](int i, int j)
        {
            // near interface cells are not considered
            if (abs(near_interface_id_[package_index][i][j]) > 1)
            {
                mesh_find_if2d<-1, 2>(
                    [&](int l, int m) -> bool
                    {
                        PackageGridPair neighbour_index =
                            NeighbourIndexShift<pkg_size>(
                                Arrayi(i + l, j + m), cell_neighborhood_[package_index]);
                        int near_interface_id =
                            near_interface_id_[neighbour_index.first][neighbour_index.second[0]][neighbour_index.second[1]];
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
#endif // MESH_LOCAL_DYNAMICS_2D_HPP