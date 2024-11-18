#include "mesh_local_dynamics.hpp"

#include "mesh_iterators.h"

namespace SPH
{
//=============================================================================================//
size_t BaseMeshLocalDynamics::SortIndexFromCellIndex(const Arrayi &cell_index)
{
    return cell_index[0] * all_cells_[1] + cell_index[1];
}
//=============================================================================================//
Arrayi BaseMeshLocalDynamics::CellIndexFromSortIndex(const size_t &sort_index)
{
    Array2i cell_index;
    cell_index[0] = sort_index / all_cells_[1];
    cell_index[1] = sort_index % all_cells_[1];
    return cell_index;
}
//=============================================================================================//
std::pair<size_t, Arrayi> BaseMeshLocalDynamics::
    NeighbourIndexShift(const Arrayi shift_index, const CellNeighborhood &neighbour)
{
    std::pair<size_t, Arrayi> result;
    Arrayi neighbour_index = (shift_index + pkg_size * Arrayi::Ones()) / pkg_size;
    result.first = neighbour[neighbour_index[0]][neighbour_index[1]];
    result.second = (shift_index + pkg_size * Arrayi::Ones()) - neighbour_index * pkg_size;

    return result;
}
//=============================================================================================//
void InitializeDataForSingularPackage::update(const size_t package_index, Real far_field_level_set)
{
    auto &phi = phi_.DataField()[package_index];
    auto &near_interface_id = near_interface_id_.DataField()[package_index];
    auto &phi_gradient = phi_gradient_.DataField()[package_index];
    auto &kernel_weight = kernel_weight_.DataField()[package_index];
    auto &kernel_gradient = kernel_gradient_.DataField()[package_index];

    for_each_cell_data(
        [&](int i, int j)
        {
            phi[i][j] = far_field_level_set;
            near_interface_id[i][j] = far_field_level_set < 0.0 ? -2 : 2;
            phi_gradient[i][j] = Vecd::Ones();
            kernel_weight[i][j] = far_field_level_set < 0.0 ? 0 : 1.0;
            kernel_gradient[i][j] = Vec2d::Zero();
        });
}
//=============================================================================================//
bool TagACellIsInnerPackage::isInnerPackage(const Arrayi &cell_index)
{
    return mesh_any_of(
        Array2i::Zero().max(cell_index - Array2i::Ones()),
        all_cells_.min(cell_index + 2 * Array2i::Ones()),
        [&](int l, int m)
        {
            return mesh_data_.isInnerDataPackage(Arrayi(l, m));    //actually a core test here, because only core pkgs are assigned
        });
}
//=============================================================================================//
void InitializeCellNeighborhood::update(const size_t &package_index)
{
    size_t sort_index = mesh_data_.occupied_data_pkgs_[package_index-2].first;
    Arrayi cell_index = Arrayi(sort_index / all_cells_[1], sort_index % all_cells_[1]); //[notion] there might be problems, 3d implementation needed
    CellNeighborhood &current = mesh_data_.cell_neighborhood_[package_index];
    std::pair<Arrayi, int> &metadata = meta_data_cell_[package_index];
    metadata.first = cell_index;
    metadata.second = mesh_data_.occupied_data_pkgs_[package_index-2].second;
    for (int l = -1; l < 2; l++)
        for (int m = -1; m < 2; m++)
        {
            current[l + 1][m + 1] = mesh_data_.PackageIndexFromCellIndex(cell_index + Arrayi(l, m));
        }
}
//=============================================================================================//
void InitializeBasicDataForAPackage::update(const size_t &package_index)
{
    auto &phi = phi_.DataField()[package_index];
    auto &near_interface_id = near_interface_id_.DataField()[package_index];
    Arrayi cell_index = meta_data_cell_[package_index].first;
    for_each_cell_data(
        [&](int i, int j)
        {
            Vec2d position = mesh_data_.DataPositionFromIndex(cell_index, Array2i(i, j));
            phi[i][j] = shape_.findSignedDistance(position);
            near_interface_id[i][j] = phi[i][j] < 0.0 ? -2 : 2;
        });
}
//=============================================================================================//
void InitializeBasicDataForAPackage::UpdateKernel::update(const size_t &package_index)
{
    auto &phi = phi_[package_index];
    auto &near_interface_id = near_interface_id_[package_index];
    Arrayi cell_index = meta_data_cell_[package_index].first;
    base_dynamics->for_each_cell_data(
        [&](int i, int j)
        {
            Vec2d position = mesh_data_->DataPositionFromIndex(cell_index, Array2i(i, j));
            phi[i][j] = shape_->findSignedDistance(position);
            near_interface_id[i][j] = phi[i][j] < 0.0 ? -2 : 2;
        });
}
//=============================================================================================//
void UpdateLevelSetGradient::UpdateKernel::update(const size_t &package_index)
{
    auto &neighborhood = cell_neighborhood_[package_index];
    auto &pkg_data = phi_gradient_[package_index];

    base_dynamics->for_each_cell_data(
        [&](int i, int j)
        {
            std::pair<size_t, Arrayi> x1 = base_dynamics->NeighbourIndexShift(Arrayi(i + 1, j), neighborhood);
            std::pair<size_t, Arrayi> x2 = base_dynamics->NeighbourIndexShift(Arrayi(i - 1, j), neighborhood);
            std::pair<size_t, Arrayi> y1 = base_dynamics->NeighbourIndexShift(Arrayi(i, j + 1), neighborhood);
            std::pair<size_t, Arrayi> y2 = base_dynamics->NeighbourIndexShift(Arrayi(i, j - 1), neighborhood);
            Real dphidx = (phi_[x1.first][x1.second[0]][x1.second[1]] -
                           phi_[x2.first][x2.second[0]][x2.second[1]]);
            Real dphidy = (phi_[y1.first][y1.second[0]][y1.second[1]] -
                           phi_[y2.first][y2.second[0]][y2.second[1]]);

            pkg_data[i][j] = 0.5 * Vecd(dphidx, dphidy) / data_spacing_;
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
        mesh_for_each2d<-3, 4>(
            [&](int i, int j)
            {
                Arrayi neighbor_index = Arrayi(global_index_[0] + i, global_index_[1] + j);
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
    return phi > threshold ? 1.0 : integral * data_spacing_ * data_spacing_;
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
        mesh_for_each2d<-3, 4>(
            [&](int i, int j)
            {
                Arrayi neighbor_index = Arrayi(global_index_[0] + i, global_index_[1] + j);
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

    return integral * data_spacing_ * data_spacing_;
}
//=============================================================================================//
void ReinitializeLevelSet::update(const size_t &package_index)
{
    auto phi_data = phi_.DataField();
    auto &phi_addrs = phi_data[package_index];
    auto &near_interface_id_addrs = near_interface_id_.DataField()[package_index];
    auto &neighborhood = mesh_data_.cell_neighborhood_[package_index];

    for_each_cell_data(
        [&](int i, int j)
        {
            // only reinitialize non cut cells
            if (near_interface_id_addrs[i][j] != 0)
            {
                Real phi_0 = phi_addrs[i][j];
                Real sign = phi_0 / sqrt(phi_0 * phi_0 + data_spacing_ * data_spacing_);
                using NeighbourIndex = std::pair<size_t, Arrayi>; /**< stores shifted neighbour info: (size_t)package index, (arrayi)local grid index. */
                NeighbourIndex x1 = NeighbourIndexShift(Arrayi(i + 1, j), neighborhood);
                NeighbourIndex x2 = NeighbourIndexShift(Arrayi(i - 1, j), neighborhood);
                NeighbourIndex y1 = NeighbourIndexShift(Arrayi(i, j + 1), neighborhood);
                NeighbourIndex y2 = NeighbourIndexShift(Arrayi(i, j - 1), neighborhood);
                Real dv_x = upwindDifference(sign, phi_data[x1.first][x1.second[0]][x1.second[1]] - phi_0,
                                              phi_0 - phi_data[x2.first][x2.second[0]][x2.second[1]]);
                Real dv_y = upwindDifference(sign, phi_data[y1.first][y1.second[0]][y1.second[1]] - phi_0,
                                              phi_0 - phi_data[y2.first][y2.second[0]][y2.second[1]]);
                phi_addrs[i][j] -= 0.5 * sign * (Vec2d(dv_x, dv_y).norm() - data_spacing_);
            }
        });
}
//=============================================================================================//
void MarkNearInterface::update(const size_t &package_index)
{
    auto &phi_addrs = phi_.DataField()[package_index];
    auto &near_interface_id_addrs = near_interface_id_.DataField()[package_index];
    auto neighborhood = mesh_data_.cell_neighborhood_[package_index];

    // corner averages, note that the first row and first column are not used
    // template <typename DataType>
    // using PackageTemporaryData = PackageDataMatrix<DataType, 5>; //[notion] the pkg_size should be replaced
    PackageDataMatrix<Real, 5> corner_averages;
    mesh_for_each2d<0, 5>( //[notion] same pkg_size problem above
        [&](int i, int j)
        {
            corner_averages[i][j] = mesh_data_.CornerAverage(phi_, Arrayi(i, j), Arrayi(-1, -1), neighborhood);
        });

    for_each_cell_data(
        [&](int i, int j)
        {
            // first assume far cells
            Real phi_0 = phi_addrs[i][j];
            int near_interface_id = phi_0 > 0.0 ? 2 : -2;
            if (fabs(phi_0) < small_shift)
            {
                near_interface_id = 0;
                Real phi_average_0 = corner_averages[i][j];
                // find outer cut cells by comparing the sign of corner averages
                mesh_for_each2d<0, 2>(
                    [&](int l, int m)
                    {
                        Real phi_average = corner_averages[i + l][j + m];
                        if ((phi_average_0 - small_shift) * (phi_average - small_shift) < 0.0)
                            near_interface_id = 1;
                        if ((phi_average_0 + small_shift) * (phi_average + small_shift) < 0.0)
                            near_interface_id = -1;
                    });
                // find zero cut cells by comparing the sign of corner averages
                mesh_for_each2d<0, 2>(
                    [&](int l, int m)
                    {
                        Real phi_average = corner_averages[i + l][j + m];
                        if (phi_average_0 * phi_average < 0.0)
                            near_interface_id = 0;
                    });
            }
            // assign this to package
            near_interface_id_addrs[i][j] = near_interface_id;
        });
}
//=============================================================================================//
void RedistanceInterface::update(const size_t &package_index)
{
    auto phi_data = phi_.DataField();
    auto near_interface_id_data = near_interface_id_.DataField();
    auto &neighborhood = mesh_data_.cell_neighborhood_[package_index];
    using NeighbourIndex = std::pair<size_t, Arrayi>; /**< stores shifted neighbour info: (size_t)package index, (arrayi)local grid index. */

    for_each_cell_data(
        [&](int i, int j)
        {
            int near_interface_id = near_interface_id_data[package_index][i][j];
            if (near_interface_id == 0)
            {
                bool positive_band = false;
                bool negative_band = false;
                mesh_for_each2d<-1, 2>(
                    [&](int r, int s)
                    {
                        NeighbourIndex neighbour_index = NeighbourIndexShift(Arrayi(i + r, j + s), neighborhood);
                        int neighbor_near_interface_id = near_interface_id_data[neighbour_index.first][neighbour_index.second[0]][neighbour_index.second[1]];
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
                            NeighbourIndex neighbour_index = NeighbourIndexShift(Arrayi(i + x, j + y), neighborhood);
                            auto &neighbor_phi = phi_.DataField()[neighbour_index.first];
                            auto &neighbor_phi_gradient = phi_gradient_.DataField()[neighbour_index.first];
                            auto &neighbor_near_interface_id = near_interface_id_.DataField()[neighbour_index.first];
                            if (neighbor_near_interface_id[neighbour_index.second[0]][neighbour_index.second[1]] >= 1)
                            {
                                Real phi_p_ = neighbor_phi[neighbour_index.second[0]][neighbour_index.second[1]];
                                Vecd norm_to_face = neighbor_phi_gradient[neighbour_index.second[0]][neighbour_index.second[1]];
                                norm_to_face /= norm_to_face.norm() + TinyReal;
                                min_distance_p = SMIN(min_distance_p, (Vecd((Real)x, (Real)y) * data_spacing_ + phi_p_ * norm_to_face).norm());
                            }
                        });
                    phi_data[package_index][i][j] = -min_distance_p;
                    // this immediate switch of near interface id
                    // does not intervening with the identification of unresolved interface
                    // based on the assumption that positive false_and negative bands are not close to each other
                    near_interface_id_data[package_index][i][j] = -1;
                }
                if (negative_band == false)
                {
                    Real min_distance_n = 5.0 * data_spacing_;
                    mesh_for_each2d<-4, 5>(
                        [&](int x, int y)
                        {
                            NeighbourIndex neighbour_index = NeighbourIndexShift(Arrayi(i + x, j + y), neighborhood);
                            auto &neighbor_phi = phi_.DataField()[neighbour_index.first];
                            auto &neighbor_phi_gradient = phi_gradient_.DataField()[neighbour_index.first];
                            auto &neighbor_near_interface_id = near_interface_id_.DataField()[neighbour_index.first];
                            if (neighbor_near_interface_id[neighbour_index.second[0]][neighbour_index.second[1]] <= -1)
                            {
                                Real phi_n_ = neighbor_phi[neighbour_index.second[0]][neighbour_index.second[1]];
                                Vecd norm_to_face = neighbor_phi_gradient[neighbour_index.second[0]][neighbour_index.second[1]];
                                norm_to_face /= norm_to_face.norm() + TinyReal;
                                min_distance_n = SMIN(min_distance_n, (Vecd((Real)x, (Real)y) * data_spacing_ - phi_n_ * norm_to_face).norm());
                            }
                        });
                    phi_data[package_index][i][j] = min_distance_n;
                    // this immediate switch of near interface id
                    // does not intervening with the identification of unresolved interface
                    // based on the assumption that positive false_and negative bands are not close to each other
                    near_interface_id_data[package_index][i][j] = 1;
                }
            }
        });
}
//=============================================================================================//
void DiffuseLevelSetSign::update(const size_t &package_index)
{
    auto phi_data = phi_.DataField();
    auto near_interface_id_data = near_interface_id_.DataField();
    auto &neighborhood = mesh_data_.cell_neighborhood_[package_index];

    for_each_cell_data(
        [&](int i, int j)
        {
            // near interface cells are not considered
            if (abs(near_interface_id_data[package_index][i][j]) > 1)
            {
                mesh_find_if2d<-1, 2>(
                    [&](int l, int m) -> bool
                    {
                        using NeighbourIndex = std::pair<size_t, Arrayi>; /**< stores shifted neighbour info: (size_t)package index, (arrayi)local grid index. */
                        NeighbourIndex neighbour_index = NeighbourIndexShift(Arrayi(i + l, j + m), neighborhood);
                        int near_interface_id = near_interface_id_data[neighbour_index.first][neighbour_index.second[0]][neighbour_index.second[1]];
                        bool is_found = abs(near_interface_id) == 1;
                        if (is_found)
                        {
                            Real phi_0 = phi_data[package_index][i][j];
                            near_interface_id_data[package_index][i][j] = near_interface_id;
                            phi_data[package_index][i][j] = near_interface_id == 1 ? fabs(phi_0) : -fabs(phi_0);
                        }
                        return is_found;
                    });
            }
        });
}
//=============================================================================================//
void WriteMeshFieldToPlt::update(std::ofstream &output_file)
{
    Arrayi number_of_operation = mesh_data_.global_mesh_.AllGridPoints();

    output_file << "\n";
    output_file << "title='View'"
                << "\n";
    output_file << "variables= "
                << "x, "
                << "y, "
                << "phi, "
                << "n_x, "
                << "n_y "
                << "near_interface_id ";
    output_file << "kernel_weight, "
                << "kernel_gradient_x, "
                << "kernel_gradient_y "
                << "\n";
    output_file << "zone i=" << number_of_operation[0] << "  j=" << number_of_operation[1] << "  k=" << 1
                << "  DATAPACKING=BLOCK  SOLUTIONTIME=" << 0 << "\n";

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            Vecd data_position = mesh_data_.global_mesh_.GridPositionFromIndex(Arrayi(i, j));
            output_file << data_position[0] << " ";
        }
        output_file << " \n";
    }

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            Vecd data_position = mesh_data_.global_mesh_.GridPositionFromIndex(Arrayi(i, j));
            output_file << data_position[1] << " ";
        }
        output_file << " \n";
    }

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            output_file << mesh_data_.DataValueFromGlobalIndex(phi_, Arrayi(i, j))
                        << " ";
        }
        output_file << " \n";
    }

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            output_file << mesh_data_.DataValueFromGlobalIndex(phi_gradient_, Arrayi(i, j))[0]
                        << " ";
        }
        output_file << " \n";
    }

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            output_file << mesh_data_.DataValueFromGlobalIndex(phi_gradient_, Arrayi(i, j))[1]
                        << " ";
        }
        output_file << " \n";
    }

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            output_file << mesh_data_.DataValueFromGlobalIndex(near_interface_id_, Arrayi(i, j))
                        << " ";
        }
        output_file << " \n";
    }

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            output_file << mesh_data_.DataValueFromGlobalIndex(kernel_weight_, Arrayi(i, j))
                        << " ";
        }
        output_file << " \n";
    }

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            output_file << mesh_data_.DataValueFromGlobalIndex(kernel_gradient_, Arrayi(i, j))[0]
                        << " ";
        }
        output_file << " \n";
    }

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            output_file << mesh_data_.DataValueFromGlobalIndex(kernel_gradient_, Arrayi(i, j))[1]
                        << " ";
        }
        output_file << " \n";
    }
}
//=============================================================================================//
} // namespace SPH
//=============================================================================================//