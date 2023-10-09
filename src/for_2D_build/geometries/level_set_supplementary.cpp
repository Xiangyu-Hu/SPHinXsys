#include "level_set.h"

#include "base_body.h"
#include "base_kernel.h"
#include "base_particle_dynamics.h"
#include "base_particles.h"
#include "mesh_iterators.hpp"

namespace SPH
{
//=================================================================================================//
LevelSet::LevelSet(BoundingBox tentative_bounds, Real data_spacing,
                   Shape &shape, SPHAdaptation &sph_adaptation)
    : LevelSet(tentative_bounds, data_spacing, 4, shape, sph_adaptation)
{
    mesh_parallel_for(MeshRange(Arrayi::Zero(), all_cells_),
                      [&](size_t i, size_t j)
                      {
                          initializeDataInACell(Arrayi(i, j));
                      });

    finishDataPackages();
}
//=================================================================================================//
void LevelSet::initializeDataForSingularPackage(LevelSetDataPackage *data_pkg, Real far_field_level_set)
{
    auto &phi = data_pkg->getPackageData(phi_);
    auto &near_interface_id = data_pkg->getPackageData(near_interface_id_);
    auto &phi_gradient = data_pkg->getPackageData(phi_gradient_);
    auto &kernel_weight = data_pkg->getPackageData(kernel_weight_);
    auto &kernel_gradient = data_pkg->getPackageData(kernel_gradient_);

    data_pkg->for_each_data(
        [&](int i, int j)
        {
            phi[i][j] = far_field_level_set;
            near_interface_id[i][j] = far_field_level_set < 0.0 ? -2 : 2;
            phi_gradient[i][j] = Vecd::Ones();
            kernel_weight[i][j] = far_field_level_set < 0.0 ? 0 : 1.0;
            kernel_gradient[i][j] = Vec2d::Zero();
        });
}
//=================================================================================================//
void LevelSet::finishDataPackages()
{
    mesh_parallel_for(MeshRange(Arrayi::Zero(), all_cells_),
                      [&](size_t i, size_t j)
                      {
                          tagACellIsInnerPackage(Arrayi(i, j));
                      });

    mesh_parallel_for(MeshRange(Arrayi::Zero(), all_cells_),
                      [&](size_t i, size_t j)
                      {
                          initializePackageAddressesInACell(Arrayi(i, j));
                      });

    updateLevelSetGradient();
    updateKernelIntegrals();
}
//=================================================================================================//
bool LevelSet::isWithinCorePackage(Vecd position)
{
    Arrayi cell_index = CellIndexFromPosition(position);
    return data_pkg_addrs_[cell_index[0]][cell_index[1]]->isCorePackage();
}
//=============================================================================================//
bool LevelSet::isInnerPackage(const Arrayi &cell_index)
{
    return mesh_any_of(
        Array2i::Zero().max(cell_index - Array2i::Ones()),
        all_cells_.min(cell_index + 2 * Array2i::Ones()),
        [&](int l, int m)
        {
            return data_pkg_addrs_[l][m]->isCorePackage();
        });
}
//=================================================================================================//
void LevelSet::diffuseLevelSetSign()
{
    package_parallel_for(
        inner_data_pkgs_,
        [&](LevelSetDataPackage *data_pkg)
        {
            auto &phi_addrs = data_pkg->getPackageDataAddress(phi_);
            auto &near_interface_id_addrs = data_pkg->getPackageDataAddress(near_interface_id_);

            data_pkg->for_each_addrs(
                [&](int i, int j)
                {
                    // near interface cells are not considered
                    if (abs(*near_interface_id_addrs[i][j]) > 1)
                    {
                        mesh_find_if2d<-1, 2>(
                            [&](int l, int m) -> bool
                            {
                                int near_interface_id = *near_interface_id_addrs[i + l][j + m];
                                bool is_found = abs(near_interface_id) == 1;
                                if (is_found)
                                {
                                    Real phi_0 = *phi_addrs[i][j];
                                    *near_interface_id_addrs[i][j] = near_interface_id;
                                    *phi_addrs[i][j] = near_interface_id == 1 ? fabs(phi_0) : -fabs(phi_0);
                                }
                                return is_found;
                            });
                    }
                });
        });
}
//=============================================================================================//
void LevelSet::reinitializeLevelSet()
{
    package_parallel_for(
        inner_data_pkgs_,
        [&](LevelSetDataPackage *data_pkg)
        {
            auto &phi_addrs = data_pkg->getPackageDataAddress(phi_);
            auto &near_interface_id_addrs = data_pkg->getPackageDataAddress(near_interface_id_);

            data_pkg->for_each_addrs(
                [&](int i, int j)
                {
                    // only reinitialize non cut cells
                    if (*near_interface_id_addrs[i][j] != 0)
                    {
                        Real phi_0 = *phi_addrs[i][j];
                        Real sign = phi_0 / sqrt(phi_0 * phi_0 + data_spacing_ * data_spacing_);
                        Real dv_x = upwindDifference(sign, *phi_addrs[i + 1][j] - phi_0, phi_0 - *phi_addrs[i - 1][j]);
                        Real dv_y = upwindDifference(sign, *phi_addrs[i][j + 1] - phi_0, phi_0 - *phi_addrs[i][j - 1]);
                        *phi_addrs[i][j] -= 0.5 * sign * (Vec2d(dv_x, dv_y).norm() - data_spacing_);
                    }
                });
        });
}
//=================================================================================================//
void LevelSet::markNearInterface(Real small_shift_factor)
{
    Real small_shift = small_shift_factor * data_spacing_;

    package_parallel_for(
        inner_data_pkgs_,
        [&](LevelSetDataPackage *data_pkg)
        {
            auto &phi_addrs = data_pkg->getPackageDataAddress(phi_);
            auto &near_interface_id_addrs = data_pkg->getPackageDataAddress(near_interface_id_);

            // corner averages, note that the first row and first column are not used
            LevelSetDataPackage::PackageTemporaryData<Real> corner_averages;
            mesh_for_each2d<1, pkg_addrs_size>(
                [&](int i, int j)
                {
                    corner_averages[i][j] = data_pkg->CornerAverage(phi_addrs, Arrayi(i, j), Arrayi(-1, -1));
                });

            data_pkg->for_each_addrs(
                [&](int i, int j)
                {
                    // first assume far cells
                    Real phi_0 = *phi_addrs[i][j];
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
                    *near_interface_id_addrs[i][j] = near_interface_id;
                });
        });
}
//=================================================================================================//
void LevelSet::initializeBasicDataForAPackage(LevelSetDataPackage *data_pkg, Shape &shape)
{
    auto &phi = data_pkg->getPackageData(phi_);
    auto &near_interface_id = data_pkg->getPackageData(near_interface_id_);
    data_pkg->for_each_data(
        [&](int i, int j)
        {
            Vec2d position = data_pkg->DataPositionFromIndex(Vec2d(i, j));
            phi[i][j] = shape.findSignedDistance(position);
            near_interface_id[i][j] = phi[i][j] < 0.0 ? -2 : 2;
        });
}
//=================================================================================================//
void LevelSet::redistanceInterfaceForAPackage(LevelSetDataPackage *core_data_pkg)
{
    int l = (int)core_data_pkg->CellIndexOnMesh()[0];
    int m = (int)core_data_pkg->CellIndexOnMesh()[1];
    auto &phi_addrs = core_data_pkg->getPackageDataAddress(phi_);
    auto &near_interface_id_addrs = core_data_pkg->getPackageDataAddress(near_interface_id_);

    core_data_pkg->for_each_addrs(
        [&](int i, int j)
        {
            int near_interface_id = *near_interface_id_addrs[i][j];
            if (near_interface_id == 0)
            {
                bool positive_band = false;
                bool negative_band = false;
                mesh_for_each2d<-1, 2>(
                    [&](int r, int s)
                    {
                        int neighbor_near_interface_id = *near_interface_id_addrs[i + r][j + s];
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
                            std::pair<int, int> x_pair = CellShiftAndDataIndex(i + x);
                            std::pair<int, int> y_pair = CellShiftAndDataIndex(j + y);
                            LevelSetDataPackage *neighbor_pkg = data_pkg_addrs_[l + x_pair.first][m + y_pair.first];
                            auto &neighbor_phi = neighbor_pkg->getPackageData(phi_);
                            auto &neighbor_phi_gradient = neighbor_pkg->getPackageData(phi_gradient_);
                            auto &neighbor_near_interface_id = neighbor_pkg->getPackageData(near_interface_id_);
                            if (neighbor_near_interface_id[x_pair.second][y_pair.second] >= 1)
                            {
                                Real phi_p_ = neighbor_phi[x_pair.second][y_pair.second];
                                Vecd norm_to_face = neighbor_phi_gradient[x_pair.second][y_pair.second];
                                norm_to_face /= norm_to_face.norm() + TinyReal;
                                min_distance_p = SMIN(min_distance_p, (Vecd((Real)x, (Real)y) * data_spacing_ + phi_p_ * norm_to_face).norm());
                            }
                        });
                    *phi_addrs[i][j] = -min_distance_p;
                    // this immediate switch of near interface id
                    // does not intervening with the identification of unresolved interface
                    // based on the assumption that positive false_and negative bands are not close to each other
                    *near_interface_id_addrs[i][j] = -1;
                }
                if (negative_band == false)
                {
                    Real min_distance_n = 5.0 * data_spacing_;
                    mesh_for_each2d<-4, 5>(
                        [&](int x, int y)
                        {
                            std::pair<int, int> x_pair = CellShiftAndDataIndex(i + x);
                            std::pair<int, int> y_pair = CellShiftAndDataIndex(j + y);
                            LevelSetDataPackage *neighbor_pkg = data_pkg_addrs_[l + x_pair.first][m + y_pair.first];
                            auto &neighbor_phi = neighbor_pkg->getPackageData(phi_);
                            auto &neighbor_phi_gradient = neighbor_pkg->getPackageData(phi_gradient_);
                            auto &neighbor_near_interface_id = neighbor_pkg->getPackageData(near_interface_id_);
                            if (neighbor_near_interface_id[x_pair.second][y_pair.second] <= -1)
                            {
                                Real phi_n_ = neighbor_phi[x_pair.second][y_pair.second];
                                Vecd norm_to_face = neighbor_phi_gradient[x_pair.second][y_pair.second];
                                norm_to_face /= norm_to_face.norm() + TinyReal;
                                min_distance_n = SMIN(min_distance_n, (Vecd((Real)x, (Real)y) * data_spacing_ - phi_n_ * norm_to_face).norm());
                            }
                        });
                    *phi_addrs[i][j] = min_distance_n;
                    // this immediate switch of near interface id
                    // does not intervening with the identification of unresolved interface
                    // based on the assumption that positive false_and negative bands are not close to each other
                    *near_interface_id_addrs[i][j] = 1;
                }
            }
        });
}
//=============================================================================================//
void LevelSet::writeMeshFieldToPlt(std::ofstream &output_file)
{
    Arrayi number_of_operation = global_mesh_.AllGridPoints();

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
            Vecd data_position = global_mesh_.GridPositionFromIndex(Arrayi(i, j));
            output_file << data_position[0] << " ";
        }
        output_file << " \n";
    }

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            Vecd data_position = global_mesh_.GridPositionFromIndex(Arrayi(i, j));
            output_file << data_position[1] << " ";
        }
        output_file << " \n";
    }

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            output_file << DataValueFromGlobalIndex(phi_, Arrayi(i, j))
                        << " ";
        }
        output_file << " \n";
    }

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            output_file << DataValueFromGlobalIndex(phi_gradient_, Arrayi(i, j))[0]
                        << " ";
        }
        output_file << " \n";
    }

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            output_file << DataValueFromGlobalIndex(phi_gradient_, Arrayi(i, j))[1]
                        << " ";
        }
        output_file << " \n";
    }

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            output_file << DataValueFromGlobalIndex(near_interface_id_, Arrayi(i, j))
                        << " ";
        }
        output_file << " \n";
    }

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            output_file << DataValueFromGlobalIndex(kernel_weight_, Arrayi(i, j))
                        << " ";
        }
        output_file << " \n";
    }

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            output_file << DataValueFromGlobalIndex(kernel_gradient_, Arrayi(i, j))[0]
                        << " ";
        }
        output_file << " \n";
    }

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            output_file << DataValueFromGlobalIndex(kernel_gradient_, Arrayi(i, j))[1]
                        << " ";
        }
        output_file << " \n";
    }
}
//=============================================================================================//
Real LevelSet::computeKernelIntegral(const Vecd &position)
{
    Real phi = probeSignedDistance(position);
    Real cutoff_radius = kernel_.CutOffRadius(global_h_ratio_);
    Real threshold = cutoff_radius + data_spacing_; // consider that interface's half width is the data spacing

    Real integral(0);
    if (fabs(phi) < threshold)
    {
        Arrayi global_index_ = global_mesh_.CellIndexFromPosition(position);
        mesh_for_each2d<-3, 4>(
            [&](int i, int j)
            {
                Arrayi neighbor_index = Arrayi(global_index_[0] + i, global_index_[1] + j);
                Real phi_neighbor = DataValueFromGlobalIndex(phi_, neighbor_index);
                if (phi_neighbor > -data_spacing_)
                {
                    Vecd phi_gradient = DataValueFromGlobalIndex(phi_gradient_, neighbor_index);
                    Vecd integral_position = global_mesh_.GridPositionFromIndex(neighbor_index);
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
Vecd LevelSet::computeKernelGradientIntegral(const Vecd &position)
{
    Real phi = probeSignedDistance(position);
    Real cutoff_radius = kernel_.CutOffRadius(global_h_ratio_);
    Real threshold = cutoff_radius + data_spacing_;

    Vecd integral = Vecd::Zero();
    if (fabs(phi) < threshold)
    {
        Arrayi global_index_ = global_mesh_.CellIndexFromPosition(position);
        mesh_for_each2d<-3, 4>(
            [&](int i, int j)
            {
                Arrayi neighbor_index = Arrayi(global_index_[0] + i, global_index_[1] + j);
                Real phi_neighbor = DataValueFromGlobalIndex(phi_, neighbor_index);
                if (phi_neighbor > -data_spacing_)
                {
                    Vecd phi_gradient = DataValueFromGlobalIndex(phi_gradient_, neighbor_index);
                    Vecd integral_position = global_mesh_.GridPositionFromIndex(neighbor_index);
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
RefinedLevelSet::RefinedLevelSet(BoundingBox tentative_bounds, LevelSet &coarse_level_set,
                                 Shape &shape, SPHAdaptation &sph_adaptation)
    : RefinedMesh(tentative_bounds, coarse_level_set, 4, shape, sph_adaptation)
{
    mesh_parallel_for(MeshRange(Arrayi::Zero(), all_cells_),
                      [&](size_t i, size_t j)
                      {
                          initializeDataInACellFromCoarse(Arrayi(i, j));
                      });

    finishDataPackages();
}
//=============================================================================================//
} // namespace SPH
//=============================================================================================//
