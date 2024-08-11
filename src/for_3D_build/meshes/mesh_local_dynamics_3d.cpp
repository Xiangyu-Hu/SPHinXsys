#include "mesh_local_dynamics.h"

namespace SPH
{
//=============================================================================================//
bool TagACellIsInnerPackage::isInnerPackage(const Arrayi &cell_index)
{
    return mesh_any_of(
        Array3i::Zero().max(cell_index - Array3i::Ones()),
        all_cells_.min(cell_index + 2 * Array3i::Ones()),
        [&](int l, int m, int n)
        {
            return mesh_data_.isInnerDataPackage(Arrayi(l, m, n));    //actually a core test here, because only core pkgs are assigned
        });
}
//=============================================================================================//
void InitializeCellNeighborhood::update(const size_t &package_index)
{
    size_t sort_index = mesh_data_.occupied_data_pkgs_[package_index-2].first;
    Arrayi cell_index = Arrayi(sort_index / all_cells_[1], sort_index % all_cells_[1]); //[notion] there might be problems, 3d implementation needed
    CellNeighborhood &current = mesh_data_.cell_neighborhood_[package_index];
    std::pair<Arrayi, int> &metadata = mesh_data_.meta_data_cell_[package_index];
    metadata.first = cell_index;
    metadata.second = mesh_data_.occupied_data_pkgs_[package_index-2].second;
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
    // auto &phi = phi_.DataField()[package_index];
    // auto &near_interface_id = near_interface_id_.DataField()[package_index];
    // Arrayi cell_index = mesh_data_.CellIndexFromPackageSortIndex(package_index);
    // mesh_data_.for_each_cell_data(
    //     [&](int i, int j, int k)
    //     {
    //         Vec3d position = mesh_data_.DataPositionFromIndex(cell_index, Array3i(i, j, k));
    //         phi[i][j][k] = shape_.findSignedDistance(position);
    //         near_interface_id[i][j][k] = phi[i][j][k] < 0.0 ? -2 : 2;
    //     });
    printf("not implemented yet");
}
//=============================================================================================//
Arrayi InitializeBasicDataForAPackage::CellIndexFromSortIndex(const size_t &sort_index)
{
    return Array3i(sort_index / all_cells_[1], sort_index % all_cells_[1], sort_index % all_cells_[1]); //[notion] not implemented yet
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

    return integral * data_spacing_ * data_spacing_;
}
//=============================================================================================//
void DiffuseLevelSetSign::update(const size_t &package_index)
{
    printf("this is the execution of DiffuseLevelSetSign\n");
}
//=============================================================================================//
} // namespace SPH
//=============================================================================================//