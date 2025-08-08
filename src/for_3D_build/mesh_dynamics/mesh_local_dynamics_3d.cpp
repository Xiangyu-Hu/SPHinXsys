#include "mesh_local_dynamics.hpp"

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
void InitializeBasicPackageData::initializeSingularPackages(
    const size_t package_index, Real far_field_level_set)
{
    auto &phi = mv_phi_.Data()[package_index];
    auto &near_interface_id = mv_near_interface_id_.Data()[package_index];
    auto &phi_gradient = mv_phi_gradient_.Data()[package_index];

    mesh_for_each3d<0, pkg_size>(
        [&](int i, int j, int k)
        {
            phi[i][j][k] = far_field_level_set;
            near_interface_id[i][j][k] = far_field_level_set < 0.0 ? -2 : 2;
            phi_gradient[i][j][k] = Vecd::Ones();
        });
}
//=============================================================================================//
bool InnerCellTagging::UpdateKernel::isInnerPackage(const Arrayi &cell_index)
{
    return mesh_any_of(
        Array3i::Zero().max(cell_index - Array3i::Ones()),
        all_cells_.min(cell_index + 2 * Array3i::Ones()),
        [&](int l, int m, int n)
        {
            return data_mesh_->isInnerDataPackage(Arrayi(l, m, n)); // actually a core test here, because only core pkgs are assigned
        });
}
//=============================================================================================//
void InitializeBasicPackageData::UpdateKernel::update(const size_t &package_index)
{
    auto &phi = phi_[package_index];
    auto &near_interface_id = near_interface_id_[package_index];
    Arrayi cell_index = meta_data_cell_[package_index].first;
    mesh_for_each3d<0, pkg_size>(
        [&](int i, int j, int k)
        {
            Vec3d position = index_handler_->DataPositionFromIndex(cell_index, Array3i(i, j, k));
            phi[i][j][k] = shape_->findSignedDistance(position);
            near_interface_id[i][j][k] = phi[i][j][k] < 0.0 ? -2 : 2;
        });
}
//=============================================================================================//
} // namespace SPH
