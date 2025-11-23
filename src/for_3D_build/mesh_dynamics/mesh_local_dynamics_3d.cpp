#include "mesh_local_dynamics.hpp"

namespace SPH
{
//=============================================================================================//
void InitializeBasicPackageData::initializeSingularPackages(
    const UnsignedInt package_index, Real far_field_level_set)
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
        data_mesh_->AllCells().min(cell_index + 2 * Array3i::Ones()),
        [&](int l, int m, int n)
        {
            return data_mesh_->isInnerDataPackage(Arrayi(l, m, n)); // actually a core test here, because only core pkgs are assigned
        });
}
//=============================================================================================//
void InitializeBasicPackageData::UpdateKernel::update(const UnsignedInt &package_index)
{
    auto &phi = phi_[package_index];
    auto &near_interface_id = near_interface_id_[package_index];
    Arrayi cell_index = index_handler_.DimensionalCellIndex(pkg_1d_cell_index_[package_index]);
    mesh_for_each3d<0, pkg_size>(
        [&](int i, int j, int k)
        {
            Vec3d position = index_handler_.DataPositionFromIndex(cell_index, Array3i(i, j, k));
            phi[i][j][k] = shape_->findSignedDistance(position);
            near_interface_id[i][j][k] = phi[i][j][k] < 0.0 ? -2 : 2;
        });
}
//=============================================================================================//
} // namespace SPH
