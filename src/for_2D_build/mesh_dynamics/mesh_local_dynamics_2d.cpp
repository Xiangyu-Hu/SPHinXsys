#include "mesh_iterators.hpp"
#include "mesh_local_dynamics.hpp"

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
void InitializeBasicPackageData::initializeSingularPackages(const size_t package_index, Real far_field_level_set)
{
    auto &phi = phi_.Data()[package_index];
    auto &near_interface_id = near_interface_id_.Data()[package_index];
    auto &phi_gradient = phi_gradient_.Data()[package_index];

    mesh_for_each2d<0, pkg_size>(
        [&](int i, int j)
        {
            phi[i][j] = far_field_level_set;
            near_interface_id[i][j] = far_field_level_set < 0.0 ? -2 : 2;
            phi_gradient[i][j] = Vecd::Ones();
        });
}
//=============================================================================================//
bool InnerCellTagging::UpdateKernel::isInnerPackage(const Arrayi &cell_index)
{
    return mesh_any_of(
        Array2i::Zero().max(cell_index - Array2i::Ones()),
        all_cells_.min(cell_index + 2 * Array2i::Ones()),
        [&](int l, int m)
        {
            return data_mesh_->isInnerDataPackage(Arrayi(l, m)); // actually a core test here, because only core pkgs are assigned
        });
}
//=============================================================================================//
void InitializeBasicPackageData::UpdateKernel::update(const size_t &package_index)
{
    auto &phi = phi_[package_index];
    auto &near_interface_id = near_interface_id_[package_index];
    Arrayi cell_index = meta_data_cell_[package_index].first;
    mesh_for_each2d<0, pkg_size>(
        [&](int i, int j)
        {
            Vec2d position = index_handler_->DataPositionFromIndex(cell_index, Array2i(i, j));
            phi[i][j] = shape_->findSignedDistance(position);
            near_interface_id[i][j] = phi[i][j] < 0.0 ? -2 : 2;
        });
}
//=============================================================================================//
} // namespace SPH
