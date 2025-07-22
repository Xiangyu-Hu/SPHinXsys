#include "mesh_local_dynamics.hpp"
#include "write_mesh_field_to_plt.h"
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

    mesh_for_each3d<0, pkg_size>(
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
bool TagACellIsInnerPackage::UpdateKernel::isInnerPackage(const Arrayi &cell_index)
{
    return mesh_any_of(
        Array3i::Zero().max(cell_index - Array3i::Ones()),
        all_cells_.min(cell_index + 2 * Array3i::Ones()),
        [&](int l, int m, int n)
        {
            return mesh_data_->isInnerDataPackage(Arrayi(l, m, n));    //actually a core test here, because only core pkgs are assigned
        });
}
//=============================================================================================//
void InitializeCellNeighborhood::UpdateKernel::update(const size_t &package_index)
{
    size_t sort_index = mesh_data_->occupied_data_pkgs_[package_index-2].first;
    Arrayi cell_index = base_dynamics->CellIndexFromSortIndex(sort_index);
    CellNeighborhood &current = cell_neighborhood_[package_index];
    std::pair<Arrayi, int> &metadata = meta_data_cell_[package_index];
    metadata.first = cell_index;
    metadata.second = mesh_data_->occupied_data_pkgs_[package_index-2].second;
    for (int l = -1; l < 2; l++)
        for (int m = -1; m < 2; m++)
            for (int n = -1; n < 2; n++)
            {
                current[l + 1][m + 1][n + 1] = mesh_data_->PackageIndexFromCellIndex(cell_package_index_, cell_index + Arrayi(l, m, n));
            }
}
//=============================================================================================//
void InitializeBasicDataForAPackage::UpdateKernel::update(const size_t &package_index)
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
void WriteMeshFieldToPlt::update(std::ofstream &output_file)
{
    StdVec<Coord3D> active_cells;
    mesh_for_each(Array3i::Zero(), mesh_data_.AllCells(), [&](const Arrayi &cell_index)
                  { 
      if (mesh_data_.isCoreDataPackage(cell_index))
      {
        active_cells.push_back({cell_index[0], cell_index[1], cell_index[2]});
      } });

    StdVec<Block3D> clustered_blocks = clusterActiveCells3D(active_cells);

    output_file << "\n";
    output_file << "title='View'"
                << "\n";
    output_file << "variables= "
                << "x, "
                << "y, "
                << "z, "
                << "phi, "
                << "n_x, "
                << "n_y, "
                << "n_z, "
                << "near_interface_id "
                << "\n";
    for (size_t l = 0; l != clustered_blocks.size(); ++l)
    {
        auto &block = clustered_blocks[l];
        Array3i lower_bound_cell_index = Array3i(block.first[0], block.first[1], block.first[2]);
        Array3i block_size = Array3i(block.second[0] - block.first[0] + 1,
                                     block.second[1] - block.first[1] + 1,
                                     block.second[2] - block.first[2] + 1);
        Array3i global_lower_bound = mesh_data_.DataPackageSize() * lower_bound_cell_index - Array3i::Ones();
        Array3i global_upper_bound = mesh_data_.DataPackageSize() * (lower_bound_cell_index + block_size) + Array3i::Ones();
        output_file << "zone i=" << block_size[0] * mesh_data_.DataPackageSize() + 2 << "  j="
                    << block_size[1] * mesh_data_.DataPackageSize() + 2 << "  k="
                    << block_size[2] * mesh_data_.DataPackageSize() + 2
                    << "  DATAPACKING=POINT \n";

        mesh_for_column_major(
            global_lower_bound, global_upper_bound,
            [&](const Array3i &global_index)
            {
                Vecd data_position = mesh_data_.global_mesh_.GridPositionFromIndex(global_index);
                output_file << data_position[0] << " ";
                output_file << data_position[1] << " ";
                output_file << data_position[2] << " ";
                output_file << DataValueFromGlobalIndex(phi_.Data(), global_index, &mesh_data_, cell_package_index_.Data()) << " ";
                output_file << DataValueFromGlobalIndex(phi_gradient_.Data(), global_index, &mesh_data_, cell_package_index_.Data())[0] << " ";
                output_file << DataValueFromGlobalIndex(phi_gradient_.Data(), global_index, &mesh_data_, cell_package_index_.Data())[1] << " ";
                output_file << DataValueFromGlobalIndex(phi_gradient_.Data(), global_index, &mesh_data_, cell_package_index_.Data())[2] << " ";
                output_file << DataValueFromGlobalIndex(near_interface_id_.Data(), global_index, &mesh_data_, cell_package_index_.Data()) << " ";
                output_file << " \n";
            });
        output_file << " \n";
    }
}
//=============================================================================================//
} // namespace SPH
//=============================================================================================//