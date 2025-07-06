#include "mesh_local_dynamics.hpp"
#include "mesh_iterators.hpp"

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
void InitializeDataForSingularPackage::update(const size_t package_index, Real far_field_level_set)
{
    auto &phi = phi_.Data()[package_index];
    auto &near_interface_id = near_interface_id_.Data()[package_index];
    auto &phi_gradient = phi_gradient_.Data()[package_index];
    auto &kernel_weight = kernel_weight_.Data()[package_index];
    auto &kernel_gradient = kernel_gradient_.Data()[package_index];

    BaseMeshLocalDynamics::for_each_cell_data(
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
bool TagACellIsInnerPackage::UpdateKernel::isInnerPackage(const Arrayi &cell_index)
{
    return mesh_any_of(
        Array2i::Zero().max(cell_index - Array2i::Ones()),
        all_cells_.min(cell_index + 2 * Array2i::Ones()),
        [&](int l, int m)
        {
            return mesh_data_->isInnerDataPackage(Arrayi(l, m));    //actually a core test here, because only core pkgs are assigned
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
        {
            current[l + 1][m + 1] = mesh_data_->PackageIndexFromCellIndex(cell_package_index_, cell_index + Arrayi(l, m));
        }
}
//=============================================================================================//
void InitializeBasicDataForAPackage::UpdateKernel::update(const size_t &package_index)
{
    auto &phi = phi_[package_index];
    auto &near_interface_id = near_interface_id_[package_index];
    Arrayi cell_index = meta_data_cell_[package_index].first;
    BaseMeshLocalDynamics::for_each_cell_data(
        [&](int i, int j)
        {
            Vec2d position = index_handler_->DataPositionFromIndex(cell_index, Array2i(i, j));
            phi[i][j] = shape_->findSignedDistance(position);
            near_interface_id[i][j] = phi[i][j] < 0.0 ? -2 : 2;
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
                << "  DATAPACKING=POINT \n";

    mesh_for_column_major(
        Arrayi::Zero(), number_of_operation,
        [&](const Array2i &global_index)
        {
            Vecd data_position = mesh_data_.global_mesh_.GridPositionFromIndex(global_index);
            output_file << data_position[0] << " ";
            output_file << data_position[1] << " ";
            output_file << BaseMeshLocalDynamics::DataValueFromGlobalIndex(phi_.Data(), global_index, &mesh_data_, cell_package_index_.Data()) << " ";
            output_file << BaseMeshLocalDynamics::DataValueFromGlobalIndex(phi_gradient_.Data(), global_index, &mesh_data_, cell_package_index_.Data())[0] << " ";
            output_file << BaseMeshLocalDynamics::DataValueFromGlobalIndex(phi_gradient_.Data(), global_index, &mesh_data_, cell_package_index_.Data())[1] << " ";
            output_file << BaseMeshLocalDynamics::DataValueFromGlobalIndex(near_interface_id_.Data(), global_index, &mesh_data_, cell_package_index_.Data()) << " ";
            output_file << BaseMeshLocalDynamics::DataValueFromGlobalIndex(kernel_weight_.Data(), global_index, &mesh_data_, cell_package_index_.Data()) << " ";
            output_file << BaseMeshLocalDynamics::DataValueFromGlobalIndex(kernel_gradient_.Data(), global_index, &mesh_data_, cell_package_index_.Data())[0] << " ";
            output_file << BaseMeshLocalDynamics::DataValueFromGlobalIndex(kernel_gradient_.Data(), global_index, &mesh_data_, cell_package_index_.Data())[1] << " ";
            output_file << " \n";
        });
    output_file << " \n";
}
//=============================================================================================//
} // namespace SPH
