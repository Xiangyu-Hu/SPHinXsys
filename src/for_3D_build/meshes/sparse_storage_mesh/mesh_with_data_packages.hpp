#ifndef MESH_WITH_DATA_PACKAGES_3D_HPP
#define MESH_WITH_DATA_PACKAGES_3D_HPP

#include "mesh_with_data_packages.h"

#include "block_cluster_3d.h"

namespace SPH
{
//=================================================================================================//
template <size_t PKG_SIZE>
void MeshWithGridDataPackages<PKG_SIZE>::writeMeshVariableToPlt(std::ofstream &output_file)
{
    StdVec<Coord3D> active_cells;
    auto meta_data = meta_data_cell_.Data();
    package_for(execution::seq, num_singular_pkgs_, num_grid_pkgs_,
                [&](size_t package_index)
                {
                    if (meta_data[package_index].second == 1)
                    {
                        auto cell_index = meta_data[package_index].first;
                        active_cells.push_back({cell_index[0], cell_index[1], cell_index[2]});
                    }
                });
    StdVec<Block3D> clustered_blocks = clusterActiveCells3D(active_cells);

    output_file << "\n"
                << "title='View'" << "\n";
    output_file << " VARIABLES = " << "x, " << "y, " << "z";

    constexpr int type_index_int = DataTypeIndex<int>::value;
    for (MeshVariable<int> *variable : std::get<type_index_int>(mesh_variable_to_write_))
    {
        output_file << ",\"" << variable->Name() << "\"";
    };

    constexpr int type_index_Vecd = DataTypeIndex<Vec3d>::value;
    for (MeshVariable<Vec3d> *variable : std::get<type_index_Vecd>(mesh_variable_to_write_))
    {
        std::string variable_name = variable->Name();
        output_file << ",\"" << variable_name << "_x\""
                    << ",\"" << variable_name << "_y\""
                    << ",\"" << variable_name << "_z\"";
    };

    constexpr int type_index_Real = DataTypeIndex<Real>::value;
    for (MeshVariable<Real> *variable : std::get<type_index_Real>(mesh_variable_to_write_))
    {
        output_file << ",\"" << variable->Name() << "\"";
    };

    output_file << "\n";

    // Write the zone header
    for (size_t l = 0; l != clustered_blocks.size(); ++l)
    {
        auto &block = clustered_blocks[l];
        Array3i lower_bound_cell_index = Array3i(block.first[0], block.first[1], block.first[2]);
        Array3i block_size = Array3i(block.second[0] - block.first[0] + 1,
                                     block.second[1] - block.first[1] + 1,
                                     block.second[2] - block.first[2] + 1);
        Array3i global_lower_bound = DataPackageSize() * lower_bound_cell_index - Array3i::Ones();
        Array3i global_upper_bound = DataPackageSize() * (lower_bound_cell_index + block_size) + Array3i::Ones();
        output_file << "zone i=" << block_size[0] * DataPackageSize() + 2 << "  j="
                    << block_size[1] * DataPackageSize() + 2 << "  k="
                    << block_size[2] * DataPackageSize() + 2
                    << "  DATAPACKING=POINT \n";

        mesh_for_column_major(
            global_lower_bound, global_upper_bound,
            [&](const Array3i &global_index)
            {
                Vec3d data_position = global_mesh_.GridPositionFromIndex(global_index);
                output_file << data_position[0] << " " << data_position[1] << " " << data_position[2] << " ";

                constexpr int type_index_int = DataTypeIndex<int>::value;
                for (MeshVariable<int> *variable : std::get<type_index_int>(mesh_variable_to_write_))
                {
                    int value = DataValueFromGlobalIndex(variable->Data(), global_index, this, cell_package_index_.Data());
                    output_file << value << " ";
                };

                constexpr int type_index_Vecd = DataTypeIndex<Vec3d>::value;
                for (MeshVariable<Vec3d> *variable : std::get<type_index_Vecd>(mesh_variable_to_write_))
                {
                    Vec3d value = DataValueFromGlobalIndex(variable->Data(), global_index, this, cell_package_index_.Data());
                    output_file << value[0] << " " << value[1] << " " << value[2] << " ";
                };

                constexpr int type_index_Real = DataTypeIndex<Real>::value;
                for (MeshVariable<Real> *variable : std::get<type_index_Real>(mesh_variable_to_write_))
                {
                    Real value = DataValueFromGlobalIndex(variable->Data(), global_index, this, cell_package_index_.Data());
                    output_file << value << " ";
                };
                output_file << " \n";
            });
        output_file << " \n";
    }
}
//=================================================================================================//
} // namespace SPH
#endif // MESH_WITH_DATA_PACKAGES_3D_HPP
