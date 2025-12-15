#ifndef MESH_WITH_DATA_PACKAGES_3D_HPP
#define MESH_WITH_DATA_PACKAGES_3D_HPP

#include "mesh_with_data_packages.hxx"

#include "block_cluster_3d.h"

namespace SPH
{
//=================================================================================================//
template <int PKG_SIZE>
void MeshWithGridDataPackages<PKG_SIZE>::writeMeshVariablesToPltByMesh(
    UnsignedInt resolution_level, std::ofstream &output_file)
{
    IndexHandler &index_handler = this->getResolutionLevel(resolution_level);

    StdVec<Coord3D> active_cells;
    auto pkg_1d_cell_index = dv_pkg_1d_cell_index_->Data();
    auto pkg_type = dv_pkg_type_->Data();
    package_for(execution::seq, num_pkgs_offsets_[resolution_level], num_pkgs_offsets_[resolution_level + 1],
                [&](UnsignedInt package_index)
                {
                    if (pkg_type[package_index] == 1)
                    {
                        auto cell_index = index_handler.DimensionalCellIndex(pkg_1d_cell_index[package_index]);
                        active_cells.push_back({cell_index[0], cell_index[1], cell_index[2]});
                    }
                });
    StdVec<Block3D> clustered_blocks = clusterActiveCells3D(active_cells);

    Mesh global_mesh = index_handler.getGlobalMesh();

    output_file << "\n"
                << "title='View'" << "\n";
    output_file << " VARIABLES = " << "x, " << "y, " << "z";

    constexpr int type_index_int = DataTypeIndex<int>::value;
    for (MeshVariable<int> *variable : std::get<type_index_int>(mesh_variables_to_write_))
    {
        output_file << ",\"" << variable->Name() << "\"";
    };

    constexpr int type_index_Vecd = DataTypeIndex<Vec3d>::value;
    for (MeshVariable<Vec3d> *variable : std::get<type_index_Vecd>(mesh_variables_to_write_))
    {
        std::string variable_name = variable->Name();
        output_file << ",\"" << variable_name << "_x\""
                    << ",\"" << variable_name << "_y\""
                    << ",\"" << variable_name << "_z\"";
    };

    constexpr int type_index_Real = DataTypeIndex<Real>::value;
    for (MeshVariable<Real> *variable : std::get<type_index_Real>(mesh_variables_to_write_))
    {
        output_file << ",\"" << variable->Name() << "\"";
    };

    output_file << "\n";

    // Write the zone header
    for (UnsignedInt l = 0; l != clustered_blocks.size(); ++l)
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

        auto cell_pkg_index = mcv_cell_pkg_index_->Data();
        mesh_for_column_major(
            global_lower_bound, global_upper_bound,
            [&](const Array3i &global_index)
            {
                Vec3d data_position = global_mesh.GridPositionFromIndex(global_index);
                output_file << data_position[0] << " " << data_position[1] << " " << data_position[2] << " ";

                constexpr int type_index_int = DataTypeIndex<int>::value;
                for (MeshVariable<int> *variable : std::get<type_index_int>(mesh_variables_to_write_))
                {
                    int value = index_handler.DataValueFromGlobalIndex(variable->Data(), global_index, cell_pkg_index);
                    output_file << value << " ";
                };

                constexpr int type_index_Vecd = DataTypeIndex<Vec3d>::value;
                for (MeshVariable<Vec3d> *variable : std::get<type_index_Vecd>(mesh_variables_to_write_))
                {
                    Vec3d value = index_handler.DataValueFromGlobalIndex(variable->Data(), global_index, cell_pkg_index);
                    output_file << value[0] << " " << value[1] << " " << value[2] << " ";
                };

                constexpr int type_index_Real = DataTypeIndex<Real>::value;
                for (MeshVariable<Real> *variable : std::get<type_index_Real>(mesh_variables_to_write_))
                {
                    Real value = index_handler.DataValueFromGlobalIndex(variable->Data(), global_index, cell_pkg_index);
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
