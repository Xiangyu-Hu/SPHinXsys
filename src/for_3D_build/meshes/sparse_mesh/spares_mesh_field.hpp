#ifndef SPARSE_MESH_FIELD_3D_HPP
#define SPARSE_MESH_FIELD_3D_HPP

#include "spares_mesh_field.hxx"

#include "block_cluster_3d.h"

namespace SPH
{
//=================================================================================================//
template <int PKG_SIZE>
void SparseMeshField<PKG_SIZE>::writePackageVariablesToPltByMesh(
    UnsignedInt resolution_level, std::ofstream &output_file)
{
    IndexHandler &index_handler = this->getMesh(resolution_level);

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

    Mesh global_mesh = index_handler.GlobalMesh();

    output_file << "\n"
                << "title='View'" << "\n";
    output_file << " VARIABLES = " << "x, " << "y, " << "z";

    constexpr int type_index_int = DataTypeIndex<int>::value;
    for (PackageVariable<int> *variable : std::get<type_index_int>(pkg_variables_to_write_))
    {
        output_file << ",\"" << variable->Name() << "\"";
    };

    constexpr int type_index_Vecd = DataTypeIndex<Vec3d>::value;
    for (PackageVariable<Vec3d> *variable : std::get<type_index_Vecd>(pkg_variables_to_write_))
    {
        std::string variable_name = variable->Name();
        output_file << ",\"" << variable_name << "_x\""
                    << ",\"" << variable_name << "_y\""
                    << ",\"" << variable_name << "_z\"";
    };

    constexpr int type_index_Real = DataTypeIndex<Real>::value;
    for (PackageVariable<Real> *variable : std::get<type_index_Real>(pkg_variables_to_write_))
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
        Array3i global_lower_bound = PackageDataSize() * lower_bound_cell_index - Array3i::Ones();
        Array3i global_upper_bound = PackageDataSize() * (lower_bound_cell_index + block_size) + Array3i::Ones();
        output_file << "zone i=" << block_size[0] * PackageDataSize() + 2 << "  j="
                    << block_size[1] * PackageDataSize() + 2 << "  k="
                    << block_size[2] * PackageDataSize() + 2
                    << "  DATAPACKING=POINT \n";

        auto cell_pkg_index = mcv_cell_pkg_index_->Data();
        mesh_for_column_major(
            global_lower_bound, global_upper_bound,
            [&](const Array3i &global_index)
            {
                Vec3d data_position = global_mesh.GridPositionFromIndex(global_index);
                output_file << data_position[0] << " " << data_position[1] << " " << data_position[2] << " ";

                constexpr int type_index_int = DataTypeIndex<int>::value;
                for (PackageVariable<int> *variable : std::get<type_index_int>(pkg_variables_to_write_))
                {
                    int value = index_handler.ValueByGlobalMesh(variable->Data(), global_index, cell_pkg_index);
                    output_file << value << " ";
                };

                constexpr int type_index_Vecd = DataTypeIndex<Vec3d>::value;
                for (PackageVariable<Vec3d> *variable : std::get<type_index_Vecd>(pkg_variables_to_write_))
                {
                    Vec3d value = index_handler.ValueByGlobalMesh(variable->Data(), global_index, cell_pkg_index);
                    output_file << value[0] << " " << value[1] << " " << value[2] << " ";
                };

                constexpr int type_index_Real = DataTypeIndex<Real>::value;
                for (PackageVariable<Real> *variable : std::get<type_index_Real>(pkg_variables_to_write_))
                {
                    Real value = index_handler.ValueByGlobalMesh(variable->Data(), global_index, cell_pkg_index);
                    output_file << value << " ";
                };
                output_file << " \n";
            });
        output_file << " \n";
    }
}
//=============================================================================================//
template <int PKG_SIZE>
template <typename DataType>
DataType SparseMeshField<PKG_SIZE>::ProbeMesh<DataType>::probePackageData(
    const IndexHandler &index_handler, UnsignedInt package_index,
    const Array3i &cell_index, const Vec3d &position)
{
    Array3i data_index = index_handler.DataIndexFromPosition(cell_index, position);
    Vec3d data_position = index_handler.DataPositionFromIndex(cell_index, data_index);
    Vec3d alpha = (position - data_position) / index_handler.DataSpacing();
    Vec3d beta = Vec3d::Ones() - alpha;

    auto &neighborhood = cell_neighborhood_[package_index];
    PackageDataPair neighbour_index_1 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(0, 0, 0), neighborhood);
    PackageDataPair neighbour_index_2 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(1, 0, 0), neighborhood);
    PackageDataPair neighbour_index_3 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(0, 1, 0), neighborhood);
    PackageDataPair neighbour_index_4 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(1, 1, 0), neighborhood);

    DataType bilinear_1 =
        pkg_data_[neighbour_index_1.first](neighbour_index_1.second) * beta[0] * beta[1] +
        pkg_data_[neighbour_index_2.first](neighbour_index_2.second) * alpha[0] * beta[1] +
        pkg_data_[neighbour_index_3.first](neighbour_index_3.second) * beta[0] * alpha[1] +
        pkg_data_[neighbour_index_4.first](neighbour_index_4.second) * alpha[0] * alpha[1];

    neighbour_index_1 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(0, 0, 1), neighborhood);
    neighbour_index_2 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(1, 0, 1), neighborhood);
    neighbour_index_3 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(0, 1, 1), neighborhood);
    neighbour_index_4 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(1, 1, 1), neighborhood);

    DataType bilinear_2 =
        pkg_data_[neighbour_index_1.first](neighbour_index_1.second) * beta[0] * beta[1] +
        pkg_data_[neighbour_index_2.first](neighbour_index_2.second) * alpha[0] * beta[1] +
        pkg_data_[neighbour_index_3.first](neighbour_index_3.second) * beta[0] * alpha[1] +
        pkg_data_[neighbour_index_4.first](neighbour_index_4.second) * alpha[0] * alpha[1];

    return bilinear_1 * beta[2] + bilinear_2 * alpha[2];
}
//=================================================================================================//
} // namespace SPH
#endif // SPARSE_MESH_FIELD_3D_HPP
