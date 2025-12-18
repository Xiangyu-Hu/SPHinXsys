#ifndef SPARSE_MESH_FIELD_2D_HPP
#define SPARSE_MESH_FIELD_2D_HPP

#include "spares_mesh_field.hxx"

namespace SPH
{
//=================================================================================================//
//=================================================================================================//
template <int PKG_SIZE>
void SparseMeshField<PKG_SIZE>::writePackageVariablesToPltByMesh(
    UnsignedInt resolution_level, std::ofstream &output_file)
{
    IndexHandler &index_handler = this->getMeshLevel(resolution_level);
    Mesh global_mesh = index_handler.getGlobalMesh();

    output_file << "\n"
                << "title='View'" << "\n";
    output_file << " VARIABLES = " << "x, " << "y";

    constexpr int type_index_unsigned = DataTypeIndex<UnsignedInt>::value;
    for (PackageVariable<UnsignedInt> *variable : std::get<type_index_unsigned>(pkg_variables_to_write_))
    {
        output_file << ",\"" << variable->Name() << "\"";
    };

    constexpr int type_index_int = DataTypeIndex<int>::value;
    for (PackageVariable<int> *variable : std::get<type_index_int>(pkg_variables_to_write_))
    {
        output_file << ",\"" << variable->Name() << "\"";
    };

    constexpr int type_index_Vecd = DataTypeIndex<Vecd>::value;
    for (PackageVariable<Vecd> *variable : std::get<type_index_Vecd>(pkg_variables_to_write_))
    {
        std::string variable_name = variable->Name();
        output_file << ",\"" << variable_name << "_x\""
                    << ",\"" << variable_name << "_y\"";
    };

    constexpr int type_index_Real = DataTypeIndex<Real>::value;
    for (PackageVariable<Real> *variable : std::get<type_index_Real>(pkg_variables_to_write_))
    {
        output_file << ",\"" << variable->Name() << "\"";
    };

    output_file << " \n";

    Arrayi number_of_operation = global_mesh.AllGridPoints();
    output_file << "zone i=" << number_of_operation[0] << "  j=" << number_of_operation[1] << "  k=" << 1
                << "  DATAPACKING=POINT \n";

    UnsignedInt *cell_pkg_index = mcv_cell_pkg_index_->Data();
    mesh_for_column_major(
        Arrayi::Zero(), number_of_operation,
        [&](const Array2i &global_index)
        {
            Vecd data_position = global_mesh.GridPositionFromIndex(global_index);
            output_file << data_position[0] << " " << data_position[1] << " ";

            for (PackageVariable<UnsignedInt> *variable : std::get<type_index_unsigned>(pkg_variables_to_write_))
            {
                UnsignedInt value = index_handler.ValueByGlobalMesh(variable->Data(), global_index, cell_pkg_index);
                output_file << value << " ";
            };

            for (PackageVariable<int> *variable : std::get<type_index_int>(pkg_variables_to_write_))
            {
                int value = index_handler.ValueByGlobalMesh(variable->Data(), global_index, cell_pkg_index);
                output_file << value << " ";
            };

            for (PackageVariable<Vecd> *variable : std::get<type_index_Vecd>(pkg_variables_to_write_))
            {
                Vecd value = index_handler.ValueByGlobalMesh(variable->Data(), global_index, cell_pkg_index);
                output_file << value[0] << " " << value[1] << " ";
            };

            for (PackageVariable<Real> *variable : std::get<type_index_Real>(pkg_variables_to_write_))
            {
                Real value = index_handler.ValueByGlobalMesh(variable->Data(), global_index, cell_pkg_index);
                output_file << value << " ";
            };
            output_file << " \n";
        });
    output_file << " \n";
}
//=================================================================================================//
} // namespace SPH
#endif // SPARSE_MESH_FIELD_2D_HPP
