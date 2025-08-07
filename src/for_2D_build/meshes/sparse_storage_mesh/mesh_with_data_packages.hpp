#ifndef MESH_WITH_DATA_PACKAGES_2D_HPP
#define MESH_WITH_DATA_PACKAGES_2D_HPP

#include "mesh_with_data_packages.h"

namespace SPH
{
//=================================================================================================//
template <size_t PKG_SIZE>
void MeshWithGridDataPackages<PKG_SIZE>::writeMeshFieldToPltByMesh(std::ofstream &output_file)
{
    output_file << "\n"
                << "title='View'" << "\n";
    output_file << " VARIABLES = " << "x, " << "y";

    constexpr int type_index_int = DataTypeIndex<int>::value;
    for (MeshVariable<int> *variable : std::get<type_index_int>(mesh_variable_to_write_))
    {
        output_file << ",\"" << variable->Name() << "\"";
    };

    constexpr int type_index_Vecd = DataTypeIndex<Vecd>::value;
    for (MeshVariable<Vecd> *variable : std::get<type_index_Vecd>(mesh_variable_to_write_))
    {
        std::string variable_name = variable->Name();
        output_file << ",\"" << variable_name << "_x\""
                    << ",\"" << variable_name << "_y\"";
    };

    constexpr int type_index_Real = DataTypeIndex<Real>::value;
    for (MeshVariable<Real> *variable : std::get<type_index_Real>(mesh_variable_to_write_))
    {
        output_file << ",\"" << variable->Name() << "\"";
    };

    output_file << " \n";

    Arrayi number_of_operation = global_mesh_.AllGridPoints();
    output_file << "zone i=" << number_of_operation[0] << "  j=" << number_of_operation[1] << "  k=" << 1
                << "  DATAPACKING=POINT \n";

    mesh_for_column_major(
        Arrayi::Zero(), number_of_operation,
        [&](const Array2i &global_index)
        {
            Vecd data_position = global_mesh_.GridPositionFromIndex(global_index);
            output_file << data_position[0] << " " << data_position[1] << " ";

            constexpr int type_index_int = DataTypeIndex<int>::value;
            for (MeshVariable<int> *variable : std::get<type_index_int>(mesh_variable_to_write_))
            {
                int value = DataValueFromGlobalIndex(variable->Data(), global_index, this, cell_package_index_.Data());
                output_file << value << " ";
            };

            constexpr int type_index_Vecd = DataTypeIndex<Vecd>::value;
            for (MeshVariable<Vecd> *variable : std::get<type_index_Vecd>(mesh_variable_to_write_))
            {
                Vecd value = DataValueFromGlobalIndex(variable->Data(), global_index, this, cell_package_index_.Data());
                output_file << value[0] << " " << value[1] << " ";
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
//=================================================================================================//
} // namespace SPH
#endif // MESH_WITH_DATA_PACKAGES_2D_HPP
