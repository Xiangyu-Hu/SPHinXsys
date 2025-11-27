#include "base_mesh.hpp"

#include "mesh_iterators.hpp"

namespace SPH
{
//=================================================================================================//
void MultiLevelMeshField::writeCellVariableToPltByMesh(const Mesh &mesh, std::ofstream &output_file)
{
    output_file << "\n"
                << "title='View'" << "\n";
    output_file << " VARIABLES = " << "x, " << "y";

    constexpr int type_index_unsigned = DataTypeIndex<UnsignedInt>::value;
    for (DiscreteVariable<UnsignedInt> *variable : std::get<type_index_unsigned>(cell_variables_to_write_))
    {
        output_file << ",\"" << variable->Name() << "\"";
    };

    constexpr int type_index_int = DataTypeIndex<int>::value;
    for (DiscreteVariable<int> *variable : std::get<type_index_int>(cell_variables_to_write_))
    {
        output_file << ",\"" << variable->Name() << "\"";
    };

    constexpr int type_index_Vecd = DataTypeIndex<Vecd>::value;
    for (DiscreteVariable<Vecd> *variable : std::get<type_index_Vecd>(cell_variables_to_write_))
    {
        std::string variable_name = variable->Name();
        output_file << ",\"" << variable_name << "_x\""
                    << ",\"" << variable_name << "_y\"";
    };

    constexpr int type_index_Real = DataTypeIndex<Real>::value;
    for (DiscreteVariable<Real> *variable : std::get<type_index_Real>(cell_variables_to_write_))
    {
        output_file << ",\"" << variable->Name() << "\"";
    };

    output_file << " \n";

    Arrayi number_of_operation = mesh.AllCells();
    output_file << "zone i=" << number_of_operation[0] << "  j=" << number_of_operation[1] << "  k=" << 1
                << "  DATAPACKING=POINT \n";

    mesh_for_column_major(
        Arrayi::Zero(), number_of_operation,
        [&](const Array2i &cell_index)
        {
            UnsignedInt linear_index = mesh.LinearCellIndex(cell_index);
            Vecd data_position = mesh.CellPositionFromIndex(cell_index);
            output_file << data_position[0] << " " << data_position[1] << " ";

            for (DiscreteVariable<UnsignedInt> *variable : std::get<type_index_unsigned>(cell_variables_to_write_))
            {
                UnsignedInt value = variable->Data()[linear_index];
                output_file << value << " ";
            };

            for (DiscreteVariable<int> *variable : std::get<type_index_int>(cell_variables_to_write_))
            {
                int value = variable->Data()[linear_index];
                output_file << value << " ";
            };

            for (DiscreteVariable<Vecd> *variable : std::get<type_index_Vecd>(cell_variables_to_write_))
            {
                Vecd value = variable->Data()[linear_index];
                output_file << value[0] << " " << value[1] << " ";
            };

            for (DiscreteVariable<Real> *variable : std::get<type_index_Real>(cell_variables_to_write_))
            {
                Real value = variable->Data()[linear_index];
                output_file << value << " ";
            };
            output_file << " \n";
        });
    output_file << " \n";
}
//=================================================================================================//
} // namespace SPH
