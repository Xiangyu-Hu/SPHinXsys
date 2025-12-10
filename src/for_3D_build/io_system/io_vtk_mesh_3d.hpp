#ifndef IO_VTK_MESH_3D_HPP
#define IO_VTK_MESH_3D_HPP

#include "io_vtk_mesh.h"

namespace SPH
{
//=============================================================================================//
template <typename OutStreamType>
void BodyStatesRecordingToTriangleMeshVtp::writeCellsToVtk(OutStreamType &output_stream, BaseParticles &particles)
{
    ParticleVariables &variables_to_write = particles.VariablesToWrite();

    // write scalars
    constexpr int type_index_Real = DataTypeIndex<Real>::value;
    for (DiscreteVariable<Real> *variable : std::get<type_index_Real>(variables_to_write))
    {
        Real *data_field = variable->Data();
        output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Float32\" Format=\"ascii\">\n";
        output_stream << "    ";
        for (size_t i = 0; i != faces_.size(); ++i)
        {
            std::array<int, 3> &face = faces_[i];
            Real sum = 0.0;
            for (size_t j = 0; j != face.size(); ++j)
            {
                sum += data_field[face[j]];
            }
            output_stream << std::fixed << std::setprecision(9) << sum / Real(face.size()) << " ";
        }
        output_stream << std::endl;
        output_stream << "    </DataArray>\n";
    }

    // write vectors
    constexpr int type_index_Vecd = DataTypeIndex<Vecd>::value;
    for (DiscreteVariable<Vecd> *variable : std::get<type_index_Vecd>(variables_to_write))
    {
        Vecd *data_field = variable->Data();
        output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
        output_stream << "    ";
        for (size_t i = 0; i != faces_.size(); ++i)
        {
            std::array<int, 3> &face = faces_[i];
            Vec3d sum = ZeroData<Vec3d>::value;
            for (size_t j = 0; j != face.size(); ++j)
            {
                sum += data_field[face[j]];
            }
            Vec3d vector_value = sum / Real(face.size());
            output_stream << std::fixed << std::setprecision(9) << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
        }
        output_stream << std::endl;
        output_stream << "    </DataArray>\n";
    }

    // write matrices
    constexpr int type_index_Matd = DataTypeIndex<Matd>::value;
    for (DiscreteVariable<Matd> *variable : std::get<type_index_Matd>(variables_to_write))
    {
        Matd *data_field = variable->Data();
        output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type= \"Float32\"  NumberOfComponents=\"9\" Format=\"ascii\">\n";
        output_stream << "    ";
        for (size_t i = 0; i != faces_.size(); ++i)
        {
            std::array<int, 3> &face = faces_[i];
            Mat3d sum = ZeroData<Mat3d>::value;
            for (size_t j = 0; j != face.size(); ++j)
            {
                sum += data_field[face[j]];
            }
            Mat3d matrix_value = sum / Real(face.size());

            for (int k = 0; k != 3; ++k)
            {
                Vec3d col_vector = matrix_value.col(k);
                output_stream << std::fixed << std::setprecision(9) << col_vector[0] << " " << col_vector[1] << " " << col_vector[2] << " ";
            }
        }
        output_stream << std::endl;
        output_stream << "    </DataArray>\n";
    }
}
//=============================================================================================//
} // namespace SPH
#endif // IO_VTK_MESH_3D_HPP