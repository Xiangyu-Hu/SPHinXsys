#include "io_vtk.hpp"

#include "io_environment.h"
#include "sph_system.h"

namespace SPH
{
//=============================================================================================//
#ifdef SPHINXSYS_USE_VTK
void BodyStatesRecordingToVtp::writeWithFileName(const std::string &sequence)
{
    for (SPHBody *body : bodies_)
    {
        if (body->checkNewlyUpdated())
        {
            BaseParticles &base_particles = body->getBaseParticles();

            if (state_recording_)
            {
                std::string filefullpath = io_environment_.OutputFolder() + "/" + body->getName() + "_" + sequence + ".vtp";

                size_t total_real_particles = base_particles.TotalRealParticles();

                // Build points
                vtkNew<vtkPoints> points;
                points->SetDataTypeToFloat();
                for (size_t i = 0; i != total_real_particles; ++i)
                {
                    Vec3d pos = upgradeToVec3d(base_particles.ParticlePositions()[i]);
                    points->InsertNextPoint(pos[0], pos[1], pos[2]);
                }

                // Build vertex cells (one vertex per particle)
                vtkNew<vtkCellArray> verts;
                for (vtkIdType i = 0; i < static_cast<vtkIdType>(total_real_particles); ++i)
                {
                    verts->InsertNextCell(1, &i);
                }

                vtkNew<vtkPolyData> polydata;
                polydata->SetPoints(points);
                polydata->SetVerts(verts);

                // Add physical time as field data
                if (sph_system_.isPhysical())
                {
                    vtkNew<vtkDoubleArray> time_array;
                    time_array->SetName("TimeValue");
                    time_array->SetNumberOfTuples(1);
                    time_array->SetValue(0, sv_physical_time_->getValue());
                    polydata->GetFieldData()->AddArray(time_array);
                }

                // Add particle data arrays to point data
                addParticlesToVtkPolyData(polydata, base_particles);

                vtkNew<vtkXMLPolyDataWriter> writer;
                writer->SetInputData(polydata);
                writer->SetFileName(filefullpath.c_str());
                if (binary_output_)
                    writer->SetDataModeToBinary();
                else
                    writer->SetDataModeToAscii();
                writer->Write();
            }
        }
        body->setNotNewlyUpdated();
    }
}
//=============================================================================================//
void BodyStatesRecordingToVtp::addParticlesToVtkPolyData(vtkPolyData *polydata, BaseParticles &particles)
{
    size_t total_real_particles = particles.TotalRealParticles();
    ParticleVariables &variables_to_write = particles.VariablesToWrite();

    // Write sorted particle IDs
    vtkNew<vtkIntArray> sorted_id_array;
    sorted_id_array->SetName("SortedParticle_ID");
    sorted_id_array->SetNumberOfValues(static_cast<vtkIdType>(total_real_particles));
    for (size_t i = 0; i < total_real_particles; ++i)
        sorted_id_array->SetValue(static_cast<vtkIdType>(i), static_cast<int>(i));
    polydata->GetPointData()->AddArray(sorted_id_array);

    // Write UnsignedInt variables
    constexpr int type_index_UnsignedInt = DataTypeIndex<UnsignedInt>::value;
    for (DiscreteVariable<UnsignedInt> *variable : std::get<type_index_UnsignedInt>(variables_to_write))
    {
        UnsignedInt *data_field = variable->Data();
        vtkNew<vtkUnsignedIntArray> arr;
        arr->SetName(variable->Name().c_str());
        arr->SetNumberOfValues(static_cast<vtkIdType>(total_real_particles));
        for (size_t i = 0; i < total_real_particles; ++i)
            arr->SetValue(static_cast<vtkIdType>(i), data_field[i]);
        polydata->GetPointData()->AddArray(arr);
    }

    // Write int variables
    constexpr int type_index_int = DataTypeIndex<int>::value;
    for (DiscreteVariable<int> *variable : std::get<type_index_int>(variables_to_write))
    {
        int *data_field = variable->Data();
        vtkNew<vtkIntArray> arr;
        arr->SetName(variable->Name().c_str());
        arr->SetNumberOfValues(static_cast<vtkIdType>(total_real_particles));
        for (size_t i = 0; i < total_real_particles; ++i)
            arr->SetValue(static_cast<vtkIdType>(i), data_field[i]);
        polydata->GetPointData()->AddArray(arr);
    }

    // Write Real (scalar) variables
    constexpr int type_index_Real = DataTypeIndex<Real>::value;
    for (DiscreteVariable<Real> *variable : std::get<type_index_Real>(variables_to_write))
    {
        Real *data_field = variable->Data();
        vtkNew<vtkFloatArray> arr;
        arr->SetName(variable->Name().c_str());
        arr->SetNumberOfValues(static_cast<vtkIdType>(total_real_particles));
        for (size_t i = 0; i < total_real_particles; ++i)
            arr->SetValue(static_cast<vtkIdType>(i), static_cast<float>(data_field[i]));
        polydata->GetPointData()->AddArray(arr);
    }

    // Write Vecd (vector) variables
    constexpr int type_index_Vecd = DataTypeIndex<Vecd>::value;
    for (DiscreteVariable<Vecd> *variable : std::get<type_index_Vecd>(variables_to_write))
    {
        Vecd *data_field = variable->Data();
        vtkNew<vtkFloatArray> arr;
        arr->SetName(variable->Name().c_str());
        arr->SetNumberOfComponents(3);
        arr->SetNumberOfTuples(static_cast<vtkIdType>(total_real_particles));
        for (size_t i = 0; i < total_real_particles; ++i)
        {
            Vec3d vec = upgradeToVec3d(data_field[i]);
            float tuple[3] = {static_cast<float>(vec[0]), static_cast<float>(vec[1]), static_cast<float>(vec[2])};
            arr->SetTuple(static_cast<vtkIdType>(i), tuple);
        }
        polydata->GetPointData()->AddArray(arr);
    }

    // Write Matd (matrix) variables
    constexpr int type_index_Matd = DataTypeIndex<Matd>::value;
    for (DiscreteVariable<Matd> *variable : std::get<type_index_Matd>(variables_to_write))
    {
        Matd *data_field = variable->Data();
        vtkNew<vtkFloatArray> arr;
        arr->SetName(variable->Name().c_str());
        arr->SetNumberOfComponents(9);
        arr->SetNumberOfTuples(static_cast<vtkIdType>(total_real_particles));
        for (size_t i = 0; i < total_real_particles; ++i)
        {
            Mat3d mat = upgradeToMat3d(data_field[i]);
            float tuple[9];
            for (int k = 0; k < 3; ++k)
            {
                Vec3d col = mat.col(k);
                tuple[3 * k + 0] = static_cast<float>(col[0]);
                tuple[3 * k + 1] = static_cast<float>(col[1]);
                tuple[3 * k + 2] = static_cast<float>(col[2]);
            }
            arr->SetTuple(static_cast<vtkIdType>(i), tuple);
        }
        polydata->GetPointData()->AddArray(arr);
    }
}
#else
void BodyStatesRecordingToVtp::writeWithFileName(const std::string &sequence)
{
    for (SPHBody *body : bodies_)
    {
        if (body->checkNewlyUpdated())
        {
            BaseParticles &base_particles = body->getBaseParticles();

            if (state_recording_)
            {
                std::string filefullpath = io_environment_.OutputFolder() + "/" + body->getName() + "_" + sequence + ".vtp";
                if (fs::exists(filefullpath))
                {
                    fs::remove(filefullpath);
                }
                std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);
                // begin of the XML file
                out_file << "<?xml version=\"1.0\"?>\n";
                out_file << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
                out_file << " <PolyData>\n";

                // physical time
                if (sph_system_.isPhysical())
                {
                    out_file << "<FieldData>\n";
                    out_file << "<DataArray type=\"Float64\"  Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">\n";
                    out_file << std::fixed << std::setprecision(9) << sv_physical_time_->getValue() << "\n";
                    out_file << " </DataArray>\n";
                    out_file << "</FieldData>\n";
                }

                size_t total_real_particles = base_particles.TotalRealParticles();
                out_file << "  <Piece Name =\"" << body->getName() << "\" NumberOfPoints=\"" << total_real_particles
                         << "\" NumberOfVerts=\"" << total_real_particles << "\">\n";

                // write current/final particle positions first
                out_file << "   <Points>\n";
                out_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" format=\"ascii\">\n";
                out_file << "    ";
                for (size_t i = 0; i != total_real_particles; ++i)
                {
                    Vec3d particle_position = upgradeToVec3d(base_particles.ParticlePositions()[i]);
                    out_file << particle_position[0] << " " << particle_position[1] << " " << particle_position[2] << " ";
                }
                out_file << std::endl;
                out_file << "    </DataArray>\n";
                out_file << "   </Points>\n";

                // write header of particles data
                out_file << "   <PointData  Vectors=\"vector\">\n";
                writeParticlesToVtk(out_file, base_particles);
                out_file << "   </PointData>\n";

                // write empty cells
                out_file << "   <Verts>\n";
                out_file << "    <DataArray type=\"Int32\"  Name=\"connectivity\"  format=\"ascii\">\n";
                out_file << "    ";
                for (size_t i = 0; i != total_real_particles; ++i)
                {
                    out_file << i << " ";
                }
                out_file << std::endl;
                out_file << "    </DataArray>\n";
                out_file << "    <DataArray type=\"Int32\"  Name=\"offsets\"  format=\"ascii\">\n";
                out_file << "    ";
                for (size_t i = 0; i != total_real_particles; ++i)
                {
                    out_file << i + 1 << " ";
                }
                out_file << std::endl;
                out_file << "    </DataArray>\n";
                out_file << "   </Verts>\n";

                out_file << "  </Piece>\n";
                out_file << " </PolyData>\n";
                out_file << "</VTKFile>\n";

                out_file.close();
            }
        }
        body->setNotNewlyUpdated();
    }
}
#endif
//=============================================================================================//
void BodyStatesRecordingToVtpString::writeWithFileName(const std::string &sequence)
{
    for (SPHBody *body : bodies_)
    {
        if (body->checkNewlyUpdated())
        {
            if (state_recording_)
            {
                const auto &vtuName = body->getName() + "_" + sequence + ".vtu";
                std::stringstream sstream;
                // begin of the XML file
                writeVtu(sstream, body);
                _vtuData[vtuName] = sstream.str();
            }
        }
        body->setNotNewlyUpdated();
    }
}
//=============================================================================================//
void BodyStatesRecordingToVtpString::writeVtu(std::ostream &stream, SPHBody *body)
{
    stream << "<?xml version=\"1.0\"?>\n";
    stream << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    stream << " <UnstructuredGrid>\n";

    BaseParticles &base_particles = body->getBaseParticles();
    size_t total_real_particles = base_particles.TotalRealParticles();
    stream << "  <Piece Name =\"" << body->getName() << "\" NumberOfPoints=\"" << total_real_particles << "\" NumberOfCells=\"0\">\n";

    stream << "   <PointData  Vectors=\"vector\">\n";
    writeParticlesToVtk(stream, base_particles);
    stream << "   </PointData>\n";

    // write empty cells
    stream << "   <Cells>\n";
    stream << "    <DataArray type=\"Int32\"  Name=\"connectivity\"  format=\"ascii\">\n";
    stream << "    </DataArray>\n";
    stream << "    <DataArray type=\"Int32\"  Name=\"offsets\"  format=\"ascii\">\n";
    stream << "    </DataArray>\n";
    stream << "    <DataArray type=\"UInt8\"  Name=\"types\"  format=\"ascii\">\n";
    stream << "    </DataArray>\n";
    stream << "   </Cells>\n";

    stream << "  </Piece>\n";

    stream << " </UnstructuredGrid>\n";
    stream << "</VTKFile>\n";
}
//=============================================================================================//
const VtuStringData &BodyStatesRecordingToVtpString::GetVtuData() const
{
    return _vtuData;
}
//=============================================================================================//
WriteToVtpIfVelocityOutOfBound::
    WriteToVtpIfVelocityOutOfBound(SPHSystem &sph_system, Real velocity_bound)
    : BodyStatesRecordingToVtp(sph_system), out_of_bound_(false)
{
    for (size_t i = 0; i < bodies_.size(); ++i)
    {
        check_bodies_.push_back(
            check_bodies_keeper_.createPtr<ReduceDynamics<VelocityBoundCheck>>(*bodies_[i], velocity_bound));
    }
}
//=============================================================================================//
void WriteToVtpIfVelocityOutOfBound::writeWithFileName(const std::string &sequence)
{
    for (auto check_body : check_bodies_)
    {
        out_of_bound_ = out_of_bound_ || check_body->exec();
    }

    if (out_of_bound_)
    {
        BodyStatesRecordingToVtp::writeWithFileName(sequence);
        std::cout << "\n Velocity is out of bound at iteration step " << sequence
                  << "\n The body states have been outputted and the simulation terminates here. \n";
    }
}
//=============================================================================================//
void ParticleGenerationRecordingToVtp::writeWithFileName(const std::string &sequence)
{

    if (state_recording_)
    {
        std::string filefullpath = io_environment_.OutputFolder() + "/" + sph_body_.getName() +
                                   "particle_generation_" + sequence + ".vtp";
        if (fs::exists(filefullpath))
        {
            fs::remove(filefullpath);
        }
        std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);

        size_t total_generated_particles = position_.size();
        // begin of the XML file
        out_file << "<?xml version=\"1.0\"?>\n";
        out_file << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
        out_file << " <PolyData>\n";

        out_file << "  <Piece Name =\"" << sph_body_.getName() << "\" NumberOfPoints=\"" << total_generated_particles
                 << "\" NumberOfVerts=\"" << total_generated_particles << "\">\n";

        // write current/final particle positions first
        out_file << "   <Points>\n";
        out_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" format=\"ascii\">\n";
        out_file << "    ";
        for (size_t i = 0; i != total_generated_particles; ++i)
        {
            Vec3d particle_position = upgradeToVec3d(position_[i]);
            out_file << particle_position[0] << " " << particle_position[1] << " " << particle_position[2] << " ";
        }
        out_file << std::endl;
        out_file << "    </DataArray>\n";
        out_file << "   </Points>\n";

        // write empty cells
        out_file << "   <Verts>\n";
        out_file << "    <DataArray type=\"Int32\"  Name=\"connectivity\"  format=\"ascii\">\n";
        out_file << "    ";
        for (size_t i = 0; i != total_generated_particles; ++i)
        {
            out_file << i << " ";
        }
        out_file << std::endl;
        out_file << "    </DataArray>\n";
        out_file << "    <DataArray type=\"Int32\"  Name=\"offsets\"  format=\"ascii\">\n";
        out_file << "    ";
        for (size_t i = 0; i != total_generated_particles; ++i)
        {
            out_file << i + 1 << " ";
        }
        out_file << std::endl;
        out_file << "    </DataArray>\n";
        out_file << "   </Verts>\n";

        out_file << "  </Piece>\n";
        out_file << " </PolyData>\n";
        out_file << "</VTKFile>\n";

        out_file.close();
    }
}
//=================================================================================================//
} // namespace SPH
