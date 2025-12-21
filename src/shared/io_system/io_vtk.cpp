#include "io_vtk.hpp"

#include "io_environment.h"

namespace SPH
{
//=============================================================================================//
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
                out_file << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
                out_file << " <PolyData>\n";

                size_t total_real_particles = base_particles.TotalRealParticles();
                out_file << "  <Piece Name =\"" << body->getName() << "\" NumberOfPoints=\"" << total_real_particles
                         << "\" NumberOfVerts=\"" << total_real_particles << "\">\n";

                // write current/final particle positions first
                out_file << "   <Points>\n";
                out_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
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
                out_file << "    <DataArray type=\"Int32\"  Name=\"connectivity\"  Format=\"ascii\">\n";
                out_file << "    ";
                for (size_t i = 0; i != total_real_particles; ++i)
                {
                    out_file << i << " ";
                }
                out_file << std::endl;
                out_file << "    </DataArray>\n";
                out_file << "    <DataArray type=\"Int32\"  Name=\"offsets\"  Format=\"ascii\">\n";
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
    stream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    stream << " <UnstructuredGrid>\n";

    BaseParticles &base_particles = body->getBaseParticles();
    size_t total_real_particles = base_particles.TotalRealParticles();
    stream << "  <Piece Name =\"" << body->getName() << "\" NumberOfPoints=\"" << total_real_particles << "\" NumberOfCells=\"0\">\n";

    writeParticlesToVtk(stream, base_particles);

    stream << "   </PointData>\n";

    // write empty cells
    stream << "   <Cells>\n";
    stream << "    <DataArray type=\"Int32\"  Name=\"connectivity\"  Format=\"ascii\">\n";
    stream << "    </DataArray>\n";
    stream << "    <DataArray type=\"Int32\"  Name=\"offsets\"  Format=\"ascii\">\n";
    stream << "    </DataArray>\n";
    stream << "    <DataArray type=\"types\"  Name=\"offsets\"  Format=\"ascii\">\n";
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
        out_file << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        out_file << " <PolyData>\n";

        out_file << "  <Piece Name =\"" << sph_body_.getName() << "\" NumberOfPoints=\"" << total_generated_particles
                 << "\" NumberOfVerts=\"" << total_generated_particles << "\">\n";

        // write current/final particle positions first
        out_file << "   <Points>\n";
        out_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
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
        out_file << "    <DataArray type=\"Int32\"  Name=\"connectivity\"  Format=\"ascii\">\n";
        out_file << "    ";
        for (size_t i = 0; i != total_generated_particles; ++i)
        {
            out_file << i << " ";
        }
        out_file << std::endl;
        out_file << "    </DataArray>\n";
        out_file << "    <DataArray type=\"Int32\"  Name=\"offsets\"  Format=\"ascii\">\n";
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
