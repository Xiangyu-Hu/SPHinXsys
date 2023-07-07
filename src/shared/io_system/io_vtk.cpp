/**
 * @file 	io_vtk.cpp
 * @author	Luhui Han, Chi Zhang and Xiangyu Hu
 */

#include "io_vtk.h"

namespace SPH
{
//=============================================================================================//
void BodyStatesRecordingToVtp::writeWithFileName(const std::string &sequence)
{
    for (SPHBody *body : bodies_)
    {
        if (body->checkNewlyUpdated())
        {
            std::string filefullpath = io_environment_.output_folder_ + "/" + body->getName() + "_" + sequence + ".vtp";
            if (fs::exists(filefullpath))
            {
                fs::remove(filefullpath);
            }
            std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);
            // begin of the XML file
            out_file << "<?xml version=\"1.0\"?>\n";
            out_file << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
            out_file << " <PolyData>\n";

            BaseParticles &base_particles = body->getBaseParticles();
            size_t total_real_particles = base_particles.total_real_particles_;
            out_file << "  <Piece Name =\"" << body->getName() << "\" NumberOfPoints=\"" << total_real_particles
                     << "\" NumberOfVerts=\"" << total_real_particles << "\">\n";

            body->writeParticlesToVtpFile(out_file);

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
            const auto &vtuName = body->getName() + "_" + sequence + ".vtu";
            std::stringstream sstream;
            // begin of the XML file
            writeVtu(sstream, body);
            _vtuData[vtuName] = sstream.str();
        }
        body->setNotNewlyUpdated();
    }
}
//=============================================================================================//
void BodyStatesRecordingToVtpString::writeVtu(std::ostream &stream, SPHBody *body) const
{
    stream << "<?xml version=\"1.0\"?>\n";
    stream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    stream << " <UnstructuredGrid>\n";

    BaseParticles &base_particles = body->getBaseParticles();
    size_t total_real_particles = base_particles.total_real_particles_;
    stream << "  <Piece Name =\"" << body->getName() << "\" NumberOfPoints=\"" << total_real_particles << "\" NumberOfCells=\"0\">\n";

    body->writeParticlesToVtuFile(stream);

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
    WriteToVtpIfVelocityOutOfBound(IOEnvironment &io_environment, SPHBodyVector bodies, Real velocity_bound)
    : BodyStatesRecordingToVtp(io_environment, bodies), out_of_bound_(false)
{
    for (size_t i = 0; i < bodies_.size(); ++i)
    {
        check_bodies_.push_back(
            check_bodies_ptr_keeper_.createPtr<ReduceDynamics<VelocityBoundCheck>>(*bodies[i], velocity_bound));
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
//=================================================================================================//
} // namespace SPH
