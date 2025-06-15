#include "io_vtk_mesh_3d.hpp"
#include "io_vtk.hpp"

#include "io_environment.h"
#include "triangle_mesh_shape.h"

namespace SPH
{
//=================================================================================================//
BodyStatesRecordingToTriangleMeshVtp::BodyStatesRecordingToTriangleMeshVtp(
    SPHBody &body, TriangleMeshShape &triangle_mesh_shape)
    : BodyStatesRecordingToVtp(body), faces_(triangle_mesh_shape.getFaces()) {}
//=================================================================================================//
void BodyStatesRecordingToTriangleMeshVtp::writeWithFileName(const std::string &sequence)
{
    for (SPHBody *body : bodies_)
    {
        if (body->checkNewlyUpdated() && state_recording_)
        {
            std::string filefullpath = io_environment_.output_folder_ + "/" + body->getName() + "_" + sequence + ".vtp";
            if (fs::exists(filefullpath))
            {
                fs::remove(filefullpath);
            }
            std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);

            BaseParticles &particles = body->getBaseParticles();

            // begin of the XML file
            out_file << "<?xml version=\"1.0\"?>\n";
            out_file << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
            out_file << "<PolyData>\n";

            // Write point data
            out_file << "<Piece NumberOfPoints=\"" << particles.TotalRealParticles()
                     << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\""
                     << faces_.size() << "\">\n";
            out_file << "<Points>\n";
            out_file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

            Vec3d *pos = particles.ParticlePositions();
            for (size_t i = 0; i != particles.TotalRealParticles(); ++i)
            {
                out_file << pos[i][0] << " " << pos[i][1] << " " << pos[i][2] << "\n";
            }

            out_file << "</DataArray>\n";
            out_file << "</Points>\n";

            // Write face data
            out_file << "<Polys>\n";
            out_file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

            for (const auto &face : faces_)
            {
                for (const auto &vertex : face)
                {
                    out_file << vertex << " ";
                }
                out_file << "\n";
            }

            out_file << "</DataArray>\n";
            out_file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";

            size_t offset = 0;
            for (const auto &face : faces_)
            {
                offset += face.size();
                out_file << offset << " ";
            }

            out_file << "\n</DataArray>\n";
            out_file << "</Polys>\n";

            // Write face attribute data
            out_file << "<CellData>\n";

            writeCellsToVtk(out_file, particles);

            out_file << "</CellData>\n";

            // Write file footer
            out_file << "</Piece>\n";
            out_file << "</PolyData>\n";
            out_file << "</VTKFile>\n";

            out_file.close();
        }
        body->setNotNewlyUpdated();
    }
}
//=================================================================================================//
} // namespace SPH
