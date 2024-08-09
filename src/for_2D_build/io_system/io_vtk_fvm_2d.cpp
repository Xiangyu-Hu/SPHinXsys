#include "io_vtk.hpp"
#include "io_vtk_fvm.h"
namespace SPH
{
//=================================================================================================//
void BodyStatesRecordingInMeshToVtp::writeWithFileName(const std::string &sequence)
{
    for (SPHBody *body : bodies_)
    {
        if (body->checkNewlyUpdated() && state_recording_)
        {
            // TODO: we can short the file name by without using SPHBody
            std::string filefullpath = io_environment_.output_folder_ + "/" + body->getName() + "_" + sequence + ".vtp";
            if (fs::exists(filefullpath))
            {
                fs::remove(filefullpath);
            }
            std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);
            // begin of the XML file
            out_file << "<?xml version=\"1.0\"?>\n";
            out_file << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
            out_file << "<PolyData>\n";

            // Write point data
            out_file << "<Piece NumberOfPoints=\"" << node_coordinates_.size()
                     << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\""
                     << elements_nodes_connection_.size() << "\">\n";
            out_file << "<Points>\n";
            out_file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

            size_t total_nodes = node_coordinates_.size();
            for (size_t node = 0; node != total_nodes; ++node)
            {
                Vec3d particle_position = upgradeToVec3d(node_coordinates_[node]);
                out_file << particle_position[0] << " " << particle_position[1] << " " << particle_position[2] << "\n";
            }

            out_file << "</DataArray>\n";
            out_file << "</Points>\n";

            // Write face data
            out_file << "<Polys>\n";
            out_file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

            for (const auto &element : elements_nodes_connection_)
            {
                for (const auto &vertex : element)
                {
                    out_file << vertex << " ";
                }
                out_file << "\n";
            }

            out_file << "</DataArray>\n";
            out_file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";

            size_t offset = 0;
            for (const auto &face : elements_nodes_connection_)
            {
                offset += face.size();
                out_file << offset << " ";
            }

            out_file << "\n</DataArray>\n";
            out_file << "</Polys>\n";

            // Write face attribute data
            out_file << "<CellData>\n";

            BaseParticles &particles = body->getBaseParticles();
            writeParticlesToVtk(out_file, particles);

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
void BodyStatesRecordingInMeshToVtu::writeWithFileName(const std::string &sequence)
{
    std::cout << "For 2D build:"
              << "The method BodyStatesRecordingInMeshToVtu::writeWithFileName not implemented yet."
              << std::endl;
    exit(1);
}
//=================================================================================================//
} // namespace SPH
