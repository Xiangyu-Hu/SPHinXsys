#include "io_vtk.hpp"
#include "io_vtk_mesh.h"
#include "mesh_helper.h"

namespace SPH
{
//=================================================================================================//
void BodyStatesRecordingToMeshVtu::writeWithFileName(const std::string &sequence)
{
    for (SPHBody *body : bodies_)
    {
        if (body->checkNewlyUpdated() && state_recording_)
        {
            std::string filefullpath = io_environment_.output_folder_ + "/SPHBody_" + body->getName() + "_" + sequence + ".vtu";
            if (fs::exists(filefullpath))
            {
                fs::remove(filefullpath);
            }
            std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);

            MeshFileHelpers::vtuFileHeader(out_file);
            Real range_max = 0.0;
            MeshFileHelpers::vtuFileNodeCoordinates(out_file, node_coordinates_, elements_nodes_connection_, bounds_, range_max);

            MeshFileHelpers::vtuFileInformationKey(out_file, range_max);

            MeshFileHelpers::vtuFileCellConnectivity(out_file, elements_nodes_connection_, node_coordinates_);

            MeshFileHelpers::vtuFileOffsets(out_file, elements_nodes_connection_);

            MeshFileHelpers::vtuFileTypeOfCell(out_file, elements_nodes_connection_);

            // write Particle data to vtu file
            out_file << "<CellData>\n";

            BaseParticles &particles = body->getBaseParticles();
            writeParticlesToVtk(out_file, particles);

            out_file << "</CellData>\n";
            // Write VTU file footer
            out_file << "</Piece>\n";
            out_file << "</UnstructuredGrid>\n";
            out_file << "</VTKFile>\n";
            out_file.close();
        }
        body->setNotNewlyUpdated();
    }
}
//=================================================================================================//
void BodyStatesRecordingToMeshVtp::writeWithFileName(const std::string &sequence)
{
    std::cout << "For 3D build:"
              << "The method BodyStatesRecordingToMeshVtp::writeWithFileName not implemented yet."
              << std::endl;
    exit(1);
}
//=================================================================================================//
} // namespace SPH
