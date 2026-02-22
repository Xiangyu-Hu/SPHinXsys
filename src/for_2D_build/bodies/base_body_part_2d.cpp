#include "base_body_part.h"

#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
void AlignedBoxPart::writeShapeProxy(SPHSystem &sph_system)
{
    std::string filefullpath = sph_system.getIOEnvironment().OutputFolder() + "/" +
                               svAlignedBox()->Name() + "Proxy.vtp";

    if (fs::exists(filefullpath))
    {
        fs::remove(filefullpath);
    }
    std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);

    Vecd halfsize = aligned_box_.HalfSize();
    Transform transform = aligned_box_.getTransform();

    // 4 corners in local frame (z=0 for 2D)
    Vecd local_corners[4] = {
        Vecd(-halfsize[0], -halfsize[1]),
        Vecd(halfsize[0], -halfsize[1]),
        Vecd(halfsize[0], halfsize[1]),
        Vecd(-halfsize[0], halfsize[1])};

    // begin of the XML file
    out_file << "<?xml version=\"1.0\"?>\n";
    out_file << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    out_file << "<PolyData>\n";

    out_file << "<Piece NumberOfPoints=\"4\" NumberOfVerts=\"0\" NumberOfLines=\"0\" "
             << "NumberOfStrips=\"0\" NumberOfPolys=\"1\">\n";
    out_file << "<Points>\n";
    out_file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (int i = 0; i < 4; ++i)
    {
        Vecd global_point = transform.shiftFrameStationToBase(local_corners[i]);
        out_file << global_point[0] << " " << global_point[1] << " 0\n";
    }

    out_file << "</DataArray>\n";
    out_file << "</Points>\n";

    out_file << "<Polys>\n";
    out_file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    out_file << "0 1 2 3\n";
    out_file << "</DataArray>\n";
    out_file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    out_file << "4\n";
    out_file << "</DataArray>\n";
    out_file << "</Polys>\n";

    out_file << "</Piece>\n";
    out_file << "</PolyData>\n";
    out_file << "</VTKFile>\n";

    out_file.close();
}
//=================================================================================================//
} // namespace SPH
