#include "triangle_mesh_shape.h"

#include "sph_system.h"
namespace SPH
{
//=================================================================================================//
TriangleMeshShape::TriangleMeshShape(const std::string &shape_name) : Shape(shape_name) {}
//=================================================================================================//
void TriangleMeshShape::writeMeshToFile(SPHSystem &sph_system, Transform transform)
{
    std::string filefullpath = sph_system.getIOEnvironment().OutputFolder() + "/" + name_ + ".vtp";

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
    out_file << "<Piece NumberOfPoints=\"" << vertices_.size()
             << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\""
             << faces_.size() << "\">\n";
    out_file << "<Points>\n";
    out_file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (const auto &vertex : vertices_)
    {
        Vecd transformed_vertex = transform.shiftFrameStationToBase(Vecd(vertex[0], vertex[1], vertex[2]));
        out_file << transformed_vertex[0] << " " << transformed_vertex[1] << " " << transformed_vertex[2] << "\n";
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

    // Write file footer
    out_file << "</Piece>\n";
    out_file << "</PolyData>\n";
    out_file << "</VTKFile>\n";

    out_file.close();
}
//=================================================================================================//
void TriangleMeshShape::initializeFromPolygonalMesh(const SimTK::PolygonalMesh &poly_mesh)
{
    SimTK::ContactGeometry::TriangleMesh triangle_mesh(poly_mesh);
    if (!SimTK::ContactGeometry::TriangleMesh::isInstance(triangle_mesh))
    {
        std::cout << "\n Error the triangle mesh is not valid" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        throw;
    }
    std::cout << "TriangleMesh:" << name_ << std::endl;
    std::cout << "num of vertices:" << triangle_mesh.getNumVertices() << std::endl;
    std::cout << "num of faces:" << triangle_mesh.getNumFaces() << std::endl;

    vertices_.reserve(triangle_mesh.getNumVertices());
    for (int i = 0; i < triangle_mesh.getNumVertices(); i++)
    {
        const auto &p = triangle_mesh.getVertexPosition(i);
        vertices_.push_back({Real(p[0]), Real(p[1]), Real(p[2])});
    }

    faces_.reserve(triangle_mesh.getNumFaces());
    for (int i = 0; i < triangle_mesh.getNumFaces(); i++)
    {
        auto f1 = triangle_mesh.getFaceVertex(i, 0);
        auto f2 = triangle_mesh.getFaceVertex(i, 1);
        auto f3 = triangle_mesh.getFaceVertex(i, 2);
        faces_.push_back({f1, f2, f3});
    }
    triangle_mesh_distance_.construct(vertices_, faces_);
}
//=================================================================================================//
void TriangleMeshShape::initializeFromSTLMesh(
    const std::string &file_path_name, Vec3d translation, Real scale_factor)
{
    std::vector<float> coords, normals;
    std::vector<size_t> tris, solids;
    stl_reader::ReadStlFile(file_path_name.c_str(), coords, normals, tris, solids);
    const size_t numVertices = coords.size() / 3;
    const size_t numTriangles = tris.size() / 3;

    vertices_.reserve(numVertices);
    for (size_t i = 0; i < numVertices; i++)
    {
        Real p1 = Real(coords[i * 3]) * scale_factor + Real(translation[0]);
        Real p2 = Real(coords[i * 3 + 1]) * scale_factor + Real(translation[1]);
        Real p3 = Real(coords[i * 3 + 2]) * scale_factor + Real(translation[2]);
        vertices_.push_back({p1, p2, p3});
    }

    faces_.reserve(numTriangles);
    for (size_t i = 0; i < numTriangles; i++)
    {
        size_t f1 = tris[i * 3];
        size_t f2 = tris[i * 3 + 1];
        size_t f3 = tris[i * 3 + 2];
        faces_.push_back({int(f1), int(f2), int(f3)});
    }
    triangle_mesh_distance_.construct(vertices_, faces_);
}
//=================================================================================================//
bool TriangleMeshShape::checkContain(const Vec3d &probe_point, bool BOUNDARY_INCLUDED)
{
    Real distance = triangle_mesh_distance_.signed_distance(probe_point).distance;

    return distance < 0.0 ? true : false;
}
//=================================================================================================//
Vecd TriangleMeshShape::findClosestPoint(const Vecd &probe_point)
{
    auto closest_pnt = triangle_mesh_distance_.signed_distance(probe_point).nearest_point;
    return Vecd(closest_pnt[0], closest_pnt[1], closest_pnt[2]);
}
//=================================================================================================//
BoundingBoxd TriangleMeshShape::findBounds()
{
    // initial reference values
    Vec3d lower_bound = SimTKToEigen(SimTKVec3(MaxReal));
    Vec3d upper_bound = SimTKToEigen(SimTKVec3(MinReal));
    for (size_t i = 0; i != vertices_.size(); ++i)
    {
        Vec3d vertex_position = Vec3d(vertices_[i][0], vertices_[i][1], vertices_[i][2]);
        for (int j = 0; j != 3; ++j)
        {
            lower_bound[j] = SMIN(lower_bound[j], vertex_position[j]);
            upper_bound[j] = SMAX(upper_bound[j], vertex_position[j]);
        }
    }
    return BoundingBoxd(lower_bound, upper_bound);
}
//=================================================================================================//
TriangleMeshShapeBrick::TriangleMeshShapeBrick(Vecd halfsize, int resolution, Vecd translation,
                                               const std::string &shape_name)
    : TriangleMeshShape(shape_name)
{
    SimTK::PolygonalMesh poly_mesh = SimTK::PolygonalMesh::createBrickMesh(
        SimTKVec3(halfsize[0], halfsize[1], halfsize[2]), resolution);
    initializeFromPolygonalMesh(
        poly_mesh.transformMesh(SimTKVec3(translation[0], translation[1], translation[2])));
}
//=================================================================================================//
TriangleMeshShapeBrick::TriangleMeshShapeBrick(
    const TriangleMeshShapeBrick::ShapeParameters &shape_parameters, const std::string &shape_name)
    : TriangleMeshShapeBrick(shape_parameters.halfsize_, shape_parameters.resolution_,
                             shape_parameters.translation_, shape_name) {}
//=================================================================================================//
TriangleMeshShapeSphere::TriangleMeshShapeSphere(Real radius, int resolution, Vec3d translation,
                                                 const std::string &shape_name)
    : TriangleMeshShape(shape_name)
{
    SimTK::PolygonalMesh poly_mesh = SimTK::PolygonalMesh::createSphereMesh(radius, resolution);
    initializeFromPolygonalMesh(
        poly_mesh.transformMesh(SimTKVec3(translation[0], translation[1], translation[2])));
}
//=================================================================================================//
TriangleMeshShapeCylinder::TriangleMeshShapeCylinder(
    Vecd axis, Real radius, Real halflength, int resolution,
    Vecd translation, const std::string &shape_name)
    : TriangleMeshShape(shape_name)
{
    SimTK::PolygonalMesh poly_mesh =
        SimTK::PolygonalMesh::createCylinderMesh(SimTK::UnitVec3(axis[0], axis[1], axis[2]), radius, halflength, resolution);
    initializeFromPolygonalMesh(
        poly_mesh.transformMesh(SimTKVec3(translation[0], translation[1], translation[2])));
}
//=================================================================================================//
TriangleMeshShapeSTL::TriangleMeshShapeSTL(const std::string &filepathname, Vec3d translation,
                                           Real scale_factor, const std::string &shape_name)
    : TriangleMeshShape(shape_name)
{
    if (!fs::exists(filepathname))
    {
        std::cout << "\n Error: the input file:" << filepathname << " is not exists" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        throw;
    }
    initializeFromSTLMesh(filepathname, translation, scale_factor);
}
//=================================================================================================//
} // namespace SPH
