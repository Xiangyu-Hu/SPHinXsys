#include "triangle_mesh_shape.h"

namespace SPH
{
//=================================================================================================//
TriangleMeshShape::TriangleMeshShape(const std::string &shape_name)
    : Shape(shape_name) {}
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
BoundingBox TriangleMeshShape::findBounds()
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
    return BoundingBox(lower_bound, upper_bound);
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
    SimTK::UnitVec3 axis, Real radius, Real halflength, int resolution,
    Vecd translation, const std::string &shape_name)
    : TriangleMeshShape(shape_name)
{
    SimTK::PolygonalMesh poly_mesh =
        SimTK::PolygonalMesh::createCylinderMesh(axis, radius, halflength, resolution);
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
