#include "triangle_mesh_shape.h"

namespace SPH
{
//=================================================================================================//
TriangleMeshShape::TriangleMeshShape(const std::string &shape_name)
    : Shape(shape_name), triangle_mesh_(nullptr) {}
//=================================================================================================//
void TriangleMeshShape::initializeFromPolygonalMesh(const SimTK::PolygonalMesh &poly_mesh)
{
    triangle_mesh_ = triangle_mesh_ptr_keeper_.createPtr<TriangleMesh>(poly_mesh);
    if (!TriangleMesh::isInstance(*triangle_mesh_))
    {
        std::cout << "\n Error the triangle mesh is not valid" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        throw;
    }
    std::cout << "TriangleMesh:" << name_ << std::endl;
    std::cout << "num of vertices:" << triangle_mesh_->getNumVertices() << std::endl;
    std::cout << "num of faces:" << triangle_mesh_->getNumFaces() << std::endl;

    std::vector<std::array<Real, 3>> vertices;
    vertices.reserve(triangle_mesh_->getNumVertices());
    for (int i = 0; i < triangle_mesh_->getNumVertices(); i++)
    {
        const auto &p = triangle_mesh_->getVertexPosition(i);
        vertices.push_back({Real(p[0]), Real(p[1]), Real(p[2])});
    }

    std::vector<std::array<int, 3>> faces;
    faces.reserve(triangle_mesh_->getNumFaces());
    for (int i = 0; i < triangle_mesh_->getNumFaces(); i++)
    {
        auto f1 = triangle_mesh_->getFaceVertex(i, 0);
        auto f2 = triangle_mesh_->getFaceVertex(i, 1);
        auto f3 = triangle_mesh_->getFaceVertex(i, 2);
        faces.push_back({f1, f2, f3});
    }
    triangle_mesh_distance_.construct(vertices, faces);
}
//=================================================================================================//
TriangleMesh *TriangleMeshShape::getTriangleMesh()
{
    if (triangle_mesh_ == nullptr)
    {
        std::cout << "\n Error: TriangleMesh not setup yet! \n";
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    return triangle_mesh_;
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
    int number_of_vertices = triangle_mesh_->getNumVertices();
    // initial reference values
    Vec3d lower_bound = SimTKToEigen(SimTKVec3(MaxReal));
    Vec3d upper_bound = SimTKToEigen(SimTKVec3(MinReal));
    for (int i = 0; i != number_of_vertices; ++i)
    {
        Vec3d vertex_position = SimTKToEigen(triangle_mesh_->getVertexPosition(i));
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
    SimTK::PolygonalMesh poly_mesh;
    poly_mesh.loadStlFile(filepathname);
    poly_mesh.scaleMesh(scale_factor);
    poly_mesh.transformMesh(SimTKVec3(translation[0], translation[1], translation[2]));
    initializeFromPolygonalMesh(poly_mesh);
}
//=================================================================================================//
} // namespace SPH
