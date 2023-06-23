#include "triangle_mesh_shape.h"

namespace SPH
{
//=================================================================================================//
SimTK::ContactGeometry::TriangleMesh *TriangleMeshShape::generateTriangleMesh(const SimTK::PolygonalMesh &poly_mesh)
{
    SimTK::ContactGeometry::TriangleMesh *triangle_mesh;
    triangle_mesh = triangle_mesh_ptr_keeper_.createPtr<SimTK::ContactGeometry::TriangleMesh>(poly_mesh);
    if (!SimTK::ContactGeometry::TriangleMesh::isInstance(*triangle_mesh))
    {
        std::cout << "\n Error the triangle mesh is not valid" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        throw;
    }
    std::cout << "num of faces:" << triangle_mesh->getNumFaces() << std::endl;

    return triangle_mesh;
}
//=================================================================================================//
bool TriangleMeshShape::checkContain(const Vec3d &probe_point, bool BOUNDARY_INCLUDED)
{
    SimTKVec2 uv_coordinate;
    bool inside = false; // note that direct prediction is not reliable sometime.
    int face_id;

    SimTKVec3 closest_pnt = triangle_mesh_->findNearestPoint(EigenToSimTK(probe_point), inside, face_id, uv_coordinate);
    Vec3d from_face_to_pnt = probe_point - SimTKToEigen(closest_pnt);
    Real distance_to_pnt = from_face_to_pnt.norm();
    Vec3d direction_to_pnt = from_face_to_pnt / (distance_to_pnt + TinyReal);
    Vec3d face_normal = SimTKToEigen(triangle_mesh_->getFaceNormal(face_id));
    Real cosine_angle = face_normal.dot(direction_to_pnt);

    int ite = 0;
    while (fabs(cosine_angle) < Eps)
    {
        Vec3d jittered = probe_point; // jittering
        for (int l = 0; l != probe_point.size(); ++l)
            jittered[l] = probe_point[l] + (((Real)rand() / (RAND_MAX)) - 0.5) * (SqrtEps + distance_to_pnt * 0.1);
        Vec3d from_face_to_jittered = jittered - SimTKToEigen(closest_pnt);
        Vec3d direction_to_jittered = from_face_to_jittered / (from_face_to_jittered.norm() + TinyReal);
        cosine_angle = face_normal.dot(direction_to_jittered);

        ite++;
        if (ite > 100)
        {
            std::cout << "\n Error: TriangleMeshShape::checkContain not able to achieve!  " << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }
    }

    return cosine_angle < 0.0 ? true : false;
}
//=================================================================================================//
Vecd TriangleMeshShape::findClosestPoint(const Vecd &probe_point)
{
    bool inside = false;
    int face_id;
    SimTKVec2 norm;
    SimTKVec3 closest_pnt = triangle_mesh_->findNearestPoint(SimTKVec3(probe_point[0], probe_point[1], probe_point[2]), inside, face_id, norm);
    if (face_id < 0 && face_id > triangle_mesh_->getNumFaces())
    {
        std::cout << "\n Error the nearest point is not valid" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        throw;
    }
    return Vecd(closest_pnt[0], closest_pnt[1], closest_pnt[2]);
}
//=================================================================================================//
BoundingBox TriangleMeshShape::findBounds()
{
    int number_of_vertices = triangle_mesh_->getNumVertices();
    // initial reference values
    Vec3d lower_bound = SimTKToEigen(SimTKVec3(Infinity));
    Vec3d upper_bound = SimTKToEigen(SimTKVec3(-Infinity));
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
TriangleMeshShapeSTL::TriangleMeshShapeSTL(const std::string &filepathname, Vecd translation, Real scale_factor,
                                           const std::string &shape_name)
    : TriangleMeshShape(shape_name)
{
    if (!fs::exists(filepathname))
    {
        std::cout << "\n Error: the input file:" << filepathname << " is not exists" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        throw;
    }
    SimTK::PolygonalMesh polymesh;
    polymesh.loadStlFile(filepathname);
    polymesh.scaleMesh(scale_factor);
    triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(SimTKVec3(translation[0], translation[1], translation[2])));
}
//=================================================================================================//
TriangleMeshShapeSTL::TriangleMeshShapeSTL(const std::string &filepathname, Mat3d rotation,
                                           Vec3d translation, Real scale_factor, const std::string &shape_name)
    : TriangleMeshShape(shape_name)
{
    if (!fs::exists(filepathname))
    {
        std::cout << "\n Error: the input file:" << filepathname << " is not exists" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        throw;
    }
    SimTK::PolygonalMesh polymesh;
    polymesh.loadStlFile(filepathname);

    polymesh.scaleMesh(scale_factor);
    SimTK::Transform_<double> transform(SimTK::Rotation_<double>(SimTKMat33((double)rotation(0, 0), (double)rotation(0, 1), (double)rotation(0, 2),
                                                                            (double)rotation(1, 0), (double)rotation(1, 1), (double)rotation(1, 2),
                                                                            (double)rotation(2, 0), (double)rotation(2, 1), (double)rotation(2, 2))),
                                        SimTKVec3((double)translation[0], (double)translation[1], (double)translation[2]));
    triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(transform));
}
//=================================================================================================//
#ifdef __EMSCRIPTEN__
TriangleMeshShapeSTL::TriangleMeshShapeSTL(const uint8_t *buffer, Vecd translation, Real scale_factor,
                                           const std::string &shape_name)
    : TriangleMeshShape(shape_name)
{
    SimTK::PolygonalMesh polymesh;
    polymesh.loadStlBuffer(buffer);
    polymesh.scaleMesh(scale_factor);
    triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(SimTKVec3(translation[0], translation[1], translation[2])));
}
#endif
//=================================================================================================//
TriangleMeshShapeBrick::TriangleMeshShapeBrick(Vecd halfsize, int resolution, Vecd translation,
                                               const std::string &shape_name)
    : TriangleMeshShape(shape_name)
{
    SimTK::PolygonalMesh polymesh = SimTK::PolygonalMesh::createBrickMesh(SimTKVec3(halfsize[0], halfsize[1], halfsize[2]), resolution);
    triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(SimTKVec3(translation[0], translation[1], translation[2])));
}
//=================================================================================================//
TriangleMeshShapeBrick::TriangleMeshShapeBrick(const TriangleMeshShapeBrick::ShapeParameters &shape_parameters,
                                               const std::string &shape_name)
    : TriangleMeshShapeBrick(shape_parameters.halfsize_, shape_parameters.resolution_,
                             shape_parameters.translation_, shape_name) {}
//=================================================================================================//
TriangleMeshShapeSphere::TriangleMeshShapeSphere(Real radius, int resolution, Vec3d translation,
                                                 const std::string &shape_name)
    : TriangleMeshShape(shape_name)
{
    SimTK::PolygonalMesh polymesh = SimTK::PolygonalMesh::createSphereMesh(radius, resolution);
    triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(SimTKVec3(translation[0], translation[1], translation[2])));
}
//=================================================================================================//
TriangleMeshShapeCylinder::TriangleMeshShapeCylinder(SimTK::UnitVec3 axis, Real radius, Real halflength, int resolution,
                                                     Vecd translation, const std::string &shape_name)
    : TriangleMeshShape(shape_name)
{
    SimTK::PolygonalMesh polymesh =
        SimTK::PolygonalMesh::createCylinderMesh(axis, radius, halflength, resolution);
    triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(SimTKVec3(translation[0], translation[1], translation[2])));
}
//=================================================================================================//
} // namespace SPH
