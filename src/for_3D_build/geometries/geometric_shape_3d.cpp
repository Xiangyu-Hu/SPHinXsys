#include "geometric_shape.h"

#include "io_environment.h"
#include "triangle_mesh_shape.h"

namespace SPH
{
//=================================================================================================//
void GeometricShapeBox::writeGeometricShapeBoxToVtp()
{
    TriangleMeshShapeBrick shape(HalfSize(), 1, Vecd::Zero(), getName());
    shape.writTriangleMeshShapeToVtp(getTransform());
}
//=================================================================================================//
} // namespace SPH