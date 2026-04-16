#include "geometric_shape.h"

#include "io_environment.h"
#include "triangle_mesh_shape.h"

namespace SPH
{
//=================================================================================================//
void GeometricShapeBox::writeGeometricShapeBoxToVtp()
{
    TriangleMeshShapeBrick shape_proxy(HalfSize(), 1, Vecd::Zero(), "Shape" + getName());
    shape_proxy.writeMeshToFile(getTransform());
}
//=================================================================================================//
} // namespace SPH