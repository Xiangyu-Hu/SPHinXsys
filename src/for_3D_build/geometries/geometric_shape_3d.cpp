#include "geometric_shape.h"

#include "triangle_mesh_shape.h"
#include "io_environment.h"

namespace SPH
{
//=================================================================================================//
void GeometricShapeBox::writeProxy()
{
    TriangleMeshShapeBrick shape_proxy(HalfSize(), 1, Vecd::Zero(), getName() + "Proxy");
    shape_proxy.writeMeshToFile(getTransform());
}
//=================================================================================================//
} // namespace SPH