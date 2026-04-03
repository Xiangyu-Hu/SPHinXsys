#include "base_body_part.h"

#include "triangle_mesh_shape.h"
#include "complex_geometry.h"

namespace SPH
{
//=================================================================================================//
void AlignedBoxPart::writeShapeProxy()
{
    TriangleMeshShapeBrick shape_proxy(
        aligned_box_.HalfSize(), 1, Vecd::Zero(), svAlignedBox()->Name() + "Proxy");
    shape_proxy.writeMeshToFile(aligned_box_.getTransform());
}
//=================================================================================================//
} // namespace SPH
