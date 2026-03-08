#include "base_body_part.h"

#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
void AlignedBoxPart::writeShapeProxy()
{
    TriangleMeshShapeBrick shape_proxy(
        aligned_box_.HalfSize(), 1, Vecd::Zero(), svAlignedBox()->Name() + "Proxy");
    shape_proxy.writeMeshToFile(sph_system_, aligned_box_.getTransform());
}
//=================================================================================================//
} // namespace SPH
