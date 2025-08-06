#include "level_set_shape.h"

namespace SPH
{
//=================================================================================================//
LevelSetShape::
    LevelSetShape(const ParallelDevicePolicy &par_device, Shape &shape,
                  SharedPtr<SPHAdaptation> sph_adaptation, Real refinement_ratio, UsageType usage_type)
    : LevelSetShape(shape.getBounds(), shape, sph_adaptation, refinement_ratio)
{
    finishInitialization(par_device, usage_type);
}
//=================================================================================================//
LevelSetShape::LevelSetShape(const ParallelDevicePolicy &par_device,
                             SPHBody &sph_body, Shape &shape, Real refinement_ratio, UsageType usage_type)
    : LevelSetShape(shape.getBounds(), sph_body, shape, refinement_ratio)
{
    finishInitialization(par_device, usage_type);
}
//=================================================================================================//
} // namespace SPH