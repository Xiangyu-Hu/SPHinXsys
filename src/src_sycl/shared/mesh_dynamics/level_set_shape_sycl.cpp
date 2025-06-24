#include "level_set_shape.h"

namespace SPH
{
//=================================================================================================//
LevelSetShape::
    LevelSetShape(const ParallelDevicePolicy &par_device, Shape &shape, SharedPtr<SPHAdaptation> sph_adaptation, Real refinement_ratio)
    : LevelSetShape(shape, sph_adaptation, refinement_ratio)
{
    finishInitialization(par_device);
}
//=================================================================================================//
LevelSetShape::LevelSetShape(const ParallelDevicePolicy &par_device, SPHBody &sph_body, Shape &shape, Real refinement_ratio)
    : LevelSetShape(par_device, sph_body, shape, refinement_ratio)
{
    finishInitialization(par_device);
}
//=================================================================================================//
} // namespace SPH