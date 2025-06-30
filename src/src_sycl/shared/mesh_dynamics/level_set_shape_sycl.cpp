#include "level_set_shape.h"

namespace SPH
{
//=================================================================================================//
LevelSetShape::
    LevelSetShape(const ParallelDevicePolicy &par_device, Shape &shape, SharedPtr<SPHAdaptation> sph_adaptation, Real refinement_ratio)
    : Shape(shape.getName()), sph_adaptation_(sph_adaptation),
      level_set_(*level_set_keeper_.movePtr(sph_adaptation->createLevelSet(shape, refinement_ratio)))
{
    bounding_box_ = shape.getBounds();
    is_bounds_found_ = true;
    finishInitialization(par_device);
}
//=================================================================================================//
LevelSetShape::LevelSetShape(const ParallelDevicePolicy &par_device, SPHBody &sph_body, Shape &shape, Real refinement_ratio)
    : Shape(shape.getName()),
      level_set_(*level_set_keeper_.movePtr(
          sph_body.getSPHAdaptation().createLevelSet(shape, refinement_ratio)))
{
    bounding_box_ = shape.getBounds();
    is_bounds_found_ = true;
    finishInitialization(par_device);
}
//=================================================================================================//
} // namespace SPH