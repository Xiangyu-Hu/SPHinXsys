#ifndef LEVEL_SET_SHAPE_SYCL_HPP
#define LEVEL_SET_SHAPE_SYCL_HPP

#include "level_set_shape.h"
#include "levelset_adaptation.hpp"
namespace SPH
{
LevelSetShape::LevelSetShape(const ParallelDevicePolicy &par_device, SPHBody &sph_body,
                              Shape &shape, Real refinement_ratio)
    : Shape(shape.getName()),
      level_set_(*level_set_keeper_.movePtr(
          sph_body.sph_adaptation_->createLevelSet(par_device, shape, refinement_ratio)))
{
    bounding_box_ = shape.getBounds();
    is_bounds_found_ = true;
}
}

#endif