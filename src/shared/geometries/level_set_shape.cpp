#include "level_set_shape.h"

#include "base_body.h"
#include "io_all.h"
#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
LevelSetShape::
    LevelSetShape(Shape &shape, SharedPtr<SPHAdaptation> sph_adaptation, Real refinement_ratio)
    : Shape(shape.getName()), sph_adaptation_(sph_adaptation),
      level_set_(*level_set_keeper_.movePtr(sph_adaptation->createLevelSet(shape, refinement_ratio)))
{
    bounding_box_ = shape.getBounds();
    is_bounds_found_ = true;
}
//=================================================================================================//
LevelSetShape::LevelSetShape(SPHBody &sph_body, Shape &shape, Real refinement_ratio)
    : Shape(shape.getName()),
      level_set_(*level_set_keeper_.movePtr(
          sph_body.sph_adaptation_->createLevelSet(shape, refinement_ratio)))
{
    bounding_box_ = shape.getBounds();
    is_bounds_found_ = true;
}
//=================================================================================================//
void LevelSetShape::writeLevelSet(IOEnvironment &io_environment)
{
    MeshRecordingToPlt write_level_set_to_plt(io_environment, level_set_);
    write_level_set_to_plt.writeToFile(0);
}
//=================================================================================================//
LevelSetShape *LevelSetShape::cleanLevelSet(Real small_shift_factor)
{
    level_set_.cleanInterface(small_shift_factor);
    return this;
}
//=================================================================================================//
LevelSetShape *LevelSetShape::correctLevelSetSign(Real small_shift_factor)
{
    level_set_.correctTopology(small_shift_factor);
    return this;
}
//=================================================================================================//
bool LevelSetShape::checkContain(const Vecd &probe_point, bool BOUNDARY_INCLUDED)
{
    return level_set_.probeSignedDistance(probe_point) < 0.0 ? true : false;
}
//=================================================================================================//
Vecd LevelSetShape::findClosestPoint(const Vecd &probe_point)
{
    Real phi = level_set_.probeSignedDistance(probe_point);
    Vecd normal = level_set_.probeNormalDirection(probe_point);
    return probe_point - phi * normal;
}
//=================================================================================================//
BoundingBox LevelSetShape::findBounds()
{
    if (!is_bounds_found_)
    {
        std::cout << "\n FAILURE: LevelSetShape bounds should be defined at construction!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    return bounding_box_;
}
//=================================================================================================//
Vecd LevelSetShape::findLevelSetGradient(const Vecd &probe_point)
{
    return level_set_.probeLevelSetGradient(probe_point);
}
//=================================================================================================//
Real LevelSetShape::computeKernelIntegral(const Vecd &probe_point, Real h_ratio)
{
    return level_set_.probeKernelIntegral(probe_point, h_ratio);
}
//=================================================================================================//
Vecd LevelSetShape::computeKernelGradientIntegral(const Vecd &probe_point, Real h_ratio)
{
    return level_set_.probeKernelGradientIntegral(probe_point, h_ratio);
}
//=================================================================================================//
} // namespace SPH
