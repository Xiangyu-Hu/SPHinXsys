#include "level_set_shape.h"

#include "all_io.h"
#include "base_body.h"
#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
LevelSetShape::
    LevelSetShape(Shape &shape, SharedPtr<SPHAdaptation> sph_adaptation,
                  Real refinement_ratio, UsageType usage_type)
    : LevelSetShape(shape.getBounds(), shape, sph_adaptation, refinement_ratio)
{
    finishInitialization(execution::par_host, usage_type);
}
//=================================================================================================//
LevelSetShape::LevelSetShape(SPHBody &sph_body, Shape &shape,
                             Real refinement_ratio, UsageType usage_type)
    : LevelSetShape(shape.getBounds(), sph_body, shape, refinement_ratio)
{
    finishInitialization(execution::par_host, usage_type);
}
//=================================================================================================//
LevelSetShape::LevelSetShape(BoundingBoxd bounding_box, Shape &shape,
                             SharedPtr<SPHAdaptation> sph_adaptation, Real refinement_ratio)
    : Shape(shape.getName()), sph_adaptation_(sph_adaptation),
      level_set_(*level_set_keeper_.movePtr(sph_adaptation->createLevelSet(shape, refinement_ratio)))
{
    bounding_box_ = shape.getBounds();
    is_bounds_found_ = true;
}
//=================================================================================================//
LevelSetShape::LevelSetShape(BoundingBoxd bounding_box, SPHBody &sph_body,
                             Shape &shape, Real refinement_ratio)
    : Shape(shape.getName()),
      level_set_(*level_set_keeper_.movePtr(
          sph_body.getSPHAdaptation().createLevelSet(shape, refinement_ratio)))
{
    bounding_box_ = shape.getBounds();
    is_bounds_found_ = true;
}
//=================================================================================================//
LevelSetShape *LevelSetShape::writeLevelSet(SPHSystem &sph_system)
{
    MeshRecordingToPlt write_level_set_to_plt(sph_system, level_set_);
    write_level_set_to_plt.writeToFile(0);
    return this;
}
//=================================================================================================//
LevelSetShape *LevelSetShape::cleanLevelSet(UnsignedInt repeat_times)
{
    level_set_.cleanInterface(repeat_times);
    return this;
}
//=================================================================================================//
LevelSetShape *LevelSetShape::correctLevelSetSign()
{
    level_set_.correctTopology();
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
BoundingBoxd LevelSetShape::findBounds()
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
Matd LevelSetShape::computeKernelSecondGradientIntegral(const Vecd &probe_point, Real h_ratio)
{
    return level_set_.probeKernelSecondGradientIntegral(probe_point, h_ratio);
}
//=================================================================================================//
} // namespace SPH
