#include "level_set_shape.h"

#include "all_io.h"
#include "sph_system.h"
#include "base_body.h"

namespace SPH
{
//=================================================================================================//
LevelSetShape::LevelSetShape(
    SPHBody &sph_body, Shape &shape, Real refinement, UsageType usage_type)
    : LevelSetShape(sph_body.getSPHSystem(), sph_body.getSPHAdaptation(), shape, refinement)
{
    finishInitialization(execution::par_host, usage_type);
}
//=================================================================================================//
LevelSetShape::LevelSetShape(
    SPHSystem &sph_system, const SPHAdaptation &sph_adaptation, Shape &shape, Real refinement)
    : Shape(shape.getName()), sph_system_(sph_system),
      level_set_(*level_set_keeper_.movePtr(sph_adaptation.createLevelSet(shape, refinement)))
{
    bounding_box_ = shape.getBounds();
    is_bounds_found_ = true;
}
//=================================================================================================//
LevelSetShape *LevelSetShape::writeLevelSet()
{
    MeshRecordingToPlt write_level_set_to_plt(sph_system_, level_set_);
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
