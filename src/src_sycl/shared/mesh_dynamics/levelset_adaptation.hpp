#ifndef LEVELSET_ADAPTATION_HPP
#define LEVELSET_ADAPTATION_HPP

#include "adaptation.h"
#include "level_set_sycl.hpp"
namespace SPH
{
//=================================================================================================//
UniquePtr<MultilevelLevelSet> SPHAdaptation::createLevelSet(const ParallelDevicePolicy &par_device,
                                                            Shape &shape, Real refinement_ratio)
{
    // estimate the required mesh levels
    int total_levels = (int)log10(MinimumDimension(shape.getBounds()) / ReferenceSpacing()) + 2;
    Real coarsest_spacing = ReferenceSpacing() * pow(2.0, total_levels - 1);
    MultilevelLevelSet coarser_level_sets(par_device, shape.getBounds(), coarsest_spacing / refinement_ratio,
    total_levels - 1, shape, *this);
    // return the finest level set only
    return makeUnique<MultilevelLevelSet>(par_device, shape.getBounds(), coarser_level_sets.getMeshLevels().back(), shape, *this);
}
//=================================================================================================//
UniquePtr<MultilevelLevelSet> ParticleWithLocalRefinement::createLevelSet(const ParallelDevicePolicy &par_device,
                                                                          Shape &shape, Real refinement_ratio)
{
    return makeUnique<MultilevelLevelSet>(par_device, shape.getBounds(), ReferenceSpacing() / refinement_ratio,
    getLevelSetTotalLevel(), shape, *this);
}
//=================================================================================================//
}

#endif