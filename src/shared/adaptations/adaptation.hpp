#ifndef ADAPTATION_HPP
#define ADAPTATION_HPP

#include "adaptation.h"
#include "base_body.h"

namespace SPH
{
//=================================================================================================//
template <class MeshType, typename... Args>
MeshType SPHAdaptation::createBackGroundMesh(SPHBody &sph_body, Args &&...args)
{
    return MeshType(sph_body.getSPHSystemBounds(), kernel_ptr_->CutOffRadius(), 2,
                    std::forward<Args>(args)...);
}
//=================================================================================================//
template <class ExecutionPolicy, class EnclosureType>
ParticleWithLocalRefinement::AdaptationKernel::
    AdaptationKernel(const ExecutionPolicy &ex_policy, EnclosureType &enclosure)
    : h_ratio_(enclosure.dv_h_ratio_->DelegatedData(ex_policy)),
      level_(enclosure.dv_level_->DelegatedData(ex_policy)) {}
//=================================================================================================//
} // namespace SPH
#endif // ADAPTATION_HPP