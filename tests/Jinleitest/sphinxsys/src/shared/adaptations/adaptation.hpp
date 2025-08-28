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
} // namespace SPH
#endif // ADAPTATION_HPP