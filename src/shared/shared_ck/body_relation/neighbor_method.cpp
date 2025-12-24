#include "neighbor_method.h"

namespace SPH
{
//=================================================================================================//
NeighborMethod<SPHAdaptation, SPHAdaptation>::NeighborMethod(
    SharedPtr<Kernel> base_kernel, Real h, Real search_increment)
    : NeighborMethod<Base>(base_kernel), inv_h_(1.0 / h),
      search_depth_(static_cast<int>(std::ceil((h - Eps) / search_increment))),
      search_box_(BoundingBoxi(Arrayi::Constant(search_depth_))) {}
//=================================================================================================//
} // namespace SPH
