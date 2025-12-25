#include "neighbor_method.hpp"

namespace SPH
{
//=================================================================================================//
NeighborMethod<Base>::NeighborMethod(
    SharedPtr<Kernel> base_kernel,
    DiscreteVariable<Vecd> *dv_src_pos, DiscreteVariable<Vecd> *dv_tar_pos)
    : base_kernel_(base_kernel), dv_src_pos_(dv_src_pos), dv_tar_pos_(dv_tar_pos) {}
//=================================================================================================//
NeighborMethod<Base>::NeighborMethod(SharedPtr<Kernel> base_kernel)
    : base_kernel_(base_kernel), dv_src_pos_(nullptr), dv_tar_pos_(nullptr) {}
//=================================================================================================//
NeighborMethod<Base>::SmoothingKernel::SmoothingKernel(NeighborMethod<Base> &encloser)
    : KernelTabulatedCK(*encloser.base_kernel_), src_pos_(nullptr), tar_pos_(nullptr) {}
//=================================================================================================//
NeighborMethod<SPHAdaptation, SPHAdaptation>::NeighborMethod(
    SharedPtr<Kernel> base_kernel, Real h, Real search_increment)
    : NeighborMethod<Base>(base_kernel), inv_h_(1.0 / h),
      search_depth_(static_cast<int>(std::ceil((h - Eps) / search_increment))),
      search_box_(BoundingBoxi(Arrayi::Constant(search_depth_))) {}
//=================================================================================================//
NeighborMethod<SPHAdaptation, SPHAdaptation>::SmoothingKernel::SmoothingKernel(
    NeighborMethod<SPHAdaptation, SPHAdaptation> &encloser)
    : BaseKernel(encloser), inv_h_(encloser.inv_h_),
      inv_h_squared_(inv_h_ * inv_h_), inv_h_cubed_(inv_h_squared_ * inv_h_),
      inv_h_fourth_(inv_h_cubed_ * inv_h_), inv_h_fifth_(inv_h_fourth_ * inv_h_) {}
//=================================================================================================//
} // namespace SPH
