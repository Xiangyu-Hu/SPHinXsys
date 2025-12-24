#ifndef NEIGHBOR_METHOD_HPP
#define NEIGHBOR_METHOD_HPP

#include "neighbor_method.h"

namespace SPH
{
//=================================================================================================//
inline Real NeighborMethod<Base>::SmoothingKernel::W(
    const Real &inv_h_squared, const Vec2d &displacement, const Real &inv_h) const
{
    return inv_h_squared * dimension_factor_2D_ * normalized_W((displacement * inv_h).norm());
}
//=================================================================================================//
inline Real NeighborMethod<Base>::SmoothingKernel::W(
    const Real &inv_h_cubed, const Vec3d &displacement, const Real &inv_h) const
{
    return inv_h_cubed * dimension_factor_3D_ * normalized_W((displacement * inv_h).norm());
}
//=================================================================================================//
inline Real NeighborMethod<Base>::SmoothingKernel::dW(
    const Real &inv_h_cubed, const Vec2d &displacement, const Real &inv_h) const
{
    return inv_h_cubed * dimension_factor_2D_ * normalized_dW((displacement * inv_h).norm());
}
//=================================================================================================//
inline Real NeighborMethod<Base>::SmoothingKernel::dW(
    const Real &inv_h_fourth, const Vec3d &displacement, const Real &inv_h) const
{
    return inv_h_fourth * dimension_factor_3D_ * normalized_dW((displacement * inv_h).norm());
}
//=================================================================================================//
inline Real NeighborMethod<Base>::SmoothingKernel::d2W(
    const Real &inv_h_fourth, const Vec2d &displacement, const Real &inv_h) const
{
    return inv_h_fourth * dimension_factor_2D_ * normalized_d2W((displacement * inv_h).norm());
}
//=================================================================================================//
inline Real NeighborMethod<Base>::SmoothingKernel::d2W(
    const Real &inv_h_fifth, const Vec3d &displacement, const Real &inv_h) const
{
    return inv_h_fifth * dimension_factor_3D_ * normalized_d2W((displacement * inv_h).norm());
}
//=================================================================================================//
template <typename SourceIdentifier, typename TargetIdentifier>
NeighborMethod<SPHAdaptation, SPHAdaptation>::NeighborMethod(
    SourceIdentifier &source_identifier, TargetIdentifier &target_identifier)
    : NeighborMethod<Base>(source_identifier.getSPHAdaptation().getKernelPtr())
{
    Real source_h = source_identifier.getSPHAdaptation().ReferenceSmoothingLength();
    Real target_h = target_identifier.getSPHAdaptation().ReferenceSmoothingLength();
    inv_h_ = 1.0 / SMAX(source_h, target_h);
    search_depth_ = static_cast<int>(std::ceil((source_h - Eps) / target_h));
    search_box_ = BoundingBoxi(Arrayi::Constant(search_depth_));
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
NeighborMethod<SPHAdaptation, SPHAdaptation>::SmoothingKernel::SmoothingKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseKernel(encloser),
      inv_h_(encloser.inv_h_), inv_h_squared_(inv_h_ * inv_h_),
      inv_h_cubed_(inv_h_squared_ * inv_h_), inv_h_fourth_(inv_h_cubed_ * inv_h_),
      inv_h_fifth_(inv_h_fourth_ * inv_h_) {}
//=================================================================================================//
inline Real NeighborMethod<SPHAdaptation, SPHAdaptation>::SmoothingKernel::W(
    const Vec2d &displacement, UnsignedInt, UnsignedInt) const
{
    return BaseKernel::W(inv_h_squared_, displacement, inv_h_);
};
//=================================================================================================//
inline Real NeighborMethod<SPHAdaptation, SPHAdaptation>::SmoothingKernel::W(
    const Vec3d &displacement, UnsignedInt, UnsignedInt) const
{
    return BaseKernel::W(inv_h_cubed_, displacement, inv_h_);
};
//=================================================================================================//
inline Real NeighborMethod<SPHAdaptation, SPHAdaptation>::SmoothingKernel::dW(
    const Vec2d &displacement, UnsignedInt, UnsignedInt) const
{
    return BaseKernel::dW(inv_h_cubed_, displacement, inv_h_);
};
//=================================================================================================//
inline Real NeighborMethod<SPHAdaptation, SPHAdaptation>::SmoothingKernel::dW(
    const Vec3d &displacement, UnsignedInt, UnsignedInt) const
{
    return BaseKernel::dW(inv_h_fourth_, displacement, inv_h_);
};
//=================================================================================================//
inline Real NeighborMethod<SPHAdaptation, SPHAdaptation>::SmoothingKernel::d2W(
    const Vec2d &displacement, UnsignedInt, UnsignedInt) const
{
    return BaseKernel::d2W(inv_h_fourth_, displacement, inv_h_);
}
//=================================================================================================//
inline Real NeighborMethod<SPHAdaptation, SPHAdaptation>::SmoothingKernel::d2W(
    const Vec3d &displacement, UnsignedInt, UnsignedInt) const
{
    return BaseKernel::d2W(inv_h_fifth_, displacement, inv_h_);
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
NeighborMethod<SPHAdaptation, SPHAdaptation>::NeighborCriterion::NeighborCriterion(
    const ExecutionPolicy &ex_policy, EncloserType &encloser,
    DiscreteVariable<Vecd> *dv_source_pos, DiscreteVariable<Vecd> *dv_target_pos)
    : source_pos_(dv_source_pos->DelegatedData(ex_policy)),
      target_pos_(dv_target_pos->DelegatedData(ex_policy)),
      kernel_size_squared_(math::pow(encloser.base_kernel_->KernelSize(), 2)),
      inv_h_(encloser.inv_h_) {}
//=================================================================================================//
inline bool NeighborMethod<SPHAdaptation, SPHAdaptation>::NeighborCriterion::operator()(
    UnsignedInt i, UnsignedInt j) const
{
    return (inv_h_ * (source_pos_[i] - target_pos_[j])).squaredNorm() < kernel_size_squared_;
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
NeighborMethod<SPHAdaptation, SPHAdaptation>::SearchBox::SearchBox(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : search_box_(encloser.search_box_) {}
//=================================================================================================//
template <typename SourceIdentifier, typename TargetIdentifier>
NeighborMethod<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::NeighborMethod(
    SourceIdentifier &source_identifier, TargetIdentifier &target_identifier)
    : NeighborMethod<Base>(*source_identifier.getSPHAdaptation().getKernel())
{
}
//=================================================================================================//
inline Real NeighborMethod<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::SmoothingKernel::W(
    const Vec2d &displacement, UnsignedInt i, UnsignedInt j) const
{
    Real inv_h = invH(i, j);
    return BaseKernel::W(math::pow(inv_h, 2), displacement, inv_h);
};
//=================================================================================================//
inline Real NeighborMethod<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::SmoothingKernel::W(
    const Vec3d &displacement, UnsignedInt i, UnsignedInt j) const
{
    Real inv_h_i = src_h_ratio_[i] * src_inv_h_ref_;
    return BaseKernel::W(math::pow(inv_h_i, 3), displacement, inv_h_i);
};
//=================================================================================================//
inline Real NeighborMethod<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::SmoothingKernel::dW(
    const Vec2d &displacement, UnsignedInt i, UnsignedInt j) const
{
    Real inv_h_i = src_h_ratio_[i] * src_inv_h_ref_;
    return BaseKernel::dW(math::pow(inv_h_i, 3), displacement, inv_h_i);
};
//=================================================================================================//
inline Real NeighborMethod<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::SmoothingKernel::dW(
    const Vec3d &displacement, UnsignedInt i, UnsignedInt j) const
{
    Real inv_h = invH(i, j);
    return BaseKernel::dW(math::pow(inv_h, 4), displacement, inv_h);
};
//=================================================================================================//
inline Real NeighborMethod<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::SmoothingKernel::d2W(
    const Vec2d &displacement, UnsignedInt i, UnsignedInt j) const
{
    Real inv_h = invH(i, j);
    return BaseKernel::d2W(math::pow(inv_h, 4), displacement, inv_h);
}
//=================================================================================================//
inline Real NeighborMethod<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::SmoothingKernel::d2W(
    const Vec3d &displacement, UnsignedInt i, UnsignedInt j) const
{
    Real inv_h = invH(i, j);
    return BaseKernel::d2W(math::pow(inv_h, 5), displacement, inv_h);
}
//=================================================================================================//
inline Real NeighborMethod<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::SmoothingKernel::invH(
    UnsignedInt i, UnsignedInt j) const
{
    return SMIN(src_h_ratio_[i] * src_inv_h_ref_, src_h_ratio_[j] * tar_inv_h_ref_);
}
//=================================================================================================//
} // namespace SPH
#endif // NEIGHBOR_METHOD_HPP
