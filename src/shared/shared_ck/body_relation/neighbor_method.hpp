#ifndef NEIGHBOR_METHOD_HPP
#define NEIGHBOR_METHOD_HPP

#include "neighbor_method.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
Neighbor<Base>::SmoothingKernel::SmoothingKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : KernelTabulatedCK(*encloser.base_kernel_),
      src_pos_(encloser.dv_src_pos_->DelegatedData(ex_policy)),
      tar_pos_(encloser.dv_tar_pos_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline Real Neighbor<Base>::SmoothingKernel::W(
    const Real &inv_h_squared, const Vec2d &displacement, const Real &inv_h) const
{
    return inv_h_squared * dimension_factor_2D_ * normalized_W((displacement * inv_h).norm());
}
//=================================================================================================//
inline Real Neighbor<Base>::SmoothingKernel::W(
    const Real &inv_h_cubed, const Vec3d &displacement, const Real &inv_h) const
{
    return inv_h_cubed * dimension_factor_3D_ * normalized_W((displacement * inv_h).norm());
}
//=================================================================================================//
inline Real Neighbor<Base>::SmoothingKernel::dW(
    const Real &inv_h_cubed, const Vec2d &displacement, const Real &inv_h) const
{
    return inv_h_cubed * dimension_factor_2D_ * normalized_dW((displacement * inv_h).norm());
}
//=================================================================================================//
inline Real Neighbor<Base>::SmoothingKernel::dW(
    const Real &inv_h_fourth, const Vec3d &displacement, const Real &inv_h) const
{
    return inv_h_fourth * dimension_factor_3D_ * normalized_dW((displacement * inv_h).norm());
}
//=================================================================================================//
inline Real Neighbor<Base>::SmoothingKernel::d2W(
    const Real &inv_h_fourth, const Vec2d &displacement, const Real &inv_h) const
{
    return inv_h_fourth * dimension_factor_2D_ * normalized_d2W((displacement * inv_h).norm());
}
//=================================================================================================//
inline Real Neighbor<Base>::SmoothingKernel::d2W(
    const Real &inv_h_fifth, const Vec3d &displacement, const Real &inv_h) const
{
    return inv_h_fifth * dimension_factor_3D_ * normalized_d2W((displacement * inv_h).norm());
}
//=================================================================================================//
template <typename SourceIdentifier, typename TargetIdentifier>
Neighbor<SPHAdaptation, SPHAdaptation>::Neighbor(
    SourceIdentifier &source_identifier, TargetIdentifier &target_identifier,
    DiscreteVariable<Vecd> *dv_source_pos, DiscreteVariable<Vecd> *dv_target_pos)
    : Neighbor<Base>(
          source_identifier.getSPHAdaptation().getKernelPtr(), dv_source_pos, dv_target_pos)
{
    Real source_h = source_identifier.getSPHAdaptation().ReferenceSmoothingLength();
    Real target_h = target_identifier.getSPHAdaptation().ReferenceSmoothingLength();
    inv_h_ = 1.0 / SMAX(source_h, target_h);
    search_depth_ = static_cast<int>(std::ceil((source_h - Eps) / target_h));
    search_box_ = BoundingBoxi(Arrayi::Constant(search_depth_));
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::SmoothingKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseKernel(ex_policy, encloser),
      inv_h_(encloser.inv_h_), inv_h_squared_(inv_h_ * inv_h_),
      inv_h_cubed_(inv_h_squared_ * inv_h_), inv_h_fourth_(inv_h_cubed_ * inv_h_),
      inv_h_fifth_(inv_h_fourth_ * inv_h_) {}
//=================================================================================================//
inline Real Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::
    W(const Vec2d &displacement) const
{
    return BaseKernel::W(inv_h_squared_, displacement, inv_h_);
};
//=================================================================================================//
inline Real Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::
    W(const Vec3d &displacement) const
{
    return BaseKernel::W(inv_h_cubed_, displacement, inv_h_);
};
//=================================================================================================//
inline Real Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::
    dW(const Vec2d &displacement) const
{
    return BaseKernel::dW(inv_h_cubed_, displacement, inv_h_);
};
//=================================================================================================//
inline Real Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::
    dW(const Vec3d &displacement) const
{
    return BaseKernel::dW(inv_h_fourth_, displacement, inv_h_);
};
//=================================================================================================//
inline Real Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::
    d2W(const Vec2d &displacement) const
{
    return BaseKernel::d2W(inv_h_fourth_, displacement, inv_h_);
}
//=================================================================================================//
inline Real Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::
    d2W(const Vec3d &displacement) const
{
    return BaseKernel::d2W(inv_h_fifth_, displacement, inv_h_);
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
Neighbor<SPHAdaptation, SPHAdaptation>::NeighborCriterion::NeighborCriterion(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : src_pos_(encloser.dv_src_pos_->DelegatedData(ex_policy)),
      tar_pos_(encloser.dv_tar_pos_->DelegatedData(ex_policy)),
      kernel_size_squared_(math::pow(encloser.base_kernel_->KernelSize(), 2)),
      inv_h_(encloser.inv_h_) {}
//=================================================================================================//
inline bool Neighbor<SPHAdaptation, SPHAdaptation>::NeighborCriterion::operator()(
    UnsignedInt j, UnsignedInt i) const
{
    return (inv_h_ * (src_pos_[i] - tar_pos_[j])).squaredNorm() < kernel_size_squared_;
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
Neighbor<SPHAdaptation, SPHAdaptation>::SearchBox::SearchBox(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : search_box_(encloser.search_box_) {}
//=================================================================================================//
template <typename SourceIdentifier, typename TargetIdentifier>
Neighbor<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::Neighbor(
    SourceIdentifier &source_identifier, TargetIdentifier &target_identifier,
    DiscreteVariable<Vecd> *dv_src_pos, DiscreteVariable<Vecd> *dv_tar_pos)
    : Neighbor<Base>(
          source_identifier.getSPHAdaptation().getKernelPtr(), dv_src_pos, dv_tar_pos)
{
    AdaptiveSmoothingLength &src_adaptation = source_identifier.getAdaptation();
    src_inv_h_ref_ = 1.0 / src_adaptation.ReferenceSmoothingLength();
    src_inv_h_min_ = 1.0 / src_adaptation.MinimumSmoothingLength();
    dv_src_h_ratio_ = src_adaptation.dvSmoothingLengthRatio();

    AdaptiveSmoothingLength &tar_adaptation = target_identifier.getAdaptation();
    tar_inv_h_ref_ = 1.0 / tar_adaptation.ReferenceSmoothingLength();
    tar_inv_h_min_ = 1.0 / tar_adaptation.MinimumSmoothingLength();
    dv_tar_h_ratio_ = tar_adaptation.dvSmoothingLengthRatio();
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
Neighbor<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::SmoothingKernel::SmoothingKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseKernel(encloser),
      src_inv_h_ref_(encloser.src_inv_h_ref_), tar_inv_h_ref_(encloser.tar_inv_h_ref_),
      src_h_ratio_(encloser.dv_src_h_ratio_->DelegatedData(ex_policy)),
      tar_h_ratio_(encloser.dv_tar_h_ratio_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline Real Neighbor<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::SmoothingKernel::
    W(const Vec2d &displacement, UnsignedInt i, UnsignedInt j) const
{
    Real inv_h = invH(i, j);
    return BaseKernel::W(math::pow(inv_h, 2), displacement, inv_h);
};
//=================================================================================================//
inline Real Neighbor<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::SmoothingKernel::
    W(const Vec3d &displacement, UnsignedInt i, UnsignedInt j) const
{
    Real inv_h_i = src_h_ratio_[i] * src_inv_h_ref_;
    return BaseKernel::W(math::pow(inv_h_i, 3), displacement, inv_h_i);
};
//=================================================================================================//
inline Real Neighbor<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::SmoothingKernel::
    dW(const Vec2d &displacement, UnsignedInt i, UnsignedInt j) const
{
    Real inv_h_i = src_h_ratio_[i] * src_inv_h_ref_;
    return BaseKernel::dW(math::pow(inv_h_i, 3), displacement, inv_h_i);
};
//=================================================================================================//
inline Real Neighbor<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::SmoothingKernel::
    dW(const Vec3d &displacement, UnsignedInt i, UnsignedInt j) const
{
    Real inv_h = invH(i, j);
    return BaseKernel::dW(math::pow(inv_h, 4), displacement, inv_h);
};
//=================================================================================================//
inline Real Neighbor<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::SmoothingKernel::
    d2W(const Vec2d &displacement, UnsignedInt i, UnsignedInt j) const
{
    Real inv_h = invH(i, j);
    return BaseKernel::d2W(math::pow(inv_h, 4), displacement, inv_h);
}
//=================================================================================================//
inline Real Neighbor<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::SmoothingKernel::
    d2W(const Vec3d &displacement, UnsignedInt i, UnsignedInt j) const
{
    Real inv_h = invH(i, j);
    return BaseKernel::d2W(math::pow(inv_h, 5), displacement, inv_h);
}
//=================================================================================================//
inline Real Neighbor<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::SmoothingKernel::
    invH(UnsignedInt i, UnsignedInt j) const
{
    return SMIN(src_h_ratio_[i] * src_inv_h_ref_, src_h_ratio_[j] * tar_inv_h_ref_);
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
Neighbor<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::NeighborCriterion::
    NeighborCriterion(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : src_pos_(encloser.dv_src_pos_->DelegatedData(ex_policy)),
      tar_pos_(encloser.dv_tar_pos_->DelegatedData(ex_policy)),
      kernel_size_squared_(math::pow(encloser.base_kernel_->KernelSize(), 2)),
      src_inv_h_ref_(encloser.src_inv_h_ref_),
      src_h_ratio_(encloser.dv_src_h_ratio_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline bool Neighbor<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::NeighborCriterion::
operator()(UnsignedInt i, UnsignedInt j) const
{
    Real inv_h_i = src_h_ratio_[i] * src_inv_h_ref_;
    return (inv_h_i * (src_pos_[i] - tar_pos_[j])).squaredNorm() < kernel_size_squared_;
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
Neighbor<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::ReverseNeighborCriterion::
    ReverseNeighborCriterion(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : src_pos_(encloser.dv_src_pos_->DelegatedData(ex_policy)),
      tar_pos_(encloser.dv_tar_pos_->DelegatedData(ex_policy)),
      kernel_size_squared_(math::pow(encloser.base_kernel_->KernelSize(), 2)),
      tar_inv_h_ref_(encloser.tar_inv_h_ref_),
      tar_h_ratio_(encloser.dv_tar_h_ratio_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline bool Neighbor<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::
    ReverseNeighborCriterion::operator()(UnsignedInt i, UnsignedInt j) const
{
    Real inv_h_j = tar_h_ratio_[j] * tar_inv_h_ref_;
    return (inv_h_j * (src_pos_[i] - tar_pos_[j])).squaredNorm() < kernel_size_squared_;
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
Neighbor<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::SearchBox::SearchBox(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : src_inv_h_ref_(encloser.src_inv_h_ref_), tar_inv_h_min_(encloser.tar_inv_h_min_),
      src_h_ratio_(encloser.dv_src_h_ratio_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline BoundingBoxi Neighbor<AdaptiveSmoothingLength, AdaptiveSmoothingLength>::
    SearchBox::operator()(UnsignedInt i) const
{
    Real search_depth =
        static_cast<int>(std::ceil((src_h_ratio_[i] * src_inv_h_ref_ - Eps) * tar_inv_h_min_));
    return BoundingBoxi(Arrayi::Constant(search_depth));
}
//=================================================================================================//
} // namespace SPH
#endif // NEIGHBOR_METHOD_HPP
