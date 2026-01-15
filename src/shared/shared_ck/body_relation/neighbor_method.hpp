#ifndef NEIGHBOR_METHOD_HPP
#define NEIGHBOR_METHOD_HPP

#include "neighbor_method.h"

#include "adaptation.h"

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
inline Real Neighbor<Base>::SmoothingKernel::W0(
    const Real &inv_h_squared, const Vec2d &) const
{
    return inv_h_squared * dimension_factor_2D_ * normalized_W(0);
}
//=================================================================================================//
inline Real Neighbor<Base>::SmoothingKernel::W0(
    const Real &inv_h_cubed, const Vec3d &) const
{
    return inv_h_cubed * dimension_factor_3D_ * normalized_W(0);
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
    Real src_h = source_identifier.getSPHAdaptation().ReferenceSmoothingLength();
    Real tar_h = target_identifier.getSPHAdaptation().ReferenceSmoothingLength();
    src_inv_h_ = 1.0 / src_h;
    inv_h_ = 1.0 / SMAX(src_h, tar_h);
    search_depth_ = static_cast<int>(std::ceil((src_h - Eps) / tar_h));
    search_box_ = BoundingBoxi(Arrayi::Constant(search_depth_));
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::SmoothingKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseKernel(ex_policy, encloser),
      src_inv_h_(encloser.src_inv_h_), src_inv_h_squared_(src_inv_h_ * src_inv_h_),
      src_inv_h_cubed_(src_inv_h_squared_ * src_inv_h_),
      inv_h_(encloser.inv_h_), inv_h_squared_(inv_h_ * inv_h_),
      inv_h_cubed_(inv_h_squared_ * inv_h_), inv_h_fourth_(inv_h_cubed_ * inv_h_),
      inv_h_fifth_(inv_h_fourth_ * inv_h_) {}
//=================================================================================================//
inline Real Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::
    W0(UnsignedInt i, const Vec2d &zero) const
{
    return BaseKernel::W0(src_inv_h_squared_, zero);
}
//=================================================================================================//
inline Real Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::
    W0(UnsignedInt i, const Vec3d &zero) const
{
    return BaseKernel::W0(src_inv_h_cubed_, zero);
}
//=================================================================================================//
inline Real Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::
    W(const Vec2d &displacement) const
{
    return BaseKernel::W(src_inv_h_squared_, displacement, src_inv_h_);
};
//=================================================================================================//
inline Real Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::
    W(const Vec3d &displacement) const
{
    return BaseKernel::W(src_inv_h_cubed_, displacement, src_inv_h_);
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
template <class SourceAdaptationType, class TargetAdaptationType>
template <typename SourceIdentifier, typename TargetIdentifier>
Neighbor<SourceAdaptationType, TargetAdaptationType>::Neighbor(
    SourceIdentifier &source_identifier, TargetIdentifier &target_identifier,
    DiscreteVariable<Vecd> *dv_src_pos, DiscreteVariable<Vecd> *dv_tar_pos)
    : Neighbor<Base>(source_identifier.getSPHAdaptation().getKernelPtr(), dv_src_pos, dv_tar_pos),
      src_adaptation_(DynamicCast<SourceAdaptationType>(this, source_identifier.getSPHAdaptation())),
      tar_adaptation_(DynamicCast<TargetAdaptationType>(this, target_identifier.getSPHAdaptation()))
{
    src_inv_h_ref_ = 1.0 / src_adaptation_.ReferenceSmoothingLength();
    src_inv_h_min_ = 1.0 / src_adaptation_.MinimumSmoothingLength();

    tar_inv_h_ref_ = 1.0 / tar_adaptation_.ReferenceSmoothingLength();
    tar_inv_h_min_ = 1.0 / tar_adaptation_.MinimumSmoothingLength();
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
template <class ExecutionPolicy, class EncloserType>
Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::SmoothingKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseKernel(ex_policy, encloser),
      src_inv_h_ref_(encloser.src_inv_h_ref_), tar_inv_h_ref_(encloser.tar_inv_h_ref_),
      src_h_ratio_(ex_policy, encloser.src_adaptation_),
      tar_h_ratio_(ex_policy, encloser.tar_adaptation_) {}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    W(const Vec2d &displacement, UnsignedInt i, UnsignedInt j) const
{
    Real inv_h_i = src_h_ratio_(i) * src_inv_h_ref_;
    return BaseKernel::W(math::pow(inv_h_i, 2), displacement, inv_h_i);
};
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    W(const Vec3d &displacement, UnsignedInt i, UnsignedInt j) const
{
    Real inv_h_i = src_h_ratio_(i) * src_inv_h_ref_;
    return BaseKernel::W(math::pow(inv_h_i, 3), displacement, inv_h_i);
};
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    dW(const Vec2d &displacement, UnsignedInt i, UnsignedInt j) const
{
    Real inv_h = invSmoothingLength(i, j);
    return BaseKernel::dW(math::pow(inv_h, 3), displacement, inv_h);
};
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    dW(const Vec3d &displacement, UnsignedInt i, UnsignedInt j) const
{
    Real inv_h = invSmoothingLength(i, j);
    return BaseKernel::dW(math::pow(inv_h, 4), displacement, inv_h);
};
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    d2W(const Vec2d &displacement, UnsignedInt i, UnsignedInt j) const
{
    Real inv_h = invSmoothingLength(i, j);
    return BaseKernel::d2W(math::pow(inv_h, 4), displacement, inv_h);
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    d2W(const Vec3d &displacement, UnsignedInt i, UnsignedInt j) const
{
    Real inv_h = invSmoothingLength(i, j);
    return BaseKernel::d2W(math::pow(inv_h, 5), displacement, inv_h);
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    invSmoothingLength(UnsignedInt i, UnsignedInt j) const
{
    return SMIN(src_h_ratio_(i) * src_inv_h_ref_, tar_h_ratio_(j) * tar_inv_h_ref_);
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
template <class ExecutionPolicy, class EncloserType>
Neighbor<SourceAdaptationType, TargetAdaptationType>::NeighborCriterion::
    NeighborCriterion(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : src_pos_(encloser.dv_src_pos_->DelegatedData(ex_policy)),
      tar_pos_(encloser.dv_tar_pos_->DelegatedData(ex_policy)),
      kernel_size_squared_(math::pow(encloser.base_kernel_->KernelSize(), 2)),
      src_inv_h_ref_(encloser.src_inv_h_ref_),
      src_h_ratio_(ex_policy, encloser.src_adaptation_) {}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
inline bool Neighbor<SourceAdaptationType, TargetAdaptationType>::NeighborCriterion::
operator()(UnsignedInt j, UnsignedInt i) const
{
    Real inv_h_i = src_h_ratio_(i) * src_inv_h_ref_;
    return (inv_h_i * (src_pos_[i] - tar_pos_[j])).squaredNorm() < kernel_size_squared_;
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
template <class ExecutionPolicy, class EncloserType>
Neighbor<SourceAdaptationType, TargetAdaptationType>::ReverseNeighborCriterion::
    ReverseNeighborCriterion(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : src_pos_(encloser.dv_src_pos_->DelegatedData(ex_policy)),
      tar_pos_(encloser.dv_tar_pos_->DelegatedData(ex_policy)),
      kernel_size_squared_(math::pow(encloser.base_kernel_->KernelSize(), 2)),
      tar_inv_h_ref_(encloser.tar_inv_h_ref_),
      tar_h_ratio_(ex_policy, encloser.tar_adaptation_) {}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
bool Neighbor<SourceAdaptationType, TargetAdaptationType>::
    ReverseNeighborCriterion::operator()(UnsignedInt i, UnsignedInt j) const
{
    Real inv_h_j = tar_h_ratio_(j) * tar_inv_h_ref_;
    return (inv_h_j * (src_pos_[i] - tar_pos_[j])).squaredNorm() < kernel_size_squared_;
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
template <class ExecutionPolicy, class EncloserType>
Neighbor<SourceAdaptationType, TargetAdaptationType>::SearchBox::SearchBox(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : src_inv_h_ref_(encloser.src_inv_h_ref_), tar_inv_h_min_(encloser.tar_inv_h_min_),
      src_h_ratio_(ex_policy, encloser.src_adaptation_) {}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
inline BoundingBoxi Neighbor<SourceAdaptationType, TargetAdaptationType>::
    SearchBox::operator()(UnsignedInt i) const
{
    Real src_h = 1.0 / (src_h_ratio_(i) * src_inv_h_ref_);
    Real search_depth = static_cast<int>(std::ceil((src_h - Eps) * tar_inv_h_min_));
    return BoundingBoxi(Arrayi::Constant(search_depth));
}
//=================================================================================================//
} // namespace SPH
#endif // NEIGHBOR_METHOD_HPP
