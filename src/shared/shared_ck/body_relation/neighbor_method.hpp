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
inline Vecd Neighbor<SPHAdaptation, SPHAdaptation>::
    SmoothingKernel::nablaW_ij(UnsignedInt i, UnsignedInt j) const
{
    Vecd disp = vec_r_ij(i, j);
    return dW(disp) * disp.normalized();
}
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
Neighbor<SPHAdaptation, SPHAdaptation>::CutOff::CutOff(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : cut_off_(encloser.base_kernel_->KernelSize() * Vecd::Ones() / encloser.src_inv_h_) {}
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
inline Vecd Neighbor<SourceAdaptationType, TargetAdaptationType>::
    SmoothingKernel::nablaW_ij(UnsignedInt i, UnsignedInt j) const
{
    auto measure = getTransformedMeasure(i, j);
    const Vecd &disp_transform = std::get<0>(measure);
    if (std::get<2>(measure))
    {
        return src_h_ratio_.KernelTransform(i) * dW(disp_transform, std::get<1>(measure)) *
               src_h_ratio_.GradientTransform(i) * disp_transform.normalized();
    };

    return tar_h_ratio_.KernelTransform(j) * dW(disp_transform, std::get<1>(measure)) *
           tar_h_ratio_.GradientTransform(j) * disp_transform.normalized();
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    W_ij(UnsignedInt i, UnsignedInt j) const
{
    auto measure = getTransformedMeasure(i, j);
    Real kernel_transform = std::get<2>(measure) ? src_h_ratio_.KernelTransform(i) : tar_h_ratio_.KernelTransform(j);
    return kernel_transform * W(std::get<0>(measure), std::get<1>(measure));
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    dW_ij(UnsignedInt i, UnsignedInt j) const
{
    auto measure = getTransformedMeasure(i, j);
    Real kernel_transform = std::get<2>(measure) ? src_h_ratio_.KernelTransform(i) : tar_h_ratio_.KernelTransform(j);
    return kernel_transform * dW(std::get<0>(measure), std::get<1>(measure));
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    W0(UnsignedInt i, const Vec2d &zero) const
{
    Real inv_h_i = src_h_ratio_(i) * src_inv_h_ref_;
    return src_h_ratio_.KernelTransform(i) * BaseKernel::W0(math::pow(inv_h_i, 2), zero);
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    W0(UnsignedInt i, const Vec3d &zero) const
{
    Real inv_h_i = src_h_ratio_(i) * src_inv_h_ref_;
    return src_h_ratio_.KernelTransform(i) * BaseKernel::W0(math::pow(inv_h_i, 3), zero);
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    W(const Vec2d &disp_transform, Real inv_h) const
{
    return BaseKernel::W(math::pow(inv_h, 2), disp_transform, inv_h);
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    W(const Vec3d &disp_transform, Real inv_h) const
{
    return BaseKernel::W(math::pow(inv_h, 3), disp_transform, inv_h);
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    dW(const Vec2d &disp_transform, Real inv_h) const
{
    return BaseKernel::dW(math::pow(inv_h, 3), disp_transform, inv_h);
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    dW(const Vec3d &disp_transform, Real inv_h) const
{
    return BaseKernel::dW(math::pow(inv_h, 4), disp_transform, inv_h);
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    d2W(const Vec2d &disp_transform, Real inv_h) const
{
    return BaseKernel::d2W(math::pow(inv_h, 4), disp_transform, inv_h);
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    d2W(const Vec3d &disp_transform, Real inv_h) const
{
    return BaseKernel::d2W(math::pow(inv_h, 5), disp_transform, inv_h);
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
std::tuple<Vecd, Real, bool> Neighbor<SourceAdaptationType, TargetAdaptationType>::
    SmoothingKernel::getTransformedMeasure(UnsignedInt i, UnsignedInt j) const
{
    Real inv_h_i = src_h_ratio_(i) * src_inv_h_ref_;
    Real inv_h_j = tar_h_ratio_(j) * tar_inv_h_ref_;
    Vecd disp = vec_r_ij(i, j);
    Vecd disp_i = src_h_ratio_.transform(disp, i);
    Vecd disp_j = tar_h_ratio_.transform(disp, j);
    return disp_i.squaredNorm() * inv_h_i * inv_h_i < disp_j.squaredNorm() * inv_h_j * inv_h_j
               ? std::tuple<Vecd, Real, bool>(disp_i, inv_h_i, true)
               : std::tuple<Vecd, Real, bool>(disp_j, inv_h_j, false);
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
    Vecd displacement = src_pos_[i] - tar_pos_[j];
    return (inv_h_i * src_h_ratio_.transform(displacement, i)).squaredNorm() < kernel_size_squared_;
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
    Vecd displacement = src_pos_[i] - tar_pos_[j];
    return (inv_h_j * tar_h_ratio_.transform(displacement, j)).squaredNorm() < kernel_size_squared_;
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
template <class ExecutionPolicy, class EncloserType>
Neighbor<SourceAdaptationType, TargetAdaptationType>::CutOff::CutOff(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : kernel_size_(encloser.base_kernel_->KernelSize()),
      src_inv_h_ref_(encloser.src_inv_h_ref_),
      src_h_ratio_(ex_policy, encloser.src_adaptation_) {}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
inline Vecd Neighbor<SourceAdaptationType, TargetAdaptationType>::
    CutOff::operator()(UnsignedInt i) const
{
    return kernel_size_ * src_h_ratio_.inverseTransform(Vecd::Ones(), i) / (src_h_ratio_(i) * src_inv_h_ref_);
}
//=================================================================================================//
} // namespace SPH
#endif // NEIGHBOR_METHOD_HPP
