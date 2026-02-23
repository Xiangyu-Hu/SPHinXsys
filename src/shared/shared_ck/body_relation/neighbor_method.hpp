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
inline Real Neighbor<Base>::SmoothingKernel::W2D(
    const Real &inv_h_squared, const Real &scaled_r) const
{
    return inv_h_squared * dimension_factor_2D_ * normalized_W(scaled_r);
}
//=================================================================================================//
inline Real Neighbor<Base>::SmoothingKernel::W3D(
    const Real &inv_h_cubed, const Real &scaled_r) const
{
    return inv_h_cubed * dimension_factor_3D_ * normalized_W(scaled_r);
}
//=================================================================================================//
inline Real Neighbor<Base>::SmoothingKernel::W02D(
    const Real &inv_h_squared) const
{
    return inv_h_squared * dimension_factor_2D_ * normalized_W(0);
}
//=================================================================================================//
inline Real Neighbor<Base>::SmoothingKernel::W03D(
    const Real &inv_h_cubed) const
{
    return inv_h_cubed * dimension_factor_3D_ * normalized_W(0);
}
//=================================================================================================//
inline Real Neighbor<Base>::SmoothingKernel::dW2D(
    const Real &inv_h_cubed, const Real &scaled_r) const
{
    return inv_h_cubed * dimension_factor_2D_ * normalized_dW(scaled_r);
}
//=================================================================================================//
inline Real Neighbor<Base>::SmoothingKernel::dW3D(
    const Real &inv_h_fourth, const Real &scaled_r) const
{
    return inv_h_fourth * dimension_factor_3D_ * normalized_dW(scaled_r);
}
//=================================================================================================//
inline Real Neighbor<Base>::SmoothingKernel::d2W2D(
    const Real &inv_h_fourth, const Real &scaled_r) const
{
    return inv_h_fourth * dimension_factor_2D_ * normalized_d2W(scaled_r);
}
//=================================================================================================//
inline Real Neighbor<Base>::SmoothingKernel::d2W3D(
    const Real &inv_h_fifth, const Real &scaled_r) const
{
    return inv_h_fifth * dimension_factor_3D_ * normalized_d2W(scaled_r);
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
inline Vecd Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::
    nablaW_ij(UnsignedInt i, UnsignedInt j) const
{
    Vecd disp = vec_r_ij(i, j);
    return dW(disp) * disp.normalized();
}
//=================================================================================================//
inline Real Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::
    W0(UnsignedInt i, const Vec2d &zero) const
{
    return BaseKernel::W02D(src_inv_h_squared_);
}
//=================================================================================================//
inline Real Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::
    W0(UnsignedInt i, const Vec3d &zero) const
{
    return BaseKernel::W03D(src_inv_h_cubed_);
}
//=================================================================================================//
inline Real Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::
    W(const Vec2d &displacement) const
{
    return BaseKernel::W2D(inv_h_squared_, displacement.norm() * inv_h_);
};
//=================================================================================================//
inline Real Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::
    W(const Vec3d &displacement) const
{
    return BaseKernel::W3D(inv_h_cubed_, displacement.norm() * inv_h_);
};
//=================================================================================================//
inline Real Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::
    dW(const Vec2d &displacement) const
{
    return BaseKernel::dW2D(inv_h_cubed_, displacement.norm() * inv_h_);
};
//=================================================================================================//
inline Real Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::
    dW(const Vec3d &displacement) const
{
    return BaseKernel::dW3D(inv_h_fourth_, displacement.norm() * inv_h_);
};
//=================================================================================================//
inline Real Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::
    d2W(const Vec2d &displacement) const
{
    return BaseKernel::d2W2D(inv_h_fourth_, displacement.norm() * inv_h_);
}
//=================================================================================================//
inline Real Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel::
    d2W(const Vec3d &displacement) const
{
    return BaseKernel::d2W3D(inv_h_fifth_, displacement.norm() * inv_h_);
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
    Vecd disp = vec_r_ij(i, j);
    if constexpr (!std::is_base_of_v<AnisotropicAdaptation, SourceAdaptationType> &&
                  !std::is_base_of_v<AnisotropicAdaptation, TargetAdaptationType>)
    {
        return dW(disp, invSmoothingLength(i, j)) * disp.normalized();
    }
    else
    {
        return selectKernelFunction(
            i, j, disp, src_h_ratio_, tar_h_ratio_,
            [&](const Vecd &disp_transform, Real inv_h, auto &src_h_ratio)->Vecd
            { return src_h_ratio.KernelTransform(i) * dW(disp_transform, inv_h) *
                     src_h_ratio_.GradientTransform(i) * disp_transform.normalized(); },
            [&](const Vecd &disp_transform, Real inv_h, auto &tar_h_ratio)->Vecd
            { return tar_h_ratio.KernelTransform(j) * dW(disp_transform, inv_h) *
                     tar_h_ratio_.GradientTransform(j) * disp_transform.normalized(); });
    }
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    W_ij(UnsignedInt i, UnsignedInt j) const
{
    Vecd disp = vec_r_ij(i, j);
    if constexpr (!std::is_base_of_v<AnisotropicAdaptation, SourceAdaptationType> &&
                  !std::is_base_of_v<AnisotropicAdaptation, TargetAdaptationType>)
    {
        return W(disp, invSmoothingLength(i, j));
    }
    else
    {
        return selectKernelFunction(
            i, j, disp, src_h_ratio_, tar_h_ratio_,
            [&](const Vecd &disp_transform, Real inv_h, auto &src_h_ratio)->Real
            { return src_h_ratio.KernelTransform(i) * W(disp_transform, inv_h); },
            [&](const Vecd &disp_transform, Real inv_h, auto &tar_h_ratio)->Real
            { return tar_h_ratio.KernelTransform(j) * W(disp_transform, inv_h); });
    }
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    dW_ij(UnsignedInt i, UnsignedInt j) const
{
    return dW(vec_r_ij(i, j), invSmoothingLength(i, j));
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    W0(UnsignedInt i, const Vec2d &zero) const
{
    Real inv_h_i = src_h_ratio_(i) * src_inv_h_ref_;
    return src_h_ratio_.KernelTransform(i) * BaseKernel::W02D(math::pow(inv_h_i, 2));
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    W0(UnsignedInt i, const Vec3d &zero) const
{
    Real inv_h_i = src_h_ratio_(i) * src_inv_h_ref_;
    return src_h_ratio_.KernelTransform(i) * BaseKernel::W03D(math::pow(inv_h_i, 3));
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    W(const Vec2d &disp_transform, Real inv_h) const
{
    return BaseKernel::W2D(math::pow(inv_h, 2), disp_transform.norm() * inv_h);
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    W(const Vec3d &disp_transform, Real inv_h) const
{
    return BaseKernel::W3D(math::pow(inv_h, 3), disp_transform.norm() * inv_h);
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    dW(const Vec2d &disp_transform, Real inv_h) const
{
    return BaseKernel::dW2D(math::pow(inv_h, 3), disp_transform.norm() * inv_h);
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    dW(const Vec3d &disp_transform, Real inv_h) const
{
    return BaseKernel::dW3D(math::pow(inv_h, 4), disp_transform.norm() * inv_h);
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    d2W(const Vec2d &disp_transform, Real inv_h) const
{
    return BaseKernel::d2W2D(math::pow(inv_h, 4), disp_transform.norm() * inv_h);
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
Real Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::
    d2W(const Vec3d &disp_transform, Real inv_h) const
{
    return BaseKernel::d2W3D(math::pow(inv_h, 5), disp_transform.norm() * inv_h);
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
template <typename FuncI, typename FuncJ>
auto Neighbor<SourceAdaptationType, TargetAdaptationType>::SmoothingKernel::selectKernelFunction(
    UnsignedInt i, UnsignedInt j, const Vecd &disp,
    const SourceSmoothingLengthRatio &src_h_ratio, const TargetSmoothingLengthRatio &tar_h_ratio,
    const FuncI &func_i, const FuncJ &func_j) const
{
    Real inv_h_i = src_h_ratio_(i) * src_inv_h_ref_;
    Real inv_h_j = tar_h_ratio_(j) * tar_inv_h_ref_;
    const Vecd &disp_i = src_h_ratio_.transform(disp, i);
    const Vecd &disp_j = tar_h_ratio_.transform(disp, j);
    return disp_i.squaredNorm() * inv_h_i * inv_h_i < disp_j.squaredNorm() * inv_h_j * inv_h_j
               ? func_i(disp_i, inv_h_i, src_h_ratio)
               : func_j(disp_j, inv_h_j, tar_h_ratio);
}
//=================================================================================================//
template <class SourceAdaptationType, class TargetAdaptationType>
template <
    typename U, typename V,
    std::enable_if_t<
        !std::is_base_of_v<AnisotropicAdaptation, U> && !std::is_base_of_v<AnisotropicAdaptation, V>, int>>
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
