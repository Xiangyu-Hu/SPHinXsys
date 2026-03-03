#ifndef ADAPTATION_HPP
#define ADAPTATION_HPP

#include "adaptation.h"

namespace SPH
{
//=================================================================================================//
inline Real AdaptiveSmoothingLength::SmoothedSpacing::operator()(
    const Real &measure, const Real &transition_thickness)
{
    Real ratio_ref = measure / (2.0 * transition_thickness);
    Real target_spacing = coarsest_spacing_bound_;
    if (ratio_ref < kernel_size_)
    {
        Real weight = smoothing_kernel_.normalized_W(ratio_ref) * inv_w0_;
        target_spacing = weight * finest_spacing_bound_ + (1.0 - weight) * coarsest_spacing_bound_;
    }
    return target_spacing;
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
AdaptiveNearSurface::LocalSpacing::ComputingKernel::
    ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : smoothed_spacing_(encloser.smoothed_spacing_),
      signed_distance_(ex_policy, encloser.level_set_, "LevelSet"),
      spacing_ref_(encloser.spacing_ref_) {}
//=================================================================================================//
template <typename... Args>
AdaptiveNearSurface::AdaptiveNearSurface(Args &&...args)
    : AdaptiveByShape(std::forward<Args>(args)...) {}
//=================================================================================================//
inline Real AdaptiveNearSurface::LocalSpacing::ComputingKernel::operator()(const Vecd &position)
{
    return smoothed_spacing_(math::fabs(signed_distance_(position)), spacing_ref_);
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
AdaptiveWithinShape::LocalSpacing::ComputingKernel::
    ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : smoothed_spacing_(encloser.smoothed_spacing_),
      signed_distance_(ex_policy, encloser.level_set_, "LevelSet"),
      spacing_ref_(encloser.spacing_ref_) {}
//=================================================================================================//
inline Real AdaptiveWithinShape::LocalSpacing::ComputingKernel::operator()(const Vecd &position)
{
    Real phi = signed_distance_(position);
    return phi < 0.0 ? smoothed_spacing_.FinestSpacingBound() : smoothed_spacing_(phi, 2.0 * spacing_ref_);
}
//=================================================================================================//
template <class ExecutionPolicy>
AnisotropicAdaptation::AnisotropicSmoothingLengthRatio::
    AnisotropicSmoothingLengthRatio(const ExecutionPolicy &ex_policy, AnisotropicAdaptation &adaptation)
    : ContinuousSmoothingLengthRatio(adaptation),
      deformation_matrix_(adaptation.dv_deformation_matrix_->DelegatedData(ex_policy)),
      deformation_det_(adaptation.dv_deformation_det_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline Vecd AnisotropicAdaptation::AnisotropicSmoothingLengthRatio::
    transform(const Vecd &original, UnsignedInt index_i) const
{
    return deformation_matrix_[index_i] * original;
}
//=================================================================================================//
inline Vecd AnisotropicAdaptation::AnisotropicSmoothingLengthRatio::
    inverseTransform(const Vecd &original, UnsignedInt index_i) const
{
    return deformation_matrix_[index_i].inverse() * original;
}
//=================================================================================================//
inline Real AnisotropicAdaptation::AnisotropicSmoothingLengthRatio::
    KernelTransform(UnsignedInt index_i) const
{
    return deformation_det_[index_i];
}
//=================================================================================================//
inline Matd AnisotropicAdaptation::AnisotropicSmoothingLengthRatio::
    GradientTransform(UnsignedInt index_i) const
{
    return deformation_matrix_[index_i];
}
//=================================================================================================//
} // namespace SPH
#endif // ADAPTATION_HPP