#ifndef ADAPTATION_HPP
#define ADAPTATION_HPP

#include "adaptation.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
AdaptiveNearSurface::LocalSpacing::ComputingKernel::
    ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : smoothing_kerel_(*encloser.kernel_ptr_),
      signed_distance_(ex_policy, encloser.level_set_, "LevelSet"),
      spacing_ref_(encloser.spacing_ref_), finest_spacing_bound_(encloser.finest_spacing_bound_),
      coarsest_spacing_bound_(encloser.coarsest_spacing_bound_),
      kernel_size_(smoothing_kerel_.KernelSize()), inv_w0_(1.0 / smoothing_kerel_.normalized_W(0)) {}
//=================================================================================================//
inline Real AdaptiveNearSurface::LocalSpacing::ComputingKernel::operator()(const Vecd &position)
{
    Real ratio_ref = math::fabs(signed_distance_(position)) / (2.0 * spacing_ref_);
    Real target_spacing = coarsest_spacing_bound_;
    if (ratio_ref < kernel_size_)
    {
        Real weight = smoothing_kerel_.normalized_W(ratio_ref) * inv_w0_;
        target_spacing = weight * finest_spacing_bound_ + (1.0 - weight) * coarsest_spacing_bound_;
    }
    return target_spacing;
}
//=================================================================================================//
} // namespace SPH
#endif // ADAPTATION_HPP