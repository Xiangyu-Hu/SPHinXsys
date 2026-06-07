#ifndef APHI_COLD_CRUCIBLE_CASE_CK_HPP
#define APHI_COLD_CRUCIBLE_CASE_CK_HPP

#include "electromagnetic_dynamics/benchmark/aphi_cold_crucible_case_ck.h"

namespace SPH
{
namespace electromagnetics
{
namespace benchmark
{

inline AssignColdCrucibleRegionMaterialsCK::AssignColdCrucibleRegionMaterialsCK(
    SPHBody &sph_body, const AphiColdCrucibleUnitBoxLayout &layout, const AphiMaterialNames &material_names)
    : LocalDynamics(sph_body),
      dv_position_(particles_->template getVariableByName<Vecd>("Position")),
      dv_sigma_(particles_->template getVariableByName<Real>(material_names.sigma)),
      dv_nu_(particles_->template getVariableByName<Real>(material_names.nu)), layout_(layout)
{
}

template <class ExecutionPolicy, class EncloserType>
inline AssignColdCrucibleRegionMaterialsCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                       EncloserType &encloser)
    : position_(encloser.dv_position_->DelegatedData(ex_policy)),
      sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)), nu_(encloser.dv_nu_->DelegatedData(ex_policy)),
      layout_(encloser.layout_)
{
}

inline void AssignColdCrucibleRegionMaterialsCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    const Vecd &position = position_[index_i];
    if (insideBoxRegion(position, layout_.melt))
    {
        sigma_[index_i] = layout_.melt_material.sigma;
        nu_[index_i] = layout_.melt_material.nu;
        return;
    }
    if (insideBoxRegion(position, layout_.crucible_wall))
    {
        sigma_[index_i] = layout_.crucible_wall_material.sigma;
        nu_[index_i] = layout_.crucible_wall_material.nu;
        return;
    }
    if (insideAnyCoilRegion(position, layout_))
    {
        sigma_[index_i] = layout_.coil_material.sigma;
        nu_[index_i] = layout_.coil_material.nu;
        return;
    }
    sigma_[index_i] = layout_.air.sigma;
    nu_[index_i] = layout_.air.nu;
}

inline AssignColdCrucibleCoilSourceCK::AssignColdCrucibleCoilSourceCK(
    SPHBody &sph_body, const AphiBlockNames &rhs_block, const AphiColdCrucibleUnitBoxLayout &layout,
    const Vecd &current_real, const Vecd &current_imag, Real peak_amplitude)
    : LocalDynamics(sph_body),
      dv_position_(particles_->template getVariableByName<Vecd>("Position")),
      dv_rhs_a_real_(particles_->template getVariableByName<Vecd>(rhs_block.a_real)),
      dv_rhs_a_imag_(particles_->template getVariableByName<Vecd>(rhs_block.a_imag)), layout_(layout),
      current_real_(current_real), current_imag_(current_imag), peak_amplitude_(peak_amplitude)
{
}

template <class ExecutionPolicy, class EncloserType>
inline AssignColdCrucibleCoilSourceCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                  EncloserType &encloser)
    : position_(encloser.dv_position_->DelegatedData(ex_policy)),
      rhs_a_real_(encloser.dv_rhs_a_real_->DelegatedData(ex_policy)),
      rhs_a_imag_(encloser.dv_rhs_a_imag_->DelegatedData(ex_policy)), layout_(encloser.layout_),
      current_real_(encloser.current_real_), current_imag_(encloser.current_imag_),
      peak_amplitude_(encloser.peak_amplitude_)
{
}

inline void AssignColdCrucibleCoilSourceCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    const Vecd &position = position_[index_i];
    const Real profile_left = smoothBoxProfile(position, layout_.coil_left, peak_amplitude_);
    const Real profile_right = smoothBoxProfile(position, layout_.coil_right, peak_amplitude_);
    const Real profile = profile_left + profile_right;
    rhs_a_real_[index_i] = profile * current_real_;
    rhs_a_imag_[index_i] = profile * current_imag_;
}

} // namespace benchmark
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_COLD_CRUCIBLE_CASE_CK_HPP
