#ifndef APHI_BENCHMARK_CASE_CK_HPP
#define APHI_BENCHMARK_CASE_CK_HPP

#include "electromagnetic_dynamics/benchmark/aphi_benchmark_case_ck.h"

namespace SPH
{
namespace electromagnetics
{
namespace benchmark
{

inline AssignPiecewiseSigmaHalfSpaceCK::AssignPiecewiseSigmaHalfSpaceCK(SPHBody &sph_body, Real x_interface,
                                                                        Real sigma_left, Real sigma_right, Real nu,
                                                                        const AphiMaterialNames &material_names)
    : LocalDynamics(sph_body),
      dv_position_(particles_->template getVariableByName<Vecd>("Position")),
      dv_sigma_(particles_->template getVariableByName<Real>(material_names.sigma)),
      dv_nu_(particles_->template getVariableByName<Real>(material_names.nu)),
      x_interface_(x_interface), sigma_left_(sigma_left), sigma_right_(sigma_right), nu_value_(nu)
{
}

template <class ExecutionPolicy, class EncloserType>
inline AssignPiecewiseSigmaHalfSpaceCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                   EncloserType &encloser)
    : position_(encloser.dv_position_->DelegatedData(ex_policy)),
      sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)), nu_(encloser.dv_nu_->DelegatedData(ex_policy)),
      x_interface_(encloser.x_interface_), sigma_left_(encloser.sigma_left_), sigma_right_(encloser.sigma_right_),
      nu_value_(encloser.nu_value_)
{
}

inline void AssignPiecewiseSigmaHalfSpaceCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    sigma_[index_i] = position_[index_i][0] <= x_interface_ ? sigma_left_ : sigma_right_;
    nu_[index_i] = nu_value_;
}

inline AssignTeam7LikeRegionMaterialsCK::AssignTeam7LikeRegionMaterialsCK(SPHBody &sph_body,
                                                                          const AphiTeam7LikeUnitBoxLayout &layout,
                                                                          const AphiMaterialNames &material_names)
    : LocalDynamics(sph_body),
      dv_position_(particles_->template getVariableByName<Vecd>("Position")),
      dv_sigma_(particles_->template getVariableByName<Real>(material_names.sigma)),
      dv_nu_(particles_->template getVariableByName<Real>(material_names.nu)), layout_(layout)
{
}

template <class ExecutionPolicy, class EncloserType>
inline AssignTeam7LikeRegionMaterialsCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                    EncloserType &encloser)
    : position_(encloser.dv_position_->DelegatedData(ex_policy)),
      sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)), nu_(encloser.dv_nu_->DelegatedData(ex_policy)),
      layout_(encloser.layout_)
{
}

inline void AssignTeam7LikeRegionMaterialsCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    const Vecd &position = position_[index_i];
    if (insideBoxRegion(position, layout_.conductor))
    {
        sigma_[index_i] = layout_.conductor_material.sigma;
        nu_[index_i] = layout_.conductor_material.nu;
        return;
    }
    if (insideBoxRegion(position, layout_.coil))
    {
        sigma_[index_i] = layout_.coil_material.sigma;
        nu_[index_i] = layout_.coil_material.nu;
        return;
    }
    sigma_[index_i] = layout_.air.sigma;
    nu_[index_i] = layout_.air.nu;
}

inline AssignImpressedCurrentRhsCK::AssignImpressedCurrentRhsCK(SPHBody &sph_body, const AphiBlockNames &rhs_block,
                                                                const AphiBoxRegion &source_region,
                                                                const Vecd &current_real, const Vecd &current_imag,
                                                                Real peak_amplitude)
    : LocalDynamics(sph_body),
      dv_position_(particles_->template getVariableByName<Vecd>("Position")),
      dv_rhs_a_real_(particles_->template getVariableByName<Vecd>(rhs_block.a_real)),
      dv_rhs_a_imag_(particles_->template getVariableByName<Vecd>(rhs_block.a_imag)), source_region_(source_region),
      current_real_(current_real), current_imag_(current_imag), peak_amplitude_(peak_amplitude)
{
}

template <class ExecutionPolicy, class EncloserType>
inline AssignImpressedCurrentRhsCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : position_(encloser.dv_position_->DelegatedData(ex_policy)),
      rhs_a_real_(encloser.dv_rhs_a_real_->DelegatedData(ex_policy)),
      rhs_a_imag_(encloser.dv_rhs_a_imag_->DelegatedData(ex_policy)), source_region_(encloser.source_region_),
      current_real_(encloser.current_real_), current_imag_(encloser.current_imag_),
      peak_amplitude_(encloser.peak_amplitude_)
{
}

inline void AssignImpressedCurrentRhsCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    const Real profile = smoothBoxProfile(position_[index_i], source_region_, peak_amplitude_);
    rhs_a_real_[index_i] = profile * current_real_;
    rhs_a_imag_[index_i] = profile * current_imag_;
}

inline AssignZeroSigmaInAnnularRegionCK::AssignZeroSigmaInAnnularRegionCK(SPHBody &sph_body,
                                                                          const AphiAnnularSourceGeometry &geometry,
                                                                          const AphiMaterialNames &material_names)
    : LocalDynamics(sph_body), geometry_(geometry),
      dv_position_(particles_->template getVariableByName<Vecd>("Position")),
      dv_sigma_(particles_->template getVariableByName<Real>(material_names.sigma))
{
}

template <class ExecutionPolicy, class EncloserType>
inline AssignZeroSigmaInAnnularRegionCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                    EncloserType &encloser)
    : position_(encloser.dv_position_->DelegatedData(ex_policy)),
      sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)), geometry_(encloser.geometry_)
{
}

inline void AssignZeroSigmaInAnnularRegionCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    if (insideAnnularSourceRegion(position_[index_i], geometry_))
    {
        sigma_[index_i] = 0.0;
    }
}

inline AssignImpressedCurrentRhsAnnularCK::AssignImpressedCurrentRhsAnnularCK(
    SPHBody &sph_body, const AphiBlockNames &rhs_block, const AphiAnnularSourceGeometry &geometry,
    const Vecd &current_real, const Vecd &current_imag, Real peak_amplitude)
    : LocalDynamics(sph_body), geometry_(geometry),
      dv_position_(particles_->template getVariableByName<Vecd>("Position")),
      dv_rhs_a_real_(particles_->template getVariableByName<Vecd>(rhs_block.a_real)),
      dv_rhs_a_imag_(particles_->template getVariableByName<Vecd>(rhs_block.a_imag)), current_real_(current_real),
      current_imag_(current_imag), peak_amplitude_(peak_amplitude)
{
}

template <class ExecutionPolicy, class EncloserType>
inline AssignImpressedCurrentRhsAnnularCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                      EncloserType &encloser)
    : position_(encloser.dv_position_->DelegatedData(ex_policy)),
      rhs_a_real_(encloser.dv_rhs_a_real_->DelegatedData(ex_policy)),
      rhs_a_imag_(encloser.dv_rhs_a_imag_->DelegatedData(ex_policy)), geometry_(encloser.geometry_),
      current_real_(encloser.current_real_), current_imag_(encloser.current_imag_),
      peak_amplitude_(encloser.peak_amplitude_)
{
}

inline void AssignImpressedCurrentRhsAnnularCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    if (!insideAnnularSourceRegion(position_[index_i], geometry_))
    {
        rhs_a_real_[index_i] = Vecd(0.0, 0.0, 0.0);
        rhs_a_imag_[index_i] = Vecd(0.0, 0.0, 0.0);
        return;
    }
    rhs_a_real_[index_i] = peak_amplitude_ * current_real_;
    rhs_a_imag_[index_i] = peak_amplitude_ * current_imag_;
}

inline AssignSolenoidalCurlCurrentRhsCK::AssignSolenoidalCurlCurrentRhsCK(
    SPHBody &sph_body, const AphiBlockNames &rhs_block, const AphiBoxRegion &source_region, Real peak_amplitude,
    Real imag_to_real_ratio)
    : LocalDynamics(sph_body),
      dv_position_(particles_->template getVariableByName<Vecd>("Position")),
      dv_rhs_a_real_(particles_->template getVariableByName<Vecd>(rhs_block.a_real)),
      dv_rhs_a_imag_(particles_->template getVariableByName<Vecd>(rhs_block.a_imag)), source_region_(source_region),
      peak_amplitude_(peak_amplitude), imag_to_real_ratio_(imag_to_real_ratio)
{
}

template <class ExecutionPolicy, class EncloserType>
inline AssignSolenoidalCurlCurrentRhsCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                   EncloserType &encloser)
    : position_(encloser.dv_position_->DelegatedData(ex_policy)),
      rhs_a_real_(encloser.dv_rhs_a_real_->DelegatedData(ex_policy)),
      rhs_a_imag_(encloser.dv_rhs_a_imag_->DelegatedData(ex_policy)), source_region_(encloser.source_region_),
      peak_amplitude_(encloser.peak_amplitude_), imag_to_real_ratio_(encloser.imag_to_real_ratio_)
{
}

inline void AssignSolenoidalCurlCurrentRhsCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    const Vecd &position = position_[index_i];
    Real g = 0.0;
    Real dgdx = 0.0;
    Real dgdy = 0.0;
    Real dgdz = 0.0;
    smoothBoxProfileWithGradients(position, source_region_, peak_amplitude_, g, dgdx, dgdy, dgdz);
    (void)dgdz;

    const Real x = position[0];
    const Real y = position[1];
    const Real radius_squared = x * x + y * y;
    const Vecd current_real(-0.5 * (dgdy * radius_squared + 2.0 * g * y),
                            0.5 * (dgdx * radius_squared + 2.0 * g * x), 0.0);
    rhs_a_real_[index_i] = current_real;
    rhs_a_imag_[index_i] = imag_to_real_ratio_ * current_real;
}

inline AssignInterfaceFluxMatchedPhiFieldsCK::AssignInterfaceFluxMatchedPhiFieldsCK(
    SPHBody &sph_body, Real x_interface, Real sigma_left, Real sigma_right, const AphiBlockNames &block_names)
    : LocalDynamics(sph_body),
      dv_position_(particles_->template getVariableByName<Vecd>("Position")),
      dv_a_real_(particles_->template getVariableByName<Vecd>(block_names.a_real)),
      dv_a_imag_(particles_->template getVariableByName<Vecd>(block_names.a_imag)),
      dv_phi_real_(particles_->template getVariableByName<Real>(block_names.phi_real)),
      dv_phi_imag_(particles_->template getVariableByName<Real>(block_names.phi_imag)), x_interface_(x_interface),
      slope_left_(Real(1)), slope_right_(sigma_left / (sigma_right + TinyReal)),
      offset_right_(slope_left_ * x_interface - slope_right_ * x_interface)
{
}

template <class ExecutionPolicy, class EncloserType>
inline AssignInterfaceFluxMatchedPhiFieldsCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                         EncloserType &encloser)
    : position_(encloser.dv_position_->DelegatedData(ex_policy)),
      a_real_(encloser.dv_a_real_->DelegatedData(ex_policy)), a_imag_(encloser.dv_a_imag_->DelegatedData(ex_policy)),
      phi_real_(encloser.dv_phi_real_->DelegatedData(ex_policy)),
      phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy)), x_interface_(encloser.x_interface_),
      slope_left_(encloser.slope_left_), slope_right_(encloser.slope_right_), offset_right_(encloser.offset_right_)
{
}

inline void AssignInterfaceFluxMatchedPhiFieldsCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    const Real x = position_[index_i][0];
    a_real_[index_i] = Vecd(0.0, 0.0, 0.0);
    a_imag_[index_i] = Vecd(0.0, 0.0, 0.0);
    phi_imag_[index_i] = Real(0);
    if (x <= x_interface_)
    {
        phi_real_[index_i] = slope_left_ * x;
        return;
    }
    phi_real_[index_i] = slope_right_ * x + offset_right_;
}

inline AssignUniformLinearPhiFieldsCK::AssignUniformLinearPhiFieldsCK(SPHBody &sph_body, Real electric_field_x,
                                                                      const AphiBlockNames &block_names)
    : LocalDynamics(sph_body),
      dv_position_(particles_->template getVariableByName<Vecd>("Position")),
      dv_a_real_(particles_->template getVariableByName<Vecd>(block_names.a_real)),
      dv_a_imag_(particles_->template getVariableByName<Vecd>(block_names.a_imag)),
      dv_phi_real_(particles_->template getVariableByName<Real>(block_names.phi_real)),
      dv_phi_imag_(particles_->template getVariableByName<Real>(block_names.phi_imag)), electric_field_x_(electric_field_x)
{
}

template <class ExecutionPolicy, class EncloserType>
inline AssignUniformLinearPhiFieldsCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                EncloserType &encloser)
    : position_(encloser.dv_position_->DelegatedData(ex_policy)),
      a_real_(encloser.dv_a_real_->DelegatedData(ex_policy)), a_imag_(encloser.dv_a_imag_->DelegatedData(ex_policy)),
      phi_real_(encloser.dv_phi_real_->DelegatedData(ex_policy)),
      phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy)), electric_field_x_(encloser.electric_field_x_)
{
}

inline void AssignUniformLinearPhiFieldsCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    a_real_[index_i] = Vecd(0.0, 0.0, 0.0);
    a_imag_[index_i] = Vecd(0.0, 0.0, 0.0);
    phi_imag_[index_i] = Real(0);
    phi_real_[index_i] = -electric_field_x_ * position_[index_i][0];
}

} // namespace benchmark
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_BENCHMARK_CASE_CK_HPP
