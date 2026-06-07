#ifndef APHI_BLOCK_ZERO_CK_HPP
#define APHI_BLOCK_ZERO_CK_HPP

#include "electromagnetic_dynamics/aphi_block_zero_ck.h"

namespace SPH
{
namespace electromagnetics
{

inline AphiZeroBlockCK::AphiZeroBlockCK(SPHBody &sph_body, const AphiBlockNames &block_names)
    : LocalDynamics(sph_body),
      dv_a_real_(particles_->template getVariableByName<Vecd>(block_names.a_real)),
      dv_a_imag_(particles_->template getVariableByName<Vecd>(block_names.a_imag)),
      dv_phi_real_(particles_->template getVariableByName<Real>(block_names.phi_real)),
      dv_phi_imag_(particles_->template getVariableByName<Real>(block_names.phi_imag))
{
}

template <class ExecutionPolicy, class EncloserType>
inline AphiZeroBlockCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : a_real_(encloser.dv_a_real_->DelegatedData(ex_policy)),
      a_imag_(encloser.dv_a_imag_->DelegatedData(ex_policy)),
      phi_real_(encloser.dv_phi_real_->DelegatedData(ex_policy)),
      phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy))
{
}

inline void AphiZeroBlockCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    a_real_[index_i] = Vecd::Zero();
    a_imag_[index_i] = Vecd::Zero();
    phi_real_[index_i] = 0.0;
    phi_imag_[index_i] = 0.0;
}

inline AphiCopyBlockCK::AphiCopyBlockCK(SPHBody &sph_body, const AphiBlockNames &dst_names, const AphiBlockNames &src_names)
    : LocalDynamics(sph_body),
      dv_dst_a_real_(particles_->template getVariableByName<Vecd>(dst_names.a_real)),
      dv_dst_a_imag_(particles_->template getVariableByName<Vecd>(dst_names.a_imag)),
      dv_dst_phi_real_(particles_->template getVariableByName<Real>(dst_names.phi_real)),
      dv_dst_phi_imag_(particles_->template getVariableByName<Real>(dst_names.phi_imag)),
      dv_src_a_real_(particles_->template getVariableByName<Vecd>(src_names.a_real)),
      dv_src_a_imag_(particles_->template getVariableByName<Vecd>(src_names.a_imag)),
      dv_src_phi_real_(particles_->template getVariableByName<Real>(src_names.phi_real)),
      dv_src_phi_imag_(particles_->template getVariableByName<Real>(src_names.phi_imag))
{
}

template <class ExecutionPolicy, class EncloserType>
inline AphiCopyBlockCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : dst_a_real_(encloser.dv_dst_a_real_->DelegatedData(ex_policy)),
      dst_a_imag_(encloser.dv_dst_a_imag_->DelegatedData(ex_policy)),
      dst_phi_real_(encloser.dv_dst_phi_real_->DelegatedData(ex_policy)),
      dst_phi_imag_(encloser.dv_dst_phi_imag_->DelegatedData(ex_policy)),
      src_a_real_(encloser.dv_src_a_real_->DelegatedData(ex_policy)),
      src_a_imag_(encloser.dv_src_a_imag_->DelegatedData(ex_policy)),
      src_phi_real_(encloser.dv_src_phi_real_->DelegatedData(ex_policy)),
      src_phi_imag_(encloser.dv_src_phi_imag_->DelegatedData(ex_policy))
{
}

inline void AphiCopyBlockCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    dst_a_real_[index_i] = src_a_real_[index_i];
    dst_a_imag_[index_i] = src_a_imag_[index_i];
    dst_phi_real_[index_i] = src_phi_real_[index_i];
    dst_phi_imag_[index_i] = src_phi_imag_[index_i];
}

inline AphiComputeResidualCK::AphiComputeResidualCK(SPHBody &sph_body, const AphiVariableNames &variable_names)
    : LocalDynamics(sph_body),
      dv_residual_a_real_(particles_->template getVariableByName<Vecd>(variable_names.residual.a_real)),
      dv_residual_a_imag_(particles_->template getVariableByName<Vecd>(variable_names.residual.a_imag)),
      dv_residual_phi_real_(particles_->template getVariableByName<Real>(variable_names.residual.phi_real)),
      dv_residual_phi_imag_(particles_->template getVariableByName<Real>(variable_names.residual.phi_imag)),
      dv_rhs_a_real_(particles_->template getVariableByName<Vecd>(variable_names.rhs.a_real)),
      dv_rhs_a_imag_(particles_->template getVariableByName<Vecd>(variable_names.rhs.a_imag)),
      dv_rhs_phi_real_(particles_->template getVariableByName<Real>(variable_names.rhs.phi_real)),
      dv_rhs_phi_imag_(particles_->template getVariableByName<Real>(variable_names.rhs.phi_imag)),
      dv_lhs_a_real_(particles_->template getVariableByName<Vecd>(variable_names.lhs.a_real)),
      dv_lhs_a_imag_(particles_->template getVariableByName<Vecd>(variable_names.lhs.a_imag)),
      dv_lhs_phi_real_(particles_->template getVariableByName<Real>(variable_names.lhs.phi_real)),
      dv_lhs_phi_imag_(particles_->template getVariableByName<Real>(variable_names.lhs.phi_imag))
{
}

template <class ExecutionPolicy, class EncloserType>
inline AphiComputeResidualCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : residual_a_real_(encloser.dv_residual_a_real_->DelegatedData(ex_policy)),
      residual_a_imag_(encloser.dv_residual_a_imag_->DelegatedData(ex_policy)),
      residual_phi_real_(encloser.dv_residual_phi_real_->DelegatedData(ex_policy)),
      residual_phi_imag_(encloser.dv_residual_phi_imag_->DelegatedData(ex_policy)),
      rhs_a_real_(encloser.dv_rhs_a_real_->DelegatedData(ex_policy)),
      rhs_a_imag_(encloser.dv_rhs_a_imag_->DelegatedData(ex_policy)),
      rhs_phi_real_(encloser.dv_rhs_phi_real_->DelegatedData(ex_policy)),
      rhs_phi_imag_(encloser.dv_rhs_phi_imag_->DelegatedData(ex_policy)),
      lhs_a_real_(encloser.dv_lhs_a_real_->DelegatedData(ex_policy)),
      lhs_a_imag_(encloser.dv_lhs_a_imag_->DelegatedData(ex_policy)),
      lhs_phi_real_(encloser.dv_lhs_phi_real_->DelegatedData(ex_policy)),
      lhs_phi_imag_(encloser.dv_lhs_phi_imag_->DelegatedData(ex_policy))
{
}

inline void AphiComputeResidualCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    residual_a_real_[index_i] = rhs_a_real_[index_i] - lhs_a_real_[index_i];
    residual_a_imag_[index_i] = rhs_a_imag_[index_i] - lhs_a_imag_[index_i];
    residual_phi_real_[index_i] = rhs_phi_real_[index_i] - lhs_phi_real_[index_i];
    residual_phi_imag_[index_i] = rhs_phi_imag_[index_i] - lhs_phi_imag_[index_i];
}

inline AphiComputeBlockResidualCK::AphiComputeBlockResidualCK(SPHBody &sph_body, const AphiBlockNames &output_names,
                                                             const AphiBlockNames &rhs_names,
                                                             const AphiBlockNames &lhs_names)
    : LocalDynamics(sph_body),
      dv_output_a_real_(particles_->template getVariableByName<Vecd>(output_names.a_real)),
      dv_output_a_imag_(particles_->template getVariableByName<Vecd>(output_names.a_imag)),
      dv_output_phi_real_(particles_->template getVariableByName<Real>(output_names.phi_real)),
      dv_output_phi_imag_(particles_->template getVariableByName<Real>(output_names.phi_imag)),
      dv_rhs_a_real_(particles_->template getVariableByName<Vecd>(rhs_names.a_real)),
      dv_rhs_a_imag_(particles_->template getVariableByName<Vecd>(rhs_names.a_imag)),
      dv_rhs_phi_real_(particles_->template getVariableByName<Real>(rhs_names.phi_real)),
      dv_rhs_phi_imag_(particles_->template getVariableByName<Real>(rhs_names.phi_imag)),
      dv_lhs_a_real_(particles_->template getVariableByName<Vecd>(lhs_names.a_real)),
      dv_lhs_a_imag_(particles_->template getVariableByName<Vecd>(lhs_names.a_imag)),
      dv_lhs_phi_real_(particles_->template getVariableByName<Real>(lhs_names.phi_real)),
      dv_lhs_phi_imag_(particles_->template getVariableByName<Real>(lhs_names.phi_imag))
{
}

template <class ExecutionPolicy, class EncloserType>
inline AphiComputeBlockResidualCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : output_a_real_(encloser.dv_output_a_real_->DelegatedData(ex_policy)),
      output_a_imag_(encloser.dv_output_a_imag_->DelegatedData(ex_policy)),
      output_phi_real_(encloser.dv_output_phi_real_->DelegatedData(ex_policy)),
      output_phi_imag_(encloser.dv_output_phi_imag_->DelegatedData(ex_policy)),
      rhs_a_real_(encloser.dv_rhs_a_real_->DelegatedData(ex_policy)),
      rhs_a_imag_(encloser.dv_rhs_a_imag_->DelegatedData(ex_policy)),
      rhs_phi_real_(encloser.dv_rhs_phi_real_->DelegatedData(ex_policy)),
      rhs_phi_imag_(encloser.dv_rhs_phi_imag_->DelegatedData(ex_policy)),
      lhs_a_real_(encloser.dv_lhs_a_real_->DelegatedData(ex_policy)),
      lhs_a_imag_(encloser.dv_lhs_a_imag_->DelegatedData(ex_policy)),
      lhs_phi_real_(encloser.dv_lhs_phi_real_->DelegatedData(ex_policy)),
      lhs_phi_imag_(encloser.dv_lhs_phi_imag_->DelegatedData(ex_policy))
{
}

inline void AphiComputeBlockResidualCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    output_a_real_[index_i] = rhs_a_real_[index_i] - lhs_a_real_[index_i];
    output_a_imag_[index_i] = rhs_a_imag_[index_i] - lhs_a_imag_[index_i];
    output_phi_real_[index_i] = rhs_phi_real_[index_i] - lhs_phi_real_[index_i];
    output_phi_imag_[index_i] = rhs_phi_imag_[index_i] - lhs_phi_imag_[index_i];
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_BLOCK_ZERO_CK_HPP
